import parameters as param
import logging
import numpy as np
import pandas as pd
import time
import tqdm
from brainpy import isotopic_variants
from parameters import MASS_H
from pyteomics import mzml
from util_functions import (unit_process_spectrum, get_averagine_mass, get_averagine,
                            isostamp_convolution, contains_precursor, dynamic_rt_filter,
                            unit_fuzzy_match, pearson_correlation, ts_reindex_v2, reduce_by_dfs)

logging.basicConfig(level=logging.INFO)
np.seterr(all='ignore')


class IsoStampEngine(object):
    def __init__(self, worker, output_path):
        self.worker = worker
        self.output_path = output_path

    @staticmethod
    def yield_spectrum(mzml_path):
        for spectrum in mzml.read(mzml_path):
            yield spectrum

    def split_mzml(self, mzml_path):
        """
        consumes a profile mzml profile iterator into MS1 and MS2 dataframes
        :param mzml_path:
        :return:
        """
        ms1 = []
        ms2 = []
        
        for output_dict in tqdm.tqdm(self.worker.imap(unit_process_spectrum, self.yield_spectrum(mzml_path))):
            if output_dict['ms level'] == 1:
                ms1.append(output_dict)
            elif output_dict['ms level'] == 2:
                ms2.append(output_dict['precursor_info'])

        
        ms1 = pd.DataFrame(ms1).sort_values('ms1_scan_number').drop('ms level', 1)
        ms2 = pd.DataFrame(ms2, columns=['obs_precursor_mz', 'obs_precursor_charge', 'ms1_scan_number',
                                         'ms2_scan_number', 'ms2_rt'])

        ms2['obs_precursor_mz'] -= param.ISOLATION_WINDOW_OFFSET

        return ms1, ms2

    @staticmethod
    def yield_brainpy_task(ms1, ms2_reduced):
        """returns the parameters for brainpy_function to use multiprocessing. Does MUCH MORE than the yield
            name implies.

            Also filters ms1 information based on ms2 precursor information. These domain knowledge rules reduce
            false positives rate and the search space.
        """
        for _, ms2_row in ms2_reduced.iterrows():
            # filter 1: ms1 is between ms2 rt time window
            ms1_subset = (ms1[ms1['ms1_rt'].between(ms2_row['rt_start'], ms2_row['rt_end'])]
                           .sort_values('ms1_rt')
                           .copy())

            # filter 2: ms1 contains the ms2 obs_precursor_mz with an intensity floor
            try:
                precursor_mask = ms1_subset.apply(lambda x: contains_precursor(x, ms2_row['obs_precursor_mz'],
                                                                        param.PRECURSOR_MZ_WINDOW,
                                                                        param.MS1_INTENSITY_FLOOR), 1)
                ms1_subset = ms1_subset.loc[precursor_mask]
            except:
                logging.warning(f'filter for precursor_mz failed, ms2_row: {ms2_row}')

            # create a dataframe of ms1
            # index: mz, columns: rt, values: intensity
            ms1_subset_df = []
            mz_start = ms2_row['obs_precursor_mz'] - param.PRECURSOR_MZ_WINDOW
            mz_end = ms2_row['obs_precursor_mz'] + param.PRECURSOR_MZ_WINDOW
            for idx, row in ms1_subset.iterrows():
                row_mz_intensity = pd.Series(row['intensity'], index=row['mz'].round(param.PRECISION))[mz_start:mz_end]
                row_mz_intensity = row_mz_intensity.groupby(row_mz_intensity.index).sum()
                if row_mz_intensity.sum() > 0:
                    ms1_subset_df.append(row_mz_intensity)
            ms1_subset_df = pd.concat(ms1_subset_df, 1).fillna(0)

            # filter 3: rt (runtime) bounds based on a constant intensity floor
            # round rt columns to precision and groupby.
            # transpose allows for groupby of the columns as an index.
            # transpose converts a float into an object for some reason. return to float so subsequent indexing works
            ms1_subset_df.columns = ms1_subset['ms1_rt'].round(param.PRECISION)
            ms1_subset_df = ms1_subset_df.T.reset_index().groupby('ms1_rt').sum().T
            ms1_subset_df.index = ms1_subset_df.index.astype(float)
            rt_floor, rt_ceil, intensity_ts = dynamic_rt_filter(ms1_subset_df, ms2_row, param.INTENSITY_FLOOR)

            # filter 4: a precursor base peak should exist and correlate in intensity_ts
            # if not, pass zeroes to return a zero isostamp score
            precursor_base_peak = (ms2_row['obs_precursor_mz'] * ms2_row['obs_precursor_charge'] - param.MASS_H * 2) / ms2_row['obs_precursor_charge']
            base_ms2_row = ms2_row.copy()
            base_ms2_row['obs_precursor_mz'] = precursor_base_peak
            base_rt_floor, base_rt_ceil, base_intensity_ts = dynamic_rt_filter(ms1_subset_df, base_ms2_row, param.INTENSITY_FLOOR * 0.1)
            rt_intensity_ts = (pd.concat([base_intensity_ts, intensity_ts], 1)
                               .fillna(0)
                               .loc[min(rt_floor, base_rt_floor): max(rt_ceil, base_rt_ceil)])
            intensity_corr = rt_intensity_ts.corr().iloc[0, 1]

            # return zero propensity if
            if (np.isnan(intensity_corr)
                    or intensity_corr < 0.2
                    # the lack of data points above intensity is likely a false positive
                    or intensity_ts[rt_floor: rt_ceil].shape[0] <= 5
                    or base_intensity_ts[base_rt_floor: base_rt_ceil].shape[0] <= 5
                    # testing on main peak mean was too restrictive. I don't want to relax the floor/ceil calc either
                    or intensity_ts[rt_floor: rt_ceil].mean() < (param.INTENSITY_FLOOR * 0.5)
                    or base_intensity_ts[base_rt_floor: base_rt_ceil].mean() < param.INTENSITY_FLOOR * 0.1):
                ms1_subset_df = ms1_subset_df * 0

            ms1_subset_df = ms1_subset_df.loc[:, min(rt_floor, base_rt_floor):max(rt_ceil, base_rt_ceil)].sum(1)

            # prepare master key
            precursor_mz = ms2_row['obs_precursor_mz']
            precursor_charge = ms2_row['obs_precursor_charge']
            task_dict = {'idx': ms2_row['isostamp_ref_index'],
                         'precursor_mz': precursor_mz,
                         'precursor_charge': precursor_charge,
                         'ms1_mz_intensity_df': ms1_subset_df,
                         'obs_precursor_mass': (precursor_mz - MASS_H) * precursor_charge,
                         'base_peak_correlation': intensity_corr}
            yield task_dict

    @staticmethod
    def brainpy_function(task_dict):
        # unpack input
        isostamp_idx = task_dict['idx']
        precursor_mz = task_dict['precursor_mz']
        precursor_charge = task_dict['precursor_charge']
        obs_precursor_mass = task_dict['obs_precursor_mass']
        base_peak_correlation = task_dict['base_peak_correlation']

        # extract info
        pd_mode = 1. / precursor_charge
        mz_lower_bound = precursor_mz - param.PRECURSOR_MZ_WINDOW
        mz_upper_bound = precursor_mz + param.PRECURSOR_MZ_WINDOW

        data_md = ts_reindex_v2(task_dict['ms1_mz_intensity_df'], precision=1)

        # conditions that if met will return zero propensity
        if data_md.empty or data_md[mz_lower_bound:mz_upper_bound].sum() == 0:
            #logging.warning(f'{isostamp_idx} mz intensity spectrum is empty')
            return [isostamp_idx, 0, 0, 0, 0]
        match_intensity = data_md[(precursor_mz - 2.5 * pd_mode):(precursor_mz + 2.5 * pd_mode)]
        if match_intensity.sum() == 0:
            #logging.warning(f'{isostamp_idx} match_intensity is empty')
            return [isostamp_idx, 0, 0, 0, 0]

        # Brainpy predicts the "regular pattern" to compare against actual data
        fold_diff = obs_precursor_mass / get_averagine_mass()
        composition_estimate = {k: int(v * fold_diff) for k, v in get_averagine().items()}
        cluster = isotopic_variants(composition_estimate, charge=precursor_charge)
        # line profiler found that the first time accessing "cluster" takes ~16 ms per hit
        # despite that it is of type 'list'. explicitly cast to list
        cluster = list(cluster)
        regular_pattern = pd.Series({peak.mz: peak.intensity for peak in cluster}).sort_index().round(param.PRECISION)

        # remove peaks that are too low in intensity_array
        regular_pattern = regular_pattern[regular_pattern.values > (regular_pattern.max() * param.INTENSITY_MIN)]
        regular_pattern.index = regular_pattern.index + precursor_mz - regular_pattern.index[0]
        if param.ISOSTAMP is not None:
            isostamp_pattern = isostamp_convolution(regular_pattern, precursor_charge, param.ISOSTAMP['element_type'],
                                                    param.ISOSTAMP['element_count'], param.ISOSTAMP['ratio']).round(param.PRECISION)
        else:
            isostamp_pattern = regular_pattern.copy()

        # isostamp pattern
        isostamp_pattern = isostamp_pattern[isostamp_pattern.values > isostamp_pattern.max() * param.INTENSITY_MIN]
        isostamp_pattern.index = isostamp_pattern.index + precursor_mz - isostamp_pattern.index[2]

        # scaling constant of the relative maximum peak
        regular_pattern = regular_pattern * (match_intensity.max() / regular_pattern.max())
        isostamp_pattern = isostamp_pattern * (match_intensity.max() / isostamp_pattern.max())

        # get the diagnostic m/z for later scoring
        regular_diag_mz = np.concatenate([isostamp_pattern.index[0:2] - isostamp_pattern.index[2] + regular_pattern.index[0], regular_pattern.index])
        isostamp_diag_mz = isostamp_pattern.index

        # test off-by-one hypothesis shift m/z array to +/- some pd_mode(s). report the best one
        shift_array = np.arange(-3, 4)
        regular_propensity_array = pd.Series(0, index=shift_array)
        isostamp_propensity_array = regular_propensity_array.copy()
        for shift in regular_propensity_array.index:
            regular_mz_array = np.array(regular_pattern.index + shift * pd_mode * MASS_H)
            regular_mask = ((regular_mz_array >= mz_lower_bound)
                            & (regular_mz_array <= mz_upper_bound)
                            & (regular_pattern.values > 0))
            regular_mz_int = pd.Series(regular_pattern.values, index=regular_mz_array)[regular_mask]
            regular_diag_mz_shifted = regular_diag_mz + shift * pd_mode * MASS_H
            regular_propensity_array.loc[shift] = pearson_correlation(ts_reindex_v2(regular_mz_int, 1),
                                                               data_md, regular_diag_mz_shifted)

            isostamp_mz_array = np.array(isostamp_pattern.index + shift * pd_mode * MASS_H)
            isostamp_mask = ((isostamp_mz_array >= mz_lower_bound)
                             & (isostamp_mz_array <= mz_upper_bound)
                             & (isostamp_pattern.values > 0))
            isostamp_mz_int = pd.Series(isostamp_pattern.values, index=isostamp_mz_array)[isostamp_mask]
            isostamp_diag_mz_shifted = isostamp_diag_mz + shift * pd_mode * MASS_H
            isostamp_propensity_array.loc[shift] = pearson_correlation(ts_reindex_v2(isostamp_mz_int, 1),
                                                              data_md, isostamp_diag_mz_shifted)

        if regular_propensity_array.empty:
            regular_propensity = 0
            regular_best_shift = 0
        elif regular_propensity_array.max() > param.PRECURSOR_PROPENSITY_FLOOR:
            regular_propensity = regular_propensity_array.max()
            regular_best_shift = regular_propensity_array.idxmax()
        else:
            regular_propensity = regular_propensity_array.max()
            regular_best_shift = 0
        regular_precursor_mass = (regular_pattern.index[0] - MASS_H + regular_best_shift * pd_mode * MASS_H) * precursor_charge

        if isostamp_propensity_array.empty:
            isostamp_propensity = 0
            isostamp_best_shift = 0
        elif isostamp_propensity_array.max() > param.PRECURSOR_PROPENSITY_FLOOR:
            isostamp_propensity = isostamp_propensity_array.max()
            isostamp_best_shift = isostamp_propensity_array.idxmax()
        else:
            isostamp_propensity = isostamp_propensity_array.max()
            isostamp_best_shift = 0
        isostamp_precursor_mass = (isostamp_pattern.index[2] - MASS_H + isostamp_best_shift * pd_mode * MASS_H) * precursor_charge

        return [isostamp_idx, isostamp_propensity, regular_propensity, isostamp_precursor_mass,
                regular_precursor_mass, base_peak_correlation]

    def isostamp_run(self, mzml_path):
        """
        perform the search algorithm, then report the final result
        """

        logging.info('step 1: split mzml file')
        time_start = time.time()
        
        try:
            ms1, ms2 = self.split_mzml(mzml_path)
        except:
            pass
        
        
        logging.info(f'splitting mzml file took {((time.time()-time_start)/60.):.2f} min')
        logging.info(f'MS1 df shape: {ms1.shape}')
        logging.info(f'MS2 df shape: {ms2.shape}')

        logging.info('step 2: transform MS2 precursor info dataframe to unique feature')
        logging.info('MS2 precursor info into adjacency list')
        adj_dict = unit_fuzzy_match(ms2)

        logging.info('perform depth-first search')
        # I don't really know what this does to filter down the adj_dict rows
        adj_list_reduced = reduce_by_dfs(adj_dict)
        logging.info(f'first depth reduces ms2 size from {ms2.shape[0]} to {len(adj_list_reduced)}')

        # consolidate redundancy
        logging.info('step 3: perform depth-first graph search to reduce search redundancy')
        asso_list = []
        for index, s in enumerate(adj_list_reduced):
            relevant_df = ms2.loc[s]
            precursor_mz = np.median(relevant_df['obs_precursor_mz']).round(param.PRECISION)
            precursor_charge = relevant_df['obs_precursor_charge'].min()
            # account for dynamic exclusion
            rt_start = relevant_df['ms2_rt'].min() - param.MS1_RT_WINDOW
            rt_end = relevant_df['ms2_rt'].max() + param.MS1_RT_WINDOW
            rt = (rt_start + rt_end) / 2.
            ms2_scan_number = '+'.join(relevant_df['ms2_scan_number'].apply(int).apply(str))
            asso_list.append([index, precursor_mz, precursor_charge, rt_start, rt_end, rt, ms2_scan_number])

        ms2_reduced = pd.DataFrame(asso_list,
                                   columns=['isostamp_ref_index', 'obs_precursor_mz',
                                            'obs_precursor_charge', 'rt_start', 'rt_end', 'rt', 'ms2_scan_number'])

        # sort by start of retention time
        ms2_reduced.sort_values(by='rt_start', ascending=True, inplace=True)
        ms2_reduced.reset_index(inplace=True, drop=True)

        logging.info('step 4: slice MS1 by runtime (RT) and pass to multiprocess brainpy fitting')
        brainpy_result = []

        # todo: remove this debugging 20191102
        # tmp = ms2_reduced.copy()
        # ms2_reduced = ms2_reduced[ms2_reduced['ms2_scan_number'].str.contains('20224')]
        # ms2_reduced = ms2_reduced.iloc[29:31]
        for result in tqdm.tqdm(self.worker.imap(self.brainpy_function,
                                                         self.yield_brainpy_task(ms1, ms2_reduced)),
                                        total=ms2_reduced.shape[0]):
            brainpy_result.append(result)

        brainpy_df = pd.DataFrame.from_records(brainpy_result,
                                               columns=['isostamp_ref_index', 'isostamp_propensity',
                                                        'regular_propensity', 'isostamp_precursor_mass',
                                                        'regular_precursor_mass', 'base_peak_correlation'])

        ms1_df = ms2_reduced.merge(brainpy_df, on='isostamp_ref_index', how='inner')
        logging.info('running IsoStamp took {:.2f} min'.format((time.time() - time_start) / 60.))
        return ms1_df

