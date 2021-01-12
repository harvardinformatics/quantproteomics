import parameters as param
import logging
import numpy as np
import pandas as pd
import re
from scipy import stats

# this will ignore all division-related numpy warning
np.seterr(all='ignore')


def calculate_q_value(tda_list, score):
    # calculate q-value
    # preprocessing
    tda_list = np.array(tda_list)
    score = np.array(score)
    sort_key = np.argsort(score)[::-1]
    reverse_key = np.argsort(sort_key)
    sorted_tda_list = tda_list[sort_key]
    # reference doi: 10.1002/pmic.200800473
    # first step: calculate cumulative count of T and D
    count_obs = len(sorted_tda_list)
    cum_count_t_array = np.zeros((count_obs,))
    cum_count_d_array = np.zeros((count_obs,))
    cum_count_t = 0
    cum_count_d = 0
    for index in range(len(sorted_tda_list)):
        if sorted_tda_list[index] == 'target':
            cum_count_t += 1
        else:
            cum_count_d += 1
        cum_count_t_array[index] = cum_count_t
        cum_count_d_array[index] = cum_count_d
    # second step: calculate estimated FDR and q-value
    est_fdr_array = np.zeros((count_obs,))
    q_value_array = np.zeros((count_obs,))
    current_min = float('inf')
    for index in range(count_obs - 1, 0, -1):
        count_t = cum_count_t_array[index]
        count_d = cum_count_d_array[index]
        est_fdr = 2 * float(count_d) / float(count_d + count_t)
        est_fdr_array[index] = est_fdr
        if est_fdr <= current_min:
            current_min = est_fdr
        q_value_array[index] = current_min
    # return
    return q_value_array[reverse_key]


def camelcase_to_snake_case(name):
    """replaces camelcase with letters to snake case"""
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()


def unit_process_spectrum(spectrum):
    """parses a mzml spectrum"""

    # unpack input
    spectrum = {camelcase_to_snake_case(k): v for k, v in spectrum.items()}
    # container
    output_dict = dict()
    # differential processing
    if spectrum['ms level'] == 1:
        # extract info
        ms1_scan_number = spectrum['index'] + 1
        rt = spectrum['scan_list']['scan'][0]['scan start time']

        # process info
        # roughly half of intensity values are zero
        output_dict['mz'] = spectrum['m/z array']
        output_dict['intensity'] = spectrum['intensity array']

        # store info
        output_dict['ms level'] = 1
        output_dict['ms1_scan_number'] = ms1_scan_number
        output_dict['ms1_rt'] = rt

    elif spectrum['ms level'] == 2:
        # extract info
        ms1_scan_number = int(spectrum['precursor_list']['precursor'][0]['spectrumRef'].split('=')[-1])
        ms2_scan_number = spectrum['index'] + 1
        ms2_rt = spectrum['scan_list']['scan'][0]['scan start time']
        key_dict = spectrum['precursor_list']['precursor'][0]['selectedIonList']['selectedIon'][0]
        obs_precursor_mz = key_dict['selected ion m/z']
        obs_precursor_charge = key_dict['charge state']

        # store info
        output_dict['ms level'] = 2
        output_dict['ms2_scan_number'] = ms2_scan_number
        output_dict['precursor_info'] = [obs_precursor_mz, obs_precursor_charge, ms1_scan_number,
                                         ms2_scan_number, ms2_rt]

    return output_dict

def local_rsquare(model_md, data_md):
    """prepare for normal Root mean square deviation calculation
       assumes that the index values are unique"""

    # align the two series to calculate sum of square differences
    tmp = pd.concat([data_md, model_md], 1, sort=True).fillna(0)
    ss_residual = ((tmp.iloc[:, 1] - tmp.iloc[:, 0]) ** 2).sum()
    ss_total = (data_md ** 2).sum()
    try:
        return np.clip(1 - ss_residual / ss_total, 0, 1)
    except ValueError:
        logging.warning(f"""local_rsquare calculation failed. set rsq to zero
                        model_md shape: {model_md.shape}, data_md shape: {data_md.shape}""")
        return 0


def local_rsquare_product(model_md, data_md, diag_mz):
    """
    Calculates two R-squared numbers "high" and "low". Then multiplies these two
    I don't understanding the logic behind the low_bound and high_bound calculation

    :param model_md: pd.Series with mz as index and intensity as values
    :param data_md: pd.Series with mz as index and intensity as values
    :param diag_mz: First two values of the array are a special calculation of isostamp and regular pattern. Values
                    and after are the regular pattern. I'm confident that these should exist as separate object, but
                    I don't understand the conceptual meaning yet
    :return:
    """
    if model_md.empty or data_md.empty:
        return 0
    # masking to include central region
    # low bound appears to have a bug
    low_bound = (diag_mz[0] - 1.5 * (diag_mz[1] - diag_mz[0]), diag_mz[-2] + 0.5 * (diag_mz[-2] - diag_mz[-3]))
    high_bound = (diag_mz[2] - 0.5 * (diag_mz[2] - diag_mz[1]), diag_mz[-1] + 0.5 * (diag_mz[-1] - diag_mz[-2]))

    data_sp_low = data_md[low_bound[0]:low_bound[1]]
    model_sp_low = model_md[low_bound[0]:low_bound[1]]
    local_rsq_low = local_rsquare(model_sp_low, data_sp_low)

    data_sp_high = data_md[high_bound[0]:high_bound[1]]
    model_sp_high = model_md[high_bound[0]:high_bound[1]]
    local_rsq_high = local_rsquare(model_sp_high, data_sp_high)

    return round(local_rsq_low * local_rsq_high, 4)


def pearson_correlation(model_md, data_md, diag_mz):
    """
    Unsatisfied with the rsquared calculation when the model and data have different variances
    """
    from scipy.stats.stats import pearsonr
    if model_md.empty or data_md.empty:
        return 0
    # masking to include central region
    low_bound = (diag_mz[0] - 1.5 * (diag_mz[1] - diag_mz[0]), diag_mz[2] + 0.5 * (diag_mz[3] - diag_mz[2]))
    high_bound = (diag_mz[2] - 0.5 * (diag_mz[2] - diag_mz[1]), diag_mz[-1] + 0.5 * (diag_mz[-1] - diag_mz[-2]))

    sp_low = pd.concat([model_md, data_md], 1).fillna(0).loc[low_bound[0]:low_bound[1]]
    pearson_low = max(pearsonr(sp_low.iloc[:, 0], sp_low.iloc[:, 1])[0], 0)

    sp_high = pd.concat([model_md, data_md], 1).fillna(0).loc[high_bound[0]:high_bound[1]]
    pearson_high = max(pearsonr(sp_high.iloc[:, 0], sp_high.iloc[:, 1])[0], 0)

    return np.round(np.mean([pearson_low, pearson_high]), 4)


def within_tol_prob(obs, model, tol, unit, prob=False):
    """
    return 0 if it's beyond tolerance
    return Gaussian prior probability otherwise
    """
    if unit == 'ppm':
        if abs(obs - model) / obs * 1e6 > tol:
            return 0
        elif not prob:
            return 1
        else:
            sigma = model * tol * 1e-6
            normalizing_factor = stats.norm.pdf(model, loc=model, scale=sigma)
            return stats.norm.pdf(obs, loc=model, scale=sigma) / normalizing_factor
    elif unit == 'Da':
        if abs(obs - model) > tol:
            return 0
        elif not prob:
            return 1
        else:
            normalizing_factor = stats.norm.pdf(model, loc=model, scale=tol)
            return stats.norm.pdf(obs, loc=model, scale=tol) / normalizing_factor
    else:
        raise NotImplementedError('unit ' + unit + ' not recognized')


def get_averagine():
    """frequency in the human proteome"""
    return {'H': 0.0710,
            'C': 0.0435,
            'N': 0.0126,
            'O': 0.0138,
            'S': 0.00037}


def get_atomic_mass():
    """
    I don't know the original source, but these look correct from a quick search online
    https://www.lenntech.com/periodic-chart-elements/atomic-mass.htm"""
    return {'H': 1.00794,
            'C': 12.0107,
            'N': 14.0067,
            'O': 15.999,
            'S': 32.065}


def get_averagine_mass():
    """ sum product of averagine abundance in the human proteome * atomic mass"""
    averagine = get_averagine()
    atomic_mass = get_atomic_mass()
    return sum({averagine[key] * atomic_mass[key] for key in atomic_mass})


def isostamp_convolution(spectrum, charge, element_type, element_count, ratio):
    if element_type == 'H':
        mass_diff = (param.MASS_D - param.MASS_H) * element_count
    elif element_type == 'C':
        mass_diff = (param.MASS_C13 - param.MASS_C12) * element_count
    else:
        mass_diff = None
        ValueError('unsupported element type ' + str(element_type))

    heavy_spectrum = pd.Series(spectrum.values * ratio, index=spectrum.index + mass_diff / charge)
    combined_spectrum = pd.concat([spectrum, heavy_spectrum]).sort_index()

    """
    doesn't appear to have any function. current_spectrum is discarded
    
    # iteratively sum the intensity of similar mz within a tolerance
    last_spectrum = None
    current_spectrum = sum_intensity_of_similar_mz(combined_spectrum)
    while not current_spectrum.equals(last_spectrum):
        last_spectrum = current_spectrum
        current_spectrum = sum_intensity_of_similar_mz(current_spectrum)
    """
    return sum_intensity_of_similar_mz(combined_spectrum)


def sum_intensity_of_similar_mz(spectrum):
    """if the mz are within tolerance, sum the intensities"""
    idx = 0
    mz_array = np.array(spectrum.index)
    int_array = np.array(spectrum.values)
    mz_array_blurred = []
    int_array_blurred = []
    while (idx + 1) < mz_array.shape[0]:
        if within_tol_prob(mz_array[idx], mz_array[idx + 1], param.MS1_TOLERANCE, 'ppm'):
            # perform merging
            mz_array_blurred.append((mz_array[idx] + mz_array[idx + 1]) / 2.)
            int_array_blurred.append(int_array[idx] + int_array[idx + 1])
            idx += 2
        else:
            # m/z array is guranteed to be strictly increasing
            # so the only case left is mz1+2*tolerance < mz2
            # push idx and advance both
            mz_array_blurred.append(mz_array[idx])
            int_array_blurred.append(int_array[idx])
            idx += 1
    # handle residual peaks by extending with original unmodified values
    if idx < mz_array.shape[0]:
        mz_array_blurred.extend(mz_array[idx::])
        int_array_blurred.extend(int_array[idx::])
    return pd.Series(np.array(int_array_blurred), index=np.array(mz_array_blurred)).sort_index()


def print_messagebox(msg):
    print(u'\u2554' + u'\u2550' * (len(msg) + 2) + u'\u2557')
    print(u'\u2551' + ' ' + msg + ' ' + u'\u2551')
    print(u'\u255a' + u'\u2550' * (len(msg) + 2) + u'\u255d')


def reducer(accumulator, element):
    for key, value in element.items():
        accumulator[key] = accumulator.get(key, 0) + value
    return accumulator


def contains_precursor(ms1_row, precursor_mz, precursor_window, intensity_floor):
    ms1_mz_int = pd.Series(ms1_row['intensity'], ms1_row['mz'])
    return ms1_mz_int.loc[
           precursor_mz - precursor_window: precursor_mz + precursor_window].sum() > intensity_floor


def dynamic_rt_filter(ms1_subset_df, ms2_row, intensity_floor):
    """
    plots the ms1 intensity over time (rt), limited to the ms2 precursor_mz.
    From the ms2 precursor time (rt), search for an upper and lower bound time (rt) that is still above the
    intensity floor. The goal is to identify the time window outputting quality data

    :param ms1_subset_df: dataframe
    :param ms2_row: pd.Series. one row of the ms2
    :return: pd.DataFrame of ms1 with a common mz as the index and rt as column
    """
    intensity_ts = (ms1_subset_df.loc[ms2_row['obs_precursor_mz'] - param.MS1_PRECURSOR_WINDOW:
                                      ms2_row['obs_precursor_mz'] + param.MS1_PRECURSOR_WINDOW].sum(0))
    rolling_intensity_ts = intensity_ts.rolling(window=3, min_periods=1, center=True).mean()
    try:
        # lower bound
        lower_half = (rolling_intensity_ts.loc[:ms2_row['rt']][::-1])
        lower_half_cumsum = (lower_half < intensity_floor).cumsum()
        lower_half_above_floor = lower_half_cumsum[lower_half_cumsum == 0]
        if lower_half_above_floor.shape[0] > 0:
            rt_floor = lower_half_above_floor.index[-1]
        else:
            rt_floor = ms2_row['rt']

        # upper bound
        upper_half = (rolling_intensity_ts.loc[ms2_row['rt']:])
        upper_half_cumsum = (upper_half < intensity_floor).cumsum()
        upper_half_above_floor = upper_half_cumsum[upper_half_cumsum == 0]
        if upper_half_above_floor.shape[0] > 0:
            rt_ceil = upper_half_above_floor.index[-1]
        else:
            rt_ceil = ms2_row['rt']

        return rt_floor, rt_ceil, rolling_intensity_ts
    except IndexError:
        logging.info(f"""no ms1 match for ms2: {ms2_row}. Return original without applying dynamic_rt_filter""")
        return ms1_subset_df.columns.min(), ms1_subset_df.columns.max(), rolling_intensity_ts


def unit_fuzzy_match(ms2):
    """
    iterates through ms2 rows and returns rows with similar mz, rt, and equal charge
    :param ms2: dataframe of ms2 observations
    :return: a dictionary of sets containing indices similar to each ms2
    """
    output = {}
    for idx, row in ms2.iterrows():
        mask_mz = abs(1 - row['obs_precursor_mz'] / ms2['obs_precursor_mz']) <= param.MS1_TOLERANCE
        mask_charge = row['obs_precursor_charge'] == ms2['obs_precursor_charge']
        mask_rt = abs(ms2['ms2_rt'] - row['ms2_rt']) < (2 * param.MS2_RT_WINDOW)
        output[idx] = set(ms2[mask_mz & mask_charge & mask_rt].drop(idx).index.unique())
    return output


def unit_fuzzy_match_v2(ms2):
    """
    optimization attemp to group the calculation by charge state
    :param ms2: dataframe of ms2 observations
    :return: a dictionary of sets containing indices similar to each ms2
    """
    output = {}
    for charge in ms2['obs_precursor_charge'].unique():
        tmp = ms2[ms2['obs_precursor_charge'] == charge]
        for idx, row in tmp.iterrows():
            mask_mz = abs(1 - row['obs_precursor_mz'] / tmp['obs_precursor_mz']) <= param.MS1_TOLERANCE
            mask_rt = abs(tmp['ms2_rt'] - row['ms2_rt']) < (2 * param.MS2_RT_WINDOW)
            output[idx] = set(tmp[mask_mz & mask_rt].drop(idx).index.unique())
    return output


def unit_fuzzy_match_v3(ms2):
    """This doesn't work due to memory issues. I'm trying to do a self join to avoid iterating through the dataframe"""
    for charge in ms2['charge'].unique():
        mdf = (ms2.loc[ms2['obs_precursor_charge'] == charge, ['obs_precursor_mz', 'ms2_rt']].assign(foo=1)
               .merge(ms2.loc[ms2['obs_precursor_charge'] == charge, ['obs_precursor_mz', 'ms2_rt']].assign(foo=1),
                      on='foo', suffixes=('', '_y'))
               .drop('foo', 1))


def ts_reindex(ts, precision):
    fill_idx = np.arange(ts.index.min(),
                         ts.index.max(), step=10 ** -precision)
    # must round to remove floating point accuracy issues. i.e. 32 vs 31.999999999
    return ts.reindex(fill_idx.round(precision), fill_value=0)


def ts_reindex_v2(ts, precision):
    """
    rounds to a precision, groups by that precision and fills in gaps
    reset_index assumes a single index. does not accept multi-index"""
    tmp = ts.reset_index()
    tmp.iloc[:, 0] = tmp.iloc[:, 0].round(precision)
    tmp = tmp.fillna(0).groupby(tmp.columns[0]).sum()
    fill_idx = np.arange(tmp.index.min(),
                         tmp.index.max() + 10 ** -precision, step=10 ** -precision)
    tmp = tmp.reindex(fill_idx.round(precision), fill_value=0)

    # reset index changes the type to pd.DataFrame. return Series if input is series
    if type(ts) == pd.core.series.Series:
        tmp = tmp.iloc[:, 0]
    return tmp


def depth_first_search(graph, start):
    """from the start index, explore and return the connected nodes (a.k.a. visited)"""
    # https://eddmann.com/posts/depth-first-search-and-breadth-first-search-in-python/
    visited = set()
    stack = [start]
    while len(stack) > 0:
        vertex = stack.pop()
        if vertex not in visited:
            visited.add(vertex)
            stack.extend(graph[vertex] - visited)
    return visited


def reduce_by_dfs(adj_dict):
    """iterates through the keys and appends the connected nodes.
       Excludes the connected nodes from the keys to iterate"""
    global_visited = set()
    group_list = []
    for key in adj_dict.keys():
        if key not in global_visited:
            local_visited = depth_first_search(adj_dict, key)
            group_list.append(local_visited)
            global_visited.update(local_visited)
    return group_list
