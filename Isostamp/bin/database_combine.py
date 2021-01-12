import logging
import numpy as np
import os
import pandas as pd
import parameters as param
import time
from pathlib import Path
from database_combine_functions import database_read_transform, propensity_from_first_scan_v2
from util_functions import camelcase_to_snake_case
logging.basicConfig(level=logging.INFO)


def run(directory, input_folder, output_folder, target_csv_filename, decoy_csv_filename, output_filename):
    # consider adding file name handling for anticipated file names. i.e. uppercase vs lowercase

    input_path = os.path.join(directory, input_folder)
    output_path = os.path.join(directory, output_folder)

    # These two files are optional. Used to merge the ms1df with these Sequest files
    if target_csv_filename:
        target_csv_path = os.path.join(input_path, target_csv_filename)
    else:
        target_csv_path = None

    if decoy_csv_filename:
        decoy_csv_path = os.path.join(input_path, decoy_csv_filename)
    else:
        decoy_csv_path = None

    # information about the input files
    time_start = time.time()
    psmdf = database_read_transform(target_csv_path, decoy_csv_path)

    mdf = []
    for spectrum_file in psmdf['spectrum_file'].unique():
        logging.info(f'running {spectrum_file}')
        psmdf_filter = psmdf[psmdf['spectrum_file'] == spectrum_file].copy()
        file_stem = Path(spectrum_file).resolve().stem
        ms1df_file = file_stem + '_ms1df.csv'
        ms1df_path = os.path.join(input_path, ms1df_file)

        try:
            ms1df = pd.read_csv(ms1df_path, index_col=0)
            # input files are generally camelcase. rename to snake case consistent with python pep8
            ms1df.columns = [camelcase_to_snake_case(col) for col in ms1df.columns]
            ms1df['ms2_scan_number'] = ms1df['ms2_scan_number'].apply(lambda x: [int(s) for s in x.split('+')])

            mdf.append(propensity_from_first_scan_v2(psmdf_filter, ms1df))

        except:
            logging.warning(f'WARNING! {spectrum_file} does not have a corresponding ms1df ({ms1df_path})')
            mdf.append(psmdf_filter.assign(isostamp_propensity=np.nan))

    mdf = pd.concat(mdf, axis=0, sort=False).sort_values('spectrum_file')

    logging.info('step 4. apply 40% normalized isostamp score and %60 normalized xcorrelation weighting')
    # mdf['normalized_isostamp'] = scale(mdf['isostamp_propensity'])
    # mdf['normalized_regular'] = scale(mdf['regular_propensity'])
    # mdf['normalized_xcorr'] = scale(mdf['xcorr'])
    # mdf['combo_score'] = (0.4 * mdf['normalized_isostamp'] + 0.6 * mdf['normalized_xcorr'])

    logging.info('step 5. re-rank by Target-Decoy Approach')
    # mdf['combo_fdr'] = calculate_q_value(mdf['tda'].values, mdf['combo_score'])

    mdf = mdf.round(param.PRECISION)
    mdf = mdf[['activation_type', 'annotated_sequence', 'apex_rt_min', 'charge',
               'checked', 'confidence', 'delta_cn', 'delta_m_ppm', 'delta_mz_da',
               'delta_score', 'file_id', 'first_scan', 'identifying_node', 'index',
               'ion_inject_time_ms',  'isolation_interference_pct', 'master_protein_accessions',
               'master_protein_descriptions', 'mhplus_da', 'modifications', 'ms_order',
               'mz_da', 'number_missed_cleavages', 'number_protein_groups',
               'number_proteins', 'precursor_abundance', 'protein_accessions',
               'psm_ambiguity', 'quan_info', 'rank', 'rt_min',
               'search_engine_rank', 'spectrum_file', 'tda', 'theo_mhplus_da', 'xcorr',
               'regular_propensity', 'isostamp_propensity', 'base_peak_correlation']]
    mdf.to_csv(os.path.join(output_path, output_filename), index=False)

    logging.info(f'target Sequest search result CSV file path: {target_csv_path}')
    logging.info(f'decoy Sequest search result CSV file path: {decoy_csv_path}')
    logging.info(f'Completed in {((time.time() - time_start) / 60.):.2f} minutes')


if __name__ == '__main__':
    directory = '/hd/data/david_manuscript/'
    input_folder = 'isostamp_test_output'
    output_folder = 'isostamp_test_output'
    target_csv_filename = 'DKM151_cleavage_PSMs.txt'
    decoy_csv_filename = ''
    output_filename = target_csv_filename.split('.')[0] + '_isostamp.csv'

    run(directory, input_folder, output_folder, target_csv_filename, decoy_csv_filename, output_filename)
