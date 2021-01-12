import logging
import numpy as np
import parameters as param
import pandas as pd
logging.basicConfig(level=logging.INFO)


def database_read_transform(target_csv_path, decoy_csv_path):
    """
    This is function is highly dependent on the input file formats. For example, xcorr column is labeled different
    names based on the search procedure in a way that I don't understand yet. So, do not be surprised if this function
    raises an exception. Debug by looking at the dataframe and adding the appropriate column mapping.
    :param target_csv_path:
    :param decoy_csv_path:
    :return:
    """

    logging.info('read Sequest target and/or decoy files')
    # Isostamp output is used as an input in other projects
    # I chose to camelcase_to_snake_case columns into a pythonic naming convention
    column_map = {'Checked': 'checked',
                  'Confidence': 'confidence',
                  'Identifying Node': 'identifying_node',
                  'PSM Ambiguity': 'psm_ambiguity',
                  'Annotated Sequence': 'annotated_sequence',
                  'Modifications': 'modifications',
                  '# Proteins': 'number_proteins',
                  'Master Protein Accessions': 'master_protein_accessions',
                  'Protein Accessions': 'protein_accessions',
                  '# Missed Cleavages': 'number_missed_cleavages',
                  'Charge': 'charge',
                  'DeltaScore': 'delta_score',
                  'DeltaCn': 'delta_cn',
                  'Rank': 'rank',
                  'Search Engine Rank': 'search_engine_rank',
                  'm/z [Da]': 'mz_da',
                  'MH+ [Da]': 'mhplus_da',
                  'Theo. MH+ [Da]': 'theo_mhplus_da',
                  'DeltaM [ppm]': 'delta_m_ppm',
                  'Deltam/z [Da]': 'delta_mz_da',
                  'Activation Type': 'activation_type',
                  'MS Order': 'ms_order',
                  'Isolation Interference [%]': 'isolation_interference_pct',
                  'Ion Inject Time [ms]': 'ion_inject_time_ms',
                  'RT [min]': 'rt_min',
                  'First Scan': 'first_scan',
                  'Spectrum File': 'spectrum_file',
                  'File ID': 'file_id',
                  '# Protein Groups': 'number_protein_groups',
                  'Master Protein Descriptions': 'master_protein_descriptions',
                  'Quan Info': 'quan_info',
                  'Precursor Abundance': 'precursor_abundance',
                  'Apex RT [min]': 'apex_rt_min',
                   'XCorr': 'xcorr'}

    # either or both files can exist
    if target_csv_path and decoy_csv_path:
        target_df = (pd.read_csv(target_csv_path, delimiter='\t')
                     .rename(columns=column_map)
                     .assign(tda='target'))
        decoy_df = (pd.read_csv(decoy_csv_path, delimiter='\t')
                    .rename(columns=column_map)
                    .assign(tda='decoy'))
        psmdf = pd.concat([target_df, decoy_df], axis=0, sort=False)
    elif target_csv_path:
        psmdf = (pd.read_csv(target_csv_path, delimiter='\t')
                 .rename(columns=column_map)
                 .assign(tda='target'))
    elif decoy_csv_path:
        psmdf = (pd.read_csv(decoy_csv_path, delimiter='\t')
                 .rename(columns=column_map)
                 .assign(tda='decoy'))
    else:
        return

    psmdf = psmdf.loc[psmdf['search_engine_rank'] == 1]
    if 'xcorr' not in psmdf.columns:
        psmdf.rename(columns={'|Log Prob|': 'xcorr'})
    psmdf.sort_values('xcorr', ascending=False, inplace=True)

    psmdf.reset_index(inplace=True)

    return psmdf


def propensity_from_first_scan_v2(psmdf, ms1df):
    """
    This is an experimental function to reformat the ms1df to enable a dataframe merge to remove a for loop
    :param ms1df:
    :return:
    """

    # ms2_scan_number is a column of lists. expand each ms2_scan_number to a separate index
    # source: https://stackoverflow.com/questions/27263805/pandas-column-of-lists-create-a-row-for-each-list-element
    lst_col = 'ms2_scan_number'
    ms1 = (pd.DataFrame({
        col: np.repeat(ms1df[col].values, ms1df[lst_col].str.len())
        for col in ms1df[['isostamp_propensity', 'regular_propensity', 'base_peak_correlation']]}
    ).assign(**{lst_col: np.concatenate(ms1df[lst_col].values)}))

    # take the first value to ensure uniqueness
    ms1 = ms1.groupby('ms2_scan_number', as_index=False).first()
    ms1 = ms1.rename(columns={'ms2_scan_number': 'first_scan'})

    df = psmdf.merge(ms1, on='first_scan', how='left')

    return df
