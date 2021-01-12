import csv
import json
import os
import sqlite3
import tqdm
from parameters import MASS_H
from pyteomics import mass
from scipy import special, interpolate
from sklearn.preprocessing import scale
from util_functions import *


"""
These are functions that were inherited from Hung Yi, but are not used in the current isostamp package.
Segregating these in preparation for deletion in a future version
"""

# average mass of a single amino acid
# reference http://bionumbers.hms.harvard.edu/bionumber.aspx?id=107786
mass_aa = 118.9

# amino acid residue mass dictionary
mass_NH3 = 17.031
mass_H2O = 18.01056
# parameter for Sequest spectrum pre-processing
sequest_mz_start = 0
sequest_mz_end = 2000
sequest_index_offset = 0.5
sequest_default_bin_size = 1.0005079
sequest_max_delay = mass_aa / 2.

# common modification table {name:monoisotopic mass change}
common_modification_table = {'Deamidated': 0.98402,
                             'Carbamidomethyl': 57.02146,
                             'Oxidation': 15.99491}


def expand_db_func(input_tuple):
    # unpack input
    db_path = input_tuple[0]
    mod_path = input_tuple[1]
    bucket_start_point = input_tuple[2]
    bucket_size = input_tuple[3]
    expanded_db_path = input_tuple[4]
    # clean up
    if os.path.isfile(expanded_db_path):
        os.remove(expanded_db_path)
    # load modification rule
    with open(mod_path, 'r') as filePointer:
        mod_rule = json.load(filePointer)
    # open database reader
    db_file_pointer = open(db_path, 'r')
    db_reader = csv.reader(db_file_pointer, delimiter=',')
    # open SQLite database
    conn = sqlite3.connect(expanded_db_path)
    cur = conn.cursor()
    cur.execute("""CREATE TABLE expandedDB
                (uniprotAccession TEXT,
                peptideSequence TEXT,
                mass_lower_bound REAL,
                mass_upper_bound REAL);""")
    # go through database
    count_insert = 0
    for rowIndex, row in enumerate(db_reader):
        if rowIndex < bucket_start_point:
            continue
        elif rowIndex > bucket_start_point + bucket_size - 1:
            break
        else:
            uniprot_accession = row[0]
            peptide_sequence = row[1]
            # calculate possible mass range
            try:
                min_possible_mass, max_possible_mass = calculate_massbound(candidate_sequence=peptide_sequence,
                                                                           mod_rule=mod_rule)
                cur.execute("""INSERT INTO expandedDB 
                        (uniprotAccession, peptideSequence, mass_lower_bound, mass_upper_bound)
                        VALUES (?, ?, ?, ?);""", [uniprot_accession, peptide_sequence, min_possible_mass, max_possible_mass])
            except:
                continue
    # create index for fast slicing
    cur.execute("""CREATE INDEX index_mass_lower_bound ON expandedDB (mass_lower_bound)""")
    cur.execute("""CREATE INDEX index_mass_upper_bound ON expandedDB (mass_upper_bound)""")
    conn.commit()
    conn.close()
    db_file_pointer.close()
    return


def calculate_massbound(candidate_sequence, mod_rule):
    """
    if observed precursor mass is larger than candidate sequence with all positive modifications
    of smaller than candidate sequence with all negative modifications
    return False to skip the search
    """
    # static
    base_mass = mass.fast_mass(sequence=candidate_sequence)
    for mod in mod_rule['negativeStatic']:
        monoisotopic_mass_change = mod[1]
        aa = mod[2]
        count = candidate_sequence.count(aa)
        if count > 0:
            base_mass += count * monoisotopic_mass_change
    for mod in mod_rule['positiveStatic']:
        monoisotopic_mass_change = mod[1]
        aa = mod[2]
        count = candidate_sequence.count(aa)
        if count > 0:
            base_mass += count * monoisotopic_mass_change
    # dynamic, with Sequest default setting
    get_mimc = lambda item: item[1]
    min_possible_mass = base_mass
    key_list = sorted(mod_rule['negativeDynamic'], key=get_mimc)
    mod_quota = 4
    pal_quota = 1
    for mod in key_list:
        mod_name = mod[0]
        monoisotopic_mass_change = mod[1]
        aa = mod[2]
        count = min(candidate_sequence.count(aa), 3)
        if '[PAL]' in mod_name and 0 < pal_quota < mod_quota:
            min_possible_mass += pal_quota * monoisotopic_mass_change
            mod_quota -= 1
            pal_quota -= 1
        elif '[PAL]' in mod_name and pal_quota > 0 and mod_quota <= pal_quota:
            min_possible_mass += mod_quota * monoisotopic_mass_change
            break
        elif '[PAL]' not in mod_name and mod_quota <= count:
            min_possible_mass += mod_quota * monoisotopic_mass_change
            break
        elif '[PAL]' not in mod_name and mod_quota > count:
            min_possible_mass += count * monoisotopic_mass_change
            mod_quota -= count
    max_possible_mass = base_mass
    key_list = sorted(mod_rule['positiveDynamic'], key=get_mimc)[::-1]
    mod_quota = 4
    pal_quota = 1
    for mod in key_list:
        mod_name = mod[0]
        monoisotopic_mass_change = mod[1]
        aa = mod[2]
        count = min(candidate_sequence.count(aa), 3)
        if '[PAL]' in mod_name and 0 < pal_quota < mod_quota:
            max_possible_mass += pal_quota * monoisotopic_mass_change
            mod_quota -= 1
            pal_quota -= 1
        elif '[PAL]' in mod_name and pal_quota > 0 and mod_quota <= pal_quota:
            max_possible_mass += mod_quota * monoisotopic_mass_change
            break
        elif '[PAL]' not in mod_name and mod_quota <= count:
            max_possible_mass += mod_quota * monoisotopic_mass_change
            break
        elif '[PAL]' not in mod_name and mod_quota > count:
            max_possible_mass += count * monoisotopic_mass_change
            mod_quota -= count
    return min_possible_mass, max_possible_mass


def expectation(input_tuple):
    # unpack input
    srg_path = input_tuple[0]
    occurrence_dict = input_tuple[1]
    prior_dict = input_tuple[2]
    mode = input_tuple[3]
    full_mod_name_list = input_tuple[4]
    # load search result of an MS2
    srg = pd.read_csv(srg_path, index_col=0)

    # maximization of likelihood
    # use prior to find the top PSM
    # use numpy array to speed up
    # columns: (count_aa, count_mod) for each mod in prior_dict.keys(), then log_posterior and posterior
    index_array = np.array(srg.index)
    peptide_sequence_array = srg['peptideSequence'].values
    mod_string_array = srg['mod_string'].values
    mod_count_string_array = srg['mod_count_string'].values
    grouped_likelihood_array = srg['groupedLikelihood'].values
    mass_prior_array = srg['massPrior'].values
    pal_prior_array = srg['PALPrior'].values
    count_mod = len(prior_dict)
    operation_array = np.zeros((index_array.shape[0], 2 * count_mod + 2))
    for i in range(srg.shape[0]):
        mod_count_string = mod_count_string_array[i]
        log_posterior = np.log(grouped_likelihood_array[i]) + np.log(mass_prior_array[i]) + np.log(pal_prior_array[i])
        for modIndex, mod_name in enumerate(prior_dict.keys()):
            # calculate maximum and actual occurrence
            aa_list = occurrence_dict[mod_name]
            count_aa = sum([peptide_sequence_array[i].count(AA) for AA in aa_list])
            if count_aa == 0:
                continue
            if mod_count_string is None or pd.isnull(mod_count_string) or len(mod_count_string) == 0:
                count_mod = 0
            else:
                relevant_mod_count_string = [m for m in mod_count_string.split(',') if m.startswith(mod_name + '(')][0]
                count_mod = int(relevant_mod_count_string.split('(')[1][0:-1])
            count_unmod = count_aa - count_mod
            # model it by multiple non-IID Bernoulli experiment
            # in order not to penalize dominant case, set a nonlinearity here
            if prior_dict[mod_name] == 1:
                # to handle numerical instability
                p = 1
                q = 1e-3
            elif prior_dict[mod_name] == 0:
                # to handle numerical instability
                p = 1e-3
                q = 1
            elif prior_dict[mod_name] == 0.5:
                p = 1
                q = 1
            elif prior_dict[mod_name] > 0.5:
                p = 1
                q = (1 - prior_dict[mod_name]) * 2
            else:
                p = prior_dict[mod_name] * 2
                q = 1
            prior = (p ** count_mod) * (q ** (count_aa - count_mod))
            log_posterior += np.log(prior)
            operation_array[i, 2 * modIndex] = count_aa
            operation_array[i, 2 * modIndex + 1] = count_mod
        operation_array[i, -2] = log_posterior
    operation_array[:, -2] -= special.logsumexp(operation_array[:, -2])
    operation_array[:, -1] = np.exp(operation_array[:, -2])
    # find most-likely PSM
    order_key = np.argsort(operation_array[:, -1])[::-1]
    index_array = index_array[order_key]
    operation_array = operation_array[order_key, :]
    if mode == 'optimize':
        # report prevalence
        prevalence = {modName: 0 for modName in occurrence_dict.keys()}
        for modIndex, mod_name in enumerate(prevalence.keys()):
            # hard clustering
            prevalence[mod_name] = [operation_array[0, 2 * modIndex], operation_array[0, 2 * modIndex + 1]]
        # report top posterior
        return prevalence, operation_array[0, -1]
    elif mode == 'evaluate':
        info_list = srg.loc[index_array[0],
                           ['ms2_scan_number', 'uniprotAccession', 'peptideSequence',
                            'mod_string', 'obs_precursor_mass', 'groupedLikelihood', 'massPrior',
                            'PALPrior']].tolist()
        srg['posterior'] = operation_array[:, -1]
        path_list = srg_path.split('/')
        path_list[-1] = 'validated' + path_list[-1][7::]
        srv_path = os.path.join(*path_list)
        srg.to_csv(srv_path)
        return info_list + [operation_array[0, -1]]


def aggregate_psm(input_tuple):
    # unpack input
    sr_path = input_tuple[0]
    mod_name_list = input_tuple[1]
    # load dataframe
    sr = pd.read_csv(sr_path, index_col=0)
    # get shared properties
    ms2_scan_number = sr.loc[0, 'ms2_scan_number']
    obs_precursor_mass = sr.loc[0, 'obs_precursor_mass']
    # prepare array for fast random access
    peptide_sequence_array = sr['peptideSequence'].values
    uniprot_accession_array = sr['uniprotAccession'].values
    mod_string_array = sr['mod_string'].values
    # camelcase_to_snake_case XCorr to probability through softmax and then operate in the probability space
    log_prob_array = sr['XCorr'].values
    log_prob_array -= special.logsumexp(log_prob_array)
    prob_array = np.exp(log_prob_array)
    # other stuff
    mass_prior_array = sr['massPrior'].values
    pal_prior_array = sr['PALPrior'].values
    record_dict = {}
    for index in range(sr.shape[0]):
        mod_pos_dict = {modName: [] for modName in mod_name_list}
        mod_string = mod_string_array[index]
        if pd.notnull(mod_string) and len(mod_string) > 0:
            for modName in mod_name_list:
                mod_pos_dict[modName] = [m.split('(')[0] for m in mod_string.split(',') if modName in m]
        mod_count_list = []
        for modName in mod_name_list:
            mod_count_list.append(modName + '({})'.format(len(mod_pos_dict[modName])))
        mod_count_string = ','.join(mod_count_list)
        peptide_sequence = peptide_sequence_array[index]
        id = (peptide_sequence, mod_count_string)
        if id not in record_dict.keys():
            record_dict[id] = [uniprot_accession_array[index],
                               peptide_sequence,
                               {k: set(mod_pos_dict[k]) for k in mod_pos_dict},
                               [prob_array[index]],
                               mass_prior_array[index],
                               pal_prior_array[index],
                               mod_count_string]
        else:
            for modName in mod_name_list:
                record_dict[id][2][modName].update(mod_pos_dict[modName])
            record_dict[id][3].append(prob_array[index])
    record = []
    for key in record_dict.keys():
        uniprot_accession = record_dict[key][0]
        peptide_sequence = record_dict[key][1]
        prob_list = record_dict[key][3]
        mass_prior = record_dict[key][4]
        pal_prior = record_dict[key][5]
        mod_count_string = record_dict[key][6]
        # prepare aggregated modification string
        mod_list = []
        for modName in record_dict[key][2].keys():
            pos_set = record_dict[key][2][modName]
            if len(pos_set) > 0:
                pos_string = '{' + '|'.join(sorted(pos_set, key=lambda i: int(i[1::]))) + '}'
                mod_list.append(pos_string + '(' + modName + ')')
        mod_string = ','.join(mod_list)
        # prepare union probability
        union_prob = 1 - np.prod(1 - np.array(prob_list))
        # record
        record.append([uniprot_accession,
                       peptide_sequence,
                       mod_string,
                       union_prob,
                       mass_prior,
                       pal_prior,
                       mod_count_string])
    srg = pd.DataFrame.from_records(record,
                                    columns=['uniprotAccession',
                                             'peptideSequence',
                                             'mod_string',
                                             'groupedLikelihood',
                                             'massPrior',
                                             'PALPrior',
                                             'modCountString'])
    # put shared properties back
    srg['ms2_scan_number'] = ms2_scan_number
    srg['obs_precursor_mass'] = obs_precursor_mass
    # export
    path_list = sr_path.split('/')
    path_list[-1] = 'grouped' + path_list[-1]
    srg_path = os.path.join(*path_list)
    srg.to_csv(srg_path)
    return srg_path


def interp_poisson_cdf(x, rate):
    """
    use scipy.interpolate to approximate Poisson CDF
    """
    if x < 0:
        return 0
    else:
        upper_bound = max(np.ceil(x) + 1, 4)
        discrete_grid = np.arange(0, upper_bound)
        discrete_cdf = stats.poisson.cdf(discrete_grid, rate)
        ecdf = interpolate.interp1d(discrete_grid, discrete_cdf, kind='cubic')
        return ecdf(x)


def modification_generator(mod_rule, seq, mass_lower_bound, mass_upper_bound):
    """
    a generator/iterable for all possible dynamic modifications
    output is a modification dictionary ready for generating spectrum
    """
    # check necessity of processing
    rule_dict = dict()
    dynamic_mod_name_list = []
    dynamic_mod_mass_list = []
    for mod in mod_rule['positiveDynamic'] + mod_rule['negativeDynamic']:
        mod_name = mod[0]
        mimc = mod[1]
        aa = mod[2]
        dynamic_mod_name_list.append(mod_name)
        dynamic_mod_mass_list.append(mimc)
        if aa in rule_dict.keys():
            rule_dict[aa][1].append((mod_name, mimc))
        else:
            rule_dict[aa] = ['dynamic', [(mod_name, mimc)]]
    for mod in mod_rule['positiveStatic'] + mod_rule['negativeStatic']:
        mod_name = mod[0]
        mimc = mod[1]
        aa = mod[2]
        if aa in rule_dict.keys():
            rule_dict[aa][1].append((mod_name, mimc))
        else:
            rule_dict[aa] = ['static', [(mod_name, mimc)]]

    if len(set(seq).intersection(set(rule_dict.keys()))) == 0:
        yield {}, ''
    else:
        # find occurrences in the sequence
        choice_list = []
        pos_list = []
        for pos, aa in enumerate(seq):
            if aa in rule_dict.keys():
                mod_type = rule_dict[aa][0]
                choice = rule_dict[aa][1]
                if mod_type == 'dynamic':
                    choice_list.append(choice + [(None, 0)])
                    pos_list.append(pos)
                elif mod_type == 'static':
                    choice_list.append(choice)
                    pos_list.append(pos)
        counter_dict = {k: [w, 0] for k, w in zip(dynamic_mod_name_list, dynamic_mod_mass_list)}
        for combo in recursive_mod_func(choice_list, counter_dict, mass_lower_bound, mass_upper_bound):
            mod_dict = dict()
            mod_string_list = []
            for pos, mod in zip(pos_list, combo):
                if mod[0] is not None:
                    aa = seq[pos]
                    mod_name = mod[0]
                    mod_mass = mod[1]
                    mod_dict[pos] = mod_mass
                    mod_string_list.append(aa + str(pos + 1) + '(' + mod_name + ')')
            mod_string = ','.join(mod_string_list)
            yield mod_dict, mod_string


def recursive_mod_func(choice_list, counter_dict, mass_lower_bound, mass_upper_bound):
    choice = choice_list[0]
    for decision in choice:
        mod_name = decision[0]
        # default Sequest setting:
        # up to 3 same dynamic modifications per peptide
        # up to 4 total dynamic modifications per peptide
        if mod_name is None:
            violate_single_mod_constraint = False
            violate_total_mod_constraint = False
            current_mass = sum([counter_dict[k][0] * counter_dict[k][1] for k in counter_dict.keys()])
            next_mass = 0
            violate_mass_constraint = current_mass + next_mass > mass_upper_bound
        else:
            violate_single_mod_constraint = counter_dict[mod_name][1] >= 3
            violate_total_mod_constraint = sum([v[1] for v in counter_dict.values()]) >= 4
            current_mass = sum([counter_dict[k][0] * counter_dict[k][1] for k in counter_dict.keys()])
            next_mass = counter_dict[mod_name][0]
            violate_mass_constraint = current_mass + next_mass > mass_upper_bound
        if (mod_name is not None) and (mod_name in counter_dict.keys()) \
                and (violate_single_mod_constraint or violate_total_mod_constraint or violate_mass_constraint):
            continue
        # green light for this modification decision
        current_counter_dict = {k: [counter_dict[k][0], counter_dict[k][1]] for k in counter_dict.keys()}
        if mod_name is not None and mod_name in current_counter_dict.keys():
            current_counter_dict[mod_name][1] += 1
        if len(choice_list) > 1:
            for nextDecision in recursive_mod_func(choice_list[1::], current_counter_dict,
                                                   mass_lower_bound, mass_upper_bound):
                yield [decision] + nextDecision
        elif current_mass + next_mass >= mass_lower_bound:
            yield [decision]


def sequest_preprocessing(int_array, mz_start=sequest_mz_start, mz_end=sequest_mz_end, step=sequest_default_bin_size):
    # step 1. sqrt transform of intensity
    int_array = np.sqrt(int_array)
    # step 2. region-normalize intensity
    zone_count = int((mz_end - mz_start) / mass_aa)
    loc = int((mz_start - sequest_mz_start) / step + sequest_index_offset * step)
    scale = int((mz_end - mz_start) / step)
    zone_size, extra = divmod(int_array.shape[0], zone_count)
    zone_size_list = extra * [zone_size + 1] + (zone_count - extra) * [zone_size]
    zone_start_list = np.cumsum([0] + zone_size_list)[0:-1]
    for zone_start, zone_size in zip(zone_start_list, zone_size_list):
        base_int = int_array[zone_start:zone_start + zone_size].max()
        if base_int > 0:
            int_array[zone_start:zone_start + zone_size] /= base_int
    return int_array


def seq_to_spectrum_sequest(seq, precursor_charge=2, isotope_array=np.arange(1), ion_type=('b', 'y'),
                            neutral_loss=('H2O', 'NH3'), mod_string='', step=sequest_default_bin_size):
    """
    camelcase_to_snake_case sequence to spectrum
    modification dictionary format: {0:+57, 3:+1}
    supports common modification table (common_modification_table) annotation
    use modification string to detect/handle pattern-generating modifications
    """
    # reconstruct mod_dict
    mod_dict = {}
    if len(mod_string) > 0:
        for s in mod_string.split(','):
            pos = int(s.split('(')[0][1::])
            mod_name = s.split('(')[1][0:-1]
            if mod_name in common_modification_table.keys():
                mod_dict[pos - 1] = common_modification_table[mod_name]
    # detect all mods
    mod_pos_array = np.zeros(len(mod_dict))
    for index, modPos in enumerate(mod_dict.keys()):
        mod_pos_array[index] = modPos
    # detect PAL mods
    pal_mod_pos_list = []
    for mod in mod_string.split(','):
        if len(mod) > 0:
            pos = int(mod.split('(')[0][1::])
            mod_name = mod.split('(')[1][0:-1]
            if '[PAL]' in mod_name:
                pal_mod_pos_list.append(pos)
    pal_mod_pos_array = np.array(pal_mod_pos_list)
    # detect neutral loss
    neutral_loss_list = [0]
    for n in neutral_loss:
        if n == 'H2O':
            neutral_loss_list.append(-mass_H2O)
        elif n == 'NH3':
            neutral_loss_list.append(-mass_NH3)
    # empty container
    mz_list = []
    int_list = []
    max_charge = min(2, int(precursor_charge - 1))
    tier1_int = 1
    tier2_int = 0.5
    tier3_int = 0.2
    padding_int = 0.25
    for pos in range(1, len(seq)):
        for charge in np.arange(1, max_charge + 1):
            for ion_type in ion_type:
                if ion_type in 'abc':
                    ion_mz = mass.fast_mass(seq[0:pos], ion_type=ion_type, charge=charge)
                    for modPos in mod_pos_array[mod_pos_array < pos]:
                        ion_mz += mod_dict[modPos] / charge
                    # handle neutral loss
                    for loss in neutral_loss_list:
                        ion_mz_n = ion_mz + loss / charge
                        # handle isotope and PAL
                        pal_count = (pal_mod_pos_array < pos).sum()
                        for p in np.arange(0, pal_count + 1):
                            shift = p * 2 * (-MASS_H) / charge
                            ion_mz_np = ion_mz_n + shift
                            ratio = (1. / 3.) ** p
                            for i in isotope_array:
                                mz_list.append(ion_mz_np + i * MASS_H / charge)
                                if loss == 0 and charge == 1:
                                    int_list.append(tier1_int * ratio)
                                elif loss == 0 and charge > 1:
                                    int_list.append(tier2_int * ratio)
                                else:
                                    int_list.append(tier3_int * ratio)
                elif ion_type in 'xyz':
                    ion_mz = mass.fast_mass(seq[pos::], ion_type=ion_type, charge=charge)
                    for modPos in mod_pos_array[mod_pos_array >= pos]:
                        ion_mz += mod_dict[modPos] / charge
                    # handle neutral loss
                    for loss in neutral_loss_list:
                        ion_mz_n = ion_mz + loss / charge
                        # handle isotope and PAL
                        pal_count = (pal_mod_pos_array >= pos).sum()
                        for p in np.arange(0, pal_count + 1):
                            shift = p * 2 * (-MASS_H) / charge
                            ion_mz_np = ion_mz_n + shift
                            ratio = (1. / 3.) ** p
                            for i in isotope_array:
                                mz_list.append(ion_mz_np + i * MASS_H / charge)
                                if loss == 0 and charge == 1:
                                    int_list.append(tier1_int * ratio)
                                elif loss == 0 and charge > 1:
                                    int_list.append(tier2_int * ratio)
                                else:
                                    int_list.append(tier3_int * ratio)
    # sort
    order_key = np.argsort(mz_list)
    mz_array = np.array(mz_list)[order_key]
    int_array = np.array(int_list)[order_key]
    # camelcase_to_snake_case to histogram
    int_array = spectrum_histogram((mz_array, int_array), step=step)
    # pad satellite
    major_peak_array = np.where(int_array == tier1_int)[0]
    satellite_pos_array = np.concatenate([major_peak_array - 1, major_peak_array + 1])
    satellite_pos_array = satellite_pos_array[(satellite_pos_array >= 0) & (satellite_pos_array < int_array.shape[0])]
    # avoid conflict
    satellite_pos_set = set(satellite_pos_array)
    all_peak_set = set(np.where(int_array > 0)[0])
    satellite_pos_array = np.array(list(satellite_pos_set.difference(all_peak_set)))
    int_array[satellite_pos_array] = padding_int
    # generate dynamic region for data normalization
    return int_array


def seq_to_spectrum(seq, mod_dict={}, precursor_charge=2, isotope_array=np.arange(3), ion_type_list=['b', 'y'],
                    neutral_loss=['H2O', 'NH3'], mod_string=''):
    """
    camelcase_to_snake_case sequence to spectrum
    modification dictionary format: {0:+57, 3:+1}
    supports common modification table (common_modification_table) annotation
    use modification string to detect/handle pattern-generating modifications
    WARNING: HARD-CODED FOR CURRENT VERSION OF ISOTAG ('C', 2, 3) AND ('H', 2, 3)
    """
    # detect all mods
    mod_pos_array = np.zeros(len(mod_dict))
    for index, modPos in enumerate(mod_dict.keys()):
        mod_pos_array[index] = modPos
    # detect PAL mods
    pal_mod_pos_list = []
    for mod in mod_string.split(','):
        if len(mod) > 0:
            pos = int(mod.split('(')[0][1::])
            mod_name = mod.split('(')[1][0:-1]
            if '[PAL]' in mod_name:
                pal_mod_pos_list.append(pos)
    pal_mod_pos_array = np.array(pal_mod_pos_list)
    # detect neutral loss
    neutral_loss_list = [0]
    for n in neutral_loss:
        if n == 'H2O':
            neutral_loss_list.append(-mass_H2O)
        elif n == 'NH3':
            neutral_loss_list.append(-mass_NH3)
    # empty container
    ion_list = []
    for pos in range(1, len(seq)):
        for charge in np.arange(1, int(precursor_charge) + 1):
            for ionType in ion_type_list:
                if ionType in 'abc':
                    ion_mz = mass.fast_mass(seq[0:pos], ion_type=ionType, charge=charge)
                    for modPos in mod_pos_array[mod_pos_array < pos]:
                        ion_mz += mod_dict[modPos] / charge
                    # handle neutral loss
                    for loss in neutral_loss_list:
                        ion_mz_n = ion_mz + loss / charge
                        # handle isotope and PAL
                        pal_count = pal_mod_pos_array[pal_mod_pos_array < pos].sum()
                        for p in np.arange(0, pal_count + 1):
                            shift = p * 2 * (-MASS_H) / charge
                            ion_mz_np = ion_mz_n + shift
                            for i in isotope_array:
                                ion_list.append(ion_mz_np + i * MASS_H / charge)
                elif ionType in 'xyz':
                    ion_mz = mass.fast_mass(seq[pos::], ion_type=ionType, charge=charge)
                    for modPos in mod_pos_array[mod_pos_array >= pos]:
                        ion_mz += mod_dict[modPos] / charge
                    # handle neutral loss
                    for loss in neutral_loss_list:
                        ion_mz_n = ion_mz + loss / charge
                        # handle isotope and PAL
                        pal_count = pal_mod_pos_array[pal_mod_pos_array >= pos].sum()
                        for p in np.arange(0, pal_count + 1):
                            shift = p * 2 * (-MASS_H) / charge
                            ion_mz_np = ion_mz_n + shift
                            for i in isotope_array:
                                ion_list.append(ion_mz_np + i * MASS_H / charge)
                else:
                    continue
    # remove duplicates at 1e-4 precision
    unique_set = set([format(value, '.4f') for value in ion_list])
    unique_list = [float(value) for value in unique_set]
    theoretical_mz = np.array(sorted(unique_list))
    return theoretical_mz, np.ones(theoretical_mz.shape)


def prob_compare_spectrum(obs_sp, ref_sp, tol, unit):
    """
    probabilistic approach
    given a sample, knowing that random values are more close to zero
    use percentile to estimate the likelihood that a value is out of random
    return number of matches between observed spectrum and reference spectrum
    weighted by the percentile of the peak intensity
    """
    # basic setup
    match_pos_list = []
    i = 0
    j = 0
    mz1 = np.array(obs_sp[0])
    mz2 = np.array(ref_sp[0])
    # loop
    while i < mz1.shape[0] and j < mz2.shape[0]:
        if within_tol_prob(mz1[i], mz2[j], tol, unit):
            # case of a match
            # add count and advance both pointer
            match_pos_list.append(i)
            # don't advance data spectrum pointer when match
            # because one data can be assigned to two theoretical peaks
            # at least PD does so
            j += 1
        elif mz1[i] > mz2[j]:
            # case of a mismatch
            # advance the lagging pointer
            j += 1
        else:
            i += 1
    # sanity check
    if len(match_pos_list) == 0:
        return 0
    else:
        match_pos_array = np.array(match_pos_list)
        # here the rank is subtracted by one
        # in order to make the percentile start at zero
        # this is the prior belief: the lowest peak is definitely noise
        percentile_array = (stats.rankdata(obs_sp[1]) - 1) / float(len(obs_sp[1]))
        match_percentile_array = percentile_array[match_pos_array]
        return match_percentile_array.sum()


def spectrum_histogram(spectrum, step=sequest_default_bin_size):
    mz_array = np.array(spectrum[0])
    int_array = np.array(spectrum[1])
    mz_start = sequest_mz_start
    mz_end = sequest_mz_end
    grid_length = int((mz_end - mz_start) / step)
    # calculate new index
    new_index_array = ((mz_array - mz_start) / step + sequest_index_offset * step).astype(int)
    # 1. restrict to in-space index
    # 2. keep only largest one upon collision
    mask = (new_index_array >= 0) & (new_index_array < grid_length)
    new_index_array = new_index_array[mask]
    int_array = int_array[mask]
    # 2. keep only largest one upon collision
    order_key = np.argsort(int_array)
    new_index_array = new_index_array[order_key]
    int_array = int_array[order_key]
    bin_int_array = np.zeros(grid_length)
    bin_int_array[new_index_array] = int_array
    return bin_int_array


def fast_xcorr_preprocessing(int_array, step):
    int_array2 = np.zeros(int_array.shape[0])
    step_max_delay = int(sequest_max_delay / step)
    for tau in range(-step_max_delay, step_max_delay + 1):
        if tau > 0:
            int_array2[tau:] += int_array[:-tau]
        elif tau < 0:
            int_array2[:tau] += int_array[-tau:]
    return int_array - int_array2 / (2 * step_max_delay)


def early_mass_check(obs_precursor_mass, candidate_bare_mass, candidate_sequence, mod_rule, aa_count, tol, unit):
    """
    if observed precursor mass is larger than candidate sequence with all positive modifications
    of smaller than candidate sequence with all negative modifications
    return False to skip the search
    """
    # check negative mass first, since most mods are positive
    min_possible_mass = candidate_bare_mass
    for mod in mod_rule['negativeDynamic'] + mod_rule['negativeStatic']:
        monoisotopic_mass_change = mod[1]
        aa = mod[2]
        count = int(aa_count[proteogenicAA.index(aa)])
        if count > 0:
            min_possible_mass += count * monoisotopic_mass_change
    if not within_tol_prob(obs_precursor_mass, min_possible_mass, tol, unit) and min_possible_mass > obs_precursor_mass:
        return False
    max_possible_mass = candidate_bare_mass
    for mod in mod_rule['positiveDynamic'] + mod_rule['positiveStatic']:
        monoisotopic_mass_change = mod[1]
        aa = mod[2]
        count = int(aa_count[proteogenicAA.index(aa)])
        if count > 0:
            max_possible_mass += count * monoisotopic_mass_change
    if not within_tol_prob(obs_precursor_mass, max_possible_mass, tol, unit) and max_possible_mass < obs_precursor_mass:
        return False
    # finally
    return True


def search_func(task_dict):
    """ basic search function for parallel processing"""
    # unpack input
    ms2_path = task_dict['ms2_path']
    db_path = task_dict['databasePath']
    mod_path = task_dict['modificationPath']
    output_folder_path = task_dict['output_path']
    prior = task_dict['prior']
    use_mass_prior = task_dict['use_mass_prior']
    # load data
    with open(ms2_path, 'r') as file_pointer:
        spectrum = json.load(file_pointer)
    with open(mod_path, 'r') as file_pointer:
        mod_rule = json.load(file_pointer)
    # extract information
    ms2_scan_number = int(ms2_path.split('/')[-1][0:-5])
    obs_spectrum = (spectrum['m/zList'], spectrum['intensityList'])
    obs_precursor_mz = spectrum['obsPrecursorMz']
    obs_precursor_charge = spectrum['obsPrecursorCharge']
    # check if prior is used
    if prior is None:
        empirical_precursor_mass = obs_precursor_mz * obs_precursor_charge - obs_precursor_charge * MASS_H
        search_dict = {'isotag': [1, empirical_precursor_mass],
                       'regular': [1, empirical_precursor_mass]}
    else:
        search_dict = {'isotag': [prior[0], prior[2]],
                       'regular': [prior[1], prior[3]]}
    # get tolerance
    if spectrum['massAnalyzer'] == 'ITMS':
        precursor_tolerance = [10, 'ppm']
        fragment_tolerance = [0.6, 'Da']
    elif spectrum['massAnalyzer'] == 'FTMS':
        precursor_tolerance = [10, 'ppm']
        fragment_tolerance = [0.02, 'Da']
    # container
    hit_dict = dict()
    hit_index = 0
    top_hit_index = float('nan')
    total_peak_db = 0
    total_match_db = 0
    # default output if no hit
    # ['ms2_scan_number',
    #   'uniprotAccessions', 'peptideSequence', 'modifications',
    #   'peptideMass', 'nPeak', 'nMatch', 'massPrior', 'PALPrior']
    default_search_result = [ms2_scan_number] + [''] * 3 + [float('nan')] * 6
    # search
    for searchType in search_dict.keys():
        obs_precursor_mass = search_dict[searchType][1]
        pal_prior = search_dict[searchType][0]
        file_pointer = open(db_path, 'r')
        db_reader = csv.reader(file_pointer, delimiter=',')
        for line in db_reader:
            peptide_sequence = line[0]
            bare_mass = float(line[1])
            uniprot_accession = line[2]
            # early checkpoint
            if not early_mass_check(obs_precursor_mass, bare_mass, peptide_sequence, mod_rule,
                                    precursor_tolerance[0], precursor_tolerance[1]):
                continue
            # handle static modifications
            for mod_dict, mod_string in modification_generator(mod_rule, peptide_sequence):
                # filter based on search type
                if (searchType == 'isotag' and mod_string.count('[PAL]') != 1) \
                        or (searchType == 'regular' and '[PAL]' in mod_string):
                    continue
                # merge modification dictionaries
                peptide_mass = bare_mass + sum(mod_dict.values())
                # get prior probability on precursor mass
                mass_prior = within_tol_prob(obs_precursor_mass, peptide_mass, precursor_tolerance[0], precursor_tolerance[1],
                                            prob=use_mass_prior)
                # get prior probability on modification
                if mass_prior > 0 and pal_prior > 0:
                    theoretical_spectrum = seq_to_spectrum(peptide_sequence, precursor_charge=obs_precursor_charge,
                                                          mod_dict=mod_dict, mod_string=mod_string)
                    n_peak = theoretical_spectrum[0].shape[0]
                    n_match = prob_compare_spectrum(obs_spectrum, theoretical_spectrum,
                                                   fragment_tolerance[0], fragment_tolerance[1])
                    # record
                    total_peak_db += n_peak
                    total_match_db += n_match
                    hit_dict[hit_index] = [ms2_scan_number, uniprot_accession, peptide_sequence,
                                           mod_string, obs_precursor_mass, n_peak, n_match, mass_prior, pal_prior]
                    hit_index += 1
        file_pointer.close()
    # scoring
    # in case there's no hit in database
    if total_match_db == 0:
        return default_search_result
    else:
        # take out info
        key_array = np.array([key for key in hit_dict.keys()])
        n_match_array = np.array([hit_dict[key][6] for key in hit_dict.keys()])
        n_peak_array = np.array([hit_dict[key][5] for key in hit_dict.keys()])
        mass_prior_array = np.array([hit_dict[key][7] for key in hit_dict.keys()])
        pal_prior_array = np.array([hit_dict[key][8] for key in hit_dict.keys()])
        # probability that a hit is not random
        score_1_p_value = np.empty(key_array.shape)
        for index, key in enumerate(key_array):
            total_match_hit = hit_dict[key][6]
            total_peak_hit = hit_dict[key][5]
            random_freq = total_peak_hit * total_match_db / total_peak_db
            score_1_p_value[index] = interp_poisson_cdf(total_match_hit - 1, random_freq)
        score_p_value = 1 - score_1_p_value
        # adapt a numerically stable way to normalize the compound probability
        score_log_p = np.log(score_p_value)
        score_log_1_p = np.log(score_1_p_value)
        # sanity check: is there perfect hit(s)?
        perfect_hit_mask = np.logical_not(np.isfinite(score_log_p))
        count_perfect_hit = np.sum(perfect_hit_mask)
        if count_perfect_hit > 0:
            best_key = key_array[perfect_hit_mask][np.argmax(n_match_array[perfect_hit_mask] / n_peak_array[perfect_hit_mask] \
                                                         * pal_prior_array[perfect_hit_mask])]
            return hit_dict[best_key] + [1]
        # take numerically stable way, through np.logsumexp
        # also separately handle those (1-p) == 0
        one_mask = np.isfinite(score_log_1_p)
        if np.sum(one_mask) == 0:
            return default_search_result
        zero_mask = np.logical_not(one_mask)
        log_likelihood_array = np.empty(n_match_array.shape)
        likelihood_array = np.empty(n_match_array.shape)
        # assign zero likelihood to those whose 1-pValue is too small
        likelihood_array[zero_mask] = 0
        # for the rest, calculate compound probability
        # np.where(oneMask)[0] returns an array of indices where
        # the boolean in oneMask is True
        log_mass_prior_array = np.log(mass_prior_array)
        log_mod_prior_array = np.log(pal_prior_array)
        for index in np.where(one_mask)[0]:
            others_mask = np.copy(one_mask)
            others_mask[index] = False
            log_likelihood_array[index] = score_log_1_p[index] \
                                        + log_mod_prior_array[index]
        log_normalizer = misc.logsumexp(log_likelihood_array[one_mask])
        likelihood_array[one_mask] = np.exp(log_likelihood_array[one_mask] - log_normalizer)
        # save whole list to disk for later EM algorithm
        output_list = []
        for index, key in enumerate(key_array):
            likelihood = likelihood_array[index]
            output_list.append(hit_dict[key] + [likelihood])
        output_df = pd.DataFrame.from_records(output_list,
                                             columns=['ms2_scan_number', 'uniprotAccession', 'peptideSequence',
                                                      'mod_string', 'obs_precursor_mass', 'nPeak', 'nMatch', 'massPrior',
                                                      'PALPrior', 'likelihood'])
        output_df.to_csv(os.path.join(output_folder_path, 'PSM_' + str(ms2_scan_number) + '.csv'))
        # now report top hit
        top_hit_pos = np.argmax(likelihood_array)
        top_hit_key = key_array[top_hit_pos]
        return hit_dict[top_hit_key] + [likelihood_array[top_hit_pos]]


def search_func_sequest(task_dict):
    """
    basic search function for parallel processing
    """

    # unpack input
    ms2_path = task_dict['MS2Path']
    db_path_list = task_dict['databasePathList']
    output_folder_path = task_dict['output_path']
    pal_prior = task_dict['prior']
    use_mass_prior = task_dict['useMassPrior']
    mod_path = task_dict['modPath']
    # load data
    with open(mod_path, 'r') as file_pointer:
        mod_rule = json.load(file_pointer)
    with open(ms2_path, 'r') as file_pointer:
        spectrum = json.load(file_pointer)
    # extract information
    ms2_scan_number = int(ms2_path.split('/')[-1][0:-5])
    obs_spectrum = (spectrum['m/zList'], spectrum['intensityList'])
    obs_precursor_mz = spectrum['obsPrecursorMz']
    obs_precursor_charge = spectrum['obsPrecursorCharge']
    # check if prior is used
    if pal_prior is None:
        empirical_precursor_mass = obs_precursor_mz * obs_precursor_charge - obs_precursor_charge * MASS_H
        search_dict = {'isotag': [1, empirical_precursor_mass],
                       'regular': [1, empirical_precursor_mass]}
    else:
        search_dict = {'isotag': [pal_prior[0], pal_prior[2]],
                       'regular': [pal_prior[1], pal_prior[3]]}
    # get tolerance
    if spectrum['massAnalyzer'] == 'ITMS':
        precursor_tolerance = [10, 'ppm']
        fragment_tolerance = [0.6, 'Da']
    elif spectrum['massAnalyzer'] == 'FTMS':
        precursor_tolerance = [10, 'ppm']
        fragment_tolerance = [0.02, 'Da']
    # container
    hit_index = 0
    hit_dict = dict()
    # default output if no hit
    # ['ms2_scan_number',
    #   'uniprotAccessions', 'peptideSequence', 'modifications',
    #   'peptideMass', 'XCorr', 'massPrior', 'PALPrior']
    default_search_result = [ms2_scan_number] + [''] * 3 + [float('nan')] * 4
    # preprocess observed spectrum
    try:
        mz_start = np.min(obs_spectrum[0])
        mz_end = np.max(obs_spectrum[0])
        step = fragment_tolerance[0]
        obs_int_array = spectrum_histogram(obs_spectrum, step=step)
        obs_int_array = sequest_preprocessing(obs_int_array, mz_start=sequest_mz_start, mz_end=sequest_mz_end, step=step)
        obs_int_array = fast_xcorr_preprocessing(obs_int_array, step=step)
    except:
        print('problem pre-processing MS2 spectrum of scan {}'.format(ms2_scan_number))
        return default_search_result
    # search
    for searchType in search_dict.keys():
        obs_precursor_mass = search_dict[searchType][1]
        pal_prior = search_dict[searchType][0]
        # get upper and lower bounds of precursor mass
        # here we use double tolerance because when creating the database we did not add tolerance
        if precursor_tolerance[1] == 'ppm':
            precursor_mass_upper_bound = obs_precursor_mass * (1 + 2 * precursor_tolerance[0] * 1e-6)
            precursor_mass_lower_bound = obs_precursor_mass * (1 - 2 * precursor_tolerance[0] * 1e-6)
        else:
            precursor_mass_upper_bound = obs_precursor_mass + 2 * precursor_tolerance[0]
            precursor_mass_lower_bound = obs_precursor_mass - 2 * precursor_tolerance[0]
        for DBPath in db_path_list:
            # connect to SQLite database
            conn = sqlite3.connect(DBPath)
            cur = conn.cursor()
            query = cur.execute('SELECT * FROM expandedDB WHERE ? <= mass_upper_bound AND ? >= mass_lower_bound',
                                [precursor_mass_lower_bound, precursor_mass_upper_bound])
            for line in query:
                # expand modification space
                uniprot_accession = line[0]
                peptide_sequence = line[1]
                base_mass = mass.fast_mass(sequence=peptide_sequence)
                for modDict, mod_string in modification_generator(mod_rule=mod_rule, seq=peptide_sequence,
                                                                  mass_lower_bound=precursor_mass_lower_bound - base_mass,
                                                                  mass_upper_bound=precursor_mass_upper_bound - base_mass):
                    # check PAL criteria
                    violate_pal = (searchType == 'isotag' and '[PAL]' not in mod_string) or (
                                searchType == 'regular' and '[PAL]' in mod_string)
                    # check mass
                    compound_mass = base_mass + sum(modDict.values())
                    mass_prior = within_tol_prob(obs_precursor_mass, compound_mass, precursor_tolerance[0],
                                                 precursor_tolerance[1],
                                                 prob=use_mass_prior)
                    if mass_prior > 0 and not violate_pal and pal_prior > 0:
                        theoretical_int_array = seq_to_spectrum_sequest(peptide_sequence,
                                                                        precursor_charge=obs_precursor_charge,
                                                                        mod_string=mod_string, step=step)
                        x_corr = theoretical_int_array.dot(obs_int_array)
                        # record
                        hit_dict[hit_index] = [ms2_scan_number, uniprot_accession, peptide_sequence,
                                               mod_string, obs_precursor_mass, x_corr, mass_prior, pal_prior]
                        hit_index += 1
            conn.close()
    # scoring
    # in case there's no hit in database
    if len(hit_dict) == 0:
        return default_search_result
    else:
        # take out info
        key_array = np.array([key for key in hit_dict.keys()])
        x_corr_array = np.array([hit_dict[key][5] for key in hit_dict.keys()])
        # save whole list to disk for later EM algorithm
        output_list = []
        for index, key in enumerate(key_array):
            output_list.append(hit_dict[key])
        output_df = pd.DataFrame.from_records(output_list,
                                             columns=['ms2_scan_number', 'uniprotAccession', 'peptideSequence',
                                                      'mod_string', 'obs_precursor_mass', 'XCorr', 'massPrior',
                                                      'PALPrior'])
        output_df.to_csv(os.path.join(output_folder_path, 'PSM_' + str(ms2_scan_number) + '.csv'))
        #        print(outputDF.sort_values('XCorr', ascending = False).head())
        #        print(outputDF.sort_values('XCorr', ascending = False).iloc[0])
        # now report top hit
        top_hit_pos = np.argmax(x_corr_array)
        top_hit_key = key_array[top_hit_pos]
        return hit_dict[top_hit_key]


def unit_fuzzy_match(input_tuple):
    """
    This multiprocessing function was used in isostamp_engine. Replaced with a single thread simpler method
    for pos, adj_set in tqdm.tqdm(self.worker.imap(self.unit_fuzzy_match,
                                  [(i, ms2) for i in ms2.index]), total=ms2.shape[0]):
    adj_dict[pos] = adj_set

    only search within +/- 1 min retention time
    then utilize the power of depth-first search to find the connection
    """
    # unpack inputs
    pos1 = input_tuple[0]
    asso_df = input_tuple[1]
    # extract info
    mz1 = asso_df.loc[pos1, 'obs_precursor_mz']
    charge1 = asso_df.loc[pos1, 'obs_precursor_charge']
    rt1 = asso_df.loc[pos1, 'ms2_rt']
    # within +/- 1 minute retention time
    relevant_df = asso_df.loc[(asso_df['ms2_rt'] - rt1).abs() < 1]

    # get mask
    mask_mz = relevant_df['obs_precursor_mz'].apply(
        lambda x: within_tol_prob(x, mz1, param.MS1_TOLERANCE, 'ppm'))
    mask_charge = relevant_df['obs_precursor_charge'] == charge1
    mask_rt = abs(relevant_df['ms2_rt'] - rt1) < (2 * param.MS2_RT_WINDOW)
    mask_overall = mask_mz & mask_charge & mask_rt

    # get duplicated index, this will definitely include itself
    matching_idx = set(relevant_df[mask_overall].drop(pos1).index.unique())
    # return in the format of a tuple
    return pos1, matching_idx


def compare_spectrum(spectrum1, spectrum2):
    if len(spectrum1[0]) != len(spectrum2[0]):
        return False
    else:
        for i in range(len(spectrum1[0])):
            if spectrum1[0][i] != spectrum2[0][i]:
                return False
            elif spectrum1[1][i] != spectrum2[1][i]:
                return False
    return True


def modeldict_to_spectrum(model_dict):
    # camelcase_to_snake_case dictionary to tuple of numpy array
    mz_array = np.array([float(string) for string in model_dict.keys()])
    int_array = np.array([model_dict[k] for k in model_dict.keys()])
    order_key = np.argsort(mz_array)
    mz_array = mz_array[order_key]
    int_array = int_array[order_key]
    return mz_array, int_array


def merge_model_dict(md1, md2):
    md1_key_set = set(md1.keys())
    md2_key_set = set(md2.keys())
    key_both = md1_key_set.intersection(md2_key_set)
    key1_only = md1_key_set.difference(md2_key_set)
    key2_only = md2_key_set.difference(md1_key_set)
    md3 = dict()
    for key in key_both:
        md3[key] = md1[key] + md2[key]
    for key in key1_only:
        md3[key] = md1[key]
    for key in key2_only:
        md3[key] = md2[key]
    return md3


def spectrum_to_modeldict(spectrum):
    mz_array = spectrum[0]
    int_array = spectrum[1]
    model_dict = dict()
    for gridMz, gridInt in zip(mz_array, int_array):
        grid_mz_string = format(gridMz, '.4f')
        if grid_mz_string not in model_dict.keys():
            model_dict[grid_mz_string] = gridInt
        else:
            model_dict[grid_mz_string] += gridInt
    return model_dict


def isostamp_database_combine(worker, psmdf, ms1df_path):
    """superceded with a v2 that doesn't use multiprocessing. one thread takes < 1 minute"""
    print_messagebox('Merging IsoStamp result with Sequest result...')
    score_dict = {}

    # psmdf['isostamp_score'] = psmdf.apply(lambda x: find_isostamp_score((x.first_scan, x.is_modified, ms1df_path)), 1)
    find_isostamp_score_param = [(ms2_scan_number, modified, ms1df_path) for ms2_scan_number, modified in
                                 zip(psmdf['first_scan'].values, psmdf['is_modified'].values)]
    for ms2_scan_number, score in tqdm.tqdm(worker.imap(find_isostamp_score, find_isostamp_score_param),
                                            total=psmdf.shape[0]):
        score_dict[ms2_scan_number] = score

    psmdf['isostamp_score'] = psmdf['first_scan'].apply(lambda x: score_dict[x])

    logging.info('step 4. apply 40% normalized isostamp score and %60 normalized xcorrelation weighting')
    psmdf['normalized_isostamp_score'] = scale(psmdf['isostamp_score'])
    psmdf['normalized_xcorr'] = scale(psmdf['xcorr'])
    psmdf['combo_score'] = 0.4 * psmdf['normalized_isostamp_score'] + 0.6 * psmdf['normalized_xcorr']

    logging.info('step 5. re-rank by Target-Decoy Approach')
    psmdf['combo_fdr'] = calculate_q_value(psmdf['tda'].values, psmdf['combo_score'])

    return psmdf


def find_isostamp_score(input_tuple):
    ms2_scan_number = input_tuple[0]
    modified = input_tuple[1]
    ms1_df_path = input_tuple[2]
    ms1_df = pd.read_csv(ms1_df_path, index_col=0)
    # hacky naming change to handle inconsistent naming
    ms1_df.columns = [camelcase_to_snake_case(col) for col in ms1_df.columns]
    for index_is in ms1_df.index:
        ms2_scan_number_list = [int(s) for s in ms1_df.loc[index_is, 'ms2_scan_number'].split('+')]
        if ms2_scan_number in ms2_scan_number_list:
            if modified:
                return ms2_scan_number, ms1_df.loc[index_is, 'isotag_propensity']
            else:
                return ms2_scan_number, ms1_df.loc[index_is, 'regular_propensity']
    return ms2_scan_number, float('nan')


def gauss_function(x, a, x0, sigma):
    """Possible source? https://stackoverflow.com/questions/19206332/gaussian-fit-for-python"""
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))


def gaussian_mapping(mz_array, intensity_array, precision=param.PRECISION, z_sigma=param.Z_SIGMA):
    """
    Empirical sum of an array of gaussian (normal) distributions.
    A normal distibution can be defined by two parameters (mean and standard deviation), but this function
    calculates the empirical results

    :param mz_array: np.array
    :param intensity_array: np.array
    :param precision: int, number of decimals to calculate the emp
    :return:
    """

    model_dict = dict()

    for mz, intensity in zip(mz_array, intensity_array):
        sigma = mz * z_sigma

        for grid_mz in np.arange(round(mz - 3 * sigma, precision), round(mz + 3 * sigma, precision), 10**-precision):
            key = round(grid_mz, precision)
            value = gauss_function(grid_mz, intensity, mz, sigma)
            # add a requirement of a non-zero value
            if key in model_dict.keys() and value:
                model_dict[key] += value
            else:
                model_dict[key] = value
    return model_dict


def gaussian_1d(ts, sigma=1):
    """begin replacing gaussian mapping with scipy implementation
    https://stackoverflow.com/questions/32900854/how-to-smooth-a-line-using-gaussian-kde-kernel-in-python-setting-a-bandwidth"""
    from scipy import ndimage
    return pd.Series(ndimage.gaussian_filter1d(ts, sigma=sigma), index=ts.index)