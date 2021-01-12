import logging
import os
import time
from isostamp_engine import IsoStampEngine
from multiprocessing import Pool
from parameters import MULTI_PROCESSES
logging.basicConfig(level=logging.INFO)


def run(directory, input_folder, output_folder, nprocesses):
    input_path = os.path.join(directory, input_folder)
    output_path = os.path.join(directory, output_folder)
    profile_files = sorted([p for p in os.listdir(input_path) if p.endswith('_profile.mzML')])

    run_start = time.time()
    worker = Pool(processes=nprocesses)
    logging.info(f'CPU detected, {worker._processes} threads specified')
    engine = IsoStampEngine(worker, output_path)

    # outputs a _ms1df.csv for each profile file
    for idx, profile_file in enumerate(profile_files):
        iteration_start = time.time()
        logging.info(f'isostamp run: {idx}, profile mzML file path: {profile_file}')

        # main script
        ms1_df = engine.isostamp_run(os.path.join(input_path, profile_file))

        # assumes an input file structure of a [unique_file_identifier]_profile.mzML
        # saves an output csv of [unique_file_indentifier]_ms1df.csv
        full_output_file = f'{profile_file.split("_profile")[0]}_ms1df.csv'
        ms1_df.to_csv(os.path.join(output_path, full_output_file), index=False)

        logging.info(f'{full_output_file} exported with {ms1_df.shape} shape')
        logging.info(f'iteration completed in {((time.time() - iteration_start) / 60.):.1f} minutes')

    logging.info(f"""{len(profile_files)} files IsoStamped. 
                     Completed in {((time.time() - run_start) / 3600.):.1f} hours""")


if __name__ == '__main__':
    directory_ = '/hd/data/hw202'
    input_folder_ = 'input'
    output_folder_ = 'output'

    run(directory_, input_folder_, output_folder_, MULTI_PROCESSES)
