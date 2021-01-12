# number of CPU cores to run. If unsure, set None to use all cores
MULTI_PROCESSES = 8
# unsure of this parameter meaning
ISOSTAMP = {'element_type': 'C', 'element_count': 2, 'ratio': 3}
# unsure of this parameter meaning
ISOLATION_WINDOW_OFFSET = 0
# precision of the mz intensity accumulation. Going from 2 to 3 precision increases runtime by 5x
PRECISION = 2
# relative to maximum intensity. exclude intensities below this threshold
INTENSITY_MIN = 0.01

# Group MS2 with similar characteristics
# dynamic exclusion rt (runtime, minutes)
MS2_RT_WINDOW = 0.5

# PPM tolerance for feature merging in MS1
MS1_TOLERANCE = 20 * 1e-6

# brainpy yield parameters. Filters MS1 based on MS2 information
MS1_RT_WINDOW = 0.5
MS1_PRECURSOR_WINDOW = 0.02
MS1_INTENSITY_FLOOR = 5000
PRECURSOR_MZ_WINDOW = 4
INTENSITY_FLOOR_COEF = 0.5  # 1.0 sets the floor as the mean intensity
INTENSITY_FLOOR = 100000

# brainpy model peak width (standard deviation) in PPM
# used when fitting failed
Z_SIGMA = 20. / 3. * 1e-6
# minimum R-squared score required for precursor mass correction in brainpy fitting
# otherwise default to no correction
PRECURSOR_PROPENSITY_FLOOR = 0.9

# amino acid residue mass dictionary. Here for reference, but don't alter unless the universe changes.
MASS_H = 1.00794
MASS_D = 2.01410
MASS_C12 = 12.0
MASS_C13 = 13.003355