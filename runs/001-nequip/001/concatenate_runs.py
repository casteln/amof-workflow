# Imports
import os
import sys
import glob

import itertools 
import logging
import configparser
import pathlib

import sadi.files.operation as sop
import sadi.files.lammps as slmp
# working_dir = os.path.dirname(sys.argv[0]) # doesn't work in terminal
working_dir = os.path.dirname(os.path.realpath(__file__))
os.chdir(working_dir) # change working dir as where the python script is

config = configparser.ConfigParser()
config.read('run_lammps.ini')

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("run_lammps.log"), 
        logging.StreamHandler()
    ]  
)

logger = logging.getLogger(__name__)

logger.warning('---------------------')

logger.info('Concatenate runs python script')

# List of paths
logger.info('working directory: ' + working_dir)


# Parameters

constant = dict(config['constant_parameters'])

variable_dict = {key : config['variable_parameters'][key].replace(" ", "").split(",") for key in config['variable_parameters']} # rm spaces and split with ,

if config['configuration']['variation.style'] == 'combination':
    product = list(itertools.product(*variable_dict.values())) 
    variable_list = [dict(zip(variable_dict.keys(), v)) for v in product]
elif config['configuration']['variation.style'] in ['linear', 'successive']:
    variable_len = set([len(value) for value in variable_dict.values()])
    if len(variable_len) != 1:
        logger.warning("linear variation style - not all parameters are provided for every run")
    N = min(variable_len)
    set_of_var = [[value[j] for value in variable_dict.values()] for j in range(N)]
    variable_list = [dict(zip(variable_dict.keys(), v)) for v in set_of_var]

def run_main(constant, variable, id):
    """
    input: 
    two dicts containing the constant and variable parameters
    id: int providing the id used for naming"""
    parameters = {**constant, **variable} # merge two dicts

    run_exp = format(id, '03d') + "_" + '_'.join([key+val for key, val in variable.items()])    
    parameters['run_exp_name'] = '_'.join([key+val for key, val in variable.items()])
    parameters['run_exp'] = id

    i = 0
    while i < 100:
        try:
            nb_restarts = int(config['configuration'][f'num.restart{i}'])
            i += 1
        except KeyError:
            break
    if i == 100:
        logger.error('No num.restart found in config.ini')

    folders = [f'restart{i}' for i in range(nb_restarts)] # sorted in run order
    files_to_zip_pattern = ['dump_vel', 'log.slurm'] 
    clean_xyz_pattern = ['*.xyz'] 

    files_to_concatenate_pattern = [*clean_xyz_pattern,'log.run.*'] # will turn into a list for every run_exp depending on what's present with this pattern
    files_to_concatenate_and_compress_pattern = files_to_concatenate_pattern
    # will subsequently be zipped as well

    # add command line options
    compress_files = False
    decompress_files = False
    try:
        if sys.argv[1] == 'compress':
            compress_files = True
    except:
        pass
    try:
        if sys.argv[1] == 'decompress':
            decompress_files = True
    except:
        pass


    def file_list_from_pattern(pattern_list):
        file_list = []
        for file in pattern_list:
            for folder in folders:
                for ext in ['.gz','']:
                    for filename in glob.glob(f'{run_exp}/{folder}/{file}{ext}'):
                        file_list.append(pathlib.Path(filename).name.replace('.gz', ''))
        return list(set(file_list))

    files_to_zip = file_list_from_pattern(files_to_zip_pattern)
    files_to_concatenate = file_list_from_pattern(files_to_concatenate_pattern)
    clean_xyz = file_list_from_pattern(clean_xyz_pattern)
    files_to_concatenate_and_compress = file_list_from_pattern(files_to_concatenate_and_compress_pattern)

    # concatenate files and decompress what's needed
    decompressed_files = []
    for file in files_to_concatenate:
        filenames = []
        for folder in folders:
            for filename in glob.glob(f'{run_exp}/{folder}/{file}.gz'):
                sop.decompress(filename[:-3], remove=decompress_files) # remove '.gz'
                decompressed_files.append(filename[:-3])
            for filename in glob.glob(f'{run_exp}/{folder}/{file}'):
                filenames.append(filename)
        logger.info("concatenate %s", file)
        sop.concatenate(filenames, f'{run_exp}/{file}')

    # clean files
    for file in clean_xyz:
        logger.info("remove_duplicate_timesteps %s", file)
        slmp.remove_duplicate_timesteps(f'{run_exp}/{file}')           

    # compress files that were concatenated (or remove uncompressed file if already compressed)
    for file in files_to_concatenate_and_compress:
        for folder in folders:
            for filename in glob.glob(f'{run_exp}/{folder}/{file}'):
                if (not decompress_files) and (compress_files or filename in decompressed_files):
                    sop.compress(filename, remove_if_exists=True)

    # compress other files
    for folder in folders:
        for file in files_to_zip:
            for filename in glob.glob(f'{run_exp}/{folder}/{file}*'):
                if compress_files and filename[-3:] != '.gz':
                    sop.compress(filename)
                if decompress_files and filename[-3:] == '.gz':
                    sop.decompress(filename[:-3])
     

logger.info('Serie number #' + config['system_info']['run_serie'])

i = 1 #id
for variable in variable_list:
    run_main(constant, variable, i)
    i += 1