# Imports

# Standard libraries
import os

import atomman as am

import json
import shutil
import itertools 
import pathlib

import pandas as pd
import numpy as np

# Logging
import logging

import configparser

# working_dir = os.path.dirname(sys.argv[0]) # doesn't work in terminal
working_dir = os.path.dirname(os.path.realpath(__file__))
os.chdir(working_dir) # change working dir as where the python script is

config = configparser.ConfigParser()
config.read('run_lammps.ini')


# create log file with same name as python file
logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(os.path.basename(__file__)[:-3]+".log"), 
        logging.StreamHandler()
    ]  
)

logger = logging.getLogger(__name__)

logger.warning('---------------------')
logger.info('Running lammps python script')



# List of paths
logger.info('working directory: ' + working_dir)

def format_time(t, dt):
    """
    input: t as float in ps, dt timestep as float in ps
    WARNING: different than for real units where dt is in fs
    output: t as str in lammps units
    """
    # return round(1000 * t / float(dt)) # for real units
    return round(t / float(dt)) # for metal units


def read_template(template_name, parameters):
    """
    template_name: name of template file used in .ini file
    parameters: dic containing all the necessary input for templating
    """
    with open(config['path']['template'] + config['filename'][template_name]) as f:
        template = f.read()
    script = am.tools.filltemplate(template, parameters, '{', '}') + '\n'  
    return script


def md_run(thermodynamic_ensemble, output_name, T_i, T_f, run_time, simulation_time, parameters, run_description, multiple_restarts = None, volume_change = None):
    """
    Args: 
        thermodynamic_ensemble: 'nve', 'nvt', or 'npt'
        run_time, simulation_time floats in ps
        multiple_restarts: None if no restarts, int representing number of restarts
        volume_change: None if no deformation, float representing relative volumic change in %
    Returns:
        code, simulation_time updated
    """
    simulation_time += run_time
    dt = float(parameters['timestep'])
    fill_var = {'output_name': output_name, 'run_time': format_time(run_time, dt), 'T_i': str(T_i), 'T_f': str(T_f)}    
    fill_var = {**fill_var, **parameters} # adds parameters such as tdamp
    template_name = thermodynamic_ensemble
    if multiple_restarts is not None:
        template_name += '_multiple_restarts'
        fill_var['freq_restart'] = format_time(run_time/multiple_restarts, dt)
    if volume_change is not None:
        template_name += '_volume_change'
        V = 1 + 0.01 * float(volume_change)
        fill_var['linear_scale'] = np.around(V**(1/3), 8)
    code = read_template(template_name, fill_var)   
    run_description.append({'output_name': output_name, 'first_step': format_time(simulation_time-run_time, dt), 'last_step': format_time(simulation_time, dt), 'thermodynamic_ensemble': thermodynamic_ensemble, 'initial_temp': T_i, 'final_temp': T_f})       
    return code, simulation_time, run_description   


def nve(output_name, run_time, simulation_time, parameters, run_description):
    """
    simple wrapper for md_run
    input: run_time, simulation_time floats in ps
    output: code, simulation_time updated
    """
    return md_run('nve', output_name, 0, 10, run_time, simulation_time, parameters, run_description)

def nvt(output_name, T_i, T_f, run_time, simulation_time, parameters, run_description, multiple_restarts = None, volume_change = None):
    """
    simple wrapper for md_run
    input: run_time, simulation_time floats in ps
    output: code, simulation_time updated
    """
    return md_run('nvt', output_name, T_i, T_f, run_time, simulation_time, parameters, run_description, multiple_restarts, volume_change)
 
def npt(output_name, T_i, T_f, run_time, simulation_time, parameters, run_description, multiple_restarts = None):
    """
    simple wrapper for md_run
    input: run_time, simulation_time floats in ps
    output: code, simulation_time updated
    """
    return md_run('npt', output_name, T_i, T_f, run_time, simulation_time, parameters, run_description, multiple_restarts)  


def write_lammps_input(folder_path, parameters, run_part, run_part_name, first_run = False, next_folder = None, previous_folder = None, first_run_exp = False):
    """
    Args:
        folder_path: path to folder
        parameters: dic containing all the necessary input for templating
        run_part, run_part name: int and string describing run_part (e.g. '2' and 'melt')
        first_run: boolean, True iff first run (no restart, reset_timestep, minimize)
        next_folder: string, name of next_folder
        previous_folder: string, name of previous_folder. Not equal to None if a new serie of lammps runs is ran starting from this. 
    """
    # copy raw_files
    os.chdir(working_dir)
    folder_path = pathlib.Path(folder_path)
    logger.info('Write folder %s', folder_path)
    rawfiles_path = pathlib.Path(config['path']['rawfiles'])
    folder_path.mkdir(parents=True, exist_ok=True)
    shutil.copy(rawfiles_path / f"{parameters['model']}.pth", folder_path)
    run_on = config['configuration']['run_on']

    parameters['code.minimize'] = ''
    
    parameters['initial_restart_file'] = 'initial.restart'
    parameters['final_restart_file'] = 'final.restart'

    if first_run:
        # Commented restart and use data for this starting point
        # for file in rawfiles_path.glob(f"equilibrate.*_300K.restart.*"):
        #     restart_file = pathlib.Path(file).name
        #     shutil.copy(file, folder_path)
        #     parameters['initial_restart_file'] = restart_file
        # parameters['code.read_file'] = read_template('read_restart', parameters)
        for file in rawfiles_path.glob(f"crystal.lmp"): 
            data_file = pathlib.Path(file).name
            shutil.copy(file, folder_path)
            parameters['initial_data_file'] = data_file
        parameters['code.read_file'] = read_template('read_data', parameters)        
        parameters['initial_restart_file'] = 'initial.restart' # set up for next_job
    else:
        parameters['code.read_file'] = read_template('read_restart', parameters)
    
    if first_run_exp or first_run:
        parameters['code.reset_timestep'] = read_template('reset_timestep', parameters)
    else:
        parameters['code.reset_timestep'] = ''      


    if previous_folder is not None:
        parameters['previous_folder'] = previous_folder
        parameters['code.previous_job'] = read_template('previous_job', parameters)
    else:
        parameters['code.previous_job'] = ''

    if next_folder is not None:
        parameters['next_folder'] = next_folder
        if run_on in ['occigen', 'zay']:
            parameters['code.next_job'] = read_template('next_job', parameters)
        elif run_on == 'curie':
            parameters['code.next_job'] = read_template('next_job_curie', parameters)
    else:
        parameters['code.next_job'] = ''

    # fill code part
    parameters['code.run_part'] = parameters['code.' + run_part_name] 
    parameters['log_run_part'] = str(run_part) + '.' + run_part_name
    parameters['code.part'] = read_template('part', parameters)

    # write final restart
    parameters['code.part'] += read_template('write_restart', {'restart_file': parameters['final_restart_file']})


    # Generate script from template and lammps_variable
    
    lmp_input =  read_template('core', parameters)
    # Save script
    with open(folder_path / 'lmp_input.in', 'w') as f:
        f.write(lmp_input)
    if run_on in ['zay']:
        parameters['nodes'] = config['configuration']['nodes']
        parameters['walltime'] = config['configuration']['walltime']
        if parameters['lmp_version'] == 'lmp_nequip_stress':
            if parameters['pytorch_version'] == '1.13.0':
                parameters['path_to_lammps'] = 'pair_nequip_stress_branch_pytorch1.13.0/lammps/build/lmp'
            elif parameters['pytorch_version'] == '1.10.1':
                parameters['path_to_lammps'] = 'pair_nequip_stress_branch/lammps/build/lmp_nequip_stress'
        slurm_input = read_template(f'job_{run_on}', parameters)
        with open(folder_path / 'job.slurm', 'w') as f:
            f.write(slurm_input)
    elif run_on == 'curie':
        parameters['ncores'] = config['configuration']['ncores']
        sh_input = read_template('job_curie', parameters)
        with open(folder_path / 'job.sh', 'w') as f:
            f.write(sh_input)


def write_description(folder_path, output_dict, run_description):
    """
    contains commented legacy code to run on curie, modify to make it an option
    """
    os.chdir(folder_path+"/")
    with open('config.json', 'w') as fp:
        json.dump(output_dict, fp)

    df = pd.DataFrame(run_description)
    df[["run_part_name","run_subpart_name"]] = df.output_name.str.split('.', n = 1, expand = True) 
    df.drop(columns =["output_name"], inplace = True) 
    df['run_part'] = df.groupby("run_part_name", sort=False).ngroup() + 1 # starts at 1 up to number of parts
    df['run_subpart'] = df.groupby("run_part_name", sort=False).cumcount()+1
    # df['run_exp'] = run_exp
    df.to_csv("run_description.csv", index=False)

def write_sub_restart(folder_path, prev_nb_restarts, nb_restarts):
    """
    """
    sub_restart = read_template('sub_restart_first', {'restart_number':prev_nb_restarts,'jid_job':1})
    jid_job = 2 
    for restart in range(prev_nb_restarts + 2, nb_restarts + 1):
        if restart < nb_restarts:
            sub_restart += read_template('sub_restart_core', {'restart_number':restart - 1,'jid_job':jid_job,'jid_job_previous':jid_job - 1})
        else:
            sub_restart += read_template('sub_restart_last', {'restart_number':restart - 1,'jid_job_previous':jid_job - 1})
        jid_job += 1
    with open(pathlib.Path(folder_path) / f'sub_restart{prev_nb_restarts}-{nb_restarts-1}.sh', 'w') as f:
        f.write(sub_restart)


def write_sub_restart_successive(folder_path, run_exp_dict, nb_restarts):
    """
    subdirs_list: ordered list of subdirs
    """
    # subdirs_list = sorted(run_exp_dict.values())
    jid_job = 2 
    sub_restart = """#!/bin/sh  \n"""
    for i, subdir in sorted(run_exp_dict.items()):
        # subdir = subdirs_list[i]
        sub_restart += 'cd ' + subdir + '\n'
        first_restart_of_run_exp = 1
        if i == 1:
            sub_restart += read_template('sub_restart_first', {'restart_number':0,'jid_job':1})
            first_restart_of_run_exp += 1
        for restart in range(first_restart_of_run_exp, nb_restarts + 1):
            if restart < nb_restarts or i < len(run_exp_dict):
                sub_restart += read_template('sub_restart_core', {'restart_number':restart - 1,'jid_job':jid_job,'jid_job_previous':jid_job - 1})
            else:
                sub_restart += read_template('sub_restart_last', {'restart_number':restart - 1,'jid_job_previous':jid_job - 1})
            jid_job += 1
        sub_restart += 'cd .. \n'
    run_exp_id = sorted(run_exp_dict.keys())
    with open(pathlib.Path(folder_path) / f'sub_run_exp{run_exp_id[0]}-{run_exp_id[-1]}.sh', 'w') as f:
        f.write(sub_restart)

def run_exp(constant, variable, id, previous_folder = None, next_folder = None):
    """
    input: 
    two dicts containing the constant and variable parameters
    id: int providing the id used for naming"""
    os.chdir(working_dir)

    parameters = {**constant, **variable, **dict(config['system_info'])} # merge two dicts + add system info

    run_exp = format(id, '03d') + "_" + '_'.join([key+val for key, val in variable.items()])    
    parameters['run_exp_name'] = '_'.join([key+val for key, val in variable.items()])
    parameters['run_exp'] = id

    try:
        parameters['gaillac_id'] = config['configuration']['gaillac_file_prefix'] + format(int(parameters['gaillac_glass_number']), '02d')
    except:
        pass

    logger.info('Write %s', run_exp)

    # create output dic
    output_dict = {**parameters}

    # change timestep, tdamp and pdamp in ps (metal units) from fs
    # this way, in the output dic it is stored in fs (for constency with other schemes)
    for var_name in ['timestep', 'tdamp', 'pdamp']:
        if var_name in parameters.keys():
            parameters[var_name] = str(1e-3 * float(parameters[var_name]))
    # parameters['run_serie'] = config['system_info']['run_serie']
    
    # create list of dict for run part description
    run_description = []



    simulation_time = 0 # in ps
    # write lammps input code from nve/nvt/npt bricks


    # no volume change, only equilibration sim
    i_restart = 0
    nb_restarts = int(config['configuration'][f'num.restart{i_restart}'])
    prev_nb_restarts = 0
    for i in range(nb_restarts):
        parameters['code.equilibrate'], simulation_time, run_description = nvt(f'equilibrate.{i+1}', int(parameters['temp']), int(parameters['temp']), float(parameters['time.equilibrate']), simulation_time, parameters, run_description, multiple_restarts=int(parameters['num.equilibrate.write']))

        folder = f'{run_exp}/restart{i}'
        first_run = (i == 0) and (previous_folder is None)
        first_run_exp = (i == 0)
        next_folder_restart = f'restart{i+1}' if i + 1 < nb_restarts else next_folder # handles next_folder=None in input
        previous_folder_restart = previous_folder if (i == 0) else None
        parameters['job-name'] = '-'.join(['7',str(int(parameters['run_serie'])), str(id), str(i)])
        write_lammps_input(folder, parameters, 1, 'equilibrate', first_run = first_run, next_folder = next_folder_restart, previous_folder = previous_folder_restart, first_run_exp = first_run_exp) # comment line for restarts if already ran
    if config['configuration']['variation.style'] != 'successive':
        write_sub_restart(run_exp, prev_nb_restarts, nb_restarts)

    # Uncomment the following block for restart, and comment write_lammps_input and write_sub_restart above

    # i_restart += 1
    # prev_nb_restarts = nb_restarts
    # nb_restarts = int(config['configuration'][f'num.restart{i_restart}'])
    # for i in range(prev_nb_restarts, nb_restarts):
    #     parameters['code.equilibrate'], simulation_time, run_description = nvt(f'equilibrate.{i+1}', int(parameters['temp']), int(parameters['temp']), float(parameters[f'time.equilibrate.restart{i_restart}']), simulation_time, parameters, run_description, multiple_restarts=int(parameters[f'num.equilibrate.write.restart{i_restart}']))

    #     folder = f'{run_exp}/restart{i}'

    #     new_restart = (i == prev_nb_restarts)
    #     previous_folder = f'restart{i-1}' if new_restart else None

    #     next_folder_restart = f'restart{i+1}' if i + 1 < nb_restarts else next_folder # handles next_folder=None in input
    #     previous_folder_restart = previous_folder if (i == prev_nb_restarts) else None
    #     parameters['job-name'] = '-'.join(['7',str(int(parameters['run_serie'])), str(id), str(i)])
    #     write_lammps_input(folder, parameters, 1, 'equilibrate', next_folder = next_folder_restart, previous_folder = previous_folder_restart) # comment line for restarts if already ran
    # if config['configuration']['variation.style'] != 'successive':
    #     write_sub_restart(run_exp, prev_nb_restarts, nb_restarts)

    write_description(run_exp, output_dict, run_description)



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


logger.info('Serie number #' + config['system_info']['run_serie'])

if config['configuration']['variation.style'] == 'successive':
    i = 1 #id
    run_exp_dict = {}
    for variable in variable_list:
        run_exp_dict[i] = format(i, '03d') + "_" + '_'.join([key+val for key, val in variable.items()])    
        i += 1
    
    nb_restarts = int(config['configuration']['num.restart1']) # 1 b/c deform + equilibrate made with two restarts

    write_sub_restart_successive('', run_exp_dict, nb_restarts)   

    i = 1 #id
    for variable in variable_list:
        # previous final.restart will be copied twice (once per sbatch) but it has no effect on the simulation
        # previous_folder = None
        previous_folder = f'../{run_exp_dict[i-1]}/restart{nb_restarts-1}' if i != 1 else None 
        next_folder = f'../{run_exp_dict[i+1]}/restart0' if i < len(variable_list) else None
        run_exp(constant, variable, i, previous_folder = previous_folder, next_folder = next_folder)
        i += 1        
else:
    i = 1 #id
    for variable in variable_list:
        run_exp(constant, variable, i)
        i += 1