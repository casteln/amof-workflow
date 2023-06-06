"""
Read thermo
---
Read lammps log files and merge them to an xarray dataset
"""

"""Imports"""
import pandas as pd
import numpy as np
import xarray as xr

import atomman.lammps as lmp

import ase.io

import re
from copy import deepcopy

import os

import json

import logging

import shared_functions as sf
from shared_variables import *

working_dir = os.path.dirname(os.path.realpath(__file__))
os.chdir(working_dir) # change working dir as where the python script is


logging.basicConfig(
    level=logging.INFO, #can be DEBUG, INFO, etc.
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("scripts.log"), 
        logging.StreamHandler()
    ]  
)
logger = logging.getLogger(__name__)


def read_log_file(file, dic_output, list_of_dict, base_dic):
    """
    Read log file and create run_id to tag it

    Args:
        dic_output: dic of output where keys are variable names ('thermo') and val are dicts representing output
                output keys are run_id and values are xarray dataframe containing the corresponding thermo data
        list of dict: each dict correspond to a run command in lammps for a given run_exp
        base_dict: serves as basis for every dict that list_of_dict will contain

    Returns:
        dic_thermo, list of dict: updated
    """
    results = lmp.Log()
    results.read(file)
    for j in range(len(results.simulations)):
        dic = deepcopy(base_dic)
        dic['run_subpart'] = j + 1
        run_id = '-'.join([format(dic[key], '03d') for key in ['run_serie','run_exp','run_part','run_subpart']])
        dic['run_id'] = run_id
        df = results.simulations[j]['thermo']
        df = df.set_index('Step')

        if 'Pxx' in df.columns:
            df_stress = df[['Pxx', 'Pxy', 'Pxz',  # use the symetry of the stress tensor
                            'Pxy', 'Pyy', 'Pyz',
                            'Pxz', 'Pyz', 'Pzz']]

            stress_ar = xr.DataArray(
                    np.array([c.reshape(3,3) for c in np.array(df_stress)]),
                    coords = [('Step', np.array(df_stress.index)), 
                                ('row', ['x', 'y', 'z']), 
                                ('col', ['x', 'y', 'z'])]
                    )
            stress_ar.attrs["units"] = "bars"
            dic_output['stress'][run_id] = stress_ar                    

            df = df.drop(columns= ['Pxx', 'Pyy', 'Pzz', 'Pxy', 'Pxz','Pyz'])
        ar = xr.DataArray(df, dims=("Step", "thermo_var"))
        dic_output['thermo'][run_id] = ar
        list_of_dict.append(dic)
    return dic_output, list_of_dict

def read_ase_output(folder, dic_output, list_of_dict, base_dic):
    """
    Read ase output file and create run_id to tag it

    Args:
        dic_output: dic of output where keys are variable names ('thermo') and val are dicts representing output
                output keys are run_id and values are xarray dataframe containing the corresponding thermo data
        list of dict: each dict correspond to a run command in lammps for a given run_exp
        base_dict: serves as basis for every dict that list_of_dict will contain

    Returns:
        dic_output, list of dict: updated
    """
    dic = deepcopy(base_dic)
    run_id = '-'.join([format(dic[key], '03d') for key in ['run_serie','run_exp','run_part','run_subpart']])
    dic['run_id'] = run_id
    
    from_ev_per_ang_3_to_bars = 1.602176634 * 1e6 # convert from ase units to bars (lmp units)

    for file in folder.glob('*.log'):
        with open(file) as f:
            first_line = f.readline()
            second_line = f.readline()
        line_list = second_line.split()
        # first two values are treated separately as they don't have units attached
        names = line_list[:4:2] + line_list[4:][::3]
        rename_var = {'Frame': 'Step', 'SimulationTime':'Time', 
            'Temperature':'Temp', 'PotentialEnergy':'PotEng', 'TotalEnergy': 'TotEng'}
        names = [s.strip(':') for s in names]
        names = [rename_var[s] if s in rename_var.keys() else s for s in names]
        usecol = [1,3] + list(np.arange(5, len(line_list), 3))
        log = pd.read_table(file, skiprows=1, delim_whitespace=True, 
            usecols=usecol, 
            names=names)       

        log['Step'] = log['Step'] + dic['first_step']
        log['Time'] = log['Time'] + dic['first_step'] * float(dic['timestep'])

        # convert units from ase to lammps metal units
        log['Time'] = log['Time'] / 1000 # from fs to ps
        log = log.set_index('Step')

        if 'xx' in log.columns:
            # multiply by -1 to have the same convention than lammps
            df_stress = (-1) * from_ev_per_ang_3_to_bars * log[['xx', 'xy', 'xz',  # use the symetry of the stress tensor
                            'xy', 'yy', 'yz',
                            'xz', 'yz', 'zz']]
            stress_ar = xr.DataArray(
                    np.array([c.reshape(3,3) for c in np.array(df_stress)]),
                    coords = [('Step', np.array(df_stress.index)), 
                                ('row', ['x', 'y', 'z']), 
                                ('col', ['x', 'y', 'z'])]
                    )
            stress_ar.attrs["units"] = "bars"
            dic_output['stress'][run_id] = stress_ar                    

            log['Press'] = df_stress[['xx', 'yy', 'zz']].mean(axis=1)   

            log = log.drop(columns= ['xx', 'yy', 'zz', 'xy', 'xz','yz'])

    thermo = log


    for file in folder.parent.glob('*.xyz'):
        if 'Press' in thermo.columns and 'Volume' in thermo.columns:
            # Load only first frame
            traj_index = slice(1)      
        else:
            endpoint = int(dic['freq.dump']) == int(dic['freq.thermo']) # special case, need to include last frame
            traj_index = sf.get_traj_index(log.index.min(), log.index.max(), 
                int(dic['freq.dump']), endpoint=endpoint)      
        traj = ase.io.read(file, index=traj_index)

    if 'Press' not in thermo.columns and 'stress' in traj[0].info:

        from_ev_per_ang_3_to_bars = 1.602176634 * 1e6 # convert from ase units to bars (lmp units)
        # multiply by -1 to have the same convention than lammps
        stress = (-1) * from_ev_per_ang_3_to_bars * np.array([atom.info['stress'] for atom in traj])
        step_ar = sf.get_interval(log.index.min(), log.index.max(), 
                int(dic['freq.dump']), endpoint=endpoint)
        stress_ar = xr.DataArray(stress,
                coords = [('Step', step_ar), 
                            ('row', ['x', 'y', 'z']), 
                            ('col', ['x', 'y', 'z'])]
                        )
        stress_ar.attrs["units"] = "bars"
        dic_output['stress'][run_id] = stress_ar

        pressure = np.trace(stress, axis1 = 1, axis2 = 2) / 3

        # create "pressure" column by merging to cover the case where the steps are sampled differently
        pressure_df = pd.DataFrame(pressure, index=step_ar, columns = ['Press'])
        pressure_df.index.name = 'Step'
        thermo = pd.merge(thermo, pressure_df, how="left", left_index=True, right_index=True)

    if 'Volume' not in thermo.columns:
        volume = [atom.cell.volume for atom in traj]

        # create "Volume" column by merging to cover the case where the steps are sampled differently
        volume_df = pd.DataFrame(volume, index=step_ar, columns = ['Volume'])
        volume_df.index.name = 'Step'
        thermo = pd.merge(thermo, volume_df, how="left", left_index=True, right_index=True)

    atomic_mass_constant =  1.66053906660 # in 1E-27 kg, 
    mass = np.sum(traj[0].get_masses())
    thermo['Density'] = atomic_mass_constant * mass / thermo['Volume']


    thermo_ar = xr.DataArray(thermo, dims=("Step", "thermo_var"))
    dic_output['thermo'][run_id] = thermo_ar

    list_of_dict.append(dic)
    return dic_output, list_of_dict

def read_run_exp(folder_path, dic_output, list_config, base_dic, log_prefix = "log.run"):
    """
    Read folder corresponding to run_exp: read each log and run_description/config files

    Args:
        dic_thermo: keys are run_id and values are xarray dataframe containing the corresponding thermo data
        list_config: list of pandas dataframes containg 'run_description' and 'config.json' content
        base_dict: serves as basis for every dict that list_of_dict will contain

    Returns:
        dic_output, list_config: updated
    """
    config_dic = json.load(open(folder_path / "config.json"))
    # handle "multiple" keyword in runserie_description
    # keep the value from config.json
    columns_to_check = ['system', 'model'] # just add a name here to handle a new column
    for column_name in columns_to_check:
        if column_name in base_dic and base_dic[column_name] == 'multiple' and column_name in config_dic:
            base_dic.pop(column_name)

    base_dic = {**config_dic, **base_dic} # base_dic takes precedence
    run_description = pd.read_csv(folder_path / "run_description.csv")
    list_of_dict = []
    if base_dic['implementation'] == 'ase':
        for folder in sorted(folder_path.glob('restart*')):
            if len(list(folder.rglob('*.xyz*'))) != 0:
                dic = deepcopy(base_dic)
                dic['run_part'] = 1
                dic['run_subpart'] = int(re.search('restart(.*)', folder.name).group(1)) + 1
                dic['first_step'] = int(run_description[run_description.run_subpart==dic['run_subpart']].first_step)
                dic_output, list_of_dict = read_ase_output(folder, dic_output, list_of_dict, dic)
    elif base_dic['implementation'] == 'lammps':
        for file in sorted(folder_path.glob(f'{log_prefix}.*.*')):
            dic = deepcopy(base_dic)
            run_part = file.name.split(log_prefix)[-1].split('.')
            dic['run_part'] = int(run_part[1])
            dic['run_part_name'] = run_part[2]  
            dic_output, list_of_dict = read_log_file(file, dic_output, list_of_dict, dic)
    if list_of_dict != []: # doesn't add this run_exp if no file found
        config_df = pd.DataFrame(list_of_dict)
        config_df = pd.merge(config_df, run_description, how = 'left', on=['run_part', 'run_subpart'], suffixes=('', '_y'))
        config_df.drop(config_df.filter(regex='_y$').columns.tolist(),axis=1, inplace=True) # remove duplicate columns, keep the version from the right-hand df
        # config_df = config_df.convert_dtypes() # automatically detect dtypes
        list_config.append(config_df)
    return dic_output, list_config

def read_run_serie(folder_path, dic_output, list_config, base_dic):
    """
    Read folder corresponding to run_serie: read each run_exp

    Args:
        dic_output: keys are run_id and values are xarray dataframe containing the corresponding thermo data
        list_config: list of pandas dataframes containg 'config.json' content
        run_exp: 3 digit string
        base_dict: serves as basis for every dict that list_of_dict will contain

    Returns:
        dic_output, list_config: updated
    """
    for folder in sorted(folder_path.glob("*")):
        run_exp = folder.name[0:3]
        if run_exp.isdigit() and len(list(folder.rglob('*.xyz'))) != 0:
            logger.info("Read run_exp %s", run_exp)
            dic = deepcopy(base_dic)
            dic['run_exp'] = int(run_exp) # may remove former run_exp content  
            dic_output, list_config = read_run_exp(folder, dic_output, list_config, dic)
    return dic_output, list_config

def main(path_to_traj, path_to_dataset, overwrite = False):
    """
    Read folder containging every run_serie desribed in runserie_description,
    create xarray dataset and write it

    Args:
        path_to_traj: path to folder containing run_series folders and runserie_description
        path_to_dataset: where to write netCFD files

    Returns:
        data: xarray dataset
    """
    logger.info("Read input path %s", path_to_traj)    
    runserie_description = pd.read_csv(path_to_traj / f"{runcategory}-runserie_description.csv")
    runserie_description = runserie_description.set_index('run_serie')

    for index, row in runserie_description.iterrows():
        run_serie = format(index, '03d')
        if row['include_in_dataset']:
            list_config = []
            dic_output = {'thermo': {}, 'stress': {}} # dic of dic output
            # dic_thermo = {}
            filename = "run_serie_" + run_serie 
            if overwrite or not (path_to_dataset / 'thermo' / f'{filename}.comp.nc').exists():
                logger.info("Read run_serie %s", run_serie)    
                base_dic = row.to_dict()
                base_dic['run_serie'] = int(run_serie)
                dic_output, list_config = read_run_serie(path_to_traj / run_serie, dic_output, list_config, base_dic)

                logger.info("Converting config to xarray")    
                df = pd.concat(list_config)
                df = df.set_index('run_id')
                config_ds = df.to_xarray()

                logger.info("Converting thermo to xarray")   
                dic_xa = {}
                for var_name, output in dic_output.items():
                    if output != {}:
                        xa = xr.Dataset(output)
                        xa = xa.to_array("run_id", "run_id")
                        dic_xa[var_name] = xa
                output_ds =  xr.Dataset(dic_xa) # one large data_array containing thermo variables
                # xa = xa.to_dataset(dim='thermo') # separate thermo variables

                logger.info("Merging xarrays")    
                data = xr.merge([output_ds, config_ds])

                sf.write_dataset(data, path_to_dataset / 'thermo' / f'{filename}.comp.nc', compress=True)
            if overwrite or not (path_to_dataset / 'thermo-flat' / f'{filename}.comp.nc').exists():
                logger.info("Flatten run_serie %s", run_serie)    
                data = xr.open_dataset(path_to_dataset / 'thermo' / f'{filename}.comp.nc')
                data = sf.flatten_dataset(data)
                sf.write_dataset(data, path_to_dataset / 'thermo-flat' / f'{filename}.comp.nc', compress=True)
        
        else:
            logger.info("Ignore run_serie %s", run_serie)

    preprocess_func = (lambda ds: ds.drop_vars('thermo').drop_dims(['thermo_var', 'Step']).reindex(run_id = np.sort(ds.run_id)))
    config_da = xr.open_mfdataset(str(path_to_dataset) + '/thermo/run_serie_*.comp.nc', parallel=True, concat_dim = 'run_id', data_vars = 'minimal', preprocess = preprocess_func)
    config_da.load()
    sf.write_dataset(config_da, path_to_dataset / f'config.comp.nc', compress=True)
    config_da_flat = sf.flatten_dataset(config_da)
    sf.write_dataset(config_da_flat, path_to_dataset / f'config-flat.comp.nc', compress=True)
        


if __name__ == "__main__":
    # Makes sure the output directories exists
    for path in [path_to_dataset / folder for folder in ['thermo', 'thermo-flat']]:
        path.mkdir(parents=True, exist_ok=True) 

    main(path_to_traj, path_to_dataset, overwrite = False)