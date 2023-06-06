import numpy as np

import xarray as xr

import logging
import os

# force numpy to use one thread to avoid msd launching on all cores (doesn't work if large memory use (>size ram))
os.environ["OMP_NUM_THREADS"] = "1" # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = "1" # export OPENBLAS_NUM_THREADS=4 
os.environ["MKL_NUM_THREADS"] = "1" # export MKL_NUM_THREADS=6
os.environ["NUMEXPR_NUM_THREADS"] = "1" # export NUMEXPR_NUM_THREADS=6
os.environ["OPENBLAS_MAIN_FREE"] = "1"

import gc
import pandas as pd
import sys

import ase.io

import amof.rdf as srdf
import amof.msd as smsd
import amof.trajectory as straj
import amof.elastic as sela
import amof.pore as spore
import amof.bad as sbad
import amof.cn as scn
import amof.coordination.reduce as sred
import amof.ring as sring

from shared_variables import *
import shared_functions as sf


working_dir = os.path.dirname(os.path.realpath(__file__))
os.chdir(working_dir)  # change working dir as where the python script is

logging.basicConfig(
    level=logging.INFO,  # can be DEBUG, INFO, etc.
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("scripts.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def get_cell(x, freq='dump', sampling_rate=1):
    """return cell array
    Args: 
        x:  run_id row of dataset
        freq: can be either 'dump' or 'thermo' and specify the spacing 
        sampling_rate: int, multiply freq 
    Returns:
        cell: xarray where each line is a cell vector"""
    ds = x.sel(thermo_var=['Cella', 'Cellb', 'Cellc',
               'CellAlpha', 'CellBeta', 'CellGamma'])
    ds = ds.sel(Step=sf.get_interval(int(x.first_step), int(
        x.last_step), int(x[f'freq.{freq}']), sampling_rate = sampling_rate))
    cell = ds.thermo
    # shorten cell
    cell = cell.where(cell.notnull(), drop=True)
    try:
        cell = cell.squeeze('run_id')
    except KeyError:
        logger.debug("run_id is not a dim of x.cell")
    # except IndexError:
    #     logger.debug("cell is empty")
    # ensure order of cell vector
    cell = cell.reindex(
        thermo_var=['Cella', 'Cellb', 'Cellc', 'CellAlpha', 'CellBeta', 'CellGamma'])
    return cell

def format_id(id):
    """return 3 digit string of id, xarray with one element"""
    return format(int(id.item()), '03d')

def construct_traj(x, path_to_traj, sampling_rate=1):
    """construct trajectory from cell info in dataset and xyz file
    Args:
        x: run_id row x of dataset
        sampling_rate: int, multiply freq 
    """
    logger.info('Construct trajectory for run_id %s', x.run_id.item())
    traj_index = sf.get_traj_index(int(x.first_step), int(
            x.last_step), int(x['freq.dump']), sampling_rate)    
    if x['implementation'] == 'ase':
        for filename in path_to_traj.glob(f'{format_id(x.run_serie)}/{format_id(x.run_exp)}*/*.xyz'):
            path_to_xyz = filename  # only one value in glob
        traj = ase.io.read(path_to_xyz, traj_index, 'extxyz')   
    elif x['implementation'] == 'lammps':
        cell = get_cell(x, sampling_rate=sampling_rate)
        for filename in path_to_traj.glob(f'{format_id(x.run_serie)}/{format_id(x.run_exp)}*/ZIF4.xyz'):
            path_to_xyz = filename  # only one value in glob
        traj = straj.read_lammps_traj(path_to_xyz, index=traj_index, cell=cell)
    return traj


def compute_properties_by_run_id(x, path_to_traj, path_to_output, compute_property, overwrite):
    """return cell array of run_id row x of dataset"""
    run_id_str = x.run_id.item()
    filename = f'run_id_{run_id_str}'
    traj, red_traj = None, None
    default_sampling_rate = 2 # 1 <-> 5 fs per frame
    reduction_sampling_rate = 20  # 200 <-> 100 ps per frame
    if compute_property['rdf'] and (overwrite or not (path_to_output['rdf'] / f'{filename}.rdf').exists()):
        if traj is None:
            traj = construct_traj(x, path_to_traj, sampling_rate=default_sampling_rate)
        logger.info("Compute rdf of %s", filename)
        rdf = srdf.Rdf.from_trajectory(traj, rmax=12)
        rdf.write_to_file(path_to_output['rdf'] / filename)
    if compute_property['cn'] and (overwrite or not (path_to_output['cn'] / f'{filename}.cn').exists()):
        if traj is None:
            traj = construct_traj(x, path_to_traj, sampling_rate=default_sampling_rate)
        logger.info("Compute coordination number of %s", filename)
        cn = scn.CoordinationNumber.from_trajectory(
            traj, {'Zn-N': 2.5}, delta_Step=int(x['freq.dump']), first_frame = int(x.first_step), parallel=True)
        cn.write_to_file(path_to_output['cn'] / filename)
    if compute_property['bad'] and (overwrite or not (path_to_output['bad'] / f'{filename}.bad').exists()):
        if traj is None:
            traj = construct_traj(x, path_to_traj, sampling_rate=default_sampling_rate)
        logger.info("Compute bad of %s", filename)
        bad = sbad.Bad.from_trajectory(traj, {'Zn-N': 2.5}, parallel=True)
        bad.write_to_file(path_to_output['bad'] / filename)       
    if compute_property['bad_by_cn'] and (overwrite or not (path_to_output['bad_by_cn'] / f'{filename}.bad').exists()):
        if traj is None:
            traj = construct_traj(x, path_to_traj)
        logger.info("Compute bad of %s", filename)
        bad = sbad.BadByCn.from_trajectory(traj, {'Zn-N': 2.5}, parallel=True)
        bad.write_to_file(path_to_output['bad_by_cn'] / filename)           
    if compute_property['msd_by_run_id'] and (overwrite or not (path_to_output['msd'] / f'{filename}.msd').exists()):
        if traj is None:
            traj = construct_traj(x, path_to_traj, sampling_rate=default_sampling_rate)
        logger.info("Compute msd of %s", filename)
        # msd = smsd.Msd.from_trajectory(traj, delta_Step=int(x['freq.dump']), first_frame = int(x.first_step))
        msd = smsd.WindowMsd.from_trajectory(traj, 
            delta_time=round(float(x.timestep.item()) * float(x['freq.dump'].item()) * default_sampling_rate),
            max_time=int(20e3), 
            timestep = round(float(x.timestep.item()) * float(x['freq.dump'].item()) * default_sampling_rate), parallel = True, 
            unwrap = True)
        msd.write_to_file(path_to_output['msd'] / filename)
    if compute_property['pore'] and (overwrite or not (path_to_output['pore'] / f'{filename}.pore').exists()):
        sampling_rate = 200  # 200 <-> 1 ps per frame
        if traj is None or sampling_rate % default_sampling_rate != 0:
            traj_sampled = construct_traj(
                x, path_to_traj, sampling_rate=sampling_rate)
        else:
            traj_sampled = traj[::sampling_rate // default_sampling_rate]  
        logger.info("Pore analysis of %s", filename)
        try:
            pore = spore.Pore.from_trajectory(
                traj_sampled, delta_Step=sampling_rate * int(x['freq.dump']), first_frame = int(x.first_step), parallel=True)
        except OSError:
            logger.warning("Unsufficient memory, reduce number of cores to 8")
            pore = spore.Pore.from_trajectory(
                traj_sampled, delta_Step=sampling_rate * int(x['freq.dump']), first_frame = int(x.first_step), parallel=8)
        pore.write_to_file(path_to_output['pore'] / filename)
    if compute_property['reduction'] and (overwrite or not (path_to_output['reduction'] / f'{filename}.xyz').exists()):
        # if x['system'].item() == 'abinitio_glass':
        if traj is None or reduction_sampling_rate % default_sampling_rate != 0:
            traj_sampled = construct_traj(
                x, path_to_traj, sampling_rate=reduction_sampling_rate)
        else:
            traj_sampled = traj[::reduction_sampling_rate // default_sampling_rate]  
        # if traj is None:
        #     traj = construct_traj(x, path_to_traj)
        logger.info("Reduce %s", filename)
        red_traj = sred.reduce_trajectory(traj_sampled, x['mof'].item(), 
            delta_Step=reduction_sampling_rate*int(x['freq.dump']), first_frame = int(x.first_step), parallel=True)
        red_traj.write_to_file(path_to_output['reduction'] / filename)
    if (compute_property['ring'] and 
            (overwrite or not (path_to_output['ring'] / f'{filename}.ring').exists()) and
            (path_to_output['reduction'] / f'{filename}.xyz').exists()):
        sampling_rate = 200  # 200 <-> 1 ps per frame
        if red_traj is None:
            red_traj = straj.ReducedTrajectory.from_file(path_to_output['reduction'] / filename, sampling=sampling_rate // reduction_sampling_rate)
        else:
            red_traj.sample(sampling_rate // reduction_sampling_rate)
        if red_traj.trajectory != [] and red_traj.trajectory[0].get_global_number_of_atoms() > 350: #333 supercell is 384
            red_traj.sample(5)
        logger.info("Compute ring statistics for %s", filename)
        ring = sring.Ring.from_reduced_trajectory(red_traj, 
            max_search_depth=32,
            discard_if_potentially_undiscovered_rings=False,
            parallel = True)
        ring.write_to_file(path_to_output['ring'] / filename)

    
    # garbage collector call
    try:
        del traj, traj_sampled, red_traj
        gc.collect()
    except:
        pass
    return xr.DataArray()  # .map() method requires this func to return a dataarray


def compute_properties_by_run_exp(x, path_to_traj, path_to_output, compute_property, overwrite):
    """return cell array of run_id row x of dataset"""
    run_id_str = x.run_id.head(1).item()
    filename = f'run_exp_{run_id_str}'
    traj = None
    if compute_property['msd'] and (overwrite or not (path_to_output['msd'] / f'{filename}.msd').exists()):
        sampling_rate = 10  # temp fix while msd not changed to handle ram
        if traj is None:
            # x = sf.flatten(x)
            traj = construct_traj(x, path_to_traj, sampling_rate=sampling_rate)
        logger.info("Compute msd of %s", filename)
        # msd = smsd.Msd.from_trajectory(traj, delta_Step=sampling_rate * int(x['freq.dump']), first_frame = int(x.first_step), parallel=2) # 2 to reduce mem use
        msd = smsd.WindowMsd.from_trajectory(traj, 
            delta_time=round(2 * float(x.timestep.item()) * float(x['freq.dump'].item()) * sampling_rate),
            max_time=int(30e3), 
            timestep = round(float(x.timestep.item()) * float(x['freq.dump'].item()) * sampling_rate), parallel = True, 
            unwrap = True)
        msd.write_to_file(path_to_output['msd'] / filename)
    if compute_property['elastic'] and x.compute_elastic_constants.head(1).item() and (overwrite or not (path_to_output['elastic'] / f'{filename}.elastic').exists()):
        # can change to thermo to compute with a higher dump frequency (gives same Cmat)
        cell = get_cell(x, freq='dump')
        logger.info("Compute elastic constant of %s", filename)
        try:
            elastic_config = pd.read_csv("elastic_config.csv").set_index("run_id")
            first_frame = elastic_config.compute_from_frame[run_id_str]
        except:
            first_frame = 2000000
        logger.info("Remove first %s ns", first_frame * float(x.timestep.item()) * 1e-6)
        cell = cell.sel(Step=slice(first_frame, None))
        ela = sela.ElasticConstant.from_cell(
            cell, int(x.temp.item()), step=cell.Step, final_value=False)
        ela.write(path_to_output['elastic'] / filename)
    return xr.DataArray()  # .map() method requires this func to return a dataarray

def compute_properties(da, path_to_traj, path_to_output, compute_property, overwrite, by='run_id'):
    """
    Compute rdf for every run_id contain in da and write it

    Args:
        da: xarray dataset
        path_to_output: dictionary to know where to write amof rdf files
        compute_property: dictionary where key is prop and value is boolean indicating whether the property should be computed
        by: str, 'run_id' or 'run_exp'

    Returns:
        rdf_array: xarray dataset containg rdf data
    """
    if by == 'run_id':
        # da['run_id'] = da['run_id'].astype('str')
        da.groupby("run_id").map(compute_properties_by_run_id, args=[
            path_to_traj, path_to_output, compute_property, overwrite])
    elif by == 'run_exp':
        da.groupby("run_id").map(compute_properties_by_run_exp, args=[
            path_to_traj, path_to_output, compute_property, overwrite])


def main(compute_property, overwrite=False):
    dataset_name = 'thermo.comp.nc'

    logger.info("Start compute properties")

    for filename in path_to_dataset.glob('thermo/run_serie_*[0-9].comp.nc'):
        da = xr.open_dataset(filename)
        compute_properties(da, path_to_traj, path_to_output, compute_property, overwrite, by = 'run_id')
    for filename in path_to_dataset.glob('thermo-flat/run_serie_*[0-9].comp.nc'):
        da = xr.open_dataset(filename)
        compute_properties(da, path_to_traj, path_to_output,
                           compute_property, overwrite, by='run_exp')
    path_to_output


if __name__ == "__main__":
    compute_property = {'rdf': True, 'msd': True,
                        'msd_by_run_id': True, 'elastic': True, 'pore': True, 'cn': True, 'bad': True,  'bad_by_cn':True,
                        'reduction': True, 'ring': True}
    try:
        if sys.argv[1] in properties_dict.keys():
            compute_property = properties_dict[sys.argv[1]]
    except:
        pass     
    main(compute_property, overwrite=False)
