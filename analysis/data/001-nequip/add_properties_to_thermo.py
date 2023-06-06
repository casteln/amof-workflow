import xarray as xr

import logging
import os
import numpy as np
import re
import sys

import amof.rdf as srdf
import amof.msd as smsd
import amof.elastic as sela
import amof.pore as spore
import amof.bad as sbad
import amof.cn as scn
import amof.ring as sring
import amof.trajectory as straj

import shared_functions as sf
from shared_variables import *


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

def slice_elastic_by_run_id(x, Cmat):
    return Cmat.sel(Step=slice(int(x.first_step), int(x.last_step)))


def read_elastic_by_run_exp(x, path_to_output, add_property, slice_by_run_id=False):
    """return cell array of run_id row x of dataset"""
    run_id_str = x.run_serie_exp.head(1).item()
    filename = f'run_exp_{run_id_str}'
    Cmat = sela.ElasticConstant().Cmat  # create empty Cmat to start with
    for file in path_to_output['elastic'].glob(f'{filename}.elastic'):
        elastic = sela.ElasticConstant.from_file(file)
        Cmat = elastic.Cmat.elastic  # load data array
    if slice_by_run_id:
        elastic = x.groupby('run_id_str').map(
            slice_elastic_by_run_id, args=[Cmat])
    else:
        elastic = Cmat
    return elastic

def reorder_dataset_by_run_id(ds):
    """
    Remove run_serie_exp as dimension so that only run_id_str remains
    """
    ds['run_id_var'] = ds.run_id_str
    ds = (
        ds
        .stack(run_id=['run_serie_exp', 'run_id_str'])
        .set_index(run_id='run_id_var')
        .dropna('run_id', how='all')
    )
    return ds

def write_with_config(ds, name, flatten = False):
    """
    write with config info, add suffix to name
    """
    config_da = xr.open_dataset(path_to_dataset / 'config.comp.nc')
    ds_merged = xr.merge([config_da, ds], join = 'left')
    sf.write_dataset(ds_merged, path_to_dataset / f'{name}.comp.nc', compress=True)
    if flatten:
        ds_merged = sf.flatten_dataset(ds_merged)
        sf.write_dataset(ds_merged, path_to_dataset / f'{name}-flat.comp.nc', compress=True)




def read_properties(run_serie, add_property, write_properties_to_file=False):
    """
    Read every prop file

    Args:
        run_serie: string indicating run_serie, can also be "*" to include every run_serie
        add_property: dictionary where key is prop and value is boolean indicating whether the property should be computed
        write_properties_to_file: write datasets for every prop containing every run_serie

    Return:
        dataset_to_merge: list of dataset of every prop
    """
    dataset_to_merge = []

    if add_property['rdf']:
        dic_rdf = {}
        for file in path_to_output['rdf'].glob(f'run_id_{run_serie}-*.rdf'):
            rdf = srdf.Rdf.from_file(file)
            run_id = file.stem.split('run_id_')[-1]
            dic_rdf[run_id] = xr.DataArray(
                rdf.data.set_index('r'), dims=['r', 'atom_pair'])

        rdf_array = xr.Dataset(dic_rdf)
        rdf_array = rdf_array.to_array("run_id", "run_id")
        # one large data_array containing rdf variables
        rdf_ds = xr.Dataset({'rdf': rdf_array})
        rdf_ds = rdf_ds.reindex(atom_pair = np.sort(rdf_ds.atom_pair.astype('str')))
        if write_properties_to_file:
            write_with_config(rdf_ds, "rdf")
        dataset_to_merge.append(rdf_ds)

    if add_property['bad']:
        dic_bad = {}
        for file in path_to_output['bad'].glob(f'run_id_{run_serie}-*.bad'):
            bad = sbad.Bad.from_file(file)
            run_id = file.stem.split('run_id_')[-1]
            dic_bad[run_id] = xr.DataArray(
                bad.data.set_index('theta'), dims=['theta', 'atom_triple'])

        bad_array = xr.Dataset(dic_bad)
        bad_array = bad_array.to_array("run_id", "run_id")
        # one large data_array containing rdf variables
        bad_ds = xr.Dataset({'bad': bad_array})
        bad_ds = bad_ds.reindex(atom_triple = np.sort(bad_ds.atom_triple.astype('str')))
        if write_properties_to_file:
            write_with_config(bad_ds, "bad")
        dataset_to_merge.append(bad_ds)

    if add_property['bad_by_cn']:
        prop = 'bad_by_cn'
        file_ext = 'bad'
        dic_data = {}
        for file in path_to_output[prop].glob(f'run_id_{run_serie}-*.{file_ext}'):
            bad = sbad.BadByCn.from_file(file)
            run_id = file.stem.split('run_id_')[-1]
            data =  bad.data
            dic_data[run_id] = data.bad
        if dic_data != {}:
            data_array = xr.Dataset(dic_data)
            data_array = data_array.to_array("run_id", "run_id")
            data_ds = xr.Dataset({prop: data_array})
            # data_ds['ring_var'] = data_ds.ring_var.astype('str')
            if write_properties_to_file:
                write_with_config(data_ds, prop)
            dataset_to_merge.append(data_ds)

    if add_property['msd_by_run_id']:
        dic_msd = {}
        for file in path_to_output['msd'].glob(f'run_id_{run_serie}-*.msd'):
            msd = smsd.WindowMsd.from_file(file)
            run_id = file.stem.split('run_id_')[-1]
            dic_msd[run_id] = xr.DataArray(msd.data.set_index(
                'Time'), dims=['Time', 'atom'])

        msd_array = xr.Dataset(dic_msd)
        msd_array = msd_array.to_array("run_id", "run_id")
        # one large data_array containing msd variables
        msd_ds = xr.Dataset({'msd': msd_array})
        msd_ds = msd_ds.reindex(atom = np.sort(msd_ds.atom.astype('str')))
        if write_properties_to_file:
            write_with_config(msd_ds, "msd_by_run_id")
        dataset_to_merge.append(msd_ds)

    if add_property['pore']:
        dic_pore = {}
        for file in path_to_output['pore'].glob(f'run_id_{run_serie}-*.pore'):
            pore = spore.Pore.from_file(file)
            run_id = file.stem.split('run_id_')[-1]
            dic_pore[run_id] = xr.DataArray(
                pore.data.set_index('Step'), dims=['Step', 'pore_var'])

        pore_array = xr.Dataset(dic_pore)
        pore_array = pore_array.to_array("run_id", "run_id")
        # one large data_array containing rdf variables
        pore_ds = xr.Dataset({'pore': pore_array})
        pore_ds = pore_ds.reindex(pore_var = np.sort(pore_ds.pore_var.astype('str')))
        if write_properties_to_file:
            write_with_config(pore_ds, "pore", flatten=True)
        dataset_to_merge.append(pore_ds)


    if add_property['cn']:
        dic_cn = {}
        for file in path_to_output['cn'].glob(f'run_id_{run_serie}-*.cn'):
            cn = scn.CoordinationNumber.from_file(file)
            run_id = file.stem.split('run_id_')[-1]
            dic_cn[run_id] = xr.DataArray(
                cn.data.set_index('Step'), dims=['Step', 'atom_pair'])
        cn_array = xr.Dataset(dic_cn)
        cn_array = cn_array.to_array("run_id", "run_id")
        # one large data_array containing rdf variables
        cn_ds = xr.Dataset({'coordination_number': cn_array})
        cn_ds = cn_ds.reindex(atom_pair = np.sort(cn_ds.atom_pair.astype('str')))
        if write_properties_to_file:
            write_with_config(cn_ds, "cn", flatten=True)
        dataset_to_merge.append(cn_ds)

    if add_property['ring']:
        prop = 'ring'
        dic_data = {}
        for file in path_to_output[prop].glob(f'run_id_{run_serie}-*.{prop}'):
            ring = sring.Ring.from_file(file)
            run_id = file.stem.split('run_id_')[-1]
            if np.sum(list(dict(ring.data.sizes).values())) != 0:
                data =  ring.data
                dic_data[run_id] = data.ring
        if dic_data != {}:
            data_array = xr.Dataset(dic_data)
            data_array = data_array.to_array("run_id", "run_id")
            data_ds = xr.Dataset({prop: data_array})
            data_ds['ring_var'] = data_ds.ring_var.astype('str')
            if write_properties_to_file:
                write_with_config(data_ds, prop, flatten=True)
            dataset_to_merge.append(data_ds)


    # read all run_exp files without a dataset as input
    if write_properties_to_file:
        if add_property['msd']:
            dic_msd = {}
            for file in path_to_output['msd'].glob(f'run_exp_{run_serie}-*.msd'):
                msd = smsd.WindowMsd.from_file(file)
                run_id = file.stem.split('run_exp_')[-1]
                dic_msd[run_id] = xr.DataArray(msd.data.set_index(
                    'Time'), dims=['Time', 'atom'])
            msd_array = xr.Dataset(dic_msd)
            msd_array = msd_array.to_array("run_id", "run_id")
            msd_ds = xr.Dataset({'msd': msd_array})
            msd_ds = msd_ds.reindex(atom = np.sort(msd_ds.atom.astype('str')))
            sf.write_dataset(msd_ds, path_to_dataset /
                             "msd.comp.nc", compress=True)
        if add_property['elastic']:
            list_elastic = []
            for file in path_to_output['elastic'].glob(f'run_exp_{run_serie}-*.elastic'):
                elastic = sela.ElasticConstant.from_file(file)
                run_id = file.stem.split('run_exp_')[-1]
                Cmat =  elastic.Cmat.elastic
                Cmat['run_id'] = run_id
                list_elastic.append(Cmat)             
            elastic_array = xr.concat(list_elastic, dim = 'run_id')
            elastic_ds = xr.Dataset({'elastic': elastic_array})
            sf.write_dataset(elastic_ds, path_to_dataset /
                             "elastic.comp.nc", compress=True)
    return dataset_to_merge

def read_properties_by_run_exp(da, run_serie, add_property):
    """
    Read every prop file for run_exp using a flatten dataframe as input

    Args:
        da: xarray dataset containing the run_exp to load
        run_serie: string indicating run_serie, can also be "*" to include every run_serie
        add_property: dictionary where key is prop and value is boolean indicating whether the property should be computed

    Return:
        dataset_to_merge: list of dataset of every prop
    """
    dataset_to_merge = []

    if add_property['msd']:
        dic_msd = {}
        for file in path_to_output['msd'].glob(f'run_exp_{run_serie}-*.msd'):
            msd = smsd.WindowMsd.from_file(file)
            run_id = file.stem.split('run_exp_')[-1]
            dic_msd[run_id] = xr.DataArray(msd.data.set_index(
                'Time'), dims=['Time', 'atom'])
        msd_array = xr.Dataset(dic_msd)
        msd_array = msd_array.to_array("run_id", "run_id")
        msd_ds = xr.Dataset({'msd': msd_array})
        msd_ds = msd_ds.reindex(atom = np.sort(msd_ds.atom.astype('str')))
        dataset_to_merge.append(msd_ds)

    return dataset_to_merge


def add_properties_to_dataset(da, run_serie, add_property, write=False, write_flatten=False):
    """
    read every prop contained in path_to_output and add it to merged dataset

    Args:
        da: xarray dataset containing the run_exp to load
        add_property: dictionary where key is prop and value is boolean indicating whether the property should be computed
        run_serie: string indicating run_serie, can also be "*" to include every run_serie

    Returns:
    """
    if write:
        logger.info("Write run_serie %s", run_serie)
        dataset_to_merge = [da]
        dataset_to_merge += read_properties(run_serie, add_property, write_properties_to_file=False)
        # merge with left to ignore rdf files that don't correspond to run_id comprised in thermo
        data = xr.merge(dataset_to_merge, join='left') # takes a lot of RAM
        sf.write_dataset(data, path_to_dataset / 'thermo-prop' / f'run_serie_{run_serie}.comp.nc', compress=True)

    if write_flatten:
        logger.info("Write flatten run_serie %s", run_serie)
        if not write:
            data = xr.open_dataset(path_to_dataset / 'thermo-prop' / f'run_serie_{run_serie}.comp.nc')
        # naive flatten: doesn't take benefit of previous flatten
        da_f = sf.flatten_dataset(data)
        dataset_to_merge = [da_f]
        dataset_to_merge += read_properties_by_run_exp(da_f, run_serie, add_property)
        data = xr.merge(dataset_to_merge, join='left') # takes a lot of RAM
        sf.write_dataset(data, path_to_dataset / 'thermo-prop-flat' / f'run_serie_{run_serie}.comp.nc', compress=True)

def main(add_property, write_run_serie_dataset, write_properties_dataset, overwrite = False):
    """
    read every prop contained in path_to_output and add it to merged dataset

    Args:
        add_property: dictionary where key is prop and value is boolean indicating whether the property should be computed

    Returns:
    """
    logger.info("Start %s", __file__)

    if write_run_serie_dataset:
        logger.info("Write run serie datasets")
        for filename in path_to_dataset.glob('thermo/run_serie_*[0-9].comp.nc'):
            run_serie = re.search(
                'run_serie_(\d{3}).comp.nc', filename.name).group(1)
            write = overwrite or not (path_to_dataset / 'thermo-prop' / f'run_serie_{run_serie}.comp.nc').exists()
            write_flatten = overwrite or not (path_to_dataset / 'thermo-prop-flat' / f'run_serie_{run_serie}.comp.nc').exists()
            if write or write_flatten:
                da = xr.open_dataset(filename)
                add_properties_to_dataset(
                    da, run_serie, add_property, write = write, write_flatten = write_flatten)
    
    if write_properties_dataset:
        logger.info("Write properties datasets")
        for key, value in add_property.items(): # don't compute something that already exists unless explicitely asked
            add_property[key] = value and (not (path_to_dataset / f"{key}.comp.nc").exists() or overwrite)
        read_properties("*", add_property, write_properties_to_file=True)

    logger.info("End of %s", __file__)


if __name__ == "__main__":
    add_property = {'rdf': True, 'msd': True,
                    'msd_by_run_id': True, 'elastic': False, 'pore': True, 'cn': True, 'bad':True, 'bad_by_cn':True, 'reduction': True, 'ring':True}
    write_run_serie_dataset = True                    
    try:
        if sys.argv[1] in properties_dict.keys():
            add_property = properties_dict[sys.argv[1]]
            if sys.argv[1] == 'elastic':   
                write_run_serie_dataset = False                       
    except:
        pass                 
    main(add_property, write_run_serie_dataset=write_run_serie_dataset, write_properties_dataset=False, overwrite = False)
    main(add_property, write_run_serie_dataset=False, write_properties_dataset=True, overwrite = True)
