import xarray as xr
import logging
import numpy as np

logger = logging.getLogger(__name__)

def write_dataset(data, filename, compress = True):
    if not compress:
        logger.info("Write netCFD file for %s", filename)
        data.to_netcdf(filename)
    else:
        logger.info("Write compressed netCFD file for %s", filename)
        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in data.data_vars}
        # potentially empty coord that are considered objects and cause netcfd writing to fail
        if 'ring_var' in data.coords: 
            data['ring_var'] = data['ring_var'].astype('str')
        data.to_netcdf(filename, encoding=encoding)


def concatenate_str(x, column_names):
    """return new columns which concatenates the columns identified by column_names for run_id row x of dataset"""
    return xr.DataArray('-'.join([str(x[name].item()) for name in column_names]))

def concatenate_run_id(x, column_names):
    """return new columns which concatenates the run_id (int, that will be formated with 3digits) identified by column_names for run_id row x of dataset"""
    return xr.DataArray('-'.join([format(int(x[name].item()), '03d') for name in column_names]))


def flatten(x):
    """
    Return x flattened by 'run_part' and 'run_subpart'
    Take the max of every variable except for 
        first_step: min
        computed properties: mean
    For uniform quanitites (e.g. config) max will take the constant value
    """
    x = x.unstack('run_id')
    var_names = list(x.data_vars.keys())
    for var in var_names:
        if var in ['first_step']:
            x[var] = x[var].min(['run_part', 'run_subpart'])
        elif var in ['rdf', 'bad']:
            x[var] = x[var].mean(['run_part', 'run_subpart'])
        elif var in ['msd']:
            x = x.drop_dims('Time') # as Time is only used for MSD so far
        else:
            x[var] = x[var].max(['run_part', 'run_subpart'])
    x = x.drop_dims(['run_part', 'run_subpart'])
    x = x.swap_dims({'run_serie_exp': 'run_id'}
                    ).set_index(run_id='run_serie_exp')
    return x

def flatten_dataset(da):
    """
    Args:
        da: xarray dataset with run_id = run_serie, exp, part and subpart
    Return:
        da: xarray dataset with run_id = run_serie, exp 
    """
    da['run_serie_exp'] = da.groupby("run_id").map(
        concatenate_run_id, args=[['run_serie', 'run_exp']])
    if 'Step' in da.coords:
        da = da.chunk({'Step': 1000000}) # turn to dask array of size 1 million steps
    da = da.set_index(run_id=['run_serie_exp', 'run_part', 'run_subpart'])
    da = da.groupby("run_serie_exp").map(flatten)
    da = da.load() # load into memory for further use as a standard xarray
    return da

def compute_bound(a, c):
    """return a/c if a divisible by c and else floor(a/c) + 1"""
    bound = a // c
    if bound * c != a:
        bound += 1
    return bound


def get_interval(a, b, c, sampling_rate=None, endpoint=False):
    """returns np arange of x in [a;b[ separated by c so that each x is divisible by c, sampled by sampling rate

    The endpoint of the interval can optionally be included.

    e.g. (0, 1000, 100) will give [  0, 100, 200, 300, 400, 500, 600, 700, 800, 900]
         (0, 1000, 100, 2) will give [  0, 200, 400, 600, 800]
         (50, 950, 100, 2) will give [100, 300, 500, 700, 900] 
    """
    if endpoint == False:
        return np.arange(compute_bound(a, c) * c, compute_bound(b, c) * c, c)[::sampling_rate]
    else:
        return np.arange(compute_bound(a, c) * c, (compute_bound(b, c) + 1) * c, c)[::sampling_rate]


def get_traj_index(a, b, c, sampling_rate = None, endpoint=False):
    """
    a, b, c expressed in Steps
    sampling_rate: int
    The endpoint of the interval can optionally be included.

    returns traj index of number of frame in traj in [a;b[ separated by c so that each x is divisible by c, sampled by sampling rate
    e.g. (0, 1000, 100) will give slice(0, 10, None) ie frames indexes [0, 1, 2, ... 9] and step [0, 100, 200, ... 900]
         (0, 1000, 100, 2) will give slice(0, 10, 2) ie frames indexes [0, 2, 4 ... 8] and step [0, 200, 400, ... 800]
         (50, 950, 100, 2) will give slice(1, 10, 2) ie frames indexes [1, 3, 5 ... 9] and step [100, 300, 500, ... 900]
    """
    if endpoint == False:
        return slice(int(compute_bound(a, c)), int(compute_bound(b, c)), sampling_rate)
    else:
        return slice(int(compute_bound(a, c)), int(compute_bound(b, c)) + 1, sampling_rate)
