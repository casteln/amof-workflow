# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # 001
#

# %%
import hvplot.xarray # noqa
import panel.widgets as pnw
import xarray as xr
import numpy as np

# %%
# To save plots
import pathlib
from amof.plot import save_hvplot
path_to_plot = pathlib.Path('figures/001') 
path_to_plot.mkdir(parents=True, exist_ok=True) # create dir

# %% [markdown]
# #### import xarray

# %%
path_to_dataset = '../../data/001-nequip/dataset/'

# run_serie_list = [1, 2, 3] 
run_serie_list = [1] 

# %%
list_da = [f"""{path_to_dataset}thermo-prop-flat/run_serie_{format(i, '03d')}.comp.nc""" for i in run_serie_list]
da = xr.open_mfdataset(list_da, parallel=True, concat_dim = 'run_id', data_vars = 'minimal')
da.load()

# %%
# Use the system name as ID, here 'crystal' (can be changed to something more relevant when comparing different runs)
da = da.set_index(run_id = 'system')
da = da.rename({"run_id":"system"})
# Convert temperature to int (not string) as it is used to compute the PMF later
da['temp'] = da['temp'].astype('int')

# %% [markdown]
# ### thermo

# %%
da.thermo.sel(Step = slice(None, None, 100)).interactive().sel(thermo_var=pnw.Select).dropna('Step'
                    ).hvplot(x = 'Step', y = 'thermo', by='system', kind='line')

# %% [markdown]
# #### density

# %%
da.thermo.sel(thermo_var='Density').mean('Step').to_series()

# %% [markdown]
# ### rdf

# %%
da.r.attrs["units"] = 'Ang'

# %%
atom_pair = "Zn-N"
xmin, xmax = 1.5, 6
plot = da.rdf.sel(atom_pair = atom_pair).hvplot.line(x = 'r', by='system', xlim=(xmin, xmax), title = atom_pair)
# save_hvplot(plot, path_to_plot / f'rdf_{atom_pair}')
plot

# %% [markdown]
# Saving plot:

# %%
# As "png"
# save_hvplot(plot, path_to_plot / f'rdf_{atom_pair}', format = 'png')
# As "svg"
# save_hvplot(plot, path_to_plot / f'rdf_{atom_pair}', format = 'svg')
# Both
save_hvplot(plot, path_to_plot / f'rdf_{atom_pair}')

# %%
atom_pair = "Zn-Zn"
xmin, xmax = 3, 7.5
plot = da.rdf.sel(atom_pair = atom_pair).hvplot.line(x = 'r', by='system', xlim=(xmin, xmax), title = atom_pair)
# save_hvplot(plot, path_to_plot / f'rdf_{atom_pair}')
plot

# %% [markdown]
# ### PMF

# %%
# Compute PMF
from scipy.constants import Boltzmann, Avogadro
k_kjpermol = Boltzmann * Avogadro * 10**-3
pmf = - k_kjpermol * da.temp * np.log(da.rdf / da.rdf.max('r')) # Normalized to be equal to 0 at the minimum
pmf.name = 'PMF'
pmf.attrs['units'] = 'kJ/mol'
pmf = pmf.where(xr.ufuncs.isfinite(pmf))

# %%
atom_pair = "Zn-N"
xmin, xmax = 1.5, 6
plot = pmf.sel(atom_pair = atom_pair).hvplot.line(x = 'r', by='system', xlim=(xmin, xmax), title = atom_pair)
# save_hvplot(plot, path_to_plot / f'pmf_{atom_pair}')
plot

# %%
atom_pair = "Zn-Zn"
xmin, xmax = 3, 7.5
plot = pmf.sel(atom_pair = atom_pair).hvplot.line(x = 'r', by='system', xlim=(xmin, xmax), title = atom_pair)
# save_hvplot(plot, path_to_plot / f'pmf_{atom_pair}')
plot

# %% [markdown]
# ### bad

# %%
da.theta.attrs["units"] = 'Deg'

# %%
prop = "N-Zn-N"
xmin, xmax = 1.5, 6
plot = da.bad.sel(atom_triple = prop).hvplot.line(x = 'theta', by='system', title = prop)
# save_hvplot(plot, path_to_plot / f'bad_{prop}')
plot

# %% [markdown]
# ### cn

# %%
# Averaged value
prop = "Zn-N"
da.coordination_number.sel(atom_pair = prop, Step = slice(None, None, 100)).groupby('system').map(lambda x: x.dropna('Step').tail(Step=40000).mean('Step')).to_series()

# %%
# As function of time
da.coordination_number.sel(atom_pair='Zn-N').dropna('Step').rolling(Step = 15, center = True).mean('Step').hvplot(x = 'Step', by='system')

# %% [markdown]
# ### pore

# %%
av = da.pore.sel(pore_var = 'AV_cm^3/g').dropna('Step', how='all') * 1000
av.attrs["units"] = "cm^3/kg"
av.name = 'Accessible volume'
prop = av.name
plot = av.hvplot(y = prop, by='system', kind='hist', bins=50, alpha = 0.7)
# save_hvplot(plot, path_to_plot / f'{prop}')
plot

# %%
tv = da.pore.sel(pore_var=['AV_cm^3/g','NAV_cm^3/g']).dropna('Step', how='all').sum('pore_var') * 1000
tv.attrs["units"] = "cm^3/kg"
tv.name = 'Total porous volume'
prop = tv.name

# %%
plot = tv.hvplot(y = tv.name, by='system', kind='hist', bins = 50, alpha = 0.5)
# save_hvplot(plot, path_to_plot / f'{prop}')
plot

# %% [markdown]
# ### ring

# %%
# as function of time
da.ring.interactive().sel(ring_var=pnw.Select).dropna('Step', how='all').hvplot(x = 'Step', y = 'ring', by='ring_size', kind = 'line')

# %%
# averaged
da.ring.interactive().sel(ring_var=pnw.Select).dropna('Step', how='all').mean('Step').hvplot(x = 'ring_size', y = 'ring', kind = 'bar')

# %% [markdown]
# ### msd

# %%
da.msd.attrs['units'] = 'Ang^2'
da.msd.dropna('Time').hvplot(x='Time', y ='msd', by = 'atom')

# %%
prop = "Zn"
plot = da.msd.sel(atom=prop).dropna('Time').hvplot(x='Time', y ='msd', by = 'system', title = prop)
# save_hvplot(plot, path_to_plot / f'msd_{prop}')
plot

# %%
