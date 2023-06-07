# amof-workflow

An example workflow to use [`amof`](https://github.com/coudertlab/amof), a python package to analyze Molecular Dynamics (MD) trajectories of amorphous Metal-Organic Frameworks (MOFs). 

It is a simplified version of the workflow I used for most on [my PhD](https://www.theses.fr/s255984). Unlike `amof`, it is not constituted of generic functions but rather of problem and hardware specific functions.
It can nonetheless be used as inspiration or to speed-up the set-up of a functional workflow.


## Structure

principles and how is it structured

TODO

## Examples

### Simple NVT run

100ps of NVT at 300K consisting of a single restart.

#### Generate inputs for MD

In `runs/001-nequip/001`

Input the desired settings in `run_lammps.ini`, change the input lammps data file `initial.lmp` and MLP model in `raw_files`.

`template/job_zay.slurm` is dependant on the cluster it is launched on, and should be updated (e.g. with the correct path to the LAMMPS installation).

Then launch `run_lammps.py`:
```
python run_lammps.py
```

When the simulation is over, concatenate the multiple restarts (only one in this example) and compress unused files.

Make sure the `include_in_dataset` column in `001-runserie_description.csv` is set to `True`.

```
python concatenate_runs.py compress
```

`concatenate_runs.py` takes as options `compress`, `decompress` or nothing which allows to compress/decompress unused files in addition to concatenating the run.

#### Launch analysis with amof

In `analysis/data/001-nequip`

```
python read_thermo.py && python compute_properties.py && python add_properties_to_thermo.py
```

`compute_properties.py` and `add_properties_to_thermo.py` can take options to only compute certain properties: `fast, fewcores, elastic, nopore`.

The list of computed properties for each option can be found in `shared_variables.py`.
TODO: detail options

#### Look at the results in jupyter notebooks

Open `analysis/run_results/001-nequip/001 - look at single MD run.ipynb` with jupyter notebook.

### Create a new run

Copy one of the desired example with a new run_serie in `runs/001-nequip`, and add the corresponding line in `001-runserie_description.csv`.
Modify the scripts as desired and run.


## Dev

### Current version

Minimal with the simplest nequip job

### Todo

Add an example input/analysis for every sort of MD sim I used

- [ ] MD scheme:
  - [ ] ReaxFF
  - [ ]  MOF-FF
  - [ ] CP2K
  - [ ]  Nequip
- [ ]  Type of run
  - [ ] Equilibration
  - [ ] Volume change
    - [ ] NVT
    - [ ] NPT
  - [ ] Melt-quench

Add a nicer doc:
- [ ] this readme
- [ ] examples

Add extended comments on the python requirements
