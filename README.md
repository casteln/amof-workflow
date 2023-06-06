# amof-workflow
Simplified workflow to use amof to compute properties of MD trajectories

# Structure

principles and how is it structured

TODO

# Example

## Generate inputs for MD

In `runs/001-nequip/001`

Input the desired settings in `run_lammps.ini`, change the input lammps data file `initial.lmp` and MLP model in `raw_files`.

`template/job_zay.slurm` is dependant on the cluster it is launched on, and should be updated (e.g. with the correct path to the LAMMPS installation).

Then launch `run_lammps.py`:
```
python run_lammps.py
```

When the simulation is over, concatenate the multiple restarts (only one in this example) and compress unused files.

```
python concatenate_runs.py compress
```

TODO: detail compress/decompress

## Launch analysis with amof

In `analysis/data/001-nequip`

```
python read_thermo.py && python compute_properties.py && python add_properties_to_thermo.py
```

TODO: detail options

## Look at the results in jupyter notebooks

Open `analysis/run_results/001-nequip/001 - look at single MD run.ipynb` with jupyter notebook.

# Dev

## Current version

Minimal with the simplest nequip job

## Todo

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
