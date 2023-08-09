# amof-workflow

An example workflow to use [`amof`](https://github.com/coudertlab/amof), a python package to analyze Molecular Dynamics (MD) trajectories of amorphous Metal-Organic Frameworks (MOFs). 

It is a simplified version of the workflow I used for most on [my PhD](https://www.theses.fr/s255984). Unlike `amof`, it is not constituted of generic functions but rather of problem specific and hardware specific functions.
Therefore it cannot be used straigh away, and require adaptation (e.g. to your supercomputing clusters for the MD runs).
It can nonetheless be used as inspiration or to speed-up the set-up of a functional workflow.


## Structure

Runs (where the MD are ran) and analyses are separated in different folders to make sure no modification is applied on the MD output when analysing them.

### Runs

Runs are stuctured by flavor of MD (e.g. `001-nequip`).
For a given flavor of MD, each folder (e.g.`001`) designate a `run_serie`, i.e. a set of MD runs that were launched simultaneously (e.g. different pressures with same parameters).
The description of all run_series is contained in the `XXX-runserie_description.csv`.

### Analyses

Analyses are separated in two folders: `data` (actual computation of the properties), and `run_results` (jupyter notebooks used to look at those properties).


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

#### Look at the results in jupyter notebooks

Open `analysis/run_results/001-nequip/001 - look at single MD run.ipynb` with jupyter notebook.

<!-- TODO:

### Parameter variation

Give an example of combination and  linear

### Volume deformation

several examples : 
a simple one
combination-successive
both-directions

### Melt-quenching

### Different MD schemes

#### NequIP with ASE

Simple example

#### AIMD

Simple

#### MOF-FF

Same as ReaxFF, one example of crystal and one of glass
One NPT and one NVT to show different ensembles
Link to citable data


#### ReaxFF

Simple example, say that can do the same that nequip for inspiration
Link to citable data -->



### Create another new run

For example with NequIP, copy one of the desired example with a new run_serie in `runs/001-nequip`, and add the corresponding line in `001-runserie_description.csv`.
Modify the scripts as desired and run.

