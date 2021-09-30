# Running oneSolver 2.0 on the Intel(R) with PBS

The folder contains example scripts for building and running
oneSolver 2.0 on a cluster equipped with Torque PBS manager.

## Building

The project binaries can be build using `mpibuild.sh` script. The script
should be run the project main folder as shown

```bash
qsub mpibuild.sh
```

If the command is executed without erros the PBS job will be created and
submitted to the job queue. After successful completion of the job the
compiled binaries would be placed in the build directory and its
subdirectories.

## Running the solvers

Files `mpirunex.sh` and `mpirunann.sh` are two examples scripts
for running oneSolver 2.0 using Torque PBS manager for exhaustive and
simulated annealing algorithms respectively.

In order to submit a job for calling the samplers to Torque PBS one has
to adapt the example scripts and execute the following command

```bash
qsub -d . mpirunex.sh
```

In these examples it is assumed that the command is called from
oneSolver directory.
