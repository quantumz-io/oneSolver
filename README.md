# oneSolver
Created by hperQ - Ewa Hendzel

## Build instructions

```
mkdir build && cd build
cmake .. 
cmake --build .
```

## Running the  applications
Binaries are ``one-solver-anneal`` and ``one-solver-exhaustive`` created
in the app folder. Binary for tests is located in ``tests``.

The binary named `one-solver-anneal` allows to solve a QUBO problem using
the simulated annealing algorithm whereas `one-solver-exhaustive` uses
the brute force search method.

In the ``tests`` folder the application called `main-tests` is build.

In order to obtain help about the applications arguments 
use the `--help` switch. For example, the command
`./one-solver-anneal --help` prints:

```
Allowed options:
  --help                                produce help message
  --input arg                           input file
  --output arg                          output file
  --num-iter arg (=100)                 number of iterations of the algorithm
  --num-tries arg (=100)                number of trajectories to try
  --schedule-type arg (=geometric)      type of beta schedule tu use, either 
                                        linear or geometric
  --beta-min arg (=0.10000000000000001) minimum value of beta in the annealing 
                                        schedule (default 0.1)
  --beta-max arg (=1)                   maximum value of beta in the annealing 
                                        schedule (default 1.0)
  --device-type arg (=host)             device type to use (cpu, gpu or host)
```

The command `./one-solver-exhaustive --help` prints:

```
Allowed options:
  --help                    produce help message
  --input arg               input file
  --output arg              output file
  --device-type arg (=host) device type to use (cpu, gpu or host)
```

The following applications call

```
./one-solver-exhaustive --input ../../examples/test1.qubo --output result1.qubo
```

executes the program for the model described in ``test1.qubo`` file and 
the obtained results are written to ``result1.qubo`` file. 

## Original goal of the project

The project\'s goal is to develop a spin-glass (i.e., Ising model)
solver based on a family of simulated annealing algorithms for large
system sizes and an exhaustive search for small Ising instances. This
solver will be implemented using an oneAPI development stack to utilize
massive parallel capabilities of present-day computing architectures.

It is well established that NP-hard problems can be solved by cleverly
mapping them to the classical Ising Hamiltonians (i.e., spin-glasses).
With the rapid development of both quantum and classical annealers, the
scientific and computer engineering communities have intertwined once
again in their efforts to reformulate many real-life combinatorial
optimization problems using the ideal of interacting spin variables
(i.e., spins). Various dedicated devices and computer chips are being
developed currently to tackle the very problem of finding the low-energy
configurations (i.e., spectrum) of the Ising models. For instance,
Fujitsu and Hitachi have produced their digital and analog classical
annealers. In stark contrast, D-Wave Inc. has manufactured its quantum
annealers, which promises to solve the Ising instances using the law of
quantum physics; the famous adiabatic theorem to be more specific.

Although there are many sophisticated methods (e.g., tensor networks,
simulated quantum annealing, path integrals, etc.) to approach the
aforementioned problems, not that many of them can benefit from the
High-Performance Computing (HPC) available today. On the other hand,
fundamental techniques such as the Brute-Force (that simply executes an
exhaustive search) or simulated annealing (that tries to mimic a famous
physical process of annealing known in the material science) can take
advantage of highly parallel architectures (graphics card, many-cores).
This project is devoted precisely for that purpose; i.e., we aim to
provide conceptually simple solvers that can be executed concurrently on
present-day, massively parallel, classical architectures in an agnostic
manner.

Literature:
1. [Brute-forcing spin-glass problems with CUDA](https://arxiv.org/abs/1904.03621)
2. [Demonstration of a scaling advantage for a quantum annealer over simulated annealing](https://arxiv.org/abs/1705.07452)
3. [Parallel in time dynamics with quantum annealers](https://www.nature.com/articles/s41598-020-70017-x)
