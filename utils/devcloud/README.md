# Script for allocation nodes on Intel DevCloud

For submitting jobs on Intel(R) DevCloud computing nodes, the _qsub_ utility of the Torque resource manager should be used. But at the time of writing this tutorial, there was some issue with using the native way of this utility.

As a workaround, scripts have been developed which allow us to use particular features of the Intel(R) DevCloud computational cluster.


For the rest of this instruction, it is assumed that the commands are  executed from the main project folder  _oneSolver_.

## Stand alone solution

The script _nodes.sh_ allows to check the availability of free computational nodes, allocate and execute the program on the Intel(R) DevCloud cluster.

The command to execute this script is:

```bash
./utils/devcloud/nodes.sh
```

The following command shows the allowed options for the script.

```bash
-n|--nodes NUM --- allows to allocate NUM nodes
-g|--gpu       --- allows to allocate nodes with graphics card
-a|--alloc     --- allows to allocate free nodes with the *qsub* utility for running the executable
-h|--help      --- shows help on using the script
```

This script has an additional feature that allows us to pass additional flags without modifying them to the _qsub_ utility.

Examples of execution:

To check if there are three free nodes, run the following command:
```bash
 ./utils/devcloud/nodes.sh -n 3
```
The output is:

```bash
To allocate run command: qsub -l nodes=s001-n201+s001-n203+s001-n205:ppn=2
```

To check if there are four free nodes with graphics cards, run the following command :

```bash
./utils/devcloud/nodes.sh -n 4 -g
```

The output is:
 
```bash
To allocate run command: qsub -l nodes=s001-n215+s001-n216+s001-n218+s001-n219:gpu:ppn=2
```

To check if there are five free nodes with gpu and allocate them in the interactive mode
with _qsub_ utility:

```bash
 ./utils/devcloud/nodes.sh -n 5 -g -a -I
```

If the required number of free nodes are available, nodes would be allocated for computations in interactive mode.

## Solution using backticks feature of Linux shell

The Linux shell has a feature that allows redirecting the output of a command in the given command chain. The second script _qsnodes.sh_ uses this feature along with _qsub_ resource manager to allocate the required number of nodes.

The command to execute this script is:

```bash
./utils/devcloud/qsnodes.sh
```
The output of this command shows the allowed options which could be passed to the script and exit to shell. Allowed options for this script are very similar to those of the _nodes.sh_ script:

```bash
-n|--nodes NUM --- allows to allocate NUM nodes
-g|--gpu       --- allows to allocate nodes with the graphics card
-h|--help      --- shows help about using the script
```

Examples of execution:

To allocate six free nodes equipped with graphics cards, run the following command:

```bash
./utils/devcloud/qsnodes.sh -n 6 -g
```
The output will be a command which could be directly passed to the _qsub_ utility:

```bash
-l nodes=s001-n218+s001-n219+s001-n220+s001-n221+s001-n222+s001-n223:gpu:ppn=2
```

So to allocate six free nodes equipped with graphics cards for computation the command should look like:

```bash
qsub `./utils/devcloud/qsnodes.sh -n 6 -g` -I
```
Take care about the "`" backticks used in the command. The rest of the options in this command is passed to the _qsub_ resource manager. If six free nodes equipped with graphics cards are available then they will be allocated for computations.