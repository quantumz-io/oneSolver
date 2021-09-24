#Script for working with Intel DevCloud allocation nodes

For using IntelDevcloud computing nodes the _qsub_ utility of Torque resource manager should be use. But at the moment of writing this tutorial
there is some troubles with using this utility in native way.

For the workaround the scripts which allows to use full features of Intel DevCloud computational cluster have been prepared.

For the rest of this instruction the assumption is taken as true,
that the commands given in this instruction are executed from the main folder of project which is named _oneSolver_

#Stand alone solution

The script _nodes.sh_ allows to check availability, allocate and
even execute computations on the DevCloud cluster.

The basic command for execution of this script is:

```bash
./utils/devcloud/nodes.sh
```

This command shows the allowed option which could be added during
script execution. These options are:

```bash
-n|--nodes NUM --- which allows to allocate NUM nodes
-g|--gpu       --- switch which allows to allocate nodes which are equipment with the graphics card, without this switch enabled the nodes its not important
-a|--alloc     --- allocate the free nodes with the *qsub* utility for running software
-h|--help      --- shows help about using the script and exit to shell
```
There is important feature which allows to pass additional flag to the script. And these flags are passed without changing them
to the _qsub_ utility.

Examples of execution:

Check if there are three nodes available the script should be executed:
```bash
 ./utils/devcloud/nodes.sh -n 3
```
The example output is:

```bash
To allocate run command: qsub -l nodes=s001-n201+s001-n203+s001-n205:ppn=2
```

Check if there are four nodes equipped with graphics cards are available:

```bash
./utils/devcloud/nodes.sh -n 4 -g
```

The output would be:
To allocate run command: 

```bash
qsub -l nodes=s001-n215+s001-n216+s001-n218+s001-n219:gpu:ppn=2
```

Check if there are five nodes with gpu available and allocate them in the interactive mode
with _qsub_ utility:

```bash
 ./utils/devcloud/nodes.sh -n 5 -g -a -I
```

If there are the required number of nodes available they would be allocated for computations.


#Solution which uses backticks feature of Linux shell

The Linux shell has feature which allows to redirect output of given command in the some commands
chain. The second script which is named _qsnodes.sh_ uses this feature and allows to use _qsub_
resource manager and only use the script for requesting and allocating required number of nodes.

The basic command for execution of this script is:

```bash
./utils/devcloud/qsnodes.sh
```
The output of this command shows the allowed options which could passed to the script and exit to shell.

Allowed options in this script are very similar to these which could be used in the _nodes.sh_:

```bash
-n|--nodes NUM --- which allows to allocate NUM nodes
-g|--gpu       --- switch which allows to allocate nodes which are equipment with the graphics card, without this switch enabled the nodes its not important
-h|--help      --- shows help about using the script and exit to shell
```

The script executed with the _qsub_ utility itself prints if the required nodes are available, for example  
to ask for six nodes equipped with graphics cards:

```bash
./utils/devcloud/qsnodes.sh -n 6 -g
```
The output will be an option which could be directly passed to the _qsub_ utility:

```bash
-l nodes=s001-n218+s001-n219+s001-n220+s001-n221+s001-n222+s001-n223:gpu:ppn=2
```

So when these nodes should be allocated for computation the command should look like:

```bash
qsub `./utils/devcloud/qsnodes.sh -n 6 -g` -I
```
Take care about the "`" backticks which are used in the script command. The rest of this command  
is typical for the using _qsub_ resource manager. If on the cluster six nodes equipped with graphics  
cards are available they will be allocated and granted for computations.

 





