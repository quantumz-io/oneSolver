usage()
{
echo "
NAME
       qsnodes.sh - check availiblity and allocate nodes in qsub chain command 

SYNOPSIS
       ./qsnodes.sh [OPTION] ...

DESCRIPTION
      Script allows to check availible computational nodes on the Intel DevCloud
      computational resources. This script could be used with qsub resource manger
      with the using backticks feature of linux shell.

      Generally the script could be called with its own arguments in backticks
      and its output is redirected to qsub utility as its arguments.

      Options for the script itself are not obligatory but the options for the
      qsub manager have to be after the script options.

      The scripts of options:
      
       -n NUM, --nodes NUM
              the NUM is the number of nodes which should be asked for allocation

       -g, --gpu
              when the flag is present the nodes have to had graphics cards

       -h, --help
             display this help and exit

SEE ALSO
      qsub utiity for submitting PBS job
      Documentation availible at:
      <http://docs.adaptivecomputing.com/torque/4-0-2/Content/topics/commands/qsub.htm>
"
}

if [[ $# -eq 0 ]]; then
	usage
	exit 1
fi

gpu=0

for arg in $@;
do
     case $arg in
    -n|--nodes)
	shift
        neednodes=$1
	shift
	;;
    -g|--gpu)
	gpu=1
	shift
	;;
     -h|--help)
         usage
	 exit
	;;
     esac
done

if [ "$gpu" = "1" ]; then
    freenodes=$(pbsnodes -l free :gpu | awk '{print $1}')
else
    freenodes=$(pbsnodes -l free | awk '{print $1}')
fi
declare -a listofnodes=( $freenodes )

nnodes=${#listofnodes[@]}
randomnode=$((RANDOM * ($nnodes + 1) /37668))

rnodename=${listofnodes[$randomnode]}

if (( $neednodes  <= 0 )); then
    echo "Bad value of nodes number"
elif (( $neednodes == 1 )) && (( $gpu == 0 )); then
    echo "-l nodes=$rnodename:ppn=2"
elif (( $neednodes == 1 )) && (( $gpu == 1 )); then
    echo "-l nodes=$rnodename:gpu:ppn=2"
elif (( $neednodes > 1 )) && (( $neednodes <= $nnodes )) && (( $gpu == 0 )); then
    if (( $randomnode - $neednodes < 0 )) || (( $randomnode + $neednodes > $nnodes )); then
	beginnode=0
	endnode=$neednodes
    else
	beginnode=$randomnode
	endnode=$(($randomnode + $neednodes - 1))
    fi
    for ((i=$beginnode; i<=$endnode; i++ )); do
	combine+=${listofnodes[$i]}+
    done
    echo "-l nodes=${combine%?}:ppn=2"
elif (( $neednodes > 1 )) &&  (( $neednodes <= $nnodes )) && (( $gpu == 1 )); then
    if (( $randomnode - $neednodes < 0 )) || (( $randomnode + $neednodes > $nnodes )); then
	beginnode=0
	endnode=$neednodes
    else
	beginnode=$randomnode
	endnode=$(($randomnode + $neednodes - 1))
    fi
    for ((i=$beginnode; i<=$endnode; i++ )); do
	combine+=${listofnodes[$i]}+
    done
    echo "-l nodes=${combine%?}:gpu:ppn=2"
elif (( $neednodes > $nnodes )); then
    echo "To many nodes. Number of availible nodes: $nnodes"
else
    echo "Something went wrong"
fi
