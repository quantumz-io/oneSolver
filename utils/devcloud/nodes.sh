usage()
{
echo "
NAME
       nodes.sh - check availiblity, show command and allocate nodes 

SYNOPSIS
       ./nodes.sh [OPTION] [OPTION FOR QSUB]...

DESCRIPTION
      Script allows to check availible computational nodes on the Intel DevCloud
      computational resources. This script has its own set of options and allows
      to redirect additional options to the qsub resource manager.

      Options for the script itself are not obligatory but the options for the
      qsub manager have to be after the script options.

      The scripts of options:
      
       -n NUM, --nodes NUM
              the NUM is the number of nodes which should be asked for allocation

       -g, --gpu
              when the flag is present the nodes have to had graphics cards

       -a, --alloc
              when the flat is present the script should try to allocate
              the required number of nodes
	      
       -h, --help
             display this help and exit

The rest options are redirected without any changes on the to the qsub resource manager.

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
allc=0

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
     -a|--allc)
         allc=1
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
elif (( $neednodes == 1 )) && (( $gpu == 0 )) && (( allc == 0 )); then
    echo "To allocate run command: qsub -l nodes=$rnodename:ppn=2"
elif (( $neednodes == 1 )) && (( $gpu == 0 )) && (( allc == 1 )); then
    qsub -l nodes=$rnodename:ppn=2 $@
elif (( $neednodes == 1 )) && (( $gpu == 1 )) && (( allc == 0 )); then
    echo "To allocate run command: qsub -l nodes=$rnodename:gpu:ppn=2"
elif (( $neednodes == 1 )) && (( $gpu == 1 )) && (( allc == 1 )); then
    qsub -l nodes=$rnodename:gpu:ppn=2 $@
elif (( $neednodes > 1 )) && (( $neednodes <= $nnodes )) && (( $gpu == 0 )) && (( allc == 0 )); then
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
    echo "To allocate run command: qsub -l nodes=${combine%?}:ppn=2"
elif (( $neednodes > 1 )) && (( $neednodes <= $nnodes )) && (( $gpu == 0 )) && (( allc == 1 )); then
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
    qsub -l nodes=${combine%?}:ppn=2 $@
 elif (( $neednodes > 1 )) &&  (( $neednodes <= $nnodes )) && (( $gpu == 1 )) && ((allc == 0 )); then
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
    echo "To allocate run command: qsub -l nodes=${combine%?}:gpu:ppn=2"
elif (( $neednodes > 1 )) &&  (( $neednodes <= $nnodes )) && (( $gpu == 1 )) && ((allc == 1 )); then
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
    echo "To allocate run command: qsub -l nodes=${combine%?}:gpu:ppn=2"
    qsub -l nodes=${combine%?}:gpu:ppn=2 $@
elif (( $neednodes > $nnodes )); then
    echo "To many nodes requested. Number of availible nodes: $nnodes"
else
    echo "Something went wrong"
fi
