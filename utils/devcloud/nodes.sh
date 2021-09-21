usage()
{
	echo "usage help"
}

if [[ $# -eq 0 ]]; then
	usage
	exit 1
fi

gpu=0
allc=0
#echo "Arguments are: $@"
#argc="$@"
remargs=()

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

echo "Remainded args: $@" 

if [ "$gpu" = "1" ]; then
    echo "Free nodes with gpu"
    freenodes=$(pbsnodes -l free :gpu | awk '{print $1}')
else
    echo "Free nodes without gpu"
    freenodes=$(pbsnodes -l free | awk '{print $1}')
fi
declare -a listofnodes=( $freenodes )

nnodes=${#listofnodes[@]}
randomnode=$((RANDOM * ($nnodes + 1) /37668))

rnodename=${listofnodes[$randomnode]}
echo "name of node: $rnodename"


if (( $neednodes  <= 0 )); then
    echo "Bad value of nodes number"
elif (( $neednodes == 1 )) && (( $gpu == 0 )) && (( allc == 0 )); then
    echo "needed one node"
    echo "availble nodee = $rnodename"
    echo "allocation command: qsub -l nodes=$rnodename:ppn=2"
elif (( $neednodes == 1 )) && (( $gpu == 0 )) && (( allc == 1 )); then
    echo "needed one node"
    echo "availble nodee = $rnodename"
    echo "allocation command: qsub -l nodes=$rnodename:ppn=2"
    qsub -I -l nodes=$rnodename:ppn=2 -$@
elif (( $neednodes == 1 )) && (( $gpu == 1 )) && (( allc == 0 )); then
    echo "One node without gpu"
    echo "allocation command: qsub -l nodes=$rnodename:gpu:ppn=2"
elif (( $neednodes == 1 )) && (( $gpu == 1 )) && (( allc == 1 )); then
    echo "One node with gpu"
    echo "allocation command: qsub -l nodes=$rnodename:gpu:ppn=2"
    qsub -I -l nodes=$rnodename:gpu:ppn=2
elif (( $neednodes > 1 )) && (( $neednodes <= $nnodes )) && (( $gpu == 0 )) && (( allc == 0 )); then
    echo "needed more than one node"
    if (( $randomnode - $neednodes < 0 )) || (( $randomnode + $neednodes > $nnodes )); then
	beginnode=0
	endnode=$neednodes
    else
	beginnode=$randomnode
	endnode=$(($randomnode + $neednodes - 1))
    fi
    echo "Begin: $beginnode, end: $endnode"
    for ((i=$beginnode; i<=$endnode; i++ )); do
	echo "${listofnodes[$i]}"
	combine+=${listofnodes[$i]}+
    done
    echo "Combined: ${combine%?}"
    echo "allocatin command: qsub -l ${combine%?}:ppn=2"
elif (( $neednodes > 1 )) && (( $neednodes <= $nnodes )) && (( $gpu == 0 )) && (( allc == 1 )); then
    echo "needed more than one node"
    if (( $randomnode - $neednodes < 0 )) || (( $randomnode + $neednodes > $nnodes )); then
	beginnode=0
	endnode=$neednodes
    else
	beginnode=$randomnode
	endnode=$(($randomnode + $neednodes - 1))
    fi
    echo "Begin: $beginnode, end: $endnode"
    for ((i=$beginnode; i<=$endnode; i++ )); do
	echo "${listofnodes[$i]}"
	combine+=${listofnodes[$i]}+
    done
    echo "Combined: ${combine%?}"
    echo "allocatin command: qsub -l nodes=${combine%?}:ppn=2"
    qsub -I -l nodes=${combine%?}:ppn=2
elif (( $neednodes > 1 )) &&  (( $neednodes <= $nnodes )) && (( $gpu == 1 )) && ((allc == 0 )); then
    echo "needed more than one node"
    if (( $randomnode - $neednodes < 0 )) || (( $randomnode + $neednodes > $nnodes )); then
	beginnode=0
	endnode=$neednodes
    else
	beginnode=$randomnode
	endnode=$(($randomnode + $neednodes - 1))
    fi
    echo "Begin: $beginnode, end: $endnode"
    for ((i=$beginnode; i<=$endnode; i++ )); do
	echo "${listofnodes[$i]}"
	combine+=${listofnodes[$i]}+
    done
    echo "Combined: ${combine%?}"
    echo "allocatin command: qsub -l nodes=${combine%?}:gpu:ppn=2"
elif (( $neednodes > 1 )) &&  (( $neednodes <= $nnodes )) && (( $gpu == 1 )) && ((allc == 1 )); then
    echo "needed more than one node"
    if (( $randomnode - $neednodes < 0 )) || (( $randomnode + $neednodes > $nnodes )); then
	beginnode=0
	endnode=$neednodes
    else
	beginnode=$randomnode
	endnode=$(($randomnode + $neednodes - 1))
    fi
    echo "Begin: $beginnode, end: $endnode"
    for ((i=$beginnode; i<=$endnode; i++ )); do
	echo "${listofnodes[$i]}"
	combine+=${listofnodes[$i]}+
    done
    echo "Combined: ${combine%?}"
    echo "allocatin command: qsub -l nodes=${combine%?}:gpu:ppn=2"
    qsub -I -l nodes=${combine%?}:gpu:ppn=2
elif (( $neednodes > $nnodes )); then
    echo "To many nodes. Number of availible nodes: $nnodes"
else
    echo "Something went wrong"
fi
