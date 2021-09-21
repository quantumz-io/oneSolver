usage()
{
	echo "usage help"
}

if [[ $# -eq 0 ]]; then
	usage
	exit 1
fi

gpu=0
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
     -h|--help)
         usage
	 exit
	;;
     esac
done

if [ "$gpu" = "1" ]; then
#    echo "Free nodes with gpu"
    freenodes=$(pbsnodes -l free :gpu | awk '{print $1}')
else
#    echo "Free nodes without gpu"
    freenodes=$(pbsnodes -l free | awk '{print $1}')
fi
declare -a listofnodes=( $freenodes )

nnodes=${#listofnodes[@]}
randomnode=$((RANDOM * ($nnodes + 1) /37668))

rnodename=${listofnodes[$randomnode]}
#echo "name of node: $rnodename"


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
#	echo "${listofnodes[$i]}"
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
#	echo "${listofnodes[$i]}"
	combine+=${listofnodes[$i]}+
    done
    echo "-l nodes=${combine%?}:gpu:ppn=2"
elif (( $neednodes > $nnodes )); then
    echo "To many nodes. Number of availible nodes: $nnodes"
else
    echo "Something went wrong"
fi
