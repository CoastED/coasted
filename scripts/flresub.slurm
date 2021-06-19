#!/bin/bash

program="flresub"
ncores="$1"
flmlarg="$2"
jobname="$3"

execution_time="12:00"
sparse_flag=""



export SLURM_ACCOUNT="XXXX"     # Replace with your account name
queue_name="normal"		# Replace with queue name

# project_group="d67"



# Always set this so bc does floating point
export BC_ENV_ARGS="-l"


if [ "$jobname" == "" ]
then
    echo "usage: $program <ncores> <flml prefix> <job name> [options]"
    echo "	    [ -t #_openmp_threads (def.=$omp_threads) ]"
    echo "	    [ -s ]"
    echo "	    [ -x execution_time (format=HH:MM) (def.=$execution_time) ]"
#    echo "	    [ -p project_group (def.=$project_group) ]"
#    echo "	    [ -q queue_name (def.=standard. | short, long, low) ]"
    echo ""
    exit 1
fi



OPTIND=4
while getopts "t:sx:" option
do
    case "${option}"
        in
        	x) execution_time=${OPTARG};;
                s) sparse_flag="sparse";;
                t) omp_threads=${OPTARG};;
    esac
done


true_omp_threads="$omp_threads"
if [ "$sparse_flag" == "sparse" ]
then
    node_width=48
    socket_width=24
else
    node_width=96
    socket_width=48
fi



# Hardware limits


# calculated 
mpi_per_socket=`echo "$socket_width/$omp_threads" | bc | sed 's/\..*$//g'`
mpi_per_node=`echo "$node_width/$omp_threads" | bc | sed 's/\..*$//g'`
total_cores=`echo "$ncores * $omp_threads" | bc | sed 's/\..*$//g'`

nremainder=`echo "$total_cores/$node_width" | bc | sed 's/.*\(\..*$\)/0\1/g'`

testzero=`echo "$nremainder==0" | bc`



if [ "$testzero" == "0" ]
then
    nnodes=`echo "($total_cores + $node_width)/$node_width" | bc | sed 's/\..*$//g'`
else
    nnodes=`echo "$total_cores/$node_width" | bc | sed 's/\..*$//g'`

fi

if [ "$sparse_flag" == "sparse" ]
then
    charged_cores=`echo "$nnodes * $node_width * 2" | bc | sed 's/\..*$//g'`
else
    charged_cores=`echo "$nnodes * $node_width" | bc | sed 's/\..*$//g'`
fi


flprefix="$flmlarg"
flfiles=`\ls -f "$flprefix"_*.flml 2> /dev/null`


if [ "$flfiles" == "" ]
then
    flfiles=`\ls -f "$flprefix".flml 2> /dev/null`
    if [ "$flfiles" == "" ]
    then
	echo "*** error: $program cannot find FLML files matching '""$flmlarg""_*.flml'"
	exit 1
    fi
fi


last_file=""
last_simtime=-1
last_finishtime=-1

for flmlfile in $flfiles
do
    simtime=`cat $flmlfile | grep -A1 "<current_time>" | tail -n 1 - | sed 's/.*<real_value.*>\(.*\)<\/real_value>.*/\1/g'`

    finishtime=`cat $flmlfile | grep -A1 "<finish_time>" | tail -n 1 - | sed 's/.*<real_value.*>\(.*\)<\/real_value>.*/\1/g'`
    
    if [ "$simtime" != "" ] && [ "$finishtime" != "" ]
    then
	if [ `echo "$simtime > $last_simtime" | bc` == 1 ]
	then
	    last_file="$flmlfile"
	    last_simtime=$simtime
	    last_finishtime=$finishtime
	fi
    fi
done



if [ `echo "$last_simtime >= $last_finishtime" | bc ` == 1 ] \
    || [ "$last_simtime" == "-1" ]
then
    echo "**** END OF SIMULATION"
    echo " current_time: $last_simtime"
    echo " finish_time : $last_finishtime"
    echo " flml_file   : $last_file" 
    echo ""
    echo "Exiting flresub"

    exit 0
fi



echo "**** RESTARTING SIMULATION"
echo "Still not reached finish_time ($last_finishtime)."
echo "Using $last_file"
echo "current_time = $last_simtime"



datfile=$flprefix.detectors.dat
ndats=`seq 0 100`

if test -f $datfile
then
    echo "*** Backing up $datfile"

    for num in $ndats
    do
        if ! test -f $datfile.$num
        then
            echo "Moving to $datfile.$num"
            mv -f $datfile $datfile.$num
            break
        fi
    done
fi




qsubfile=/tmp/.qsub.$USER.$$
# jobname="lill_t=$last_simtime"

fluidityExec=$FLUIDITY_HOME/bin/fluidity

if [ "$FLUIDITY_INTEL" == "yes" ]
then
    fluidityModule="fluidity-intel"
else
    fluidityModule="local-fluidity-gcc"
fi



wkdir="`pwd`"
rm -f slurm-*.{out,err}

echo "#!/bin/bash --login
#SBATCH -J $jobname
#SBATCH -D $wkdir
#SBATCH --error=slurm-%j.err
#SBATCH --time=0-$execution_time
#SBATCH --nodes=$nnodes
#SBATCH --ntasks=$ncores

# #SBATCH --threads-per-core=$true_omp_threads
#SBATCH --threads-per-core=2
#SBATCH --nvram-options=1LM:1000

# This may or may not do something. Set to 50 Mb for unformatted

# export GFORTRAN_UNFORMATTED_BUFFER_SIZE=50000000
# export GFORTRAN_FORMATTED_BUFFER_SIZE=1000000

# OpenMP

export OMP_NUM_THREADS=$true_omp_threads
export OMP_PLACES=threads
export OMP_PROC_BIND=close

# MVAPICH tweaks

export MV2_ON_DEMAND_THRESHOLD=1
export MV2_CPU_BINDING_POLICY=hybrid
export MV2_HYBRID_BINDING_POLICY=linear
export MV2_ENABLE_AFFINITY=0
export MV2_USE_SHARED_MEM=1
export MV2_ALLGATHER_REVERSE_RANKING=1

# # PSM2 stuff

export PSM2_MULTI_EP=1
export PSM2_MULTIRAIL=1
export PSM2_MULTIRAIL_MAP=0:1,1:1

export FI_PROVIDER=psm2
export FI_PSM2_LAZY_CONN=0
# export PSM2_DEVICES=self,shm,hfi
export PSM2_DEVICES=self,hfi,shm

#export I_MPI_HYDRA_TOPOLIB=
#export I_MPI_DEBUG_OUTPUT=debug_output.txt
#export I_MPI_DEBUG=100


# export PSM2_MEMORY=large
# export PSM2_RANKS_PER_CONTEXT=4
# export HFI_DISABLE_MMAP_MALLOC=1

export FI_PSM2_LOCK_LEVEL=0
export PSM2_RCVTHREAD=0

module load local-fluidity-gcc


# Override system malloc for faster memory alloc/dealloc
export LD_PRELOAD=\$CONTRIB_DIR/base/1.0/lib/libtcmalloc.so

# The following take a copy of the Fluidity Python directory and 
# put it in the current directory. If we don't do this, we get import errors. 
export WORKING_DIR=\$(pwd -P)

rm -rf .python
mkdir -p .python/fluidity; mkdir .python/tidetools; mkdir .python/windtools

cp -rf \$FLUIDITY_HOME/python/* .python/fluidity
cp -f \$HOME/windtools/*.py .python/windtools
cp -f \$HOME/tidetools/*.py .python/tidetools

# copy turbine previous state to backup
# cp -f turbines.state turbines.state.bak 2>/dev/null

cd $wkdir

export PYTHONPATH=.python/fluidity:.python/windtools:.python/tidetools:python:\$PYTHONPATH


# MVAPICH2
# srun --distribution=block:block \$HOME/fluidity/bin/fluidity -l -v1 $last_file

srun --distribution=block:block \$HOME/fluidity/bin/fluidity -l $last_file

## clean up the python directory
# rm -rf .python

" > $qsubfile

# Make a copy for debugging

rm -f qsub.[0123456789]*
cp $qsubfile qsub.$$

echo ""
echo "--=----=----=----=----=----=----=----=----=----=----=----=----=----=----=----=--"
echo "qsub: submitting job $jobname -"
echo "   cores: $total_cores"
echo " charged: $charged_cores"
echo ""
echo " runtime: $run_hours""h"
echo " srun \$HOME/fluidity/bin/fluidity -l -v1  $last_file"

if [ "$sparse_flag" == "sparse" ]
then
    echo "*** sparsely populated"
fi
echo ""

rm -f $jobname.{e,o}*
rm -f slurm-*.out

# Submit job
sbatch $qsubfile


echo "--=----=----=----=----=----=----=----=----=----=----=----=----=----=----=----=--"

# If anything goes amiss, halt submissions
if [ $? -eq 1 ]
then
    echo "**** error: sbatch $qsubfile failed. Aborting."
    exit 1
fi

rm -f $qsubfile
