
#PBS -lselect=1:ncpus=16:mem=20gb
#PBS -lwalltime=72:00:00

# Load modules for any applications

module load anaconda3/personal
export MESA_DIR=$HOME/MESA/mesa
export MESASDK_ROOT=$HOME/MESA/mesasdk

source $MESASDK_ROOT/bin/mesasdk_init.sh

export OMP_NUM_THREADS=16

# Change to the temporary submission directory

cd $PBS_O_WORKDIR

# Run program
# must point directly to the location of the run-script

python $EPHEMERAL/Young_planets_2/Model_0/pampas_run.py 0 0 39 1

# copy the files to WORK directory for storage

#mkdir $WORK/$PBS_JOBID
#cp -r * $WORK/$PBS_JOBID
 
#output files now in WORK with JOBID as folder name
