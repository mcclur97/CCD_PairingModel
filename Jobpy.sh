#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########

#SBATCH --time=119:59:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1                 # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=1           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name julia         # you can give your job a name for easier identification (same as -J)
#SBATCH -q normal

######### Command Lines to Run #########

for pnum in 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40
do
	echo $pnum $hnum
	python3 Test.py $pnum $hnum
done

scontrol show job $SLURM_JOB_ID     ### write job information to output file#

