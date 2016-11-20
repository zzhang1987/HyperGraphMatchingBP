#!/bin/sh
#SBATCH --job-name=house_mmodes    # Job name for easy identification in queue
#SBATCH --time=72:00:00      # Time requested (HH:MM:SS, or DAYS-HH:MM:SS)
#SBATCH --nodes=1           # Request a single node only
#SBATCH --ntasks=32         # Job will consist of only one process
#SBATCH --cpus-per-task=1   # Number of CPUs per process (for multithreaded programs)
#SBATCH --mem=48GB           # Memory requested per nod#!/bin/sh
#python RunMModesCar.py 2> Error.log
python RunMModesHouse.py 2> HouseMModesError.log
