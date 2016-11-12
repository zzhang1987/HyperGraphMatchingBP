#!/bin/sh
#SBATCH --job-name=demo_modes    # Job name for easy identification in queue
#SBATCH --time=48:00:00      # Time requested (HH:MM:SS, or DAYS-HH:MM:SS)
#SBATCH --nodes=1           # Request a single node only
#SBATCH --ntasks=4         # Job will consist of only one process
#SBATCH --cpus-per-task=1   # Number of CPUs per process (for multithreaded programs)
#SBATCH --mem=32GB           # Memory requested per nod#!/bin/sh
LD_PRELOAD=/usr/lib64/libstdc++.so.6 python demoCarHyper.py > run.log 2> run.err
