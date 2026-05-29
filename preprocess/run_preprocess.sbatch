#!/bin/bash
#SBATCH --job-name=wrf_preprocess   # Name showing up in squeue
#SBATCH --account=strong-kp    	    # Court's node (kingspeak11)
#SBATCH --partition=strong-kp
#SBATCH --nodes=1                   # Single node is best for single script multi-threading
#SBATCH --ntasks=16         
#SBATCH --output=wrf_pre_%j.out     # Standard output log file (%j dynamically inserts Job ID)
#SBATCH --error=wrf_pre_%j.err      # Error log file

# Exit script immediately if any command fails
set -e

# Clear system modules to avoid library conflicts
module purge

# Initialize shell support for your Miniforge/Conda installation
source /uufs/chpc.utah.edu/common/home/strong-group7/sydney/miniforge3_envs/karakoram

# Activate your exact working environment
conda activate karakoram

# Run the python script natively on the allocated compute node
python preprocess.py
