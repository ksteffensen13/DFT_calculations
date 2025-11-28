#!/bin/bash

# for running optimized_unit_cell_analyses.sh

############################################## MANUAL INPUT ##############################################
  # analysis you want to run: "elastic_constant_deformation", "elastic_constant_calculation"
  analysis="elastic_constant_deformation"
  system="Al"
  symmetry="HCP"

  # directory on cluster where you want calculations to run/results to be saved
  compute_directory="/home/ksteffen/scratch/QUEENS"
  # path to folder containing all of our scripts for generating/running calculations
  script_directory="/home/ksteffen/scratch/QUEENS/DFT_scripts"
  # path to folder containing index files (atomic_weights.txt, list_of_pseudopotentials.txt etc.)
  index_directory="/home/ksteffen/scratch/QUEENS/DFT_index_files"
  # pseudopotential folder
  pseudo_directory="/home/ksteffen/scratch/QUEENS/PSEUDOPOTENTIALS"
  # where is .txt file with list of python wheels for easy installation when calling python
  python_wheels="/home/ksteffen/python_wheels.txt"

source optimized_unit_cell_analyses.sh
