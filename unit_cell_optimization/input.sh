#!/bin/bash

# After making necessary changes (ie calculation type), run as executable: ./input.sh

############################################## MANUAL INPUT ##############################################

  # state calculation type from: kpoints, ecutwfc, ecutrho, lattice, geometry, relax, vc-relax
  calculation='kpoints'

# ============================================================= #
#### IF CONVERGENCE PARAMETERS HAVE BEEN DETERMINED ALREADY ####
# ============================================================= #
  # kpoint grid (if already determined from convergence). In format "2 2 2".
  # 3 numbers with a space in between, all in quotes
  kpoint_grid=
  # ecutwfc (if already determined from convergence)
  ecutwfc=
  # ecutrho (if already determined from convergence)
  ecutrho=

# ============================================================= #
        #### INFORMATION ON THE ELEMENTS PRESENT ####
# ============================================================= #
  # system we're studying. Ex: Al or AlCe or Al11Ce3 (include stoichiometry)
  system="Al"
  # list of unique elements. Must be in format like "Al Ce Cu" with a space between, all within quotations
  elements="Al"

  # directory on cluster where you want calculations to run/results to be saved
  compute_directory="/path/to/QUEENS"
  # path to folder containing all of our scripts for generating/running calculations
  script_directory="/path/to/DFT_scripts"
  # path to folder containing index files (atomic_weights.txt, list_of_pseudopotentials.txt etc.)
  index_directory="/path/to/DFT_index_files"
  # pseudopotential folder
  pseudo_directory="/path/to/PSEUDOPOTENTIALS"
  # where is .txt file with list of python wheels for easy installation when calling python
  python_wheels="/path/to/python_wheels.txt"

# ============================================================= #
                #### LATTICE INFORMATION ####
# ============================================================= #
  # FCC, BCC, HCP, Orthorhombic
  symmetry="FCC"
  # Primitive or Unit
  cell_type="Primitive"

  # If you have a .cif file, we can extract data using cif2cell output. Set to true. Otherwise, false.
  # NOTE: cif_file has to be named after the system. Ex: Al.cif, AlFe.cif, Al11Ce3.cif etc.
  # whatever the you have for "system" variable above is the name we need for .cif file
  # NOTE: IF NOT USING A CIF FILE, YOU MUST ENTER ATOMIC POSITIONS MANUALLY IN SECTION BELOW
  # true or false
  using_cif_file=false

# ============================================================= #
            #### LATTICE INFORMATION (NO .CIF) ####
# ============================================================= #
  # number of atoms in cell. If using .cif, leave blank
  natoms=1

  # INITIAL lattice parameter (pre-optimization), in angstroms. If using a .cif, leave blank
  starting_lattice=3.9

  # if orthorhombic or hcp and NO .cif file, need to manually input b/a or c/a ratios.
  # For orthorhombic, start with: b_a=1.633, c_a=1.73205081 <-- sqrt(3)
  # for HCP, start with: b_a=1, c_a=1.633
  # doesn't REALLY matter what these values are, because we'll run vc-relax which will vary these for us. We just want
  # to reduce computational expense/time of vc-relax by getting as close as possible through manual calculations
  b_a=1
  c_a=1

  # If no .cif file, enter atomic positions manually
  # here, just modify atomic positions based on the symmetry/elements you're using
  # remember, for primitive cells put 1 atom at origin
  # for primitive hcp, 1 at origin and 1 at (2/3, 1/3, 1/2)
  manual_atomic_positions="ATOMIC_POSITIONS {crystal}
Al 0.000000000000000 0.000000000000000 0.000000000000000"

# ============================================================= #
    #### INFORMATION FOR HOW TO HANDLE THE GIVEN SYSTEM ####
# ============================================================= #
  # PAW or USPP
  pseudopotential_type="PAW"
  # scalar ("Scalar") or full relativistic ("Full")
  pseudo_rel_type="Scalar"
  # PBE, PBESOL, LDA, or R2SCAN
  functional_type="PBE"

# ============================================================= #
                  #### KPOINT PARAMETERS ####
# ============================================================= #
  # set to true to define manual_kpoint_array instead of calculating it.
  # true or false. Generally do false unless doing slab calculation or custom cell shape
  use_manual_kpoints=false

  # kpoint array. For HCP and orthorhombic where a=b, kpoints are the same along a/b, and different over c
  # but if a =/ b =/ c, then Kpoints will be different in each direction. Hav eto make them roughly proportional to
  # b/a and c/a ratios. Script calculates kpoints based on these ratios.
  # Can't automate for all possibilities though, so allow for option to enter manually
  # If converging kpoints for orthorhombic , enter the kpoint arrays you want to sample over
  # for example, if b=c you could do ("2 1 1" "4 2 2" "6 4 4" "8 6 6" "10 6 6" "12 8 8") etc.
  manual_kpoint_array=

  # true for Γ-centered grid (default); false for shifted grid
  gamma_centered=true

############################################ END MANUAL INPUT ############################################

# run the calculation
case "$calculation" in
  "kpoints"|"ecutwfc"|"ecutrho"|"lattice"|"geometry")
    source scf_calculations.sh
  ;;
  "relax"|"vc-relax")
    source relax_calculations.sh
  ;;
esac
