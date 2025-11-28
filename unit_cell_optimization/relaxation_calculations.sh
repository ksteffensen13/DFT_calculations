#!/bin/bash

# first, calculate elastic constants
# take optimized unit cell, run scf, make sure there's no forces acting on the atoms, then run 'relax' for each deformation
# from the deformations, extract the stress tensors (invert signs) and convert to GPa
# note: inverse sign of stress tensor because the code reports external stress (positive when pointing in, negative when pointing out)
# but we want internal stress (the opposite)
# then convert to Voigt notation


# ============================================
# DEFINE SYSTEM NAME
# ============================================
file_naming="${system}_${symmetry}"
system_directory="${compute_directory}/${file_naming}"
# if the compute directory doesn't exist, then create it
mkdir -p "${system_directory}"
mkdir -p "${system_directory}"/optimized_unit_analyses


# function for running matrix operations to calculate new cell_parameters
matrix_operations() {
  # grab deformation type, as well as strain amounts
  local deformation=$1
  echo "deformation type is $deformation"
  e1=$2
  # convert to 15 decimals for cell_parameter section later
  local e1=$(printf "%.15f" "$e1")
  echo "strain1 is $e1"
  e2=$3
  local e2=$(printf "%.15f" "$e2")
  echo "strain2 is $e2"
  # declare cell parameter 3x3 matrix elements to be global variables, so they can be called later
  declare -g cell_param11
  declare -g cell_param12
  declare -g cell_param13
  declare -g cell_param21
  declare -g cell_param22
  declare -g cell_param23
  declare -g cell_param31
  declare -g cell_param32
  declare -g cell_param33


  # Original lattice vectors (root tensor components)
  # ax, ay, az etc. must be set globally before calling this function
  local vectors=("a" "b" "c")
  local i j

  # Deformation matrix elements (initialize with zeros)
  local F=(0 0 0 0 0 0 0 0 0)  # F11 F12 F13 F21 F22 F23 F31 F32 F33
  # so, in QE the cell_parameter 3x3 matrix is the transposed version of the root tensor
  # so, if 3 lattice vectors a=(ax, ay, az) b=(bx, by, bz) c=(cx, cy, cz),
  # cell parameters:    Root tensor:
  #    ax ay az           ax bx cx
  #    bx by bz           ay by cy
  #    cx cy cz           az bz cz

  # when applying deformation matrix, it goes to the root tensor not the cell_parameter matrix
  # so if deformation matrix F:
  # 1.03 0 0 * ax = 1.03ax+0ay+0az = 1.03*ax
  # 0.00 1 0   ay = 0.00ax+1ay+0az = 1.00*ay
  # 0.00 0 1   az = 0.00ax+0ay+1az = 1.00*az
  # repeat for vectors b and c,
  # 1.03*ax 1.03*bx 1.03*cx
  # 1.00*ay 1.00*by 1.00*cy
  # 1.00*az 1.00*bz 1.00*cz
  # then transpose into cell parameters:
  # 1.03*ax 1.00*ay 1.00*az
  # 1.03*bx 1.00*by 1.00*bz
  # 1.03*cx 1.00*cy 1.00*cz

  # for each deformation possiblity, define the components of deformation matrix in form:
  # F11 F12 F13
  # F21 F22 F23
  # F31 F32 F33
  local zero=$(printf "%.15f" "0")
  local one=$(printf "%.15f" "1")
  local def=$(printf "%.15f" "$(bc -l <<< "1 + $e1")")

  # Assign deformation matrix based on type
  # Deformation types:
  # F1=         F2=       F3=         F4=      F5=    F6=
  # 1+e1 0 0   1  0   0   1 0  0     1 e2 0   1 0 e2  1 0 0
  #  0   1 0   0 1+e1 0   0 1  0     0 1  0   0 1 0   0 1 e2
  #  0   0 1   0  0   1   0 0 1+e1   0 0  1   0 0 1   0 0 1
  case "$deformation" in
      "deformation1") F=($def $zero $zero $zero $one $zero $zero $zero $one) ;;
      "deformation2") F=($one $zero $zero $zero $def $zero $zero $zero $one) ;;
      "deformation3") F=($one $zero $zero $zero $one $zero $zero $zero $def) ;;
      "deformation4") F=($one $e2 $zero $zero $one $zero $zero $zero $one) ;;
      "deformation5") F=($one $zero $e2 $zero $one $zero $zero $zero $one) ;;
      "deformation6") F=($one $zero $zero $zero $one $e2 $zero $zero $one) ;;
  esac
  echo "deformation matrix is
  ${F[0]} ${F[1]} ${F[2]}
  ${F[3]} ${F[4]} ${F[5]}
  ${F[6]} ${F[7]} ${F[8]}" >> matrix_operations.txt
  # now calculate the new cell parameters based on lattice vectors a, b, c
  # new root tensor will be:
  # new_ax new_bx new_cx
  # new_ay new_by new_cy
  # new_az new_bz new_cz
  # use a loop to calculate new values for each

  # Compute new root tensor (F * original vectors)
  echo "calculating deformed lattice vectors"
  # loop over a, b, c to calculate ax/ay/az, then bx/by/bz, then cx/cy/cz
  for vec in "${vectors[@]}"; do
      local x y z
      eval x=\$${vec}x
      eval y=\$${vec}y
      eval z=\$${vec}z
      echo "calculating new vectors ${vec}x ${vec}y ${vec}z"
      # When calculate, make sure to convert to 15 decimal points for consistent CELL_PARAMETER formatting
      eval new_${vec}x=$(printf "%.15f" "$(bc -l <<< "(${F[0]} * $x) + (${F[1]} * $y) + (${F[2]} * $z)")")
      eval new_${vec}y=$(printf "%.15f" "$(bc -l <<< "(${F[3]} * $x) + (${F[4]} * $y) + (${F[5]} * $z)")")
      eval new_${vec}z=$(printf "%.15f" "$(bc -l <<< "(${F[6]} * $x) + (${F[7]} * $y) + (${F[8]} * $z)")")
  done

# transposing new root tensor:
# new_ax new_ay new_az
# new_bx new_by new_bz
# new_cx new_cy new_cz

# which means new cell_parameters will be:
# cell_param11 cell_param12 cell_param13
# cell_param21 cell_param22 cell_param23
# cell_param31 cell_param32 cell_param33

cell_param11=${new_ax}
cell_param12=${new_ay}
cell_param13=${new_az}

cell_param21=${new_bx}
cell_param22=${new_by}
cell_param23=${new_bz}

cell_param31=${new_cx}
cell_param32=${new_cy}
cell_param33=${new_cz}

echo "old parameters were
${ax} ${ay} ${az}
${bx} ${by} ${bz}
${cx} ${cy} ${cz}" >> matrix_operations.txt

echo "Cell Parameters will be
${cell_param11} ${cell_param12} ${cell_param13}
${cell_param21} ${cell_param22} ${cell_param23}
${cell_param31} ${cell_param32} ${cell_param33} " >> matrix_operations.txt

}


deform_cell_for_elastic_constants() {

declare -g ax 
declare -g ay 
declare -g az
declare -g bx 
declare -g by 
declare -g bz
declare -g cx 
declare -g cy 
declare -g cz


echo "making elastic_constants directory"
mkdir -p "${system_directory}"/optimized_unit_analyses/elastic_constants
cd "${system_directory}"/optimized_unit_analyses/elastic_constants

# first, pull most of the .in information from given file
# in this case, use .out from vc-relax since it is the most optimized
output_file_to_use="${system_directory}/vc-relax/${file_naming}.vc-relax.out"
input_file_to_use="${system_directory}/vc-relax/${file_naming}.vc-relax.in"

echo "extracting parameters for .in file"
# for each of the below lines, using tail -n1 gets us the last instance
# (aka the optimized cell, rather than the input parameters in case vc-relax made changes)
# grab number of atoms
nat=$(grep "number of atoms/cell" ${output_file_to_use} | tail -n1 | cut -d= -f2 | tr -d ' ')
# number of elements
ntyp=$(grep "number of atomic types" ${output_file_to_use} | tail -n1 | cut -d= -f2 | tr -d ' ')
# ecutwfc
ecutwfc=$(grep "kinetic-energy cutoff" ${output_file_to_use} | tail -n1 | cut -d= -f2 | tr -d ' ' | cut -d. -f1)
# ecutrho
ecutrho=$(grep "charge density cutoff" ${output_file_to_use} | tail -n1 | cut -d= -f2 | tr -d ' ' | cut -d. -f1)
# lattice parameters
# NOTE: Don't need to re-calculate, QE will do it for us
# we supply CELL_PARAMETERS, which define the cell dimensions for the deformed cells
# CELL_PARAMETERS are in terms of lattice parameter.
# So if we wanted to increase x by 20%, we could calculate new lattice, or just put cell parameter of 1.2 in x
# which will tell QE that the cell will be 20% longer than the previous lattice parameter in that direction
output_celldm1=$(grep "CELL_PARAMETERS (alat=" ${output_file_to_use} | tail -n1 | sed -E 's/.*alat= *([0-9.]+).*/\1/')
lattice_parameter="celldm(1)=${output_celldm1},"

# each element of the lattice vectors/cell parameters can be found from the section after "End final coordinates"
# where it lists lattice vectors a(1) a(2) and a(3), which we'll refer to as a, b, and c respectively
read ax ay az bx by bz cx cy cz <<< $(awk '/End final coordinates/{flag=1; next} flag && /a\(1\)/ {a1=$4; a2=$5; a3=$6}
     flag && /a\(2\)/ {b1=$4; b2=$5; b3=$6}
     flag && /a\(3\)/ {c1=$4; c2=$5; c3=$6}
     END {print a1,a2,a3,b1,b2,b3,c1,c2,c3}' ${output_file_to_use})

# then convert lattice vectors to 15 decimals to fill in CELL_PARAMETERS section later on
ax=$(printf "%.15f" "$ax")
ay=$(printf "%.15f" "$ay")
az=$(printf "%.15f" "$az")
bx=$(printf "%.15f" "$bx")
by=$(printf "%.15f" "$by")
bz=$(printf "%.15f" "$bz")
cx=$(printf "%.15f" "$cx")
cy=$(printf "%.15f" "$cy")
cz=$(printf "%.15f" "$cz")

# grab atomic species from input file (since this won't change at all)
atomic_species=$(awk '
  /^ATOMIC_SPECIES/ {capture=1; print; next}
  capture {
      # stop if we reach a line starting with an uppercase QE section
      if ($1 ~ /^[A-Z_]+$/) { exit }
      print
  }
' ${input_file_to_use})

# now grab final atomic positions from output file (vc-relax) in case the atoms have shifted
# this will be between the lines "Begin final coordinates" and "End final coordinates"
# this command is setup to capture the lines right after it, regardless of number of spaces between columns/decimals
atomic_positions=$(
awk '
  /Begin final coordinates/ { in_final=1; next }
  /End final coordinates/   { in_final=0; capture=0 }

  # Detect the start of the block
  in_final && /^ATOMIC_POSITIONS/ {
    capture=1
    print
    next
  }

  # Match an atomic line:
  #   element_symbol   num   num   num
  # The numeric regex matches integers, floats, scientific notation.
  capture &&
  $1 ~ /^[A-Za-z][A-Za-z0-9]*$/ &&
  $(NF-2) ~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/ &&
  $(NF-1) ~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/ &&
  $NF   ~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/ {
    print
    next
  }

  # Any other line ends capture
  capture { capture=0 }
' ${output_file_to_use}
)

# from input file, grab the kpoint grid since this won't change
kpoints=$(awk '
  /^K_POINTS/ {capture=1; print; next}
  capture {
      if ($0 ~ /^$/) exit
      print
  }
' ${input_file_to_use}
)

# now we want to deform the cell
# we have 6 deformations to perform,
# For each deformation, we build a 3x3 deformation matrix and apply it to the ROOT TENSOR of our cell
# NOTE: in QE, cell parameters section of .in is 3x3 matrix, which is transposed root tensor
# cell parameters:    Root tensor:
#    ax ay az           ax bx cx
#    bx by bz           ay by cy
#    cx cy cz           az bz cz

# then, create list of deformation types so the script can determine how to deform the cell
deformation_list=("deformation1" "deformation2" "deformation3" "deformation4" "deformation5" "deformation6")
# and a list of different combinations of e1 and e2
# Do multiple values, and opposing signs, so that we have multiple values to fit later for more accuracy
strain_list=("-0.01 -0.03" "-0.005 -0.015" "0.005 0.015" "0.01 0.03")

# now, for each pair of strains, run through the 6 deformation types and generate a .in file
for strain in "${strain_list[@]}"; do
  read -r strain1 strain2 <<< "$strain"
  for def_type in "${deformation_list[@]}"; do
    # run the function to calculate cell parameters by applying deformation to root tensor
    matrix_operations ${def_type} ${strain1} ${strain2}

cat > "${file_naming}".tensors_strain${strain1}_${def_type}.in << EOF2
&CONTROL
calculation='relax',
outdir='.',
prefix='${file_naming}_tensors_strain${strain1}_${def_type}',
pseudo_dir='${pseudo_directory}',
verbosity='low',
tprnfor=.true.,
tstress=.true.,
/

&SYSTEM
ibrav = 0,
${lattice_parameter}
nat=${nat},
ntyp=${ntyp},
ecutwfc=${ecutwfc},
ecutrho=${ecutrho},
occupations='smearing',
smearing='mp',
degauss=0.03,
/

&ELECTRONS
conv_thr=1d-08,
mixing_beta=0.7d0,
diagonalization='cg',
/

&IONS
ion_dynamics='bfgs',
/

CELL_PARAMETERS {alat}
${cell_param11} ${cell_param12} ${cell_param13}
${cell_param21} ${cell_param22} ${cell_param23}
${cell_param31} ${cell_param32} ${cell_param33}

${atomic_species}

${atomic_positions}

${kpoints}

EOF2
  done
done

cat > "${file_naming}".run_tensors.sh << EOF1
#!/bin/bash
#SBATCH --time=0-02:00:00
#SBATCH --account=def-diak01
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000
#SBATCH --output=${file_naming}_tensor_deformation.out
#SBATCH --job-name=${file_naming}_tensor_deformation
#SBATCH --mail-user=karl.steffensen@queensu.ca
#SBATCH --mail-type=ALL

module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 intel/2023.2.1 gcccore/.12.3
module load quantumespresso/7.3.1

directory="${system_directory}/optimized_unit_analyses/elastic_constants"

cd \${directory}

for input_file in \$(ls \${directory}/*.in)
do
    srun pw.x -in "\$input_file" > "\${input_file%.*}.out"
done
EOF1

sbatch "${file_naming}".run_tensors.sh

}


build_tensors() {
# To calculate elastic constants, we need stress and strain tensors
# stress tensors just come from .out file of each deformation calculation
# strain tensors (E) are calculated via Green-lagrange: E=0.5*(F.T*F - I)
# Here, F is deformation matrix, F.T is transposed deformation matrix, and I is unit matrix (100 010 001)
# since we have 4 combinations of strains, and 6 deformation types, we will have 24 of each tensor
# these will be used in python to fit to a final matrix

analysis_directory="${system_directory}/optimized_unit_analyses/elastic_constants"
cd "${analysis_directory}"


# for the list of deformations
deformation_list=("deformation1" "deformation2" "deformation3" "deformation4" "deformation5" "deformation6")
# and strains
strain_list=("-0.01 -0.03" "-0.005 -0.015" "0.005 0.015" "0.01 0.03")

# set counters for the loops
a=1
# then loop over all combinations of strains/deformation type
for strain in "${strain_list[@]}"; do
  read -r strain1 strain2 <<< "$strain"
  b=1
  for def_type in "${deformation_list[@]}"; do
    # identify the .out file to use
    output_file_name="${file_naming}.tensors_strain${strain1}_${def_type}.out"
    output_file_to_use="${analysis_directory}/${output_file_name}"

    # extract the stress tensor from the lines after "total stress"
    stress_tensor=($(awk '
    /total[[:space:]]+stress/ {
        getline
        printf "%.10f %.10f %.10f ", $1, $2, $3
        getline
        printf "%.10f %.10f %.10f ", $1, $2, $3
        getline
        printf "%.10f %.10f %.10f\n", $1, $2, $3
    }' "${output_file_to_use}"))

    # save each stress tensor in the larger array
    # each entry will be categorized by strain and deformation type
    # ie stress tensor for strain of 0.03 and deformation4 will be stored under stress_tensors[0.03,deformation4]
    # so you can look into the array and tell it to give you the values stored under "0.03,deformation4"
    # makes it easy to pull the stress tensor we're interested in

    # Now, quantum espresso lists stress tensor as external stress, but we want internal stress
    # so we need to reverse the signs of all stress values
    # while we're at it, convert from Ry/bohr^3 to GPa (1 Ry/bohr^3=14710.5076 GPa)
    for i in "${!stress_tensor[@]}"; do
      stress_tensor[$i]="$(awk -v x="${stress_tensor[$i]}" 'BEGIN{print -x*14710.5076}')"
    done

    # Store in .txt file
    echo "${stress_tensor[*]}" >> "stress_tensor_strain${a}_def${b}.txt"

    # now do it for strain tensor.
    # Strain tensors will be calculated in Python because it's easier
    # So we will just build deformation matrices F, which will feed into Python to find F.T and calculate E

    # pull out the current strain
    local e1=$(printf "%.15f" "$strain1")
    local e2=$(printf "%.15f" "$strain2")

    # build deformation matrix F
    local F=(0 0 0 0 0 0 0 0 0)

    local zero=$(printf "%.15f" "0")
    local one=$(printf "%.15f" "1")
    local def=$(printf "%.15f" "$(bc -l <<< "1 + $e1")")

    # Assign deformation matrix based on type
    # in format F=F11 F12 F13 F21 F22 F23 F31 F32 F33 from F=
    # F11 F12 F13
    # F21 F22 F23
    # F31 F32 F33

    # Deformation types:
    # F1=       F2=       F3=         F4=       F5=     F6=
    # 1+e1 0 0   1  0   0   1 0  0     1 e2 0   1 0 e2  1 0 0
    #  0   1 0   0 1+e1 0   0 1  0     0 1  0   0 1 0   0 1 e2
    #  0   0 1   0  0   1   0 0 1+e1   0 0  1   0 0 1   0 0 1
    case "$def_type" in
        "deformation1") F=($def $zero $zero $zero $one $zero $zero $zero $one) ;;
        "deformation2") F=($one $zero $zero $zero $def $zero $zero $zero $one) ;;
        "deformation3") F=($one $zero $zero $zero $one $zero $zero $zero $def) ;;
        "deformation4") F=($one $e2 $zero $zero $one $zero $zero $zero $one) ;;
        "deformation5") F=($one $zero $e2 $zero $one $zero $zero $zero $one) ;;
        "deformation6") F=($one $zero $zero $zero $one $e2 $zero $zero $one) ;;
    esac
    # store in .txt file
    echo "${F[*]}" >> "deformation_matrix_strain${a}_def${b}.txt"

    b=$(( b + 1 ))
  done
  a=$(( a + 1 ))
done

}

calculate_elastic_constants() {

  # now that we have our stress/strain tensors, run python script to convert them to 3x3 matrices, modify into
  # Voigt format, then calculate the elastic constant 6x6 matrix via matrix operations
  # run cif2cell
  module load StdEnv/2023 python/3.13.2
  ENVDIR=/tmp/$RANDOM
  virtualenv --no-download $ENVDIR
  source $ENVDIR/bin/activate
  pip install --no-index --upgrade pip
  pip install -r ${python_wheels}
  python ${script_directory}/calculate_tensors.py "${analysis_directory}"
  deactivate
  rm -rf $ENVDIR

}



case "$analysis" in
  "elastic_constant_deformation")
    deform_cell_for_elastic_constants ;;
  "elastic_constant_calculation")
    build_tensors
    calculate_elastic_constants ;;
esac

