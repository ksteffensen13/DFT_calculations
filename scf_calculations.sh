#!/bin/bash

# ============================================
# DEFINE SYSTEM NAME
# ============================================
file_naming="${system}_${symmetry}"
system_directory="${compute_directory}/${file_naming}"
# if the compute directory doesn't exist, then create it
mkdir -p "${system_directory}"

# ============================================
# LOOK UP ATOMIC WEIGHTS AND PP INFORMATION
# ============================================

set_pseudo_atomic_info() {
  declare -g max_ecut
  declare -g max_ecutrho
  declare -g max_ecutrho_ratio
  declare -gA pseudopotential_file
  declare -gA atomic_weight

  # split element list into array
  read -r -a element_array <<< "$elements"

  # for each element, find the matching line, and print column 2 (weight)
  # store this value in the array, under an entry matching the element
  for el in "${element_array[@]}"; do
    weight=$(awk -v el="$el" '$1 == el {print $2}' "${index_directory}/atomic_weights.txt")
    atomic_weight["$el"]="$weight"

    # then for everything in the element array, look in the .txt file, find the lines where
    # entries match our inputs (element name, PP type, XCF type, scalar vs full relativistic)
    # for the line where all entries match, pull out the PP file name from the last column
    pseudo_file=$(awk -v el="$el" -v pt="$pseudopotential_type" -v rel="$pseudo_rel_type" -v xc="$functional_type" '
    NR > 1 && toupper($1)==toupper(el) && toupper($2)==toupper(pt) && \
              toupper($3)==toupper(rel) && toupper($4)==toupper(xc) { print $5; exit }
  ' "${index_directory}/list_of_pseudopotentials.txt")

    # save the PP file name in the pseudopotential array, under the element name
    # ex: AL_USPP.UPF would be saved in the array under "Al", so that if say "look in the array and tell me what you find
    # saved under element Al" it will pull up "AL_USPP.UPF". Want this for all elements in the system, hence the loop
    pseudopotential_file["$el"]="$pseudo_file"
  done


  # DETERMINE MAX ECUT AND ECUTRHO FROM PP FILES

  # remember, each PP file gives values for ecutwfc and ecutrho. The minimum values we should use are the LARGEST of all
  # recommended values for all PP's we're using. For example, if ecutwfc is 40 for Al and 60 for Ce, we want to use
  # ecutwfc of AT LEAST 60. Same process for ecutrho

  # first, set the variables for maximum to 0
  max_ecut=0
  max_ecutrho=0

  # Then, for each element pull up the entry in pseudopotential array that matches the current element
  # Then for that file name, look inside and pull out the ecut and ecutrho values
  # store within max_ecut and max_ecutrho
  # repeat for all elements, and compare to the already stored values
  # if greater than the currently stored max values, update with current values. If not, continue with stored values
  for el in "${element_array[@]}"; do
      pseudo_file="${pseudopotential_file[$el]}"

      echo "Reading pseudopotential for $el → $pseudo_file"

      # Extract ecutwfc
      ecut=$(awk '/wavefunctions/ {
          gsub(/\.$/, "", $6);
          print $6
      }' "${pseudo_directory}/${pseudo_file}")

      # Extract ecutrho
      ecutrho_PP=$(awk '/charge density/ {
          gsub(/\.$/, "", $7);
          print $7
      }' "${pseudo_directory}/${pseudo_file}")

      # If greater than current max_ecut, update the values to store the new largest values
      (( $(echo "$ecut > $max_ecut" | bc -l) )) && max_ecut=$ecut
      (( $(echo "$ecutrho_PP > $max_ecutrho" | bc -l) )) && max_ecutrho=$ecutrho_PP
  done

  # based on the max values, calculate the ratio of ecutrho to ecutwfc and round to nearest whole number
  ecutrho_ratio=$(awk -v x="$max_ecut" -v y="$max_ecutrho" 'BEGIN { printf "%.0f", y/x }')

  # based on the largest values from pseudopotentials
  # list them to ensure our final convergence choice is above this threshold
  echo "minimum ecutwfc from PP: $max_ecut"
  echo "minimum ecutrho from PP: $max_ecutrho"
  echo "ecut_rho ratio from PP: $ecutrho_ratio"
}

# ============================================
# GENERATE ALL .in INFORMATION
# ============================================

# KPOINTS GRID TYPE
# Γ-centered or non-centered?
if [ "$gamma_centered" = true ]; then
  kshift="0 0 0"
else
  kshift="1 1 1"
fi

# define function to set variables that can be used to fill in .in file
# based on symmetry, will have different ibrav, format for lattice parameter etc.
# will also have different number of atoms based on symmetry AND whether we want primitive or unit cell
# if ibrav=0, we also have to define cell parameters in terms of lattice vectors (alat)
set_input_parameters_no_cif() {
  local celltype=$1
  local sym=$2

  echo "generating .in parameters from symmetry"
  case "$sym" in
    "FCC")
      ibrav=2
      nat="nat = $([ "$celltype" = "Unit" ] && echo 4 || echo 1)"
      ;;
    "BCC")
      ibrav=3
      nat="nat = $([ "$celltype" = "Unit" ] && echo 2 || echo 1)"
      ;;
    "HCP")
      ibrav=4
      nat="nat = $([ "$celltype" = "Unit" ] && echo 6 || echo 2)"
      ;;
    "Orthorhombic")
      ibrav=0
      nat="nat = ${natoms}"
      alat_b=$(printf "%.15f" "$b_a")
      alat_c=$(printf "%.15f" "$c_a")
      cell_parameters="CELL_PARAMETERS {alat}
  1.000000000000000   0.000000000000000   0.000000000000000
  0.000000000000000   ${alat_b}   0.000000000000000
  0.000000000000000   0.000000000000000   ${alat_c}"
      ;;
  esac

  # calculate number of unique elements
  # split into positional parameters
  set -- $elements
  # count them
  ntyp="ntyp = $#"

  echo "building ATOMIC_POSITIONS"
  # build atomic_positions block from manual entry
  atomic_positions="${manual_atomic_positions}"

  echo "building ATOMIC_SPECIES"
  # Build ATOMIC_SPECIES section
  atomic_species="ATOMIC_SPECIES"
for el in "${element_array[@]}"; do
  w="${atomic_weight[$el]}"
  pseudo="${pseudopotential_file[$el]}"
  atomic_species+="
${el} ${w} ${pseudo}"
done
}

# if working from .cif file, extract the information
set_input_parameters_cif() {
  local celltype=$1
  local sym=$2

  echo "generating .in parameters from .cif file"

  # cif2cell output file to extract positions from
  cif_file="${compute_directory}/DFT_cif_files/${system}.cif"
  cif2cell_output="${system_directory}/${system}_converted_cif.in"

  echo "converting .cif file to .in"
  # run cif2cell
  module load StdEnv/2023 python/3.13.2
  ENVDIR=/tmp/$RANDOM
  virtualenv --no-download $ENVDIR
  source $ENVDIR/bin/activate
  pip install --no-index --upgrade pip
  pip install -r /home/ksteffen/python_wheels.txt
  cif2cell -f "${cif_file}" -p quantum-espresso -o "${cif2cell_output}"
  deactivate
  rm -rf $ENVDIR

  echo "extracting from converted .cif"
  case "$sym" in
    "FCC")
      ibrav=2
      nat=$(awk '/^[[:space:]]*nat[[:space:]]*=/ {print; exit}' "$cif2cell_output")
      ntyp=$(awk '/^[[:space:]]*ntyp[[:space:]]*=/ {print; exit}' "$cif2cell_output")
      ;;
    "BCC")
      ibrav=3
      nat=$(awk '/^[[:space:]]*nat[[:space:]]*=/ {print; exit}' "$cif2cell_output")
      ntyp=$(awk '/^[[:space:]]*ntyp[[:space:]]*=/ {print; exit}' "$cif2cell_output")
      ;;
    "HCP")
      ibrav=4
      nat=$(awk '/^[[:space:]]*nat[[:space:]]*=/ {print; exit}' "$cif2cell_output")
      ntyp=$(awk '/^[[:space:]]*ntyp[[:space:]]*=/ {print; exit}' "$cif2cell_output")
      ;;
    "Orthorhombic")
      ibrav=0
      nat=$(awk '/^[[:space:]]*nat[[:space:]]*=/ {print; exit}' "$cif2cell_output")
      ntyp=$(awk '/^[[:space:]]*ntyp[[:space:]]*=/ {print; exit}' "$cif2cell_output")
      cell_parameters=$(awk '/^CELL_PARAMETERS/{flag=1; print; next} /^ATOMIC_SPECIES/{flag=0} flag' "$cif2cell_output")
      ;;
  esac

  echo "building ATOMIC_SPECIES"
  # build ATOMIC_SPECIES block
  atomic_species=$(awk '/^ATOMIC_SPECIES/{flag=1; print; next} /^ATOMIC_POSITIONS/{flag=0} flag' "$cif2cell_output")

  # Rebuild ATOMIC_SPECIES block with correct pseudopotential filenames (cif2cell just has element_PSEUDO)
  updated_atomic_species=""
  while IFS= read -r line; do
    if [[ "$line" =~ ^[[:space:]]*([A-Za-z]+)[[:space:]]+([0-9.]+)[[:space:]]+ ]]; then
      el="${BASH_REMATCH[1]}"
      w="${BASH_REMATCH[2]}"
      pseudo="${pseudopotential_file[$el]}"
      if [[ -n "$pseudo" ]]; then
        updated_atomic_species+="${el} ${w} ${pseudo}\n"
      else
        updated_atomic_species+="${line}\n"
      fi
    else
      # Keep non-element lines (e.g., the "ATOMIC_SPECIES" header)
      updated_atomic_species+="${line}\n"
    fi
  done <<< "$atomic_species"

  # Convert \n to actual newlines
  atomic_species="$(echo -e "$updated_atomic_species")"

  echo "building ATOMIC_POSITIONS"
  atomic_positions=$(awk '/^ATOMIC_POSITIONS/{flag=1; print; next} /^$/{flag=0} flag' "$cif2cell_output")
}

set_lattice_parameter_no_cif() {
  local calc=$1
  local sym=$2
  declare -g lattice_parameter

  case "$calc" in
  "kpoints"|"ecutwfc"|"ecutrho")
    case "$sym" in
    "FCC"|"BCC"|"Orthorhombic")
      lattice_parameter="A = ${starting_lattice}"
      ;;
    "HCP")
      celldm1=$(bc -l <<< "${starting_lattice} / 0.529177525830478")
      lattice_parameter="celldm(1) = ${celldm1}
celldm(3) = ${c_a}"
      ;;
    esac
    ;;
  "lattice"|"geometry")
    lattice_parameter='${lattice_parameter}'
    ;;
  esac
}

set_lattice_parameter_cif() {
  local calc=$1
  local sym=$2
  declare -g lattice_parameter

  # cif2cell output file to extract positions from
  cif_file="${compute_directory}/DFT_cif_files/${system}.cif"
  cif2cell_output="${system_directory}/${system}_converted_cif.in"

  case "$calc" in
  "kpoints"|"ecutwfc"|"ecutrho")
    case "$sym" in
    "FCC"|"BCC"|"Orthorhombic")
      a=$(awk '/_cell_length_a/{print $2}' "$cif_file")
      lattice_parameter="A=${a}"
      b=$(awk '/_cell_length_b/{print $2}' "$cif_file")
      c=$(awk '/_cell_length_c/{print $2}' "$cif_file")
      b_a=$(awk -v a="$lattice_parameter" -v b="$b" 'BEGIN { printf "%.0f", b/a }')
      c_a=$(awk -v a="$lattice_parameter" -v c="$c" 'BEGIN { printf "%.0f", c/a }')
      ;;
    "HCP")
      starting_lattice=$(awk '/_cell_length_a/{print $2}' "$cif_file")
      celldm1=$(bc -l <<< "${starting_lattice} / 0.529177525830478")
      b=$(awk '/_cell_length_b/{print $2}' "$cif_file")
      c=$(awk '/_cell_length_c/{print $2}' "$cif_file")
      b_a=$(awk -v a="$lattice_parameter" -v b="$b" 'BEGIN { printf "%.0f", b/a }')
      c_a=$(awk -v a="$lattice_parameter" -v c="$c" 'BEGIN { printf "%.0f", c/a }')

      lattice_parameter="celldm(1) = ${celldm1}
celldm(3) = ${c_a}"
      ;;
    esac
    ;;
  "lattice"|"geometry")
    lattice_parameter='${lattice_parameter}'
    ;;
  esac
}


# ============================================
# FUNCTIONS TO CALL BASED ON TYPE OF CALCULATION
# ============================================

set_kpoints() {
  local celltype=$1
  local sym=$2

  echo "running kpoint function"
  # first check if using manual kpoint array
  if [ "$use_manual_kpoints" = true ]; then
      echo "using manual kpoint array"
      convergence_array=("${manual_kpoint_array[@]}")

  # if not, generate automatically based on symmetry
  else
    echo "generating kpoint array from symmetry"
    case "$sym" in
    "FCC"|"BCC")
      # For cubic systems, uniform in all directions
      # so set k1 to loop from 2-->40 in incremenets of 2, then k2 and k3=k1, and create array of k1 k2 k3
        convergence_array=()
        for k1 in {2..40..2}; do
          k2=$k1
          k3=$k1
          convergence_array+=("$k1 $k2 $k3")
        done
      ;;

    "HCP"|"Orthorhombic")
      # create an empty array. For values 2 --> 40 in increments of 2,
      # calculate kpoints along b and c based on b/a and c/a ratios (since kpoints scale with length in real space)
      convergence_array=()
      for k1 in {2..40..2}; do
        k2=$(awk -v k1="$k1" -v ba="$b_a" 'BEGIN { printf "%.0f", k1/ba }')
        k3=$(awk -v k1="$k1" -v ca="$c_a" 'BEGIN { printf "%.0f", k1/ca }')

        # round b and c kpoint numbers to nearest EVEN number (odd numbers can increase compute time due to symmetry)
        for var in k2 k3; do
          val=$(eval echo \$$var)
          if (( val % 2 != 0 )); then
            (( val++ ))
          fi
          eval $var=$val
        done

        # take rounded numbers and input into array
        convergence_array+=("$k1 $k2 $k3")
      done
      ;;
    esac
  fi

  # for kpoints, just use largest recommended ecut values from PP files
  ecutwfc=$max_ecut
  ecutrho=$max_ecutrho

  # then, create a variable called 'loop'
  # each type of calculation will have a different loop to be inserted into the script file
  # this one iterates over the kpoint array, extracts k1/k2/k3, then creates kpoint variable
  echo "generating kpoint loop to insert into script"
  loop=$(cat << 'EOF'
for k in "${convergence_array[@]}"; do
  read -r k1 k2 k3 <<< "$k"
  value=${k1}
  kpoints="${k1} ${k2} ${k3}"
EOF
)

}

set_ecutwfc() {
  local celltype=$1
  local sym=$2

  echo "running ecutwfc function"
  echo "generating ecutwfc array"
  convergence_array=()
  for ec in {20..120..5}; do
    convergence_array+=("${ec} $(( ec * ecutrho_ratio ))")
  done

  echo "generating ecutwfc loop to insert into script"
  loop=$(cat << 'EOF'
for e in "${convergence_array[@]}"; do
  read -r ecutwfc ecutrho <<< "$e"
  value=${ecutwfc}
EOF
)

}

set_ecutrho() {
  local celltype=$1
  local sym=$2

  echo "running ecutrho function"

  echo "check ecutwfc ( $ecutwfc ) vs recommended ( $max_ecut )"
  # check if converged ecutwfc is less than the recommended from PP. If that's the case, just use the PP value
  # this is because we may find that convergence happens before reaching the recommended minimum value from PP
  # increasing ecutwfc above convergence in this case will just increase accuracy, so can just up it to the recommended
  if [[ $ecutwfc -lt $max_ecut ]]; then
      ecutwfc=$max_ecut
      echo "Below recommended. Updating ecutwfc ( $ecutwfc ) --> recommended ( $max_ecut )"
  fi

  echo "generating ecutrho array"
  convergence_array=()
  for ec in {4..14..1}; do
    convergence_array+=("${ecutwfc} $(( ecutwfc * ec ))")
  done

  echo "generating ecutrho loop to insert into script"
  loop=$(cat << 'EOF'
for e in "${convergence_array[@]}"; do
  read -r ecutwfc ecutrho <<< "$e"
  value=${ecutrho}
EOF
)
}

set_lattice() {
  local celltype=$1
  local sym=$2

  echo "running lattice convergence function"
  # check if converged ecutwfc is less than the recommended from PP. If that's the case, just use the PP value
  # this is because we may find that convergence happens before reaching the recommended minimum value from PP
  # increasing ecutwfc above convergence in this case will just increase accuracy, so can just up it to the recommended
  echo "check ecutwfc ( $ecutwfc ) vs recommended ( $max_ecut )"
  if [[ $ecutwfc -lt $max_ecut ]]; then
      ecutwfc=$max_ecut
      echo "Below recommended. Updating ecutwfc ( $ecutwfc ) --> recommended ( $max_ecut )"
  fi

  echo "check ecutrho ( $ecutrho ) vs recommended ( $max_ecutrho )"
  if [[ $ecutrho -lt $max_ecutrho ]]; then
      ecutrho=$max_ecutrho
      echo "Below recommended. Updating ecutrho ( $ecutrho ) --> recommended ( $max_ecutrho )"
  fi

  case "$sym" in
    "FCC"|"BCC"|"Orthorhombic")
     convergence_array=()
      while IFS= read -r a; do
        convergence_array+=("$a")
      done < <(seq 2.5 0.02 5)
    ;;
    "HCP")
      convergence_array=()
      while IFS= read -r a; do
        celldm1=$(bc -l <<< "${a} / 0.529177525830478")
        convergence_array+=("$celldm1")
      done < <(seq 2.5 0.02 5)
    ;;
  esac

  echo "generating lattice parameter loop to insert into script"
  loop=$(cat << 'EOF'
for a in "${convergence_array[@]}"; do
  if [ "$lattice_symmetry" = "HCP" ]; then
    lattice_parameter="celldm(1) = ${a}
celldm(3) = ${c_a}"
  else
    lattice_parameter="A=${a}"
  fi
  value=${a}
EOF
)

}

set_volume_optimize_birch_murnaghan() {
  local sym=$1
  local calc=$2
  declare -g output

  # extract lattice parameter from file name and energy (in Ry) from output files

  # hcp, volume is (sqrt(3)/2) * a^2 * c
  if [ "$sym" = "HCP" ]; then
    # first, extract all lattice parameters, calculate volume, and store volume-energies
    awk -v c_a="$c_a" '
    /kinetic-energy/{ecut=$4}
    /^!.*total/{
      match(FILENAME, /lattice([0-9.]+)\.out/, arr)

      # convert a from Bohr → Ang
      a = arr[1] * 0.529177525830478
      c = a * c_a

      # compute volume
      volume = (sqrt(3)/2) * a * a * c

      print volume, $5
    }' ${system_directory}/${calc}/*.out \
      > ${system_directory}/${file_naming}_volume_energy.txt

  elif [ "$sym" = "Orthorhombic" ]; then
    awk -v c_a="$c_a" -v b_a="$b_a" '
    /kinetic-energy/{ecut=$4}
    /^!.*total/{
      match(FILENAME, /lattice([0-9.]+)\.out/, arr)

      a = arr[1]
      b = a * b_a
      c = a * c_a

      # compute volume
      volume = a * b * c

      print volume, $5
    }' ${system_directory}/${calc}/*.out \
      > ${system_directory}/${file_naming}_volume_energy.txt

  else
    # cubic systems feed lattice parameter into birch-murnaghan not volume
    awk '
    /kinetic-energy/{ecut=$4}
    /^!.*total/{
      match(FILENAME, /lattice([0-9.]+)\.out/, arr)
      print arr[1], $5
    }' ${system_directory}/${calc}/*.out \
      > ${system_directory}/${file_naming}_volume_energy.txt
  fi

  # set variables to fit to birch-murnaghan
  case "$sym" in
    "HCP"|"Orthorhombic")
      type="noncubic"
    ;;
    "FCC")
      type="fcc"
    ;;
    "BCC")
      type="bcc"
    ;;
  esac

  # now, birch-murnaghan is only valid for volumes +-10% of ideal. So, we want to filter our .txt files
  # first, calculate the volume with the lowest energy in the .txt files, then calculate +10%/-10%
  # then filter out any lines beyond that range
  datafile="${system_directory}/${file_naming}_volume_energy.txt"
  # Find volume with minimum energy
  v0=$(awk 'NR==1 || $2 < min_energy {min_energy=$2; min_vol=$1} END{print min_vol}' "$datafile")

  # Compute range of volumes
  vmin=$(awk -v v0="$v0" 'BEGIN{printf "%.6f", v0*0.9}')
  vmax=$(awk -v v0="$v0" 'BEGIN{printf "%.6f", v0*1.1}')

  echo "Equilibrium volume: $v0"
  echo "Lower bound: $vmin"
  echo "Upper bound: $vmax"

  # Filter original file into a new trimmed file
  awk -v vmin="$vmin" -v vmax="$vmax" '$1 >= vmin && $1 <= vmax' "$datafile" \
    > "${datafile%.txt}_trimmed.txt"

  unit="Ang"
  input="${datafile%.txt}_trimmed.txt"
  output="${system_directory}/${file_naming}_birch_murnaghan.txt"
  echo -e "${unit} \n${type} \n2 \n${input} \n${output}" | ev.x
}

set_extract_BM_data() {
  declare -g V0

  # from the birch-murnaghan output file, extract ideal volume.
  # this accounts for unit's shown as either Ang^3 OR A^3. Also accounts for scientific notation
  V0=$(grep -oP 'V0\s*=\s*\K[0-9.E+-]+(?=\s*A(?:ng)?\^?3)' "$output")
  A0=$(grep -oP 'a0\s*=\s*\K[0-9.E+-]+(?=\s*A(?:ng)?)' "$output")
  B0=$(grep -oP 'k0\s*=\s*\K[0-9.E+-]+(?=\s*GPa)' "$output")

  # save a log of BM findings
  output_file="${compute_directory}/log_of_DFT.txt"
  echo "$(date '+%Y-%m-%d %H:%M:%S') BM-lattice opt System = $file_naming V0 = $V0 Ang^3 a0 = $A0 Ang B0 = $B0 GPa" >> "$output_file"
}

set_geometry() {
  echo "running geometry optimization for HCP or Orthorhombic"

  echo "check ecutwfc ( $ecutwfc ) vs recommended ( $max_ecut )"
  if [[ $ecutwfc -lt $max_ecut ]]; then
      ecutwfc=$max_ecut
      echo "Below recommended. Updating ecutwfc ( $ecutwfc ) --> recommended ( $max_ecut )"
  fi

  echo "check ecutrho ( $ecutrho ) vs recommended ( $max_ecutrho )"
  if [[ $ecutrho -lt $max_ecutrho ]]; then
      ecutrho=$max_ecutrho
      echo "Below recommended. Updating ecutrho ( $ecutrho ) --> recommended ( $max_ecutrho )"
  fi

  # now for HCP and orthorhombic, have to vary c/a and b/a, respectively, and run a bunch of scf calculations at V0
  # then, based on which ratio has the lowest energy, that is the ideal ratio
  case "$symmetry" in
  "HCP")
    # for HCP, we have ideal volume V0, a=b, and vary c/a.
    # So, for each c/a value, back-calculate a from V0, then convert to bohr and store in array
    # for HCP, V=(sqrt(3)/2)*a^2*c and c=a*c_a --> V=(sqrt(3)/2)*a^3*c_a --> a=(2*V/sqrt(3)*c_a)^(1/3)
    echo "back calculating a from ideal volume ( $V0 ) and c_a ratios --> creating array"
    convergence_array=()
    while IFS= read -r c_a; do
      a=$(bc -l <<< "e(l((2 * $V0) / (sqrt(3) * $c_a)) / 3)")
      a_bohr=$(bc -l <<< "${a} / 0.529177525830478")
      convergence_array+=("$a_bohr $c_a")
    done < <(seq 1.55 0.01 1.75)

    loop=$(cat << 'EOF'
for a in "${convergence_array[@]}"; do
    read -r celldm1 celldm3 <<< "$a"

    lattice_parameter="celldm(1) = ${celldm1}
celldm(3) = ${celldm3}"
    value=${celldm3}
EOF
)
    ;;
  "Orthorhombic")
    echo "back calculating a from ideal volume ( $V0 ) and b_a ratios --> creating array"
    # for orthorhombic, same process as HCP, but volume is easier
    # V=a*b*c and b=b_a*a and c=c_a*a --> V=a^3*b_a*c_a --> a=(V/b_a*c_a)^(1/3)
    convergence_array=()
    while IFS= read -r b_a; do
      a=$(bc -l <<< "e(l($V0 / ($b_a * $c_a)) / 3)")
      convergence_array+=("$a $b_a $c_a")
    done < <(seq 1.55 0.01 1.75)

    loop=$(cat << 'EOF'
for a in "${convergence_array[@]}"; do
    read -r lattice b_a c_a <<< "$a"
    lattice_parameter="A=${lattice}"
    value=${b_a}
    alat_b=$(printf "%.15f" "$b_a")
    alat_c=$(printf "%.15f" "$c_a")
EOF
)

    cell_parameters=$(cat << 'EOF'
"CELL_PARAMETERS {alat}
  1.000000000000000   0.000000000000000   0.000000000000000
  0.000000000000000   ${alat_b}   0.000000000000000
  0.000000000000000   0.000000000000000   ${alat_c}"
EOF
)
    ;;
  esac

}

# ============================================
# CALL DESIRED FUNCTIONS
# ============================================

set_pseudo_atomic_info

# if using a .cif file, need to extract information like starting_lattice, nat, and atomic_positions
# will use cif2cell to convert from .cif to .in, which we can then use to fill in our expanded .in file
if [ "$using_cif_file" = true ]; then
  echo "running .cif function"
  set_input_parameters_cif "$cell_type" "$symmetry"
  set_lattice_parameter_cif "$calculation" "$symmetry"
else
  echo "running no .cif function"
  set_input_parameters_no_cif "$cell_type" "$symmetry"
  set_lattice_parameter_no_cif "$calculation" "$symmetry"
fi

if [ "$calculation" = "geometry" ]; then
  # no matter symmetry, fit to birch-murnaghan then extract the data
  module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 intel/2023.2.1 gcccore/.12.3
  module load quantumespresso/7.3.1

  echo "Before optimizing geomtry, fit to B-M..."
  set_volume_optimize_birch_murnaghan "$symmetry" lattice
  echo "extracting B-M results..."
  set_extract_BM_data
fi

case "$calculation" in
  "kpoints"|"ecutwfc"|"ecutrho"|"lattice")
    set_${calculation} "$cell_type" "$symmetry"
  ;;
  "geometry")
    if [[ "$symmetry" == "HCP" || "$symmetry" == "Orthorhombic" ]]; then
  set_${calculation} "$cell_type" "$symmetry"
  else
    echo "no need for geometry optimization for symmetry $symmetry"
  fi
  ;;
esac

########################################## GENERATE SCRIPT ##########################################

# go to directory
cd ${system_directory}

# create script that makes all .in files, then runs them and extracts energy info
cat > "${file_naming}"."${calculation}".sh << EOF2
#!/bin/bash
#SBATCH --time=0-01:00:00
#SBATCH --account=def-diak01
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000
#SBATCH --output=${file_naming}_${calculation}.out
#SBATCH --job-name=${file_naming}_${calculation}
#SBATCH --mail-user=karl.steffensen@queensu.ca
#SBATCH --mail-type=ALL

module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 intel/2023.2.1 gcccore/.12.3
module load quantumespresso/7.3.1

lattice_symmetry="${symmetry}"

# pull values from input script
# in the case of convergence, each value will be overwritten during respective calculation
kpoints=("${kpoint_grid}")
ecutwfc=${ecutwfc}
ecutrho=${ecutrho}
b_a=${b_a}
c_a=${c_a}

cd ${system_directory}
mkdir -p "${calculation}"
cd "${calculation}"

# now, insert a loop that will be used to create different .in files
EOF2

# Write the exact array contents into the generated file
declare -p convergence_array >> "${file_naming}"."${calculation}".sh

# Append the rest of the script logic
cat >> "${file_naming}"."${calculation}".sh << EOF2
${loop}

  cat > "${file_naming}"."${calculation}""\${value}".in <<EOF1
&CONTROL
calculation='scf',
outdir='.',
prefix='${file_naming}_${calculation}\${value}',
pseudo_dir='${pseudo_directory}',
verbosity='low',
tprnfor=.true.,
tstress=.true.,
/

&SYSTEM
ibrav = ${ibrav}
${lattice_parameter}
${nat}
${ntyp}
ecutwfc=\${ecutwfc},
ecutrho=\${ecutrho},
occupations='smearing',
smearing='mp',
degauss=0.03,
/

&ELECTRONS
conv_thr=1d-08,
mixing_beta=0.7d0,
diagonalization='cg',
/

${cell_parameters}

${atomic_species}

${atomic_positions}

K_POINTS {automatic}
\${kpoints} ${kshift}

EOF1
done

# run the calculations
directory="${system_directory}/${calculation}"

for input_file in \$(ls \${directory}/*.in)
do
    srun pw.x -in "\$input_file" > "\${input_file%.*}.out"
done

cd ${system_directory}

calculation=${calculation}

# here, based on the type of calculation, extract the relevant results

case "\$calculation" in
  "lattice"|"geometry")
    awk -v calc="\$calculation" '
/^!.*total/ {
  regex = "\\\." calc "([0-9.]+)\\\.out\$"
  match(FILENAME, regex, arr)
  val = arr[1]
  print val, \$5
}' \${directory}/*.out \
| sort -n \
> ${file_naming}_\${calculation}_output.txt
    ;;
  *)
    awk -v calc="\$calculation" '
/^!.*total/ {
  regex = "\\\." calc "([0-9]+)\\\.out\$"
  match(FILENAME, regex, arr)
  val = arr[1]
  print val, \$5
}' \${directory}/*.out \
| sort -n \
> ${file_naming}_\${calculation}_output.txt

  awk '
NR==1 { prev=\$2; print \$1, \$2, 0; next }
{ diff = (\$2 - prev >= 0) ? \$2 - prev : prev - \$2; print \$1, \$2, diff; prev=\$2 }
' ${file_naming}_\${calculation}_output.txt > tmp && mv tmp ${file_naming}_\${calculation}_output.txt
    ;;
esac


EOF2

case "$calculation" in
  "geometry")
    if [[ "$symmetry" == "HCP" || "$symmetry" == "Orthorhombic" ]]; then
      sbatch "${file_naming}"."${calculation}".sh
    else
      echo "All done, no need for geometry optimization for symmetry $symmetry"
    fi
    ;;
  *)
    sbatch "${file_naming}"."${calculation}".sh
    ;;
esac
