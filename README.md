# this file will continually be updated as new scripts are incorporated into the workflow

####################################################################################################

THIS IS AN INSTRUCTION MANUAL FOR RUNNING DFT CONVERGENCE AND CALCULATIONS FOR ALUMINUM ALlOYS USING THE SCRIPTS THAT I CREATED

####################################################################################################

# General workflow of script:

**scf calculations:**
1. Run convergence calculations (kpoints, ecutwfc, ecutrho).
2. Using convergence parameters, optimize lattice parameter (test large range of potential values). For non-cubic systems, b/a and c/a ratios kept constant.
3. Take lattice parameter results and fit to Birch-Murnaghan equation of state (for cubic systems, feed lattice parameter+energy, for non-cubic feed volume+energy). **Note:** For non-cubic, the script calculates the volume of the cell for the given lattice parameter, then generates a volume-energy.txt file. For cubic, the script extracts lattice parameter/energy and builds a lattice parameter-energy.txt file. The script automatically filters these files, identifying the lattice parameter/volume with the minimum energy, and keeps all other lattice parameter/volume lines within +-10% (Birch-Murnaghan is valid within 10% of ideal volume). Filtered .txt files are fit to Birch-Murnaghan.
4. Ideal volume/lattice parameter are extracted from Birch-Murnaghan fit. For non-cubic, ideal lattice parameter a is back-calculated (because b/a and c/a have been kept constant).
5. **For non-cubic only:** Because b/a (orthorhombic) and c/a (HCP) can vary (constant in cubic), there is an extra degree of freedom to optimize. For a **constant volume (equal to ideal volume from B-M fit)**, run geometry optimization testing a wide range of b/a (orthorhombic, with constant c/a) or c/a (HCP, where a=b so b/a=1).
6. Isolate ideal unit cell parameters. For cubic systems, this is ideal lattice parameter from Birch-Murnaghan fit. For HCP, this is c/a ratio with the lowest energy from geometry optimization. For Orthorhombic, this is b/a ratio with the lowest energy from geometry optimization.

**Non-scf calculations**
1.  Relax atomic positions for ideal unit cell parameters ('relax' calculation, fix unit cell shape/volume+allow atoms to move).
2.  Run variable cell relaxation ('vc-relax', allow atoms, unit cell volume, and unit cell shape to vary).

**Output here will be a unit cell that has been fully optimized in every degree of freedom**

**Note:** Because vc-relax allows all dimensions to vary at once, you could just use that from your starting structure (after convergence) and skip the optimization steps (lattice parameter, geometry). The only problem is that if your starting structure is far away from the ideal structure, this can take a long time. By running optimization steps, the calculations are much faster because we're only varying 1 degree of freedom at a time. This gets us much closer to the final, optimized structure, so vc-relax doesn't have as much that it needs to vary (resulting in a net decrease in calculation time, even if there are more steps involved).


# PREPARE FOR RUNNING DFT:

1. Download all pseudopotential files for the elements you want to use. Include PAW, USPP, and any pseudopotential-XCF combination (PAW-PBE, PAW-PBESol, USPP-PBE etc.). Include multiple pseudopotentials per element in case of needing to test for the most effective pseudopotential. 

2. Put all pseudopotentials into a folder and copy to cluster.

3. In cluster, create a folder for all DFT scripts and copy them over (should be input.sh, scf_combined.sh, relaxation_combined.sh).

4. In cluster, create a folder for ‘index’ files (more on that **below**).

5. Make sure pseudo potential file is copied over.

6. In ‘index file’ folder, create index file for atomic weights (atomic_weights.txt). Inside have 2 columns, first with the elemental symbol, second with its atomic weight. The scripts look up this information to help generate input files. Any element you want to calculate should have an entry. 

**Example**:
```
Al 26.9815
Ce 140.116
Cu 63.546
Fe 55.845
```

7. Create an index file for pseudo potentials (“list_of_pseudopotentials.txt”).Inside, have 5 columns: Element, pseudopotential type (PAW vs USPP), relativism (Full vs Scalar), exchange correlation functional (PBE, PBESol etc.), and the pseudopotential file name corresponding to the choices in the other 4 columns (each has a unique name). This information is used by the script to generate input parameters. 

**Example:**

```
Element PP_type Relativism	XCF     PP_name
Al	PAW     Scalar  PBE     Al.pbe-n-kjpaw_psl.0.1.UPF
```


####################################################################################################

# CREATE AND MODIFY INPUT FILE:

Copy input.sh to cluster (or create manually via nano). Inside the file:

1. Convergence results
	•	**If convergence has already been run**, enter the associated results (ie kpoints, ecutwfc, ecutrho).
	•	If you haven’t run convergence tests yet, leave the variables blank for the ones you haven’t converged (ie if you just converged kpoints, leave ecutwfc and ecutrho blank).
	•	Kpoints should be in format “12 12 12”, where the numbers are all in quotations and separated by a space.

**Example**:
```
# kpoint grid taken from kpoint convergence results
kpoint_grid="16 16 10"
# ecutwfc value taken from ecutwfc convergence results
ecutwfc=90
# ecutrho value taken from ecutrho convergence results
ecutrho=450
```

2. System information
* First, enter system name within the quotation marks (chemical formula, ie “Al”, “AlCe”, “AlCe3”, “Al11Ce3” etc.).
* Second, enter the unique elements in each system in quotations, separated by a space (ie “Al”, “Al Ce”, “Al Si Ce” etc.).
* Third, enter paths to general directory containing ALL subfolders for your index files, pseudopotential files, calculations etc (this is ‘compute_directory’). Enter path to location of all your scripts (‘script_directory’), path to location of index files from above (‘index_directory’), and a path to pseudopotential files (‘pseudo_directory’). 
* NOTE: In all directory paths, make sure you don’t include a slash at the end
    * ie, use "/home/username/scratch" NOT "/home/username/scratch/“, because the extra slash at the end will cause the system to be unable to find the files its looking for. For example, if you need to find “Al.cif”, using the first option it will look for "/home/username/scratch/Al.cif”, which it will find as long as "/home/ksteffen/scratch/QUEENS" is the correct path. Using the slash at the end will cause it to look for "/home/username/scratch//Al.cif”. That extra slash will mess it up, it won’t be able to find the correct .cif file, and the script will fail. 

**Example:**
```
# chemical formula of system (ie "Al", "Al11Ce3" "AlFe" etc.)
system="Al"
# each unique element, formatted like "Al Fe Ce" etc. Each element separated by a space
elements="Al"
# overall directory containing ALL subfolders
compute_directory="/home/username/scratch/"
# path to all DFT scripts (.sh and .py)
script_directory="/home/username/scratch/DFT_scripts"
# path to directory with index files (atomic_weights.txt and list_of_pseudopotentials.txt)
index_directory="/home/username/scratch/DFT_index_files"
# path to pseudopotential files
pseudo_directory="/home/username/scratch/pseudopotentials"
```


3. Lattice information
* First, enter symmetry type (FCC, BCC, HCP, Orthorhombic).
    * Note: Have not setup script for other types of cells yet (like triclinic).
* State if you want to use unit or primitive cell (typically primitive cell unless orthorhombic).
* Then enter “true” or “false” if using a .cif file (if true, the script will extract parameters from the .cif file, if false you have to supply additional information about the lattice).

**Example:**
```
# HCP, FCC, BCC, Orthorhombic  
symmetry="HCP"
# Unit or Primitive
cell_type="Primitive"

# If you have a .cif file, we can extract data using cif2cell output. Set to true. Otherwise, false.
# NOTE: cif_file has to be named after the system. Ex: Al.cif, AlFe.cif, Al11Ce3.cif etc.
# whatever the you have for "system" variable above is the name we need for .cif file
# NOTE: IF NOT USING A CIF FILE, YOU MUST ENTER ATOMIC POSITIONS MANUALLY IN SECTION BELOW
# true or false
using_cif_file=false
```


4. Lattice information (if NO .cif)
* NOTE: If using a .cif file, you can leave these variables blank. 
* atoms is the number of atoms in the cell (based on unit vs primitive cell, type of symmetry etc.).
* starting_lattice is just an initial guess of lattice parameter that will be used for convergence calculations. Doesn’t really matter what it is, as long as it’s not wildly wrong (like 1 Angstrom or 20 Angstroms or something).
* b_a is the b/a ratio (usually 1 unless orthorhombic cell, because a=b for all cubic and HCP systems. Good initial value maybe 1.633 for Orthorhombic).
* c_a is the c/a ratio (1 for all cubic, variable for HCP systems. Good initial value is 1.633 for HCP, maybe 1.73205081 (sqrt(3)) for Orthorhombic).
* **NOTE:** Final step of structure optimization will vary b/a and c/a for us, so our initial guesses don’t have to be accurate, we just want to reduce the computational time of that final calculation by honing in on an optimized structure. 
* Finally, enter manual atomic positions. I prefer crystal format (positions along lattice vectors a,b,c. Other options are xyz coordinates in terms of Angstroms). For most primitive cells, just 1 atom at origin (0, 0, 0). For HCP primitive, 1 atom at origin (0, 0, 0) and 2nd atom at (2/3, 1/3, 1/2). Input files I’ve worked with so far have had 15 decimal points for these coordinates, so try to keep that just in case. The left-most column for each line of atomic coordinates has the element for which those coordinates apply. 
    * These are automatically extracted from .cif file if using one 

Example:
```
  # number of atoms in cell. If using .cif, leave blank
  natoms=2

  # INITIAL lattice parameter (pre-optimization), in angstroms. If using a .cif, leave blank
  starting_lattice=3.789

  # if orthorhombic or hcp and NO .cif file, need to manually input b/a or c/a ratios.
  # For orthorhombic, start with: b_a=1.633, c_a=1.73205081 <-- sqrt(3)
  # for HCP, start with: b_a=1, c_a=1.633
  # doesn't REALLY matter what these values are, because we'll run vc-relax which will vary these for us. We just want
  # to reduce computational expense/time of vc-relax by getting as close as possible through manual calculations
  b_a=1
  c_a=1.633

  # If no .cif file, enter atomic positions manually
  # here, just modify atomic positions based on the symmetry/elements you're using
  # remember, for primitive cells put 1 atom at origin
  # for primitive hcp, 1 at origin and 1 at (2/3, 1/3, 1/2)
  manual_atomic_positions="ATOMIC_POSITIONS {crystal}
Al 0.000000000000000 0.000000000000000 0.000000000000000
Al 0.666666666666666 0.333333333333333 0.500000000000000"

```

5. Calculation information
* Need to tell the system information about our pseudopotential type (PAW or USPP), whether we are using full relativistic (“Full”) or scalar relativistic (“Scalar”), as well as the exchange correlation functional (PBE, PBESOL, LDA, R2SCAN).
* This is important because based on these inputs, the script will pull the corresponding pseudopotential file from our folder containing all of our pseudo potentials (first, it ‘looks’ in our index file to find the line matching all of our choices here for each atom, then it pulls the pseudopotential file name on that line, then it pulls that corresponding file for the calculations).

**Example:**
```
# PAW or USPP
pseudopotential_type="PAW"
# scalar ("Scalar") or full relativistic ("Full")
pseudo_rel_type="Scalar"
# PBE, PBESOL, LDA, or R2SCAN
functional_type="PBE"
```


6. kpoint parameters
* The main things here are telling it whether you are using a manual kpoint array (“true”) or not (“false”), as well as whether you want gamma_centered kpoint array (“true”), which is the default.
* If set gamma_centered to “false”, it will do 1 1 1 shift in kpoints.
* If entering a manual kpoint array, you have to have sets of 3 numbers for what you want to test. Each set of 3 should be in quotations with a space in between each number, and ALL SETS of 3 numbers should be enclosed in round brackets. Have to do this even if they’re all the same like in a cubic system (but you shouldn’t use manual in that case because we sample a huge range of kpoints for that). The main reasons to use manual kpoint would be if you’re doing slab calculations or supercells or something and need different numbers of kpoints in each direction, or if you want to sample odd numbers of kpoints (my script only uses even because using odd numbers can sometimes decrease computing efficiency).

**Example:**
```
# set to true to define manual_kpoint_array instead of calculating it.
# true or false. Generally do false unless doing slab calculation or custom cell shape
use_manual_kpoints=true

# kpoint array. For HCP and orthorhombic where a=b, kpoints are the same along a/b, and different over c
# but if a =/ b =/ c, then Kpoints will be different in each direction. Hav eto make them roughly proportional to
# b/a and c/a ratios. Script calculates kpoints based on these ratios.
# Can't automate for all possibilities though, so allow for option to enter manually
# If converging kpoints for orthorhombic , enter the kpoint arrays you want to sample over
# for example, if b=c you could do ("2 1 1" "4 2 2" "6 4 4" "8 6 6" "10 6 6" "12 8 8") etc.
manual_kpoint_array=(“4 2 1” “6 4 2” “8 6 4” “10 8 6” “12 10 8” “14 12 10” “16 14 12” “18 16 14” “20 18 16”)

# true for Γ-centered grid (default); false for shifted grid
gamma_centered=true
```
 
####################################################################################################

# RUNNING THE DFT CALCULATIONS:

To run, just navigate to the folder containing all of your DFT scripts.

Then, for each script use chmod +x to turn them into executable files. Then, once your input file has been fully modified as necessary, run it as an executable file with ./

**Example:**
```
cd /home/ksteffen/scratch/DFT_scripts

chmod +x input.sh
chmod +x scf_combined.sh
chmod +x relaxation_combined.sh
chmod +x optimization_analyses.sh

./input.sh
```


Now, the script will submit a job that runs the specified type of calculation in your input file.

####################################################################################################

# WHAT TO DO AFTER EACH STEP:

Some manual checking is still required after each step of the process (ie checking the best parameters to use for convergence). Each step of the process outputs a .txt file that you can easily read with ‘nano’. 

For all steps, go to your main directory containing all of your results. For example, if you’re looking at FCC Aluminum with no other elements, this would be a folder called “Al_FCC”. Within that folder will be subfolders for each calculation type (called ‘kpoints’, ‘ecutwfc’, ‘ecutrho’, ‘lattice’, ‘geometry’, ‘relax’, and ‘vc-relax’). In this case, “Al_FCC” will be created wherever you specify the “compute_directory” variable to be in your input file. Below is what you have to do after each step:

	•	kpoints: 
	⁃	Once kpoint convergence has been run, check the “kpoints_output.txt” file. Each line shows the results of an scf calculation for a given kpoint grid. 
	⁃	The 1st column lists the first kpoint number in that grid (ie if the grid is “12 10 8”, then this will show 12). 
	⁃	The 2nd column lists the final energy of the system. 
	⁃	The 3rd column lists the difference (absolute value) between the energy of that line and the energy of the previous line.
	⁃	What we want is to find a k-point grid that results in little to no difference in energy (that is, the system can’t improve the calculation any further, so using more kpoints past that point is just a waste of computing infrastructure). 
	⁃	What you have to do is decide on a suitable kpoint array based on the energy difference column. Choose an array that balances energy difference for computational time, then input that for “kpoint_grid” variable at the top of the input.sh file. 

Example:

For an output file containing:

2 -78.87907422 0
4 -78.80528285 0.0737914
6 -78.81289438 0.00761153
8 -78.81200825 0.00088613
10 -78.81306237 0.00105412
12 -78.81233947 0.0007229
14 -78.81253245 0.00019298
16 -78.81249354 3.891e-05
18 -78.81251273 1.919e-05
20 -78.81251090 1.83e-06
22 -78.81250434 6.56e-06
24 -78.81250673 2.39e-06
26 -78.81250652 2.1e-07
28 -78.81250639 1.3e-07
30 -78.81250636 3e-08
32 -78.81250632 4e-08
34 -78.81250634 2e-08
36 -78.81250634 0
38 -78.81250634 0
40 -78.81250633 9.99999e-09

We can see that beyond the kpoint grid “16”, the differences between energies are incredibly small, but every increase in kpoints requires more computing time and power. In this example, it was for an HCP structure, so if I go to the kpoints directory and look into the “kpoints16.in” file, I can see that the kpoint array for this is “16 16 10”. Now I’ll use “16 16 10” as my kpoint grid for any future calculations in this system. 



	•	ecutwfc:
	⁃	After running ecutwfc convergence, look in “ecutwfc_output.txt” file. 
	⁃	Like with kpoints, first column is the value tested, 2nd is the energy, 3rd is the difference in energy. 
	⁃	Again, we want to choose a value for ecutwfc after which the energy differences are consistently very small and  input value for ‘ecutwfc’ variable at the top of input.sh. 
	⁃	NOTE: Each pseudopotential file has a value for ecutwfc and ecutrho inside of them. These are the minimum recommended values to use for that element. We must always use those values or higher. In the case of multiple elements, we should always use whatever values are the highest (ie if Al pseudopotential suggests ecutwfc of 40, and Fe pseudopotential suggests ecutwfc of 60, we should use 60 or higher). If our convergence values end up being below these pseudopotential recommendations, the script will automatically update to the pseudopotential recommendations, so you don’t have to worry about checking to make sure your convergence values are above the recommended, just choose values at which point the system is converged. 

Example:

20 -78.81129574 0
30 -78.81193593 0.00064019
40 -78.81233275 0.00039682
50 -78.81232970 3.05e-06
60 -78.81219853 0.00013117
70 -78.81228469 8.616e-05
80 -78.81217917 0.00010552
90 -78.81223206 5.289e-05
100 -78.81218015 5.191e-05
110 -78.81218613 5.98e-06
120 -78.81218584 2.9e-07

Here, the energy difference is consistently 10^-5 or less after an ecutwfc of 90, so going forward use ecutwfc=90 for all calculations in this system. As stated above, if this value is less than the recommended value in the pseudopotential file, it will automatically be updated for you (because if 90 represents convergence, and the pseudopotential file recommends a minimum of 100, then we know 100 is converged because it comes after 90 where the system is converged, but use the 100 instead because its the minimum recommended value).



	•	ecutrho:
	⁃	Same process as before, look in “ecutrho_output.txt”
	⁃	With ecutrho convergence, we are testing different values based on a ratio of ecutrho to ecutwfc. For PAW pseudopotentials, a ratio of ecutrho:ecutwfc of 4 or 5 is usually good, while USPP recommends a ratio of 10-12. 
	⁃	For kpoint and ecutwfc convergence, the script just calculates the ratio of ecutrho:ecutwfc from the recommended pseudopotential values, and uses that to calculate ecutrho based on the ecutwfc values we’re trying. Here, we test a range from 4 to 14.
	⁃	Again, choose the ecutrho value representing convergence and input value for ‘ecutrho’ variable at the top of input.sh. 

Example:

360 -78.81223206 0
450 -78.81219334 3.872e-05
540 -78.81216520 2.814e-05
630 -78.81214941 1.579e-05
720 -78.81215318 3.77e-06
810 -78.81215506 1.88e-06
900 -78.81214358 1.148e-05
990 -78.81214688 3.3e-06
1080 -78.81214498 1.9e-06
1170 -78.81214130 3.68e-06
1260 -78.81214414 2.84e-06

In this example, convergence occurs at 450 and above (which, since we decided on an ecutwfc of 90 in the previous step, is a ratio of ecutrho:ecutwfc of 5:1, in line with the recommendations for PAW pseudopotentials). 



At this point, the only manual work you have to do is checking that everything is running properly (ie after updating input.sh file, then ./input.sh to run the next calculation type, checking the script file it generates, check the .in/.out files to make sure everything is formatted properly etc.). 

From here: 
	•	Running ‘lattice’ as calculation type will test a huge range of lattice parameters. The script will pull out energies from each lattice parameter .out file.
	•	Running ‘geometry’ as calculation type will run full geometry optimization. First, however, we have to use the previous results to figure out what the ideal unit cell volume is.
	⁃	The script takes the lattice parameters from the previous step, calculates unit cell volume, determine which volume has the lowest associated energy (from the lattice parameter calculation), filter the data to only include volumes/energies where the volume is +-10% of the volume with the lowest energy (explained in a second), then will fit the filtered data to BIrch-Murnaghan equation of state
	⁃	Birch-Murnaghan fitting is valid within +-10% of the ideal volume, hence the filtering of the data above to get a better fit. 
	⁃	The script automatically adjusts the input file for B-M fitting based on the symmetry type (for cubic, it wants lattice parameter vs energy, for non cubic it wants volume vs energy).
	⁃	So, to prepare for fitting to B-M, the script will convert to Angstroms if needed (HCP systems have lattice parameter provided in Bohr units), calculate the volume, filter the data, and supply either lattice parameter+energy or volume+energy. 
	⁃	Then, from the B-M fit the script will extract the ideal volume (V0), bulk modulus estimate, and (for cubic systems only) the ideal lattice parameter A0.
	⁃	For cubic systems, we have the ideal volume/lattice parameter so we are done FOR NOW. For non cubic systems, we have the ideal volume from volume-energy relationship, but now need to figure out the best combination of a/b/c dimensions of the unit cell that give that volume. 
	⁃	So, the script takes the ideal volume and back-calculates the lattice parameter ‘a’ based on the b/a and c/a ratios. 
	⁃	For HCP, b=a so b/a=1, so want to test a range of c/a ratios to find the best one (and therefore the best a and b values) that produces the ideal volume. By varying c/a, the script back calculates ‘a’ from the ideal volume V0 for each possible c/a ratio, runs a calculation, and extracts the energy associated with each c/a ratio.
	⁃	For Orthorhombic, a=/b=/c. Too much to vary BOTH b/a and c/a (because as change 1, have to change the other), so keep c/a constant (if using .cif file keep it at the initial c/a ratio from the .cif file, if no .cif file use the manually input c_a variable). Then, using the same process as HCP, vary b/a ratio and back calculate lattice parameter ‘a’ from ideal volume V0, calculating energy associated with each b/a ratio.
	•	Running ‘relax’ calculation allows atomic relaxation (that is, atoms can now move) while keeping the unit cell in tact (at the ideal volume). 
	⁃	The script takes the ideal volume from before, then the ideal lattice parameter (a0 from B-M fit for cubic, calculated from ideal c/a for HCP based on which produces the lowest energy, or calculated from ideal b/a for Orthorhombic based on which produces the lowest energy), and allows the atoms to move in case there are any forces acting on the atoms. 
	⁃	NOW YOU HAVE TO MANUALLY CHECK THE ‘RELAX’ OUTPUT (.out) FILE. Within the .out file is a line that says “Begin final coordinates”. Find this and make sure the atomic positions haven’t changed (compare to manual ‘atomic_position’ input or initial .cif file coordinates). 
	⁃	If atomic positions haven’t changed, then you can continue onto the next calculation. 
	⁃	If atomic positions HAVE changed, then update your atomic positions by altering the input file manual)atomic_positions variable (even if you have been using a .cif file).


For the final step, run ‘vc-relax’ calculation. It is possible to bypass the lattice parameter, geometry optimization, and relax calculation steps by running vc-relax from the beginning. However, vc-relax allows atom positions AND unit cell volume/shape to all vary at the same time in order to find the ideal setup (the ideal cell shape/volume and atom spacing within that shape/volume). Because of the number of degrees of freedom, this can take a really long time depending on your system. Running the steps that we have up until now dramatically increases the computing speed/efficiency because each step gets us closer and closer to the ‘finish line’, so we just use vc-relax to get us the last little bit. By running a bunch of fast, low-computing power scf calculations, we shave off the heaviest computing steps. 

Now, once vc-relax has finished, go to the vc-relax.out file and scroll down to the line with “Beginning final coordinates”. Compare the “new unit-cell volume’, ‘cell parameters’, and ‘atomic positions’ to your previous ones (which can be found at the top of the .out file). Has the shape, size, or atomic positions changed substantially? If so, this is your ideal unit cell that you should use going forward. If not, then we had already found the ideal unit cell/atomic arrangements from our optimization process. Either way, we now have our ideal unit cell for running analyses or generating supercells. 

Example:

This is from the vc-relax.out file of our system above. We found that ‘relax’ did not change any atomic coordinates, so we moved onto vc-relax. At the top of the .out file is our input parameters:

     bravais-lattice index     =            4
     lattice parameter (alat)  =       5.4038  a.u.
     unit-cell volume          =     224.1118 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            1
     number of electrons       =         6.00
     number of Kohn-Sham states=            7
     kinetic-energy cutoff     =      90.0000  Ry
     charge density cutoff     =     450.0000  Ry
     scf convergence threshold =      1.0E-08
     mixing beta               =       0.7000
     number of iterations used =            8  plain     mixing
     energy convergence thresh.=      1.0E-04
     force convergence thresh. =      1.0E-03
     press convergence thresh. =      5.0E-01
     Exchange-correlation= SLA  PW   PBX  PBC
                           (   1   4   3   4   0   0   0)
     nstep                     =           50


     celldm(1)=   5.403769  celldm(2)=   0.000000  celldm(3)=   1.640000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (   1.000000   0.000000   0.000000 )
               a(2) = (  -0.500000   0.866025   0.000000 )
               a(3) = (   0.000000   0.000000   1.640000 )

And further down is the calculated unit cell:

Begin final coordinates
     new unit-cell volume =    224.12617 a.u.^3 (    33.21207 Ang^3 )
     density =      2.69805 g/cm^3

CELL_PARAMETERS (alat=  5.40376856)
   0.999382691  -0.000000000   0.000000000
  -0.499691346   0.865490799  -0.000000000
   0.000000000   0.000000000   1.642132016

ATOMIC_POSITIONS (crystal)
Al               0.0000000000        0.0000000000        0.0000000000
Al               0.6666666667        0.3333333333        0.5000000000
End final coordinates

As we can see, initial unit cell volume (which if you remember is our V0 from B-M fit) is 224.1118 a.u.^3, while our new volume is 224.12617 a.u.^3. Lattice parameter (which is celldm(1) in Bohr units since this is an HCP structure), went from 5.403769 to 5.40376856 (so unchanged). Atomic positions were unchanged. C/a (found in the bottom right of CELL_PARAMETERS) went from 1.64 to 1.642132016, so a very slight increase. Then CELL_PARAMETERS show a very slight change in the overall shape, as this lists xyz coordinates of the cell in terms of lattice parameter (so 1 lattice vector goes to (1,0,0)*a, 2nd lattice vector goes to (-0.5,0.866025,0)*a, and 3rd goes to (0,0,1.64)*a in our input structure). 

To know if these minuscule differences are meaningful, compare energy of input (from ‘geometry’) which is -79.00133768 Ry vs output of -79.00133929 Ry. Almost no difference. Then looking at stress/pressure tensors:


     Computing stress (Cartesian axis) and pressure

          total   stress  (Ry/bohr**3)                   (kbar)     P=        0.04
  -0.00000281  -0.00000000   0.00000000           -0.41       -0.00        0.00
  -0.00000000  -0.00000281   0.00000000           -0.00       -0.41        0.00
   0.00000000   0.00000000   0.00000640            0.00        0.00        0.94


     Computing stress (Cartesian axis) and pressure

          total   stress  (Ry/bohr**3)                   (kbar)     P=       -0.03
  -0.00000130  -0.00000000   0.00000000           -0.19       -0.00        0.00
  -0.00000000  -0.00000130   0.00000000           -0.00       -0.19        0.00
   0.00000000   0.00000000   0.00000207            0.00        0.00        0.30

Given our input parameter had a convergence threshold of 0.5 for pressure, and the energy only decreased by 2.2*10^-5 eV, we can consider these the same structures and our original from geometry optimization is fine to use if we wanted to. 



##################################################################################################################################

OPTIMIZED UNIT CELL ANALYSES:

First, extract elastic constants. 

	•	A GOOD RULE OF THUMB FOR KPOINTS FOR METALS:
	⁃	As per literature, when extracting elastic constants for metals, a good number of kpoints is approximately 7000/atom (so, for a 1 atom primitive cell something like 20x20x20)

