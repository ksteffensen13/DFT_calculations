import numpy as np
import sys

analysis_directory = sys.argv[1]

# create deformation matrix F, composed of 24 3x3 matrices (from our 24 deformation .txt files)
F = np.zeros((24, 3, 3))
# create a matrix composed of 24, 3x3 stress tensors (from our 24 stress tensor .txt files)
stress_mats = np.zeros((24, 3, 3))

# Also store which (strain, deformation) each entry belongs to
strain_index = np.zeros(24, dtype=int)
def_index = np.zeros(24, dtype=int)

# load our files in order (strain1_deformation1, strain1_deformation2, strain1_deformation3 etc.)
i = 0
for a in range(1, 5):
    for b in range(1, 7):
        deformation_file = f"{analysis_directory}/deformation_matrix_strain{a}_def{b}.txt"
        F[i] = np.loadtxt(deformation_file).reshape(3, 3)

        stress_file = f"{analysis_directory}/stress_tensor_strain{a}_def{b}.txt"
        stress_raw = np.loadtxt(stress_file)
        stress_mats[i] = np.loadtxt(stress_file).reshape(3, 3)

        # record mapping
        strain_index[i] = a - 1
        def_index[i] = b - 1

        i += 1

# create an empty stack of matrices "strain_mats" the same length as F
strain_mats = np.zeros_like(F)
F_T = np.zeros_like(F)
# create unit matrix I (100 010 001)
I = np.eye(3)
# then for each of the N 3x3 matrices, calculate the strain tensor via Green-Lagrange
for i in range(24):
    F_T[i]=F[i].T
    strain_mats[i] = 0.5 * (F_T[i] @ F[i] - I)

# Make functions to convert to Voigt notation
def strain_to_voigt_eps(e):
    """
    Convert strain tensor to Voigt form, which is (e11, e22, e33, 2e23, 2e13, 2e12)
    because matrix element indexes are 1 less, e[0,0] is e11, e[1,1] is e22 etc.
    also, strain matrix is symmetric so e13=e31, e21=e12, and e32=e23
    """
    return np.array([
        e[0,0], e[1,1], e[2,2],
        2*e[1,2], 2*e[0,2], 2*e[0,1]
    ], dtype=float)

def stress_to_voigt_eps(s):
    """
    Convert symmetric stress tensor to Voigt (s11, s22, s33, s23, s13, s12)
    This is the same as strain tensor, except it doesn't have 2*the last 3 elements
    """
    return np.array([
        s[0,0], s[1,1], s[2,2],
        s[1,2], s[0,2], s[0,1]
    ], dtype=float)

# Build E (strain) and S (stress) matrices from Voigt

# create a final E and S matrix
# these are both 6x6 matrices, and there should be 4 of them since there are 4 strain amounts
# each column of the E/S matrices is the voigt notation of 1 deformation type
# for example, column1 is voigt notation of deformation1, column2 is voigt notation of deformation2 etc.
# so first create 4, 6x6 matrices full of 0s
E = np.zeros((4, 6, 6))
S = np.zeros((4, 6, 6))

# for each of our 24 3x3 matrices
for n in range(24):
    # pull the strain index (going to be 0-3 representing strain1-4)
    si = strain_index[n]
    # and pull our deformation index (0-5 representing deformations1-6)
    dj = def_index[n]
    # basically, our first of 24 matrices (matrix index 0) is from strain1_deformation1
    # so, for matrix 1 (index0), pull strain index (0) and deformation index (0)
    # and fill

    E[si][:, dj] = strain_to_voigt_eps(strain_mats[n])
    S[si][:, dj] = stress_to_voigt_eps(stress_mats[n])



# Solve for C from equation S = C*E --> C = S*E^-1
# calculate inverse of E
E_pinv = np.linalg.pinv(E)
# then solve for C
C = E_pinv @ S
C1 = C[0]
C2 = C[1]
C3 = C[2]
C4 = C[3]
C_average = (C1 + C2 + C3 + C4) / 4
C_average_T = C_average.T

# make C symmetric
C = 0.5 * (C_average + C_average_T)




# compliance matrix (for compliance tensor), S = C^(-1)

Smat = np.linalg.inv(C)
print("\nCompliance matrix S = C^{-1}:\n")
print(Smat)

# calculate elastic moduli from C
# Bulk modulus (Voigt & Reuss)
KV = (C[0,0] + C[1,1] + C[2,2] + 2*(C[0,1] + C[1,2] + C[0,2])) / 9.0
KR = 1.0 / (Smat[0,0] + Smat[1,1] + Smat[2,2] + 2*(Smat[0,1] + Smat[1,2] + Smat[0,2]))
KH = 0.5*(KV + KR)

# Shear modulus (Voigt & Reuss)
GV = (
    C[0,0] + C[1,1] + C[2,2]
    - (C[0,1] + C[1,2] + C[0,2])
    + 3*(C[3,3] + C[4,4] + C[5,5])
) / 15.0

GR = 15.0 / (
    4*(Smat[0,0] + Smat[1,1] + Smat[2,2])
    - 4*(Smat[0,1] + Smat[1,2] + Smat[0,2])
    + 3*(Smat[3,3] + Smat[4,4] + Smat[5,5])
)

GH = 0.5 * (GV + GR)

E_young = 9*KH*GH / (3*KH + GH)
ν_poisson = (3*KH - 2*GH) / (2*(3*KH + GH))

print("Calculated Mechanical Properties")
print(f"Bulk modulus (Voigt)      KV = {KV:.4f}")
print(f"Bulk modulus (Reuss)      KR = {KR:.4f}")
print(f"Bulk modulus (Hill avg)   KH = {KH:.4f}\n")

print(f"Shear modulus (Voigt)     GV = {GV:.4f}")
print(f"Shear modulus (Reuss)     GR = {GR:.4f}")
print(f"Shear modulus (Hill avg)  GH = {GH:.4f}\n")

print(f"Young's modulus           E  = {E_young:.4f}")
print(f"Poisson ratio             ν  = {ν_poisson:.4f}")

# per-direction Young's Modulus

print("\nDirectional Young’s moduli from compliance:")
Ex = 1.0 / Smat[0,0]
Ey = 1.0 / Smat[1,1]
Ez = 1.0 / Smat[2,2]

print(f"E_x = {Ex:.4f}")
print(f"E_y = {Ey:.4f}")
print(f"E_z = {Ez:.4f}")

print("\nAll done.")




output_file = f"{analysis_directory}/elastic_constants_results.txt"

with open(output_file, "w") as f:
    f.write("=== Deformation matrices F (24x9) ===\n")
    # flatten each 3x3 matrix into 9 values per row
    F_flat = F.reshape(24, 9)
    np.savetxt(f, F_flat, fmt="%.10f")
    f.write("\n")

    f.write("=== Transposed deformation matrices F_T (24x9) ===\n")
    F_T_flat = F.transpose(0,2,1).reshape(24, 9)
    np.savetxt(f, F_T_flat, fmt="%.10f")
    f.write("\n")

    f.write("=== Green-Lagrange strain tensors (24x9) ===\n")
    strain_flat = strain_mats.reshape(24, 9)
    np.savetxt(f, strain_flat, fmt="%.10f")
    f.write("\n")

    f.write("=== Stress tensors (24x9) ===\n")
    stress_flat = stress_mats.reshape(24, 9)
    np.savetxt(f, stress_flat, fmt="%.10f")
    f.write("\n")

    f.write("=== Voigt notation E (strain, 24x6) ===\n")
    f.write(f"Voigt E: {E}\n")
    f.write("\n")

    f.write("=== Voigt notation S (stress, 24x6) ===\n")
    f.write(f"Voigt S: {S}\n")
    f.write("\n")

    f.write("=== Pseudoinverse of E (6x24) ===\n")
    f.write(f"Voigt E_pinv: {E_pinv}\n")
    f.write("\n")

    f.write("=== Elastic stiffness C (6x6) ===\n")
    np.savetxt(f, C, fmt="%.10f")
    f.write("\n")

    f.write("=== Compliance matrix Smat (6x6) ===\n")
    np.savetxt(f, Smat, fmt="%.10f")
    f.write("\n")

    f.write("=== Derived Elastic Moduli ===\n")
    f.write(f"Bulk modulus (Voigt)      KV = {KV:.10f}\n")
    f.write(f"Bulk modulus (Reuss)      KR = {KR:.10f}\n")
    f.write(f"Bulk modulus (Hill avg)   KH = {KH:.10f}\n")
    f.write(f"Shear modulus (Voigt)     GV = {GV:.10f}\n")
    f.write(f"Shear modulus (Reuss)     GR = {GR:.10f}\n")
    f.write(f"Shear modulus (Hill avg)  GH = {GH:.10f}\n")
    f.write(f"Young's modulus           E  = {E_young:.10f}\n")
    f.write(f"Poisson ratio             ν  = {ν_poisson:.10f}\n")
    f.write("\nDirectional Young's moduli:\n")
    f.write(f"E_x = {Ex:.10f}\n")
    f.write(f"E_y = {Ey:.10f}\n")
    f.write(f"E_z = {Ez:.10f}\n")

print(f"All results saved to {output_file}")
