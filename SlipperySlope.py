"""Slippage Calculator"""
import sys
import numpy as np

__version__ = "0.1.3"


class Molecule:
    """Molecule loading"""
    def __init__(self):
        try:
            with open(sys.argv[1]) as file:
                lines = file.readlines()
            self.natoms = int(lines[0])
            self.element = []
            self.xyz = []
            for line in lines[2:]:
                line = line.split()
                self.element.append(line[0])
                self.xyz.append(
                    np.array([float(line[i]) for i in range(1, 4)])
                )

        except:
            print('Error during importing molecule. Check XMol format!')
            quit()


def dist(atom_a, atom_b):
    """Calulate the distance between two atoms"""
    return np.linalg.norm(atom_a-atom_b)


def ang(atom_a, atom_b, atom_c):
    """Calulate the angle given three atoms"""
    b_a = atom_a - atom_b
    b_c = atom_c - atom_b
    rad_angle = np.arccos(
        np.dot(b_a, b_c) / (dist(atom_a, atom_b) * dist(atom_b, atom_c))
    )
    return rad_angle


def lisp(ring_coor, m_coor):
    """Calculate the Label Independednt Slippage Parameter"""
    ring_ctd = np.mean(ring_coor, axis=0)
    d_ring_ctd_m = dist(ring_ctd, m_coor)
    d_list_norm = []
    for idx, val in enumerate(ring_coor):
        if idx != len(ring_coor)-1:
            d_ring_pair_mid = np.mean(
                [ring_coor[idx], ring_coor[idx+1]], axis=0)
            print("<{}-{}>...\t\t{:.6f}  {:.6f}  {:.6f}".format(
                RING_IDX[idx] + 1,
                RING_IDX[idx+1] + 1,
                d_ring_pair_mid[0],
                d_ring_pair_mid[1],
                d_ring_pair_mid[2]
                ))
        else:
            d_ring_pair_mid = np.mean([ring_coor[idx], ring_coor[0]], axis=0)
            print("<{}-{}>...\t\t{:.6f}  {:.6f}  {:.6f}".format(
                RING_IDX[idx] + 1,
                RING_IDX[0] + 1,
                d_ring_pair_mid[0],
                d_ring_pair_mid[1],
                d_ring_pair_mid[2]
                ))

        theta = ang(d_ring_pair_mid, ring_ctd, m_coor)
        d_list_norm.append(abs(d_ring_ctd_m*np.sin(theta-np.pi/2)))

    return np.mean(d_list_norm), ring_ctd, d_ring_ctd_m

print("""
    ================================
             SlipperySlope
    ====================== Mk. {}
""".format(__version__))

MOL = Molecule()
print("Load molecule...\tDONE!\n")

print("--> INPUT\n")
try:
    RING_IDX = [int(i)-1 for i in input("Index for the ring...\t").split()]
    M_IDX = int(input("Index for the M atom...\t"))-1
except ValueError as val_err:
    print("The input index is not valid: ", val_err)
    quit()

print("\n--> SANITY CHECK\n")

try:
    print("Ring definition...\t\t", end="")
    if len(RING_IDX) >= 3:
        print("OK!")
    else:
        print("ERROR!")
        raise ValueError()

    D_LIST = []
    for idx, val in enumerate(RING_IDX):
        if idx != len(RING_IDX)-1:
            D_LIST.append(
                dist(MOL.xyz[val], MOL.xyz[RING_IDX[idx+1]])
            )
        else:
            D_LIST.append(
                dist(MOL.xyz[val], MOL.xyz[RING_IDX[0]])
            )
    D_LIST_VAR = np.var(D_LIST)
    D_LIST_STD = np.std(D_LIST)
    print("Variance...\t\t{:.3f}".format(D_LIST_VAR))
    print("Std. deviation...\t", end="")
    if D_LIST_STD < 0.1:
        print("{:.3f}\tOK!".format(D_LIST_STD))
    else:
        print("{:.3f}\tTOO HIGH!".format(D_LIST_STD))
        raise ValueError()

except ValueError:
    print("\nWARNING: CALCULATION ABORTED DUE TO UNRELIABLE RESULTS")
    print("""
    ================================
           SANITY CHECK FAILED
    ================================
    """)
    quit()

print("\n--> LISP CALCULATION\n")
RING_COOR = [MOL.xyz[i] for i in RING_IDX]
LISP, RING_CTD, D_RING_CTD_M = lisp(RING_COOR, MOL.xyz[M_IDX])

print("\nRing centroid...\t{:.6f}  {:.6f}  {:.6f}".format(RING_CTD[0], RING_CTD[1], RING_CTD[2]))
print("Dist. centroid-M...\t{:.3f}".format(D_RING_CTD_M))
print("LISP...\t\t\t{:.3f}".format(LISP))

print("\n--> BASOLO CALCULATION\n")
D_LIST_M_RING = [dist(MOL.xyz[M_IDX], i) for i in RING_COOR]

L_SORT = sorted((e, i) for i, e in enumerate(D_LIST_M_RING))
LONG_1 = L_SORT[-1]
LONG_2 = L_SORT[-2]

print("Longer C-M vicinity check...\t", end="")

if LONG_2[1] == LONG_1[1] + 1:
    if LONG_1[1] - 1 == -1:
        ADJ_1_IDX = len(L_SORT) - 1
    else:
        ADJ_1_IDX = LONG_1[1] - 1

    if LONG_2[1] + 1 == len(L_SORT):
        ADJ_2_IDX = 0
    else:
        ADJ_2_IDX = LONG_2[1] + 1
    print("OK!")
elif LONG_2[1] == LONG_1[1] - 1:
    if LONG_2[1] - 1 == -1:
        ADJ_2_IDX = len(L_SORT) - 1
    else:
        ADJ_2_IDX = LONG_2[1] - 1

    if LONG_1[1] + 1 == len(L_SORT):
        ADJ_1_IDX = 0
    else:
        ADJ_1_IDX = LONG_1[1] + 1
    print("OK!")
elif LONG_1[1] + 1 == len(L_SORT) and LONG_2[1] == 0:
    ADJ_1_IDX = LONG_1[1] - 1
    ADJ_2_IDX = 1
    print("OK!")
elif LONG_2[1] + 1 == len(L_SORT) and LONG_1[1] == 0:
    ADJ_1_IDX = LONG_1[1] + 1
    ADJ_2_IDX = LONG_2[1] - 1
    print("OK!")
else:
    print("FAILED!\nThe two longer C-M carbon atoms are not consecutive.\n Basolo Delta cannot be calcualted.")
    BASOLO_DELTA = "N/D"

BASOLO_DELTA = ((LONG_1[0]+LONG_2[0])-(D_LIST_M_RING[ADJ_1_IDX]+D_LIST_M_RING[ADJ_2_IDX]))/2

if BASOLO_DELTA == "N/D":
    print("Basolo Delta...\t\t\t{}".format(BASOLO_DELTA))
else:
    print("Long Bond:\t\t\tC{}-M, C{}-M".format(RING_IDX[LONG_1[1]]+1, RING_IDX[LONG_2[1]]+1))
    print("Neighboring Bond:\t\tC{}-M, C{}-M".format(RING_IDX[ADJ_1_IDX]+1, RING_IDX[ADJ_2_IDX]+1))
    print("Basolo Delta...\t\t\t{:.3f}".format(BASOLO_DELTA))

print("""
    ================================
           NORMAL TERMINATION
    ================================
""")
