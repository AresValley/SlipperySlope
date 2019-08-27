"""LISP Calculator"""
import sys
import numpy as np

__version__ = "0.1.1"


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
    ba = atom_a - atom_b
    bc = atom_c - atom_b
    rad_angle = np.arccos(
        np.dot(ba, bc) / (dist(atom_a, atom_b) * dist(atom_b, atom_c))
    )
    return rad_angle


def lisp(ring_coor, M_coor):
    """Calculate the Label Independednt Slippage Parameter"""
    ring_ctd = np.mean(ring_coor, axis=0)
    d_ring_ctd_M = dist(ring_ctd, M_coor)
    d_list_norm = []
    for idx, val in enumerate(ring_coor):
        if idx != len(ring_coor)-1:
            d_ring_pair_mid = np.mean(
                [ring_coor[idx], ring_coor[idx+1]], axis=0)
            print("<{}-{}>...\t\t{:.6f}  {:.6f}  {:.6f}".format(
                ring_idx[idx] + 1,
                ring_idx[idx+1] + 1,
                d_ring_pair_mid[0],
                d_ring_pair_mid[1],
                d_ring_pair_mid[2]
                ))
        else:
            d_ring_pair_mid = np.mean([ring_coor[idx], ring_coor[0]], axis=0)
            print("<{}-{}>...\t\t{:.6f}  {:.6f}  {:.6f}".format(
                ring_idx[idx] + 1,
                ring_idx[0] + 1,
                d_ring_pair_mid[0],
                d_ring_pair_mid[1],
                d_ring_pair_mid[2]
                ))

        theta = ang(d_ring_pair_mid, ring_ctd, M_coor)
        d_list_norm.append(abs(d_ring_ctd_M*np.sin(theta-np.pi/2)))

    return np.mean(d_list_norm), ring_ctd, d_ring_ctd_M

print("""
    ================================
             SlipperySlope
    ====================== Mk. {}
""".format(__version__))

MOL = Molecule()
print("Load molecule...\tDONE!\n")

print("--> INPUT\n")
try:
    ring_idx = [int(i)-1 for i in input("Index for the ring...\t").split()]
    M_idx = int(input("Index for the M atom...\t"))-1
except ValueError as val_err:
    print("The input index is not valid: ", val_err)
    quit()

print("\n--> SANITY CHECK\n")

try:
    print("Ring definition...\t\t", end="")
    if len(ring_idx) >= 3:
        print("OK!")
    else:
        print("ERROR!")
        raise ValueError()

    d_list = []
    for idx, val in enumerate(ring_idx):
        if idx != len(ring_idx)-1:
            d_list.append(
                dist(MOL.xyz[val], MOL.xyz[ring_idx[idx+1]])
            )
        else:
            d_list.append(
                dist(MOL.xyz[val], MOL.xyz[ring_idx[0]])
            )
    d_list_var = np.var(d_list)
    d_list_std = np.std(d_list)
    print("Variance...\t\t{:.3f}".format(d_list_var))
    print("Std. deviation...\t", end="")
    if d_list_std < 0.1:
        print("{:.3f}\tOK!".format(d_list_std))
    else:
        print("{:.3f}\tTOO HIGH!".format(d_list_std))
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
ring_coor = [MOL.xyz[i] for i in ring_idx]
LISP, ring_ctd, d_ring_ctd_M = lisp(ring_coor, MOL.xyz[M_idx])

print("\nRing centroid...\t{:.6f}  {:.6f}  {:.6f}".format(ring_ctd[0], ring_ctd[1], ring_ctd[2]))
print("Dist. centroid-M...\t{:.3f}".format(d_ring_ctd_M))
print("LISP...\t\t\t{:.3f}".format(LISP))

print("""
    ================================
           NORMAL TERMINATION
    ================================
""")
