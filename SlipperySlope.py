"""LISP Calculator"""
import sys
import numpy as np

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
                    np.array([float(line[i]) for i in range(1,4)])
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
    return np.degrees(rad_angle)


print("""
    ================================
                 LISP   
    ================================
""")

mol=Molecule()
print("Load molecule...\tDONE!\n")

print("--> INPUT\n")
ring_idx=[int(i)-1 for i in input("Index for the ring...\t").split()]
M_idx=int(input("Index for the M atom...\t"))

##sanity check
print("\n--> SANITY CHECK\n")

try:
    print("Ring definition:\t\t", end="")
    if len(ring_idx)>=3:
        print("OK!")
    else:
        print("ERROR!")
        raise ValueError()
    d_list=[]

    for idx, val in enumerate(ring_idx):
        if val != ring_idx[-1]:
            d_list.append(
                dist(mol.xyz[val], mol.xyz[ring_idx[idx+1]])
                )
        else:
            d_list.append(
                dist(mol.xyz[val], mol.xyz[ring_idx[0]])
                )
    d_list_var=np.var(d_list)
    d_list_std=np.std(d_list)
    print("Variance:\t\t{:.3f}".format(d_list_var))

    print("Std. deviation:\t\t", end="")
    if d_list_std < 0.1:
        print("{:.3f} => OK!".format(d_list_std))
    else:
        print("{:.3f} => TOO HIGH!".format(d_list_std))
        raise ValueError()


except ValueError:
    print("\nWARNING: CALCULATION ABORTED DUE TO UNRELIABLE RESULTS")
    print("""
    ================================
           SANITY CHECK FAILED
    ================================
    """)
    quit()

###