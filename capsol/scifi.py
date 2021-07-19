
import numpy as np
from collections import OrderedDict, defaultdict


def read_file(fname):
    with open(fname, "r") as f:
        return f.read().splitlines()

def write_file(fname, string):
    with open(fname, 'w') as fh:
        if isinstance(string, str):
            fh.write(string)
        else:
            try:
                fh.write("\n".join(string))
            except:
                print(f"Cannot iterate over {string}")
                raise


# file = read_file("mgo_mt.INN")

# lines = """

# .<Toleranc>
# printing medium
# shells moveall
# xyz-file write
# traj_fmt 2
# movie
# g98-file write
# max_opt   200
# optim bfgs
# shake 0.100000E-01
# en_tol 0.100000E-05
# es_tol 0.100000E-02
# frc_tol 0.100000E-02
# crd_tol 0.100000E-02
# do_lateral
# .<End>

# """.split("\n")


# expected_output = dict(Toleranc={"printing": "medium", "shells": "moveall"})

def process_scifi_input(lines):
    # Remove blank lines
    lines_no_blank = [line for line in lines if line != ""]

    out = dict()
    in_block = False
    key = None
    for line in lines_no_blank:
        # This signals that we should start a block...
        if line[0] == "." and not in_block:
            in_block = True
            key = line[2:-1] # Strip leading . and < and trailing >
            out[key] = []
            continue
        
        # End a block
        if line == ".<End>" and in_block:
            key = None
            in_block = False
            continue

        # With a block, process things on whitespace...
        if in_block:
            out[key].append(line)

    return out

def process_coord(lines):
    group = "all"
    out = defaultdict(lambda : [])
    for line in lines:
        if line in ["frozen-tip", "tip-atoms", "tip-buffer", "cluster-atoms", "cluster-buffer"]:
            group = line
        else:
            out[group].append(line)
    
    return out

def xyz_file(atom_group, comment=""):
    atom_xyz = [" ".join(line.strip().split()[:4]) for line in atom_group]
    N_atoms = len(atom_xyz)
    all_lines = [str(N_atoms), comment]
    all_lines.extend(atom_xyz)
    return "\n".join(all_lines)

def get_atoms(lines):
    return [line.strip().split()[0] for line in lines]

def get_coords(lines):
    return np.array([[float(x) for x in line.strip().split()[1:4]] for line in lines])

from dataclasses import dataclass




@dataclass
class XYZFile:
    atoms: list
    coords: np.ndarray

    def write_str(self, comment="") -> str:
        N_atoms = len(self.atoms)
        atom_xyz = [atom+" "+" ".join([str(x) for x in coord]) for atom, coord in zip(self.atoms, self.coords)]
        all_lines = [str(N_atoms), comment]
        all_lines.extend(atom_xyz)
        return "\n".join(all_lines)



def get_xyz(fname):
    file = read_file(fname)
    atoms = []
    coords = []
    for line in file[2:]:
        items = line.split()
        if line != "":
            atoms.append(items[0])
            coords.append([float(x) for x in items[1:]])
        
    coords = np.array(coords)
    return XYZFile(atoms, coords)
            



# s = process_scifi_input(read_file("mgo_mt-Si.INN"))

# out = process_coord(s['Coord'])

# write_file("mgo_mt-Si.frozen-tip.xyz", xyz_file(out['frozen-tip']))

# write_file("mgo_mt-Si.tip-atoms.xyz", xyz_file(out['tip-atoms']))

# Si81 #81 -2.482 3.415 -0.573



# cluster_atom_labels = [line.strip().split()[0] for line in cluster_atoms]

# cluster_data = np.array([[float(x) for x in line.strip().split()[1:7]] for line in cluster_atoms])

# cluster_xyz = [" ".join(line.strip().split()[:4]) for line in cluster_atoms]

# cluster_coords = cluster_data[:, :3]
# cluster_frozen = np.array(cluster_data[:, 3:], dtype=int)

# print("\n".join(cluster_xyz))





# # Output the file
# open("mgo-cluster-atoms-new.xyz", 'w').write(xyz_file(cluster_atoms))

# new_lines = open("mgo-cluster-atoms-new-cut.xyz", 'r').read().splitlines()[2:]


# n_atom_labels = get_atoms(new_lines)
# n_atom_coords = get_coords(new_lines)


# inds = [np.argmin(np.sum(abs(coords - cluster_coords), axis=1)) for coords in n_atom_coords]

# n_atom_frozen = np.array([cluster_frozen[i] for i in inds])

# new_lines = []
# for atom, coords, frozen in zip(n_atom_labels, n_atom_coords, n_atom_frozen):
#     line = [atom]
#     line.extend([str(x) for x in coords])
#     line.extend([str(x) for x in frozen])
#     new_lines.append(line)


# out['cluster-atoms'] = [" ".join(line) for line in new_lines]

# cluster_lines = [".<Coord>"]
# for key, val in out.items():
#     cluster_lines.append(key)
#     cluster_lines.extend(val)

# cluster_lines.append(".<End>")

# cluster_str = "\n".join(cluster_lines)

# open("mgo-cluster.Coord", 'w').write(cluster_str)

# print(n_atom_coords)
# print(cluster_coords)




# print(out)
# print(out.keys())
# print("\n".join(out['frozen-tip']))
# print("\n".join(out['tip-atoms']))