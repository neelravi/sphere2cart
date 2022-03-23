import argparse
import numpy as np
np.set_printoptions(threshold=np.inf)

# Instantiate the parser
parser = argparse.ArgumentParser(description='Python Converter for conversion of mixed type MOs to all-cartesian MOs for CHAMP code')

# Required positional argument
parser.add_argument("--lcao", "-s", "--orb", dest='filename_lcao', type=str,
                    help='Required: Filename (including extension) of .lcao or .orb file in the old format')

# Required positional argument
parser.add_argument("--bfinfo", "-r", "--bf", dest='filename_bfinfo', type=str,
                    help='Required: Filename (including extension) of .bfinfo file in the old format')

# Required positional argument
parser.add_argument("--geom", "-g", "--xyz", dest='filename_geom', type=str,
                    help='Required: Filename of .geom file in the old format')

args = parser.parse_args()

print("Filenames parsed are :")
print(" lcao file old   :: ", args.filename_lcao)
print(" bfinfo file old :: ", args.filename_bfinfo)
print(" geom file old   :: ", args.filename_geom)


### Read the lcao file first

with open(args.filename_lcao) as f1:
    lines = f1.readlines()

    iorb = 0; mocoeffs = []
    for line in lines:
        # Skip the comments
        if line.startswith('#'):
            continue
        if line.startswith('end'):
            continue
        # Read the number of basis and number of coeffs
        if line.startswith('lcao'):
            nbasis = int(line.split()[1])
            ncoeff = int(line.split()[2])
            continue

        temp = [float(i) for i in list(line.split()) ]
        mocoeffs.append(temp)
        iorb += 1

# print (mocoeffs)


### Read the geometry file
with open(args.filename_geom) as f2:
    lines = f2.readlines()

    dict_atom_type = {}
    atom_type = []
    for line_num, line in enumerate(lines):
        # Skip the comments
        if line.startswith('&atoms'):
            nctype = int(line.split()[2])
            natoms = int(line.split()[4])
            continue
        if line.startswith('&atom_types'):
            ntokens = len(line.split()) - 1
            for i in range(0,ntokens,2):
                dict_atom_type[int(line.split()[i+1])] = line.split()[i+2]
            continue

        if line.startswith('geometry'):
            coord_block_start = line_num + 1

    print ("nstoms ntypes and other data")
    print (natoms, nctype, dict_atom_type)

    for i in range(coord_block_start, coord_block_start+natoms):
        atom_type.append(lines[i].split()[3])

    print (atom_type)
    atom_type_symbol = []
    for i in atom_type:
        atom_type_symbol.append(dict_atom_type[int(i)])

    print (atom_type_symbol)

## utility function to get the patterns in bfinfo file
def guess_seq_len(seq):
    guess = 1
    max_len = int(len(seq) / 2)
    for x in range(2, max_len):
        print ("all ", seq[0:x])
        if seq[0:x] == seq[x:2*x] :
            print ("squence ", seq[0:x])
            return x
    return guess


### Read the bfinfo file
with open(args.filename_bfinfo) as f3:
    lines = f3.readlines()

    for line_num, line in enumerate(lines):
        print (line, line_num)
        # Skip the comments
        # Skip the comments
        if line.startswith('qmc_bf_info'):
            coord_block_start = line_num + 1
        if line.startswith('end'):
            continue

    unique_atoms, indices = np.unique(atom_type_symbol, return_index=True)

    print ("unique elements", unique_atoms)
    print ("indices ", indices)
    num_unique_atoms = len(unique_atoms)

    dict_num_per_shell = {} #only the odd numbered rows of data
    for i in range(coord_block_start, coord_block_start+2*num_unique_atoms,2):
        dict_num_per_shell[i] = lines[i].split()

    print ("dict num per shell", dict_num_per_shell)


    dict_radial_pointers = {} #only the even numbered rows of data
    for i in range(coord_block_start+1, coord_block_start+2*num_unique_atoms,2):
        dict_radial_pointers[i] = lines[i].split()

    print ("radial pointers", dict_radial_pointers)



### Writing part starts here!

## Get the bfinfo file with cartesian ordering of shell pointers
# Cartesian ordering
order = [[0],
         [0, 1, 2],
         [0, 1, 2, 3, 4, 5],
         [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]]

dict_l_of_shells = {1:0, 3:1, 6:2, 10:3}



for vals in dict_radial_pointers.values():
    bfcounter = 1; bf_representation = {}
    counter0 = 0; basis_per_atom_counter = 0
    pointers = [int(i) for i in vals]
    sorted_pointers = sorted(pointers)
    unique_shells, shell_count = np.unique(sorted_pointers, return_counts=True)
    print ("unique inds   ", shell_count)

    # Counter of how many p,d,f shells for the given atom
    _, counter = np.unique(shell_count, return_counts=True)
    print ("counter of shells per atom ", counter)
    if len(counter) == 4:
        local_ind_p = np.zeros((3,counter[1]),dtype=int)
        local_ind_d = np.zeros((3,counter[2]),dtype=int)
        local_ind_f = np.zeros((3,counter[3]),dtype=int)
    elif len(counter) == 3:
        local_ind_p = np.zeros((3,counter[1]),dtype=int)
        local_ind_d = np.zeros((3,counter[2]),dtype=int)
    elif len(counter) == 2:
        local_ind_p = np.zeros((3,counter[1]),dtype=int)


    for j in shell_count:
        # local np array for reshuffling later
        print ("the shell count", j)
        l = dict_l_of_shells[j]
        print ("the l", l)
        print (counter0)
        for k in order[l]:
            bf_representation[counter0] = bfcounter
            counter0 += 1
            basis_per_atom_counter += 1
        bfcounter += 1

    print ("bf repre ", bf_representation)