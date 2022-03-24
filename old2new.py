import argparse
import numpy as np
np.set_printoptions(threshold=np.inf)
from collections import Counter

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

    print ("atom symbols ", atom_type_symbol)

## utility functions
def get_key(dictionary, val):
    for key, value in dictionary.items():
        if val == value:
            return key


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
    _, count_each_type_atoms = np.unique(atom_type_symbol, return_counts=True)
    print ("unique elements", unique_atoms)
    print ("indices ", indices)
    print ("count each atom type", count_each_type_atoms)
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

shells = {}
shells[0] = ['S']
shells[1] = ['X','Y','Z']
shells[2] = ['XX','XY','XZ','YY','YZ','ZZ']
shells[3] = ['XXX','XXY','XXZ','XYY','XYZ','XZZ','YYY','YYZ','YZZ','ZZZ']

dict_l_of_shells = {1:0, 3:1, 6:2, 10:3}


index_unique_atom = 0
basis_shell_ang_mom_unique = {}
dict_shell_counter = {}

for vals in dict_radial_pointers.values():
    pointers = [int(i) for i in vals]
    sorted_pointers = sorted(pointers)
    unique_shells, shell_count = np.unique(sorted_pointers, return_counts=True)
    print ("unique inds   ",   shell_count)
    print ("sorted pointers   ", sorted_pointers)

    # Counter of how many s,p,d,f shells for the given atom
    counter = np.zeros(4, dtype=int)
    _, temp_counter = np.unique(shell_count, return_counts=True)
    for i in range(len(temp_counter)):
        counter[i] = temp_counter[i]

    dict_shell_counter[unique_atoms[index_unique_atom]] = list(counter)
    print ("dict shell counter", dict_shell_counter)

    # Get and save the first line of the new bfinfo file
    first_line = []
    for ind, count in enumerate(counter):
        multiplier = get_key(dict_l_of_shells,ind)
        first_line.extend([count for i in range(multiplier)])

    first_line.extend([0 for i in range(len(first_line),20) ])
    # first line of bfinfo file obtained here.

    ### trial of something
    i = 0
    # get the length of individual sublist and stack if length matched with earlier sublist
    earlier_len = 1; temp_list = []; temp_list2 = []
    for group_size in shell_count:
        temp_list2.append(dict_l_of_shells[group_size])
        print ("grouping", sorted_pointers[i:i + group_size] )
        if len(sorted_pointers[i:i + group_size]) == earlier_len:
            temp_list.append(sorted_pointers[i:i + group_size])

        earlier_len = len(sorted_pointers[i:i + group_size])
        i += group_size
        # print ("grouping", [sorted_pointers[i:i + group_size] for i in range(0, len(sorted_pointers), group_size)] )
    print ("outside loop ", temp_list)
    basis_shell_ang_mom_unique[unique_atoms[index_unique_atom]] = temp_list2


    index_unique_atom += 1
print ("basis shell ang mom ", basis_shell_ang_mom_unique)

## Come outside the unique atoms loop

basis_shell_ang_mom = []; shell_counter_all_atoms = []
for i in atom_type_symbol:
    basis_shell_ang_mom.append(basis_shell_ang_mom_unique[i])
    shell_counter_all_atoms.append(dict_shell_counter[i])

print ("full list ang mom ", basis_shell_ang_mom)
print ("full shell counter list ", shell_counter_all_atoms)

# This part is for reshuffling to make the AO basis in the CHAMP's new own ordering
index_dict = {}; shell_representation = {}; bf_representation = {}
new_shell_representation = []
counter = 0; basis_per_atom = []
ind = 0; champ_ao_ordering = []
for atom_index in range(len(atom_type_symbol)):
    bfcounter = 1; basis_per_atom_counter = 0
    pindex = 0; dindex = 0; findex = 0
    for l in basis_shell_ang_mom[atom_index]:
        # run a small loop to reshuffle the shell ordering
        if l == 0:
            new_shell_representation.append(shells[l][0])
            champ_ao_ordering.append(ind)
            ind += 1

        local_p = np.zeros((3,shell_counter_all_atoms[atom_index][1]),dtype='U1')
        local_ind_p = np.zeros((3,shell_counter_all_atoms[atom_index][1]),dtype=int)
        if l == 1:
            pindex += 1; ind = champ_ao_ordering[-1] + 1
            for j in range(shell_counter_all_atoms[atom_index][1]):
                #loop over all 3 p orbitals
                for k in order[l]:
                    local_p[k,j] = shells[l][k]
                    local_ind_p[k,j] = ind
                    ind += 1

            if pindex == shell_counter_all_atoms[atom_index][1]:
                new_shell_representation.extend(list(local_p.flatten()))
                champ_ao_ordering.extend(list(local_ind_p.flatten()))

        local_d = np.zeros((6,shell_counter_all_atoms[atom_index][2]),dtype='U2')
        local_ind_d = np.zeros((6,shell_counter_all_atoms[atom_index][2]),dtype=int)
        if l == 2:
            dindex += 1; ind = champ_ao_ordering[-1] + 1
            for j in range(shell_counter_all_atoms[atom_index][2]):
                #loop over all 6 d orbitals
                for k in order[l]:
                    local_d[k,j] = shells[l][k]
                    local_ind_d[k,j] = ind
                    ind += 1

            if dindex == shell_counter_all_atoms[atom_index][2]:
                new_shell_representation.extend(list(local_d.flatten()))
                champ_ao_ordering.extend(list(local_ind_d.flatten()))


        local_f = np.zeros((10,shell_counter_all_atoms[atom_index][3]),dtype='U3')
        local_ind_f = np.zeros((10,shell_counter_all_atoms[atom_index][3]),dtype=int)
        if l == 3:
            findex += 1; ind = champ_ao_ordering[-1] + 1
            for j in range(shell_counter_all_atoms[atom_index][3]):
                #loop over all 10 f orbitals
                for k in order[l]:
                    local_f[k,j] = shells[l][k]
                    local_ind_f[k,j] = ind
                    ind += 1

            if findex == shell_counter_all_atoms[atom_index][3]:
                new_shell_representation.extend(list(local_f.flatten()))
                champ_ao_ordering.extend(list(local_ind_f.flatten()))

        # Get number of AO basis per atom
        for k in order[l]:
            shell_representation[counter] = shells[l][k]
            index_dict[counter] =  counter
            bf_representation[counter] = bfcounter
            counter += 1
            basis_per_atom_counter += 1
        bfcounter += 1
    basis_per_atom.append(basis_per_atom_counter)

print ("champ ao ordering: ", champ_ao_ordering)
print ("new_shell_representation: ", new_shell_representation)
print ("old_shell_representation: ", shell_representation.values())
print ("basis per atom: ", basis_per_atom)
print ("BF representation: ", bf_representation.values())


# ## Reorder orbitals according to the ordering of the CHAMP ordering
# reordered_mo_array = dict_mo["coefficient"][:,champ_ao_ordering]

# The next two arrays are needed for bfinfo file
reordered_bf_array = {k: bf_representation[k] for k in champ_ao_ordering}
reordered_bf_array_values = list(reordered_bf_array.values())
shell_representation_values = list(shell_representation.values())

print( "bf   ", reordered_bf_array_values)

accumumulated_basis_per_atom = np.cumsum(basis_per_atom)

start_index = 0
basis_pointer_per_atom = []
shell_reprensentation_per_atom = []
for i in range(len(basis_per_atom)):
    end_index = accumumulated_basis_per_atom[i]
    basis_pointer_per_atom.append(reordered_bf_array_values[start_index:end_index])
    shell_reprensentation_per_atom.append(shell_representation_values[start_index:end_index])
    start_index = end_index

print ("basis pointer per atom ", basis_pointer_per_atom)
print ("shell prepresentation per atom ", shell_reprensentation_per_atom)


## Write the new bfinfo file begins here ----------------
new_filename_bfinfo = "champ_v2_new_cartesian_" + args.filename_bfinfo
if new_filename_bfinfo is not None:
    if isinstance(new_filename_bfinfo, str):
        ## Write down a symmetry file in the new champ v2.0 format
        with open(new_filename_bfinfo, 'w') as file:

            # qmc bfinfo line printed below
            file.write("qmc_bf_info 1 \n")

            # pointers to the basis functions
            for i in indices:
                count_shells_per_atom = list(Counter(shell_reprensentation_per_atom[i]).values())
                # Write the number of types of shells for each unique atom
                for num in count_shells_per_atom:
                    file.write(f"{num} ")
                # Write down zeros for shells that are not present. Total shells supported are S(1) + P(3) + D(6) + F(10) = 20
                for rem in range(len(count_shells_per_atom), 20):
                    file.write(f"0 ")
                file.write(f"\n")

                # Write the pointers to the basis functions
                for pointer in basis_pointer_per_atom[i]:
                    file.write(f"{pointer} ")
                file.write(f"\n")
            file.write("end\n")
        file.close()

    else:
        raise ValueError
# all the bfinfo file information written to the file
