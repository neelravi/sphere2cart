import argparse
import numpy as np
np.set_printoptions(threshold=np.inf)

# Instantiate the parser
parser = argparse.ArgumentParser(description='Python Converter for conversion of mixed type MOs to all-cartesian MOs for CHAMP code')

# Required positional argument
parser.add_argument("--lcao", dest='filename_lcao', type=str,
                    help='Required: Filename (including extension) of .lcao or .orb file in the old format')

# Required positional argument
parser.add_argument("--bfinfo", dest='filename_bfinfo', type=str,
                    help='Required: Filename (including extension) of .bfinfo file in the old format')

# Required positional argument
parser.add_argument("--geom", dest='filename_geom', type=str,
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
    for line in lines:
        # Skip the comments
        if line.startswith('&atoms'):
            nctype = int(line.split()[2])
            natoms = int(line.split()[4])
            continue
        if line.startswith('&atom_types'):
            ntokens = len(line.split()) - 1
            for i in range(1,ntokens+1):
                dict_atom_type[int(line.split()[i*2])] = line.split()[i*2]
            continue
        if line.startswith('geometry'):

            continue
        # Read the number of basis and number of coeffs
        if line.startswith('lcao'):
            nbasis = int(line.split()[1])
            ncoeff = int(line.split()[2])
            continue

        temp = [float(i) for i in list(line.split()) ]
        mocoeffs.append(temp)
        iorb += 1
