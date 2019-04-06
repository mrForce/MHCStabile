from Bio.PDB import *
import os
import collections
import math
cwd = os.getcwd()
peptide_list = []
bad = open('bad.txt', 'w')
with open('structures.txt', 'r') as f:
    for line in f:
        line = line.strip()
        if not line.startswith('#')  and len(line) > 1:
            pdbl = PDBList()
            pdbl.retrieve_pdb_file(line, pdir='structures', file_format = 'pdb')
            p = PDBParser()
            structure = p.get_structure('X', os.path.join('structures', 'pdb' + line.lower() + '.ent'))
            ppb = PPBuilder()
            peptides = sorted([str(x.get_sequence()) for x in ppb.build_peptides(structure)], key=lambda x:abs(len(x) -9))
            peptide_list.append(peptides[0])
            if len(peptides[0]) < 8:
                bad.write('peptide: %s, id: %s\n' % (peptides[0], line.lower()))
                bad.write(', '.join(peptides))
c = collections.Counter([len(x) for x in peptide_list])
for k,v in c.items():
    print('length: %d, number: %d' % (k, v))
bad.close()
