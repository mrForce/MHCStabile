from Bio.PDB import *
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio import SeqIO
import os
import os.path
import urllib.request
from collections import *

def get_uniprot(uniprot_id):
    if ':' in uniprot_id:
        #then of the form UNP:P01891
        uniprot_id = uniprot_id.split(':')[1]

    filename = uniprot_id + '.xml'
    file_location = os.path.join('uniprot', filename)
    if not os.path.isfile(file_location):
        handler = urllib.request.urlretrieve('http://www.uniprot.org/uniprot/' + filename, file_location)
    record = SeqIO.read(file_location, 'uniprot-xml')
    return str(record.seq)

def get_identity_score(seq_one, seq_two):
    assert(len(seq_one) == len(seq_two))
    num_overlap = 0
    num_matches = 0
    for i in range(0, len(seq_one)):
        if seq_one[i] != '-' and seq_two[i] != '-':
            num_overlap += 1
            if seq_one[i] == seq_two[i]:
                num_matches += 1
    return (num_overlap, num_matches*1.0/num_overlap)

def get_sequence_

hla_sequence = 'GALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQKMEPRAPWIEQEGPEYWDQETRNMKAHSQTDRANLGTLRGYYNQSEDGSHTIQIMYGCDVGPDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAVHAAEQRRVYLEGRCVDGLRRYLENGKETLQRTDPPKTHMTHHPISDHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWELSSQPTIPIVGIIAGLVLLGAVITGAVVAAVMWRRKSSDRKGGSYTQAASSDSAQGSDVSLTACKV'

matrix = matlist.blosum62

position_dict = defaultdict(list)


with open('structures.txt', 'r') as f:
    for line in f:
        line = line.strip()
        if not line.startswith('#') and len(line) > 1:
            p = PDBParser()
            structure = p.get_structure('X', os.path.join('structures', 'pdb' + line.lower() + '.ent'))
            chains = structure.get_chains()

            best_chain = None
            best_identity_score = 0
            best_alignment = None
            chain_list = list(chains)
            peptide_chain = None
            peptide_length = 0
            chain_dict = {x.id: x for x in chain_list}
            manual_peptide_chain = False
            for chain in chain_list:
                sequence = PPBuilder().build_peptides(chain)[0].get_sequence()
                if len(sequence) in [8, 9, 10, 11, 12, 13] and manual_peptide_chain is False:
                    peptide_length = len(sequence)
                    if peptide_chain is not None:
                        print('pdb id: ' + line.lower())
                        peptide_chain_id = input('What\'s the peptide chain? Here\'s a list: ' + ', '.join([x.id for x in chain_list]))
                        peptide_chain = chain_dict[peptide_chain_id]
                        manual_peptide_chain = True
                    else:
                        peptide_chain = chain
                alignment = pairwise2.align.localds(hla_sequence, sequence, matrix, -10, -0.5, one_alignment_only = True)
                num_overlap, identity_score = get_identity_score(alignment[0][0], alignment[0][1])
                print('identity score: %f' % identity_score)
                print('num overlap: %d', num_overlap)
                if num_overlap >= 200 and identity_score > best_identity_score:                    
                    best_chain = chain
                    best_identity_score = identity_score
                    best_alignment = alignment
            print(best_alignment)
            if best_identity_score < 0.80:
                print('PDB ID: ' + line.lower())
                print('best identity score: %f' % best_identity_score)
                while True:
                    chain_id = input('What\'s the best chain? Here is a list: ' + ', '.join([x.id for x in chain_list]) + ': ')                                    
                    if chain_id in chain_dict:
                        best_chain = chain_dict[chain_id]
                        break
            """
            Now that we have both the peptide and HLA chains, and the alignment, iterate 
            """
            peptide_residues = peptide_chain.get_residues()
            hla_residues = best_chain.get_residues()
            for peptide_residue in peptide_residues:
                for hla_residue in hla_residues:
                    if 'CA' in peptide_residue and 'CA' in hla_residue:
                        print(peptide_residue['CA'] - hla_residue['CA'])
                
                
                
                
            """
            best_alignment_score = 0
            best_alignment = None
            for record in SeqIO.parse(os.path.join('structures', 'pdb' + line.lower() + '.ent'), 'pdb-seqres'):
                if 'UNP' in record.dbxrefs[0]:
                    sequence = get_uniprot(record.dbxrefs[0])
                    
                    
                    
                if record.annotations['chain'] == 'A':
                    print(record.id)
                    print('database refs: ')
                    print(record.dbxrefs)
                    has_a = True
                    assert('UNP' in record.dbxrefs[0])
            assert(has_a)
            """
            
            #print('line: %s' % line)
            #for chain in chains:
            #    sequence = PPBuilder().build_peptides(chain)[0].get_sequence()
            #    alignment = pairwise2.align.localxx(sequence, hla_sequence)[0]
                
                
            #print('line: %s, chain lengths: %s' % (line, ', '.join([str(len(list(x.get_residues()))) for x in chains])))
            
            #ppb = PPBuilder()
            #peptides = sorted([str(x.get_sequence()) for x in ppb.build_peptides(structure)], key=lambda x:abs(len(x) -9))
            
