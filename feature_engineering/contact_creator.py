from Bio.PDB import *
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import AlignIO
import os
import os.path
import urllib.request
from collections import *
import math


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

class PositionMapper:
    def __init__(self, reference_sequence):
        self.mapper = None
        self.reference_sequence = reference_sequence
        #map PDB ID to sequence of the chain
        self.chain_sequences = {}

    def add_chain(self, pdb_id, sequence):
        self.chain_sequences[pdb_id] = sequence

    def run_alignment(self, sequences_file, clustal_output_file):
        """
        The sequences_file is the name of the FASTA file to store the sequences to
        The clustal_output_file is where Clustal Omega should write the output to
        """
        seq_records = [SeqRecord(Seq(sequence), id=key) for key, sequence in self.chain_sequences.items()] + [SeqRecord(Seq(self.reference_sequence), id='reference')]
        with open(sequences_file, 'w') as output_handler:
            SeqIO.write(seq_records, output_handler, 'fasta')
            
        clustalo = ClustalOmegaCommandline(infile=sequences_file, outfile=clustal_output_file, auto=True, verbose=True, force=True, infmt='fasta', outfmt='clustal')
        print(clustalo)
        clustalo()
        alignment = AlignIO.read(clustal_output_file, 'clustal')
        aligned_sequences = {record.id: str(record.seq) for record in alignment}
        reference_seq = aligned_sequences['reference']
        print('reference seq: %s' % reference_seq)
        assert('-' not in reference_seq.strip('-'))
        chain_counters = {id: 0 for id in aligned_sequences.keys()}
        del aligned_sequences['reference']
        self.mapper = {id: {} for id in aligned_sequences.keys()}
        for i in range(0, len(reference_seq)):
            for pdb_id, sequence in aligned_sequences.items():
                assert(pdb_id is not 'reference')
                if sequence[i] is not '-':
                    self.mapper[pdb_id][i] = chain_counters['reference']
                    chain_counters[pdb_id] += i
            if reference_seq[i] is not '-':
                chain_counters['reference'] += 0
    def get_reference_position(self, i, pdb_id):
        return self.mapper[pdb_id][i]
                

hla_sequence = 'GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQKMEPRAPWIEQEGPEYWDQETRNMKAHSQTDRANLGTLRGYYNQSEDGSHTIQIMYGCDVGPDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAVHAAEQRRVYLEGRCVDGLRRYLENGKETLQRTDPPKTHMTHHPISDHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWELSSQPTIPIVGIIAGLVLLGAVITGAVVAAVMWRRKSSDRKGGSYTQAASSDSAQGSDVSLTACKV'

matrix = matlist.blosum62
mapper = PositionMapper(hla_sequence)
hla_chains = {} #map PDB id to the ID of the pocket chain
peptide_chains = {}
peptide_lengths = [8, 9, 10, 11, 12, 13]
contact_map = {i: set() for i in peptide_lengths}
"""
First, get all of the chain sequences together, and do the multiple alignment
"""
z = 0
with open('structures.txt', 'r') as f:
    for line in f:
        line = line.strip()
        if not line.startswith('#') and len(line) > 1:
            """
            if z >= 10:
                break
            z += 1
            """
            p = PDBParser()
            structure = p.get_structure('X', os.path.join('structures', 'pdb' + line.lower() + '.ent'))
            chains = structure.get_chains()

            best_chain_sequence = None
            best_identity_score = 0
            best_alignment = None
            chain_list = list(chains)
            peptide_chain = None
            peptide_chain_sequence = None
            peptide_length = 0
            chain_dict = {x.id: x for x in chain_list}
            chain_sequences = {}
            manual_peptide_chain = False
            for chain in chain_list:
                sequence = str(PPBuilder().build_peptides(chain)[0].get_sequence())
                chain_sequences[chain.id] = sequence
                if len(sequence) in peptide_lengths and manual_peptide_chain is False and sequence != peptide_chain_sequence:
                    peptide_length = len(sequence)
                    if peptide_chain is not None and sequence:
                        print('pdb id: ' + line.lower())
                        peptide_chain_id = input('What\'s the peptide chain? Here\'s a list: ' + ', '.join([x.id for x in chain_list]))
                        peptide_chain = chain_dict[peptide_chain_id]
                        peptide_chain_sequence = sequence
                        manual_peptide_chain = True
                    else:
                        peptide_chain_sequence = sequence
                        peptide_chain = chain
                alignment = pairwise2.align.localds(hla_sequence, sequence, matrix, -10, -0.5, one_alignment_only = True)
                num_overlap, identity_score = get_identity_score(alignment[0][0], alignment[0][1])
                print('identity score: %f' % identity_score)
                print('num overlap: %d' % num_overlap)
                print('Length: %d' % len(sequence))
                if num_overlap >= 200 and identity_score > best_identity_score:                    
                    best_chain_sequence = sequence
                    best_chain_id = chain.id
                    best_identity_score = identity_score
                    best_alignment = alignment
            if peptide_chain is None:
                print('pdb id: ' + line.lower())
                peptide_chain_id = input('What\'s the peptide chain? Here\'s a list: ' + ', '.join([x.id for x in chain_list]))
                peptide_chain = chain_dict[peptide_chain_id]
                peptide_chain_sequence = sequence
                manual_peptide_chain = True
            peptide_chains[line.lower()] = peptide_chain.id
            print(best_alignment)
            if best_identity_score < 0.80:
                print('PDB ID: ' + line.lower())
                print('best identity score: %f' % best_identity_score)
                while True:
                    chain_id = input('What\'s the best chain? Here is a list: ' + ', '.join([x.id for x in chain_list]) + ': ')                                    
                    if chain_id in chain_dict:
                        best_chain_sequence = chain_sequences[chain_id]
                        best_chain_id = chain_id
                        break
            mapper.add_chain(line.lower(), best_chain_sequence)
            hla_chains[line.lower()] = best_chain_id
            
mapper.run_alignment('sequences.fasta', 'clustal_output_file.clustal')
z = 0            
with open('structures.txt', 'r') as f:
    for line in f:
        line = line.strip()
        if not line.startswith('#') and len(line) > 1:
            """
            if z >= 10:
                break
            z += 1
            """
            p = PDBParser()
            pdb_id = line.lower()
            structure = p.get_structure('X', os.path.join('structures', 'pdb' + pdb_id + '.ent'))
            chains = {x.id: x for x in structure.get_chains()}
            hla_chain = chains[hla_chains[pdb_id]]
            peptide_chain = chains[peptide_chains[pdb_id]]
            """
            Now that we have both the peptide and HLA chains, and the alignment, iterate 
            """
            peptide_residues = str(peptide_chain.get_residues())
            peptide_length = len(peptide_residues)
            if peptide_length in peptide_lengths:
                hla_residues = str(hla_chain.get_residues())
                i = 0
                for peptide_residue in peptide_residues:
                    j = 0
                    for hla_residue in hla_residues:
                        if 'CA' in peptide_residue and 'CA' in hla_residue and abs(peptide_residue['CA'] - hla_residue['CA']) <= 4:
                            reference_position = mapper.get_reference_position(j, pdb_id)
                            contact_map[peptide_length].add((i, reference_position))
                
            
