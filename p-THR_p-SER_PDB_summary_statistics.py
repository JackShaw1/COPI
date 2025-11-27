from Bio.PDB import MMCIFParser, NeighborSearch
import os
from collections import Counter
import csv

"""

This script extracts data from each phosphoserine and phosphothreonine deposited to the PDB.
Each of these phosphosites is within 6 Angstroms of exactly one foreign protein chain. In the output
file, all_SEP_TPO_in_PDB_nov_26_2025.csv, the pdb id, chain id for phosphorylated protein, residue index of phosphorylated
amino acid (w/ the first amino acid in the chain starting as index=1), the chain id for the lone binding
protein, and binding amino acid index (again starting with index=1 at N terminus) are written across one
line per phosphosite. 

These phosphosites were used to train The Phospho-Interaction Library's COPI model. Known stimulatory sites 
from the set were also identified using these data.

Requirement: pdb_structures_sep_tpo_nov_26_2025 data folder must be present in cwd where this script is executed!

"""

amino_acids = [
    "ALA",  # Alanine
    "ARG",  # Arginine
    "ASN",  # Asparagine
    "ASP",  # Aspartic acid
    "CYS",  # Cysteine
    "GLN",  # Glutamine
    "GLU",  # Glutamic acid
    "GLY",  # Glycine
    "HIS",  # Histidine
    "ILE",  # Isoleucine
    "LEU",  # Leucine
    "LYS",  # Lysine
    "MET",  # Methionine
    "PHE",  # Phenylalanine
    "PRO",  # Proline
    "SER",  # Serine
    "THR",  # Threonine
    "TRP",  # Tryptophan
    "TYR",  # Tyrosine
    "VAL",  # Valine

    # Phosphorylation
    "SEP", # Phosphoserine
    "TPO", # Phosphothreonine
    "PTR", # Phosphotyrosine

    # Methylation
    "MLY",  # Methylated Lysine → Lysine
    "MLZ",  # Dimethylated Lysine → Lysine
    "M3L",  # Trimethylated Lysine → Lysine
    "MRA",  # Methylated Arginine → Arginine

    # Acetylation
    "ALY",  # N-acetylated Lysine → Lysine
    "ACE",  # N-acetylated Serine → Serine

    # Glycosylation
    "NAG",  # N-acetylglucosaminyl-Asparagine → Asparagine
    "GAL",  # O-linked N-acetylgalactosamine (Ser/Thr) → Serine

    # Ubiquitination & SUMOylation
    "UBI",  # Ubiquitinated Lysine → Lysine
    "SUM",  # SUMOylated Lysine → Lysine

    # Hydroxylation
    "HYP",  # Hydroxyproline → Proline
    "HYL",  # Hydroxylysine → Lysine

    # Other PTMs
    "SEC",  # Selenocysteine → Cysteine
    "PYL",  # Pyrrolysine → Lysine
    "CGU",  # Carboxylated Glutamate → Glutamate
    "FME",
    "MSE"
]

aa_map = {
    # Canonical 20 AA
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",

    # Phosphorylation
    "SEP": "S",   # phosphoserine
    "TPO": "T",   # phosphothreonine
    "PTR": "Y",   # phosphotyrosine

    # Methylation (→ Lys or Arg)
    "MLY": "K",
    "MLZ": "K",
    "M3L": "K",
    "MRA": "R",

    # Acetylation
    "ALY": "K",
    "ACE": "S",

    # Glycosylation
    "NAG": "N",
    "GAL": "S",   # O-GalNAc on Ser/Thr → map to Ser (your choice)

    # Ubiquitination / SUMOylation
    "UBI": "K",
    "SUM": "K",

    # Hydroxylation
    "HYP": "P",
    "HYL": "K",

    # Other PTMs / uncommon residues
    "SEC": "C",   # Selenocysteine
    "PYL": "K",   # Pyrrolysine → Lysine functional group
    "CGU": "E",   # γ-carboxy-glutamate → Glu
    "FME": "M",   # N-formylmethionine → Met
    "MSE": "M"
}

def get_index(chain, res_index):
    counter = 0
    for residue in chain:
        for atom in residue:
            if atom.get_name() == "CA":
                counter += 1
                if residue.get_id()[1] == res_index:
                    return counter
                
def get_sequence(chain):
    seq = ""
    for residue in chain:
        for atom in residue:
            if atom.get_name() == "CA" and residue.get_resname() in amino_acids:
                seq += aa_map[residue.get_resname()]
                break
    return seq

def non_bio_check(chain):
    for residue in chain:
        CA_check = False
        C_check = False
        O_check = False
        for atom in residue:
            if atom.get_name() == "CA":
                CA_check = True
            if atom.get_name() == "C":
                C_check = True
            if atom.get_name() == "O":
                O_check = True
        if CA_check == True and C_check == True and O_check == True and residue.get_resname() not in amino_acids:
            return False    
    return True

def non_bio_check_2(structure, chain_id):
    for model in structure:
        for chain in model:
            if chain.get_id() == chain_id:
                for residue in chain:
                    N_check = False
                    CA_check = False
                    C_check = False
                    O_check = False
                    for atom in residue:
                        if atom.get_name() == "N":
                            N_check = True
                        if atom.get_name() == "CA":
                            CA_check = True
                        if atom.get_name() == "C":
                            C_check = True
                        if atom.get_name() == "O":
                            O_check = True
                    if N_check == True and CA_check == True and C_check == True and O_check == True and residue.get_resname() not in amino_acids:
                        return False   
        break
    return True

def analyze_struct(structure, ns):
    total_num_phos_holder = 0
    total_num_inter_holder = 0
    contacts_holder = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() in ['SEP', 'TPO']:
                    for atom in residue:
                        if atom.get_name() == 'P':
                            bio_check = non_bio_check(chain)
                            if bio_check == True:
                                phos_seq = ""
                                bind_seq = ""
                                phos_seq = get_sequence(chain)
                                total_num_phos_holder += 1
                                neighbors = ns.search(atom.get_coord(), 6.0)
                                taken = []
                                for atom2 in neighbors:
                                    if f"{atom2.get_parent().get_parent().get_id()} - {atom2.get_parent().get_id()[1]}" not in taken and atom2.get_parent().get_parent().get_id() != chain.get_id() and atom2.get_parent().get_resname() in amino_acids:
                                        taken.append(f"{atom2.get_parent().get_parent().get_id()} - {atom2.get_parent().get_id()[1]}")
                                        contacts_holder.append(atom2.get_parent().get_resname())
                                chain_neighbor_counter = []
                                inter_distance = 100
                                inter_index = None
                                for atom2 in neighbors:
                                    if atom2.get_parent().get_parent().get_id() != chain.get_id() and atom2.get_parent().get_resname() in amino_acids and atom2.get_parent().get_parent().get_id() not in chain_neighbor_counter:
                                        chain_neighbor_counter.append(atom2.get_parent().get_parent().get_id())
                                for atom2 in neighbors:
                                    if atom2.get_parent().get_parent().get_id() != chain.get_id() and atom2.get_parent().get_resname() in amino_acids and non_bio_check_2(structure, chain_neighbor_counter[0]) == True:
                                        if atom2 - atom < inter_distance:
                                            inter_distance = atom2 - atom
                                            inter_index = get_index(atom2.get_parent().get_parent(), atom2.get_parent().get_id()[1])
                                            bind_seq = get_sequence(atom2.get_parent().get_parent())
                                if len(chain_neighbor_counter) == 1 and bind_seq != "":
                                    with open('all_SEP_TPO_in_PDB_nov_26_2025_1.csv', 'a', newline='') as outfile:
                                        writer = csv.writer(outfile)
                                        writer.writerow([structure.get_id(), chain.get_id(), get_index(chain, residue.get_id()[1]), chain_neighbor_counter[0], inter_index, phos_seq, bind_seq])
                                    total_num_inter_holder += 1
        break
    return [total_num_phos_holder, total_num_inter_holder]

if __name__ == "__main__":
    counter = 0
    total_num_inter = 0
    total_num_phos = 0
    parser = MMCIFParser(QUIET = True)
    print(f'\nTotal number of files to analyze: {len(os.listdir("pdb_structures_sep_tpo_nov_26_2025"))}\n')
    file_counter = 0
    for file in sorted(os.listdir('pdb_structures_sep_tpo_nov_26_2025')):
        print(file_counter)
        file_counter += 1
        structure = parser.get_structure(file.split('.')[0], f"pdb_structures_sep_tpo_nov_26_2025/{file}")
        ns = NeighborSearch(list(structure.get_atoms()))
        results = analyze_struct(structure, ns)
        total_num_phos += results[0]
        total_num_inter += results[1]

    print(f"Total number of p-THR and p-SER: {total_num_phos}")
    print(f"Total number of p-THR and p-SER interacting with one other protein: {total_num_inter}")
