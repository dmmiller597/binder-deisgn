import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.esmfold import EsmFoldv1

# from Bio.PDB import PDBParser, PDBIO
# from Bio.SeqUtils import seq1

# # Parse the PDB file
# parser = PDBParser()
# structure = parser.get_structure('6Y97', '6Y97.pdb')

# # Assert that there are two chains
# assert len(list(structure.get_chains())) == 2, "There should be exactly two chains in the structure."

# for chain in structure.get_chains():
#     print(f">Chain {chain.id}")
#     sequence = ""
#     for residue in chain.get_residues():
#         resname = residue.get_resname()
#         sequence += seq1(resname)
#     print(sequence)

if __name__ == "__main__":
    esmfold = EsmFoldv1()
    esmfold.load("cuda:0")

    CD20_A = "FMRESKTLGAVQIMNGLFHIALGGLLMIPAGIYAPICVTVWYPLWGGIMYIISGSLLAATEKNSRKCLVKGKMIMNSLSLFAAISGMILSIMDILNIKISHFLKMESLNFIRAHTPYINIYNCEPANPSEKNSPSTQYCYSIQSLFLGILSVMLIFAFFQELVIAGIVE"

    CD20_B = "FMRESKTLGAVQIMNGLFHIALGGLLMIPAGIYAPICVTVWYPLWGGIMYIISGSLLAATEKNSRKCLVKGKMIMNSLSLFAAISGMILSIMDILNIKISHFLKMESLNFIRAHTPYINIYNCEPANPSEKNSPSTQYCYSIQSLFLGILSVMLIFAFFQELVIAGIVE"
    

    with open("data/input_sequence.fasta", "r") as file:
        binder = file.readlines()[1].strip()
    
        # strip the last part of the sequence separated by :
        binder = binder.split(":")[0]

    multimer_gap = 256

    linker_sequence = "G" * 100 # need to make at least 100 so that it can circle all the way back and around

    # they also suggest changing the order of the sequences
    sequence = CD20_A + linker_sequence + CD20_B + linker_sequence + binder

    multimer_gap = 256

    residue_indices = list(range(1, len(sequence) + 1))


    residue_indices[len(CD20_A) + len(linker_sequence):] = [
        i + multimer_gap for i in residue_indices[len(CD20_A) + len(linker_sequence):]
    ]

    residue_indices[len(CD20_A) + len(linker_sequence) + len(CD20_B) + len(linker_sequence):] = [
        i + multimer_gap for i in residue_indices[len(CD20_A) + len(linker_sequence) + len(CD20_B) + len(linker_sequence):]
    ]


    esmfold.fold(
        sequence=sequence,
        residue_indices=residue_indices,
        print_pdb=True,
        pdb_file="outputs/esmfold_output.pdb",
    )

    print("Done")




