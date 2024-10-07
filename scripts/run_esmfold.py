import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.esmfold import EsmFoldv1

if __name__ == "__main__":
    device = 'cuda'

    esmfold = EsmFoldv1()
    esmfold.load(
        device=device,
        memory_light=True,
    )

    filepath = "../outputs/dimer5/SolubleMPNN/dimer5-solMPNN.fasta"

    # iterate through and get all sequences
    from Bio.SeqIO.FastaIO import SimpleFastaParser

    sequences = []
    with open(filepath, "r") as file:
        for title, sequence in SimpleFastaParser(file):
            sequences.append(sequence)


    for i, sequence in enumerate(sequences):
        print(f"Running sequence {i} / {len(sequences)}: {sequence}")

        # split by :
        sequences = [sequence.split(":")[0] for sequence in sequences]

        CD20_A = sequences[0]
        CD20_B = sequences[1]
        binder = sequences[2]

        # CD20_A = "FMRESKTLGAVQIMNGLFHIALGGLLMIPAGIYAPICVTVWYPLWGGIMYIISGSLLAATEKNSRKCLVKGKMIMNSLSLFAAISGMILSIMDILNIKISHFLKMESLNFIRAHTPYINIYNCEPANPSEKNSPSTQYCYSIQSLFLGILSVMLIFAFFQELVIAGIVE"

        # CD20_B = "FMRESKTLGAVQIMNGLFHIALGGLLMIPAGIYAPICVTVWYPLWGGIMYIISGSLLAATEKNSRKCLVKGKMIMNSLSLFAAISGMILSIMDILNIKISHFLKMESLNFIRAHTPYINIYNCEPANPSEKNSPSTQYCYSIQSLFLGILSVMLIFAFFQELVIAGIVE"
        

        # # with open("data/test.fasta", "r") as file:
        # binder = 'FLIFLHHSNSYLYNDAYHNRMMQKPDDEREGDPEYALYAAQRNAHMGMMYWMYEWHNWLGFKATPNEYPGLTEVVPELTR'
        
        # # strip the last part of the sequence separated by :
        # binder = binder.split(":")[0]

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


        print('Running sequence')

        result = esmfold.fold(
            sequence=sequence,
            residue_indices=residue_indices,
            print_pdb=True,
            pdb_file=f"outputs/dimer5/EsmFold/fold_{i}.pdb",
            device=device,
        )
        # save also result json, delete result.atoms
        result_json = result.model_dump()
        del result_json["atoms"]
        with open(f"outputs/dimer5/EsmFold/fold_{i}.json", "w") as file:
            json.dump(result_json, file)

        print("Done")




