import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.esmfold import EsmFoldv1

from Bio.SeqIO.FastaIO import SimpleFastaParser
import torch

if __name__ == "__main__":
    device = 'cuda'

    esmfold = EsmFoldv1()
    esmfold.load(
        device=device,
        memory_light=True,
    )

    design = 'simpleAA'

    for design in ['simpleAA', 'simpleBB', 'simpleAB']:

        for j in range(0, 10):

            filepath = f"outputs/{design}/SolubleMPNN/seqs/_{j}.fa"


            sequences = []
            with open(filepath, "r") as file:
                for title, sequence in SimpleFastaParser(file):
                    sequences.append(sequence)


            for i, sequence in enumerate(sequences):
                print(f"Running sequence {i} / {len(sequences)}: {sequence}")

                CD20_A = "FMRESKTLGAVQIMNGLFHIALGGLLMIPAGIYAPICVTVWYPLWGGIMYIISGSLLAATEKNSRKCLVKGKMIMNSLSLFAAISGMILSIMDILNIKISHFLKMESLNFIRAHTPYINIYNCEPANPSEKNSPSTQYCYSIQSLFLGILSVMLIFAFFQELVIAGIVE"

                CD20_B = "FMRESKTLGAVQIMNGLFHIALGGLLMIPAGIYAPICVTVWYPLWGGIMYIISGSLLAATEKNSRKCLVKGKMIMNSLSLFAAISGMILSIMDILNIKISHFLKMESLNFIRAHTPYINIYNCEPANPSEKNSPSTQYCYSIQSLFLGILSVMLIFAFFQELVIAGIVE"
                
                # remove CD20_A from sequence
                sequence = sequence.replace(CD20_A, "")
                # remove :
                sequence = sequence.replace(":", "")

                binder = sequence


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

                if not os.path.exists(f"outputs/{design}/EsmFold/"):
                    os.makedirs(f"outputs/{design}/EsmFold/")

                result = esmfold.fold(
                    sequence=sequence,
                    residue_indices=residue_indices,
                    print_pdb=True,
                    pdb_file=f"outputs/{design}/EsmFold/fold_{j}_{i}.pdb",
                    device=device,
                )
                # save also result json, delete result.atoms
                # turn a dataclass into a dictionary
                from dataclasses import asdict
                import json
                import numpy as np
                result_dict = asdict(result)
                del result_dict["atoms"]
                # turn all attributes into strings (if it is a numpy array turn it into lists of lists)
                for key, value in result_dict.items():
                    if isinstance(value, (np.ndarray, torch.Tensor)):
                        result_dict[key] = value.tolist()
                    elif isinstance(value, dict):
                        for k, v in value.items():
                            if isinstance(v, (np.ndarray, torch.Tensor)):
                                value[k] = v.tolist()
                with open(f"outputs/{design}/EsmFold/fold_{j}_{i}.json", "w") as file:
                    json.dump(result_dict, file)

                print(f"Done {design} {j}")


