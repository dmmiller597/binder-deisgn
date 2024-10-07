from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Chain import Chain


def fixer(pdb_path):
    print('Splitting chains...')
    # SEE https://stackoverflow.com/questions/74735845/splitting-and-renaming-protein-chain-with-biopythons-biopdb
    structure = PDBParser().get_structure('dimer1', pdb_path)
    # get residue id where to break
    last_id = 0
    res_to_chain_B = []
    res_to_chain_C = []
    chainB = False
    chainC = False
    for model in structure:
        for chain in model:
            for residue in chain:
                if chainC:
                    res_to_chain_C.append(residue)
                    continue
                elif chainB:
                    res_to_chain_B.append(residue)
                    if abs(residue.id[1] - last_id) > 2:
                        chainC = True
                        res_to_chain_C.append(residue)
                    continue
                elif abs(residue.id[1] - last_id) > 2:
                    chainB = True
                    res_to_chain_B.append(residue)
                last_id = residue.id[1]

    ### SEE https://stackoverflow.com/questions/25884758/deleteing-residue-from-pdb-using-biopython-library
    for model in structure:
        for chain in model:
            for res in res_to_chain_B:
                if res.get_id() in chain.child_dict:
                    chain.detach_child(res.get_id())
            for res in res_to_chain_C:
                if res.get_id() in chain.child_dict:
                    chain.detach_child(res.get_id())

    ### SEE https://stackoverflow.com/questions/33364370/how-to-add-chain-id-in-pdb
    my_chain_B = Chain("B")
    my_chain_C = Chain("C")
    model.add(my_chain_B)
    model.add(my_chain_C)
    for index, res in enumerate(res_to_chain_B):
        # reindex
        res.id = (' ', index+1, ' ')
        my_chain_B.add(res)
    for index, res in enumerate(res_to_chain_C):
        # reindex
        res.id = (' ', index+1, ' ')
        my_chain_C.add(res)

    # Create a new structure to avoid PDBConstructionException
    from Bio.PDB import Structure, Model
    new_structure = Structure.Structure('new_dimer1')
    new_model = Model.Model(0)
    new_structure.add(new_model)
    new_model.add(my_chain_B)
    new_model.add(my_chain_C)

    io = PDBIO()
    io.set_structure(new_structure)
    io.save('../data/designs/dimer1/original.pdb')


if __name__ == "__main__":
    fixer("/home/jakub/phd/binder-design/data/stefano/dimer/CARTCD20__best_123.pdb")