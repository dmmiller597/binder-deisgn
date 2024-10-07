from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Chain import Chain


def fixer(pdb_path, design_name):
    print('Splitting chains...')        
    # SEE https://stackoverflow.com/questions/74735845/splitting-and-renaming-protein-chain-with-biopythons-biopdb
    structure = PDBParser().get_structure(design_name, pdb_path)
    # get residue ids where to break
    all_ids = []
    for model in structure:
        for chain in model:
            for residue in chain:
                all_ids.append(residue.id[1])
    # find the indices of the residues where to break
    indices_to_break = []
    for i in range(len(all_ids) - 1):
        if abs(all_ids[i] - all_ids[i + 1]) > 2:
            indices_to_break.append(i)

    chains = {'A': [], 'B': [], 'C': []}

    from Bio.PDB import Structure, Model
    new_structure = Structure.Structure(design_name)
    new_model = Model.Model(0)
    new_structure.add(new_model)

    for model in structure:
        for chain in model:
            if 'simple' not in design_name:
                # get residues 0 to first break
                residues = [i for i in chain.get_residues()]

                # first chain
                chains['A'] = residues[0:indices_to_break[0]+1]
                # second chain
                chains['B'] = residues[indices_to_break[0]+1:indices_to_break[1]+1]
                # third chain
                chains['C'] = residues[indices_to_break[1]+1:]
            else:
                # get residues 0 to first break
                residues = [i for i in chain.get_residues()]
                # first chain
                chains['A'] = residues[0:indices_to_break[0]+1]
                # second chain
                chains['B'] = residues[indices_to_break[0]+1:]
                
            
    # add these to a new structure
    for chain_id, residues in chains.items():
        i = 1
        new_chain = Chain(chain_id)
        for residue in residues:
            residue.id = (residue.id[0], i, residue.id[2])
            new_chain.add(residue)
            i += 1
        new_model.add(new_chain)

    io = PDBIO()
    io.set_structure(new_structure)
    if 'simple' in design_name:
        if 'originalA' in pdb_path:
            io.save(f'/home/jakub/phd/binder-design/data/designs/{design_name}/fixedA.pdb')
        elif 'originalB' in pdb_path:
            io.save(f'/home/jakub/phd/binder-design/data/designs/{design_name}/fixedB.pdb')
    else:
        io.save(f'/home/jakub/phd/binder-design/data/designs/{design_name}/fixed.pdb')


if __name__ == "__main__":
    fixer("/home/jakub/phd/binder-design/data/designs/dimer1/original.pdb", 'dimer1')
    fixer("/home/jakub/phd/binder-design/data/designs/dimer5/original.pdb", 'dimer5')
    fixer("/home/jakub/phd/binder-design/data/designs/simple/originalA.pdb", 'simple')
    fixer("/home/jakub/phd/binder-design/data/designs/simple/originalB.pdb", 'simple')