from Bio.PDB import PDBParser

from Bio import PDB
from Bio.PDB import Superimposer


def linker(pdb_path_1, pdb_path_2, design_name):
    # load the structures
    CD20_dimer_structure = PDBParser().get_structure('6Y97', '../../data/cd20/6Y97_reindexed.pdb')
    # this one has the two monomers and its epitope
    # just get chain A from this structure

    filepath_1 = pdb_path_1
    filepath_2 = pdb_path_2
    simple_design_1 = PDBParser().get_structure('simple', filepath_1)
    simple_design_2 = PDBParser().get_structure('simple', filepath_2)

    # align first of the designs to CD20 dimer chain A
    cd20_atoms_chain_A = []
    simple_design_atoms = []

    for model in CD20_dimer_structure:
        for chain in model:
            if chain.id == 'A':
                residues = [i for i in chain.get_residues()]
                for residue in residues:
                    if residue.has_id("CA"):
                        cd20_atoms_chain_A.append(residue["CA"])

    for model in simple_design_1:
        for chain in model:
            if chain.id == 'A':
                residues = [i for i in chain.get_residues()]
                for residue in residues:
                    if residue.has_id("CA"):
                        simple_design_atoms.append(residue["CA"])


    super_imposer = Superimposer()
    super_imposer.set_atoms(cd20_atoms_chain_A, simple_design_atoms)
    super_imposer.apply(simple_design_1.get_atoms()) # epitope 1 is aligned to CD20

    # now do the same with simple_design_2
    cd20_atoms_chain_B = []
    simple_design_atoms = []

    for model in CD20_dimer_structure:
        for chain in model:
            if chain.id == 'B':
                residues = [i for i in chain.get_residues()]
                for residue in residues:
                    if residue.has_id("CA"):
                        cd20_atoms_chain_B.append(residue["CA"])

    for model in simple_design_2:
        for chain in model:
            if chain.id == 'A':
                residues = [i for i in chain.get_residues()]
                for residue in residues:
                    if residue.has_id("CA"):
                        simple_design_atoms.append(residue["CA"])

    super_imposer = Superimposer()
    super_imposer.set_atoms(cd20_atoms_chain_B, simple_design_atoms)
    super_imposer.apply(simple_design_2.get_atoms()) # epitope 2 is aligned to CD20

    # do it differently, rather combine the two simple_designs into one structure
    # and then add the CD20 dimer
    new_structure = PDB.Structure.Structure('binder_cd20_dimer')
    new_model = PDB.Model.Model(0)
    new_structure.add(new_model)

    # add the binders
    cd20_monomer_1, binder_1 = simple_design_1.get_chains()
    cd20_monomer_2, binder_2 = simple_design_2.get_chains()

    cd20_monomer_1.id = 'A'
    cd20_monomer_2.id = 'B'
    new_model.add(cd20_monomer_1)
    new_model.add(cd20_monomer_2)
    
    # combine the binders into a single chain and make the indices be separated by a 40 gap
    from Bio.PDB.Chain import Chain
    new_chain = Chain('C')
    # Add residues from binder_1
    for residue in binder_1:
        new_chain.add(residue.copy())

    # add 40 random residues
    from Bio.PDB.Residue import Residue
    for i in range(21, 21 + 40):
        residue.id = (' ', i, ' ')
        new_chain.add(residue.copy())

    # Add residues from binder_2 with a gap of 40
    last_resid = max(int(res.id[1]) for res in new_chain)
    for residue in binder_2:
        new_residue = residue.copy()
        new_resid = int(residue.id[1]) + last_resid
        new_residue.id = (' ', new_resid, ' ')
        new_chain.add(new_residue)

    # Renumber atoms to ensure unique atom serial numbers
    atom_count = 1
    for residue in new_chain:
        for atom in residue:
            atom.serial_number = atom_count
            atom_count += 1
    
    # Add the new combined chain to the model
    new_model.add(new_chain)

    io = PDB.PDBIO()
    io.set_structure(new_structure)
    io.save(f"../../data/designs/simple/{design_name}.pdb")

if __name__ == "__main__":
    linker("../../data/designs/simple/fixedA.pdb", "../../data/designs/simple/fixedA.pdb", "fixedAA")
    linker("../../data/designs/simple/fixedA.pdb", "../../data/designs/simple/fixedB.pdb", "fixedAB")
    linker("../../data/designs/simple/fixedB.pdb", "../../data/designs/simple/fixedB.pdb", "fixedBB")