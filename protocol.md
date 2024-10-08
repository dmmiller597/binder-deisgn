# Protocol for Designing an 80 Amino Acid De Novo Binder for CAR T-cell Construct

## Overview

This protocol provides a detailed, step-by-step guide to designing an 80 amino acid de novo binder targeting the CD20 antigen for incorporation into the CAR T-cell construct specified in `data/car-construct.txt`. The objective is to develop a binder with high specificity and affinity, ensuring compatibility and optimal functionality within the CAR framework.

## Step-by-Step Plan

### 1. Define the Target Antigen

- **Identify the Target Antigen**
  - **Target:** CD20 antigen.
  - **Function:** CD20 is a B-cell surface molecule involved in B-cell proliferation and differentiation.

- **Obtain High-Resolution 3D Structure**
  - **Tool:** Retrieve the crystal structure from the Protein Data Bank (PDB).
    - **Example PDB ID:** `P04206` (Replace with actual CD20 structure ID).
  - **Procedure:**
    1. Visit [RCSB PDB](https://www.rcsb.org/).
    2. Search for "CD20" and download the relevant PDB file.
    3. Verify the resolution and quality of the structure.

### 2. Analyze Antigen Structure

- **Inspect Antigen Structure**
  - **Tool:** PyMOL or ChimeraX.
  - **Procedure:**
    1. Load the PDB file into PyMOL/ChimeraX.
    2. Visualize secondary and tertiary structures.
  
- **Identify Potential Epitopes and Binding Sites**
  - **Tool:** Epitope mapping software such as **IEDB Analysis Resource**.
  - **Procedure:**
    1. Input the CD20 structure into IEDB.
    2. Analyze surface accessibility and antigenicity to determine potential epitopes.
  
- **Map Surface Properties**
  - **Tool:** APBS (Adaptive Poisson-Boltzmann Solver) for electrostatic surface mapping.
  - **Procedure:**
    1. Calculate electrostatic potential.
    2. Identify positively and negatively charged regions conducive to binding.

### 3. Generate Initial Binder Scaffolds

- **Use De Novo Protein Design Software**
  - **Tool:** RosettaRemodel from the Rosetta Suite.
  - **Procedure:**
    1. Set up Rosetta environment.
    2. Define backbone conformations targeting identified epitopes.
    3. Generate a diverse set of backbone scaffolds (e.g., 1000 variants).

### 4. Dock Scaffolds to Antigen

- **Perform Protein-Protein Docking Simulations**
  - **Tool:** RosettaDock.
  - **Procedure:**
    1. Align binder scaffolds near target epitopes.
    2. Execute docking simulations for each scaffold.
  
- **Score and Rank Docking Poses**
  - **Tool:** Rosetta scoring functions.
  - **Procedure:**
    1. Evaluate binding affinity scores.
    2. Rank poses based on lowest energy and highest binding probability.

### 5. Design Binder Interface

- **Optimize Interface Residues**
  - **Tool:** RosettaDesign.
  - **Procedure:**
    1. Identify key interacting residues at the interface.
    2. Perform side-chain optimization to enhance binding contacts.

- **Generate Multiple Sequence Variants**
  - **Tool:** RosettaScripts with sequence optimization protocols.
  - **Procedure:**
    1. Introduce mutations to interface residues.
    2. Generate variants with improved binding characteristics.

### 6. Sequence Optimization

- **Enhance Sequences Using Machine Learning Models**
  - **Tool:** ProteinMPNN.
  - **Procedure:**
    1. Input optimized binder scaffolds into ProteinMPNN.
    2. Generate sequence variants with enhanced stability and affinity.

### 7. Filter and Rank Candidates

- **Evaluate Candidates**
  - **Criteria:** Binding affinity, structural stability, and specificity.
  - **Tool:** Rosetta binding energy calculations and FoldX for stability.
  - **Procedure:**
    1. Calculate binding free energies.
    2. Assess thermodynamic stability of each variant.
    3. Filter out candidates with suboptimal scores.

- **Rank Candidates**
  - **Tool:** Custom ranking algorithm or Excel.
  - **Procedure:**
    1. Assign scores based on predefined criteria.
    2. Select top 50 candidates for further analysis.

### 8. Predict Binder Structures

- **Use AlphaFold2 for 3D Structure Prediction**
  - **Tool:** AlphaFold2.
  - **Procedure:**
    1. Input amino acid sequences of top candidates into AlphaFold2.
    2. Generate predicted 3D structures.

- **Compare with Designed Models**
  - **Tool:** PyMOL or ChimeraX.
  - **Procedure:**
    1. Align predicted structures with initial designs.
    2. Identify discrepancies and refine models accordingly.

### 9. Validate Binding In Silico

- **Redock Predicted Structures**
  - **Tool:** RosettaDock.
  - **Procedure:**
    1. Dock AlphaFold2 predicted structures back to CD20.
    2. Assess consistency with initial docking results.

- **Perform Molecular Dynamics Simulations**
  - **Tool:** GROMACS or AMBER.
  - **Procedure:**
    1. Set up simulation parameters for binder-antigen complexes.
    2. Run simulations to assess stability and dynamics.

- **Calculate Binding Free Energies**
  - **Tool:** MM-PBSA or MM-GBSA within GROMACS/AMBER.
  - **Procedure:**
    1. Extract energy components from simulations.
    2. Compute binding free energies to quantify affinity.

### 10. Assess Immunogenicity

- **Analyze for Potential Immunogenic Epitopes**
  - **Tool:** NetMHCpan.
  - **Procedure:**
    1. Input binder sequences into NetMHCpan.
    2. Predict binding to MHC class I and II molecules.

- **Modify Sequences to Reduce Immunogenicity Risks**
  - **Tool:** Epitope prediction feedback.
  - **Procedure:**
    1. Identify high-risk epitopes.
    2. Introduce mutations to eliminate or reduce epitope presentation.

### 11. Finalize Binder Sequences

- **Compile Top Candidate Sequences**
  - **Tool:** Excel or database management system.
  - **Procedure:**
    1. Select sequences with best binding, stability, and low immunogenicity.
  
- **Validate Compatibility with CAR Construct**
  - **Tool:** In silico cloning tools (e.g., SnapGene).
  - **Procedure:**
    1. Insert binder sequences into `data/car-construct.txt`.
    2. Ensure proper reading frame and linker compatibility.

### 12. Insert Binders into CAR Construct

- **Replace `<de_novo_binder>` Placeholder**
  - **Tool:** Text editor or script (e.g., Python).
  - **Procedure:**
    1. Open `data/car-construct.txt`.
    2. Replace `<de_novo_binder>` with selected binder sequence.

- **Ensure Proper Linkage with Flanking Regions**
  - **Tool:** Sequence alignment tools.
  - **Procedure:**
    1. Verify N-terminal and C-terminal linkers are correctly connected.
    2. Confirm absence of frame shifts or premature stop codons.

### 13. Documentation and Data Management

- **Record Design Parameters and Methodologies**
  - **Tool:** Laboratory Information Management System (LIMS) or documentation software.
  - **Procedure:**
    1. Document each step, parameters used, and decisions made.
  
- **Maintain a Database of All Candidate Binders**
  - **Tool:** SQL database or cloud-based storage (e.g., GitHub).
  - **Procedure:**
    1. Store sequences, structures, and evaluation scores.
    2. Ensure data is backed up and version-controlled.

## Key Software and Tools

- **Structural Analysis & Visualization**
  - PyMOL, ChimeraX
- **De Novo Protein Design**
  - Rosetta Suite (RosettaRemodel, RosettaDock, RosettaDesign)
- **Machine Learning for Sequence Optimization**
  - ProteinMPNN
- **3D Structure Prediction**
  - AlphaFold2
- **Molecular Dynamics Simulations**
  - GROMACS, AMBER
- **Immunogenicity Prediction**
  - NetMHCpan
- **Sequence Alignment & Cloning**
  - SnapGene
- **Energy Calculations**
  - FoldX, MM-PBSA, MM-GBSA
- **Epitope Mapping**
  - IEDB Analysis Resource
- **Electrostatic Surface Mapping**
  - APBS
- **Data Management**
  - SQL Databases, GitHub

## Notes

- **Diversity of Candidates:** Generate a wide array of binder variants to enhance the probability of identifying high-affinity, specific binders.
- **Complementary Validation:** Pair in silico findings with experimental validations to confirm binding efficacy and specificity.
- **Iterative Refinement:** Utilize feedback from simulations and predictions to iteratively improve binder designs.

## Updating the CAR Construct

- **Replace `<de_novo_binder>` Placeholder:**
  - Insert the finalized binder sequence into the `data/car-construct.txt` at the `<de_novo_binder>` position.
  
- **Verify Correct Alignment:**
  - Ensure that the inserted binder sequence aligns seamlessly with upstream and downstream regions.
  - Confirm that linkers facilitate proper expression and functionality of the CAR construct.
