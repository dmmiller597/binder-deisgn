# Protein Design Competition: EGFR Binding Protein

## Overview

This competition challenges protein designers to create a novel binding protein targeting the extracellular domain of EGFR (Epidermal Growth Factor Receptor), a significant cancer-associated drug target. The most promising designs will undergo experimental validation in our state-of-the-art automated wet lab, with all results being open-sourced to benefit the scientific community.

## Target Specifications

- **Protein**: EGFR (Epidermal Growth Factor Receptor)
- **Organism**: Human
- **UniProt Accession ID**: P00533

## Submission Guidelines

Each participant may submit up to 10 unique designs adhering to the following criteria:

1. Utilize only natural amino acids
2. Maximum length of 200 amino acids
3. Single-chain design
4. Maintain at least 10 amino acid edit distance from known EGFR binders

## Evaluation Process

1. Initial ranking based on PAE_interaction scores using AlphaFold2 (AF2)
2. Top 100 designs selected for experimental testing.
3. Additional designs may be chosen based on novelty, creative design approach, or other noteworthy factors.
4. The top 5 designers achieving the best binding affinity (lowest KD value) will receive awards.


## Approach from leading candidate

1. Set up the computational environment
   - Install and configure ColabFold (https://github.com/sokrypton/ColabFold)
     - Use Google Colab or local installation with GPU support
   - Set up a Modal app with minimalaf (https://github.com/whitead/minimalaf)
     - Create a Modal account and install the Modal CLI
     - Set up a new Modal app using the minimalaf template
     - Configure the app with necessary compute resources (CPU, GPU, memory)
   - Prepare the MSA for EGFR
     - Use UniRef90 database for diverse sequences
     - Apply MMseqs2 for sequence search and clustering
     - Align sequences using MAFFT or Clustal Omega
     - Save MSA in A3M format for use with AlphaFold

2. Implement the initial experiments
   - Write scripts to analyze EGF subsequences
     - Extract EGF sequence from UniProt (P01133)
     - Create sliding window function (window sizes: 20, 30, 40 amino acids)
     - Evaluate each subsequence with AlphaFold Multimer
   - Create an alanine scan function
     - Implement systematic alanine substitution for each residue
     - Use AlphaFold to predict structure and calculate iPAE for each mutant
   - Develop a random mutation generator
     - Use numpy for random number generation
     - Implement functions for point mutations, insertions, and deletions
     - Ensure mutations maintain sequence length within 200 amino acids

3. Develop the proposal engine
   - Integrate PepMLM (https://github.com/programmablebio/pepmlm)
     - Install PepMLM and its dependencies
     - Fine-tune PepMLM on EGFR-binding sequences if available
   - Create functions to interface with the AlphaFold Multimer oracle
     - Implement sequence input and iPAE score output functions
     - Use Modal for distributed computation of AlphaFold predictions

4. Build the dataset generation pipeline
   - Implement functions to feed sequences to the oracle and collect iPAE scores
     - Create batching system for efficient processing
     - Implement error handling and retries for failed predictions
   - Create a database or file system to store the sequence-iPAE pairs
     - Use SQLite for local storage or PostgreSQL for distributed setup
     - Implement functions for efficient data insertion and retrieval

5. Develop the surrogate model
   - Implement the 1D convolutional neural network ensemble
     - Use PyTorch or TensorFlow for model implementation
     - Create 3-5 models with varying architectures (e.g., different filter sizes, layer depths)
     - Implement ensemble averaging for predictions
   - Create training and evaluation scripts for the surrogate model
     - Use k-fold cross-validation for robust evaluation
     - Implement early stopping and learning rate scheduling
     - Use MSE or custom loss function based on iPAE score distribution

6. Integrate EvoProtGrad
   - Set up EvoProtGrad for guided sampling
     - Install EvoProtGrad and its dependencies
     - Configure EvoProtGrad parameters (population size, mutation rate, etc.)
   - Create an interface between the surrogate model and EvoProtGrad
     - Implement custom fitness function using surrogate model predictions
     - Ensure proper sequence encoding/decoding between systems

7. Implement the iterative improvement loop
   - Develop scripts to continuously generate, evaluate, and refine sequences
     - Use EvoProtGrad for sequence generation
     - Evaluate top candidates with AlphaFold Multimer oracle
     - Update dataset with new sequence-iPAE pairs
   - Implement the retraining process for the surrogate model
     - Retrain model every N iterations or when performance degrades
     - Implement checkpointing for model weights and training state

8. Set up scaling with Modal Labs
   - Configure the Modal app for distributed computing
     - Set up worker pools for different tasks (AlphaFold prediction, surrogate model training)
     - Implement job queuing system for efficient resource utilization
   - Adapt the iterative improvement loop for scaled execution
     - Parallelize sequence generation and evaluation steps
     - Implement distributed data storage and synchronization

9. Create the final selection process
   - Implement the iPAE cutoff filter
     - Use iPAE threshold of 8.5 as mentioned in the approach
   - Develop a greedy optimization algorithm for maximizing Hamming distance
     - Implement Hamming distance calculation function
     - Create iterative selection process to maximize diversity
   - Create a script to sample and prepare the final submission set
     - Ensure all submission guidelines are met (length, natural amino acids, etc.)
     - Generate necessary metadata for each sequence