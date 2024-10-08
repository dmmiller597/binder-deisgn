# Install



RFDiffusion
```
git submodule add https://github.com/RosettaCommons/RFdiffusion ./modules/RFdiffusion/
```

Install all submodules
```
git submodule update --init --recursive
```


```
conda create -n bits-in-binder
conda activate bits-in-binder
```

RFDiffusion
```
cd modules/RFdiffusion && mkdir models && cd models
curl -O http://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt
curl -O http://files.ipd.uw.edu/pub/RFdiffusion/12fc204edeae5b57713c5ad7dcb97d39/Base_epoch8_ckpt.pt
curl -O http://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt
curl -O http://files.ipd.uw.edu/pub/RFdiffusion/60f09a193fb5e5ccdc4980417708dbab/Complex_Fold_base_ckpt.pt
curl -O http://files.ipd.uw.edu/pub/RFdiffusion/74f51cfb8b440f50d70878e05361d8f0/InpaintSeq_ckpt.pt
curl -O http://files.ipd.uw.edu/pub/RFdiffusion/76d00716416567174cdb7ca96e208296/InpaintSeq_Fold_ckpt.pt
curl -O http://files.ipd.uw.edu/pub/RFdiffusion/5532d2e1f3a4738decd58b19d633b3c3/ActiveSite_ckpt.pt
```

```
mamba install python=3.9
mamba install pytorch pytorch-cuda=11.8 -c pytorch -c nvidia
mamba install -c dglteam/label/th24_cu118 dgl
mamba install pandas
pip install hydra-core pyrsistent

cd env/SE3Transformer
pip install --no-cache-dir -r requirements.txt
python setup.py install
cd ../.. # change into the root directory of the repository
pip install -e . # install the rfdiffusion module from the root of the repository
```
Should be installed.



LigandMPNN
```
bash get_model_params_ligandmpnn.sh "./model_params"
cd modules/LigandMPNN
pip install -r <(sed 's/==.*//' requirements.txt)
```

For ESMFold, I am using an environment built for our other pipeline with Stefano. Can provide details later, but it's a bit messy as the `transformers` package has some weird error.


RossettaFold-All-Atom
```
module load CUDA/11.2
mamba env create -f environment.yaml
conda activate RFAA

```

ESMFold
```
conda create -n esmfold
conda activate esmfold
mamba install -c "nvidia/label/cuda-11.8.0" cuda-toolkit
mamba install pytorch=2.0.0 pytorch-cuda=11.8 -c pytorch -c nvidia
pip install numpy scipy biotite einops modelcif omegaconf tree dm-tree
pip install Bio ml_collections
pip install transformers[torch]
```
Then need to edit the edit line:
`with ContextManagers([] if False else [torch.no_grad()]):`
on line 1956 in `transformers/src/transformers/models/esm/modeling_esmfold.py`

in file `/home/jakub/anaconda3/envs/esm/lib/python3.12/site-packages/transformers/models/esm/modeling_esmfold.py`


# Current Protocol

1. Use RFdiffusion on the epitopes to generate X backbones.
  - for the entire CD20 dimer, it takes about 2 minutes to generate a single backbone
2. Use ProteinMPNN to generate Y sequences to fold in X backbones (in total having X * Y designs)
3. Validate by co-folding with ESMFold ???
  - it takes about 14 seconds to fold the CD20 dimer + the binder
4. Surrogate modeling of protein-protein interactions to validate?? Also, maybe optimise for something else here?
5. Full-scale MD? (entire CAR-T complex? + CD20 + membrane)


For the scripts that I use, look into the `job_` files.

# CD20
FMRESKTLGAVQIMNGLFHIALGGLLMIPAGIYAPICVTVWYPLWGGIMYIISGSLLAATEKNSRKCLVKGKMIMNSLSLFAAISGMILSIMDILNIKISHFLKMESLNFIRAHTPYINIYNCEPANPSEKNSPSTQYCYSIQSLFLGILSVMLIFAFFQELVIAGIVE


# CAR-T construct
MALPVTALLLPLALLLHAARP<>GSTSGSGKPGSGEGSTKGLEPKSCDKTHTCPPCPDPKFWVLVVVGGVLACYSLLVTVAFIIFWVRSKRSRLLHSDYMNMTPRRPGPTRKHYQPYAPPRDFAAYRSRVKFSRSADAPAYKQGQNQLYNELNLGRREEYDVLDKRRGRDPEMGGKPRRKNPQEGLYNELQKDKMAEAYSEIGMKGERRRGKGHDGLYQGLSTATKDTYDALHMQALPPR
 
CD8 leader – MALPVTALLLPLALLLHAARP
 
<> – de novo binder
 
218 linker – GSTSGSGKPGSGEGSTKG
 
CD28 hinge – LEPKSCDKTHTCPPCPDPK
 
CD28 transmembrane – FWVLVVVGGVLACYSLLVTVAFIIFWVR
 
CD28 co-stim – SKRSRLLHSDYMNMTPRRPGPTRKHYQPYAPPRDFAAYRS
 
CD3z stim – RVKFSRSADAPAYKQGQNQLYNELNLGRREEYDVLDKRRGRDPEMGGKPRRKNPQEGLYNELQKDKMAEAYSEIGMKGERRRGKGHDGLYQGLSTATKDTYDALHMQALPPR



# Resources for the BioML Challenge 2024: Bits to Binders

## Table of Contents
- [Install](#install)
- [Current Protocol](#current-protocol)
- [CD20](#cd20)
- [CAR-T construct](#car-t-construct)
- [Resources for the BioML Challenge 2024: Bits to Binders](#resources-for-the-bioml-challenge-2024-bits-to-binders)
  - [Table of Contents](#table-of-contents)
  - [Protein Design](#protein-design)
  - [Protein Structure Prediction](#protein-structure-prediction)
  - [Inverse Folding](#inverse-folding)
  - [Generative Models](#generative-models)
  - [Language Models](#language-models)
  - [Other Useful Tools](#other-useful-tools)
  - [Relevant Literature](#relevant-literature)
  - [Miscellaneous](#miscellaneous)
  - [Biology](#biology)
    - [CAR-T Literature](#car-t-literature)
    - [CD20 Literature](#cd20-literature)
  - [Compute Resources](#compute-resources)
    - [TACC](#tacc)
    - [Other Compute](#other-compute)

## Protein Design
Here is a compilation of many tools used by the UT Austin BioML Society. This is not a comprehensive list so we are likely missing some great tools. For a comprehensive list, check out [this great repository](https://github.com/yangkky/Machine-learning-for-proteins)!

## Protein Structure Prediction
- AlphaFold2 / ColabFold [[paper](https://www.nature.com/articles/s41586-021-03819-2)] [[code](https://github.com/deepmind/alphafold)] / [[colab](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)] [[Docker](https://github.com/kalininalab/alphafold_non_docker)]
- ESMFold [[paper](https://www.nature.com/articles/s41586-022-05543-x)] [[code](https://github.com/facebookresearch/esm)] [[server](https://esmatlas.com/resources?action=fold)]
- RoseTTAFold [[paper](https://www.science.org/doi/10.1126/science.abj8754)] [[code](https://github.com/RosettaCommons/RoseTTAFold)]
- RoseTTAFoldAllAtom [[paper](https://www.biorxiv.org/content/10.1101/2023.05.24.542179v1)] [[code](https://github.com/RosettaCommons/RoseTTAFold)]
- UMol [[paper](https://www.biorxiv.org/content/10.1101/2023.08.01.551553v1)] [[code](https://github.com/patrickbryant1/Umol)]
- AlphaFold2-RAVE [[paper](https://www.biorxiv.org/content/10.1101/2023.09.01.555891v1)] [[code](https://github.com/bjornwallner/alphafold-rave)]

## Inverse Folding
- ProteinMPNN [[paper](https://www.science.org/doi/10.1126/science.add2187)] [[repo](https://github.com/dauparas/ProteinMPNN)]
- LigandMPNN [[paper](https://www.biorxiv.org/content/10.1101/2023.08.01.551615v1)] [[repo](https://github.com/dauparas/LigandMPNN)]
- ESM-IF1 [[paper](https://www.nature.com/articles/s41586-022-05543-x)] [[repo](https://github.com/facebookresearch/esm)] [[colab](https://colab.research.google.com/github/facebookresearch/esm/blob/main/examples/inverse_folding/notebook.ipynb)]

## Generative Models
- RFDiffusion [[paper](https://www.nature.com/articles/s41586-023-06415-8)] [[code](https://github.com/RosettaCommons/RFdiffusion)] [[colab](https://colab.research.google.com/github/RosettaCommons/RFdiffusion/blob/main/notebooks/diffusion.ipynb)] [[Docker](https://github.com/RosettaCommons/RFdiffusion/tree/main/docker)]
- RFDiffusionAllAtom [[paper](https://www.biorxiv.org/content/10.1101/2023.06.22.546080v1)] [[code](https://github.com/RosettaCommons/RFdiffusion)]
- Protpardelle [[paper](https://www.biorxiv.org/content/10.1101/2023.05.24.542189v1)] [[code](https://github.com/RosettaCommons/protpardelle)]
- Chroma [[paper](https://www.biorxiv.org/content/10.1101/2022.12.01.518682v3)] [[code](https://github.com/generatebio/chroma)]
- ProGen [[paper](https://www.nature.com/articles/s41587-022-01618-2)] [[code](https://github.com/salesforce/progen)]
- Framediff [[paper](https://www.biorxiv.org/content/10.1101/2023.05.16.541040v1)] [[code](https://github.com/jasonkyuyim/se3_diffusion)]
- Genie [[paper](https://www.biorxiv.org/content/10.1101/2023.05.29.542705v2)] [[code](https://github.com/LPDI-EPFL/Genie)]
- FoldFlow [[paper](https://www.biorxiv.org/content/10.1101/2023.10.09.561603v1)] [[code](https://github.com/Profluent-Internships/MMDiff)]
- FrameFlow [[paper](https://www.biorxiv.org/content/10.1101/2023.12.22.573103v1)] [[code](https://github.com/Profluent-Internships/FrameFlow)]
- Proteus [[paper](https://www.biorxiv.org/content/10.1101/2023.12.15.571823v1)] [[code](https://github.com/OATML-Markslab/Proteus)]
- Multiflow [[paper](https://www.biorxiv.org/content/10.1101/2024.01.22.576444v1)] [[code](https://github.com/Profluent-Internships/MMDiff)]
- PepFlow [[paper](https://www.biorxiv.org/content/10.1101/2024.01.22.576444v1)] [[code](https://github.com/Profluent-Internships/PepFlow)]

## Language Models
- ESM-2 [[paper](https://www.nature.com/articles/s41586-022-05543-x)] [[code](https://github.com/facebookresearch/esm)]
- SaProt [[paper](https://www.nature.com/articles/s41467-023-38054-y)] [[code](https://github.com/westlake-repl/SaProt)]
- Prost-T5 [[paper](https://www.biorxiv.org/content/10.1101/2023.07.23.550085v1)] [[code](https://github.com/HannesStark/protein-language-models)]
- ProtBERT [[paper](https://ieeexplore.ieee.org/document/9477085)] [[huggingface](https://huggingface.co/Rostlab/prot_bert)]
- ProtTrans [[paper](https://ieeexplore.ieee.org/document/9477085)] [[code](https://github.com/agemagician/ProtTrans)]

## Other Useful Tools
- FoldSeek [[paper](https://www.nature.com/articles/s41587-023-01773-0)] [[code](https://github.com/steineggerlab/foldseek)] [[server](https://search.foldseek.com/search)]
- MMseqs2 [[paper](https://www.nature.com/articles/nbt.3988)] [[code](https://github.com/soedinglab/MMseqs2)]
- MutCompute [[paper](https://www.nature.com/articles/s41587-022-01625-3)] [[server](https://mutcompute.com/)]

## Relevant Literature
- [Improving de novo protein binder design with deep learning](https://www.nature.com/articles/s41586-023-06758-2), Bennett et al., 2023
- [Improving Protein Expression, Stability, and Function with ProteinMPNN](https://www.biorxiv.org/content/10.1101/2024.01.08.574661v1), Sumida et al., 2024

## Miscellaneous
- Lecture on Protein Structure Prediction from the BioML Society [[video](https://www.youtube.com/watch?v=uAIuA1O7iE8)] [[slides](https://docs.google.com/presentation/d/1Wy-ePQqDHxEFDQqLHBhXJlKlLLZxZJxf/edit?usp=sharing&ouid=116893289863885570912&rtpof=true&sd=true)]
- Lecture on Protein Design from the BioML Society [[video](https://www.youtube.com/watch?v=uAIuA1O7iE8)] [[slides](https://docs.google.com/presentation/d/1Wy-ePQqDHxEFDQqLHBhXJlKlLLZxZJxf/edit?usp=sharing&ouid=116893289863885570912&rtpof=true&sd=true)]

## Biology

### CAR-T Literature
- [CAR T cell immunotherapy for human cancer](https://www.science.org/doi/10.1126/science.aar6711), June et al., 2018
- [Recent advances and discoveries in the mechanisms and functions of CAR T cells](https://www.nature.com/articles/s41577-020-00483-x), Larson et al., 2021
- [CAR-T design: Elements and their synergistic function](https://www.nature.com/articles/s41392-020-00262-z), Jayaraman et al., 2020

### CD20 Literature
- [The Development of a Recombinant scFv Monoclonal Antibody Targeting Canine CD20 for Use in Comparative Medicine](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4962021/), Jain et al., 2016
- [CD20 as a target for therapy](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2442443/), Winiarska et al., 2007
- [Mechanisms of action of CD20 antibodies](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3133688/), Boross et al., 2012
- [Rituximab (monoclonal anti-CD20 antibody): mechanisms of action and resistance](https://www.nature.com/articles/1206939), Smith 2003

## Compute Resources

### TACC
- Brief Overview
- Vista Supercomputer
- Containers on TACC
- Login
- [Create account](https://portal.tacc.utexas.edu/account-request)

### Other Compute
- [Synthia](https://www.synthia.ai/) - LLM that uses protein design tools
- Tools suggested by Adaptyv Bio