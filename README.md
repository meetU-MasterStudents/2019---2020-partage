# 2019---2020-partage
For exchanging material and doc


**** 
# Bibliography (Biblio)


## For profile-profile alignments and (remote) homology detection

ORION: fold recognition method based on profile-profile alignments

scoringProfiles: methods to score profile-profile alignments

HHsearch: remote homology detection using profiles (defined from hidden Markov models, not from PSSM!)

CLADE: remote homology detection using a multi-profile strategy

## For evaluating the quality of a 3D model

DOPE: "historical" statistical potential for evaluating the quality of 3D models (does not perform very well but is a good start)

Rosetta: a reference function for estimating the energy of a 3D protein conformation (be aware it's all-atom!)

SBROD: coarse-grained statistical potential to evaluate 3D models quality (can be applied to Calpha-only or backbone-only structures, recently performed well in CASP13, can be completely re-trained!) 

## For structural annotations (secondary structure, solvent accessibility, contacts... predicted from sequence)

ReviewSA: book chapter reviewing methods to predict secondary structure, solvent accessibility, torsional angles and contact maps, from sequenc information) 

CCMPRED: method to predict protein-protein contact by extracting coevolution signals (co-occurring patterns of mutations across sequences)

****
# Helper data and tools (Codes)

Params: directory containing values for the DOPE statistical potential. 

Tools: directory containing a PDB parser (in Python) and a script to weight sequences based on their similarity (in Perl).

****
# Data for training and testing (Data)


## 1009 directories corresponding to 1009 families from HOMSTRAD.

For each family:

- FASTA file containing the master (reference) sequence,

- MAP file containing a multiple sequence alignment with the master sequence and homologous (or related) sequences, 

- SCOP_ID file containing the SCOP identifier of the family,

- PDB file containing the 3D coordinates of the master sequence. Please note that the structure may contain "holes" (missing residues that either could not be resolved, or were modified/non-canonical). 

### Please note that the master sequence is named by its PDB code in the MAP file. It may not be the first sequence appearing in the file!

## 398 test sequences ("queries") to validate your program.

They are contained in the file queries398.multifasta. The name of each query sequence is as follows:

>gluts | Q09596.1 | NP_001254267.1 | 98.0%

- gluts: name of the family you should find!

- Q09596.1: UNIPROT code

- NP_001254267.1: sequence identifier (you do not care about it)

- 98.0%: max percentage of identity between the query and all sequences from the family you should find! It gives an idea of the level of difficulty associated to the query, but in principle you do not need to use this value in your calculations 

