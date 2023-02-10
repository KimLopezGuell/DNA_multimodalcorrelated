# Multimodal correlated DNA model

#### Coarse-grain model of DNA flexibility with multimodality and correlations between base pairs.

The folder contains all scripts, additional functions and inputs. 

"main_di_tetra.c" is the main file, and "dna_flex.h" is the file that contains the majority of the main code functions. Other additional utilities are "analysis.c", "energy.h", "nrutil.c", "rebuild.c", "cmn_fncs.c". \
Inputs used for the discussion of the paper can be found in he folder Inputs. 


## Dependencies

-----------------
The program SCHNArP [1] is needed for reconstruction of the DNA sequence to cartesian coordinates. 

The program X3DNA [2] is needed for the analysis and rebuilding of the 3D structures.


## Instructions for running

-----------------

If any script is changed, running makefile will compile the program again. 

To run the program, just type: 

`export PATH="[PATH]/SCHNArP/bin":'[PATH]/X3DNA/bin/':$PATH` \
`export SCHNA_AR_P="[PATH]/SCHNArP"` \
`export X3DNA="[PATH]/X3DNA/bin"` 

where [PATH] refers to whichever folder you have the programs saved in your computer. 

`cd [PATH where you have the scripts]` 

`mkdir Result_folder` \
`mkdir Result_folder/analysis`  \
`mkdir Result_folder/output`  \
`mkdir Result_folder/output_acc`  \
`mkdir Result_folder/output_cart`  \
`mkdir Result_folder/output_dnaflex`  \
`mkdir Result_folder/output_schnarp`  \
`mkdir Result_folder/output_tables_helpar` 

where "Result_folder" can have any name you desire. 

Finally, the actual program is run by: \
`./MC_Chromatin sequence.dat 5000 [PATH]/Result_folder 6` 

where sequence.dat has the sequence of interest (for example, "AATTAGCTA"), 5000 refers to the number of runs desired for the ensemble of outputs, and the last integer is chosen to select the type of model desired (6 stands for multimodal correlated, other types can be chosen**). 

** 0 = individual, 1 = dimer, 2 = tetramer, 3 = tetramer with dimer end; 4 = tetramer with dimer end (multimodal); 5 = tetramer with dimer end (unimodal, with correlations between bp i and i+1, i and i-1); 6 = tetramer with dimer end (multimodal, with correlations whenever possible (bp 3 to itot-4), computed in the MC algorithm with 12x12 matrices and in selection process before MC with an Ising approach)

------------------

## References

[1] https://doi.org/10.1006/jmbi.1997.1345. \
[2] https://doi.org/10.1038/nprot.2008.104


------------------
## Zenodo DOI


[![DOI](https://zenodo.org/badge/590646462.svg)](https://zenodo.org/badge/latestdoi/590646462)




