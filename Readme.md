# Readme
This repository contains source-code for simulations and analyses of transcription-looping, associated with our manuscript:

## How subtle changes in 3D structure can create large changes in transcription 
[![DOI](https://zenodo.org/badge/299465350.svg)](https://zenodo.org/badge/latestdoi/299465350)

Jordan Xiao, Alistair Boettiger
Stanford University, Program in Biophysics, Department of Developmental Biology

### Abstract: 
Animal genomes are organized into topologically associated domains (TADs), which exhibit more intra-domain than inter-domain contact. However, the absolute difference in contact is usually no more than twofold, even though disruptions to TAD boundaries can change gene expression by 8-10 fold. Existing models fail to explain the superlinear transcriptional response to genomic contact. We propose a futile cycle model where an enzyme stimulated by association with its products can exhibit bistability and hysteresis, allowing a small increase  in enhancer-promoter contact to produce a large change in expression without obvious correlation between E-P contact and promoter activity. Through mathematical analysis and stochastic simulation, we show that this system can create an illusion of enhancer-promoter specificity through hypersensitivity. It also potentially reconciles global cohesin loop disruption and TAD boundary deletion experiments which have led to conflicting conclusions on whether E-P contact affects gene regulation.

## Repository Organization
Scripts for all simulations described in the manuscript can be found here, the script titles match the corresponding figure names.

Stochastic simulations were written in Matlab(™) 2020a. Executable source code to run these simulations is available here. A detailed table of parameters values and parameter descriptions is provided at the start of each script. The associated figure is indicated in the filename and in the header of the corresponding script for each simulation. In brief, these simulations model the chemical master equation in explicit time and discrete time steps, with discrete transition probabilities for the addition and removal of promoter tags. For simplicity, we equate promoter modification and transcription, reducing the total parameter space without loss of generality. 

Caution: some of the simulations are computationally expensive. While some care has been taken for  efficiency, our simulations were run on a Dell T630 with 500 GB of RAM and 28 2.4GHz cores.