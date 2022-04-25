# speciation-by-extinction

This repository contains data and scripts to reproduce analyses in "Speciation-by-Extinction".


######################################################################

## Cranioleuca

Data and scripts used for species delimitation analyses with simulated extinction in *Cranioleuca antisiensis*. Phenotypic and genetic data from Seeholzer and Brumfield 2018. 

Principal Directories

### bfd - Bayes Factor species delimitation of phenotypic data

Includes all raw data and files necessary to create .xml and .job scripts for BFD runs (ran on a cluster).

##### data
- cran.snp.data.genlight.Rdata

##### scripts
- SCRIPT1.genind2snapp.R - converts raw genetic data in genlight format (cran.snp.data.genlight.Rdata)
- SCRIPT2.xml.editing.R - edits base .job and .xml files for each BFD run
- SCRIPT3.collate.analysis.files.R - collates output files from subdirectories in cluster
- SCRIPT4.process.analysis.files.R - summarizes results of BFD analyses (used in SCRIPT.Fig.2.plot.R)
- SCRIPT4b.restart.failed.runs.R - restarts failed runs (for use on cluster)

### mclust - Normal Mixture Model species delimitation of phenotypic data

##### data
- cran.morph.raw.txt - raw Cranioleuca measurements
- cran.plumage.raw.txt - raw Cranioleuca plumage coordinates
- cran.meta.txt - meta data for each Cranioleuca specimen

##### scripts
- SCRIPT1.cran.morph.processing.R - process raw morphological measurements for analysis
- SCRIPT2.cran.plumage.processing.R - process raw plumage coordinates for analysis
- SCRIPT3.mclust.analysis.R - mclust species delimitation analysis


- Seeholzer GF, Brumfield RT. 2017. Isolation-by-distance, not incipient ecological speciation, explains genetic differentiation in an Andean songbird (Aves: Furnariidae: *Cranioleuca antisiensis*, Line-cheeked Spinetail) despite near three-fold body size change across an environmental. Molecular Ecology 1–18.


## Geospiza

Data and scripts used for species delimitation analysis of simulated extinction in the morphological dataset for the Darwin's Finch used by McKay & Zink 2015 and Cadena et al. 2018. Code adapted from supplementary material in Cadena et al. 2018. 

Data Geospiza.data.csv from Cadena et al. 2018. Institution codes as follows
RM = Rotheschild Mus.
Harris = Harris
Hull = G. W. Hull
CAS = California Academy of Sciences
SU = Stanford University


- Cadena, C. D., Zapata, F. & Jiménez, I. Issues and perspectives in species delimitation using phenotypic data: Atlantean evolution in Darwin’s Finches. Syst. Biol. 67, 181–194 (2018).
- Mckay, B. D. & Zink, R. M. Sisyphean evolution in Darwin’s finches. Biol. Rev. 90, 689–698 (2015).


## Bird Subspecies

Data and script to calculate the number of bird species with three or more subspecies based on the Clements/eBird taxonomy.

- Clements, J. F. et al. The eBird/Clements checklist of birds of the world: v2019. (2019)
