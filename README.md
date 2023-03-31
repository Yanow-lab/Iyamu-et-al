# Tutorial

## Plotepiscan: Plots epitope mapping data from peptide arrays onto protein sequence alignments.

Peptide arrays of overlapping sequences can be used to identify the epitope of antibodies on a protein antigen at the amino acid level. Scanning with monoclonal antibodies for molecular recognition using peptide arrays is particularly useful for epitope identification. The `plotepiscan.md` tool visualizes the data from epitope mapping screens using a color-coded sequence alignment representation of the antigens screened. In this study, we performed epitope mapping using the monoclonal antibody 3D10 against immobilized arrays of overlapping peptides. We targeted the extracellular domain of two alleles of the *Plasmodiun falciparum* virulence factor VAR2CSA. Arbritary units (AU) of fluorescence intensity quantified the antibody recognition for each peptide on the peptide array. 
* `plotepiscan.md` normalizes the intensity data using either linear or non-linear data transformations. Several data transformations are available through the function `ArrayTools.data_transform()`. The power law with cubic exponent was the transformation implemented in this study.
* `plotepiscan.md` maps the array data onto a global alignment between the FCR3 and NF54 alleles highlighting the antibody recognition on each peptide array by a colour code, from red to white for high to low intensity, respectively. When plotting the alignment, the colored boxes on the background of the symbols corresponds to the flourescence intensity score for the corresponding 20-mer peptide that ends at that residue.
To quantify the color map representing the intensity data, a color bar is generated according to the transformation applied. 

## Input files

`plotepiscan.md` takes as input two files containing the peptide array data. The data should be a table in comma-separated value (csv) format, e.g. a spreadsheet tab exported in csv format. Each row should represent data for one peptide. Header or footer lines should be deleted. The files used in this study are "FCR3_10ug.csv" and "NF54_10ug.csv". Peptides were 20 amino acids in length with an overlap of 19 and 18 amino acids for the FCR3 or NF54 arrays, respectively. Signal intensities for each peptide were based on the average of the median foreground intensities of each spot. 

`plotepiscan.md` requires a sequence file containing the sequences of the two VAR2CSA alleles screened as peptide arrays. The sequence file should be a txt file in FASTA format with header lines indicating the sequence label. The file used in this study: "Array_Seq.txt"

## Required dependencies

* Python package [Biotite](https://www.biotite-python.org) , version 0.35.0.
* [Pandas](https://pandas.pydata.org), version 1.5.2
* Python modules `ArrayTools.py` and `SignalArray.py`, available on this repository.

## Graphical abstract

<p align="center">
 <img src="https://user-images.githubusercontent.com/123925494/228916143-e61c398f-4578-4516-9691-124c7294185c.png" width="390" height="150">
</p>


