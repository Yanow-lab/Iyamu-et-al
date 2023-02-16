---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.4
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

```python tags=[]
import biotite.sequence as seq
import biotite.sequence.align as align
import biotite.sequence.graphics as graphics
import biotite.sequence.io.fasta as fasta
```

```python
# Read avidin(avi) and streptavidin(savi) sequences from FASTA format txt-file:
fasta_file = fasta.FastaFile.read('Avi-sAvi.txt')
```

```python
#Parse protein sequences of Avi and sAvi
for name, sequence in fasta_file.items():
    if "P02701" in name:
        avidin_seq = seq.ProteinSequence(sequence) # create a protein seq obj
    elif "B8YQ01" in name:
        streptavidin_seq = seq.ProteinSequence(sequence)
```

```python
# Get BLOSUM62 matrix
matrix = align.SubstitutionMatrix.std_protein_matrix()
# Perform pairwise sequence alignment with affine gap penalty
# Terminal gaps are not penalized
alignments = align.align_optimal(avidin_seq, streptavidin_seq, matrix, # returns list of ali
                    gap_penalty=(-10, -1), local=False, terminal_penalty=False)
```

```python
A = alignments[0]
```

```python
print(A)
```

#### Get the of symbols the aligned sequences

```python
traceA = align.get_symbols(A)[0]
#traceA
```

```python
traceB = align.get_symbols(A)[1]
#traceB
```

## Load(create) pep-array data

```python
from itertools import islice
import pandas as pd
```

```python
##  Reader for FCR3_like sequences: 10 mer peptides with 9 overlap
##  formula: No of peptides = {[len(seq)-len(pep)]/step} + 1
with open("Avi-sAvi.txt", 'r+') as f:
    sec = ''
    a = 0
    b = 10
    for line in islice(f, 1, 4):  # READ lines of Seq0 only
        sec += line
    sec=sec.replace('\n','')    # remove '\n' character 
    pep_s0 = [sec[i:i+10] for i in range(len(sec)-9)] # list comprehension
    f.close
```

```python
len(pep_s0)
```

```python
##  Reader for NF54_like sequences: 10 mer peptides with 8 overlap 
##  formula: No of peptides = {[len(seq)-len(pep)]/step} + 1 
with open("Avi-sAvi.txt", 'r+') as f:
    sec = ''
    a = 0
    b = 10
    for line in islice(f, 5, 8): # READ lines of Seq1 only
        sec += line
    sec=sec.replace('\n','')
    pep_s1 = [sec[i:i+10] for i in range(0,(len(sec)-9),2)]
    f.close
```

```python
len(pep_s1)
```

#### Array data for MAb reacognition follows a long-tailed distrubution. Shape: pareto

```python
import numpy as np
```

```python
a, m = 0.5, 2. 
np.random.seed(42) # seed the randomizer method to get reproducible results
```

```python
ag1_scan = pd.DataFrame()
ag1_scan['Seq'] = pep_s0
ag1_scan['r1'] = (np.random.pareto(a, ag1_scan.shape[0]) + 1) * m
ag1_scan['r2'] = (np.random.pareto(a, ag1_scan.shape[0]) + 1.1) * m # replicate_variability
```

```python
score_res = ag1_scan['Seq'].str[-1]
ag1_scan['s_res'] = score_res
```

```python
ag2_scan = pd.DataFrame()
ag2_scan['Seq'] = pep_s1
ag2_scan['r1'] = (np.random.pareto(a, ag2_scan.shape[0]) + 1) * m
ag2_scan['r2'] = (np.random.pareto(a, ag2_scan.shape[0]) + 1.1) * m # replicate_variability
```

```python
score_res = ag2_scan['Seq'].str[-1]
ag2_scan['s_res'] = score_res
```

## Process array data

```python
import SignalArray as sa
```

```python
dfa = sa.compute_params(ag1_scan, combine = 'max', flag_noisy = True)
dfb = sa.compute_params(ag2_scan, combine = 'max', flag_noisy = True)
```

```python
sa.data_describe(dfa) # Define a threshold 
```

```python
sa.data_transform(dfa, method ='cubic', threshold = 0)
sa.data_transform(dfb, method ='cubic', threshold = 0)
```

## Convert a list of score residues from the epitope </br>scan data into a aligment-like gapped sequences 

```python
gapd_s1 = sa.gapped_seq(dfa, traceA,10, 1)
gapd_s2 = sa.gapped_seq(dfb, traceB,10, 2) # here the overlap step is 2 
                                        # (peptide = 20-mer with 18 overlap)
```

## Disparo de control

```python
len(gapd_s1) == len(gapd_s2)
```

## Create a signal_map (ndarray)

```python
score = sa.signal_map(gapd_s1, gapd_s2,)
```

## Plot

```python
import ArrayTools as at
import matplotlib.pyplot as plt
import matplotlib as mpl
```

```python
fig = plt.figure(figsize=(8, 2.5))
ax1 = fig.add_subplot(111)
at.plot_alignment_array(ax1, alignments[0], fl_score= score,
     labels=["Avi", "sAvi"],show_numbers=True,
    symbols_per_line= 50, show_line_position=True) 

# add a 2nd axes and a colorbar

ax2 = fig.add_axes([0.1,-0.15,0.8,0.1])
ax2.set_frame_on(False)
cmp = at.get_cmap(ax2, score)
cbar = at.get_colorbar(ax2, dfa, dfb, cmp, transform = 'cubic', 
                       orient = 'horizontal', title = 'signal intensity')

# to improve readability, tilt ticklabels on the colorbar

labls = cbar.ax.get_xticklabels()
plt.setp(labls, rotation=45, horizontalalignment='center')

plt.show()
```

```python

```
