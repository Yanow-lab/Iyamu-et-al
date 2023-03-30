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

```python
import matplotlib.pyplot as plt
import biotite.sequence as seq
import biotite.sequence.align as align
import biotite.sequence.io.fasta as fasta
```

```python
fasta_file = fasta.FastaFile.read('Array_Seq.txt')
```

```python
# Parse protein sequences of FCR3 and NF54
for name, sequence in fasta_file.items():
    if "AAQ73926" in name:
        FCR3_seq = seq.ProteinSequence(sequence)
    elif "EWC87419" in name:
        NF54_seq = seq.ProteinSequence(sequence)
```

```python
# Get BLOSUM62 matrix
matrix = align.SubstitutionMatrix.std_protein_matrix()
# Perform pairwise sequence alignment with affine gap penalty
# Terminal gaps are not penalized
alignments = align.align_optimal(FCR3_seq, NF54_seq, matrix,
                                 gap_penalty = (-10, -1), 
                                 terminal_penalty = False)
```

## Load epitope scan data

```python
import SignalArray as sa
```

```python
filenames = ['FCR3_10ug.csv','NF54_10ug.csv']
```

```python
d = 0
for f in filenames:
    if f == filenames[0]:
        ag1_scan = sa.read_scan(filenames[d], 20, 20)
    elif f == filenames[1]:
        ag2_scan = sa.read_scan(filenames[d], 20, 20)
    d = d + 1
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
A = alignments[0]
traceA = align.get_symbols(A)[0]
traceB = align.get_symbols(A)[1]
```

```python
gapd_s1 = sa.gapped_seq(dfa, traceA, 20, 1)
gapd_s2 = sa.gapped_seq(dfb, traceB, 20, 2) # overlap_step: 2 (pep = 20-mer with 18 overlap)
```

## Checkpoint

```python
len(gapd_s1) == len(gapd_s2)
```

## Create a signal_map (ndarray)

```python
score = sa.signal_map(gapd_s1, gapd_s2)
```

## Plot

```python
import ArrayTools as at
```

```python
fig = plt.figure(figsize=(20, 16))
ax = fig.add_subplot(111)
at.plot_alignment_array(
    ax, alignments[0], fl_score = score, labels = ["FCR3", "NF54"],
    show_numbers = True, symbols_per_line = 120,
    show_line_position = True) 

# add a 2nd axes and a colorbar

ax2 = fig.add_axes([0.1,-0.005, 0.8, 0.03])
ax2.set_frame_on(False)
cmp = at.get_cmap(ax2, score)
cbar = at.get_colorbar(ax2, dfa, dfb, cmp, transform = 'cubic', 
                       orient = 'horizontal', 
                       title = 'Fluorescence Intensity [AU]')

plt.show()
```

```python
fig.savefig('cubic_cbar.svg', transparent=True, dpi=80, bbox_inches="tight")
```

```python

```
