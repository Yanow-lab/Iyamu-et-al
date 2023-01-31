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
import matplotlib.pyplot as plt
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
alignments = align.align_optimal(avidin_seq, streptavidin_seq, matrix, # returns a list of alig
                                 gap_penalty=(-10, -1), local=False, terminal_penalty=False)
```

```python
A = alignments[0]
```

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

```python
import numpy as np
```

```python
dfa = pd.DataFrame()
dfa['seq'] = pep_s0
dfa['score'] = np.random.uniform(0.1,1.0, dfa.shape[0])
```

```python
dfb = pd.DataFrame()
dfb['seq'] = pep_s1
dfb['score'] = np.random.uniform(0.1,1.0, dfb.shape[0])
```

```python
dfa
```

## Process array data

```python
score_res = dfa['seq'].str[-1]
dfa['s_res'] = score_res
```

```python
score_res = dfb['seq'].str[-1]
dfb['s_res'] = score_res
```

```python
dfb.head()
```

```python
dfb.head()
```

```python
gapd_s1 = list(zip(dfa.s_res , dfa.score)) # convert the df columns into a list of tuples for Seq_0
lk1 = dfa["s_res"].values.tolist()     # convert the df column into a list of keys for Seq_0
```

```python
gapd_s2 = list(zip(dfb.s_res , dfb.score)) 
lk2 = dfb["s_res"].values.tolist()
```

```python
len(A.trace)
```

```python
print(A)
```

```python
######    Seq_0 create a list of tuples same length as trace[seq] #######
x=0
b=0
c=0 #cyclic counter up to the peptide length :10
p=0 #peptide counter
for b in range(len(lk1)):
    for a in traceA[x:]:
        if c < 9 : #and x < (len(traceA)-1):
            if a==None:
                gapd_s1.insert(x,(traceA[x],0)) 
                x=x+1
            elif a != lk1[b] :
                gapd_s1.insert(x,(traceA[x],0))         
                x=x+1
                c=c+1
            elif p==0:
                gapd_s1.insert(x,(traceA[x],0)) 
                x=x+1
                c=c+1 
            else:
                x=x+1
                c=c+1 
                break
        else:
            c = 0 # reset the counter        
            p=p+1
            x=x+1
            break

if len(gapd_s1) < len(traceA) and traceA[len(gapd_s1)+1]== None:
    gapd_s1_tail=[]
    for n in range(len(traceA)-len(gapd_s1)):
        gapd_s1_tail.append(('None', 0))
    gapd_s1 = gapd_s1+gapd_s1_tail


```

```python
####     Seq_1 create a list of tuples same length as trace[seq]  #####
x=0
b=0
c=0 #cyclic counter up to the peptide length :10
p=0 #peptide counter
for b in range(len(lk2)):
    for a in traceB[x:]:
        if c < 9 and p==0:            
            if a==None :
                gapd_s2.insert(x,(traceB[x],0)) 
                x=x+1
            else:
                gapd_s2.insert(x,(traceB[x],0))         
                x=x+1
                c=c+1
        elif p==0 :
            c = 0 # reset the counter        
            p=p+1
            x=x+1
            break
        if p!=0: 
            if a==None and c == 0:
                gapd_s2.insert(x,(traceB[x],0)) 
                x=x+1
            elif c % 2 == 0 : 
                if a==None:
                    gapd_s2.insert(x,(traceB[x],0)) 
                    x=x+1
                else:
                    gapd_s2.insert(x,(traceB[x],0)) 
                    x=x+1
                    c=c+1
            elif c % 2 != 0 : 
                if a==None:
                    gapd_s2.insert(x,(traceB[x],0)) 
                    x=x+1
                elif a != lk2[b]:
                    gapd_s2.insert(x,(traceB[x],0))         
                    x=x+1
                    c=c+1
                else:        
                    x=x+1
                    c=c+1
                    break
         


        
        
        
        
        
if len(gapd_s2) < len(traceB) and traceB[len(gapd_s2)+1] == None:
    gapd_s2_tail=[]
    for n in range(len(traceB)-len(gapd_s2)):
        gapd_s2_tail.append(('None', 0)) 

    gapd_s2 = gapd_s2+gapd_s2_tail               
```

## Disparo de control

```python
len(gapd_s1) == len(gapd_s2)
```

##

```python
gapd_s1[:10]
```

```python
import numpy as np
```

```python
#alignment.trace.shape[0] = seq_len
fl_score=np.zeros((len(gapd_s1),2))
for v1 in range(len(gapd_s1)):
    fl_score[v1,0]=gapd_s1[v1][1]    
    fl_score[v1,1]=gapd_s2[v1][1]
```

```python
score = fl_score

```

```python
def _get_signal( score, column_i, seq_i):
    if fl_score is None:
        signal = 0.01
    else:
        signal = fl_score[column_i, seq_i]
    return signal
```

```python
_get_signal(score, 3,0)
```

## Plot

```python
symbol_plotter = ap.ArrayPlotter(ax, score)
```

```python
import ArrayPlotter as ap

""" my module ArrayPlotter2 has a class and a function
    that uses objects from that class
"""
fig = plt.figure(figsize=(8, 2.5))
ax = fig.add_subplot(111)
ap.plot_alignment_array(ax, alignments[0], fl_score= score,
     labels=["Avi", "sAvi"],show_numbers=True,
    symbols_per_line= 50, show_line_position=True) 

fig.colorbar(symbol_plotter,ax=ax)
fig.tight_layout()

plt.show()
```

```python

```

```python
stop
```

```python
# DEBUG_Seq_0 create a list of tuples same length as trace[seq]
x=0
b=0
c=0 #cyclic counter up to the peptide length :10
p=0 #peptide counter
for b in range(len(lk1)):
    print('b-',b)
    for a in traceA[x:]:
        print('entry a-',a)
        if c < 9 : #and x < (len(traceA)-1):
            if a==None:
                gapd_s1.insert(x,(traceA[x],0)) 
                print('a==',a,a==None)
                print('x',x)
                print('<-insert01')
                x=x+1
            elif a != lk1[b] :
                gapd_s1.insert(x,(traceA[x],0))         
                print('a',a,'b',b)                   
                print('x',x)
                print(a != lk1[b], '<-insert02')
                x=x+1
                c=c+1
                print('c',c)
            elif p==0:
                gapd_s1.insert(x,(traceA[x],0)) 
                print('a',a,'b',b)                   
                print('x',x)
                print('<-insert03')
                x=x+1
                c=c+1 
                print('c',c)
            else:
                print('a',a,'b',b)                   
                print('x',x)
                print(a != lk1[b],'<-val-0')
                x=x+1
                c=c+1 
                print('c',c)
                break
        else:
            c = 0 # reset the counter        
            print('a',a,'b',b)                   
            print('x',x)
            print(a != lk1[b],'<-val-2')
            p=p+1
            print('P',p)
            x=x+1
                #c=c+1
            break

if len(gapd_s1) < len(traceA) and traceA[len(gapd_s1)+1]== None:
    gapd_s1_tail=[]
    for n in range(len(traceA)-len(gapd_s1)):
        gapd_s1_tail.append(('None', 0))
    gapd_s1 = gapd_s1+gapd_s1_tail


```

```python
#  DEBUG_2_Seq_1 (step:2)create a list of tuples same length as trace[seq]
x=0
b=0
c=0 #cyclic counter up to the peptide length :10
p=0 #peptide counter
for b in range(len(lk2)):
    print('b-',b)
    for a in traceB[x:]:
        print('entry a-',a) 
        if c < 9 and p==0:            
            if a==None :
                gapd_s2.insert(x,(traceB[x],0)) 
                print('a==',a,a==None)
                print('x',x)
                print('<-insert01')
                x=x+1
            else:
                gapd_s2.insert(x,(traceB[x],0))         
                print('a',a,'b',b)                   
                print('x',x)
                print(a != lk2[b], '<-insert02')
                x=x+1
                c=c+1
                print('c',c)
        elif p==0 :
            c = 0 # reset the counter        
            print('a',a,'b',b)                   
            print('x',x)
            print(a != lk2[b],'<-val-1')
            p=p+1
            print('P',p)
            x=x+1
                #c=c+1
            break
        if p!=0: 
            if a==None and c == 0:
                gapd_s2.insert(x,(traceB[x],0)) 
                print('a==',a,a==None)
                print('x',x)
                print('<-insert03')
                x=x+1
            elif c % 2 == 0 : 
                if a==None:
                    gapd_s2.insert(x,(traceB[x],0)) 
                    print('a',a)
                    print('x',x)
                    print('<-insert04')
                    x=x+1
                else:
                    gapd_s2.insert(x,(traceB[x],0)) 
                    print('a',a)
                    print('x',x)
                    print('<-insert05')
                    x=x+1
                    c=c+1
            elif c % 2 != 0 : 
                if a==None:
                    gapd_s2.insert(x,(traceB[x],0)) 
                    print('a==',a,a==None)
                    print('x',x)
                    print('<-insert06')
                    x=x+1
                elif a != lk2[b]:
                    gapd_s2.insert(x,(traceB[x],0))         
                    print('a',a,'b',b)                   
                    print('x',x)
                    print(a != lk2[b], '<-insert07')
                    x=x+1
                    c=c+1
                    print('c',c)
                else:        
                    print('a',a,'b',b)                   
                    print('x',x)
                    print(a != lk2[b], '<-insert07')
                    x=x+1
                    c=c+1
                    print('c',c)
                    break
         


        
        
        
        
        
if len(gapd_s2) < len(traceB) and traceB[len(gapd_s2)+1] == None:
    gapd_s2_tail=[]
    for n in range(len(traceB)-len(gapd_s2)):
        gapd_s2_tail.append(('None', 0)) 

gapd_s2 = gapd_s2+gapd_s2_tail               
```
