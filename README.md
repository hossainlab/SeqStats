# `SeqStats` - A Bioinformatics Package for DNA Sequence Analysis 

## Features 
- Reading Files 
    - Text Format 
    - FASTA 
    - FASTQ 
- Summary of Sequence 
    - Sequence Length 
    - GC Content 
    - AT Content 
    - Frequency 
    - Percentage of Nucleotide 
- DNA Sequence Frequency 
    - Bar Charts 
    - Pie Charts 
- Protein Sequence Frequency
    - Bar Chart 
    - Pie Charts 
- Calculating GC Content
- Calculating AT Content
- SLiding Window Analysis
    - GC Sliding Window Analysis 
    - AT Sliding Window Analysis 
- Protein Sequence Frequency 
- Nucleotide Density 
    - Nucleotide Density one
    - Nucleotide Density Two

## Install Package 
```bash 
pip install SeqStats
```
- [SEE MORE](https://pypi.org/project/SeqStats/)

## Import 
```python
import SeqStats as ss  
```

## Reading File 
```python
ss.readFASTA('file_path')
```

## Sequence Informations 
```python
ss.seqInfo(seq)
```

## To know more about `SeqStats` llease check out the tutorial notebook 

<h3>About the Author</h3>
This repo was created by <a href="https://datatutor.github.io/" target="_blank">Jubayer Hossain</a> <br>
<a href="https://datatutor.github.io/" target="_blank">Jubayer Hossain</a> is a student of Microbiology at Jagannath University and the founder of <a href="https://hdro.github.io/" target="_blank">Health Data Research Organization</a>. He is also a team member of a bioinformatics research group known as Bio-Bio-1. 

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.