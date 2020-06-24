# `seqstats` - A Bioinformatics Package for Sequence Statistics

## Load `seqstats` Packages


```python
# generic import
import seqstats
```

```python
# import specific function
from SeqStats import readText, readFASTA, readFASTQ
```


```python
# import as a nickname
import SeqStats as ss # here `ss` is a nickname
```

## Getting Help


```python
help(ss.readFASTA)
```

    Help on function readFASTA in module SeqStats:

    readFASTA(inputfile)
        Reads a FASTA file and skips the sequence information line with special character removed"
        Args:
            inputfile ([type]): [FASTA File]
        Returns:
            [type]: [String]



## Files Handing
- Text File
- FASTA File
- FASTQ File


```python
# reading text file
text = ss.readFASTA("./data/dna.txt")
```


```python
# reading FASTA file
fasta = ss.readFASTA("./data/Haemophilus_influenzae.fasta")
```


```python
# reading FASTQ
seqs, quals = ss.readFASTQ("./data/example.fastq")
```

## DNA Sequence Statistics


```python
# read a FASTA sequence
seq = readFASTA('./data/Haemophilus_influenzae.fasta')
```

### Head  


```python
# examine first few base pairs(1000bp)
ss.head(seq)
```




    'AACCGAAATTACAGTGCATGGACGCACAAAATCTGATGGTTATCGTGCTGATAGAATTAATTGGAAAAAAATTGGTAAAGTCCGAGAGCGTTTATCCATTCCTGTTATTGCTAACGGAGAAATTTGGCATTGGCAAGATGGTCAAGATTGCTTATCTCAAACAGGTTGTCAGGATTTAATGGTGGGACGAGGTGCATTGAATATTCCGAACTTAAGCCATGTTCTGAAATCAAATGCAGAAAAAATGCCTTGGAATGAGATTCAAAAAATCTTGCAAAAATATGCGAATGTTGAAAATGAATATGGCAGCGGTTTTTACCATGTGGCACGAATTAAACAATGGTTACGTTATTTGAATAAGGAATATGATGAGGCGAACCAAGAGTTTGATAAGATTAAGACTTGCCAAACTGCTGAAGATTTGAAATTACGGTTAAATGATAAATAAAAAACCTGCTAATCAGCAGGTTTTCTTTTTCTAAATTATTTAAAAATTCACCGCACTTTTTAATTTTTTCAGTGCATTAGTTTCTAATTGACGAATACGTTCCGCAGATACATTATATTTTGCTGCCAAGTCGTGCAATGTAGCTTTGTTATCATCTAACCAACGAGCTTTGATAATATCCTGACTACGCGCATCTAGACTTTGTAGTGCCGCACCTAATTGCTCGGTTGCTTGGCTTTCAAAGTTTTCATTTTCAAGCGCTGCTGCGAAATTAGAACTTTTATCTTCCAAATAAAGGGCTGGAGAATAGGTCTCTGTTTCTGCATCATCAGTAGGTAAATCAAAACCGACATCGGCACCGCTCATACGCGATTCCATTTCGATCACATCTTCTTTGGAAACACCCAATTCATTTGCGACCATATCAACTTCATTTTCATTGAACCAACCTAAACGTTGTTTAGTTTTACGTAGGTTGAAAAATAATTTACGTTGAGCTTTCGTGGTTGCGACTTTTACGATACGCCAGTTACGAAGAACATATTCGTGGAT'



### Tail


```python
# examine last few base pairs(1000bp)
ss.tail(seq)
```




    'CTAATCAATTAATGATTCAGCTTAAACTGGATGATGCGGGGCGTAAATTAGCGCAGGATGCGTTTCGTCGTGGTAAAGAATCTGATTTCCCAATTCGCCAAGTTATCCGTGAATTTCGCATTGGTTGTGGGCAACGTGCTGATTTATTGCGAATGTTCTTACAAGTGCAAGTTCAAGCGGCTTTTGCCGATTCAGAACTTCACGAAAATGAGAAAGAAGTGCTTTATGTGATTGCTGAAGAACTTGGTCTTTCTCGTATGCAATTTGAACAAATGATTGCGATGGAAATGGCTGCTCGAGCTTTTACTCAAGGCAGTTTTTATCAGCAATATCAACAAGGTGCCTATCAAGGTGGCTACCAATATCAGCAACAAAATAGTGGTGGTTACCAACATGCTTCTGGCCCAACTTTAAATGATGCCTATAAAGTATTGGGTGTAACTGAGTCCGACGAGCAAAATACGGTTAAGCGTGCTTATCGTCGTCTAATGAATGAACACCATCCAGATAAACTCGTGGCGAAAGGTTTACCGCCAGAAATGATGGAAATGGCAAAAGAAAAAACGCAACAGATTCAAGCTGCTTACGATTTAATTTGTAAAGCAAAAGGCTGGAAATAGTGCGTGTTATTCTTGCGCCTATGCAAGGGGTTCTTGATCCCTTTGTACGCCAACTTCTCACTGAAGTGAATGACTACGATTTATGTATAACAGAATTTGTCCGCGTAGTTGATCAACTTCTTCCTGAAAAAGTATTTTATCGTTTATGCCCTGAATTAAAAAATCAGGGCTTTACTTCTTCTAGCACGCCTGTGCGAGTGCAGTTGCTAGGGCAGCATCCAGAATACCTTGCTGAAAATGCAATTCGTGCAATCGAGCTTGGTTCTCATGGCATTGATTTAAATTGTGGTTGCCCTTCTAAAACAGTGAATGGCAGCAATGGCGGTGCGGCATTATTGAAACAGCCTGAATTGATTTATCGTGCAACTCAAGCCTTACGC'



### Summary


```python
# basic info
ss.seqInfo(seq)
```

    Sequence Length: 1856176 bp
    Frequency: {'G': 349480, 'C': 356835, 'A': 576078, 'T': 573783}
    Percentage: {'G': 0.19, 'C': 0.19, 'A': 0.31, 'T': 0.31}
    GC Content: 38.05
    AT Content: 61.95


### Frequency


```python
# frequency: look at the four bases A, T, C, G only
ss.frequency(seq)
```




    {'A': 576078, 'T': 573783, 'G': 349480, 'C': 356835}




```python
# frequency: look at every letter that appears in the sequence.
ss.frequency(seq, useall=True)
```




    {'G': 349480, 'C': 356835, 'A': 576078, 'T': 573783}




```python
# calculate percentage
ss.frequency(seq, calpc=True)
```




    {'A': 0.31, 'T': 0.31, 'G': 0.19, 'C': 0.19}



### Nucleotide Frequency Plotting


```python
# plotting frequency
ss.plotFrequency(seq, kind='bar')
```


![png](output_24_0.png)



```python
# plotting frequency
ss.plotFrequency(seq, kind='pie')
```


![png](output_25_0.png)


## GC & AT Calculation


```python
# calculate GC
ss.calculateGC(seq)
```




    38.05




```python
# calculate AT
ss.calculateAT(seq)
```




    61.95



## Sliding Window Analysis


```python
help(ss.subSeqGC)
```

    Help on function subSeqGC in module SeqStats:

    subSeqGC(seq, window=1000)
        # 7. Function for Calculating Sub-Sequence GC Content




```python
# subseq GC
gc = ss.subSeqGC(seq, 9000)
gc[100:121]
```




    [36.59,
     34.62,
     37.92,
     37.42,
     39.32,
     36.89,
     36.16,
     37.9,
     36.57,
     37.9,
     39.17,
     39.6,
     40.2,
     42.0,
     37.97,
     38.34,
     37.01,
     40.12,
     39.09,
     38.98,
     41.54]




```python
# subSeq AT
at = ss.subSeqAT(seq)
at[1:20]
```




    [63.4,
     61.7,
     62.4,
     63.9,
     62.7,
     64.9,
     65.9,
     62.9,
     63.6,
     62.3,
     62.3,
     66.1,
     63.0,
     60.2,
     74.1,
     70.5,
     66.3,
     63.1,
     66.9]



## Sliding window Plot


```python
# sling plot of GC
ss.plotSubSeqGC(seq, window=9000)
```


![png](output_34_0.png)



```python
# sling plot of AT
ss.plotSubSeqAT(seq, window=9000)
```


![png](output_35_0.png)


## K-mer Analysis

### Building K-mers


```python
# k-mer building
kmers = ss.buildKmers(seq, 2)
kmers[1:20]
```




    ['AC',
     'CC',
     'CG',
     'GA',
     'AA',
     'AA',
     'AT',
     'TT',
     'TA',
     'AC',
     'CA',
     'AG',
     'GT',
     'TG',
     'GC',
     'CA',
     'AT',
     'TG',
     'GG']



### K-mers Frequency


```python
# count k-mers
ss.countKmers(kmers)
```




    Counter({'AA': 223369,
             'AC': 93850,
             'CC': 69297,
             'CG': 73364,
             'GA': 94221,
             'AT': 169804,
             'TT': 221014,
             'TA': 134844,
             'CA': 123643,
             'AG': 89055,
             'GT': 92435,
             'TG': 120812,
             'GC': 96575,
             'GG': 66249,
             'TC': 97113,
             'CT': 90530})



### K-mer Frequency Plotting


```python
# plotting kmers
ss.plotKmerFrequency(kmers, kind='bar')
```


![png](output_42_0.png)



```python
# plotting kmers frequency
ss.plotKmerFrequency(kmers, kind='pie')
```


![png](output_43_0.png)


## Protein Sequence Statistics


```python
# reading protein sequence  using readFASTA
prt_seq = ss.readFASTA('./data/Oncorhynchus_mykiss.fasta')
```


```python
len(ss.head(prt_seq))
```




    129




```python
# basic info
ss.seqInfo(prt_seq)
```

    Sequence Length: 129 bp
    Frequency: {'E': 2, 'K': 5, 'G': 12, 'V': 9, 'Q': 5, 'A': 11, 'M': 1, 'F': 1, 'I': 5, 'P': 3, 'R': 10, 'D': 10, 'N': 9, 'W': 5, 'T': 7, 'S': 9, 'L': 10, 'H': 1, 'C': 8, 'Y': 6}
    Percentage: {'E': 0.02, 'K': 0.04, 'G': 0.09, 'V': 0.07, 'Q': 0.04, 'A': 0.09, 'M': 0.01, 'F': 0.01, 'I': 0.04, 'P': 0.02, 'R': 0.08, 'D': 0.08, 'N': 0.07, 'W': 0.04, 'T': 0.05, 'S': 0.07, 'L': 0.08, 'H': 0.01, 'C': 0.06, 'Y': 0.05}
    GC Content: 15.5
    AT Content: 13.95


## Protein Frequency


```python
ss.proteinFrequency(prt_seq)
```




    Counter({'K': 5,
             'V': 9,
             'Y': 6,
             'D': 10,
             'R': 10,
             'C': 8,
             'E': 2,
             'L': 10,
             'A': 11,
             'S': 9,
             'G': 12,
             'M': 1,
             'N': 9,
             'P': 3,
             'W': 5,
             'T': 7,
             'Q': 5,
             'I': 5,
             'F': 1,
             'H': 1})




```python
ss.plotProteinFrequency(prt_seq, kind='bar')
```


![png](output_50_0.png)


## Nucleotide Density


```python
ss.ntDensityOne(seq, windowsize=9000)
```


![png](output_52_0.png)



```python
ss.ntDensityTwo(seq, windowsize=9000)
```


![png](output_53_0.png)
