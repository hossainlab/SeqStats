# Imports 
from collections import Counter, deque 
import matplotlib.pyplot as plt 
plt.style.use('seaborn-whitegrid')
plt.figure(figsize=(10,6))

#####################################################################################
# File Handling Functions  
#####################################################################################

# Reading Text File 
def readText(inputfile): 
    """Reads a sequence file in text format and returns the file with special character removed.
    Args:
        inputfile ([Text]):
    Returns:
        [type]: [String]
    """
    with open(inputfile, "r") as seqfile: 
        # read data 
        seq = seqfile.read()
        # remove special characters \n and \t 
        seq = seq.replace("\n", "")
        seq = seq.replace("\t", "") 
    return seq

# Reading FASTA File 
def readFASTA(inputfile): 
    """Reads a FASTA file and skips the sequence information line with special character removed"
    Args:
        inputfile ([type]): [FASTA File]
    Returns:
        [type]: [String]
    """
    with open(inputfile, "r") as f: 
        # remove name line / info line 
        seq = f.readline()
        # read data 
        seq = f.read() 
        # remove special character 
        seq = seq.replace("\n", "")
        seq = seq.replace("\t", "") 
    return seq 

# Reading FASTQ File 
def readFASTQ(inputfile):
    """Reads FASTQ file and returns with special character removed"
    Args:
        filename ([type]): [FASTQ]
    Returns:
        [type]: [String]
    """
    sequences = []
    qualities = []
    with open(inputfile) as fh:
        while True:
            fh.readline() # skip name line
            seq = fh.readline().rstrip() # read base sequence
            fh.readline() # skip placeholder line
            qual = fh.readline().rstrip() #base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities 


#####################################################################################
# Head and Tail of Sequence  
#####################################################################################

# Examine first few base pairs 
def head(seq):
    """Show first 1000bp of a sequence. 
    Args:
        seq ([type]): [Text File]
    Returns:
        [type]: [String]
    """
    h = seq[0:1000] 
    return h  
# Examine last few base pairs 
def tail(seq):
    """Show last 1000bp of a sequence. 
    Args:
        seq ([type]): [Text File]
    Returns:
        [type]: [String]
    """
    t = seq[-1000:] 
    return t 

#####################################################################################
# Frequency and Percentage with Plotting 
#####################################################################################

def frequency(seq, useall=False, calpc=False): 
    """Count the frequencies of each bases in sequence including every letter"""
    length = len(seq) 
    if calpc: 
        # Make a dictionary "freqs" to contain the frequency(in % ) of each base. 
        freqs = {} 
    else: 
        # Make a dictionary "base_counts" to contain the frequency(whole number) of each base. 
        base_counts = {} 
        
    if useall: 
        # If we want to look at every letter that appears in the sequence. 
        seqset = set(seq) 
    else: 
        # If we just want to look at the four bases A, T, C, G 
        seqset = ("A", "T", "G", "C")

    for letter in seqset: 
        num = seq.count(letter)
        if calpc: 
            # The frequency is calculated out of the total sequence length, even though some bases are not A, T, C, G
            freq = round(num/length, 2) 
            freqs[letter] = freq 
        else: 
            # Contain the actual number of bases. 
            base_counts[letter] = num 
    if calpc: 
        return freqs  
    else: 
        return base_counts 

#####################################################################################
# SeqInfo 
#####################################################################################
def seqInfo(seq): 
    L = len(seq) 
    print(f"Sequence Length: {L} bp") 
    F = frequency(seq, useall=True)
    print(f"Frequency: {F}") 
    P = frequency(seq, useall=True,calpc=True)
    print(f"Percentage: {P}") 
    GC = calculateGC(seq)
    print(f"GC Content: {GC}")
    AT = calculateAT(seq) 
    print(f"AT Content: {AT}")

#####################################################################################
# Plotting 
#####################################################################################
def plotFrequency(seq, title="Nucleotide Frquency Distribution", xlab="Bases",ylab="Frequency", kind=None):
    """Makes a scatterplot of y-variable(s) against an x-variable"""
    if kind == 'bar': 
        freq = frequency(seq) 
        plt.figure(figsize=(10,6))
        plt.bar(freq.keys(), freq.values()) 
        plt.title(title,  fontsize=14) 
        plt.xlabel(xlab,  fontsize=14)
        plt.ylabel(ylab,  fontsize=14)
        plt.tight_layout() 
        plt.show() 
    elif kind == 'pie': 
        plt.figure(figsize=(10,6))
        freq = frequency(seq) 
        plt.pie(freq.values(), labels=freq.keys(), autopct='%1.1f%%', shadow=True)
        plt.tight_layout() 
        plt.show() 
    else: 
        print("Please select your visualization type either bar or pie chart!")


#####################################################################################
# GC Content and AT Content, GC Plotting, AT Plotting 
#####################################################################################
def calculateGC(seq): 
    """
    Take DNA sequence as input and calculate the GC content.
    """
    no_of_g = seq.count("G")
    no_of_c = seq.count("C")
    total = no_of_g + no_of_c 
    gc = round(total/len(seq) * 100, 2) 
    return gc 
# 7. Function for Calculating Sub-Sequence GC Content 
def subSeqGC(seq, window=1000):
    res = [0] *  100 
    for i in range(0, len(seq)-window+1, window):
        subseq = seq[i:i+window] 
        gc = calculateGC(subseq) 
        res.append(gc) 
    return res  

# 8. Plot GC Content 
def plotSubSeqGC(seq, window=1000, title="GC Content of Subsequence", xlab="Base-Pair Position(bp)", ylab="% GC"): 
    gc = subSeqGC(seq) 
    gc_kb = [x/10000 for x in gc]
    plt.figure(figsize=(10,6)) 
    plt.plot(range(len(gc_kb)), gc_kb)
    plt.title(title, fontsize=14) 
    plt.xlabel(xlab, fontsize=14) 
    plt.ylabel(ylab, fontsize=14) 
    plt.tight_layout() 
    plt.show()  

# 9. Function for Calculating AT Content 
def calculateAT(seq): 
    """
    Take DNA sequence as input and calculate the AT content.
    """
    no_of_a = seq.count("A")
    no_of_t = seq.count("T")
    total = no_of_a + no_of_t
    at = round(total/len(seq) * 100, 2) 
    return at 

def subSeqAT(seq, window=1000):
    res = [] 
    for i in range(0, len(seq)-window+1, window): 
        subseq = seq[i:i+window] 
        at = calculateAT(subseq) 
        res.append(at) 
    return res  

def plotSubSeqAT(seq, window=1000, title="AT Content of Subsequence", xlab="Base-Pair Position(bp)", ylab="% AT"): 
    at = subSeqGC(seq) 
    at_kb = [x/10000 for x in at]
    plt.figure(figsize=(10,6)) 
    plt.plot(range(len(at_kb)), at_kb)
    plt.title(title, fontsize=14) 
    plt.xlabel(xlab, fontsize=14) 
    plt.ylabel(ylab, fontsize=14) 
    plt.tight_layout() 
    plt.show() 

#####################################################################################
# K Mer Analysis, Plotting 
#####################################################################################
def buildKmers(sequence, ksize):
    """Take a sequence as input and build k-mers 
    Args:
        sequence ([type]): [String]
        ksize ([type]): [Integer]
    Returns:
        [type]: [List]
    """
    kmers = []
    n_kmers = len(sequence) - ksize + 1
    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        kmers.append(kmer)
    return kmers

# Counts K-mers 
from collections import Counter 
def countKmers(seq):
    """Counting frequencies of K-mers 
    Args:
        seq ([type]): [String]
    """
    counts = Counter(seq) 
    return counts 

# Plot K-mer Frequency 
def plotKmerFrequency(seq, title=False, xlab=False, ylab=False, kind=None): 
    counts = countKmers(seq) 
    if kind == 'bar': 
        plt.figure(figsize=(10,6))
        plt.bar(counts.keys(), counts.values())  
        plt.title(title, fontsize=14) 
        plt.xlabel(xlab, fontsize=14)
        plt.ylabel(ylab, fontsize=14) 
        plt.tight_layout() 
        plt.show() 
    elif kind == 'pie': 
        plt.figure(figsize=(10,6))
        plt.pie(counts.values(), labels=counts.keys(), autopct='%1.1f%%', shadow=True)
        plt.tight_layout() 
        plt.show()
    else: 
        print("Please select your visualization type either bar or pie chart!") 



#####################################################################################
# Protein Analysis 
#####################################################################################
def proteinFrequency(seq): 
    prt_freq = Counter(seq) 
    return prt_freq 

def plotProteinFrequency(seq, title="Protein Frquency Distribution", xlab="Proteins",ylab="Frequency", kind=None):
    """Makes a scatterplot of y-variable(s) against an x-variable"""
    if kind == 'bar': 
        freq = proteinFrequency(seq) 
        plt.figure(figsize=(10,6))
        plt.bar(freq.keys(), freq.values()) 
        plt.title(title,  fontsize=14) 
        plt.xlabel(xlab,  fontsize=14)
        plt.ylabel(ylab,  fontsize=14)
        plt.tight_layout() 
        plt.show() 
    elif kind == 'pie': 
        plt.figure(figsize=(10,6))
        freq = frequency(seq) 
        plt.pie(freq.values(), labels=freq.keys(), autopct='%1.1f%%', shadow=True)
        plt.tight_layout() 
        plt.show() 
    else: 
        print("Please select your visualization type either bar or pie chart!")


#####################################################################################
# Nucleotide Density 
#####################################################################################

def myscatterplot(xvar,ydict,xlab,ylab, title):
    """Makes a scatterplot of y-variable(s) against an x-variable"""
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(111)
    yvarnames = []
    for yvar in ydict:
        yvarnames.append(yvar)
        ax.plot(xvar,ydict[yvar])
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    ax.legend(yvarnames, loc="upper right")
    plt.show()



def ntDensityOne(seq,windowsize,verbose=False,jumpsize=1000,makePlot=True):
    """Plots the base frequency along a sequence using a sliding window
    """
    length = len(seq)
    # make a dictionary to contain four empty lists
    freqs = { "A": [], "T": [], "G": [], "C": [] }
    myset = ("A", "C", "T", "G")
    midpoints = []
    # Move the window by jumpsize bp at a time.
    for i in range(0,length-windowsize+1,jumpsize):
        subseq = seq[i:i+windowsize]
        if verbose:
            start = i
            end = i+windowsize
            print("start %d end %d subseq is %s length %d windowsize %d" % (start,end,subseq,length,windowsize))
        assert len(subseq)==windowsize, "ERROR: ntdensity2: length of subseq is not windowsize"
        for letter in myset:
            num = subseq.count(letter)
            pc = 100 * num/windowsize
            freqs[letter].append(pc)
        # Find the mid-point of the window:
        # For example, if the window is from i=1000 to i=11000,
        # midpoint = 12000/2 = 6000
        midpoint = (i + i + windowsize)/2
        midpoints.append(midpoint)
    if makePlot:
        # Call the plotting function
        midpoints2 = [x/1000 for x in midpoints] # list comprehension
        myscatterplot(midpoints2,freqs,'Base-Pair Position (kb)','% of Nucleotide', "Nucleotide Density") # Convert to kb for plotting
    # Return the results:
    # return(midpoints,freqs)
   

def ntDensityTwo(seq,windowsize,verbose=False,jumpsize=1000,makePlot=True):
    """Plots the G+C content along a sequence using a sliding window
    """
    length = len(seq)
    # Make a dictionary to contain two empty lists
    freqs = { "G+C": [], "A+T": [] }
    myset = ("A+T", "G+C")
    # Note: instead of storing G+C, A+T in a hash, could just have coded G+C=0, A+T=1, and stored these values in arrays.
    midpoints = []
    # Move the window by jumpsize bp at a time.
    # The first window we look at is from i=0 to i=windowsize.
    # For example, if the sequence is 30000 bases long, windowsize=10000.
    # In the first loop, we look at i=0 to i=10000.
    # In the second loop, we look at i=1000 to i=11000. ...
    # In the last loop, we look at i=29000 to i=30000.
    # Note: for i = range(0,10) goes from i=0...9. 
    for i in range(0,length-windowsize+1,jumpsize):
        subseq = seq[i:i+windowsize]
        if verbose:
            start = i
            end = i+windowsize
            print("start %d end %d subseq is %s length %d windowsize %d" % (start,end,subseq,length,windowsize))
        assert len(subseq)==windowsize, "ERROR: ntdensity1: length of subseq is not windowsize"
        for dimer in myset:
            letter1 = dimer[0:1]
            letter2 = dimer[2:3]
            num1 = subseq.count(letter1)
            num2 = subseq.count(letter2)
            num = num1 + num2
            pc = (100 * num)/windowsize
            freqs[dimer].append(pc)
        # Find the mid-point of the window:
        # For example, if the window is from i=1000 to i=11000,
        # midpoint = 12000/2 = 6000
        midpoint = (i + i + windowsize)/2
        midpoints.append(midpoint)
    if makePlot:
        # Call the plotting function
        midpoints2 = [x/1000 for x in midpoints] # list comprehension
        myscatterplot(midpoints2,freqs,'Base-Pair Position (kb)','%  of Nucleotide',"Di-Nucleotide Density") # Convert to kb for plotting
    # Return the results
    # return(midpoints,freqs)