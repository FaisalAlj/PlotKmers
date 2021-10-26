import itertools
from Bio import SeqIO
from matplotlib import pyplot as plt
import pandas as pd
from collections import Counter

fasta_file = open("D:\SRR003265.filt.fastq")

# read and filter the Fastq file
# return Sequences only
def reads(file):
    records = SeqIO.parse(file, format="fastq")
    seq_list = []
    for record in records:
        seq = str(record.seq)
        if "N" in seq:
            continue
        if len(seq_list) >= 20000:
            break
        seq_list.append(seq)
    return seq_list


# Histogram
# All Possible k k-mers

def find_Kmers(k, seq_list):
    kmers_dict = {}
    for seq in seq_list:
        for i in range(0, ((len(seq) - k) + 1)):
            sequence = seq[i: i + k]
            if sequence not in kmers_dict:
                kmers_dict[sequence] = 1
            else:
                kmers_dict[sequence] += 1
    return kmers_dict

# Generate all possible neighbors (kmers) with hamming distance 1


def HammingDistance(kmer, alphabet=None):
    neighbors = []
    if alphabet is None:
        alphabet = {'A', 'C', 'G', 'T'}
    for lists in itertools.combinations(range(len(kmer) - 1, -1, -1), 1):  # lists [4],[3],[2],[1] lexicographical order
        for characters in itertools.product(alphabet, repeat=1):  # Characters [A,C,G,T]
            newkmer = list(kmer)
            for index in lists:
                for char in characters:
                    if char != newkmer[index]:
                        newkmer[index] = char
                        neighbors.append("".join(newkmer))
    return neighbors  

# Corrected reads .... Error Free

def newReads(kmersFrequence, seqs, t, k):
    new_Reads = []
    for seq in seqs:
        n_of_kmers = (len(seq) - k) + 1
        for i in range(n_of_kmers):

            kmer = seq[i:i + k]
            if kmer in kmersFrequence:
                if kmersFrequence[kmer] < t:
                    all_neighbors = HammingDistance(kmer)
                    for nk in all_neighbors:
                        if nk in kmersFrequence:
                            if kmersFrequence[nk] >= t:
                                seq = seq[:i] + nk + seq[i + k:]
                                break
        new_Reads.append(seq)
    return new_Reads


def plot_2_in_One(dict1, dict2):
    kmers_dict1 = sorted(
        Counter(dict1.values()).items())  # Sortiertes Dictionary ....  {"K-mer Count":"Distinct k-mer"}....
    kmers_dict2 = sorted(Counter(dict2.values()).items())
    df1 = pd.DataFrame(kmers_dict1, columns=["K-mer Count",
                                             "Distinct k-mer"])  # Data Frame mit 2 Spalten ("K-mer Count" , "Distinct k-mer")
    df2 = pd.DataFrame(kmers_dict2, columns=["K-mer Count", "Distinct k-mer"])
    plt.plot(df1["K-mer Count"], df1["Distinct k-mer"], 'b.-', label="Error")
    plt.plot(df2["K-mer Count"], df2["Distinct k-mer"], 'r.-', label="Error Free")
    plt.title("Kmers Distribution")
    plt.xlabel("K-mer Counts")
    plt.ylabel("# Distinct K-mers with that Count")
    plt.legend()
    plt.xlim([0, 50])
    plt.show()


seqs = reads(fasta_file)
Histogram = find_Kmers(8, seqs)

newRd = newReads(Histogram, seqs, 6, 8)
newHisto = find_Kmers(8, newRd)

plot_2_in_One(Histogram, newHisto)