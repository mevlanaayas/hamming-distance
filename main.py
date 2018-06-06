#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author: mevlana ayas - mevlanaayas@gmail.com
"""
reference for nj algorithm
http://www.wiki-zero.net/index.php?q=aHR0cHM6Ly9lbi53aWtpcGVkaWEub3JnL3dpa2kvTmVpZ2hib3Jfam9pbmluZw

##### NJ Algorithm from wikipedia
Neighbor joining takes as input a distance matrix specifying the distance between each pair of taxa.
The algorithm starts with a completely unresolved tree, whose topology corresponds to that of a star network,
and iterates over the following steps until the tree is completely resolved and all branch lengths are known:

1-) Based on the current distance matrix calculate the matrix Q (defined below).
2-) Find the pair of distinct taxa i and j (i.e. with i != j) for which
    Q(i,j) has its lowest value. These taxa are joined to a newly created node,
    which is connected to the central node. In the figure at right, f and g are joined to the new node u.
3-) Calculate the distance from each of the taxa in the pair to this new node.
4-) Calculate the distance from each of the taxa outside of this pair to the new node.
5-) Start the algorithm again, replacing the pair of joined neighbors with the new node and using the distances
    calculated in the previous step.


##### The Q-matrix
Based on a distance matrix relating the n taxa, calculate Q as follows:
Q(i,j)=(n-2)d(i,j)-sum[k=1 to n](d(i,k))-sum[k=1 to n](d(j,k))
where d(i,j) is the distance between taxa i and j.


#### taxa means source of DNA sequence

#### reference for hamming distance
http://www.wiki-zero.net/index.php?q=aHR0cHM6Ly9lbi53aWtpcGVkaWEub3JnL3dpa2kvSGFtbWluZ19kaXN0YW5jZQ

def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))

The Hamming distance between:
    "karolin" and "kathrin" is 3.
    "karolin" and "kerstin" is 3.
     1011101  and  1001001  is 2.
     2173896  and  2233796  is 3.
"""

import string
from Bio import SeqIO


def hamming_distance(seq1, seq2):
    """ Calculate the Hamming distance between two bit strings """
    assert len(seq1) == len(seq2)
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def print_result(distance_data, current_index, alpha):
    """ print result to txt file """
    file = open("dm.txt", "a")
    file.write(alpha[current_index])
    for distance in distance_data:
        file.write(", " + str(distance))
    file.write("\n")
    file.close()


def print_header(cnt, alpha):
    """ print header (A, B, C..) to file """
    file = open("dm.txt", "a")
    for index in range(cnt):
        file.write(", " + alpha[index])
    file.write("\n")
    file.close()


if __name__ == '__main__':
    filename = input("please give input file with extension: ")
    extension = filename.split(".")[-1]
    if extension == 'gb':
        extension = "genbank"
    elif extension == 'fasta':
        pass
    else:
        exit(0)
    seq_input = SeqIO.parse(filename, extension)
    count = 0
    alphabet = list(string.ascii_uppercase)
    sequences = []
    result_j = []
    for seq_record in seq_input:
        sequences.append(seq_record)
        count += 1
    print_header(count, alphabet)
    for i in range(count):
        for j in range(count):
            hd = hamming_distance(sequences[i], sequences[j])
            result_j.append(hd)
        print_result(result_j, i, alphabet)
        result_j = []
