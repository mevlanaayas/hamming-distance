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
from utils.distance_matrix import distance_matrix_func
from utils.q_matrix import q_matrix, find_min
from utils.tree import calculate_branch_length_est

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
    complete_tree = {}
    for seq_record in seq_input:
        sequences.append(seq_record)
        count += 1
    # dna sequence elimizde (array olarak)
    distance_matrix = distance_matrix_func(cnt=count, alpha=alphabet, seq=sequences)
    # distance matrix elimizde

    temp_distance_matrix = [
        [0, 5, 9, 9, 8],
        [5, 0, 10, 10, 9],
        [9, 10, 0, 8, 7],
        [9, 10, 8, 0, 3],
        [8, 9, 7, 3, 0]
    ]
    temp_count = 5

    q = q_matrix(d_matrix=temp_distance_matrix, taxa_count=len(temp_distance_matrix))
    # q matrix elimizde
    index_of_a, index_of_b = find_min(q_mat=q, cnt=temp_count)
    # minimum değerlerin indexleri elimizde
    calculate_branch_length_est(d_matrix=temp_distance_matrix,
                                taxa_count=temp_count,
                                node_one=index_of_a,
                                node_two=index_of_b,
                                comp_tree=complete_tree)
    # calculated new node and distance to connection

    c = 4
