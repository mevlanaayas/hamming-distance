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
from math import sqrt

from Bio import SeqIO
from utils.distance_matrix import distance_matrix_func, update_distance_matrix, clear_distance_matrix, \
    find_new_distances
from utils.q_matrix import q_matrix, find_min
from utils.tree import calculate_branch_length_est, build_tree

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
    distance_matrix_length = 0
    alphabet = list(string.ascii_lowercase)
    sequences = []
    distance_from_joined_node_to_joint = {}
    for seq_record in seq_input:
        sequences.append(seq_record)
        distance_matrix_length += 1
    # dna sequence elimizde (array olarak)
    distance_matrix = distance_matrix_func(cnt=distance_matrix_length, seq=sequences, alpha=alphabet)
    # distance matrix elimizde

    temp_distance_dict = {
        'aa': 0, 'ab': 5, 'ac': 9, 'ad': 9, 'ae': 8,
        'ba': 5, 'bb': 0, 'bc': 10, 'bd': 10, 'be': 9,
        'ca': 9, 'cb': 10, 'cc': 0, 'cd': 8, 'ce': 7,
        'da': 9, 'db': 10, 'dc': 8, 'dd': 0, 'de': 3,
        'ea': 8, 'eb': 9, 'ec': 7, 'ed': 3, 'ee': 0
    }
    step = -1
    while int(sqrt(len(temp_distance_dict))) > 2:
        temp_count = int(sqrt(len(temp_distance_dict)))
        # while distance_matrix has nodes to join
        # do

        q = q_matrix(d_matrix=temp_distance_dict, taxa_count=temp_count, alpha=alphabet)
        # q matrix elimizde

        first_node_to_join, second_node_to_join = find_min(q_mat=q)
        # minimum değerlerin indexleri elimizde

        backward_step = calculate_branch_length_est(d_matrix=temp_distance_dict,
                                                    taxa_count=temp_count,
                                                    node_one=first_node_to_join,
                                                    node_two=second_node_to_join,
                                                    comp_tree=distance_from_joined_node_to_joint,
                                                    tree_step=step,
                                                    alpha=alphabet
                                                    )
        step -= 1
        # calculated new node and distance to connection

        clean_d_m, remaining_nodes = clear_distance_matrix(temp_distance_dict,
                                                           node_one=first_node_to_join,
                                                           node_two=second_node_to_join
                                                           )
        # delete items if key include connected nodes

        new_distances = find_new_distances(d_matrix=temp_distance_dict,
                                           nodes=remaining_nodes,
                                           node_one=first_node_to_join,
                                           node_two=second_node_to_join,
                                           step=backward_step
                                           )
        # distance of remaining nodes to connected nodes

        temp_distance_dict = update_distance_matrix(cleaned_d_matrix=clean_d_m,
                                                    alpha=alphabet,
                                                    new_dist=new_distances
                                                    )
        # update distance matrix. and repeat q matrix so on...

    # result in distance_from_joined_node_to_joint
    # and we need to add last connection that is stored in new_distances

    build_tree(complete_tree=distance_from_joined_node_to_joint, remaining_distances=new_distances)
    # DONE result in distance_from_joined_node_to_joint now.
