import copy

from utils.printer import print_header, print_result


def hamming_distance(seq1, seq2):
    """ Calculate the Hamming distance between two bit strings """
    assert len(seq1) == len(seq2)
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def distance_matrix_func(cnt, alpha, seq):
    result_j = []
    distance_matrix = []
    print_header(cnt, alpha)
    for i in range(cnt):
        for j in range(cnt):
            hd = hamming_distance(seq[i], seq[j])
            result_j.append(hd)
        print_result(result_j, i, alpha)
        temp = copy.deepcopy(result_j)
        distance_matrix.append(temp)
        result_j = []
    return distance_matrix


# def update_distance_matrix(d_matrix, node_one, node_two, cnt):
#     for i in range(cnt):
#         for j in range(cnt):
#             if i == j:
#                 print("diagonal")
#                 d_matrix[i][j] = 0
#             elif i > j:
#                 d_matrix[i][j] = d_matrix[j][i]
#             else:
#                 d_matrix[i][j] = 1/2*(d_matrix[node_one][node_two] + d_matrix[][] - d_matrix[][])
#
#
