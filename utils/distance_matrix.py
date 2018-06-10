import copy


def hamming_distance(seq1, seq2):
    """ Calculate the Hamming distance between two bit strings """
    assert len(seq1) == len(seq2)
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def distance_matrix_func(cnt, seq, alpha):
    distance_matrix = {}
    for i in range(cnt):
        for j in range(cnt):
            key = alpha[i] + alpha[j]
            reverse_key = alpha[j] + alpha[i]
            if i == j:
                distance_matrix[key] = 0
            elif i > j:
                distance_matrix[key] = distance_matrix[reverse_key]
            else:
                hd = hamming_distance(seq[i], seq[j])
                distance_matrix[key] = hd
    return distance_matrix


def clear_distance_matrix(d_matrix, node_one, node_two):
    temp = copy.deepcopy(d_matrix)
    remaining_nodes = []
    for key in d_matrix.keys():
        if node_one in key or node_two in key:
            del temp[key]
        elif key[0] == key[1]:
            remaining_nodes.append(key[0])
    return temp, remaining_nodes


def find_new_distances(d_matrix, nodes, node_one, node_two, step):
    result = {}
    for item in nodes:
        result_key = step + item
        key = node_one + node_two
        key1 = node_one + item
        key2 = node_two + item
        result[result_key] = 1/2*(d_matrix[key1] + d_matrix[key2] - d_matrix[key])
    return result


def update_distance_matrix(d_matrix, node_one, node_two, alpha, created_dist, rm_nodes, cnt):
    pass
    # result_matrix = {}
    # for item in rm_nodes:  # c, d, e
    #     if alpha.index(node_two) < alpha.index(item):
    #         for i in range(cnt):
    #             for j in range(cnt):
    #                 key = alpha[i] + alpha[j]
    #                 reverse_key = alpha[j] + alpha[i]
    #                 if i == j:
    #                     result_matrix[key] = 0
    #                 elif i > j:
    #                     result_matrix[key] = result_matrix[reverse_key]
    #                 else:
    #                     for key_item in created_dist.keys():
    #                         if item in key_item:
    #                             result_matrix[key] = created_dist[key_item]
    #     else:
    #         continue
    # return result_matrix
