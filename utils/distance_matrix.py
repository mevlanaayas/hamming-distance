import copy


def hamming_distance(seq1, seq2):
    """ Calculate the Hamming distance between two bit strings """
    assert len(seq1) == len(seq2)
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def distance_matrix_func(cnt, seq, alpha):
    distance_matrix = {}
    for alpha1 in alpha:
        for alpha2 in alpha:
            if alpha.index(alpha1) < cnt and alpha.index(alpha2) < cnt:
                key = alpha1 + alpha2
                reverse_key = alpha2 + alpha1
                if alpha1 == alpha2:
                    distance_matrix[key] = 0
                elif alpha1 > alpha1:
                    distance_matrix[key] = distance_matrix[reverse_key]
                else:
                    hd = hamming_distance(seq[alpha.index(alpha1)], seq[alpha.index(alpha2)])
                    distance_matrix[key] = hd
            else:
                break
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
        reverse_result_key = item + step
        key = node_one + node_two
        key1 = node_one + item
        key2 = node_two + item
        result[result_key] = 1/2*(d_matrix[key1] + d_matrix[key2] - d_matrix[key])
        result[reverse_result_key] = 1/2*(d_matrix[key1] + d_matrix[key2] - d_matrix[key])
    result[step + step] = 0
    return result


def update_distance_matrix(cleaned_d_matrix, alpha, new_dist):
    new_distance_matrix = {}
    for alpha1 in alpha:
        for alpha2 in alpha:
            key = alpha1 + alpha2
            reverse_key = alpha2 + alpha1
            if key in cleaned_d_matrix.keys():
                if alpha1 == alpha2:
                    new_distance_matrix[key] = 0
                elif alpha1 > alpha1:
                    new_distance_matrix[key] = new_distance_matrix[reverse_key]
                else:
                    new_distance_matrix[key] = cleaned_d_matrix[key]
            elif key in new_dist.keys():
                new_distance_matrix[key] = new_dist[key]
            else:
                continue
    return new_distance_matrix

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
