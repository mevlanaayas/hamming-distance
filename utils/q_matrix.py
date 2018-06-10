from utils.sigma import sigma_two


def q_matrix(d_matrix, taxa_count, alpha):
    q_matrix_result = {}
    for alpha1 in alpha:
        for alpha2 in alpha:
            key = alpha1 + alpha2
            reverse_key = alpha2 + alpha1
            if key in d_matrix.keys():
                if alpha1 == alpha2:
                    q_matrix_result[key] = 0
                elif alpha1 > alpha2:
                    q_matrix_result[key] = q_matrix_result[reverse_key]
                else:
                    sum_one = sigma_two(d_matrix=d_matrix, node=alpha1, alph=alpha)
                    sum_two = sigma_two(d_matrix=d_matrix, node=alpha2, alph=alpha)
                    q_matrix_result[key] = (taxa_count - 2) * d_matrix[key] - sum_one - sum_two
            else:
                continue
    return q_matrix_result


def find_min(q_mat, cnt, alpha):
    min_of_q = 0
    nodes_to_join = []
    for alpha1 in alpha:
        for alpha2 in alpha:
            if alpha.index(alpha1) < cnt and alpha.index(alpha2) < cnt:
                if alpha1 < alpha2:
                    key = alpha1 + alpha2
                    if q_mat[key] < min_of_q:
                        min_of_q = q_mat[key]
                        nodes_to_join = alpha1, alpha2
    return nodes_to_join
