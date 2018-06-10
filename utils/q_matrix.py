from utils.sigma import sigma_operator


def q_matrix(d_matrix, taxa_count, alpha):
    q_matrix_result = {}
    for i in range(taxa_count):
        for j in range(taxa_count):
            key = alpha[i] + alpha[j]
            reverse_key = alpha[j] + alpha[i]
            if i == j:
                print("diagonal")
                q_matrix_result[key] = 0
            elif i > j:
                q_matrix_result[key] = q_matrix_result[reverse_key]
            else:
                sum_one = sigma_operator(k_from=0, k_to=taxa_count, d_matrix=d_matrix, d_ij=alpha[i], alph=alpha)
                sum_two = sigma_operator(k_from=0, k_to=taxa_count, d_matrix=d_matrix, d_ij=alpha[j], alph=alpha)
                q_matrix_result[key] = (taxa_count - 2) * d_matrix[key] - sum_one - sum_two

    return q_matrix_result


def find_min(q_mat, cnt, alpha):
    min_of_q = 0
    nodes_to_join = []
    for i in range(cnt):
        for j in range(cnt):
            if i < j:
                key = alpha[i] + alpha[j]
                if q_mat[key] < min_of_q:
                    min_of_q = q_mat[key]
                    nodes_to_join = alpha[i], alpha[j]
    return nodes_to_join
