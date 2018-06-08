from utils.sigma import sigma_operator


def q_matrix(d_matrix, taxa_count):
    q_matrix_result = [
        [0] * 5,
        [0] * 5,
        [0] * 5,
        [0] * 5,
        [0] * 5
    ]
    for a in range(taxa_count):
        for b in range(taxa_count):
            if a == b:
                print("diagonal")
                q_matrix_result[a][b] = 0
            elif a > b:
                q_matrix_result[a][b] = q_matrix_result[b][a]
            else:
                sum_one = sigma_operator(k_from=0, k_to=taxa_count, d_matrix=d_matrix, d_ij=a)
                sum_two = sigma_operator(k_from=0, k_to=taxa_count, d_matrix=d_matrix, d_ij=b)
                q_matrix_result[a][b] = (taxa_count - 2) * d_matrix[a][b] - sum_one - sum_two

    return q_matrix_result


def find_min(q_mat, cnt):
    min_of_q = 0
    a, b = 0, 0
    for i in range(cnt):
        for j in range(cnt):
            if i < j:
                if q_mat[i][j] < min_of_q:
                    min_of_q = q_mat[i][j]
                    a, b = i, j
    return a, b
