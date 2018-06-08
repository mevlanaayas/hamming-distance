

def sigma_operator(k_from, k_to, d_matrix, d_ij):
    local_sum = 0
    for index in range(k_from, k_to, 1):
        local_sum += d_matrix[d_ij][index]
    return local_sum
