

def sigma_operator(k_from, k_to, d_matrix, d_ij, alph):
    local_sum = 0
    for index in range(k_from, k_to, 1):
        key = d_ij + alph[index]
        local_sum += d_matrix[key]
    return local_sum


def sigma_two(node, alph, d_matrix):
    local_sum = 0
    for alp in alph:
        key = node + alp
        if key in d_matrix.keys():
            local_sum += d_matrix[key]
    return local_sum
