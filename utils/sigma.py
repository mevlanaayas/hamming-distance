
def sigma_two(node, alph, d_matrix):
    local_sum = 0
    for alp in alph:
        key = node + alp
        if key in d_matrix.keys():
            local_sum += d_matrix[key]
    return local_sum
