from utils.sigma import sigma_operator


def calculate_branch_length_est(d_matrix, taxa_count, node_one, node_two, comp_tree, tree_step, alpha):
    key = node_one + node_two
    key1 = node_one + alpha[tree_step]
    part_one = 1/2 * d_matrix[key]
    part_two = 1/(2*(taxa_count-2))
    part_three = sigma_operator(k_from=0, k_to=taxa_count, d_matrix=d_matrix, d_ij=node_one, alph=alpha)
    part_four = sigma_operator(k_from=0, k_to=taxa_count, d_matrix=d_matrix, d_ij=node_two, alph=alpha)
    value = part_one + (part_two*(part_three - part_four))
    comp_tree[key1] = value
    key2 = node_two + alpha[tree_step]
    value = d_matrix[key] - comp_tree[key1]
    comp_tree[key2] = value
    return alpha[tree_step]
