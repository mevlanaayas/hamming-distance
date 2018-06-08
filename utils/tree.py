from utils.sigma import sigma_operator


def calculate_branch_length_est(d_matrix, taxa_count, node_one, node_two, comp_tree):
    key1 = 'd(' + str(node_one) + ', new_conn)'
    part_one = 1/2 * d_matrix[node_one][node_two]
    part_two = 1/(2*(taxa_count-2))
    part_three = sigma_operator(k_from=0, k_to=taxa_count, d_matrix=d_matrix, d_ij=node_one)
    part_four = sigma_operator(k_from=0, k_to=taxa_count, d_matrix=d_matrix, d_ij=node_two)
    value = part_one + (part_two*(part_three - part_four))
    comp_tree[key1] = value
    key2 = 'd(' + str(node_two) + ', new_conn)'
    value = d_matrix[node_one][node_two] - comp_tree[key1]
    comp_tree[key2] = value
