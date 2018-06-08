def print_result(distance_data, current_index, alpha):
    """ print result to txt file """
    file = open("dm.txt", "a")
    file.write(alpha[current_index])
    for distance in distance_data:
        file.write(", " + str(distance))
    file.write("\n")
    file.close()


def print_header(cnt, alpha):
    """ print header (A, B, C..) to file """
    file = open("dm.txt", "a")
    for index in range(cnt):
        file.write(", " + alpha[index])
    file.write("\n")
    file.close()