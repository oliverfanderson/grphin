def get_three_node_graphlet_dist_adj_list_v4(G: nx.MultiDiGraph, G_prime: nx.Graph):
    three_node_graphlet_dict = {}
    graphlet_mapper = {}
    start_time = time.time()

    # create all the binary edge vectors
    adj_list_vector = [{} for _ in range(len(G.nodes()))]
    for i, j, data in G.edges(data=True):
        label = data.get("label")
        if label == "ppi":
            if j not in adj_list_vector[i]:
                adj_list_vector[i][j] = [1, 0, 0]
            # else:
            #     adj_list_vector[i][j][0] += 1
            if i not in adj_list_vector[j]:
                adj_list_vector[j][i] = [1, 0, 0]
            # else:
            #     adj_list_vector[j][i][0] += 1
        elif label == "reg":
            if j not in adj_list_vector[i]:
                adj_list_vector[i][j] = [0, 1, 0]
            else:
                adj_list_vector[i][j][1] += 1
            if i not in adj_list_vector[j]:
                adj_list_vector[j][i] = [0, 0, 1]
            else:
                adj_list_vector[j][i][2] += 1

    # find all combinations of potential 3 node graphlets
    # pick an edge between A and B
    # for each edge pair, find the union of neighbors between A and B
    three_node_combination = set()  # Use a set for fast triplet lookups

    # Preprocess neighbors for each node once
    neighbors_dict = {i: set(G_prime.neighbors(i)) for i in G_prime.nodes()}
    completed_i = set()

    for i in G_prime.nodes():
        print(f"Node: {i}", end="\r")
        for j in neighbors_dict[i]:
            for k in neighbors_dict[j].difference(completed_i):
                if (
                    (i < k) and (i != j) and (j != k)
                ):  # Ensure no duplicates by enforcing i < k and i != j
                    triplet = tuple(sorted([i, j, k]))
                    if triplet not in three_node_combination:
                        three_node_combination.add(triplet)
                        # print(f"Triplet: {i}, {j}, {k}")

                        a = i
                        b = j
                        c = k

                        ab = ac = ba = bc = ca = cb = 0
                        if b in adj_list_vector[a]:
                            ab = adj_list_vector[a][b]
                        else:
                            ab = [0, 0, 0]
                        if c in adj_list_vector[a]:
                            ac = adj_list_vector[a][c]
                        else:
                            ac = [0, 0, 0]
                        if a in adj_list_vector[b]:
                            ba = adj_list_vector[b][a]
                        else:
                            ba = [0, 0, 0]
                        if c in adj_list_vector[b]:
                            bc = adj_list_vector[b][c]
                        else:
                            bc = [0, 0, 0]
                        if a in adj_list_vector[c]:
                            ca = adj_list_vector[c][a]
                        else:
                            ca = [0, 0, 0]
                        if b in adj_list_vector[c]:
                            cb = adj_list_vector[c][b]
                        else:
                            cb = [0, 0, 0]
                        a_b, a_c, b_a, b_c, c_a, c_b = get_three_node_graphlet_dict(
                            ab, ac, ba, bc, ca, cb
                        )

                        # get unsorted edge values for A, B, C
                        a_unsort = [a_b, a_c]
                        b_unsort = [b_a, b_c]
                        c_unsort = [c_a, c_b]

                        # order A, B, C edge values internally
                        a_sort = tuple(sorted(a_unsort))
                        b_sort = tuple(sorted(b_unsort))
                        c_sort = tuple(sorted(c_unsort))

                        # Create a list of tuples in order [A, B, C]
                        tuples_list = [a_sort, b_sort, c_sort]

                        unsorted_tuples = tuple(tuples_list)

                        # Sort the tuples first by the first index, then by the second index
                        sorted_tuples = tuple(
                            sorted(tuples_list, key=lambda x: (x[0], x[1]))
                        )

                        # Add the graphlet if it has not been seen yet and update the count for the graphlet
                        if hash(sorted_tuples) not in three_node_graphlet_dict:
                            three_node_graphlet_dict[hash(sorted_tuples)] = 0
                            graphlet_mapper[hash(sorted_tuples)] = sorted_tuples
                        three_node_graphlet_dict[hash(sorted_tuples)] += 1

        # Once we're done processing i, mark it as completed
        completed_i.add(i)
        # print(f"Completed node: {i}")

    run_time = time.time() - start_time
    print("run time : %.3f seconds" % run_time)
    return three_node_graphlet_dict, graphlet_mapper


#    # def get_orbit_per_graphlet(orbit_dict, orbit_mapper, sorted_tuples, a_edges, b_edges, c_edges, i, j, k):

#     # # 3-node line basic only ppi example 1
#     # we can check which graphlet we are looking at via the sorted tuple
#     if sorted_tuples == ((0, 1), (0, 1), (1, 1)):
#         a_expected = (0, 1)
#         b_expected = (0, 1)
#         c_expected = (1, 1)
#         # since we arbitrarily chose what i, j, and k is,
#         # we need to check which of them match the a, b, c in the sorted tuple
#         orbit_change = get_orbit_position_change(
#             a_edges, b_edges, c_edges, a_expected, b_expected, c_expected, i, j, k
#         )

#         # store orbits with the key:
#         # graphlet_sorted_tuple + orbit_position
#         # i.e
#         # (('0', '1'), ('0', '1'), ('1', '1') ('0','0')) this specific graphlet at its 0th orbit
#         # (('0', '1'), ('0', '1'), ('1', '1') ('1','1')) this specific graphlet at its 1th orbit
#         zero_orbit = ("G1", 0)
#         first_orbit = ("G1", 1)

#         a_orbit = sorted_tuples + zero_orbit
#         b_orbit = sorted_tuples + first_orbit
#         # c_orbit = graphlet_sorted_tuple + second_orbit

#         if hash(a_orbit) not in orbit_dict:
#             orbit_dict[hash(a_orbit)] = []
#             orbit_mapper[hash(a_orbit)] = "Graphlet 1, Orbit 0"
#         orbit_dict[hash(a_orbit)] += [orbit_change[0], orbit_change[1]]
#         if hash(b_orbit) not in orbit_dict:
#             orbit_dict[hash(b_orbit)] = []
#             orbit_mapper[hash(b_orbit)] = "Graphlet 1, Orbit 1"
#         orbit_dict[hash(b_orbit)] += [orbit_change[2]]

#     elif sorted_tuples == ((0, 1), (0, 4), (1, 5)):
#         a_expected = (0, 1)
#         b_expected = (0, 4)
#         c_expected = (1, 5)
#         # since we arbitrarily chose what i, j, and k is,
#         # we need to check which of them match the a, b, c in the sorted tuple
#         orbit_change = get_orbit_position_change(
#             a_edges, b_edges, c_edges, a_expected, b_expected, c_expected, i, j, k
#         )

#         # store orbits with the key:
#         zero_orbit = ("G2", 0)
#         first_orbit = ("G2", 1)
#         second_orbit = ("G2", 2)

#         a_orbit = sorted_tuples + zero_orbit
#         b_orbit = sorted_tuples + first_orbit
#         c_orbit = sorted_tuples + second_orbit

#         if hash(a_orbit) not in orbit_dict:
#             orbit_dict[hash(a_orbit)] = []
#             orbit_mapper[hash(a_orbit)] = "Graphlet 2, Orbit 0"
#         orbit_dict[hash(a_orbit)] += [orbit_change[0]]
#         if hash(b_orbit) not in orbit_dict:
#             orbit_dict[hash(b_orbit)] = []
#             orbit_mapper[hash(b_orbit)] = "Graphlet 2, Orbit 1"
#         orbit_dict[hash(b_orbit)] += [orbit_change[1]]
#         if hash(c_orbit) not in orbit_dict:
#             orbit_dict[hash(c_orbit)] = []
#             orbit_mapper[hash(c_orbit)] = "Graphlet 2, Orbit 2"
#         orbit_dict[hash(c_orbit)] += [orbit_change[2]]

#     elif sorted_tuples == ((0, 1), (0, 5), (1, 4)):
#         a_expected = (0, 1)
#         b_expected = (0, 5)
#         c_expected = (1, 4)
#         # since we arbitrarily chose what i, j, and k is,
#         # we need to check which of them match the a, b, c in the sorted tuple
#         orbit_change = get_orbit_position_change(
#             a_edges, b_edges, c_edges, a_expected, b_expected, c_expected, i, j, k
#         )

#         # store orbits with the key:
#         zero_orbit = ("G3", 0)
#         first_orbit = ("G3", 1)
#         second_orbit = ("G3", 2)

#         a_orbit = sorted_tuples + zero_orbit
#         b_orbit = sorted_tuples + first_orbit
#         c_orbit = sorted_tuples + second_orbit

#         if hash(a_orbit) not in orbit_dict:
#             orbit_dict[hash(a_orbit)] = []
#             orbit_mapper[hash(a_orbit)] = "Graphlet 3, Orbit 0"
#         orbit_dict[hash(a_orbit)] += [orbit_change[0]]
#         if hash(b_orbit) not in orbit_dict:
#             orbit_dict[hash(b_orbit)] = []
#             orbit_mapper[hash(b_orbit)] = "Graphlet 3, Orbit 1"
#         orbit_dict[hash(b_orbit)] += [orbit_change[1]]
#         if hash(c_orbit) not in orbit_dict:
#             orbit_dict[hash(c_orbit)] = []
#             orbit_mapper[hash(c_orbit)] = "Graphlet 3, Orbit 2"
#         orbit_dict[hash(c_orbit)] += [orbit_change[2]]

#     elif sorted_tuples == ((0, 4), (0, 4), (5, 5)):
#         a_expected = (0, 4)
#         b_expected = (0, 4)
#         c_expected = (5, 5)
#         # since we arbitrarily chose what i, j, and k is,
#         # we need to check which of them match the a, b, c in the sorted tuple
#         orbit_change = get_orbit_position_change(
#             a_edges, b_edges, c_edges, a_expected, b_expected, c_expected, i, j, k
#         )

#         # store orbits with the key:
#         zero_orbit = ("G4", 0)
#         first_orbit = ("G4", 1)
#         # second_orbit = ('G3', 2)

#         a_orbit = sorted_tuples + zero_orbit
#         b_orbit = sorted_tuples + first_orbit
#         # c_orbit = sorted_tuples + second_orbit

#         if hash(a_orbit) not in orbit_dict:
#             orbit_dict[hash(a_orbit)] = []
#             orbit_mapper[hash(a_orbit)] = "Graphlet 4, Orbit 0"
#         orbit_dict[hash(a_orbit)] += [orbit_change[0], orbit_change[1]]
#         if hash(b_orbit) not in orbit_dict:
#             orbit_dict[hash(b_orbit)] = []
#             orbit_mapper[hash(b_orbit)] = "Graphlet 4, Orbit 1"
#         orbit_dict[hash(b_orbit)] += [orbit_change[2]]

#     return orbit_dict




def get_adjacency_matrix(G):
    """Get the adjacency matrix for a MultiDiGraph"""

    G_adj_matrix = nx.adjacency_matrix(G)
    adj_matrix = [
        [[0, 0, 0] for _ in range(len(G.nodes()))] for _ in range(len(G.nodes()))
    ]

    for edge in G.edges(data=True):
        i = edge[0]
        j = edge[1]
        edge_type = edge[2]["label"]

        if edge_type == "ppi":
            adj_matrix[i][j][0] = 1
            adj_matrix[j][i][0] = 1
        elif edge_type == "reg":
            adj_matrix[i][j][1] = 1
            adj_matrix[j][i][2] = 1

    return G_adj_matrix



def get_two_node_graphlet_dist_adj_matrix(G, two_node_graphlet_dict):
    G_adj_matrix = get_adjacency_matrix(G)
    for i in range(len(G_adj_matrix)):
        for j in range(len(G_adj_matrix[0])):
            vector = G_adj_matrix[i][j] + G_adj_matrix[j][i]
            if hash(tuple(vector)) in two_node_graphlet_dict:
                two_node_graphlet_dict[hash(tuple(vector))] += 1
    return two_node_graphlet_dict
