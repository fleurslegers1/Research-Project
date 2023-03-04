import matplotlib.pyplot as plt
import networkx as nx
import pandas
import numpy as np
import warnings
import matplotlib.cm as cmx
import math
import matplotlib.patches as mpatches

warnings.filterwarnings("ignore")

#list_of_col = [1, 2, 3, 4, 5, 6, 7, 8, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 292, 293, 294, 295]
list_of_col = [34]

for n_col in list_of_col:

    #### Load Population Sizes ####
    global sizes
    sizes = pandas.read_csv("inwoners per rwzi.csv",
                            sep = ";",
                            header = None,
                            names = ["RWZI", "size"],
                            decimal = ",")

    #### Load STP coordinates ####
    global coordinates
    coordinates = pandas.read_csv("RWZI Coordinates.csv",
                                  sep = ";",
                                  header = 0,
                                  decimal = ",")

    coordinates["RWZI"] = coordinates["RWZI"].str.upper()
    coordinates["RWZI"] = coordinates["RWZI"].str.replace(" ", "")
    coordinates["RWZI"] = coordinates["RWZI"].str.replace("-", "")

    #### Load coefficient matrix ####
    global W

    if(n_col == 34):
        W = pandas.read_csv("W for PCA_ns, p = 5, MA = 3, n_col = 34.csv",
                                    sep = " ",
                                    header = 0,
                                    decimal = ",")
    else:
        W = pandas.read_csv("W for PCA_ns, p = 5, MA = 3, n_col = "+ str(n_col) +".csv",
                                    sep = ",",
                                    header = 0,
                                    decimal = ",")

    p = 5

    # Make sure all row sums are 1
    W = W.astype(float)
    row_sums = (np.abs(W)).sum(axis = 1)
    for row in range(0,311):
        row_sum = row_sums.iloc[row]
        for col in range(0,1556):
            W.iloc[row, col] = W.iloc[row, col]/row_sum

    # For each row, for each STP choose the lag with highest coefficient to add to the edge list
    edges = []
    p_of_edges = []
    RWZIs = W.columns[0:311]

    # Remove .L1 from RWZIs
    RWZIs_2 = []
    for i in range(0,311):
        rwzi = RWZIs[i][:-3]
        rwzi = rwzi.replace("-", "")
        rwzi = rwzi.replace(" ", "")
        RWZIs_2.append(rwzi)

    # Find largest coefficient that is not "constant"
    for row in range(0,311):
        row_name = RWZIs_2[row]
        for STP in range(0, 311):
            col_name = RWZIs_2[STP]
            largest_coeff = 0
            p_of_largest_coeff = None
            for i in range(0, p):
                if np.abs(W.iloc[row][STP + (i * 311)]) > largest_coeff:
                    largest_coeff = np.abs(W.iloc[row][STP + (i * 311)])
                    p_of_largest_coeff = i + 1
            if(row_name != col_name and largest_coeff > 0):
                tuple_ = (row_name, col_name, largest_coeff)
                edges.append(tuple_)
                p_of_edges.append(p_of_largest_coeff)

    # Delete the nodes that are not in correlation
    coordinates = coordinates[coordinates["RWZI"].isin(RWZIs_2)]

    # Make dictionary with coordinates
    coordinates["pos"] = coordinates[["Latitude", "Longitude"]].apply(tuple, axis = 1)
    coordinates = dict(zip(coordinates.RWZI, coordinates.pos))

    myKeys = list(coordinates.keys())
    myKeys.sort()
    coordinates = {i: coordinates[i] for i in myKeys}

    # Make node-size array
    node_sizes = []

    nodes = list(coordinates.keys())
    for RWZI in nodes:
        if(RWZI not in sizes.RWZI.values):
            node_sizes.append(75)
            print("No size known for " + RWZI)
        else:
            size = sizes[sizes.RWZI == RWZI]["size"].iloc[0]
            node_sizes.append(math.log(6*size,2))

    # Make Graph
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_weighted_edges_from(edges)

    i = 0
    for u,v,d in G.edges(data = True):
        d['color'] = p_of_edges[i]
        i = i + 1

    edges, colors = zip(*nx.get_edge_attributes(G, 'color').items())
    edges, weights = zip(*nx.get_edge_attributes(G, 'weight').items())

    plt.figure(1, figsize=(7,8.5))

    # Calculate node sizes
    node_sizes  = {}
    node_colors = {}

    col_sums = (np.abs(W)).sum(axis=0)
    max_col_sum = 0

    # Give nodes color and size based on whether or not the node
    # was used as variable in the model
    for i in range(len(RWZIs_2)):
        col_sum = 0
        for j in range(0, 4):
            col_sum += col_sums.values[i + j * 311]
        if col_sum > max_col_sum:
            max_col_sum = col_sum
        if (col_sum == 0):
            value_ = 0.0
            node_colors[RWZIs_2[i]] = 'tab:gray'
        else:
            value_ = col_sum
            node_colors[RWZIs_2[i]] = 'black'
        node_sizes[RWZIs_2[i]] = value_

    # Change range of node sizes to lie within (1,20)
    for key, value in node_sizes.items():
        node_sizes[key] = int(1.0 + (29.0 * node_sizes[key]) / max_col_sum)

    # Sort node_sizes and node_colors so that they appear in same order as nodes
    nx.draw_networkx_nodes(G, pos = coordinates,
            nodelist = list(node_sizes.keys()),
            node_size = np.array(list(node_sizes.values())),
            node_color= np.array(list(node_colors.values())))

    sorted_weights = list(np.abs(weights))
    sorted_weights = sorted(sorted_weights, key=int)
    index_ = int(np.ceil(len(sorted_weights) * 0.75))

    min_ = sorted_weights[-index_] # should be min of top 25%
    print(min_)

    max_ = max(np.abs(weights))
    print(max_)

    # For each lag, use a different colour map
    cmap_1 = cmx.get_cmap('Blues')
    cmap_2 = cmx.get_cmap('Oranges')
    cmap_3 = cmx.get_cmap('Greens')
    cmap_4 = cmx.get_cmap('Reds')
    cmap_5 = cmx.get_cmap('RdPu')

    # for each edge, add color of edge based on lag and weight
    rgb_colors = []
    for i in range(len(colors)):
        if colors[i] == 1:
            rgb_color = cmap_1((abs(weights[i]) - min_)/(max_ - min_))
        elif colors[i] == 2:
            rgb_color = cmap_2((abs(weights[i]) - min_)/(max_ - min_))
        elif colors[i] == 3:
            rgb_color = cmap_3((abs(weights[i]) - min_)/(max_ - min_))
        elif colors[i] == 4:
            rgb_color = cmap_4((abs(weights[i]) - min_)/(max_ - min_))
        elif colors[i] == 5:
            rgb_color = cmap_5((abs(weights[i]) - min_)/(max_ - min_))
        rgb_colors.append(rgb_color)

    # make dictionaries to be able to sort
    dict_rgb   = dict(zip(rgb_colors, np.abs(weights)))
    dict_weights = dict(zip(weights, np.abs(weights)))
    dict_edges = dict(zip(edges, np.abs(weights)))

    # sort dictionaries
    rgb_colors = sorted(rgb_colors, key = lambda ele: dict_rgb[ele])
    weights = sorted(weights, key = lambda ele: dict_weights[ele])
    edges = sorted(edges, key = lambda ele: dict_edges[ele])

    # make edges thicker
    for i in range(len(weights)):
        weights[i] = weights[i] * 10.0

    # Select only top 25% to plot
    number = int(len(edges)*.1)
    edges = edges[-number:]
    weights = weights[-number:]
    rgb_colors = rgb_colors[-number:]

    # Draw edges in network
    nx.draw_networkx_edges(G,
            pos = coordinates,
            edgelist = edges,
            width = weights,
            edge_color = rgb_colors
            )

    # Add labels
    label_1 = mpatches.Patch(color = 'blue', label = 'lag 1')
    label_2 = mpatches.Patch(color = 'orange', label = 'lag 2')
    label_3 = mpatches.Patch(color = 'green', label = 'lag 3')
    label_4 = mpatches.Patch(color = 'red', label = 'lag 4')
    label_5 = mpatches.Patch(color = 'purple', label = 'lag 5')

    # Make plot
    plt.legend(handles = [label_1, label_2, label_3, label_4, label_5])
    plt.title("number of columns: " + str(n_col))
    plt.savefig("p = 5, MA = 3, PCA_ns, n_col = " + str(n_col) + ".png")
    plt.show()

    ##########################
    #### PRINT STATEMENTS ####
    ##########################

    # Print most important edges
    n = len(edges) - 1
    print(edges[(n - 10):])

    # Print biggest STPs
    nodes_dict = {}
    for key in nodes:
        for value in node_sizes:
            nodes_dict[key] = value

    myKeys = list(nodes_dict.keys())
    myKeys.sort()
    sorted_nodes_dict = {i : nodes_dict[i] for i in myKeys}

    keys = list(sorted_nodes_dict.keys())
    print(keys[0:9])










