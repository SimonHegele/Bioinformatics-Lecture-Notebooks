# Imports for calculations
from itertools  import permutations

# Imports for visualization
from matplotlib import lines
from networkx   import circular_layout, DiGraph, draw_networkx_nodes, draw_networkx_edges, draw_networkx_labels, draw_networkx_edge_labels, get_edge_attributes
from numpy      import amax
from numpy      import log2

def maximum_overlap(string_1: str, string_2: str):
    '''
    Returns: length of the longest suffix of string_1 that is also prefix of string_2
    '''
    n = min(len(string_1),len(string_2))
    x = 0
    for i in range(n):
        if string_1[n-i:]==string_2[:i]:
            x=i
    return x

def maximum_overlap_matrix(strings: list):
    '''
    Returns the matrix of maximum overlaps for all (ordered) pairs of strings in strings (implies a directed graph)
    '''
    return [[maximum_overlap(string_1, string_2) for string_2 in strings] for string_1 in strings]

def draw_graph_from_matrix(strings, scoring_matrix, axes,  label=False):
    # Creating the Graph object
    G = DiGraph()
    nodes     = strings
    edge_list = [(strings[i],strings[j]) for i in range(len(strings)) for j in range(len(strings)) if not strings[i]==strings[j] if scoring_matrix[i][j]]
    G.add_nodes_from(nodes)
    G.add_edges_from(edge_list)
    # Drawing the Graph object
    # 1. Nodes
    draw_networkx_nodes(G, pos=circular_layout(G), node_size=2500*log2(len(strings[0])), ax=axes)
    draw_networkx_labels(G, pos=circular_layout(G), font_size=25, ax=axes)
    # 2. Edges
    colors = ['#000000', '#ff0000','#00ff00','#0000ff']
    l  = []
    t  = []
    for i in reversed(range(1, amax(scoring_matrix)+1)):
        edge_list_i   = [edge for edge in edge_list if scoring_matrix[strings.index(edge[0])][strings.index(edge[1])]==i]
        edge_labels_i = {edge: i for edge in edge_list_i}
        draw_networkx_edges(G, pos=circular_layout(G), edgelist=edge_list_i, node_size=2500*log2(len(strings[0])), connectionstyle="arc3,rad=0.1", edge_color=colors[i-1], width=i*2, ax=axes)
        l.append(lines.Line2D(range(1), range(1), color="white", marker="_", mec=colors[i-1], markersize=30, markeredgewidth=3*i))
        t.append(f'Max overlap = {i}')
        if label:
            draw_networkx_edge_labels(G, pos=circular_layout(G), edge_labels=edge_labels_i, font_size=20, label_pos=0.25, font_color=colors[i-1], ax=axes)
    # Legend
    axes.legend(tuple(l),tuple(t), fontsize=15)

def path_weight(path: list, scoring_matrix: list):
    '''
    Returns: Sum of weight of all edges in path in a directed graph implied by the scoring_matrix
    '''
    weight = 0
    for i in range(len(path)-1):
        weight += scoring_matrix[path[i]][path[i+1]]
    return weight

def brute_force_best_hamilton_paths(scoring_matrix):
    '''
    Returns a list with the hamilton paths with the highest sum of edge weights by trying every permutation of the nodes
    '''
    candidates    = list(permutations(range(len(scoring_matrix))))
    best_score    = 0
    current_bests = []
    for candidate in candidates:
        candidate_score = path_weight(candidate, scoring_matrix)
        if candidate_score > best_score:
            best_score = candidate_score
            current_bests = [candidate]
        elif candidate_score == best_score:
            current_bests.append(candidate)
    return current_bests, best_score

def get_corresponding_strings(path: list, strings: list):
    return [strings[index] for index in path]

def merge_strings(string_1: str, string_2: str):
    '''
    Returns the merged string of string_1 and string_2 with maximum overlap (pairwise merging)
    '''
    for i in range(1,len(string_1)+1):
        for j in reversed(range(1,len(string_2))):
            if string_1[i:]==string_2[:j]:
                return string_1 + string_2[j:]
    return string_1+string_2

def merge_strings_from_list(strings: list):
    s = strings[0]
    for i in range(1,len(strings)):
        s  = merge_strings(s, strings[i])
    return s

def k_mer_graph(strings: list, k: int):
    G = DiGraph()
    for string in strings:
        for i in range(0, len(string)-k+1):
            if not string[i:i+k] in G.nodes():
                G.add_node(string[i:i+k])
            if i > 0:
                if (string[i-1:i+k-1], string[i:i+k]) not in G.edges():
                    G.add_edge(string[i-1:i+k-1], string[i:i+k], weight=1)
                else:
                    G[string[i-1:i+k-1]][string[i:i+k]]['weight'] += 1
    return G

def draw_k_mer_graph(G, axes):
    pos = circular_layout(G)
    # 1. Nodes
    draw_networkx_nodes(G, pos=pos, node_size=2500*log2(len(list(G.nodes())[0])), ax=axes)
    # 2. Node labels
    draw_networkx_labels(G, pos, font_size=25, ax=axes)
    # 3. Edges
    draw_networkx_edges(G, pos=pos, node_size=2500*log2(len(list(G.nodes())[0])), width=3, ax=axes)
    # 4. Edge labels
    draw_networkx_edge_labels(G, pos, edge_labels=get_edge_attributes(G, 'weight'), font_size=25, ax=axes)

def draw_k_mer_graph(G, axes):
    pos = circular_layout(G)
    # 1. Nodes
    draw_networkx_nodes(G, pos=pos, node_size=2500*log2(len(list(G.nodes())[0])), ax=axes)
    # 2. Node labels
    draw_networkx_labels(G, pos, font_size=25, ax=axes)
    # 3. Edges
    draw_networkx_edges(G, pos=pos, node_size=2500*log2(len(list(G.nodes())[0])), width=3, ax=axes)
    # 4. Edge labels
    draw_networkx_edge_labels(G, pos, edge_labels=get_edge_attributes(G, 'weight'), font_size=25, ax=axes)

def generate_rotations(string: str):
    '''
    Parameters:
        string (str)
    Returns:
        A list of all rotations of string+'$'
    '''
    return [(string+'$')[i:]+(string+'$')[:i] for i in range(0,len(string+'$'))]

def burrows_wheeler_indices(strings: list):
    return {strings[i][0]: i for i in range(1, len(strings)) if strings[i][0]!=strings[i-1][0]}