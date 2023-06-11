# Imports for calculations
from itertools  import permutations

# Imports for visualization
from matplotlib import lines
from networkx   import circular_layout, DiGraph, draw_networkx_nodes, draw_networkx_edges, draw_networkx_labels, draw_networkx_edge_labels, get_edge_attributes
from numpy      import amax
from numpy      import log2

### 1. Genome Assembly ###

### 1.1. Hamilton

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

### 1.2. Euler ###

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

### 2. Read mapping ###

### 2.1. Read mapping . Bowtie ###

def generate_rotations(string: str):
    '''
    Parameters:
        string (str)
    Returns:
        A list of all rotations of string+'$'
    '''
    return [(string)[i:]+(string)[:i] for i in range(0,len(string))]

def count_character_occurence(string: str, alphabet: dict):
    ''' 
    Parameters:
        string   (string)
        alphabet (dict), key:    characters of the alphabet
                         values: 0 (a counter for the occurences of the characters in )
    Returns:
        counts (list), counts[i] is the i-th occurence of the character in the string at the position i
    '''
    counts   = [0 for i in range(len(string))]
    for i, character in enumerate(string):
        alphabet[character] += 1
        counts[i] = alphabet[character]
    return counts

class Bowtie:

    def __init__(self, genome:str):
        
        genome      = genome if genome.endswith('$') else genome+'$'
        alphabet    = sorted(list(set([character for character in genome])))
        bw_matrix   = sorted(generate_rotations(genome))
        first_col   = ''.join([string[ 0] for string in bw_matrix])
        last_col    = ''.join([string[-1] for string in bw_matrix])
        count       = {character: -1 for character in alphabet}
        genome_i    = count_character_occurence(genome, count.copy())
        first_col_i = count_character_occurence(first_col, count.copy())
        last_col_i  = count_character_occurence(last_col, count.copy())
        bw_indices  = {first_col[i]: i for i in range(len(first_col)) if first_col[i]!=first_col[i-1]}

        self.alphabet                = alphabet
        self.genome                  = genome
        self.genome_i                = genome_i
        self.first_col_i             = first_col_i    
        self.bw_transformed_genome   = last_col
        self.bw_transformed_genome_i = last_col_i
        self.bw_indices              = bw_indices
        self.bw_matrix               = bw_matrix 

    def bw_matrix_formatted(self, indices: list, alignment_length=0):
        '''
        Parameters:
            indices (list),
            alignment_length (int)
        Returns:
            Nicely formatted bw_matrix for displaying purposes
            The indeces tell in which rows an alignment step is to be highlighted
        '''
        return [f'{i if i>9 else " "+str(i)} {"->" if i in indices else "  "} {self.first_col_i[i]} ' +
                f'{row[:alignment_length].upper() + row[alignment_length:] if (i in indices and alignment_length>0) else row} '+
                f'{self.bw_transformed_genome_i[i]} {"<-" if i in indices else "  "}'
                   for i, row in enumerate(self.bw_matrix)]

    def match_next_character(self, last_aligned_character, character_to_align, allowed_i):
        start = self.bw_indices[last_aligned_character]
        end   = min([index for index in self.bw_indices.values() if index>self.bw_indices[last_aligned_character]])
        rows  = [row for row in range(start,end) if character_to_align==self.bw_transformed_genome[row] and self.first_col_i[row] in allowed_i]
        print(f'matched previously: {last_aligned_character}')
        print(f'matches now: {character_to_align}')
        return  [self.bw_transformed_genome_i[i] for i in rows]
    
    def map_read(self, read):
        # Start of the procedure: Mapping of the first character
        first_character = read[-1]
        first   = read[len(read)-1]
        start   = self.bw_indices[first]
        end     = min([index for index in self.bw_indices.values() if index>self.bw_indices[first]])
        allowed_character_indices = list(range(end-start))
        allowed_first_row_indices = list(range(start, end))
        out = [print(row) for row in self.bw_matrix_formatted(allowed_first_row_indices, alignment_length=1)]
        
        # Iterating remaining steps
        reversed_read = read[::-1]
        for i, character in enumerate(reversed_read):
            if i+1 < len(reversed_read):
                print()
                allowed_character_indices = self.match_next_character(reversed_read[i], reversed_read[i+1], allowed_character_indices)
                allowed_first_row_indices = [self.bw_indices[read[i]]+j for j in allowed_character_indices]
                out = [print(row) for row in self.bw_matrix_formatted(allowed_first_row_indices, alignment_length=i+2)]
