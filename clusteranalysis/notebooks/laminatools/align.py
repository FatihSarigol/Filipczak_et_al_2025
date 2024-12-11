import pandas as pd
import numpy as np
import networkx as nx

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.spatial.distance import cdist
from pybedtools import BedTool


def get_majority_assignment(group, assignment_column):
    return group.groupby(assignment_column).agg({'olap': 'sum'}).idxmax()[0]


def annotate_hmm_states(genes, states, assignment_column):
    genes = BedTool().from_dataframe(
        genes.loc[:, ['Chr', 'Start', 'End', 'Geneid']]
    )

    annotated_genes = genes.intersect(
        states,
        wao = True
    ).to_dataframe(
        names = ['chrom', 'start', 'end', 'name', 'chrom2', 'start2', 'end2', assignment_column, 'olap']
    )
    processed_assignments = []
    for _, group in annotated_genes.groupby('name', sort = False):
        if len(group) == 1:
            if group[assignment_column].values[0] == '.':
                group.loc[:, assignment_column] = 'None'

            processed_assignments.append(group)
            continue

        return_row = group.iloc[0].copy()
        return_row[assignment_column] = get_majority_assignment(
            group, 
            assignment_column
        )
        return_row.start2 = -1
        return_row.end2 = -1
        return_row.olap = -1

        processed_assignments.append(pd.DataFrame(return_row).T)

    return pd.concat(processed_assignments)


def align_columns(data_matrix_this, col_names_this, data_matrix_that, col_names_that):
    col_sort_idx_this = np.argsort(col_names_this)
    col_sort_idx_that = np.argsort(col_names_that)
    return (
        data_matrix_this[:, col_sort_idx_this],
        col_names_this[col_sort_idx_this],
        data_matrix_that[:, col_sort_idx_that],
        col_names_that[col_sort_idx_that]
    )


def build_nn_graph(p1, p2, metric = 'cityblock'):
    n1, n2 = p1.shape[0], p2.shape[0]
    n_nodes = n1 + p2.shape[0]
    adjacency_matrix = np.zeros(shape = (n_nodes, n_nodes), dtype = int)
    distances_p1_p2 = cdist(p1, p2, metric)
    for i in range(n1):
        nearest_neighbour = distances_p1_p2[i, :].argmin()
        adjacency_matrix[i, nearest_neighbour + n1] = 1
        adjacency_matrix[nearest_neighbour + n1, i] = 1
    
    # doing the same from perspective of p2 to avoid orphan clusters
    distances_p2_p1 = distances_p1_p2.T
    for i in range(n2):
        nearest_neighbour = distances_p2_p1[i, :].argmin()
        adjacency_matrix[n1 + i, nearest_neighbour] = 1
        adjacency_matrix[nearest_neighbour, n1 + i] = 1
        
    nn_graph = nx.from_numpy_array(adjacency_matrix)
    return nn_graph


def cluster_correspondences(nn_graph, p1, p2):
    n1 = p1.shape[0]
    components = nx.connected_components(nn_graph)
    correspondences = {}
    for new_cluster_label, component in enumerate(components):
        p1_idx, p2_idx = [], []
        for i in component:
            # if i is larger equal to n1 we are already in p2
            if i > n1 - 1:
                i = i - n1
                p2_idx.append(i)
            
            else:
                p1_idx.append(i)
        
        points = np.concatenate([p1[p1_idx, :], p2[p2_idx, :]])
        correspondences[new_cluster_label] = {
            'points': points,
            'p1_idx': p1_idx,
            'p2_idx': p2_idx
        }
    
    return correspondences


def pca(
    matrix_this, 
    matrix_that, 
    scale = True, 
    n_components = 2
):
    pca = PCA(n_components = n_components)

    if scale:
        scaler = StandardScaler()
        matrix_this = scaler.fit_transform(matrix_this)
        matrix_that = scaler.fit_transform(matrix_that)

    n1 = matrix_this.shape[0]
    m = np.concatenate([matrix_this, matrix_that])
    p = pca.fit_transform(m)
    p1 = p[:n1, :]
    p2 = p[n1:, :]
    # p1 = pca.fit_transform(matrix_this)
    # p2 = pca.fit_transform(matrix_that)
    return p1, p2
    

def nearest_neighbour_correspondence(
    matrix_this, 
    matrix_that,
    scale = True, 
    n_components = 2,
    metric = 'cityblock'
):
    p1, p2 = pca(
        matrix_this, 
        matrix_that, 
        scale, 
        n_components
    )
        
    graph = build_nn_graph(p1, p2, metric)
    correspondences = cluster_correspondences(graph, p1, p2)
    
    return p1, p2, graph, correspondences