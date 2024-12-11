import numpy as np

from scipy.stats import spearmanr


def cluster_correlation_vector(clusters, cluster):
    correlation_vector = np.ones(shape = len(clusters))
    correlation_vector[clusters != cluster] = -1
    return correlation_vector


def cluster_data_correlation(df, cluster_column, data_column):
    clusters = df[cluster_column].values
    cluster_correlations = {}
    for cluster in np.unique(clusters):
        correlation_vector = cluster_correlation_vector(clusters, cluster)
        r, _ = spearmanr(correlation_vector, df[data_column].values)
        cluster_correlations[cluster] = r
    
    return cluster_correlations



def cluster_data_enrichment(df, cluster_column, data_column):
    cluster_enrichments = {}
    background_signal = df[data_column].mean()
    for cluster, cluster_group in df.groupby(cluster_column):
        cluster_signal = cluster_group[data_column].mean()
        cluster_enrichments[cluster] = cluster_signal / background_signal
    
    return cluster_enrichments


def order_vector(vector, sorted_clusters):
    return np.array([vector[cluster] for cluster in sorted_clusters])


def compute_metric_over_clusters(func, df, cluster_column, data_columns):
    sorted_clusters = sorted(df[cluster_column].unique())
    matrix = np.zeros(
        shape = (
            len(sorted_clusters),
            len(data_columns)
        )
    )
    for j, data_column in enumerate(data_columns):
        cluster_results = func(
            df,
            cluster_column,
            data_column
        )
        matrix[:, j] = order_vector(
            cluster_results,
            sorted_clusters
        )
    
    return matrix, sorted_clusters


def data_metric_over_clusterings(
    metric_func,
    df, 
    cluster_columns, 
    data_columns, 
    log_results = True
):
    cluster_metrics = {}
    for cluster_column in cluster_columns:
        matrix, row_names = compute_metric_over_clusters(
            metric_func,
            df, 
            cluster_column, 
            data_columns
        )
        
        if log_results:
            matrix[matrix == 0] = 1
            matrix = np.log2(matrix)
            
        cluster_metrics[cluster_column] = {
            'matrix': matrix,
            'row_names': row_names,
            'col_names': data_columns.copy()
        }
    return cluster_metrics


def data_enrichment_over_clusterings(
    df, 
    cluster_columns, 
    data_columns, 
    log_results = True
):
    cluster_enrichments = data_metric_over_clusterings(
        cluster_data_enrichment,
        df, 
        cluster_columns, 
        data_columns, 
        log_results = log_results
    )
    return cluster_enrichments


def data_correlation_over_clusterings(
    df, 
    cluster_columns, 
    data_columns
):
    cluster_correlations = data_metric_over_clusterings(
        cluster_data_correlation,
        df, 
        cluster_columns, 
        data_columns, 
        log_results = False
    )
    return cluster_correlations