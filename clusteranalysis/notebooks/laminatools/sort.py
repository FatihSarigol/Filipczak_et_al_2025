import numpy as np
from sklearn.decomposition import PCA


def column_names_to_mask(data_columns, columns_of_interest):
    column_mask = np.ones(shape = len(data_columns))
    
    if isinstance(columns_of_interest, type(None)):
        return column_mask
    
    column_set = set(columns_of_interest)
    print(data_columns)
    for i, col in enumerate(data_columns):
        if not col in column_set:
            column_mask[i] = 0
    
    return column_mask
    
    
def transform_to_sort_data(
    data, 
    use_column_mask, 
    negate_column_mask
):
    transformed_data = []
    for value, use_value, negate_value in zip(
        data, 
        use_column_mask, 
        negate_column_mask
    ):
        if negate_value:
            value = -value
        
        if use_value:
            transformed_data.append(value)
    
    return tuple(transformed_data)


def remove_none_cluster(signals_matrix, row_names):
    filtered_signals_matrix, filtered_row_names = [], []
    for i, row_name in enumerate(row_names):
        if row_name == 'None':
            continue
        
        filtered_signals_matrix.append(signals_matrix[i, :])
        filtered_row_names.append(row_name)
    
    return np.vstack(filtered_signals_matrix), np.array(filtered_row_names)


def sort_clusters_on_data_signals(
    signals_matrix, 
    row_labels, 
    column_labels, 
    use_columns = None, 
    negate_columns = None,
    consider_none = False
):
    use_col_mask = column_names_to_mask(column_labels, use_columns)
    negate_col_mask = column_names_to_mask(column_labels, negate_columns)
    
    if not consider_none:
        signals_matrix, row_labels = remove_none_cluster(
            signals_matrix,
            row_labels
        )
        
    sort_input = [
        (i, row_label, row_data)
        for i, (row_label, row_data)
        in enumerate(zip(row_labels, signals_matrix))
    ]
    
    sorted_row_labels, sorted_signal_matrix, sorted_index = [], [], []
    for i, row_label, row in sorted(
        sort_input, 
        key = lambda x: transform_to_sort_data(x[-1], use_col_mask, negate_col_mask),
        reverse = True
    ):
        sorted_row_labels.append(row_label)
        sorted_signal_matrix.append(row)
        sorted_index.append(i)
    
    return np.array(sorted_index), np.array(sorted_row_labels), np.vstack(sorted_signal_matrix)
        

def sort_clusters_on_data_signals_pca(
    signals_matrix, 
    row_labels, 
    column_labels, 
    use_columns = None, 
    negate_columns = None,
    consider_none = False
):
    use_col_mask = column_names_to_mask(column_labels, use_columns)
    negate_col_mask = column_names_to_mask(column_labels, negate_columns)
    
    if not consider_none:
        signals_matrix, row_labels = remove_none_cluster(
            signals_matrix,
            row_labels
        )
        
    pca_input = np.vstack(
        [
            np.array(transform_to_sort_data(data, use_col_mask, negate_col_mask))
            for data 
            in signals_matrix
        ]
    )
    
    pca = PCA(n_components = 1)
    pc = pca.fit_transform(
        pca_input
    )
    
    sort_idx = np.argsort(pc.flatten())
    return sort_idx, np.array(row_labels)[sort_idx], signals_matrix[sort_idx, :]