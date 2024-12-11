import umap

import pandas as pd
import networkx as nx

from sklearn.preprocessing import StandardScaler


def umap_embedding(df, **kwargs):
    reducer = umap.UMAP(
        **kwargs
    )
    embedding = pd.DataFrame(
        reducer.fit_transform(
            StandardScaler().fit_transform(df)
        ),
        columns = ['UMAP1', 'UMAP2'],
        index = df.index
    )
    return embedding, reducer.graph_


def louvain_communities(
    umap_graph, 
    resolution,
    weight_name = 'weight'
):
    g = nx.from_numpy_array(umap_graph)
    return nx.community.louvain_communities(
        g, 
        weight = weight_name,
        resolution = resolution
    )


def add_community_annotation(
    embedding, 
    communities, 
    annotation_column
):
    embedding[annotation_column] = ''
    
    for i, indices in enumerate(communities):
        embedding.iloc[list(indices), -1] = str(i)


def annotate_multiple_resolutions(
    embedding,
    umap_graph,
    resolutions
):
    louvain_columns = []
    for resolution in resolutions:
        communities = louvain_communities(
            umap_graph,
            resolution,
            'weight'
        )

        add_community_annotation(
            embedding, 
            communities, 
            f'louvain_{resolution}'
        )
        louvain_columns.append(f'louvain_{resolution}')
    
    return louvain_columns


def rename_column(column, to_replace):
    for replace in to_replace:
        column = column.replace(f'_{replace}', '')
        
    return column


def concatenate_condition_frames(condition_frames, add_key_to_index = False):
    frames = []
    to_replace = condition_frames.keys()
    for key, df in condition_frames.items():
        df = df.rename(
            columns = {
                c: rename_column(c, to_replace) for c in df.columns
            }
        )

        if add_key_to_index:
            df.index = [i + f'_{key}' for i in df.index]

        frames.append(df)
    
    return pd.concat(frames)