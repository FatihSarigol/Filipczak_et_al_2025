import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import pylluvial as pa
import pandas as pd

from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sklearn.preprocessing import StandardScaler

from . import sort


def add_margin(xmin, xmax, ymin, ymax):
    margin = ((xmax - xmin) + (ymax - ymin)) / 2 * 0.1
    return xmin - margin, xmax + margin, ymin - margin, ymax + margin


def get_categorical_colors(categories):
    color_palette = sns.color_palette(
        'husl', 
        len(categories)
    ).as_hex()
    return {cat: color_palette[i] for i, cat in enumerate(categories)}


def plot_embedding_categorical(
    embedding, 
    label_column, 
    ax, 
    color_palette = None,
    consider_none = True
):
    plot_data = embedding.loc[:, ['UMAP1', 'UMAP2', label_column]]

    if not consider_none:
        plot_data = plot_data.loc[plot_data[label_column] != 'None', :]

    categories = sorted(plot_data[label_column].unique())
    if not color_palette:
        color_palette = get_categorical_colors(categories)
    
    for cat in categories:
        data = plot_data.loc[embedding[label_column] == cat, :]
        ax.scatter(
            data['UMAP1'],
            data['UMAP2'],
            c = color_palette[cat],
            label = cat
        )
    
    ax.legend()


def plot_embedding_scalar(
    embedding, 
    label_column, 
    fig, 
    ax
):
    color_palette = sns.color_palette('viridis', as_cmap = True)
    q05, q95 = np.percentile(embedding[label_column], [5, 95])
    norm = Normalize(q05, q95, clip = True)
    ax.scatter(
        embedding['UMAP1'],
        embedding['UMAP2'],
        c = embedding[label_column],
        cmap = color_palette,
        norm = norm
    )
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    smap = plt.cm.ScalarMappable(cmap=color_palette, norm=norm)
    fig.colorbar(smap, cax=cax, orientation='vertical')


def plot_embedding_per_category(
    embedding, 
    label_column, 
    savefile,
    color_palette = None, 
    consider_none = True
):
    plot_data = embedding.loc[:, ['UMAP1', 'UMAP2', label_column]]

    if not consider_none:
        plot_data = plot_data.loc[plot_data[label_column] != 'None', :]

    categories = sorted(plot_data[label_column].unique())

    fig, axs = plt.subplots(1, len(categories))
    xmin, xmax, ymin, ymax = add_margin(
        embedding['UMAP1'].min(),
        embedding['UMAP1'].max(),
        embedding['UMAP2'].min(),
        embedding['UMAP2'].max()
    )

    if not color_palette:
        color_palette = get_categorical_colors(categories)
    
    for ax, cat in zip(axs, categories):
        data = plot_data.loc[plot_data[label_column] == cat, :]
        ax.scatter(
            data['UMAP1'],
            data['UMAP2'],
            c = color_palette[cat]
        )
        
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_title(cat)
        ax.set_xticks([])
        ax.set_yticks([])
    
    fig.set_figwidth(len(categories) * 5)
    fig.set_figheight(5)
    fig.tight_layout()
    fig.savefig(savefile, dpi = 300)
    plt.close(fig)


def plot_and_save_annotated_embedding(
    embedding,
    annotation_columns,
    plotfile,
    plot_per_category = False,
    palettes = None,
    consider_none = True
):
    n_panels_per_row = 4
    n_rows = int(np.ceil(len(annotation_columns) / n_panels_per_row))
    fig, axs = plt.subplots(n_rows, n_panels_per_row)
    for ax, annotation_column in zip(
        axs.reshape(n_rows * n_panels_per_row), 
        annotation_columns
    ):
        if any(annotation_column.startswith(s) for s in ['chromhmm', 'louvain']):
            if not palettes:
                categories = sorted(embedding[annotation_column].unique())
                color_palette = get_categorical_colors(categories)

            else:
                color_palette = palettes[annotation_column]

            plot_embedding_categorical(
                embedding, 
                annotation_column, 
                ax, 
                color_palette,
                consider_none = consider_none
            )          

            if plot_per_category:
                plot_embedding_per_category(
                    embedding,
                    annotation_column,
                    plotfile.replace('.png', f'{annotation_column}.individual.png'),
                    color_palette,
                    consider_none = consider_none
                )
            
        else:
            plot_embedding_scalar(
                embedding, 
                annotation_column, 
                fig, 
                ax
            )
            
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(annotation_column)

        
    fig.set_figwidth(7 * n_panels_per_row)
    fig.set_figheight(7 * n_rows)
    fig.tight_layout()
    # can't increase dpi too much with png since it has an intrinsic limit of 2^16 pixels per dimension
    # see https://github.com/matplotlib/matplotlib/issues/16234
    fig.savefig(plotfile, dpi = 300)
    plt.close(fig)


def make_diagnostic_plot(
    p1,
    p2,
    nn_graph,
    correspondences,
    row_names_this, 
    row_names_that,
    label_this,
    label_that,
    axs
):        
    ax1, ax2 = axs
    X = np.concatenate([p1[:, :2], p2[:, :2]])
    x_axis_size = X[:, 0].max() - X[:, 0].min()
    
    for p, row_names, label in zip(
        [p1, p2], 
        [row_names_this, row_names_that], 
        [label_this, label_that]
    ):
        ax1.scatter(p[:, 0], p[:, 1], label = label)
        for name, (x, y) in zip(row_names, p[:, :2]):
            ax1.text(
                x + x_axis_size * 0.03, 
                y, 
                str(name), 
                ha = 'center', 
                va = 'center'
            )
    
    ax1.legend()
    
    nx.draw_networkx_edges(
        nn_graph,
        X,
        ax = ax2
    )
    
    for c, correspondence_dict in correspondences.items():
        points = correspondence_dict['points']
        ax2.scatter(points[:, 0], points[:, 1], label = c)
        
    ax2.legend()
    
    for ax in [ax1, ax2]:
        ax.set_yticks([])
        ax.set_xticks([])
        ax.set_ylabel('PC2')
        ax.set_xlabel('PC1')


def plot_recluster_data(
    correspondences, 
    matrix_this, 
    matrix_that, 
    row_names_this, 
    row_names_that, 
    label_this, 
    label_that, 
    gs,
    fig
):
    scaler = StandardScaler()
    m_this_scaled = scaler.fit_transform(matrix_this)
    m_that_scaled = scaler.fit_transform(matrix_that)
    for i, (c_label, correspondence) in enumerate(correspondences.items()):
        ax = fig.add_subplot(gs[i + 1, :])
        this_idx = correspondence['p1_idx']
        that_idx = correspondence['p2_idx']
        m = np.concatenate(
            [
                m_this_scaled[this_idx, :],
                m_that_scaled[that_idx, :]
            ]
        )
        
        ax.imshow(m, aspect = 'auto', vmin = -1, vmax = 1)
        ax.set_ylabel(i)
        ax.set_xticks([])
        ax.set_yticks(range(m.shape[0]))
        
        yticklabels_this = [f'{label_this}_{rname}' for rname in row_names_this[this_idx]]
        yticklabels_that = [f'{label_that}_{rname}' for rname in row_names_that[that_idx]]
        yticklabels = np.concatenate(
            [
                yticklabels_this,
                yticklabels_that
            ]
        )
        ax.set_yticklabels(yticklabels)
        ax.set_ylabel(c_label)


def plot_matrix(matrix, row_names, col_names, ax, vmin, vmax):
    ax.imshow(
        matrix,
        vmin = vmin,
        vmax = vmax,
        cmap = 'coolwarm',
        interpolation = None,
        aspect = 'auto'
    )
    ytick_positions = list(range(len(row_names)))
    xtick_positions = list(range(len(col_names)))
    
    for xpos in xtick_positions:
        for ypos in ytick_positions:
            ax.text(
                xpos, 
                ypos, 
                '{:.2f}'.format(matrix[ypos, xpos]),
                ha = 'center',
                va = 'center'
            )
            
    ax.set_yticks(ytick_positions)
    ax.set_yticklabels(row_names)
    ax.set_xticks(xtick_positions)
    ax.set_xticklabels(col_names, rotation = 90)


def generate_correspondence_plot(
    binsize_correlations_this, 
    binsize_correlations_that, 
    data, 
    clustering, 
    this_condition, 
    that_condition, 
    plotfile,
    vmin,
    vmax,
    consider_none = False,
    palette = 'husl'
):
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    
    matrix_this = binsize_correlations_this[clustering]['matrix']
    data_labels_this = binsize_correlations_this[clustering]['col_names']
    row_names_this = binsize_correlations_this[clustering]['row_names']
    matrix_that = binsize_correlations_that[clustering]['matrix']
    data_labels_that = binsize_correlations_that[clustering]['col_names']
    row_names_that = binsize_correlations_that[clustering]['row_names']
    
    if not consider_none:
        matrix_this, row_names_this = sort.remove_none_cluster(matrix_this, row_names_this)
        matrix_that, row_names_that = sort.remove_none_cluster(matrix_that, row_names_that)
        data = data.loc[~(data[clustering] == 'None'), :]
    
    plot_matrix(
        matrix_this,
        row_names_this,
        data_labels_this,
        ax1,
        vmin = vmin,
        vmax = vmax
    )
    ax1.set_title(this_condition)

    pa.alluvial(
        data = data,
        x = 'condition',
        stratum = clustering,
        alluvium = 'bin',
        ax = ax2,    
        palette = palette,
        stratum_gap = 2,
        stratum_width = 2,
        show_labels = True
    )

    plot_matrix(
        matrix_that,
        row_names_that,
        data_labels_that,
        ax3,
        vmin = vmin,
        vmax = vmax
    )
    ax3.set_title(that_condition)

    fig.set_figwidth(30)
    fig.set_figheight(7)
    fig.tight_layout()
    fig.savefig(plotfile)
    plt.close(fig)


def plot_cluster_correspondence(
    binsize_correlations_this, 
    binsize_correlations_that, 
    this_data, 
    that_data, 
    clustering, 
    this_condition, 
    that_condition,
    plotfile,
    vmin,
    vmax,
    subset = None,
    consider_none = False,
    palette = 'husl'
):
    if not subset:
        this_cluster_data = this_data[[clustering]].reset_index(names = 'bin')
        that_cluster_data = that_data[[clustering]].reset_index(names = 'bin')

    else:
        this_cluster_data = this_data.loc[subset, [clustering]].reset_index(names = 'bin')
        that_cluster_data = that_data.loc[subset, [clustering]].reset_index(names = 'bin')
    
    this_cluster_data['condition'] = f'1_{this_condition}'
    that_cluster_data['condition'] = f'2_{that_condition}'
    data = pd.concat([this_cluster_data, that_cluster_data])
    
    generate_correspondence_plot(
        binsize_correlations_this,
        binsize_correlations_that,
        data,
        clustering,
        this_condition,
        that_condition,
        plotfile,
        vmin = vmin,
        vmax = vmax,
        consider_none = consider_none,
        palette = palette
    )


