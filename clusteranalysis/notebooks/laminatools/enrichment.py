import gc
import logging

import gseapy as gp
import pandas as pd

from scipy.stats import t
from statsmodels.stats.multitest import fdrcorrection


def get_enriched_bins(data, alpha = 0.05):
    # background = df.loc[df[data_column] < df[data_column].mean(), :]
    # background = df.loc[np.random.choice(df.index, int(len(df) * 0.25), replace = False), :]
    
    # we compute the mean separately from the full data frame 
    # to avoid biasing the distribution to higher values
    # fitting it to the full data frame results in nonsense parameters
    # thus we fit only to non zero values
    # std/scale is similar in both cases (full and nonzero) so we use the fitted one here
    loc = data.mean()
    data = data[data > 0]
    degrees, _, scale = t.fit(data)
    
    pval = lambda x: t.sf(x, degrees, loc = loc, scale = scale)
    
    pvals = data.apply(pval)
    reject, _ = fdrcorrection(pvals, alpha = alpha)
    return data[reject].index


    
# concurrent web service use is tricky due to request limitations etc.
# we use offline instead
def gp_prerank(expression, gene_sets, fdr = 0.1, identifier = ''):
    try:
        results = gp.prerank(
            rnk = expression,
            gene_sets=gene_sets,
            outdir = None,
            min_size = 5,
            max_size = 2000
        )
        
    except Exception as e:
        logging.info('prerank {identifier} {exception}'.format(identifier = identifier, exception = str(e)))
        return pd.DataFrame()
        
    results = pd.DataFrame() \
        .from_dict(results.results, orient = 'index') \
        .reset_index(names = 'Term')

    results = results.loc[results.fdr < fdr, :]
    return results.drop(columns = ['RES'])
    

def gp_enrichment(gene_list, gene_sets, background, fdr = 0.1, identifier = ''):
    try:
        enrichment = gp.enrichr(
            gene_list=gene_list,
            gene_sets=gene_sets,
            background=background,
            outdir = None
        ).results
    
    except ValueError:
        logging.warning(f'no results for {identifier}! gene list was {str(gene_list)}')
        return pd.DataFrame()
    
    return enrichment.loc[enrichment['Adjusted P-value'] < fdr, :]
    
    
def do_enrichment_and_prerank(
    gene_list, 
    ranked_genes, 
    gene_sets,
    background,
    identifier, 
    gene_set,
    cluster_column,
    cluster,
    condition
    
):
    # we first wait for a random amount of time to not spam the server
    logging.info(f'{identifier} enrichr {gene_set}')
    cluster_enrichment = gp_enrichment(gene_list, gene_sets, background, 0.1, identifier = identifier)

    logging.info(f'{identifier} prerank {gene_set}')
    cluster_prerank = gp_prerank(
        ranked_genes,
        gene_sets,
        0.1,
        identifier
    )

    if not cluster_enrichment.empty:
        cluster_enrichment.to_csv(
            f'../counts/gene_clusterings_{cluster_column}_cluster_{cluster}_enrichment_{gene_set}_{condition}.tsv',
            sep = '\t',
            index = False
        )

    if not cluster_prerank.empty:
        cluster_prerank.to_csv(
            f'../counts/gene_clusterings_{cluster_column}_cluster_{cluster}_prerank_{gene_set}_{condition}.tsv',
            sep = '\t',
            index = False
        )
        
    del cluster_enrichment, cluster_prerank
    gc.collect()
        
    
def get_gene_ranking(gene_list, expressed_genes):
    ranked_genes = expressed_genes \
        .loc[expressed_genes.index.isin(gene_list), :] \
        .groupby('gene_name') \
        .max() \
        .sort_values(by = 'rnaseq', ascending = False)
    
    return ranked_genes
    
    
def compute_cluster_enrichments(
    group_tuple, 
    expressed_genes, 
    id_name_map, 
    differential_genes, 
    lap2_enriched_genes,
    gene_sets,
    condition,
    cluster_column
):
    cluster, group = group_tuple
    identifier = ' '.join([condition, cluster_column, cluster])
    logging.warning(f'{identifier}')
    group = group.merge(
        id_name_map,
        how = 'inner',
        left_index = True,
        right_index = True
    )
    
    expressed_group = group.loc[group.index.isin(expressed_genes.index), :]
    do_enrichment_and_prerank(
        list(expressed_group.gene_name.unique()),
        get_gene_ranking(group.index, expressed_genes),
        gene_sets,
        list(expressed_genes.gene_name.unique()),
        identifier,
        'full',
        cluster_column,
        cluster,
        condition
    )

    lap2_group = expressed_group.loc[expressed_group.index.isin(lap2_enriched_genes), :]
    do_enrichment_and_prerank(
        list(lap2_group.gene_name.unique()),
        get_gene_ranking(lap2_group.index, expressed_genes),
        gene_sets,
        list(expressed_genes.gene_name.unique()),
        identifier,
        'lap2',
        cluster_column,
        cluster,
        condition
    )
    
    del expressed_group, lap2_group
    gc.collect()
    
    for regulation, diff_genes in differential_genes.items():
        regulated_group = group.loc[group.gene_name.isin(diff_genes), :]        
        do_enrichment_and_prerank(
            list(regulated_group.gene_name.unique()),
            get_gene_ranking(regulated_group.index, expressed_genes),
            gene_sets,
            list(expressed_genes.gene_name.unique()),
            identifier,
            regulation,
            cluster_column,
            cluster,
            condition
        )
        
        lap2_regulated_group = regulated_group.loc[regulated_group.index.isin(lap2_enriched_genes), :]
        do_enrichment_and_prerank(
            list(lap2_regulated_group.gene_name.unique()),
            get_gene_ranking(lap2_regulated_group, expressed_genes),
            gene_sets,
            list(expressed_genes.gene_name.unique()),
            identifier,
            f'lap2_{regulation}',
            cluster_column,
            cluster,
            condition
        )
        
        del regulated_group, lap2_regulated_group
        gc.collect()
        