import os
import gc

import pandas as pd
import pyBigWig as pb

from . import align
from pybedtools import BedTool


def make_bin_frame_from_chrom_sizes(chrom_sizes, binsize):
    bins = []
    for _, chrom in chrom_sizes.iterrows():
        chrom_name = chrom.chrom
        chrom_size = chrom.length
        bins.extend(
            [(chrom_name, start, start + binsize if start + binsize < chrom_size else chrom_size) for start in range(0, chrom_size, binsize)]
        )
    
    return pd.DataFrame(bins, columns = ['chrom', 'start', 'end'])


def read_bigwig(bigwigfile, bin_frame, binsize):
    print(bigwigfile)
    bin_frame = bin_frame.copy()
    column_name = '_'.join(os.path.basename(bigwigfile).split('_')[:-1])
    
    bw = pb.open(bigwigfile)
    
    intervals = []
    for chrom in bin_frame.chrom.unique():
        for interval in bw.intervals(chrom):
            start, end, value = interval
            if not value:
                continue
                
            if end - start > binsize:
                expanded_interval = [
                    (
                        chrom, start, start + binsize 
                        if start + binsize < end 
                        else end, value
                    ) 
                    for start 
                    in range(start, end, binsize)
                ]
            
            else:
                expanded_interval = [(chrom, start, end, value)]
            
            intervals.extend(expanded_interval)
    
    intervals = pd.DataFrame(
        intervals, 
        columns = [
            'chrom', 
            'start', 
            'end', 
            column_name
        ]
    )
    
    bw.close()
    
    bin_frame = bin_frame.merge(
        intervals,
        on = ['chrom', 'start', 'end'],
        how = 'left'
    )
    bin_frame.index = bin_frame.apply(
        lambda chrom_bin: f'{chrom_bin.chrom}:{chrom_bin.start}-{chrom_bin.end}',
        axis = 1
    )
    return bin_frame.drop(columns = ['chrom', 'start', 'end'])


def read_segmentation(states, segment_string, binsize, condition):
    segmentationfile = segment_string.format(
        condition = condition,
        binsize = binsize,
        states = states
    ) 
    state_column = f'chromhmm_{states}'
    
    print(segmentationfile, state_column)
    
    segmentation = pd.read_csv(
        segmentationfile,
        sep = '\t',
        header = None,
        names = ['chrom', 'start', 'end', state_column]
    )
    state_per_bin = []
    for _, segment in segmentation.iterrows():
        if segment.end - segment.start > binsize:
            state_per_bin.extend(
                [
                    (
                        segment.chrom, 
                        start, 
                        start + binsize if start + binsize < segment.end else segment.end, 
                        segment[state_column]
                    ) 
                    for start 
                    in range(segment.start, segment.end, binsize)

                ]
            )
        else:
            state_per_bin.append(
                (
                    segment.chrom, 
                    segment.start, 
                    segment.end, 
                    segment[state_column]
                )
            )

    state_per_bin = pd.DataFrame(
        state_per_bin, 
        columns = [
            'chrom', 
            'start', 
            'end', 
            state_column
        ]
    )

    del segmentation
    gc.collect()

    state_per_bin.index = state_per_bin.apply(
        lambda x: f'{x.chrom}:{x.start}-{x.end}',
        axis = 1
    )
    return state_per_bin.drop(columns = ['chrom', 'start', 'end'])


def write_categorical_colors(palette, key, filename):
    with open(filename, 'a') as colorfile:
        line = '{key};{color_list}\n'
        colorfile.write(
            line.format(
                key = key,
                color_list = ','.join(
                    [f'{k}:{v}' for k, v in palette.items()]
                )
            )
        )


def read_chromhmm_and_transfer_to_genes(
    file_format_string,
    state_range,
    genes
):
    state_frames = []
    for n in range(*state_range):
        annotation_column = f'chromhmm_{n}'
        states = BedTool(file_format_string.format(n_states = n))
        annotated_genes = align.annotate_hmm_states(
            genes,
            states,
            annotation_column
        )
        state_frames.append(
            annotated_genes.loc[:, ['name', annotation_column]].set_index('name')
        )

    return pd.concat(state_frames, axis = 1)
