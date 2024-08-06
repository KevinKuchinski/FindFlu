import argparse as arg
import subprocess as sp
import os as os
import numpy as np
import pandas as pd

def main():
    version = '1.0.1'
    print(f'\nFindFlu v{version}')
    print(f'https://github.com/KevinKuchinski/FindFlu/\n')
    args = get_args()
    print(f'Query fragment end sequences: {args.input_file}')
    print(f'Reference sequence database: {args.db_file}')
    check_files(args.db_file, args.input_file, args.output_dir)
    blast_results = align_frag_end_sequences(args.db_file, args.input_file, args.output_dir,
                                             args.output_name, args.blast_threads)
    blast_results = filter_alignments_by_id_and_cov(blast_results, args.min_ID, args.min_cov)
    blast_results = filter_alignments_for_uniq_aligners(blast_results)
    blast_results = annotate_frag_end_seqs(blast_results)
    blast_results = filter_alignments_for_both_ends_aligned(blast_results)
    blast_results = filter_alignments_plus_minus_strands(blast_results)
    blast_results = get_best_ref_seqs_matches(blast_results, args.num_ranks)
    blast_results = annotate_ref_seqs(blast_results)
    blast_results = filter_alignments_for_uniq_segments_subtypes(blast_results)
    merge_frag_ends(blast_results, args.input_file, args.output_dir, args.output_name)
    blast_results = align_merged_frag_sequences(args.db_file, args.output_dir, args.output_name,
                                                args.blast_threads)
    blast_results = annotate_merged_frag_seqs(blast_results)
    blast_results = get_best_ref_seqs_matches(blast_results, args.num_ranks)
    blast_results = annotate_ref_seqs(blast_results)
    blast_results = filter_alignments_for_uniq_segments_subtypes(blast_results)
    blast_results = compare_frag_and_ref_seqs(blast_results, args.output_dir, args.output_name)
    write_reports(blast_results, args.output_dir, args.output_name)
    print('\nDone.\n')


def get_args():
    parser = arg.ArgumentParser()
    parser.add_argument('-i', '--input_file', required=True, type=str)
    parser.add_argument('-o', '--output_dir', required=True, type=str)
    parser.add_argument('-n', '--output_name', required=True, type=str)
    parser.add_argument('-d', '--db_file', required=True, type=str)
    parser.add_argument('-I', '--min_ID', required=False, default=90, type=float)
    parser.add_argument('-c', '--min_cov', required=False, default=95, type=float)
    parser.add_argument('-r', '--num_ranks', required=False, default=2, type=int)
    parser.add_argument('-t', '--blast_threads', required=False, default=1, type=int)
    args = parser.parse_args()
    return args


def check_files(db_file_path, query_file_path, output_dir):
    for file in [db_file_path, query_file_path]:
        if os.path.isfile(file) == False:
            print(f'\nERROR: {file} does not exist!\n')
            exit()
    if os.path.exists(output_dir) == False:
        print(f'Creating output directory {output_dir}...')
        os.mkdir(output_dir)
    else:
        if os.path.isdir(output_dir) == False:
            print(f'\nERROR: {output_dir} exists but is not a directory!\n')
            exit()


def align_frag_end_sequences(db_file_path, query_file_path, output_dir, output_name, blast_threads):
    suffixes = 'nhr nin nsq'.split(' ')
    if any(os.path.isfile(f'{db_file_path}.{suffix}') == False for suffix in suffixes):
        print(f'Making blastn database for {db_file_path}...')
        terminal_command = f'makeblastdb -in {db_file_path} -dbtype nucl'
        sp.run(terminal_command, shell=True, stderr=sp.DEVNULL, stdout=sp.DEVNULL)
    print('Aligning fragment end sequences to reference database...')
    blast_results_path = os.path.join(output_dir, f'{output_name}_blast_results.tsv')
    cols = 'qseqid sseqid pident qcovhsp bitscore qstart qend qlen sstart send slen qseq sseq sstrand'
    terminal_command = f'blastn -db {db_file_path} -query {query_file_path}'
    terminal_command += f' -outfmt "6 {cols}" -num_threads {blast_threads} > {blast_results_path}'
    sp.run(terminal_command, shell=True)
    cols = cols.split(' ')
    blast_results = pd.read_csv(blast_results_path, sep='\t', names=cols)
    print(f' {len(blast_results)} alignments')
    if len(blast_results) == 0:
        print(f'\nNo alignments.\nDone.\n')
        exit()
    blast_results = blast_results.drop_duplicates()
    os.remove(blast_results_path)
    return blast_results.drop_duplicates()


def filter_alignments_by_id_and_cov(blast_results, min_id, min_cov):
    print('Filtering fragment end alignments by nucleotide identity and query coverage...')
    blast_results = blast_results[blast_results['pident'] >= min_id]
    blast_results = blast_results[blast_results['qcovhsp'] >= min_cov]
    print(f' {len(blast_results)} alignments')
    if len(blast_results) == 0:
        print(f'\nNo alignments.\nDone.\n')
        exit()
    return blast_results.drop_duplicates()


def filter_alignments_for_uniq_aligners(blast_results):  
    print('Discarding alignments between fragment ends and reference sequences involving fragment ends that aligned to multiple locations in those reference sequences...')
    cols = 'qseqid sseqid qseq'.split(' ')
    group_cols = 'qseqid sseqid'.split(' ')
    uniq_aligned_queries = blast_results[cols].drop_duplicates().groupby(group_cols).count().reset_index()
    uniq_aligned_queries = uniq_aligned_queries[uniq_aligned_queries['qseq']==1]
    merge_cols = 'qseqid sseqid'.split(' ')
    blast_results = pd.merge(blast_results, uniq_aligned_queries[merge_cols], on=merge_cols)
    print(f' {len(blast_results)} alignments')
    if len(blast_results) == 0:
        print(f'\nNo alignments.\nDone.\n')
        exit()
    return blast_results.drop_duplicates()


def annotate_frag_end_seqs(blast_results):
    print('Annotating fragment sequences...')
    q_annots = blast_results[['qseqid']].drop_duplicates()
    fields = 'experiment lib_name frag_name frag_end frag_copies UMI_pair'.split(' ')
    for i, name in enumerate(fields):
        q_annots[name] = q_annots.apply(lambda row: row['qseqid'].split('|')[i], axis=1)
    blast_results = pd.merge(blast_results, q_annots, on='qseqid')
    print(f' {len(blast_results)} alignments')
    cols = 'experiment lib_name frag_name UMI_pair'.split(' ')
    print(f' {len(blast_results[cols].drop_duplicates())} fragments')
    if len(blast_results) == 0:
        print(f'\nNo alignments.\nDone.\n')
        exit()
    return blast_results.drop_duplicates()


def filter_alignments_for_both_ends_aligned(blast_results):   
    print('Discarding alignments between fragments and reference sequences if both ends of those fragments did not'
          ' align to those reference sequence...')    
    cols = 'experiment lib_name frag_name UMI_pair sseqid frag_end'.split(' ')
    group_cols = 'experiment lib_name frag_name UMI_pair sseqid'.split(' ')
    frag_ends_aligned = blast_results[cols].drop_duplicates().groupby(group_cols).count().reset_index()
    frag_ends_aligned = frag_ends_aligned[frag_ends_aligned['frag_end']==2]
    merge_cols = 'experiment lib_name frag_name UMI_pair sseqid'.split(' ')
    blast_results = pd.merge(blast_results, frag_ends_aligned[merge_cols], on=merge_cols)
    print(f' {len(blast_results)} alignments')
    cols = 'experiment lib_name frag_name UMI_pair'.split(' ')
    print(f' {len(blast_results[cols].drop_duplicates())} fragments')
    if len(blast_results) == 0:
        print(f'\nNo alignments.\nDone.\n')
        exit()
    return blast_results.drop_duplicates()


def filter_alignments_plus_minus_strands(blast_results):
    print('Discarding alignments between fragments and reference sequences if those fragments did not have exactly one plus sense alignment to those reference sequences...')
    cols = ['experiment', 'lib_name', 'frag_name', 'UMI_pair', 'sseqid', 'sseq']
    group_cols = ['experiment', 'lib_name', 'frag_name', 'UMI_pair', 'sseqid']
    num_plus_strands = blast_results[blast_results['sstrand']=='plus'][cols].drop_duplicates().groupby(group_cols).count().reset_index()
    num_plus_strands = num_plus_strands[num_plus_strands['sseq']==1]
    merge_cols = ['experiment', 'lib_name', 'frag_name', 'UMI_pair', 'sseqid']
    blast_results = pd.merge(blast_results, num_plus_strands[merge_cols], on=merge_cols)
    print(f' {len(blast_results)} alignments')
    cols = 'experiment lib_name frag_name UMI_pair'.split(' ')
    print(f' {len(blast_results[cols].drop_duplicates())} fragments')
    if len(blast_results) == 0:
        print(f'\nNo alignments.\nDone.\n')
        exit()
    print('Discarding alignments between fragments and reference sequences if those fragments did not have exactly one minus sense alignment to those reference sequences...')
    cols = ['experiment', 'lib_name', 'frag_name', 'UMI_pair', 'sseqid', 'sseq']
    group_cols = ['experiment', 'lib_name', 'frag_name', 'UMI_pair', 'sseqid']
    num_minus_strands = blast_results[blast_results['sstrand']=='minus'][cols].drop_duplicates().groupby(group_cols).count().reset_index()
    num_minus_strands = num_minus_strands[num_minus_strands['sseq']==1]
    merge_cols = ['experiment', 'lib_name', 'frag_name', 'UMI_pair', 'sseqid']
    blast_results = pd.merge(blast_results, num_minus_strands[merge_cols], on=merge_cols)
    print(f' {len(blast_results)} alignments')
    cols = 'experiment lib_name frag_name UMI_pair'.split(' ')
    print(f' {len(blast_results[cols].drop_duplicates())} fragments')
    if len(blast_results) == 0:
        print(f'\nNo alignments.\nDone.\n')
        exit()
    print('Discarding alignments between fragments and reference sequences if those fragments did not have both plus and minus sense alignments to those reference sequences...')
    cols = ['experiment', 'lib_name', 'frag_name', 'UMI_pair', 'sseqid', 'sstrand']
    group_cols = ['experiment', 'lib_name', 'frag_name', 'UMI_pair', 'sseqid']
    num_strands = blast_results[cols].drop_duplicates().groupby(group_cols).count().reset_index()
    num_strands = num_strands[num_strands['sstrand']==2]
    merge_cols = ['experiment', 'lib_name', 'frag_name', 'UMI_pair', 'sseqid']
    blast_results = pd.merge(blast_results, num_minus_strands[merge_cols], on=merge_cols)
    print(f' {len(blast_results)} alignments')
    cols = 'experiment lib_name frag_name UMI_pair'.split(' ')
    print(f' {len(blast_results[cols].drop_duplicates())} fragments')
    if len(blast_results) == 0:
        print(f'\nNo alignments.\nDone.\n')
        exit()
    return blast_results


def get_best_ref_seqs_matches(blast_results, num_ranks):   
    print('Calcuating combined bitscores for each fragment against each reference sequences...')  
    cols = 'experiment lib_name frag_name UMI_pair sseqid bitscore'.split(' ')
    group_cols = 'experiment lib_name frag_name UMI_pair sseqid'.split(' ')
    combined_bitscores = blast_results[cols].groupby(group_cols).sum().reset_index()
    print(f'Identifying the top {num_ranks} bitscores for each fragment...')
    cols = 'experiment lib_name frag_name UMI_pair bitscore'.split(' ')
    group_cols = 'experiment lib_name frag_name UMI_pair'.split(' ')
    top_combined_bitscores = combined_bitscores[cols].groupby(group_cols)['bitscore'].nlargest(num_ranks, keep='all').reset_index()
    merge_cols = 'experiment lib_name frag_name UMI_pair bitscore'.split(' ')
    top_combined_bitscores = pd.merge(top_combined_bitscores, combined_bitscores, on=merge_cols)
    print('Restricting fragment end alignments to reference sequences with top combined bitscores...')
    merge_cols = 'experiment lib_name frag_name UMI_pair sseqid'.split(' ')
    blast_results = pd.merge(blast_results, top_combined_bitscores[merge_cols], on=merge_cols)
    print(f' {len(blast_results)} alignments')
    cols = 'experiment lib_name frag_name UMI_pair'.split(' ')
    print(f' {len(blast_results[cols].drop_duplicates())} fragments')
    if len(blast_results) == 0:
        print(f'\nNo alignments.\nDone.\n')
        exit()
    return blast_results.drop_duplicates()


def annotate_ref_seqs(blast_results):
    print('Annotating reference sequences...')
    s_annots = blast_results[['sseqid', 'slen']].drop_duplicates()
    s_annots.columns = ['sseqid', 'ref_seq_length']
    fields = 'ref_seq_accession ref_seq_name segment subtype'.split(' ')
    for i, name in enumerate(fields):
        s_annots[name] = blast_results.apply(lambda row: row['sseqid'].split('|')[i], axis=1)
    blast_results = pd.merge(blast_results, s_annots, on='sseqid')
    print(f' {len(blast_results)} alignments')
    cols = 'experiment lib_name frag_name UMI_pair'.split(' ')
    print(f' {len(blast_results[cols].drop_duplicates())} fragments')
    if len(blast_results) == 0:
        print(f'\nNo alignments.\nDone.\n')
        exit()
    return blast_results.drop_duplicates()


def filter_alignments_for_uniq_segments_subtypes(blast_results):
    print("Ensuring each fragment's best reference sequence matches are the same segment...")
    cols = 'experiment lib_name frag_name UMI_pair segment'.split(' ')
    group_cols = 'experiment lib_name frag_name UMI_pair'.split(' ')
    num_segments = blast_results[cols].drop_duplicates().groupby(group_cols).count().reset_index()
    num_segments = num_segments[num_segments['segment']==1]
    merge_cols = 'experiment lib_name frag_name UMI_pair'.split(' ')
    blast_results = pd.merge(blast_results, num_segments[merge_cols], on=merge_cols)
    print(f' {len(blast_results)} alignments')
    cols = 'experiment lib_name frag_name UMI_pair'.split(' ')
    print(f' {len(blast_results[cols].drop_duplicates())} fragments')
    if len(blast_results) == 0:
        print(f'\nNo alignments.\nDone.\n')
        exit()  
    print("Ensuring each fragment's best reference sequence matches are the same subtype...")
    cols = 'experiment lib_name frag_name UMI_pair subtype'.split(' ')
    group_cols = 'experiment lib_name frag_name UMI_pair'.split(' ')
    num_subtypes = blast_results[cols].drop_duplicates().groupby(group_cols).count().reset_index()
    num_subtypes = num_subtypes[num_subtypes['subtype']==1]
    merge_cols = 'experiment lib_name frag_name UMI_pair'.split(' ')
    blast_results = pd.merge(blast_results, num_subtypes[merge_cols], on=merge_cols)
    print(f' {len(blast_results)} alignments')
    cols = 'experiment lib_name frag_name UMI_pair'.split(' ')
    print(f' {len(blast_results[cols].drop_duplicates())} fragments')
    if len(blast_results) == 0:
        print(f'\nNo alignments.\nDone.\n')
        exit()
    return blast_results.drop_duplicates()


def merge_frag_ends(blast_results, query_file_path, output_dir, output_name):
    print('Merging fragment ends into fragment sequences...')
    # Analyze plus strand
    cols = ['experiment', 'lib_name', 'frag_name', 'UMI_pair', 'qseqid', 'frag_copies', 'sseqid', 
            'segment', 'subtype', 'qstart', 'qend', 'qlen', 'sstart', 'send', 'qseq']
    plus_strand = blast_results[blast_results['sstrand']=='plus'][cols].drop_duplicates()
    plus_strand.columns = ['experiment', 'lib_name', 'frag_name', 'UMI_pair', 'plus_strand_qseqid', 
                           'frag_copies', 'sseqid', 'segment', 'subtype', 'qstart', 'qend', 'qlen', 
                           'plus_strand_start', 'plus_strand_end', 'plus_strand_qseq']
    plus_strand['plus_strand_head'] = plus_strand['qstart'] - 1
    plus_strand['plus_strand_tail'] = plus_strand['qlen'] - plus_strand['qend']
    cols = ['experiment', 'lib_name', 'frag_name', 'UMI_pair', 'plus_strand_qseqid', 'frag_copies', 
            'sseqid', 'segment', 'subtype', 'plus_strand_head', 'plus_strand_tail', 'plus_strand_start', 
            'plus_strand_end', 'plus_strand_qseq']
    plus_strand = plus_strand[cols]
    # Analyze minus strand
    cols = ['experiment', 'lib_name', 'frag_name', 'UMI_pair', 'qseqid', 'frag_copies', 'sseqid', 
            'segment', 'subtype', 'qstart', 'qend', 'qlen', 'send', 'sstart', 'qseq']
    minus_strand = blast_results[blast_results['sstrand']=='minus'][cols].drop_duplicates()
    minus_strand.columns = ['experiment', 'lib_name', 'frag_name', 'UMI_pair', 'minus_strand_qseqid', 
                            'frag_copies', 'sseqid', 'segment', 'subtype', 'qstart', 'qend', 'qlen', 
                            'minus_strand_start', 'minus_strand_end', 'minus_strand_qseq']
    minus_strand['minus_strand_head'] = minus_strand['qstart'] - 1
    minus_strand['minus_strand_tail'] = minus_strand['qlen'] - minus_strand['qend']
    rev_comp = lambda seq: ''.join({'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N', '-': '-'}[base] for base in seq[::-1])
    minus_strand['minus_strand_qseq'] = minus_strand.apply(lambda row: rev_comp(row['minus_strand_qseq']), axis=1)
    cols = ['experiment', 'lib_name', 'frag_name', 'UMI_pair', 'minus_strand_qseqid', 'frag_copies', 
            'sseqid', 'segment', 'subtype', 'minus_strand_head', 'minus_strand_tail', 'minus_strand_start', 
            'minus_strand_end', 'minus_strand_qseq']
    minus_strand = minus_strand[cols]
    # Merge plus and minus strand analysis
    merge_cols = ['experiment', 'lib_name', 'frag_name', 'UMI_pair', 'frag_copies', 'sseqid', 'segment', 'subtype']
    merged_strands = pd.merge(plus_strand, minus_strand, on=merge_cols).drop_duplicates()
    # Estimate potential lengths of each fragment based on alignment coordinates to each top-matching reference sequence
    def get_frag_length(row):
        coordinates = sorted((row['plus_strand_start'], row['plus_strand_end'], row['minus_strand_start'], row['minus_strand_end']))
        frag_start, frag_end = coordinates[0], coordinates[-1]
        frag_length = frag_end - frag_start + 1
        return frag_length
    merged_strands['frag_length'] = merged_strands.apply(get_frag_length, axis=1)
    # Get median potential length of each fragment
    cols = 'experiment lib_name frag_name UMI_pair frag_length'.split(' ')
    group_cols = 'experiment lib_name frag_name UMI_pair'.split(' ')
    median_frag_length = merged_strands[cols].groupby(group_cols).quantile(0.5, interpolation='nearest').reset_index()
    median_frag_length['frag_length'] = median_frag_length['frag_length'].astype(int)
    # Keep only fragment end pairings forming fragments of the median length 
    merge_cols = 'experiment lib_name frag_name UMI_pair frag_length'.split(' ')
    merged_strands = pd.merge(merged_strands, median_frag_length, on=merge_cols).drop_duplicates()
    # Load query seqs
    query_seqs = {}
    with open(query_file_path, 'r') as input_file:
        for line in input_file:
            if line[0] == '>':
                header = line.strip().lstrip('>')
                query_seqs[header] = ''
            else:
                query_seqs[header] += line.strip().upper()
    # Merge fragment end sequences based on 'meet in the middle' algorithm
    def meet_in_the_middle(row):
        fwd_seq = ''
        rev_seq = ''
        plus_strand_seq = row['plus_strand_qseq']
        minus_strand_seq = rev_comp(row['minus_strand_qseq'])
        i = 0
        while (i * 2) + 1 < row['frag_length']:
            if i < len(plus_strand_seq) and i < len(minus_strand_seq):
                fwd_seq += plus_strand_seq[i]
                rev_seq += minus_strand_seq[i]
            elif i < len(plus_strand_seq) and i >= len(minus_strand_seq):
                fwd_seq += plus_strand_seq[i]
            elif i >= len(plus_strand_seq) and i < len(minus_strand_seq):
                rev_seq += minus_strand_seq[i]
            elif i >= len(plus_strand_seq) and i >= len(minus_strand_seq):        
                fwd_seq += 'N'
                rev_seq += 'N'
            i += 1
        if (i * 2) < row['frag_length']:
            if i < len(plus_strand_seq):
                fwd_seq += plus_strand_seq[i]
            elif i < len(minus_strand_seq):
                rev_seq += minus_strand_seq[i]
            else:
                fwd_seq += 'N'
        fwd_seq = query_seqs[row['plus_strand_qseqid']][:row['plus_strand_head']] + fwd_seq
        rev_seq = query_seqs[row['minus_strand_qseqid']][:row['minus_strand_head']] + rev_seq
        seq = fwd_seq + rev_comp(rev_seq)
        seq = seq.replace('-', '')
        return seq
    merged_strands['merged_frag_seq'] = merged_strands.apply(meet_in_the_middle, axis=1)
    cols = 'experiment lib_name frag_name UMI_pair frag_copies segment subtype merged_frag_seq'.split(' ')
    merged_strands = merged_strands[cols].drop_duplicates()
    # If multiple sequences exist for fragments, choose the one(s) with fewest Ns
    count_Ns = lambda row: row['merged_frag_seq'].count('N')
    merged_strands['merged_frag_seq_Ns'] = merged_strands.apply(count_Ns, axis=1)
    cols = 'experiment lib_name frag_name UMI_pair merged_frag_seq_Ns'.split(' ')
    group_cols = 'experiment lib_name frag_name UMI_pair'.split(' ')
    min_Ns = merged_strands[cols].drop_duplicates().groupby(group_cols).min().reset_index()
    merge_cols = 'experiment lib_name frag_name UMI_pair merged_frag_seq_Ns'.split(' ')
    merged_strands = pd.merge(merged_strands, min_Ns, on=merge_cols).drop_duplicates()
    # If multiple sequences exist for fragments, choose the one that appears first alphabetically
    cols = 'experiment lib_name frag_name UMI_pair merged_frag_seq'.split(' ')
    group_cols = 'experiment lib_name frag_name UMI_pair'.split(' ')
    first_seqs = merged_strands[cols].drop_duplicates().groupby(group_cols).min().reset_index()
    merge_cols = 'experiment lib_name frag_name UMI_pair merged_frag_seq'.split(' ')
    merged_strands = pd.merge(merged_strands, first_seqs, on=merge_cols).drop_duplicates()
    # Write out merged fragment sequences
    frag_seqs_path = os.path.join(output_dir, f'{output_name}_flu_frag_seqs.fa')
    with open(frag_seqs_path, 'w') as output_file:
        for index, row in merged_strands.iterrows():
            seq = row['merged_frag_seq']
            header = (row['experiment'], row['lib_name'], row['frag_name'], f'{len(seq)}_bases',
                      row['frag_copies'], row['UMI_pair'], row['segment'], row['subtype'])
            header = '|'.join(header)
            header = '>' + header + '|'
            output_file.write(header + '\n')
            output_file.write(seq + '\n')
    print(f' Wrote {len(merged_strands)} merged fragment sequences')
    if len(merged_strands) == 0:
        print(f'\nNo fragments.\nDone.\n')
        exit()


def align_merged_frag_sequences(db_file_path, output_dir, output_name, blast_threads):
    suffixes = 'nhr nin nsq'.split(' ')
    if any(os.path.isfile(f'{db_file_path}.{suffix}') == False for suffix in suffixes):
        print(f'Making blastn database for {db_file_path}...')
        terminal_command = f'makeblastdb -in {db_file_path} -dbtype nucl'
        sp.run(terminal_command, shell=True, stderr=sp.DEVNULL, stdout=sp.DEVNULL)
    print('Aligning full fragment sequences to reference database...')
    query_file_path = os.path.join(output_dir, f'{output_name}_flu_frag_seqs.fa')
    blast_results_path = os.path.join(output_dir, f'{output_name}_blast_results.tsv')
    cols = 'qseqid sseqid pident qcovhsp bitscore sstart send qseq sseq slen'
    terminal_command = f'blastn -db {db_file_path} -query {query_file_path}'
    terminal_command += f' -outfmt "6 {cols}" -num_threads {blast_threads} > {blast_results_path}'
    sp.run(terminal_command, shell=True)
    cols = cols.split(' ')
    blast_results = pd.read_csv(blast_results_path, sep='\t', names=cols)
    print(f' {len(blast_results)} alignments')
    if len(blast_results) == 0:
        print(f'\nNo alignments.\nDone.\n')
        exit()
    blast_results = blast_results.drop_duplicates()
    os.remove(blast_results_path)
    return blast_results.drop_duplicates()


def annotate_merged_frag_seqs(blast_results):
    print('Annotating merged fragment sequence alignments...')
    q_annots = blast_results[['qseqid']].drop_duplicates()
    fields = 'experiment lib_name frag_name frag_length frag_copies UMI_pair'.split(' ')
    for i, name in enumerate(fields):
        q_annots[name] = q_annots.apply(lambda row: row['qseqid'].split('|')[i], axis=1)
    q_annots['frag_length'] = q_annots.apply(lambda row: int(row['frag_length'].split('_')[0]), axis=1)
    q_annots['frag_copies'] = q_annots.apply(lambda row: int(row['frag_copies'].split('_')[0]), axis=1)
    blast_results = pd.merge(blast_results, q_annots, on='qseqid')
    print(f' {len(blast_results)} alignments')
    cols = 'experiment lib_name frag_name UMI_pair'.split(' ')
    print(f' {len(blast_results[cols].drop_duplicates())} fragments')
    if len(blast_results) == 0:
        print(f'\nNo alignments.\nDone.\n')
        exit()
    return blast_results.drop_duplicates()


def compare_frag_and_ref_seqs(blast_results, output_dir, output_name):
    #
    frag_seqs = {}
    frag_seqs_path = os.path.join(output_dir, f'{output_name}_flu_frag_seqs.fa')
    with open(frag_seqs_path, 'r') as input_file:
        for line in input_file:
            if line[0] == '>':
                header = line.strip().lstrip('>')
                frag_seqs[header] = ''
            else:
                frag_seqs[header] += line.strip()
    count_sequenced_bases = lambda row: sum(frag_seqs[row['qseqid']].count(base) for base in 'ATGC')
    get_perc_sequenced = lambda row: count_sequenced_bases(row) * 100 / len(frag_seqs[row['qseqid']])
    blast_results['perc_sequenced'] = blast_results.apply(get_perc_sequenced, axis=1)
    #
    def get_aligned(data_frame):
        aligned_positions = set()
        for index, row in data_frame.iterrows():
            aligned_positions = aligned_positions.union(set(range(row['sstart'], row['send']+1)))
        return len(aligned_positions)
    cols = ['qseqid', 'sseqid', 'sstart', 'send']
    group_cols = ['qseqid', 'sseqid']
    aligned = blast_results[cols].drop_duplicates().groupby(group_cols)[cols].apply(get_aligned).reset_index()
    aligned.columns = ['qseqid', 'sseqid', 'aligned']
    blast_results = pd.merge(blast_results, aligned, on=['qseqid', 'sseqid'])
    #
    def get_identical(data_frame):
        identical_positions = set()
        for index, row in data_frame.iterrows():
            position = row['sstart']
            for qbase, sbase in zip(row['qseq'], row['sseq']):
                if qbase not in 'N-' and qbase == sbase:
                    identical_positions.add(position)
                if sbase != '-':
                    position += 1
        return len(identical_positions)
    cols = ['qseqid', 'sseqid', 'sstart', 'qseq', 'sseq']
    group_cols = ['qseqid', 'sseqid']
    identical = blast_results[cols].drop_duplicates().groupby(group_cols)[cols].apply(get_identical).reset_index()
    identical.columns = ['qseqid', 'sseqid', 'identical']
    blast_results = pd.merge(blast_results, identical, on=['qseqid', 'sseqid'])
    #
    def get_mismatches(data_frame):
        mismatched_positions = set()
        for index, row in data_frame.iterrows():
            position = row['sstart']
            for qbase, sbase in zip(row['qseq'], row['sseq']):
                if qbase not in 'N-' and sbase not in 'N-' and qbase != sbase:
                    mismatched_positions.add(position)
                if sbase != '-':
                    position += 1
        return len(mismatched_positions)
    cols = ['qseqid', 'sseqid', 'sstart', 'qseq', 'sseq']
    group_cols = ['qseqid', 'sseqid']
    mismatches = blast_results[cols].drop_duplicates().groupby(group_cols)[cols].apply(get_mismatches).reset_index()
    mismatches.columns = ['qseqid', 'sseqid', 'mismatches']
    blast_results = pd.merge(blast_results, mismatches, on=['qseqid', 'sseqid'])
    #
    def get_deletions(data_frame):
        deleted_positions = set()
        for index, row in data_frame.iterrows():
            position = row['sstart']
            last_qbase = ''
            for qbase, sbase in zip(row['qseq'], row['sseq']):
                if qbase == '-' and last_qbase != '-':
                    deleted_positions.add(position)
                if sbase != '-':
                    position += 1
                last_qbase = qbase
        return len(deleted_positions)
    cols = ['qseqid', 'sseqid', 'sstart', 'qseq', 'sseq']
    group_cols = ['qseqid', 'sseqid']
    deletions = blast_results[cols].drop_duplicates().groupby(group_cols)[cols].apply(get_deletions).reset_index()
    deletions.columns = ['qseqid', 'sseqid', 'deletions']
    blast_results = pd.merge(blast_results, deletions, on=['qseqid', 'sseqid'])
    #
    def get_insertions(data_frame):
        inserted_positions = set()
        for index, row in data_frame.iterrows():
            position = row['sstart']
            last_sbase = ''
            for qbase, sbase in zip(row['qseq'], row['sseq']):
                if sbase == '-' and last_sbase != '-':
                    inserted_positions.add(position)
                if sbase != '-':
                    position += 1
                last_sbase = sbase
        return len(inserted_positions)
    cols = ['qseqid', 'sseqid', 'sstart', 'qseq', 'sseq']
    group_cols = ['qseqid', 'sseqid']
    insertions = blast_results[cols].drop_duplicates().groupby(group_cols)[cols].apply(get_insertions).reset_index()
    insertions.columns = ['qseqid', 'sseqid', 'insertions']
    blast_results = pd.merge(blast_results, insertions, on=['qseqid', 'sseqid'])
    #
    def get_ambiguous(data_frame):
        ambiguous_positions = set()
        for index, row in data_frame.iterrows():
            position = row['sstart']
            for qbase, sbase in zip(row['qseq'], row['sseq']):
                if (qbase == 'N' and sbase in 'ATGCN') or (sbase == 'N' and qbase in 'ATGCN'):
                    ambiguous_positions.add(position)
                if sbase != '-':
                    position += 1
        return len(ambiguous_positions)
    cols = ['qseqid', 'sseqid', 'sstart', 'qseq', 'sseq']
    group_cols = ['qseqid', 'sseqid']
    ambiguous = blast_results[cols].drop_duplicates().groupby(group_cols)[cols].apply(get_ambiguous).reset_index()
    ambiguous.columns = ['qseqid', 'sseqid', 'ambiguous']
    blast_results = pd.merge(blast_results, ambiguous, on=['qseqid', 'sseqid'])
    #
    blast_results['perc_identity'] = round(blast_results['identical'] * 100 / blast_results['aligned'], 1)
    blast_results['perc_ambiguous'] = round(blast_results['ambiguous'] * 100 / blast_results['aligned'], 1)
    blast_results['perc_aligned'] = round(blast_results['aligned'] * 100 / blast_results['frag_length'], 1)
    blast_results['ref_seq_cov_sequenced'] = round(blast_results['aligned'] * 100 / blast_results['ref_seq_length'], 1)
    blast_results['ref_seq_cov_total'] = round(blast_results['frag_length'] * 100 / blast_results['ref_seq_length'], 1)
    return blast_results


def write_reports(blast_results, output_dir, output_name):
    print('Writing reports...')
    # Top ref seqs report
    cols = ['experiment', 'lib_name', 'frag_name', 'UMI_pair', 'frag_copies', 'frag_length', 'perc_sequenced',
            'ref_seq_name', 'ref_seq_length', 'segment', 'subtype', 'aligned', 'identical', 'mismatches', 'deletions',
            'insertions', 'ambiguous', 'perc_aligned', 'perc_identity', 'perc_ambiguous', 'ref_seq_cov_sequenced',
            'ref_seq_cov_total']
    report = blast_results[cols].drop_duplicates()
    report['frag_number'] = report.apply(lambda row: int(row['frag_name'].split('_')[1]), axis=1)
    sort_cols = 'experiment lib_name frag_number ref_seq_name'.split(' ')
    report_path = os.path.join(output_dir, f'{output_name}_top_ref_seqs_report.csv')
    report.sort_values(by=sort_cols)[cols].drop_duplicates().to_csv(report_path, index=False)
    # Frag report
    cols = ['experiment', 'lib_name', 'frag_name', 'UMI_pair', 'frag_copies', 'segment', 'subtype', 'frag_length',
            'perc_sequenced', 'perc_aligned', 'perc_identity', 'perc_ambiguous', 'ref_seq_cov_sequenced',
            'ref_seq_cov_total']
    group_cols = ['experiment', 'lib_name', 'frag_name', 'UMI_pair', 'frag_copies', 'segment', 'subtype',
                  'frag_length', 'perc_sequenced']
    report = blast_results[cols].groupby(group_cols).median().reset_index()
    for col in ['perc_sequenced', 'perc_aligned', 'perc_identity', 'perc_ambiguous', 'ref_seq_cov_sequenced',
                'ref_seq_cov_total']:
        report[col] = round(report[col], 1)
    report['frag_number'] = report.apply(lambda row: int(row['frag_name'].split('_')[1]), axis=1)
    sort_cols = 'experiment lib_name frag_number'.split(' ')
    report_path = os.path.join(output_dir, f'{output_name}_frag_report.csv')
    report.sort_values(by=sort_cols)[cols].drop_duplicates().to_csv(report_path, index=False)


if __name__ == '__main__':
    main()
