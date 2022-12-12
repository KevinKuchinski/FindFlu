import argparse as arg
import numpy as np
import pandas as pd
import subprocess as sp
import os as os


def main():
    version = '0.0.8'
    print(f'\nFindFlu v{version}')
    print(f'https://github.com/KevinKuchinski/FindFlu/\n')
    args = get_args()
    print(f'\nReference sequence database: {args.db_file}\n')
    blast_results = blast_sequences(args.db_file, args.input_file, args.output_dir, args.output_name, args.blast_threads)
    delete_blast_results_file(args.output_dir, args.output_name)
    blast_results = filter_frag_ends_by_id_and_cov(blast_results, args.min_ID, args.min_cov)
    blast_results = filter_fragments_by_best_matches(blast_results)
    blast_results = annotate_blast_results(blast_results)
    blast_results = calc_combined_id_and_cov(blast_results)
    blast_results = calc_frag_length(blast_results)
    blast_results = calc_subject_coverage(blast_results)
    write_frag_seqs(args.input_file, blast_results, args.output_dir, args.output_name)
    write_best_ref_seqs_report(blast_results, args.output_dir, args.output_name)
    write_frag_report(blast_results, args.output_dir, args.output_name)
    write_country_host_clade_reports(blast_results, args.output_dir, args.output_name)
    print('\nDone.\n')


def get_args():
    parser = arg.ArgumentParser()
    parser.add_argument('-i', '--input_file', required=True, type=str)
    parser.add_argument('-o', '--output_dir', required=True, type=str)
    parser.add_argument('-n', '--output_name', required=True, type=str)
    parser.add_argument('-d', '--db_file', required=True, type=str)
    parser.add_argument('-I', '--min_ID', required=False, default=95, type=float)
    parser.add_argument('-c', '--min_cov', required=False, default=95, type=float)
    parser.add_argument('-t', '--blast_threads', required=False, default=1, type=int)
    args = parser.parse_args()
    return args


def blast_sequences(db_file_path, query_file_path, output_dir, output_name, blast_threads):
    for file in [db_file_path, query_file_path]:
        if os.path.isfile(file) == False:
            print(f'\nERROR: {file} does not exist!\n')
            exit()
    suffixes = ['ndb', 'nhr', 'nin', 'njs', 'not', 'nsq', 'ntf', 'nto']
    if any(os.path.isfile(f'{db_file_path}.{suffix}') == False for suffix in suffixes):
        print(f'Making blastn database for {db_file_path}...')
        terminal_command = f'makeblastdb -in {db_file_path} -dbtype nucl'
        sp.run(terminal_command, shell=True)
    print('Aligning fragment end sequences to reference database...')
    blast_results_path = os.path.join(output_dir, f'{output_name}_blast_results.tsv')
    cols = 'qseqid sseqid bitscore qlen qstart qend qseq slen sstart send sseq'
    terminal_command = f'blastn -db {db_file_path} -query {query_file_path}'
    terminal_command += f' -outfmt "6 {cols}" -num_threads {blast_threads} > {blast_results_path}'
    sp.run(terminal_command, shell=True)
    cols = cols.split(' ')
    blast_results = pd.read_csv(blast_results_path, sep='\t', names=cols)
    if len(blast_results) == 0:
        print(f'\nNo alignments.\nDone.\n')
        exit()
    blast_results = blast_results.drop_duplicates()
    return blast_results


def delete_blast_results_file(output_dir, output_name):
    blast_results_path = os.path.join(output_dir, f'{output_name}_blast_results.tsv')
    os.remove(blast_results_path)


def count_valid_positions(qseq, sseq):
    valid_positions = 0
    for qbase, sbase in zip(qseq, sseq):
        if {qbase, sbase} != {'N'}:
            valid_positions += 1
    return valid_positions


def count_identical_bases(qseq, sseq):
    identical_bases = 0
    for qbase, sbase in zip(qseq, sseq):
        if qbase == sbase and {qbase, sbase} != {'N'}:
            identical_bases += 1
    return identical_bases


def count_covered_query_bases(qstart, qend):
    qstart, qend = sorted([qstart, qend])
    covered_query_bases = qend - qstart + 1
    return covered_query_bases
    

def filter_frag_ends_by_id_and_cov(blast_results, min_id, min_cov):
    print('Filtering fragment end alignments by nucleotide identity and query coverage...')
    if len(blast_results) == 0:
        return pd.DataFrame()
    alignment_annots = blast_results[['qseqid', 'sseqid', 'qstart', 'qend', 'qseq', 'sseq']].drop_duplicates()
    alignment_annots['valid_positions'] = alignment_annots.apply(lambda row: count_valid_positions(row['qseq'], row['sseq']), axis=1)
    alignment_annots['identical_bases'] = alignment_annots.apply(lambda row: count_identical_bases(row['qseq'], row['sseq']), axis=1)
    alignment_annots['covered_query_bases'] = alignment_annots.apply(lambda row: count_covered_query_bases(row['qstart'], row['qend']), axis=1)
    cols = ['qseqid', 'sseqid', 'valid_positions', 'identical_bases', 'covered_query_bases']
    blast_results = pd.merge(blast_results, alignment_annots[cols], on=['qseqid', 'sseqid'])
    blast_results['frag_end_id'] = blast_results['identical_bases'] * 100 / blast_results['valid_positions']
    blast_results['frag_end_cov'] = blast_results['covered_query_bases'] * 100 / blast_results['qlen']
    blast_results = blast_results[blast_results['frag_end_id'] >= min_id]
    blast_results = blast_results[blast_results['frag_end_cov'] >= min_cov]
    return blast_results


def filter_fragments_by_best_matches(blast_results):
    print('Finding best alignments...')
    if len(blast_results) == 0:
        return pd.DataFrame()    
    query_annots = blast_results[['qseqid']].drop_duplicates()
    query_annots['experiment'] = query_annots.apply(lambda row: row['qseqid'].split('|')[0], axis=1)
    query_annots['library'] = query_annots.apply(lambda row: row['qseqid'].split('|')[1], axis=1)
    query_annots['fragment'] = query_annots.apply(lambda row: row['qseqid'].split('|')[2], axis=1)
    query_annots['frag_end'] = query_annots.apply(lambda row: row['qseqid'].split('|')[3], axis=1)
    blast_results = pd.merge(blast_results, query_annots, on='qseqid')
    #
    cols = ['experiment', 'library', 'fragment', 'frag_end', 'sseqid']
    group_cols = ['experiment', 'library', 'fragment', 'sseqid']
    num_frag_ends_aligned_to_subject = blast_results[cols].drop_duplicates().groupby(group_cols).count().reset_index()
    num_frag_ends_aligned_to_subject = num_frag_ends_aligned_to_subject[num_frag_ends_aligned_to_subject['frag_end']==2]
    merge_cols = ['experiment', 'library', 'fragment', 'sseqid']
    blast_results = pd.merge(blast_results, num_frag_ends_aligned_to_subject[merge_cols], on=merge_cols)
    #
    cols = ['experiment', 'library', 'fragment', 'frag_end', 'sseqid', 'bitscore']
    combined_bitscores = blast_results[cols].drop_duplicates()
    cols = ['experiment', 'library', 'fragment', 'sseqid', 'bitscore']
    group_cols = ['experiment', 'library', 'fragment', 'sseqid']
    combined_bitscores = combined_bitscores[cols].groupby(group_cols).sum().reset_index()
    cols = ['experiment', 'library', 'fragment', 'bitscore']
    group_cols = ['experiment', 'library', 'fragment']
    max_combined_bitscores = combined_bitscores[cols].drop_duplicates().groupby(group_cols).max().reset_index()
    merge_cols = ['experiment', 'library', 'fragment', 'bitscore']
    combined_bitscores = pd.merge(combined_bitscores, max_combined_bitscores[merge_cols], on=merge_cols)
    merge_cols = ['experiment', 'library', 'fragment', 'sseqid']
    blast_results = pd.merge(blast_results, combined_bitscores[merge_cols], on=merge_cols)
    #
    subject_annots = blast_results[['sseqid']].drop_duplicates()
    subject_annots['segment'] = subject_annots.apply(lambda row: row['sseqid'].split('|')[2], axis=1)
    subject_annots['subtype'] = subject_annots.apply(lambda row: row['sseqid'].split('|')[3], axis=1)
    blast_results = pd.merge(blast_results, subject_annots, on='sseqid')
    #
    cols = ['experiment', 'library', 'fragment', 'segment']
    group_cols = ['experiment', 'library', 'fragment']
    num_segments_per_frag = blast_results[cols].drop_duplicates().groupby(group_cols).count().reset_index()
    num_segments_per_frag = num_segments_per_frag[num_segments_per_frag['segment']==1]
    merge_cols = ['experiment', 'library', 'fragment']
    blast_results = pd.merge(blast_results, num_segments_per_frag[merge_cols], on=merge_cols)
    #
    cols = ['experiment', 'library', 'fragment', 'subtype']
    group_cols = ['experiment', 'library', 'fragment']
    num_subtypes_per_frag = blast_results[cols].drop_duplicates().groupby(group_cols).count().reset_index()
    num_subtypes_per_frag = num_subtypes_per_frag[num_subtypes_per_frag['subtype']==1]
    merge_cols = ['experiment', 'library', 'fragment']
    blast_results = pd.merge(blast_results, num_subtypes_per_frag[merge_cols], on=merge_cols)
    return blast_results


def annotate_blast_results(blast_results):
    print('Annotating alignments...')
    if len(blast_results) == 0:
        return pd.DataFrame()
    query_annots = blast_results[['qseqid']].drop_duplicates()
    query_annots['frag_copies'] = query_annots.apply(lambda row: int(row['qseqid'].split('|')[4].split('_')[0]), axis=1)
    query_annots['UMI_pair'] = query_annots.apply(lambda row: row['qseqid'].split('|')[5].split('UMI_pair_')[1], axis=1)
    blast_results = pd.merge(blast_results, query_annots, on='qseqid')
    subject_annots = blast_results[['sseqid']].drop_duplicates()
    subject_annots['ref_seq_ID'] = subject_annots.apply(lambda row: row['sseqid'].split('|')[0], axis=1)
    subject_annots['ref_seq'] = subject_annots.apply(lambda row: row['sseqid'].split('|')[1], axis=1)
    subject_annots['country'] = subject_annots.apply(lambda row: row['sseqid'].split('|')[4], axis=1)
    subject_annots['host'] = subject_annots.apply(lambda row: row['sseqid'].split('|')[5], axis=1)
    subject_annots['H5_clade'] = subject_annots.apply(lambda row: row['sseqid'].split('|')[6], axis=1)
    blast_results = pd.merge(blast_results, subject_annots, on='sseqid')
    return blast_results


def calc_combined_id_and_cov(blast_results):
    print('Calculating combined nucleotide identity and query coverage for each fragment...')
    if len(blast_results) == 0:
        return pd.DataFrame()
    #
    cols = ['experiment', 'library', 'fragment', 'sseqid', 'valid_positions', 'identical_bases', 'qlen', 'covered_query_bases']
    group_cols = ['experiment', 'library', 'fragment', 'sseqid']
    combined_id_and_cov = blast_results[cols].groupby(group_cols).sum().reset_index()
    combined_id_and_cov['combined_query_id'] = combined_id_and_cov['identical_bases'] * 100 / combined_id_and_cov['valid_positions']
    combined_id_and_cov['combined_query_cov'] = combined_id_and_cov['covered_query_bases'] * 100 / combined_id_and_cov['qlen']
    cols = ['experiment', 'library', 'fragment', 'sseqid', 'combined_query_id', 'combined_query_cov']
    merge_cols = ['experiment', 'library', 'fragment', 'sseqid']
    blast_results = pd.merge(blast_results, combined_id_and_cov[cols], on=merge_cols)
    return blast_results


def calc_frag_length(blast_results):
    print('Calculating fragment lengths...')
    if len(blast_results) == 0:
        return pd.DataFrame()
    cols = ['experiment', 'library', 'fragment', 'sseqid', 'sstart', 'send']
    alignment_coordinates = blast_results[cols].drop_duplicates()
    #
    index_cols = ['experiment', 'library', 'fragment', 'sseqid']
    melted_cols = ['sstart', 'send']
    alignment_coordinates = pd.melt(alignment_coordinates, index_cols, melted_cols, 'blast_col', 'coordinate')
    #
    cols = ['experiment', 'library', 'fragment', 'sseqid', 'coordinate']
    group_cols = ['experiment', 'library', 'fragment', 'sseqid']
    alignment_starts = alignment_coordinates[cols].groupby(group_cols).min().reset_index()
    alignment_starts.columns = group_cols + ['start_coordinate']
    #
    alignment_ends = alignment_coordinates[cols].groupby(group_cols).max().reset_index()
    alignment_ends.columns = group_cols + ['end_coordinate']
    #
    alignment_coordinates = pd.merge(alignment_starts, alignment_ends, on=group_cols)
    alignment_coordinates['frag_length'] = alignment_coordinates['end_coordinate'] - alignment_coordinates['start_coordinate']
    alignment_coordinates = alignment_coordinates[group_cols + ['frag_length']]
    #
    merge_cols = ['experiment', 'library', 'fragment', 'sseqid']
    blast_results = pd.merge(blast_results, alignment_coordinates, on=merge_cols)
    return blast_results
    

def calc_subject_coverage(blast_results):
    print('Calculating subject coverage...')
    if len(blast_results) == 0:
        return pd.DataFrame()
    cols = ['experiment', 'library', 'fragment', 'sseqid', 'frag_length', 'slen']
    subject_coverage = blast_results[cols].drop_duplicates()
    subject_coverage['frag_subject_cov'] = subject_coverage['frag_length'] * 100 / subject_coverage['slen']
    blast_results = pd.merge(blast_results, subject_coverage, on=cols)
    return blast_results


def write_frag_seqs(query_file_path, blast_results, output_dir, output_name):
    print('Writing flu frag seqs to FASTA...')
    query_seqs = {}
    with open(query_file_path, 'r') as input_file:
        for line in input_file:
            if line[0] == '>':
                header = line.strip()
                query_seqs[header] = ''
            else:
                query_seqs[header] += line.strip()
    flu_seqs = {}
    for index in blast_results.index:
        experiment = blast_results.loc[index]['experiment']
        library = blast_results.loc[index]['library']
        fragment = blast_results.loc[index]['fragment']
        frag_copies = blast_results.loc[index]['frag_copies']
        UMI_pair = blast_results.loc[index]['UMI_pair']
        segment = blast_results.loc[index]['segment']
        subtype = blast_results.loc[index]['subtype']
        for frag_end in ['end_1', 'end_2']:
            header = f'>{experiment}|{library}|{fragment}|{frag_end}|{frag_copies}_copies|UMI_pair_{UMI_pair}'
            seq = query_seqs[header]
            header += f'|{segment}|{subtype}'
            flu_seqs[header] = seq
    output_fasta_path = os.path.join(output_dir, f'{output_name}_flu_frag_seqs.fa')
    with open(output_fasta_path, 'w') as output_file:
        for header, seq in flu_seqs.items():
            output_file.write(header + '\n')
            output_file.write(seq + '\n')


def write_best_ref_seqs_report(blast_results, output_dir, output_name):
    print('Writing best ref seqs report...')
    blast_results['ref_seq_length'] = blast_results['slen']
    cols = ['experiment', 'library', 'fragment', 'UMI_pair', 'ref_seq_ID', 'ref_seq', 'segment', 'subtype', 'host', 'country',
            'H5_clade', 'combined_query_id', 'combined_query_cov', 'ref_seq_length', 'frag_length', 'frag_subject_cov']
    report = blast_results[cols].drop_duplicates()
    for col in ['combined_query_id', 'combined_query_cov', 'frag_subject_cov']:
        report[col] = round(report[col], 1)
    for col in ['frag_length', 'ref_seq_length']:
        report[col] = report[col].astype(int)
    for col in ['subtype', 'H5_clade']:
        report[col] = report[col].replace('none', np.nan)
    report = report.sort_values(by=['experiment', 'library', 'fragment', 'UMI_pair', 'ref_seq_ID'])
    report_path = os.path.join(output_dir, f'{output_name}_best_ref_seqs_report.csv')
    report.to_csv(report_path, index=False)


def write_frag_report(blast_results, output_dir, output_name):
    print('Writing frag report...')
    cols = ['experiment', 'library', 'fragment', 'UMI_pair', 'frag_copies', 'segment', 'subtype',
            'combined_query_id', 'combined_query_cov', 'frag_length', 'frag_subject_cov']
    group_cols = ['experiment', 'library', 'fragment', 'UMI_pair', 'frag_copies', 'segment', 'subtype']
    report = blast_results[cols].drop_duplicates().groupby(group_cols).quantile(0.5, interpolation='midpoint').reset_index()
    for col in ['combined_query_id', 'combined_query_cov', 'frag_subject_cov']:
        report[col] = round(report[col], 1)
    for col in ['frag_length', 'frag_copies']:
        report[col] = report[col].astype(int)
    for col in ['subtype']:
        report[col] = report[col].replace('none', np.nan)
    report = report.sort_values(by=['experiment', 'library', 'fragment', 'UMI_pair', 'segment', 'subtype'])
    report_path = os.path.join(output_dir, f'{output_name}_frag_report.csv')
    report.to_csv(report_path, index=False)


def write_country_host_clade_reports(blast_results, output_dir, output_name):
    for variable in ['country', 'host', 'H5_clade']:
        print(f'Writing {variable} report...')
        cols = ['experiment', 'library', 'fragment', 'UMI_pair', 'segment', 'subtype', variable]
        report = blast_results[cols]
        if variable == 'H5_clade':
            report = report[report['subtype']=='H5']
            report = report[report['H5_clade']!='NA']
        report = report.groupby(cols).size().reset_index()
        report.columns = cols + [f'{variable}_count']
        cols = ['experiment', 'library', 'fragment', 'UMI_pair', f'{variable}_count']
        total_best_alignments = report[cols].groupby(['experiment', 'library', 'fragment', 'UMI_pair']).sum().reset_index()
        total_best_alignments.columns = ['experiment', 'library', 'fragment', 'UMI_pair', 'total_frag_matches']
        report = pd.merge(report, total_best_alignments, on=['experiment', 'library', 'fragment', 'UMI_pair'])
        report[f'{variable}_perc_of_total_matches'] = round(report[f'{variable}_count'] * 100 / report['total_frag_matches'], 1)
        report['subtype'] = report['subtype'].replace('none', np.nan)
        sort_by = ['experiment', 'library', 'fragment', 'UMI_pair', 'segment', 'subtype', f'{variable}_perc_of_total_matches']
        report = report.sort_values(by=sort_by, ascending=[True, True, True, False, True, True, False])
        report_path = os.path.join(output_dir, f'{output_name}_{variable}_report.csv')
        report.to_csv(report_path, index=False)


if __name__ == '__main__':
    main()
