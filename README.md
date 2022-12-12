# FindFlu
FindFlu is a bioinformatic tool designed to identify and characterize fragments of influenza A virus (IAV) genome from user-provided sequences. For each IAV genome fragment identified, FindFlu will report the segment of the IAV genome from which it originated, the subtype of that segment (for fragments of haemagglutinin and neuraminidase segments), and the length of the fragment. FindFlu identifies and characterizes fragments by aligning them to a database of IAV reference sequences and finding each fragment’s best-matching reference sequence. These reference sequences are annotated with the host species from which they were collected as well as the country in which they were collected. FindFlu reports these annotations for assessing the potential host range and geographical origin of IAV fragments. Additionally, haemagglutinin segment reference sequences belonging to the H5 subtype are annotated with their lineage and clade. FindFlu reports these annotations for identifying fragments of virus from recognized high-pathogenicity lineages. 

## FindFlu input
As input, FindFlu expects a FASTA file containing pairs of consensus sequences describing both ends of distinct genomic fragments. The entries in this FASTA file must have headers that use the HopDropper (v1.0.0) format (https://github.com/KevinKuchinski/HopDropper).

## FindFlu database
FindFlu aligns input sequences to a database of annotated IAV reference sequences. Reference sequences must be in a FASTA file with the following header format:

```>ref_seq_identifier|strain_name(strain_subtype)|segment|subtype|country|host|H5_clade```

Ref_seq_identifier is a unique identifier for the reference sequence, typically its GenBank accession number. Strain_name is the name of the strain (e.g. A/goose/Guangdong/1/1996), and strain_subtype is the HxNx subtype of the strain (regardless of which IAV genome segment the reference sequence represents). Segment is the IAV genome segment the reference sequence represents, and it must be encoded as one of the following: PB2, PB1, PA, HA, NP, NA, M, or NS. Subtype is the HA subtype (for HA segment sequences) or the NA subtype (for NA segment sequences). For internal segments, ‘none’ should be entered in the subtype field. Host is the host species from which the reference sequence was collected. Country is the country in which the reference sequence was collected. H5_clade is the H5 lineage or clade of the reference sequence. For reference sequences that are not HA segments belonging to subtype H5, ‘none’ should be entered in this field.

## Alignment and filtering
FindFlu aligns input sequences to the reference sequence database with blastn (v2.13.0). The following steps are performed on the blastn alignments:
1.	A nucleotide identity percentage is calculated for each alignment. This is done by dividing the number of identical positions in the alignment by the number of valid positions in the alignment. Identical alignment positions are those where both the query and subject sequences contain the same non-ambiguous nucleotide. Valid alignment positions are those where both the query and subject sequences contain non-ambiguous nucleotide bases. Alignments are discarded if their nucleotide identity is below a minimum threshold (default is 95%).

2.	A query coverage percentage is calculated for each alignment. This is done by dividing the covered query positions by the length of the query sequence. Covered query positions are calculated by subtracting the query alignment end coordinate from the query alignment start coordinate. Alignments are discarded if their query coverage is below a minimum threshold (default is 95%).

## Identifying the best-matching reference sequences for each genomic fragment
The following steps are performed on the filtered blastn alignments to identify each fragment's best-matching reference sequences:
1.	A combined bitscore is calculated for each fragment-reference sequence pairing. The combined bitscore is determined by adding together the individual bitscores for both fragment ends when aligned to the reference sequence.

2.	The best reference sequences for each fragment are chosen based on the highest combined bitscores calculated in the previous step. Best reference sequences are only retained if both ends of the fragment aligned. A fragment is discarded if all of its best reference sequences do not share the same segment and subtype annotations.

## FindFlu output
FindFlu outputs 6 files to the directory specified at runtime by the -o parameter. File names are prepended with the analysis name provided at runtime by the -n parameter.

•	<b>flu_frags_seqs.fa</b>: A FASTA file containing the IAV fragments identified by FindFlu. As in the input file, each fragment is described by a pair on consensus sequences (one for each end of the fragment).

•	<b>frag_report.csv</b>: A line list indicating segment and subtype for each fragment. Additional metrics are provided for each fragment (described below).

•	<b>best_ref_seqs_report.csv</b>: A line list indicating the best-matching reference sequences for each fragment.

•	<b>host_report.csv</b>: A table summarizing the host species annotations of the best-matching reference sequences for each fragment.

•	<b>country_report.csv</b>: A table summarizing the collection country annotations of the best-matching reference sequences for each fragment.

•	<b>H5_clade_report.csv</b>: A table summarizing the H5 clade annotations of the best-matching reference sequences for each H5 fragment.

## Fragment report
A line list that provides various metrics about each IAV fragment detected. It contains the following columns:

•	<b>experiment</b>: the name of the experiment to which the fragment belongs (extracted from the HopDropper header).

•	<b>library</b>: the name of the library in which the fragment was detected (extracted from the HopDropper header).

•	<b>fragment</b>: the fragment identifier assigned by HopDropper (extracted from the HopDropper header).

•	<b>UMI_pair</b>: the UMI pair describing the fragment (extracted from the HopDropper header). This can be used to identify the same fragment across different experiments (when it might have been assigned a different fragment identifier within its library).

•	<b>frag_copies</b>: the number of times the fragment was sequenced (extracted from the HopDropper header).

•	<b>segment</b>: the segment of the IAV genome from which the fragment originated (one of PB2, PB1, PA, HA, NP, NA, M, or NS).

•	<b>subtype</b>: the HA or NA subtype of the HA or NA genome segment from which the fragment originated (for internal segment, subtype is left blank).

•	<b>combined_query_id</b>: the combined alignment identity percentage of the fragment when it is aligned to its best-matching reference sequence, i.e. the percentage nucleotide similarity between the fragment and its best-matching reference sequence. The combined alignment identity percentage is calculated as follows: the number of identical positions in both fragment end alignments are summed, and the number of valid positions in both fragment end alignments are summed. These sums are divided to obtain the combined alignment identity percentage. Identical alignment positions are those where both the query and subject sequences contain the same non-ambiguous nucleotide. Valid alignment positions are those where both the query and subject sequences contain non-ambiguous nucleotide bases. For fragments with multiple best-matching reference sequences, a combined alignment identity percentage is calculated for each reference sequence, then the median value is reported.

•	<b>combined_query_cov</b>: the combined query coverage percentage of the fragment when it is aligned to its best-matching reference sequence, i.e. the extent of the fragment end sequences that aligned to the best-matching reference sequence. The combined query coverage percentage is calculated as follows: The number of covered query positions from both fragment end alignments are summed, and the length of both fragment end sequences are summed. These sums are divided to obtain the combined query coverage percentage. Covered query positions are calculated by subtracting the query alignment end coordinate from the query alignment start coordinate. For fragments with multiple best-matching reference sequences, a combined query coverage percentage is calculated for each reference sequence, then the median value is reported.

•	<b>frag_length</b>: the estimated length of the fragment based on its alignment to its best-matching reference sequence. Fragment start and end coordinates are determined from the minimum and maximum subject coordinates from both fragment end alignments. Fragment length is calculated by subtracting the fragment start coordinate from the fragment end coordinate. For fragments with multiple best-matching reference sequences, a fragment length is estimated for each reference sequence, then the median value is reported.

•	<b>frag_subject_cov</b>: the percentage of the fragment’s best-matching reference sequence covered by the fragment. This is calculated by dividing the estimated length of the fragment by the length of the best-matching reference sequence. For fragments with multiple best-matching reference sequences, coverage is estimated for each reference sequence, then the median value is reported.

## Best reference sequences report
A line list that reports the best-matching reference sequences for each fragment and associated alignment metrics. Each line describes one fragment-reference sequence pairing. It contains the following columns:

•	<b>experiment</b>: the name of the experiment to which the fragment belongs (extracted from the HopDropper header).

•	<b>library</b>: the name of the library in which the fragment was detected (extracted from the HopDropper header).

•	<b>fragment</b>: the fragment identifier assigned by HopDropper (extracted from the HopDropper header).

•	<b>UMI_pair</b>: the UMI pair describing the fragment (extracted from the HopDropper header). This can be used to identify the same fragment across different experiments (when it might have been assigned a different fragment identifier within its library).

•	<b>frag_copies</b>: the number of times the fragment was sequenced (extracted from the HopDropper header).

•	<b>ref_seq_ID</b>: the unique identifier of the reference sequence.

•	<b>ref_seq</b>: the strain name of the reference sequence and its HxNx subtype in parentheses.

•	<b>segment</b>: the IAV genome segment described by the reference sequence (one of PB2, PB1, PA, HA, NP, NA, M, or NS).

•	<b>subtype</b>: the HA or NA subtype of the HA or NA genome segment described by the reference sequence (for internal segment, subtype is left blank).

•	<b>combined_query_id</b>: the combined alignment identity percentage of the fragment when it is aligned to the reference sequence, i.e. the percentage nucleotide similarity between the fragment and its best-matching reference sequence. The combined alignment identity percentage is calculated as follows: the number of identical positions in both fragment end alignments are summed, and the number of valid positions in both fragment end alignments are summed. These sums are divided to obtain the combined alignment identity percentage. Identical alignment positions are those where both the query and subject sequences contain the same non-ambiguous nucleotide. Valid alignment positions are those where both the query and subject sequences contain non-ambiguous nucleotide bases.

•	<b>combined_query_cov</b>: the combined query coverage percentage of the fragment when it is aligned to the reference sequence, i.e. the extent of the fragment end sequences that aligned to the best-matching reference sequence. The combined query coverage percentage is calculated as follows: The number of covered query positions from both fragment end alignments are summed, and the length of both fragment end sequences are summed. These sums are divided to obtain the combined query coverage percentage. Covered query positions are calculated by subtracting the query alignment end coordinate from the query alignment start coordinate.

•	<b>frag_subject_cov</b>: the percentage of the reference sequence covered by the fragment. This is calculated by dividing the estimated length of the fragment by the length of the best-matching reference sequence.

## Host, country, and H5 clade reports
A line list that reports the host, country, and H5 clade annotations of fragments’ best-matching reference sequences. Each line describes the number of times a particular host/country/H5 clade annotation was observed among a fragment’s best-matching reference sequences. Only fragments originating from HA segments belonging to the H5 subtype are reported on the H5 clade report. It contains the following columns:

•	<b>experiment</b>: the name of the experiment to which the fragment belongs (extracted from the HopDropper header).

•	<b>library</b>: the name of the library in which the fragment was detected (extracted from the HopDropper header).

•	<b>fragment</b>: the fragment identifier assigned by HopDropper (extracted from the HopDropper header).

•	<b>UMI_pair</b>: the UMI pair describing the fragment (extracted from the HopDropper header). This can be used to identify the same fragment across different experiments (when it might have been assigned a different fragment identifier within its library).

•	<b>frag_copies</b>: the number of times the fragment was sequenced (extracted from the HopDropper header).

•	<b>segment</b>: the segment of the IAV genome from which the fragment originated (one of PB2, PB1, PA, HA, NP, NA, M, or NS).

•	<b>subtype</b>: the HA or NA subtype of the HA or NA genome segment from which the fragment originated (for internal segment, subtype is left blank).

•	<b>host/country/H5_clade</b>: one of the host species/collection country/H5 clade annotations observed among the best-matching reference sequences for this fragment. 

•	<b>host_count/country_count/H5_clade_count</b>: the number of times this line’s host species/collection country/H5 clade annotation appeared among this line’s fragment’s best-matching reference sequences.

•	<b>total_frag_matches</b>: the number of best-matching reference sequences associated with this line’s fragment.

•	<b>host_perc_of_total_matches/country_perc_of_total_matches/H5_clade_perc_of_total_matches</b>: the percentage of this line’s fragment’s best-matching reference sequences that had this line’s host/country/H5 clade annotation. This is calculated by dividing this line’s host_count/country_count/H5_clade_count by this line’s total_frag_matches.

## FindFlu usage
```python find_flu.py -i <input FASTA file> -o <output dir> -n <analysis name> -d <database FASTA file> [-I <min alignment identity> -c <min alignment cov> -t <blastn threads>]```

<b>Required arguments:</b>
```
  -i : input FASTA file
  -o : output directory
  -n : analysis name
  -d : database FASTA file
```
<b>Optional arguments:</b>
```
  -I : minimum fragment end alignment identity percentage (default = 95)
  -c : minimum fragment end alignment coverage percentage (default = 95)
  -t : number of threads to use for blastn (default = 1)
```
