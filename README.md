# FindFlu
FindFlu is a bioinformatic tool designed to identify and characterize fragments of influenza A virus (IAV) genome from user-provided sequences. For each IAV genome fragment identified, FindFlu will report the segment of the IAV genome from which it originated, the subtype of that segment (for fragments of haemagglutinin and neuraminidase segments), and the length of the fragment. FindFlu identifies and characterizes fragments by aligning them to a database of IAV reference sequences and finding each fragment’s best-matching reference sequence.

## FindFlu input
As input, FindFlu expects a FASTA file containing pairs of consensus sequences describing both ends of distinct genomic fragments. The entries in this FASTA file must have headers that use the HopDropper (v1.0.0) format (https://github.com/KevinKuchinski/HopDropper).

## FindFlu database
FindFlu aligns input sequences to a database of annotated IAV reference sequences. Reference sequences must be in a FASTA file with the following header format:

```>ref_seq_identifier|strain_name(strain_subtype)|segment|subtype|```

Ref_seq_identifier is a unique identifier for the reference sequence, typically its GenBank accession number. Strain_name is the name of the strain (e.g. A/goose/Guangdong/1/1996), and strain_subtype is the HxNx subtype of the strain (regardless of which IAV genome segment the reference sequence represents). Segment is the IAV genome segment the reference sequence represents, and it must be encoded as one of the following: PB2, PB1, PA, HA, NP, NA, M, or NS. Subtype is the HA subtype (for HA segment sequences) or the NA subtype (for NA segment sequences). For internal segments, ‘none’ should be entered in the subtype field.

## Alignment and filtering
FindFlu aligns input sequences to the reference sequence database with blastn (v2.13.0). The following filtering is performed  on the blastn alignments:
1.  Alignments are discarded if their percent identity is less than the minimum identity threshold (-I, default=90%).
2.  Alignments are discarded if their query coverage is less than the minimum coverage threshold (-c, default=95%).
3.  Each fragment end sequence can only align once to a reference sequence, otherwise all alignments between that fragment end and that reference sequence are discarded.
4.  Both ends of a fragment must align to a reference sequence, otherwise all alignments between either end of that fragment and that reference sequence are discarded.
5.  For each fragment-reference sequence pairing, one end must align in the plus sense and the other end must align in the minus sense.

## Identifying the best-matching reference sequences for each genomic fragment
The following steps are performed on the filtered blastn alignments to identify each fragment's best-matching reference sequences:
1.	A combined bitscore is calculated for each fragment-reference sequence pairing. The combined bitscore is determined by summing the bitscores for both fragment ends alignments to the reference sequence.
2.	The best reference sequences for each fragment are chosen based on the top <i>n</i>-ranked combined bitscores calculated in the previous step (-r, default=2).
3.	All best reference sequences for a fragment must share the same segment and subtype annotation, otherwise that fragment is discarded.
4.	The length of the fragment is estimated based on the distance between alignment coordinates for both fragment ends. If a fragment has multiple best reference sequences, a fragment length is estimated for each reference sequence, then the median fragment length is chosen.

Once the length of each fragment has been estimated, the consensus sequences for both ends of a fragment are merged into a single consensus sequence. If the fragment end sequences do not overlap each other, FindFlu inserts the appropriate number of N characters between them based on the distance between the two fragment ends and the number of positions sequenced at each fragment end. The merged fragment sequences are then re-aligned against the FindFlu database so that the number mismatches, insertions, and deletions between fragments and each of their best-matching reference sequences can be calculated.

## FindFlu output
FindFlu outputs 3 files to the directory specified at runtime by the -o parameter. File names are prepended with the analysis name provided at runtime by the -n parameter.

•	<b>flu_frags_seqs.fa</b>: A FASTA file containing the IAV fragments identified by FindFlu. Each fragment is described by a single consensus sequence. Fragments inherit their FASTA headers from the input file, but with segment and subtype appended.

•	<b>top_ref_seqs_report.csv</b>: This is FindFlu’s most detailed report. Each line describes the alignment of one fragment to one of its best-matching reference sequences. Columns are provided that indicate fragment length and the segment and subtype of the fragment/reference sequence. Alignment metrics (identity, coverage, mismatches, gaps, etc) are also provided.

•	<b>frag_report.csv</b>: This report is a simplified version of the best reference sequences report; it is a line list wherein information about each fragment is summarized on one line.

## Best reference sequences report
Each line of this report describes the alignment of one fragment to one of its best-matching reference sequences. It contains the following columns:

•	<b>experiment</b>: the name of the experiment to which this line's fragment belongs (extracted from the HopDropper header).

•	<b>lib_name</b>: the name of the library in which this line's fragment was detected (extracted from the HopDropper header).

•	<b>frag_name</b>: the fragment identifier assigned to this line's fragment by HopDropper (extracted from the HopDropper header).

•	<b>UMI_pair</b>: the UMI pair describing this line's fragment (extracted from the HopDropper header). This can be used to identify the same fragment across different experiments (when it might have been assigned a different fragment identifier within its library).

•	<b>frag_copies</b>: the number of times this line's fragment was sequenced (extracted from the HopDropper header).

•	<b>frag_length</b>: The estimated length (in nucleotide positions) of this line's fragment based on this line's reference sequence.

•	<b>perc_sequenced</b>: The estimated percentage of positions in the fragment that were sequenced (based on the fragment length estimate, the length of both fragment end sequences, and the length of any overlap between the fragment end sequences).

•	<b>ref_seq_name</b>: the strain name of this line's reference sequence and its HxNx subtype in parentheses.

•	<b>ref_seq_length</b>: the length of this line's reference sequence.

•	<b>segment</b>: the IAV genome segment of this line's reference sequence (one of PB2, PB1, PA, HA, NP, NA, M, or NS).

•	<b>subtype</b>: the HA or NA subtype of this line's reference sequence ("none" for internal segments).

•	<b>aligned</b>: the number of positions in the fragment that aligned to this line's reference sequence.

•	<b>identical</b>: the number of identical positions in the alignment between this line's fragment and reference sequence.

•	<b>mismatches</b>: the number of mismatched positions in the alignment between this line's fragment and reference sequence (positions with Ns do NOT count as mismatches).

•	<b>deletions</b>: the number of deletions in the this line's fragment when aligned this line's reference sequence.

•	<b>insertions</b>: the number of insertions in the this line's fragment when aligned this line's reference sequence.

•	<b>ambiguous</b>: the number of ambiguous positions in the alignment between this line's fragment and reference sequence.

•	<b>ambiguous</b>: the number of ambiguous positions in the alignment between this line's fragment and reference sequence.

•	<b>perc_aligned</b>: the percentage of positions in this line's fragment that aligned to this line's reference sequence.

•	<b>perc_identity</b>: the percentage of positions in this line's fragment that aligned to this line's reference sequence and were identical.

•	<b>perc_ambiguous</b>: the percentage of positions in this line's fragment that aligned to this line's reference sequence and were ambiguous.

•	<b>ref_seq_cov_sequenced</b>: the percentage of positions in this line's reference sequence that were covered by this line's fragment and sequenced.

•	<b>ref_seq_cov_total</b>: the percentage of positions in this line's reference sequence that were covered by this line's fragment including ambiguous positions and unsequenced positions between the fragment end sequences.

## Fragment report
This report is a line list wherein information about each fragment is summarized on one line. It contains the following columns:

•	<b>experiment</b>: the name of the experiment to which this line's fragment belongs (extracted from the HopDropper header).

•	<b>lib_name</b>: the name of the library in which this line's fragment was detected (extracted from the HopDropper header).

•	<b>frag_name</b>: the fragment identifier assigned to this line's fragment by HopDropper (extracted from the HopDropper header).

•	<b>UMI_pair</b>: the UMI pair describing this line's fragment (extracted from the HopDropper header). This can be used to identify the same fragment across different experiments (when it might have been assigned a different fragment identifier within its library).

•	<b>frag_copies</b>: the number of times this line's fragment was sequenced (extracted from the HopDropper header).

•	<b>segment</b>: the IAV genome segment of this line's fragment (one of PB2, PB1, PA, HA, NP, NA, M, or NS).

•	<b>subtype</b>: the HA or NA subtype of this line's fragment ("none" for internal segments).

•	<b>frag_length</b>: The estimated length (in nucleotide positions) of this line's fragment (median of frag_length estimates based on all best reference sequences).

•	<b>perc_sequenced</b>: The estimated percentage of positions in the fragment that were sequenced based on the fragment length estimate, the length of both fragment end sequences, and the length of any overlap between the fragment end sequences (median of perc_sequenced based on all best reference sequences).

•	<b>perc_aligned</b>: the percentage of positions in this line's fragment that aligned to its best reference sequences (median of perc_aligned based on all best reference sequences).

•	<b>perc_identity</b>: the percentage of positions in this line's fragment that aligned to its best reference sequences and were identical (median of perc_identity based on all best reference sequences).

•	<b>perc_ambiguous</b>: the percentage of positions in this line's fragment that aligned to its best reference sequences and were ambiguous (median of perc_ambiguous based on all best reference sequences).

•	<b>ref_seq_cov_sequenced</b>: the percentage of positions in the genome segment covered by this line's fragment and sequenced (median of ref_seq_cov_sequenced based on all best reference sequences).

•	<b>ref_seq_cov_total</b>: the percentage of positions in the genome segment covered by this line's fragment including ambiguous positions and unsequenced positions between the fragment end sequences (median of ref_seq_cov_total based on all best reference sequences).

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
  -I : minimum fragment end alignment identity percentage (default = 90)
  -c : minimum fragment end alignment coverage percentage (default = 95)
  -t : number of threads to use for blastn (default = 1)
  -r : ranks of combined bitscore to consider as best matches (default=2)
```
