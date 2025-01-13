import os
import csv
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio.Align import AlignInfo
import tempfile

HUMAN_COX1 = """ATGTTCGCCGACCGTTGACTATTCTCTACAAACCACAAAGACATTGGAACACTATACCTATTATTCGGCG
CATGAGCTGGAGTCCTAGGCACAGCTCTAAGCCTCCTTATTCGAGCCGAGCTGGGCCAGCCAGGCAACCT
TCTAGGTAACGACCACATCTACAACGTTATCGTCACAGCCCATGCATTTGTAATAATCTTCTTCATAGTA
ATACCCATCATAATCGGAGGCTTTGGCAACTGACTAGTTCCCCTAATAATCGGTGCCCCCGATATGGCGT
TTCCCCGCATAAACAACATAAGCTTCTGACTCTTACCTCCCTCTCTCCTACTCCTGCTCGCATCTGCTAT
AGTGGAGGCCGGAGCAGGAACAGGTTGAACAGTCTACCCTCCCTTAGCAGGGAACTACTCCCACCCTGGA
GCCTCCGTAGACCTAACCATCTTCTCCTTACACCTAGCAGGTGTCTCCTCTATCTTAGGGGCCATCAATT
TCATCACAACAATTATCAATATAAAACCCCCTGCCATAACCCAATACCAAACGCCCCTCTTCGTCTGATC
CGTCCTAATCACAGCAGTCCTACTTCTCCTATCTCTCCCAGTCCTAGCTGCTGGCATCACTATACTACTA
ACAGACCGCAACCTCAACACCACCTTCTTCGACCCCGCCGGAGGAGGAGACCCCATTCTATACCAACACC
TATTCTGATTTTTCGGTCACCCTGAAGTTTATATTCTTATCCTACCAGGCTTCGGAATAATCTCCCATAT
TGTAACTTACTACTCCGGAAAAAAAGAACCATTTGGATACATAGGTATGGTCTGAGCTATGATATCAATT
GGCTTCCTAGGGTTTATCGTGTGAGCACACCATATATTTACAGTAGGAATAGACGTAGACACACGAGCAT
ATTTCACCTCCGCTACCATAATCATCGCTATCCCCACCGGCGTCAAAGTATTTAGCTGACTCGCCACACT
CCACGGAAGCAATATGAAATGATCTGCTGCAGTGCTCTGAGCCCTAGGATTCATCTTTCTTTTCACCGTA
GGTGGCCTGACTGGCATTGTATTAGCAAACTCATCACTAGACATCGTACTACACGACACGTACTACGTTG
TAGCCCACTTCCACTATGTCCTATCAATAGGAGCTGTATTTGCCATCATAGGAGGCTTCATTCACTGATT
TCCCCTATTCTCAGGCTACACCCTAGACCAAACCTACGCCAAAATCCATTTCACTATCATATTCATCGGC
GTAAATCTAACTTTCTTCCCACAACACTTTCTCGGCCTATCCGGAATGCCCCGACGTTACTCGGACTACC
CCGATGCATACACCACATGAAACATCCTATCATCTGTAGGCTCATTCATTTCTCTAACAGCAGTAATATT
AATAATTTTCATGATTTGAGAAGCCTTCGCTTCGAAGCGAAAAGTCCTAATAGTAGAAGAACCCTCCATA
AACCTGGAGTGACTATATGGATGCCCCCCACCCTACCACACATTCGAAGAACCCGTATACATAAAATCTA
GA""".replace("\n", "")

def parse_arguments():
    """Set up command line argument parsing with clear help text and validation."""
    parser = argparse.ArgumentParser(
        description='Process FASTA files to filter sequences based on AT content difference from consensus and human COX1 similarity.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '-i', '--input_dir',
        required=True,
        help='Directory containing input FASTA files'
    )
    
    parser.add_argument(
        '-o', '--output_dir',
        required=True,
        help='Directory where output files will be saved'
    )
    
    parser.add_argument(
        '-d', '--at_difference',
        type=float,
        default=0.1,
        help='Maximum allowed AT content difference from consensus (0.0-1.0). Sequences with larger differences are removed'
    )
    
    parser.add_argument(
        '-c', '--consensus_threshold',
        type=float,
        default=0.7,
        help='Consensus threshold (0.0-1.0). Minimum frequency required for a base to be included in consensus'
    )
    
    parser.add_argument(
        '-u', '--human_threshold',
        type=float,
        default=0.95,
        help='Human COX1 similarity threshold (0.0-1.0). Sequences with similarity >= threshold are removed'
    )
    
    args = parser.parse_args()
    
    if not os.path.isdir(args.input_dir):
        parser.error(f"Input directory does not exist: {args.input_dir}")
    
    if not 0 <= args.at_difference <= 1:
        parser.error("AT difference threshold must be between 0.0 and 1.0")
    
    if not 0 <= args.human_threshold <= 1:
        parser.error("Human threshold must be between 0.0 and 1.0")
    
    return args

def calculate_at_content(sequence):
    """Calculate AT content for a sequence, ignoring gaps."""
    sequence_no_gaps = sequence.replace('-', '').upper()
    
    if not sequence_no_gaps:
        return 0.0
    
    at_count = sequence_no_gaps.count('A') + sequence_no_gaps.count('T')
    return at_count / len(sequence_no_gaps)

def compare_at_content(consensus_seq, query_seq):
    """
    Compare AT content between consensus and query sequence in overlapping regions.
    Returns tuple of (query AT content, consensus AT content, difference).
    """
    # Convert sequences to uppercase for comparison
    consensus_seq = consensus_seq.upper()
    query_seq = query_seq.upper()
    
    # Find positions where both sequences have bases (not gaps)
    overlapping_positions = [
        (c, q) for c, q in zip(consensus_seq, query_seq)
        if c != '-' and q != '-'
    ]
    
    if not overlapping_positions:
        return 0.0, 0.0, 0.0
    
    # Split into consensus and query sequences
    consensus_bases = ''.join(c for c, _ in overlapping_positions)
    query_bases = ''.join(q for _, q in overlapping_positions)
    
    # Calculate AT content for overlapping regions
    consensus_at = calculate_at_content(consensus_bases)
    query_at = calculate_at_content(query_bases)
    
    # Calculate absolute difference
    at_difference = abs(consensus_at - query_at)
    
    return query_at, consensus_at, at_difference

def analyze_sequence(record, consensus_seq=None):
    """Calculate sequence metrics including length and AT content comparison with consensus."""
    sequence = str(record.seq).upper()
    sequence_no_gaps = sequence.replace('-', '')
    length_no_gaps = len(sequence_no_gaps)
    
    if length_no_gaps == 0:
        return length_no_gaps, 0.0, 0.0, 0.0
    
    # Calculate overall AT content
    at_content = calculate_at_content(sequence_no_gaps)
    
    # If consensus is provided, compare AT contents
    consensus_at = 0.0
    at_difference = 0.0
    if consensus_seq is not None:
        at_content, consensus_at, at_difference = compare_at_content(
            consensus_seq, sequence
        )
    
    return length_no_gaps, at_content, consensus_at, at_difference

def calculate_sequence_similarity(seq1, seq2):
    """
    Calculate sequence similarity between two sequences.
    Returns similarity ratio between 0 and 1.
    Ignores gaps in the comparison.
    """
    # Remove gaps from both sequences
    seq1_nogaps = seq1.replace('-', '').upper()
    seq2_nogaps = seq2.replace('-', '').upper()
    
    # Get the length of the shorter sequence
    min_length = min(len(seq1_nogaps), len(seq2_nogaps))
    
    if min_length == 0:
        return 0.0
    
    # Count matches
    matches = sum(1 for i in range(min_length) if seq1_nogaps[i] == seq2_nogaps[i])
    
    # Calculate similarity ratio
    return matches / min_length

def write_sequence_report(sequences, consensus_seq, at_difference_threshold, human_threshold, output_path):
    """Write sequence metrics to CSV file."""
    with open(output_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            'sequence_name',
            'sequence_length',
            'sequence_AT_ratio',
            'consensus_AT_ratio',
            'AT_ratio_diff',
            'human_similarity',
            'sequence_status'
        ])
        
        for record in sequences:
            length_no_gaps, at_content, consensus_at, at_difference = analyze_sequence(
                record, consensus_seq
            )
            human_similarity = calculate_sequence_similarity(str(record.seq), HUMAN_COX1)
            
            # Determine status based on both AT difference and human similarity
            if human_similarity >= human_threshold:
                status = 'removed (human)'
            elif at_difference > at_difference_threshold:
                status = 'removed (AT)'
            else:
                status = 'kept'
                
            writer.writerow([
                record.id,
                length_no_gaps,
                f"{at_content:.3f}",
                f"{consensus_at:.3f}",
                f"{at_difference:.3f}",
                f"{human_similarity:.3f}",
                status
            ])

def write_fasta_single_line(records, filename):
    """Write sequences to FASTA file without line wrapping."""
    with open(filename, 'w') as handle:
        for record in records:
            handle.write(f">{record.id}\n")
            handle.write(f"{str(record.seq)}\n")

def get_base_filename(filepath):
    """Extract base filename without extension."""
    basename = os.path.basename(filepath)
    for ext in ['.fasta', '.fas', '.fa']:
        if basename.lower().endswith(ext):
            return basename[:-len(ext)]
    return basename

def process_fasta_file(input_file, at_difference_threshold=0.1, consensus_threshold=0.7, human_threshold=0.95):
    """Process a single FASTA file and return alignments and metrics."""
    alignment = AlignIO.read(input_file, "fasta")
    
    # First generate consensus sequence
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.dumb_consensus(threshold=consensus_threshold)
    consensus_seq = str(consensus)
    
    kept_records = []
    removed_records = []
    
    for record in alignment:
        _, _, _, at_difference = analyze_sequence(record, consensus_seq)
        human_similarity = calculate_sequence_similarity(str(record.seq), HUMAN_COX1)
        
        # First check human similarity, then AT difference
        if human_similarity >= human_threshold:
            removed_records.append(record)
        elif at_difference > at_difference_threshold:
            removed_records.append(record)
        else:
            kept_records.append(record)
    
    kept_alignment = None
    removed_alignment = None
    consensus_record = None
    
    if kept_records:
        kept_alignment = MultipleSeqAlignment(kept_records)
        # Generate new consensus from kept sequences
        summary_align = AlignInfo.SummaryInfo(kept_alignment)
        consensus = summary_align.dumb_consensus(threshold=consensus_threshold)
        
        base_name = get_base_filename(os.path.basename(input_file))
        consensus_id = f"{base_name}_d{at_difference_threshold}_c{consensus_threshold}"
        
        consensus_record = SeqRecord(
            Seq(str(consensus)),
            id=consensus_id,
            description=""
        )
    
    if removed_records:
        removed_alignment = MultipleSeqAlignment(removed_records)
    
    return kept_alignment, removed_alignment, consensus_record, alignment, consensus_seq

def process_directory(input_dir, output_dir, at_difference_threshold=0.1, consensus_threshold=0.7, human_threshold=0.95):
    """Process all FASTA files in a directory."""
    os.makedirs(output_dir, exist_ok=True)
    
    processed_count = 0
    
    for filename in os.listdir(input_dir):
        if filename.lower().endswith(('.fasta', '.fas', '.fa')):
            input_path = os.path.join(input_dir, filename)
            base_name = get_base_filename(filename)
            
            kept_alignment, removed_alignment, consensus, all_sequences, consensus_seq = process_fasta_file(
                input_path, at_difference_threshold, consensus_threshold, human_threshold
            )
            
            cleaned_path = os.path.join(output_dir, f"{base_name}_cleaned.fasta")
            removed_path = os.path.join(output_dir, f"{base_name}_removed.fasta")
            consensus_path = os.path.join(output_dir, f"{base_name}_cleaned_consensus.fasta")
            report_path = os.path.join(output_dir, f"{base_name}_report.csv")
            
            original_count = len(all_sequences)
            kept_count = len(kept_alignment) if kept_alignment else 0
            removed_count = len(removed_alignment) if removed_alignment else 0
            
            if kept_alignment:
                write_fasta_single_line(kept_alignment, cleaned_path)
            if removed_alignment:
                write_fasta_single_line(removed_alignment, removed_path)
            if consensus:
                write_fasta_single_line([consensus], consensus_path)
            
            write_sequence_report(
                all_sequences,
                consensus_seq,
                at_difference_threshold,
                human_threshold,
                report_path
            )
            
            print(f"Processed {filename}:")
            print(f"  - Original sequences: {original_count}")
            print(f"  - Kept sequences: {kept_count}")
            print(f"  - Removed sequences: {removed_count}")
            print(f"  - Saved cleaned alignment to: {cleaned_path}")
            print(f"  - Saved removed sequences to: {removed_path}")
            print(f"  - Saved consensus to: {consensus_path}")
            print(f"  - Saved sequence report to: {report_path}")
            print()
            
            processed_count += 1
    
    if processed_count == 0:
        print(f"No FASTA files found in {input_dir}")
    else:
        print(f"Finished processing {processed_count} FASTA files")

if __name__ == "__main__":
    args = parse_arguments()
    process_directory(
        args.input_dir,
        args.output_dir,
        args.at_difference,
        args.consensus_threshold,
        args.human_threshold
    )
