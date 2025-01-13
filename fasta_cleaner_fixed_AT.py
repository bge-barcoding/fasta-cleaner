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

def parse_arguments():
    """Set up command line argument parsing with clear help text and validation."""
    parser = argparse.ArgumentParser(
        description='Process FASTA files to filter sequences based on AT content.',
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
        '-t', '--threshold',
        type=float,
        default=0.5,
        help='AT content threshold (0.0-1.0). Sequences with AT content >= threshold are kept'
    )
    
    parser.add_argument(
        '-c', '--consensus_threshold',
        type=float,
        default=0.7,
        help='Consensus threshold (0.0-1.0). Minimum frequency required for a base to be included in consensus'
    )
    
    args = parser.parse_args()
    
    if not os.path.isdir(args.input_dir):
        parser.error(f"Input directory does not exist: {args.input_dir}")
    
    if not 0 <= args.threshold <= 1:
        parser.error("Threshold must be between 0.0 and 1.0")
    
    return args

def analyze_sequence(record):
    """Calculate sequence metrics (length excluding gaps and AT content)."""
    sequence = str(record.seq).upper()
    sequence_no_gaps = sequence.replace('-', '')
    length_no_gaps = len(sequence_no_gaps)
    
    if length_no_gaps == 0:
        return length_no_gaps, 0.0
    
    a_count = sequence_no_gaps.count('A')
    t_count = sequence_no_gaps.count('T')
    at_content = (a_count + t_count) / length_no_gaps
    
    return length_no_gaps, at_content

def write_sequence_report(sequences, at_threshold, output_path):
    """Write sequence metrics to CSV file."""
    with open(output_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Sequence Name', 'Length (excluding gaps)', 'AT ratio', 'Status'])
        
        for record in sequences:
            length_no_gaps, at_content = analyze_sequence(record)
            status = 'kept' if at_content >= at_threshold else 'removed'
            writer.writerow([
                record.id,
                length_no_gaps,
                f"{at_content:.3f}",
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

def process_fasta_file(input_file, at_threshold=0.5, consensus_threshold=0.7):
    """Process a single FASTA file and return alignments and metrics."""
    alignment = AlignIO.read(input_file, "fasta")
    
    kept_records = []
    removed_records = []
    
    for record in alignment:
        _, at_content = analyze_sequence(record)
        if at_content >= at_threshold:
            kept_records.append(record)
        else:
            removed_records.append(record)
    
    kept_alignment = None
    removed_alignment = None
    consensus_record = None
    
    if kept_records:
        kept_alignment = MultipleSeqAlignment(kept_records)
        summary_align = AlignInfo.SummaryInfo(kept_alignment)
        consensus = summary_align.dumb_consensus(threshold=consensus_threshold)
        
        base_name = get_base_filename(os.path.basename(input_file))
        consensus_id = f"{base_name}_t{at_threshold}_c{consensus_threshold}"
        
        consensus_record = SeqRecord(
            Seq(str(consensus)),
            id=consensus_id,
            description=""
        )
    
    if removed_records:
        removed_alignment = MultipleSeqAlignment(removed_records)
    
    return kept_alignment, removed_alignment, consensus_record, alignment

def process_directory(input_dir, output_dir, at_threshold=0.5, consensus_threshold=0.7):
    """Process all FASTA files in a directory."""
    os.makedirs(output_dir, exist_ok=True)
    
    processed_count = 0
    
    for filename in os.listdir(input_dir):
        if filename.lower().endswith(('.fasta', '.fas', '.fa')):
            input_path = os.path.join(input_dir, filename)
            base_name = get_base_filename(filename)
            
            kept_alignment, removed_alignment, consensus, all_sequences = process_fasta_file(
                input_path, at_threshold, consensus_threshold
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
            
            write_sequence_report(all_sequences, at_threshold, report_path)
            
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
    process_directory(args.input_dir, args.output_dir, args.threshold, args.consensus_threshold)
