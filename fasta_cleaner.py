import os
import csv
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio.Align import AlignInfo
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
import time

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
    
    parser.add_argument(
        '-T', '--threads',
        type=int,
        default=multiprocessing.cpu_count(),
        help='Number of parallel threads to use'
    )
    
    args = parser.parse_args()
    
    if not os.path.isdir(args.input_dir):
        parser.error(f"Input directory does not exist: {args.input_dir}")
    
    if not 0 <= args.threshold <= 1:
        parser.error("Threshold must be between 0.0 and 1.0")
        
    if args.threads < 1:
        parser.error("Number of threads must be at least 1")
    
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

def process_single_file(input_file, output_dir, at_threshold, consensus_threshold):
    """Process a single FASTA file with the given parameters."""
    base_name = get_base_filename(os.path.basename(input_file))
    
    # First check if file is empty
    if os.path.getsize(input_file) == 0:
        print(f"Warning: {input_file} is empty. Skipping...")
        return {"status": "empty", "filename": input_file}
        
    try:
        alignment = AlignIO.read(input_file, "fasta")
    except ValueError as e:
        print(f"Warning: Could not read {input_file}: {str(e)}")
        print("File may be corrupted or malformed. Skipping...")
        return {"status": "corrupted", "filename": input_file}
    except Exception as e:
        print(f"Error processing {input_file}: {str(e)}")
        return {"status": "error", "filename": input_file}
    
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
        
        # Convert consensus string to list for manipulation
        consensus_seq = list(str(consensus))
        
        # Get alignment length
        align_len = kept_alignment.get_alignment_length()
        
        # For each position in the alignment
        for i in range(align_len):
            # Get all bases at this position (excluding gaps)
            bases = [record.seq[i] for record in kept_records if record.seq[i] != '-']
            
            if consensus_seq[i] == 'X':
                if len(bases) == 0:  # No bases at this position (all gaps)
                    consensus_seq[i] = '-'
                else:  # Has bases but doesn't meet threshold
                    consensus_seq[i] = 'N'
        
        # Convert back to string
        consensus_seq = ''.join(consensus_seq)
        
        consensus_id = f"{base_name}_t{at_threshold}_c{consensus_threshold}"
        consensus_record = SeqRecord(
            Seq(consensus_seq),
            id=consensus_id,
            description=""
        )
    
    if removed_records:
        removed_alignment = MultipleSeqAlignment(removed_records)
    
    # Write output files
    cleaned_path = os.path.join(output_dir, f"{base_name}_cleaned.fasta")
    removed_path = os.path.join(output_dir, f"{base_name}_removed.fasta")
    consensus_path = os.path.join(output_dir, "cleaned_consensus", f"{base_name}_cleaned_consensus.fasta")
    report_path = os.path.join(output_dir, f"{base_name}_report.csv")
    
    if kept_alignment:
        write_fasta_single_line(kept_alignment, cleaned_path)
    if removed_alignment:
        write_fasta_single_line(removed_alignment, removed_path)
    if consensus_record:
        write_fasta_single_line([consensus_record], consensus_path)
    
    write_sequence_report(alignment, at_threshold, report_path)
    
    result = {
        "status": "success",
        "filename": input_file,
        "original_count": len(alignment),
        "kept_count": len(kept_alignment) if kept_alignment else 0,
        "removed_count": len(removed_alignment) if removed_alignment else 0,
        "paths": {
            "cleaned": cleaned_path,
            "removed": removed_path,
            "consensus": consensus_path,
            "report": report_path
        }
    }
    
    return result

def process_directory(input_dir, output_dir, at_threshold=0.5, consensus_threshold=0.7, num_threads=None):
    """Process all FASTA files in a directory using parallel processing."""
    os.makedirs(output_dir, exist_ok=True)
    # Create consensus subdirectory
    consensus_dir = os.path.join(output_dir, "cleaned_consensus")
    os.makedirs(consensus_dir, exist_ok=True)
    
    # Create a log file for skipped files
    log_file = os.path.join(output_dir, "processing_log.txt")
    with open(log_file, 'w') as f:
        f.write("FASTA Processing Log\n")
        f.write("==================\n\n")
    
    # Get list of FASTA files
    fasta_files = [
        os.path.join(input_dir, f) for f in os.listdir(input_dir)
        if f.lower().endswith(('.fasta', '.fas', '.fa'))
    ]
    
    if not fasta_files:
        print(f"No FASTA files found in {input_dir}")
        return
    
    start_time = time.time()
    
    # Process files in parallel
    with ProcessPoolExecutor(max_workers=num_threads) as executor:
        # Create a list of futures
        futures = [
            executor.submit(process_single_file, f, output_dir, at_threshold, consensus_threshold)
            for f in fasta_files
        ]
        
        # Process results as they complete
        results = []
        for future in futures:
            result = future.result()
            results.append(result)
            
            # Print progress for successful processing
            if result['status'] == 'success':
                print(f"Processed {os.path.basename(result['filename'])}:")
                print(f"  - Original sequences: {result['original_count']}")
                print(f"  - Kept sequences: {result['kept_count']}")
                print(f"  - Removed sequences: {result['removed_count']}")
                print(f"  - Saved cleaned alignment to: {result['paths']['cleaned']}")
                print(f"  - Saved removed sequences to: {result['paths']['removed']}")
                print(f"  - Saved consensus to: {result['paths']['consensus']}")
                print(f"  - Saved sequence report to: {result['paths']['report']}")
                print()
    
    # Compile statistics
    processed_count = sum(1 for r in results if r['status'] == 'success')
    empty_count = sum(1 for r in results if r['status'] == 'empty')
    corrupted_count = sum(1 for r in results if r['status'] == 'corrupted')
    error_count = sum(1 for r in results if r['status'] == 'error')
    
    end_time = time.time()
    total_time = end_time - start_time
    
    # Write final summary
    summary = f"""
Processing Summary:
==================
Total files processed successfully: {processed_count}
Empty files skipped: {empty_count}
Corrupted files skipped: {corrupted_count}
Files with processing errors: {error_count}
Total processing time: {total_time:.2f} seconds
Average time per file: {total_time/len(fasta_files):.2f} seconds

Detailed log written to: {log_file}
"""
    print(summary)
    
    # Write summary to log file
    with open(log_file, 'a') as f:
        f.write(summary)
        
        # Write detailed results for each file
        f.write("\nDetailed Results:\n")
        f.write("================\n")
        for result in results:
            f.write(f"\nFile: {os.path.basename(result['filename'])}\n")
            f.write(f"Status: {result['status']}\n")
            if result['status'] == 'success':
                f.write(f"Original sequences: {result['original_count']}\n")
                f.write(f"Kept sequences: {result['kept_count']}\n")
                f.write(f"Removed sequences: {result['removed_count']}\n")
    
    # Concatenate all consensus files from the cleaned_consensus subdirectory
    consensus_dir = os.path.join(output_dir, "cleaned_consensus")
    concat_consensus_path = os.path.join(output_dir, "concat_cleaned_consensus.fasta")
    
    with open(concat_consensus_path, 'w') as outfile:
        for consensus_file in sorted(os.listdir(consensus_dir)):
            if consensus_file.endswith('.fasta'):
                with open(os.path.join(consensus_dir, consensus_file)) as infile:
                    outfile.write(infile.read())
    
    print(f"\nConcatenated consensus sequences saved to: {concat_consensus_path}")

if __name__ == "__main__":
    args = parse_arguments()
    process_directory(
        args.input_dir,
        args.output_dir,
        args.threshold,
        args.consensus_threshold,
        args.threads
    )
