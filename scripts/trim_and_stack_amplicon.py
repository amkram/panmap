import sys
import pysam
import argparse
from collections import defaultdict

class Primer:
  def __init__(self, chrom, start, end, name, strand):
    self.chrom = chrom
    self.start = start
    self.end = end
    self.name = name
    self.strand = strand

def load_primers(bed_file, offset):
  primers = []
  with open(bed_file, 'r') as f:
    for line in f:
      if line.startswith('#') or not line.strip():
        continue
      fields = line.strip().split('\t')
      if len(fields) < 6:
        continue
      
      chrom = fields[0]
      start = int(fields[1]) - offset
      end = int(fields[2]) - 1 + offset
      name = fields[3]
      strand = fields[5]
      
      primers.append(Primer(chrom, start, end, name, strand))
  
  primers.sort(key=lambda p: (p.chrom, p.start))
  return primers

def binary_search_primers(primers, chrom, target_pos):
  chrom_primers = [p for p in primers if p.chrom == chrom]
  if not chrom_primers:
    return []
  
  low, high = 0, len(chrom_primers) - 1
  matches = []
  
  while low <= high:
    mid = low + (high - low) // 2
    primer = chrom_primers[mid]
    
    if target_pos >= primer.start and target_pos <= primer.end:
      matches.append(primer)
      break
    elif target_pos > primer.end:
      low = mid + 1
    else:
      high = mid - 1
  
  if not matches:
    return []
  
  mid_idx = chrom_primers.index(matches[0])
  
  i = 1
  while True:
    found_any = False
    
    right_idx = mid_idx + i
    if right_idx < len(chrom_primers):
      p = chrom_primers[right_idx]
      if target_pos >= p.start and target_pos <= p.end:
        matches.append(p)
        found_any = True
      elif target_pos < p.start:
        pass
    
    left_idx = mid_idx - i
    if left_idx >= 0:
      p = chrom_primers[left_idx]
      if target_pos >= p.start and target_pos <= p.end:
        matches.append(p)
        found_any = True
      elif target_pos > p.end:
        pass
    
    if not found_any and (right_idx >= len(chrom_primers) and left_idx < 0):
      break
    
    i += 1
    if i > len(chrom_primers):
      break
  
  return matches

def get_read_start_position(read):
  if read.is_unmapped:
    return None, None, None
  
  chrom = read.reference_name
  
  if read.is_reverse:
    start_pos = read.reference_end - 1
    strand = '-'
  else:
    start_pos = read.reference_start
    strand = '+'
  
  return chrom, start_pos, strand

def select_candidate_primer(overlapping_primers, read_strand):
  if not overlapping_primers:
    return None
  
  strand_filtered = [p for p in overlapping_primers 
                     if p.strand == read_strand or p.strand == '0']
  
  if not strand_filtered:
    return None
  
  if read_strand == '-':
    return min(strand_filtered, key=lambda p: p.start)
  else:
    return max(strand_filtered, key=lambda p: p.end)

def get_trimmed_sequence_and_quality(read):
  if read.is_unmapped or not read.cigartuples:
    return read.query_sequence, read.query_qualities
  
  query_start = read.query_alignment_start
  query_end = read.query_alignment_end
  
  trimmed_seq = read.query_sequence[query_start:query_end]
  
  if read.query_qualities is not None:
    trimmed_qual = read.query_qualities[query_start:query_end]
  else:
    trimmed_qual = None
  
  return trimmed_seq, trimmed_qual

def quality_to_phred_string(quality_array):
  if quality_array is None:
    return '*'
  return ''.join(chr(q + 33) for q in quality_array)

def write_fastq_record(output_file, read_name, sequence, quality):
  output_file.write(f"@{read_name}\n")
  output_file.write(f"{sequence}\n")
  output_file.write("+\n")
  output_file.write(f"{quality}\n")

def trim_command(args):
  bam = pysam.AlignmentFile(args.bam_file, 'rb')
  output_file = open(args.output, 'w')
  
  read_count = 0
  unmapped_count = 0
  written_count = 0
  
  for read in bam:
    read_count += 1
    
    query_name = ''
    if read.is_read1:
      query_name = read.query_name + "/1"
    elif read.is_read2:
      query_name = read.query_name + "/2"
    else:
      query_name = read.query_name
    
    
    if read.is_unmapped:
      unmapped_count += 1
      if read.query_sequence:
        qual_string = quality_to_phred_string(read.query_qualities)
        write_fastq_record(output_file, query_name, read.query_sequence, qual_string)
        written_count += 1
      continue
    
    trimmed_seq, trimmed_qual = get_trimmed_sequence_and_quality(read)
    
    if trimmed_seq:
      qual_string = quality_to_phred_string(trimmed_qual)
      write_fastq_record(output_file, query_name, trimmed_seq, qual_string)
      written_count += 1
    
    if read_count % 100000 == 0:
      print(f"Processed {read_count} reads...", file=sys.stderr)
  
  output_file.close()
  bam.close()
  
  print(f"\n-------", file=sys.stderr)
  print(f"Results:", file=sys.stderr)
  print(f"Total reads processed: {read_count}", file=sys.stderr)
  print(f"Unmapped reads: {unmapped_count}", file=sys.stderr)
  print(f"Total reads written to FASTQ: {written_count}", file=sys.stderr)
  print(f"Output written to: {args.output}", file=sys.stderr)

def stack_command(args):
  primers = load_primers(args.primer_bed, args.offset)
  print(f"Loaded {len(primers)} primers from BED file", file=sys.stderr)
  
  bam = pysam.AlignmentFile(args.bam_file, 'rb')
  tsv_file = open(args.output, 'w')
  
  read_count = 0
  mapped_count = 0
  unmapped_count = 0
  no_primer_count = 0
  
  for read in bam:
    read_count += 1
    
    query_name = ''
    if read.is_read1:
      query_name = read.query_name + "/1"
    elif read.is_read2:
      query_name = read.query_name + "/2"
    else:
      query_name = read.query_name
    
    if read.is_unmapped:
      unmapped_count += 1
      continue
    
    result = get_read_start_position(read)
    if result is None or len(result) != 3:
      continue
    
    chrom, start_pos, read_strand = result
    
    overlapping_primers = binary_search_primers(primers, chrom, start_pos)
    candidate_primer = select_candidate_primer(overlapping_primers, read_strand)
    
    if candidate_primer:
      mapped_count += 1
      tsv_file.write(f"{query_name}\t{candidate_primer.name}\n")
    else:
      no_primer_count += 1
    
    if read_count % 100000 == 0:
      print(f"Processed {read_count} reads...", file=sys.stderr)
  
  tsv_file.close()
  bam.close()
  
  print(f"\n-------", file=sys.stderr)
  print(f"Results:", file=sys.stderr)
  print(f"Total reads processed: {read_count}", file=sys.stderr)
  print(f"Unmapped reads: {unmapped_count}", file=sys.stderr)
  print(f"Reads mapped to primers: {mapped_count}", file=sys.stderr)
  print(f"Reads without matching primers: {no_primer_count}", file=sys.stderr)
  print(f"Output written to: {args.output}", file=sys.stderr)

def main():
  parser = argparse.ArgumentParser(description='Map reads to primers and output trimmed FASTQ or stack depth.')
  subparsers = parser.add_subparsers(dest='command', help='Subcommands')
  
  trim_parser = subparsers.add_parser('trim', help='Trim reads and output FASTQ')
  trim_parser.add_argument('bam_file', type=str, help='Path to the BAM file')
  trim_parser.add_argument('-o', '--output', type=str, required=True, help='Output FASTQ file')
  
  stack_parser = subparsers.add_parser('stack', help='Count stack depth and output TSV')
  stack_parser.add_argument('bam_file', type=str, help='Path to the BAM file')
  stack_parser.add_argument('primer_bed', type=str, help='Path to the primer BED file')
  stack_parser.add_argument('--offset', type=int, default=3, help='Offset value (default: 3)')
  stack_parser.add_argument('-o', '--output', type=str, required=True, help='Output TSV file')
  
  args = parser.parse_args()
  
  if args.command == 'trim':
    trim_command(args)
  elif args.command == 'stack':
    stack_command(args)
  else:
    parser.print_help()
    sys.exit(1)

if __name__ == '__main__':
  main()