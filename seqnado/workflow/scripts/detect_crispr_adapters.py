"""
Detect CRISPR adapter sequences from FASTQ files.

This script analyzes the first N reads from FASTQ files to identify common
adapter sequences at the 5' end of reads. It uses k-mer frequency analysis
to find the most common sequences that appear at the beginning of reads.

For CRISPR libraries (like TKOv3):
- R1 typically has a short adapter (e.g., ACACCG from U6 promoter) followed by guide RNA
- R2 typically has the tracrRNA scaffold at the start
"""

import gzip
from pathlib import Path
from collections import Counter, defaultdict
from typing import Dict, Optional, Tuple, List
import json
from loguru import logger


def open_fastq(file_path: str):
    """Open FASTQ file, handling both gzipped and uncompressed files."""
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    return open(file_path, 'r')


def parse_fastq(file_path: str, n_reads: int = 10000) -> List[str]:
    """Parse first n_reads sequences from a FASTQ file."""
    sequences = []
    try:
        with open_fastq(file_path) as f:
            line_count = 0
            for line in f:
                if line_count % 4 == 1:  # Sequence line
                    sequences.append(line.strip())
                    if len(sequences) >= n_reads:
                        break
                line_count += 1
        logger.debug(f"Parsed {len(sequences)} sequences from {file_path}")
    except FileNotFoundError:
        logger.error(f"FASTQ file not found: {file_path}")
        raise
    except Exception as e:
        logger.error(f"Error parsing FASTQ file {file_path}: {e}")
        raise
    return sequences


def find_common_prefix(sequences: List[str], min_length: int = 6, max_length: int = 100) -> Dict[str, int]:
    """
    Find the most common k-mers at the start of sequences.

    Returns a dictionary of k-mer -> count for various k-mer lengths.
    """
    prefix_counts = defaultdict(Counter)

    for seq in sequences:
        # Try different k-mer lengths
        for k in range(min_length, min(max_length + 1, len(seq) + 1)):
            prefix = seq[:k]
            prefix_counts[k][prefix] += 1

    return prefix_counts


def detect_adapter(sequences: List[str],
                  min_length: int = 6,
                  max_length: int = 100,
                  min_frequency: float = 0.7) -> Optional[str]:
    """
    Detect the most likely adapter sequence.

    Strategy:
    1. Find common prefixes of various lengths
    2. Identify the longest prefix that appears in at least min_frequency of reads
    3. Return the consensus adapter sequence

    Args:
        sequences: List of sequence strings
        min_length: Minimum adapter length to consider
        max_length: Maximum adapter length to consider
        min_frequency: Minimum fraction of reads that must contain the adapter

    Returns:
        The detected adapter sequence, or None if no clear adapter found
    """
    if not sequences:
        logger.warning("No sequences provided for adapter detection")
        return None

    logger.info(f"Analyzing {len(sequences)} sequences for adapter detection")
    logger.info(f"Parameters: min_length={min_length}, max_length={max_length}, min_frequency={min_frequency}")

    prefix_counts = find_common_prefix(sequences, min_length, max_length)
    n_sequences = len(sequences)

    # Find the longest prefix that meets the frequency threshold
    best_adapter = None
    best_length = 0

    for k in sorted(prefix_counts.keys(), reverse=True):
        # Get the most common prefix of this length
        most_common = prefix_counts[k].most_common(1)
        if not most_common:
            continue

        prefix, count = most_common[0]
        frequency = count / n_sequences

        logger.debug(f"Length {k}: most common prefix '{prefix}' appears in {count}/{n_sequences} reads ({frequency:.2%})")

        if frequency >= min_frequency and k > best_length:
            best_adapter = prefix
            best_length = k
            logger.info(f"Found adapter candidate of length {k}: {prefix} (frequency: {frequency:.2%})")

    if best_adapter:
        logger.info(f"Best adapter detected: {best_adapter} (length: {best_length})")
    else:
        logger.warning(f"No adapter found meeting frequency threshold of {min_frequency:.2%}")
        # Log top candidates for troubleshooting
        logger.info("Top prefix candidates for troubleshooting:")
        for k in sorted(prefix_counts.keys(), reverse=True)[:5]:
            most_common = prefix_counts[k].most_common(3)
            for prefix, count in most_common:
                frequency = count / n_sequences
                logger.info(f"  Length {k}: '{prefix}' - {frequency:.2%} ({count}/{n_sequences} reads)")

    return best_adapter


def refine_adapter_boundary(sequences: List[str], adapter: str, expected_guide_length: int = 20) -> str:
    """
    Refine the adapter boundary by looking at sequence diversity after the adapter.

    Guide RNA sequences should be highly diverse, while adapter sequences are constant.
    This function finds where diversity increases, indicating the adapter/guide boundary.

    Args:
        sequences: List of sequence strings
        adapter: Initial adapter sequence
        expected_guide_length: Expected length of guide RNA

    Returns:
        Refined adapter sequence
    """
    if not adapter:
        return adapter

    # Calculate nucleotide diversity at each position after the initial adapter
    max_check_length = len(adapter) + 10  # Check a bit beyond the initial adapter

    position_entropy = []
    for pos in range(len(adapter) - 10, min(max_check_length, len(sequences[0]))):
        base_counts = Counter()
        for seq in sequences:
            if pos < len(seq):
                base_counts[seq[pos]] += 1

        # Calculate Shannon entropy as a measure of diversity
        total = sum(base_counts.values())
        entropy = 0
        for count in base_counts.values():
            if count > 0:
                p = count / total
                entropy -= p * (p ** 0.5)  # Simplified entropy measure

        position_entropy.append((pos, entropy))

    # Find where entropy increases sharply (transition from adapter to guide)
    if len(position_entropy) > 1:
        max_increase = 0
        best_pos = len(adapter)

        for i in range(1, len(position_entropy)):
            increase = position_entropy[i][1] - position_entropy[i-1][1]
            if increase > max_increase:
                max_increase = increase
                best_pos = position_entropy[i][0]

        if best_pos != len(adapter):
            logger.info(f"Refined adapter boundary from {len(adapter)} to {best_pos}")
            return adapter[:best_pos]

    return adapter


def detect_adapters_from_fastq(fq1: str,
                               fq2: Optional[str] = None,
                               n_reads: int = 10000,
                               min_length: int = 6,
                               max_length: int = 100,
                               min_frequency: float = 0.7) -> Dict[str, Optional[str]]:
    """
    Detect adapter sequences from FASTQ files.

    Args:
        fq1: Path to R1 FASTQ file
        fq2: Path to R2 FASTQ file (optional, for paired-end)
        n_reads: Number of reads to analyze
        min_length: Minimum adapter length
        max_length: Maximum adapter length
        min_frequency: Minimum frequency of adapter in reads

    Returns:
        Dictionary with 'adapter_r1' and 'adapter_r2' keys
    """
    logger.info(f"Analyzing {n_reads} reads from {fq1}")
    sequences_r1 = parse_fastq(fq1, n_reads)
    adapter_r1 = detect_adapter(sequences_r1, min_length, max_length, min_frequency)

    result = {
        'adapter_r1': adapter_r1,
        'adapter_r2': None
    }

    if fq2:
        logger.info(f"Analyzing {n_reads} reads from {fq2}")
        sequences_r2 = parse_fastq(fq2, n_reads)
        adapter_r2 = detect_adapter(sequences_r2, min_length, max_length, min_frequency)
        result['adapter_r2'] = adapter_r2

    return result


def main():
    """Main function for Snakemake integration."""
    try:
        # Setup logging
        logger.add(snakemake.log[0], level="INFO")
        logger.info("=" * 80)
        logger.info("CRISPR Adapter Detection Script")
        logger.info("=" * 80)

        # Get input files
        fq1 = snakemake.input.get('fq1', snakemake.input[0])
        fq2 = snakemake.input.get('fq2', None) if len(snakemake.input) > 1 else None

        # Get parameters
        n_reads = snakemake.params.get('n_reads', 10000)
        min_length = snakemake.params.get('min_length', 6)
        max_length = snakemake.params.get('max_length', 100)
        min_frequency = snakemake.params.get('min_frequency', 0.7)

        logger.info("Configuration:")
        logger.info(f"  Input R1: {fq1}")
        if fq2:
            logger.info(f"  Input R2: {fq2}")
        logger.info(f"  Number of reads to analyze: {n_reads}")
        logger.info(f"  Min adapter length: {min_length}")
        logger.info(f"  Max adapter length: {max_length}")
        logger.info(f"  Min frequency threshold: {min_frequency}")
        logger.info(f"  Output file: {snakemake.output[0]}")
        logger.info("")

        # Check input files exist
        from pathlib import Path
        if not Path(fq1).exists():
            logger.error(f"Input file R1 does not exist: {fq1}")
            raise FileNotFoundError(f"Input file R1 not found: {fq1}")
        if fq2 and not Path(fq2).exists():
            logger.error(f"Input file R2 does not exist: {fq2}")
            raise FileNotFoundError(f"Input file R2 not found: {fq2}")

        logger.info("Starting adapter detection...")
        logger.info("")

        # Detect adapters
        adapters = detect_adapters_from_fastq(
            fq1=fq1,
            fq2=fq2,
            n_reads=n_reads,
            min_length=min_length,
            max_length=max_length,
            min_frequency=min_frequency
        )

        logger.info("")
        logger.info("=" * 80)
        logger.info("RESULTS:")
        logger.info("=" * 80)
        logger.info(f"R1 adapter: {adapters['adapter_r1'] if adapters['adapter_r1'] else 'None detected'}")
        if fq2:
            logger.info(f"R2 adapter: {adapters['adapter_r2'] if adapters['adapter_r2'] else 'None detected'}")
        logger.info("=" * 80)
        logger.info("")

        # Write results
        output_path = snakemake.output[0]
        with open(output_path, 'w') as f:
            json.dump(adapters, f, indent=2)

        logger.info(f"Results successfully written to: {output_path}")
        logger.info("Adapter detection complete!")
        logger.info("=" * 80)

    except Exception as e:
        logger.error("=" * 80)
        logger.error("FATAL ERROR")
        logger.error("=" * 80)
        logger.error(f"An error occurred during adapter detection: {e}")
        logger.exception("Full traceback:")
        logger.error("=" * 80)
        raise


if __name__ == "__main__":
    main()
