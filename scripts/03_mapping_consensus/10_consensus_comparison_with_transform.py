#!/usr/bin/env python3
"""
Consensus-to-Reference Comparison with Coordinate Transform

Compares consensus sequences with reference, detects variants, and
transforms coordinates from rotated reference back to original coordinates.

Usage:
    python 10_consensus_comparison_with_transform.py [--test SAMPLE_ID] [--jobs N] [--output-format {csv|tsv|json}]
"""

import os
import sys
import argparse
import csv
import json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
import logging

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError:
    print("Error: biopython required. Install: pip install biopython")
    sys.exit(1)

# Path configuration
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT_DIR = os.path.dirname(os.path.dirname(_SCRIPT_DIR))

CONSENSUS_DIR = os.environ.get(
    "HBV_CONSENSUS_DIR", os.path.join(_PROJECT_DIR, "Medaka_consensus_8/r2")
)
REFERENCE_FILE = os.environ.get("HBV_REF", "/path/to/reference.fasta")
OUTPUT_DIR = os.environ.get("HBV_OUTPUT_DIR", os.path.join(_PROJECT_DIR, "consensus_ref_viewa_9"))

# HBV genome parameters
GENOME_LENGTH = 3221
DEFAULT_OFFSET = 1823  # DR1 position corresponds to original position 1824

os.makedirs(OUTPUT_DIR, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(os.path.join(OUTPUT_DIR, "integrated_analysis.log")),
    ],
)
logger = logging.getLogger(__name__)


def load_reference_sequence(ref_file):
    """Load reference sequence file, return (seq_id, sequence)"""
    logger.info(f"Loading reference: {ref_file}")
    try:
        for record in SeqIO.parse(ref_file, "fasta"):
            ref_id = record.id
            ref_seq = str(record.seq).upper()
            logger.info(f"Loaded reference: {ref_id}, length: {len(ref_seq)}")
            return ref_id, ref_seq
        logger.error("Reference file is empty")
        return None, None
    except Exception as e:
        logger.error(f"Failed to load reference: {e}")
        return None, None


def get_sample_list():
    """Get list of all sample IDs"""
    samples = []
    consensus_path = Path(CONSENSUS_DIR)

    if not consensus_path.exists():
        logger.error(f"Consensus directory not found: {CONSENSUS_DIR}")
        return samples

    for sample_dir in consensus_path.iterdir():
        if sample_dir.is_dir() and sample_dir.name.isdigit():
            consensus_file = sample_dir / "consensus.fasta"
            if consensus_file.exists():
                samples.append(sample_dir.name)

    samples.sort()
    logger.info(f"Found {len(samples)} samples")
    return samples


def load_consensus_sequence(sample_id):
    """Load consensus sequence for a sample"""
    consensus_file = Path(CONSENSUS_DIR) / sample_id / "consensus.fasta"

    if not consensus_file.exists():
        logger.error(f"Sample {sample_id}: consensus file not found: {consensus_file}")
        return None, None

    try:
        for record in SeqIO.parse(consensus_file, "fasta"):
            return record.id, str(record.seq).upper()
    except Exception as e:
        logger.error(f"Sample {sample_id}: failed to read consensus: {e}")
        return None, None


def transform_coordinate(current_pos, offset=DEFAULT_OFFSET, genome_length=GENOME_LENGTH):
    """Transform coordinate from rotated to original coordinate system"""
    original_pos = ((current_pos - 1) + offset) % genome_length + 1
    return original_pos


def compare_sequences_with_transform(seq1, seq2, offset=DEFAULT_OFFSET):
    """
    Compare sequences and transform coordinates.

    Returns list of differences with both rotated and original positions.
    """
    differences = []

    ref_len = len(seq1)
    con_len = len(seq2)
    max_len = max(ref_len, con_len)

    # Pad shorter sequence
    if ref_len < max_len:
        seq1 += "N" * (max_len - ref_len)
    if con_len < max_len:
        seq2 += "N" * (max_len - con_len)

    i = 0
    while i < max_len:
        if seq1[i] != seq2[i]:
            rotated_pos = i + 1
            original_pos = transform_coordinate(rotated_pos, offset)

            if seq1[i] != "-" and seq2[i] != "-":
                differences.append(
                    {
                        "rotated_position": rotated_pos,
                        "original_position": original_pos,
                        "type": "SNP",
                        "reference": seq1[i],
                        "consensus": seq2[i],
                        "length": 1,
                    }
                )
            else:
                if seq1[i] == "-":
                    differences.append(
                        {
                            "rotated_position": rotated_pos,
                            "original_position": original_pos,
                            "type": "insertion",
                            "reference": "-",
                            "consensus": seq2[i],
                            "length": 1,
                        }
                    )
                elif seq2[i] == "-":
                    differences.append(
                        {
                            "rotated_position": rotated_pos,
                            "original_position": original_pos,
                            "type": "deletion",
                            "reference": seq1[i],
                            "consensus": "-",
                            "length": 1,
                        }
                    )
        i += 1

    return differences


def process_single_sample(sample_id, ref_id, reference_seq, output_format="csv"):
    """Process single sample comparison and coordinate transform"""
    logger.info(f"Processing sample: {sample_id}")

    consensus_header, consensus_seq = load_consensus_sequence(sample_id)

    if not consensus_header or not consensus_seq:
        logger.error(f"Sample {sample_id}: cannot read consensus")
        return None

    logger.info(
        f"Sample {sample_id}: ref={ref_id} (len={len(reference_seq)}), consensus len={len(consensus_seq)}"
    )

    differences = compare_sequences_with_transform(reference_seq, consensus_seq, DEFAULT_OFFSET)

    result = {
        "sample_id": sample_id,
        "reference_id": ref_id,
        "reference_length": len(reference_seq),
        "consensus_length": len(consensus_seq),
        "offset_used": DEFAULT_OFFSET,
        "total_differences": len(differences),
        "differences": differences,
        "processing_time": datetime.now().isoformat(),
    }

    output_file = Path(OUTPUT_DIR) / f"{sample_id}_comparison_transformed.{output_format}"
    save_result(result, output_file, output_format)

    logger.info(f"Sample {sample_id}: {len(differences)} differences, saved to {output_file}")
    return result


def save_result(result, output_file, output_format="csv"):
    """Save comparison result to file"""
    try:
        if output_format == "csv":
            with open(output_file, "w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                writer.writerow(["pos", "mutat", "ref", "alt", "length"])
                for diff in result["differences"]:
                    writer.writerow(
                        [
                            diff["original_position"],
                            diff["type"],
                            diff["reference"],
                            diff["consensus"],
                            diff["length"],
                        ]
                    )

        elif output_format == "tsv":
            with open(output_file, "w", encoding="utf-8") as f:
                f.write(f"# Sample Comparison Result\n")
                f.write(f"sample_id\t{result['sample_id']}\n")
                f.write(f"reference_id\t{result['reference_id']}\n")
                f.write(f"reference_length\t{result['reference_length']}\n")
                f.write(f"consensus_length\t{result['consensus_length']}\n")
                f.write(f"offset\t{result['offset_used']}\n")
                f.write(f"total_differences\t{result['total_differences']}\n")
                f.write(f"processing_time\t{result['processing_time']}\n\n")

                f.write("original_pos\trotated_pos\ttype\tref\talt\tlength\n")
                for diff in result["differences"]:
                    f.write(
                        f"{diff['original_position']}\t{diff['rotated_position']}\t{diff['type']}\t{diff['reference']}\t{diff['consensus']}\t{diff['length']}\n"
                    )

        elif output_format == "json":
            with open(output_file.with_suffix(".json"), "w", encoding="utf-8") as f:
                json.dump(result, f, indent=2, ensure_ascii=False)

    except Exception as e:
        logger.error(f"Failed to save result: {e}")


def generate_summary_report(results, output_format="csv"):
    """Generate summary report"""
    logger.info("Generating summary report...")

    summary_file = Path(OUTPUT_DIR) / f"integrated_summary.{output_format}"

    try:
        if output_format == "csv":
            with open(summary_file, "w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                writer.writerow(
                    [
                        "sample_id",
                        "reference_id",
                        "ref_length",
                        "consensus_length",
                        "offset",
                        "differences",
                        "time",
                    ]
                )
                for result in results:
                    if result:
                        writer.writerow(
                            [
                                result["sample_id"],
                                result["reference_id"],
                                result["reference_length"],
                                result["consensus_length"],
                                result["offset_used"],
                                result["total_differences"],
                                result["processing_time"],
                            ]
                        )

        elif output_format == "tsv":
            with open(summary_file, "w", encoding="utf-8") as f:
                f.write(
                    "sample_id\treference_id\tref_length\tconsensus_length\toffset\tdifferences\ttime\n"
                )
                for result in results:
                    if result:
                        f.write(
                            f"{result['sample_id']}\t{result['reference_id']}\t{result['reference_length']}\t{result['consensus_length']}\t{result['offset_used']}\t{result['total_differences']}\t{result['processing_time']}\n"
                        )

        logger.info(f"Summary saved to: {summary_file}")
        generate_variant_distribution_report(results)

    except Exception as e:
        logger.error(f"Failed to generate summary: {e}")


def generate_variant_distribution_report(results):
    """Generate variant position distribution report"""
    try:
        position_counts = {}
        variant_types = {"SNP": 0, "insertion": 0, "deletion": 0}

        for result in results:
            if result:
                for diff in result["differences"]:
                    pos = diff["original_position"]
                    if pos not in position_counts:
                        position_counts[pos] = {"count": 0, "samples": []}
                    position_counts[pos]["count"] += 1
                    position_counts[pos]["samples"].append(result["sample_id"])
                    variant_types[diff["type"]] += 1

        dist_file = Path(OUTPUT_DIR) / "variant_position_distribution.txt"
        with open(dist_file, "w", encoding="utf-8") as f:
            f.write("# Variant Position Distribution Report\n")
            f.write(f"# Generated: {datetime.now().isoformat()}\n\n")

            f.write("## Variant Type Statistics\n")
            for vtype, count in variant_types.items():
                f.write(f"{vtype}: {count}\n")
            f.write(f"Total: {sum(variant_types.values())}\n\n")

            f.write("## High-frequency Variant Positions (>=5 samples)\n")
            f.write("position\tcount\tsamples\n")

            sorted_positions = sorted(
                position_counts.items(), key=lambda x: x[1]["count"], reverse=True
            )

            for pos, data in sorted_positions:
                if data["count"] >= 5:
                    sample_list = ",".join(data["samples"][:10])
                    if len(data["samples"]) > 10:
                        sample_list += f"...({len(data['samples'])} total)"
                    f.write(f"{pos}\t{data['count']}\t{sample_list}\n")

        logger.info(f"Distribution report saved to: {dist_file}")

    except Exception as e:
        logger.error(f"Failed to generate distribution report: {e}")


def main():
    global DEFAULT_OFFSET

    parser = argparse.ArgumentParser(
        description="Consensus vs Reference Comparison with Coordinate Transform"
    )
    parser.add_argument("--test", type=str, help="Test single sample (sample ID)")
    parser.add_argument("--jobs", type=int, default=4, help="Parallel workers (default: 4)")
    parser.add_argument(
        "--output-format", choices=["csv", "tsv", "json"], default="csv", help="Output format"
    )
    parser.add_argument("--list", action="store_true", help="List all available samples")
    parser.add_argument(
        "--offset",
        type=int,
        default=DEFAULT_OFFSET,
        help=f"Coordinate offset (default: {DEFAULT_OFFSET})",
    )
    parser.add_argument("--no-confirm", action="store_true", help="Skip confirmation prompt")

    args = parser.parse_args()

    if args.offset != DEFAULT_OFFSET:
        DEFAULT_OFFSET = args.offset
        logger.info(f"Using custom offset: {DEFAULT_OFFSET}")

    logger.info("=" * 60)
    logger.info("Consensus vs Reference Comparison v1.1")
    logger.info(f"Consensus dir: {CONSENSUS_DIR}")
    logger.info(f"Reference: {REFERENCE_FILE}")
    logger.info(f"Offset: {DEFAULT_OFFSET}")
    logger.info("=" * 60)

    samples = get_sample_list()

    if args.list:
        print(f"Found {len(samples)} samples:")
        for i, sample in enumerate(samples, 1):
            print(f"{i:3d}. {sample}")
        return

    if not samples:
        logger.error("No samples found")
        return

    ref_id, reference_seq = load_reference_sequence(REFERENCE_FILE)
    if not ref_id or not reference_seq:
        logger.error("Cannot load reference, exiting")
        return

    if args.test:
        if args.test not in samples:
            logger.error(f"Sample {args.test} not found")
            return

        result = process_single_sample(args.test, ref_id, reference_seq, args.output_format)
        if result:
            print(f"\nSample {args.test} complete:")
            print(f"  Reference: {result['reference_id']}")
            print(f"  Ref length: {result['reference_length']}")
            print(f"  Consensus length: {result['consensus_length']}")
            print(f"  Offset: {result['offset_used']}")
            print(f"  Differences: {result['total_differences']}")

            if result["differences"]:
                print(f"  Examples:")
                for diff in result["differences"][:5]:
                    print(
                        f"    pos {diff['original_position']} (rotated {diff['rotated_position']}): "
                        f"{diff['type']} {diff['reference']}->{diff['consensus']}"
                    )
                if len(result["differences"]) > 5:
                    print(f"    ... {len(result['differences'])} total")
    else:
        print(f"\nProcessing {len(samples)} samples (workers: {args.jobs})")
        print(f"Reference: {ref_id}")
        print(f"Offset: {DEFAULT_OFFSET}")

        if not args.no_confirm:
            confirm = input("Continue? [y/N]: ")
            if confirm.lower() != "y":
                logger.info("User cancelled")
                return

        results = []

        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            future_to_sample = {
                executor.submit(
                    process_single_sample, sample, ref_id, reference_seq, args.output_format
                ): sample
                for sample in samples
            }

            completed = 0
            for future in as_completed(future_to_sample):
                sample = future_to_sample[future]
                try:
                    result = future.result()
                    results.append(result)
                    completed += 1
                    logger.info(
                        f"Progress: {completed}/{len(samples)} ({100*completed/len(samples):.1f}%)"
                    )
                except Exception as e:
                    logger.error(f"Error processing {sample}: {e}")

        generate_summary_report(results, args.output_format)

        successful_results = [r for r in results if r is not None]
        total_differences = sum(r["total_differences"] for r in successful_results)

        print(f"\nComplete!")
        print(f"  Processed: {len(successful_results)}/{len(samples)} samples")
        print(f"  Total differences: {total_differences}")
        print(
            f"  Average: {total_differences/len(successful_results):.1f}"
            if successful_results
            else "  Average: 0"
        )
        print(f"  Output: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
