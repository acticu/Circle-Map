import argparse
import sys
from extractor import ReadExtractor
from realigner import ProbabilisticRealigner
from graph_solver import EccDNAGraphSolver
from config import OUTPUT_DIR
import pandas as pd

def run_pipeline(ref: str, fastq_r1: str, fastq_r2: str):
    print(f"=== Starting EccDNA Detection Pipeline (Paired-End) ===")
    print(f"Reference: {ref}")
    print(f"Read 1: {fastq_r1}")
    print(f"Read 2: {fastq_r2}")
    
    # 1. Extraction (Handles PE alignment internally)
    extractor = ReadExtractor(ref, fastq_r1, fastq_r2)
    split_df, discord_df = extractor.align_and_extract()
    
    if split_df.empty:
        print("No split reads found. Exiting.")
        return

    # 2. Probabilistic Realignment (Matrix-Based HMM)
    print(">> Running Matrix-Based Probabilistic Realignment (with Gap Penalties)...")
    realigner = ProbabilisticRealigner(ref)
    validated_df = realigner.realign_junctions(split_df)
    
    valid_count = validated_df['is_valid'].sum()
    print(f"Validated {valid_count} high-confidence junctions.")
    
    if valid_count == 0:
        print("No high-confidence junctions found.")
        return

    # 3. Graph Construction & Cycle Detection
    print(">> Building Breakpoint Graph and Detecting Cycles...")
    solver = EccDNAGraphSolver()
    solver.build_graph(validated_df, discord_df)
    circles = solver.detect_cycles()
    
    # 4. Output
    output_file = f"{OUTPUT_DIR}/eccdna_calls.csv"
    circles.to_csv(output_file, index=False)
    print(f"=== Detection Complete. Results saved to {output_file} ===")
    if not circles.empty:
        print(circles.to_string())
    else:
        print("No circular DNA structures detected.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="High-Performance Multi-Segment ecDNA Detector")
    parser.add_argument("-r", "--reference", required=True, help="Path to reference FASTA")
    parser.add_argument("-1", "--fastq1", required=True, help="Path to Read 1 (FASTQ)")
    parser.add_argument("-2", "--fastq2", required=True, help="Path to Read 2 (FASTQ)")
    
    args = parser.parse_args()
    run_pipeline(args.reference, args.fastq1, args.fastq2)