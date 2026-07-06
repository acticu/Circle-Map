import subprocess
import pysam
import pandas as pd
import os
from typing import Tuple, List, Dict
from config import (
    BWA_MEM2_PATH, SAMBLASTER_PATH, SAMTOOLS_PATH,
    MIN_MAPQ, MIN_CLIP_LEN, MAX_INSERT_SIZE, OUTPUT_DIR
)

class ReadExtractor:
    def __init__(self, ref_fasta: str, fastq_r1: str, fastq_r2: str):
        self.ref = ref_fasta
        self.fastq_r1 = fastq_r1
        self.fastq_r2 = fastq_r2
        self.bam_path = os.path.join(OUTPUT_DIR, "aligned_sorted.bam")
        
    def align_and_extract(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Pipeline:
        1. bwa-mem2 (Paired-End mode)
        2. samtools sort (Name sort for samblaster)
        3. samblaster (Split discordant/split)
        4. pysam parsing (In-memory, no bedtools)
        """
        print(">> Step 1: Aligning Paired-End reads with bwa-mem2...")
        # BWA-MEM2 accepts two FASTQ files for PE data
        bwa_cmd = [
            BWA_MEM2_PATH, "mem", "-t", "8", "-M", "-Y", "-k", "19",
            self.ref, self.fastq_r1, self.fastq_r2
        ]
        
        print(">> Step 2: Sorting and Extracting with samblaster...")
        split_out = os.path.join(OUTPUT_DIR, "split_reads.bam")
        discord_out = os.path.join(OUTPUT_DIR, "discordant_reads.bam")
        
        # Pipe: bwa -> samtools sort (name) -> samblaster
        full_cmd = (
            f"{' '.join(bwa_cmd)} | "
            f"{SAMTOOLS_PATH} sort -n -@ 8 - | "
            f"{SAMBLASTER_PATH} "
            f"--splitFile {split_out} "
            f"--discordantFile {discord_out} "
            f"-o /dev/null"
        )
        
        # Execute
        subprocess.run(full_cmd, shell=True, check=True, executable='/bin/bash')
        
        print(">> Step 3: Parsing Split Reads (Junction Candidates)...")
        split_df = self._parse_split_reads(split_out)
        
        print(">> Step 4: Parsing Discordant Reads (Supporters)...")
        discord_df = self._parse_discordant_reads(discord_out)
        
        return split_df, discord_df

    def _parse_split_reads(self, bam_path: str) -> pd.DataFrame:
        """
        Parses BAM to find soft-clipped reads.
        Preserves Sequence and Quality for Matrix Realignment.
        """
        records = []
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if read.is_unmapped or read.mapping_quality < MIN_MAPQ:
                    continue
                
                clips = self._get_soft_clips(read)
                if not clips:
                    continue
                
                for clip_pos, clip_seq, clip_qual, clip_type in clips:
                    if len(clip_seq) < MIN_CLIP_LEN:
                        continue
                        
                    records.append({
                        'read_name': read.query_name,
                        'chr': read.reference_name,
                        'pos': clip_pos,
                        'strand': '-' if read.is_reverse else '+',
                        'seq': clip_seq,
                        'qual': clip_qual,
                        'type': clip_type,
                        'mate_chr': read.next_reference_name,
                        'mate_pos': read.next_reference_start,
                        'is_proper_pair': not read.is_paired or (read.is_proper_pair)
                    })
        
        df = pd.DataFrame(records)
        if df.empty:
            return df
            
        # Tag same-chromosome mates for single-segment hypothesis
        df['is_same_chr'] = df['chr'] == df['mate_chr']
        return df

    def _parse_discordant_reads(self, bam_path: str) -> pd.DataFrame:
        """
        Parses discordant pairs that support long-range connections.
        """
        records = []
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if read.is_unmapped or read.mate_is_unmapped:
                    continue
                if read.reference_name == read.next_reference_name:
                    dist = abs(read.reference_start - read.next_reference_start)
                    if dist > MAX_INSERT_SIZE:
                         records.append({
                            'read_name': read.query_name,
                            'chr1': read.reference_name,
                            'pos1': read.reference_start,
                            'chr2': read.next_reference_name,
                            'pos2': read.next_reference_start,
                            'strand1': '-' if read.is_reverse else '+',
                            'strand2': '-' if read.mate_is_reverse else '+'
                        })
                else:
                    # Inter-chromosomal also counts as discordant
                    records.append({
                            'read_name': read.query_name,
                            'chr1': read.reference_name,
                            'pos1': read.reference_start,
                            'chr2': read.next_reference_name,
                            'pos2': read.next_reference_start,
                            'strand1': '-' if read.is_reverse else '+',
                            'strand2': '-' if read.mate_is_reverse else '+'
                        })
        return pd.DataFrame(records)

    def _get_soft_clips(self, read) -> List[Tuple[int, str, str, str]]:
        """Extracts soft-clipped sequences and qualities."""
        clips = []
        cigartypes = dict(zip(range(8), "MIDNSHP=X"))
        seq = read.query_sequence
        qual = read.query_qualities
        
        if seq is None or qual is None:
            return []
            
        pos = read.reference_start
        s_pos = 0 
        
        for op, length in read.cigartuples:
            op_char = cigartypes[op]
            if op_char == 'S':
                clip_seq = seq[s_pos:s_pos+length]
                # Convert quality ints to Phred string
                clip_qual = "".join([chr(q+33) for q in qual[s_pos:s_pos+length]])
                
                # Determine coordinate: 
                # Head clip: position is start of read (reference start if mapped, else ambiguous)
                # Tail clip: position is end of aligned block
                if s_pos == 0: # Head
                    coord = pos # Approximate junction at start
                    clip_type = 'head'
                else: # Tail
                    coord = pos # Junction at end of alignment
                    clip_type = 'tail'
                
                clips.append((coord, clip_seq, clip_qual, clip_type))
            
            # Update positions based on CIGAR
            if op_char in "M=X":
                pos += length
            if op_char in "MIS=X":
                s_pos += length
                
        return clips