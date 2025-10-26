"""
Modern MHC class II binding predictor with multi-matrix support.
"""

import os
import glob
import re
from pathlib import Path
from typing import List, Tuple, Union, Optional
import numpy as np
import csv
import argparse


class HLAmatrix:
    """Class for loading and scoring a single MHC binding matrix."""

    def __init__(self, csv_file: Union[str, Path], threshold: int = 5):
        self.csv_file = Path(csv_file)
        self.threshold = threshold
        self.threshold_score = 0.0
        self.allele = self.csv_file.stem
        self.matrix_data = {}  # Dict[str, np.ndarray]
        self.load_matrix()

    def load_matrix(self) -> None:
        """Load CSV matrix into memory as NumPy arrays."""
        percentage_threshold_tmp = []
        scoring_tmp = []

        with self.csv_file.open() as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue  # skip empty lines
                parts = line.split(',')
                key = parts[0].strip()
                if not key:
                    continue  # skip lines with empty first column

                # Rows with amino acids
                if len(key) <= 2:
                    values = []
                    for x in parts[1:]:
                        x = x.strip()
                        if x == '-' or x == '':
                            values.append(0.0)
                        else:
                            values.append(float(x))
                    self.matrix_data[key[0]] = np.array(values)
                elif key == 'Percent Threshold':
                    percentage_threshold_tmp = parts[1:]
                elif key == 'Numerical Score':
                    scoring_tmp = parts[1:]

        # Determine threshold score
        for i, threshold_str in enumerate(percentage_threshold_tmp):
            if int(threshold_str.strip('%')) == self.threshold:
                self.threshold_score = float(scoring_tmp[i])



    def score_nonamer(self, nonamer: str) -> float:
        """Score a single 9-mer."""
        score = 0.0
        for i, residue in enumerate(nonamer):
            score += self.matrix_data.get(residue, np.zeros(9))[i]
        return score

    def score_sequence(self, sequence: str) -> List[Tuple[str, float, int, int]]:
        """Score a full sequence and return matching epitopes above threshold."""
        results = []
        for i in range(len(sequence) - 8):
            nonamer = sequence[i:i + 9]
            score = self.score_nonamer(nonamer)
            if score >= self.threshold_score:
                results.append((nonamer, score, i, i + 8))
        return sorted(results, key=lambda x: x[1], reverse=True)

    @staticmethod
    def write_results_csv(results: List[Tuple], filename: str, header: Optional[Tuple] = None) -> None:
        """Write results to CSV."""
        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            if header:
                writer.writerow(header)
            writer.writerows(results)


class HLAmatrixCollection:
    """Manage multiple HLAmatrix objects and score sequences against all of them."""

    def __init__(self, matrix_dir: Union[str, Path], threshold: int = 5):
        self.matrix_dir = Path(matrix_dir)
        self.threshold = threshold
        self.matrices: List[HLAmatrix] = []
        self.load_all_matrices()

    def load_all_matrices(self) -> None:
        """Load all CSV matrices in the directory."""
        for matrix_file in self.matrix_dir.glob("*.csv"):
            self.matrices.append(HLAmatrix(matrix_file, self.threshold))

    def score_sequence(self, sequence: str, base_index: int = 1) -> List[Tuple[float, int, int, str, str]]:
        """
        Score a sequence against all matrices.
        Returns list of tuples: (score, start, end, epitope, allele)
        """
        results_set = set()
        for matrix in self.matrices:
            for epitope, score, begin, end in matrix.score_sequence(sequence):
                results_set.add((score, begin + base_index, end + base_index, epitope, matrix.allele))
        return sorted(results_set, key=lambda x: x[0], reverse=True)


def parse_sequence(args) -> str:
    """Extract the peptide sequence either from --sequence or from a file."""
    if args.sequence and args.filename:
        raise ValueError('Specify only one of --sequence or -f/--file.')
    if not args.sequence and not args.filename:
        raise ValueError('You must provide a sequence or a file.')

    if args.sequence:
        # Literal peptide sequence
        return re.sub(r'\s', '', args.sequence)
    else:
        # Read sequence from file
        seq_file = Path(args.filename)
        if not seq_file.is_file():
            raise FileNotFoundError(f"Sequence file not found: {seq_file}")
        return re.sub(r'\s', '', seq_file.read_text())




def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-t', "--threshold", type=int, default=5)
    parser.add_argument('-c', '--csv', dest='csv_filename', default=None, help='CSV output filename')
    parser.add_argument('--no-header', dest='csv_header', action='store_false', help='Omit header row in CSV')
    parser.add_argument('-f', '--file', dest='filename', default=None, help='File with sequence')
    parser.add_argument('-s', '--sequence', dest='sequence', default=None, help='Peptide sequence string')
    parser.add_argument('-i', '--starting-index', dest='base_index', default=1, type=int,
                    help='Index of first peptide in sequence')
    parser.add_argument('-m', '--matrices', dest='matrix_dir', default='matrices', help='Directory with matrices')

    args = parser.parse_args()

    sequence = parse_sequence(args)

    collection = HLAmatrixCollection(args.matrix_dir, args.threshold)
    results = collection.score_sequence(sequence, base_index=args.base_index)

    if args.csv_filename:
        header = ('score', 'begin', 'end', 'epitope', 'allele') if args.csv_header else None
        HLAmatrix.write_results_csv(results, args.csv_filename, header)
    else:
        for row in results:
            print(f"{row[0]:.2f}\t{row[1]}\t{row[2]}\t{row[3]}\t{row[4]}")


if __name__ == '__main__':
    main()
