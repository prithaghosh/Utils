#! /usr/bin/python3.8

"""
Author: Pritha Ghosh
Date: 30 August 2020
Place: Warsaw, Poland

Script for parsing BED files
Usage:
python3 bedfileParser.py --bed bed_file --out output_file

For more help use:
./bedfileParser.py -h
"""

import argparse

des = "Script for parsing BED files"
bed_help_text = "BED file format:\n>Multi-line BED\nFor more details refer to: https://genome.ucsc.edu/FAQ/FAQformat.html#format1\n\nFollow these steps for generating BED files:\n\n1. Align FASTA sequences to the genome using popular tools like BLAT or Gmap\n\n2. In case of BLAT, the output format is psl, which can be converted to BED using psl2bed\n(https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/psl2bed.html)\n\n3. In case of Gmap, the output format is BAM, which can be converted to BED using bamtobed in the bedtools package\n(https://bedtools.readthedocs.io/en/latest/content/tools/bamtobed.html)\n\n"
output_help_text = "..."

parser = argparse.ArgumentParser(
    description=f"{des}\n\n{bed_help_text}\n{output_help_text}",
    formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument("-b",
                    "--bed",
                    help="Name of the bed file")
parser.add_argument("-o",
                    "--out",
                    help="Name of the output file")
args = parser.parse_args()

class BedParse: # defining the BedParse class
    """
    This class reads the contents of the input BED file.
    """
    def __init__(self, line): # initialising the BedParse class
        self.line = line
        operation = self.line.strip().split("\t")
        self.start = int(operation[1])
        self.end = int(operation[2])
        self.name = operation[3]

        self.strand = operation[5]

        self.block_sizes = operation[10].split(",")
        self.block_sizes.pop()
        self.block_sizes = [int(x) for x in self.block_sizes]

        self.block_starts = operation[11].split(",")
        self.block_starts.pop()
        self.block_starts = [int(x) for x in self.block_starts]

        self.start_norm = None
        self.end_norm = None
    
    def end_coords(self):
        """
        This function calculates the end coordinates of blocks,
        based on block start coordinates and block sizes.
        """
        block_end_coords = []
        for start, length in zip(self.block_starts, self.block_sizes):
            end = start + length - 1
            block_end_coords.append(end)
        return block_end_coords

def extract_boundaries(bed_file):
    """
    This function stores the details for each RNA, and also
    calculates the global start and end coordinates.
    """
    my_dict = {}
    all_coordinates = []
    with open(bed_file) as f:
        for line in f:
            bed = BedParse(line)
            my_dict[bed.name] = bed

            all_coordinates.append(bed.start)
            all_coordinates.append(bed.end)

        first = min(all_coordinates)
        last = max(all_coordinates)

    return my_dict, first, last

boundaries, global_start_coord, global_end_coord = extract_boundaries(args.bed)

block_end_coords = []
full_string_index = ""
coord_dict = {}

with open(args.out, "w") as f:
    for lncrna in sorted(boundaries.keys()):
        isoform = boundaries[lncrna]

        isoform.start_norm = isoform.start - global_start_coord
        isoform.end_norm = isoform.start_norm + isoform.end - isoform.start

        for i in range(len(isoform.block_sizes)):
            if i == 0:
                block_1 = "*" * isoform.block_sizes[0]
            elif i > 0:
                gap_index = "_" * (isoform.block_starts[i] - isoform.end_coords()[i-1] - 1)
                b_start = isoform.block_starts[i]
                b_end = b_start + isoform.block_sizes[i]
                block_index = "*" * (b_end - b_start)
                full_string_index += "".join([gap_index, block_index])

        gap_start = "_" * isoform.start_norm
        gap_end = "_" * (int(global_end_coord) - isoform.end)
        
        full_string = "".join([gap_start, block_1, full_string_index, gap_end])
        print(f">{lncrna}\n{full_string}", file=f)

        block_end_coords = []
        full_string_index = ""