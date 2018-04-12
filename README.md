# WooPositionFinder
This is a small program to find information of a particular nucleotide from a gbk form files.

Usage: ptyhon Woo PositionFinder.py -i postion -g file/name.gbk

-i: Input nucleotide position you want to search for

-g: Input gbk format file path for searching

Output: 

if position in CDS:

position; nucleotide, CDS, nucleotide position in CDS, amino acid residue position in CDS; codon sequence; amino acid, residue, gene name, and annotation;

if postion in RNA:

position, RNA type, nucleotide position in RNA/RNA, length, RNA name, annotation;

otherwise:

postion, nucleotide, distance to upflow gene/distance to downflow gene, upflow gene name/downflow gene name, upflow gene annotation/downflow gene annotation.
