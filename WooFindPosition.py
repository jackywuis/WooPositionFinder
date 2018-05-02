import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='WooFindPosition - A tool to find information of a particular nucleotide '
                                             'from a gbk form file. '
                                             'Usage: python WooFindPosition.py -i position -g file/name.gbk '
                                             'Output: if position in CDS: position, nucleotide, CDS, nucleotide position'
                                             ' in CDS, amino acid residue position in CDS, codon sequence, amino acid '
                                             'residue, gene name, and annotation; '
                                             'if postion in RNA: position, RNA type, nucleotide position in RNA/RNA '
                                             'length, RNA name, annotation. '
                                             'otherwise: postion, nucleotide, distance to upflow gene/distance to '
                                             'downflow gene, upflow gene name/downflow gene name, upflow gene '
                                             'annotation/downflow gene annotation. '
                                             'Copyleft: Jacky Woo from ZHK Research team, iSynBio, SIAT.')
parser.add_argument('-i', '--input', type=int, help='Input nucleotide position you want to search for')
parser.add_argument('-g', '--gbk', type=str, help='Input gbk file path for searching')
parser.add_argument('-r', '--record', type=str, required=False, help='Input your chromosome name if specific')
parser.add_argument('-a', '--allele', type=str, required=False, help='Input an allele nt for checking')

args = parser.parse_args()
item = int(args.input)
readic = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '[': ']', ']': '['}
aad = {'G': 'Gly', 'A': 'Ala', 'V': 'Val', 'L': 'Leu', 'I': 'Ile', 'F': 'Phe', 'W': 'Trp', 'Y': 'Tyr', 'D': 'Asp',
       'N': 'Asn', 'E': 'Glu', 'K': 'Lys', 'Q': 'Gln', 'M': 'Met', 'S': 'Ser', 'T': 'Thr', 'C': 'Cys', 'P': 'Pro',
       'H': 'His', 'R': 'Arg', '*': 'Stop'}
codons = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
          "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
          "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
          "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
          "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
          "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
          "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
          "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"}
gbk = open(args.gbk)
out_w = False
out_be = ['', 0, '', '']
for record in SeqIO.parse(gbk, 'genbank'):
    if (args.record is not None) and (record.name != args.record):
        continue
    else:
        for feature in record.features:
            if (item in feature.location) and (feature.type == 'CDS') and ('gene' in feature.qualifiers):
                if (args.allele is not None) and (args.allele != str(record.seq[item - 1])):
                    new = args.allele
                    if feature.location.strand == 1:
                        ntpos = int(item) - feature.location.start
                        aapos = (ntpos + 2) / 3
                        copos = str(record.seq[(feature.location.start + aapos * 3 - 3): (item - 1)]) + '[' + str(
                            record.seq[(item - 1): item]) + ']' + str(
                            record.seq[item: (feature.location.start + aapos * 3)])
                        n_copos = str(
                            record.seq[(feature.location.start + aapos * 3 - 3): (item - 1)]) + '[' + new + ']' + str(
                            record.seq[item: (feature.location.start + aapos * 3)])
                    else:
                        ntpos = feature.location.end + 1 - item
                        aapos = (ntpos + 2) / 3
                        o_copos = str(record.seq[(feature.location.end - aapos * 3): (item - 1)]) + '[' + str(
                            record.seq[(item - 1): item]) + ']' + str(
                            record.seq[item: (feature.location.end - aapos * 3 + 3)])
                        on_copos = str(
                            record.seq[(feature.location.end - aapos * 3): (item - 1)]) + '[' + new + ']' + str(
                            record.seq[item: (feature.location.end - aapos * 3 + 3)])
                        copos, n_copos = '', ''
                        for i in o_copos:
                            copos = readic[i] + copos
                        for i in on_copos:
                            n_copos = readic[i] + n_copos
                    if len(new) == 1:
                        outfh = str(item) + '\t' + str(
                            record.seq[item - 1]) + '->' + new + '\t' + feature.type + '\t' + str(ntpos) + '(' + str(
                            aapos) + ')\t' + copos + '->' + n_copos + '\t' + aad[
                                    codons[copos.replace('[', '').replace(']', '')]] + '->' + aad[
                                    codons[n_copos.replace('[', '').replace(']', '')]] + '\t' + \
                                feature.qualifiers['gene'][0] + '\t' + feature.qualifiers['gene'][0] + ':' + codons[
                                    copos.replace('[', '').replace(']', '')] + str(aapos) + codons[
                                    n_copos.replace('[', '').replace(']', '')] + '\t'
                    else:
                        outfh = str(item) + '\t' + str(
                            record.seq[item - 1]) + '->' + new + '\t' + feature.type + '\t' + str(ntpos) + '(' + str(
                            aapos) + ')\t' + copos + '->' + n_copos + '\t' + aad[
                                    codons[copos.replace('[', '').replace(']', '')]] + '\t' + \
                                feature.qualifiers['gene'][0] + '\t' + feature.qualifiers['gene'][0] + ':' + codons[
                                    copos.replace('[', '').replace(']', '')] + str(aapos) + 'fs\t'
                else:
                    if feature.location.strand == 1:
                        ntpos = int(item) - feature.location.start
                        aapos = (ntpos + 2) / 3
                        copos = record.seq[
                                (feature.location.start + aapos * 3 - 3): (feature.location.start + aapos * 3)]
                    else:
                        ntpos = feature.location.end + 1 - item
                        aapos = (ntpos + 2) / 3
                        o_copos = record.seq[(feature.location.end - aapos * 3): (feature.location.end - aapos * 3 + 3)]
                        copos = ''
                        for i in o_copos:
                            copos = readic[i] + copos
                    outfh = str(item) + '\t' + str(record.seq[item - 1]) + '\t' + feature.type + '\t' + str(
                        ntpos) + '(' + str(aapos) + ')' + str(copos) + '\t' + aad[codons[
                        copos]] + '\t' + feature.qualifiers['gene'][0] + '\t'
                if 'product' in feature.qualifiers:
                    outfh += feature.qualifiers['product'][0]
                elif 'function' in feature.qualifiers:
                    outfh += feature.qualifiers['function'][0]
                elif 'note' in feature.qualifiers:
                    outfh += feature.qualifiers['note'][0]
                out_w = True
                break
            elif (item in feature.location) and (feature.type not in ['source', 'gene', 'misc_feature']) and (
                    'gene' in feature.qualifiers):
                if feature.location.strand == 1:
                    ntpos = item - feature.location.start
                else:
                    ntpos = feature.location.end + 1 - item
                if (args.allele is not None) and (args.allele != str(record.seq[item - 1])):
                    outfh = str(item) + '\t' + str(
                        record.seq[item - 1]) + '->' + args.allele + '\t' + feature.type + '\t' + str(
                        ntpos) + '/' + str(len(feature.extract(record.seq))) + 'nt\t' + feature.qualifiers['gene'][
                        0] + '\t'
                else:
                    outfh = str(item) + '\t' + str(record.seq[item - 1]) + '\t' + feature.type + '\t' + str(
                        ntpos) + '/' + str(len(feature.extract(
                        record.seq))) + 'nt\t' + feature.qualifiers['gene'][0] + '\t'
                if 'product' in feature.qualifiers:
                    outfh += feature.qualifiers['product'][0]
                elif 'function' in feature.qualifiers:
                    outfh += feature.qualifiers['function'][0]
                elif 'note' in feature.qualifiers:
                    outfh += feature.qualifiers['note'][0]
                out_w = True
                break
            elif (item > feature.location.start) and (item > feature.location.end) and (
                    feature.type not in ['source', 'gene', 'misc_feature']) and ('gene' in feature.qualifiers):
                out_be[0] = feature.qualifiers['gene'][0]
                out_be[1] = item - feature.location.end
                if feature.location.strand == 1:
                    out_be[2] = '+'
                else:
                    out_be[2] = '-'
                if 'product' in feature.qualifiers:
                    out_be[3] = feature.qualifiers['product'][0]
                elif 'function' in feature.qualifiers:
                    out_be[3] = feature.qualifiers['function'][0]
                elif 'note' in feature.qualifiers:
                    out_be[3] = feature.qualifiers['note'][0]
            elif (int(item) < feature.location.start) and (item < feature.location.end) and (
                    feature.type not in ['source', 'gene', 'misc_feature']) and ('gene' in feature.qualifiers):
                if not out_w:
                    if (args.allele is not None) and (args.allele != str(record.seq[item - 1])):
                        outfh = str(item) + '\t' + str(record.seq[item - 1]) + '->' + args.allele + '\t' + out_be[
                            2] + str(out_be[1]) + '/'
                    else:
                        outfh = str(item) + '\t' + str(record.seq[item - 1]) + '\t-\t' + out_be[2] + str(out_be[
                            1]) + '/'
                    if feature.location.strand == 1:
                        outfh += '-'
                    elif feature.location.strand == -1:
                        outfh += '+'
                    outfh += str(feature.location.start - item + 1) + '\t' + out_be[0] + '/' + feature.qualifiers[
                        'gene'][0] + '\t' + out_be[3] + '/'
                    if 'product' in feature.qualifiers:
                        outfh += feature.qualifiers['product'][0]
                    elif 'function' in feature.qualifiers:
                        outfh += feature.qualifiers['function'][0]
                    elif 'note' in feature.qualifiers:
                        outfh += feature.qualifiers['note'][0]
                    break
    print outfh
gbk.close()
