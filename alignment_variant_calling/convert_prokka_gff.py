#!/usr/bin/env python
# genbank created GFF files are not compatible with gff3 parser - roll our own

import sys
from collections import defaultdict

reserved_keys = ['ID', 'Name', 'biotype', 'Parent']

gene_counts = defaultdict(int)
trns_counts = defaultdict(int)


def parse_gff_annots(s):
    return dict(x.split('=') for x in s.split(';'))
'''
2019-nCoV_WH01	Prodigal:002006	CDS	27177	27362	.	+	0	ID=MHEPPOEF_00135;inference=ab initio prediction:Prodigal:002006;locus_tag=MHEPPOEF_00135;product=hypothetical protein
###
gnl|charite|EBIIEENH_1	Prodigal:002006	CDS	27187	27372	.	+	0	ID=EBIIEENH_00007;inference=ab initio prediction:Prodigal:002006;locus_tag=EBIIEENH_00007;product=hypothetical protein;protein_id=gnl|charite|EBIIEENH_00007
###
GWHABKP00000001	Prodigal:002006	CDS	27169	27354	.	+	0	ID=IINNKKGC_00483;inference=ab initio prediction:Prodigal:002006;locus_tag=IINNKKGC_00483;product=hypothetical protein
###
MG772933	Prodigal:002006	CDS	27107	27292	.	+	0	ID=IINNKKGC_01234;inference=ab initio prediction:Prodigal:002006;locus_tag=IINNKKGC_01234;product=hypothetical protein

'''

def parse_line(l):
    cols = l.split('\t')
    if cols[2] == 'CDS':
        info = parse_gff_annots(cols[8])
        output_gene_line(cols, info)
        transcript = output_transcript_line(cols, info)
        output_5utr_line(cols, info, transcript)
        output_exon_line(cols, info, transcript)
        output_cds_line(cols, info, transcript)
        output_3utr_line(cols, info, transcript)
    elif cols[2] == 'stem_loop':
        cols[8] += ';biotype=misc_RNA'
        print("\t".join(cols))
    elif cols[2] == 'five_prime_UTR' or cols[2] == 'three_prime_UTR':
        return
    else:
        print("\t".join(cols))


def output_new_cols(feature, cols, new_info, old_info):
    new_info.extend("{}={}".format(k, v) for k, v in old_info.items() if k not
                    in reserved_keys)
    print("\t".join(cols[:2] + [feature] + cols[3:8] + [";".join(new_info)]))


def output_gene_line(cols, info):
    outinfo = ["ID=gene:{}".format(info['ID']),
               "Name={}".format(info['locus_tag']),
               "biotype=protein_coding"]
    output_new_cols('gene', cols, outinfo, info)


def output_transcript_line(cols, info):
    trns_counts[info['ID']] += 1
    transcript = "{}-{}".format(info['ID'], trns_counts[info['ID']])
    outinfo = ["ID=transcript:{}".format(transcript),
               "Parent=gene:{}".format(info['ID']),
               "Name={}".format(info['locus_tag']),
               "biotype=protein_coding"]
    output_new_cols('mRNA', cols, outinfo, info)
    return transcript


def output_cds_line(cols, info, parent):
    outinfo = ["ID=CDS:{}_cds".format(info['ID']),
               "Parent=transcript:{}".format(parent),
               "protein_id={}_protein".format(info['ID'])]
    output_new_cols('CDS', cols, outinfo, info)


def output_exon_line(cols, info, parent):
    outinfo = ["ID={}-exon_1".format(parent),
               "Parent=transcript:{}".format(parent),
               "Name={}".format(info['locus_tag'] + "_exon"),
               "rank=1"]
    output_new_cols('exon', cols, outinfo, info)


def output_3utr_line(cols, info, parent):
    outinfo = ["Parent=transcript:{}".format(parent)]
    start = str(int(cols[4]) + 1)
    end = str(int(cols[4]) + 2)
    adjusted_cols = cols[:3] + [start, end] + cols[5:]
    output_new_cols('three_prime_UTR', adjusted_cols, outinfo, info)


def output_5utr_line(cols, info, parent):
    outinfo = ["Parent=transcript:{}".format(parent)]
    start = str(int(cols[3]) - 2)
    end = str(int(cols[3]) - 1)
    adjusted_cols = cols[:3] + [start, end] + cols[5:]
    output_new_cols('five_prime_UTR', adjusted_cols, outinfo, info)


def main(f):
    with open(f, 'rt') as infile:
        for line in infile:
            line = line.rstrip()
            if line.startswith('#'):
                print(line)
            elif not line.rstrip():
                continue
            else:
                parse_line(line)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit("Usage: {} input.gff".format(sys.argv[0]))
    main(sys.argv[1])
