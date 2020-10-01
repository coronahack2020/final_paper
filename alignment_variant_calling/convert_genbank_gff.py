#!/usr/bin/env python
# genbank created GFF files are not compatible with gff3 parser - roll our own

import sys
from collections import defaultdict

reserved_keys = ['ID', 'Name', 'biotype', 'Parent']
skip_categories = ['five_prime_UTR', 'three_prime_UTR', 'region', 'exon']

gene_counts = defaultdict(int)
trns_counts = defaultdict(int)


def parse_gff_annots(s):
    return dict(x.split('=') for x in s.split(';'))


def parse_line(l):
    cols = l.split('\t')
    if cols[2] == 'CDS':
        info = parse_gff_annots(cols[8])
        gene_id = output_gene_line(cols, info)
        transcript = output_transcript_line(cols, info, gene_id)
        output_5utr_line(cols, info, transcript)
        output_exon_line(cols, info, transcript)
        output_cds_line(cols, info, transcript)
        output_3utr_line(cols, info, transcript)
    elif cols[2] == 'stem_loop':
        cols[8] += ';biotype=misc_RNA'
        print("\t".join(cols))
    elif cols[2] in skip_categories:
        return
    else:
        print("\t".join(cols))


def output_new_cols(feature, cols, new_info, old_info, phase='.'):
    new_info.extend("{}={}".format(k, v) for k, v in old_info.items() if k not
                    in reserved_keys)
    print("\t".join(cols[:2] + [feature] + cols[3:7] +
                    [phase, ";".join(new_info)]))


def output_gene_line(cols, info):
    try:
        gene_info = dict(x.split(':') for x in info['Dbxref'].split(','))
        gene_id = gene_info['GeneID']
    except KeyError:
        gene_id = info['ID']
    gene_counts[gene_id] += 1
    gene_id += "_{}".format(gene_counts[gene_id])
    outinfo = ["ID=gene:{}".format(gene_id),
               "Name={}-{}".format(gene_id, info['Name']),
               "biotype=protein_coding"]
    output_new_cols('gene', cols, outinfo, info)
    return gene_id


def output_transcript_line(cols, info, gene_id):
    # TODO - grok whether line comes from website downloaded GFF or converted
    # for now assume comes from bp_genbank2gff3
    #try:
    #    transcript = parent2transcript(info['Parent'])
    #except KeyError:
    #    transcript = gene2transcript(info['ID'])
    transcript = gene2transcript(info['ID'])
    outinfo = ["ID={}".format(transcript),
               "Parent=gene:{}".format(gene_id),
               "Name={}".format(info['Name']),
               "biotype=protein_coding"]
    output_new_cols('mRNA', cols, outinfo, info)
    return transcript


def output_cds_line(cols, info, parent):
    outinfo = ["ID=CDS:{}_cds".format(info['Name']),
               "Parent={}".format(parent),
               "protein_id={}_protein".format(info['Name'])]
    output_new_cols('CDS', cols, outinfo, info, '0')


def output_exon_line(cols, info, parent):
    outinfo = ["ID={}-exon_1".format(parent),
               "Parent={}".format(parent),
               "Name={}".format(info['Name'] + "_exon"),
               "rank=1"]
    output_new_cols('exon', cols, outinfo, info)


def output_3utr_line(cols, info, parent):
    outinfo = ["Parent={}".format(parent)]
    start = str(int(cols[4]) + 1)
    end = str(int(cols[4]) + 2)
    adjusted_cols = cols[:3] + [start, end] + cols[5:]
    output_new_cols('three_prime_UTR', adjusted_cols, outinfo, info)


def output_5utr_line(cols, info, parent):
    outinfo = ["Parent={}".format(parent)]
    start = str(int(cols[3]) - 2)
    end = str(int(cols[3]) - 1)
    adjusted_cols = cols[:3] + [start, end] + cols[5:]
    output_new_cols('five_prime_UTR', adjusted_cols, outinfo, info)


def parent2transcript(s):
    t = s.replace('-', ':').replace('gene', 'transcript')
    trns_counts[t] += 1
    return t + "_{}".format(trns_counts[t])


def gene2transcript(s):
    t = "transcript:{}".format(s)
    trns_counts[t] += 1
    return t + "_{}".format(trns_counts[t])

def check_file_format(f):
    features = ['gene', 'mRNA', 'exon', 'CDS']
    is_ensembl = False
    with open(f, 'rt') as infile:
        for line in infile:
            line = line.rstrip()
            if line.startswith('#'):
                continue
            cols = line.split()
            if cols[2] in features:
                is_ensembl = cols[1] == 'ensembl'
    return is_ensembl


def main(f):
    is_ensembl = check_file_format(f)
    with open(f, 'rt') as infile:
        for line in infile:
            line = line.rstrip()
            if line.startswith('#'):
                print(line)
            elif not line.rstrip():
                continue
            elif is_ensembl:
                print(line)
            else:
                parse_line(line)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit("Usage: {} input.gff".format(sys.argv[0]))
    main(sys.argv[1])
