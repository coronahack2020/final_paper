#!/usr/bin/env python3
import sys
import os
import pandas as pd
import re
from parse_vcf import VcfReader
from collections import defaultdict


def main(v, failures, out):
    vcf = VcfReader(v)
    csq_re = re.compile(r'''.*Format:\s*\'\[\*\](\S+\|*\S+)\[(\S+\|*\S+)\]''')
    match = csq_re.match(vcf.header.metadata['INFO']['BCSQ'][0]['Description'])
    csq_order = [x for x in match.group(1).split('|') + match.group(2).split('|') if x != '']
    df = defaultdict(list)
    for record in vcf:
        if 'BCSQ' not in record.INFO_FIELDS:
            continue
        csqs = [dict(zip(csq_order, x.split('|'))) for x in record.parsed_info_fields()['BCSQ']]
        i = 0
        for alt in record.ALLELES[1:]:
            for k in csq_order:
                for c in csqs:
                    if 'dna_change' in c:
                        if c['dna_change'] != "{}{}>{}".format(record.POS, record.REF, alt):
                            continue
                        for f in ['CHROM', 'POS', 'REF']:
                            df[f].append(getattr(record, f))
                        df['ALT'].append(alt)
                        df['VAC'].append(record.parsed_info_fields()['VAC'][i])
                        df['VariantSamples'].append(record.parsed_info_fields()['VariantSamples'][i])
                        for k in csq_order:
                            df[k].append(c[k])
                    else:
                        for f in ['CHROM', 'POS', 'REF']:
                            df[f].append(getattr(record, f))
                        df['ALT'].append(alt)
                        df['VAC'].append(record.parsed_info_fields()['VAC'][i])
                        df['VariantSamples'].append(record.parsed_info_fields()['VariantSamples'][i])
                        for k in csq_order:
                            df[k].append(c.get(k, '.'))
            i += 1
    df = pd.DataFrame.from_dict(df)
    unaligned = []
    with open(failures, 'rt') as infile:
        for line in infile:
            line = line.rstrip()
            u = os.path.splitext(os.path.basename(line))[0]
            unaligned.append(u)
    samples = [x for x in vcf.header.samples if x not in unaligned]
    df['AF'] = df.VAC/len(samples)
    df.to_csv(out, index=False)


if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.exit("Usage: {} annotated.vcf failures.txt out.csv".format(
            sys.argv[0]))
    main(*sys.argv[1:])
