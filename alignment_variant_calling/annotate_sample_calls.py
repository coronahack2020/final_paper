#!/usr/bin/env python
import sys
import os
import pysam


def main(f):
    vcf = pysam.VariantFile(f)
    out = os.path.splitext(f)[0]
    if out.endswith(".vcf"):
        out = os.path.splitext(out)[0]
    out += '.annot.vcf.gz'
    vcf.header.info.add(id="VAC",
                        number="A",
                        type="Integer",
                        description="Number of virus samples with ALT allele")
    vcf.header.info.add(id="VariantSamples",
                        number="A",
                        type="String",
                        description="ID of samples with ALT allele")
    vcf_out = pysam.VariantFile(out, 'w', header=vcf.header)
    for record in vcf:
        samps = list()
        for i in range(1, len(record.alleles)):
            samps.append([])
            for k, v in record.samples.items():
                if i in v.allele_indices:
                    samps[i-1].append(k)
        samp_counts = [len(samps[i]) for i in range(len(samps))]
        record.info['VAC'] = samp_counts
        record.info['VariantSamples'] = ["|".join(samps[i])
                                         for i in range(len(samps))]
        vcf_out.write(record)
    vcf.close()
    vcf_out.close()


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit("Usage: {} input.vcf".format(sys.argv[0]))
    main(sys.argv[1])
