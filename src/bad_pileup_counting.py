    outfile.write(">INS\n")
    for key in ins:
        max_allele_count == 1
        final_ins = ''
        for allele in ins[key]:
            if (ins[key][allele] > max_allele_count): # change threshold for QC-ing errors
                max_allele_count = ins[key][allele]
                final_ins = allele
        if (final_ins != ''):
            outfile.write(key + ',' + final_ins + '\n')
    outfile.write(">DEL\n")
    for key in dels:
        max_allele_count == 1
        final_del = ''
        for allele in dels[key]:
            if (dels[key][allele] > max_allele_count): # change threshold for QC-ing errors
                max_allele_count = dels[key][allele]
                final_del = allele
        if (final_del != ''):
            outfile.write(key + ',' + final_del + '\n')
    outfile.write(">SNP\n")
    for key in snps:
        max_snp_count = 1
        final_snp = ''
        for snp in snps[key]:
            if (snps[key][snp] > max_snp_count):
                max_snp_count = snps[key][snp]
                final_snp = snp
        if (final_snp != ''):
            outfile.write(key + ',' + ref_snp[key] + ',' + final_snp + '\n')

