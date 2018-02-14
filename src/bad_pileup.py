# pseudo pile-up




def pileup(reads, outfile):
    snps = {}
    ins = {}
    dels = {}

    for line in reads:
        line = line.split('\t')
        pos_aligned = line[3]
        if (line[4] != '.'):
            read_snps = line[4].split(',')
            for val in read_snps:
                actual_pos = pos_aligned + val
                if ((str(actual_pos),line[0][val]) in snps.keys()):
                    snps[(str(actual_pos),line[0][val])] = snps[(str(actual_pos),line[0][val])] + 1
                else:
                    snps[(str(actual_pos),line[0][val])] = 1
        if (line[5] != '.'):
            read_ins = line[5].split(',')
            for val in read_ins:
                actual_pos = pos_aligned + val
                if ((str(actual_pos),line[0][val:val+line[6]]) in ins.keys()):
                    ins[(str(actual_pos),line[0][val:val+line[6]])] = ins[(str(actual_pos),line[0][val:val+line[6]])] + 1
                else:
                    ins[(str(actual_pos),line[0][val:val+line[6]])] = 1
    outfile.write(">SNP\n")
    for key in snps:
        if snps[key] > 5:
            outfile.write(key[1] + "," + key[0] + "\n")



reads = open('reads_hw2grad_M_1_chr_1.txt.aln', 'r')
outfile = open('out', 'w')
pileup(reads, outfile)
