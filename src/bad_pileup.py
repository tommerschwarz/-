# pseudo pile-up




def pileup(reads, outfile):
    snps = {}
    dels = {}
    ins = {}

    for line in reads:
        line = line.split('\t')
        pos_aligned = line[3]
        if (line[4] != '.'): ## SNP
            read_snps_pos = line[4].split(',')
            read_snps_ref = line[5].split(',')
            for i in range(len(read_snps_pos)):
                read_snp_p = int(read_snps_pos[i])
                read_snp_r = read_snps_ref[i]
                actual_pos = pos_aligned + read_snp_p
                if ((str(actual_pos),line[0][val]) in snps.keys()):
                    snps[(str(actual_pos),line[0][val])] = snps[(str(actual_pos),line[0][val])] + 1
                else:
                    snps[(str(actual_pos),line[0][val])] = 1
        if (line[6] != '.'): ## DEL
            read_del = line[6].split(',')
            for val in read_del:
                actual_pos = pos_aligned + val
                if ((str(actual_pos),line[0][val:val+line[6]]) in dels.keys()):
                    dels[(str(actual_pos),line[0][val:val+line[6]])] = dels[(str(actual_pos),line[0][val:val+line[6]])] + 1
                else:
                    dels[(str(actual_pos),line[0][val:val+line[6]])] = 1
    outfile.write(">SNP\n")
    for key in snps:
        if snps[key] > 5:
            outfile.write(key[1] + "," + key[0] + "\n")



reads = open('reads_hw2grad_M_1_chr_1.txt.aln', 'r')
outfile = open('out', 'w')
pileup(reads, outfile)
