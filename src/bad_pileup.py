#! /usr/bin/env python
import sys


def pileup(reads, outfile):
    snps = {}
    dels = {}
    ins = {}

    for line in reads:
        line = line.split('\t')
        if (line[3] == "."):
            continue
        pos_aligned = int(line[3]) + 1
        if (line[4] != '.'): ## SNP
            read_snps_pos = line[4].split(',')
            read_snps_ref = line[5].split(',')
            for i in range(len(read_snps_pos)):
                read_snp_p = int(read_snps_pos[i])
                read_snp_r = read_snps_ref[i]
                actual_pos = pos_aligned + read_snp_p - 1
                if ((str(actual_pos), read_snp_r, line[0][read_snp_p]) in snps.keys()):
                    # print (str(actual_pos),read_snp_r, line[0][read_snp_p],line[2])
                    snps[(str(actual_pos), read_snp_r, line[0][read_snp_p])] = snps[(str(actual_pos), read_snp_r, line[0][read_snp_p])] + 1
                else:
                    snps[(str(actual_pos), read_snp_r, line[0][read_snp_p])] = 1
        if (line[6] != '.'): ## DEL
            read_dels_pos = line[6].split(',')
            read_dels_ref = line[7].split(',')
            for j in range(len(read_dels_pos)):
                read_del_p = int(read_dels_pos[j])
                read_del_r = read_dels_ref[j]
                actual_pos = pos_aligned + read_del_p
                if ((str(actual_pos),read_del_r) in dels.keys()):
                    # print (str(actual_pos),read_del_r,line[2])
                    dels[(str(actual_pos),read_del_r)] = dels[(str(actual_pos),read_del_r)] + 1
                else:
                    dels[(str(actual_pos),read_del_r)] = 1
        if (line[8] != '.'): ## INS
            read_inss_pos = line[8].split(',')
            read_inss_ref = line[9].split(',')
            for k in range(len(read_inss_pos)):
                read_ins_p = int(read_inss_pos[k])
                read_ins_r = read_inss_ref[k]
                actual_pos = pos_aligned + read_ins_p
                if ((str(actual_pos),read_ins_r) in ins.keys()):
                    # print (str(actual_pos),read_ins_r,line[2])
                    ins[(str(actual_pos),read_ins_r)] = ins[(str(actual_pos),read_ins_r)] + 1
                else:
                    ins[(str(actual_pos),read_ins_r)] = 1
    outfile.write(">INS\n")
    for key in ins:
        if ins[key] > 1: # change threshold for QC-ing errors
            outfile.write(key[1] + ',' + key[0] + '\n')
    outfile.write(">DEL\n")
    for key in dels:
        if dels[key] > 1: # change threshold for QC-ing errors
            outfile.write(key[1] + ',' + key[0] + '\n')
    outfile.write(">SNP\n")
    for key in snps:
        if snps[key] > 1: # change threshold for QC-ing errors
            outfile.write(key[1] + ',' + key[2] + ',' + key[0] + '\n')



reads = open(sys.argv[1], 'r')
outfile = open(sys.argv[2], 'w')
pileup(reads, outfile)





