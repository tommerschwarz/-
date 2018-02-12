import sys

def opt(s1, s2):
    a = len(s1) + 1 #rows
    b = len(s2) + 1 #cols
    max_val = 0
    x = [[0]*(b) for i in range(a)]
    bt = [[0]*(b) for i in range(a)]
    match = 1
    mm = -1
    gap = -2
    ins = -2
    dels = -1

    for i in range(1, a):
        for j in range(1, b):
            if (s1[i-1] == s2[j-1]):
                t = match
            else:
                t = mm
            arr = [x[i-1][j] + ins, x[i-1][j-1] + t, x[i][j-1] + dels]
            x[i][j] = max(arr)
            if (arr.index(x[i][j]) == 1):
                if (t == match):
                    bt[i][j] = 1 #match
                else:
                    bt[i][j] = 3 #mm
            else:
                bt[i][j] = arr.index(x[i][j])
            if ((x[i][j] > max_val) & (j == b-1) & (i > b)):
                max_val = x[i][j]
                max_i = i
                max_j = j
    
    i = max_i
    j = max_j
    if (x[i][j] < 43):
        return 0

    seq1 = ''
    seq2 = ''
    ins = {}
    dels = {}
    snps = {}
    num_match = 0
    num_ins = 0
    num_del = 0
    num_mm = 0


    while ((j > 0) & (i > 0)):
        if (bt[i][j] == 1): # match
            #bt[i][j] = '*'
            seq1 = s1[i-1] + seq1
            seq2 = s2[j-1] + seq2
            i = i - 1
            j = j - 1
            num_match = num_match + 1
            #print(s1[i]+' '+s2[j])
            #print(j)
        elif (bt[i][j] == 3): # SNP
            #bt[i][j] = '*'
            snps[i-1] = s1[j]
            seq1 = s1[i-1].lower() + seq1
            seq2 = s2[j-1].lower() + seq2
            i = i - 1
            j = j - 1
            num_mm = num_mm + 1
            #print(s1[i]+' != '+s2[j])
        elif (bt[i][j] == 2): # insertion
            #bt[i][j] = '*'
            ins[i-1] = s2[j-1]
            seq1 = '-' + seq1
            seq2 = s2[j-1] + seq2
            j = j - 1
            num_ins = num_ins + 1
        elif (bt[i][j] == 0): # deletion
            #bt[i][j] = '*'
            dels[i-1] = s2[j-1]
            seq1 = s1[i-1] + seq1
            seq2 = '-' + seq2
            i = i - 1
            num_del = num_del + 1
    bt[max_i][max_j] = 88
    
    return snps, ins, dels

#s1 = 'GGGAAAGGCAGACAGCGAC CTACACGTAAACGCCCTCAGAACAGGGGTTCCTATGTGGTTAAGCTAC GGCAAATCTTTCAATGGGCAAAAGCTTACACAGTTGCCTGTATAAATCCCTAATCGAGCACAAATCAGAACCTTAGGGGTATGTAGCTACCAAGCAGCGTTCCACACCTTATCAGCGCCGTGGGTTGGCGCGGCCCGTCTACCTCTTCCCCAGAAAGAAGCGGACCAAAAGTCAAAACTTTGCCATGTCTCCTCAGTCGTTACAAGGAGGAAATACACGGGAGCACGGGGTTGATACTAGCGCGGACAAGGTGGAGCATGTCGTAACCAAACAACCTGCTCGCATCCTGTGGCCTGTCTTGTACGGTGCTGTTATGATCTCAATTTCCTATCTGGAGGGCCCCACTTCATACGCGAAACCTAGGTGCCATAATTGCCGACCTAGATACGGTGGGGCAGTTTATGTGATATGGCAAACATGGTCTGTAGACGTCGATGCCTCTCCAGCTGCTTTTAGACGCCCTGTCAATCAATGTCCTCTCTGACCATTGACAGCTCCCCATACGCCGCTACCGTACATACCGCCCATTTGTAACTTAACAACTCTTTACTCGCCCGTCTACCGCGTTCTTTGCTACAAGATGCTCACATCTTTGCCGCAGAAGAGATAGACTGTCCGCCTACCAGATATAATTGATTTCAGTGCTTCGCTTCTGAGCTTTCAAAAAGCATACAAATAGCTTTTCTATAACCACTAGCGAGCCTCACCAATTAACCTGGTGTGGTCACATCGTTTCGTAGACGAGAATTCGGTTGGACTGGCAGATGAATCGGTGGGGGAAGAGAGAGTTCCCGCGACTGACGTTGAATACGGTGTGTCTCGAAATGTTTCGTTAACACTAGACTTGACGGGAGCTGACCCTAATTGGATGGAGTCGACTTCTGTC'
s1 = 'ATGGTAGTCGGCTACGTAACGCTACGTAACTGATCGTAGCTAGCTGACGATTGTGATGCTAGCTGATGCTAATCGTACGATCGTAGCTAG' #ref
s2 = 'GCTACGTAACTGATCGTAGCTAGCTGACGATTGTGATGCTAGCTGATGCT'  #read

print(opt(s1, s2))

