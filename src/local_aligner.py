import sys

def local_align(s1, s2):
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
    del_open = False
    ins_open = False

    for i in range(1, a):
        for j in range(1, b):
            if (s1[i-1] == s2[j-1]):
                t = match
            else:
                t = mm
            if (bt[i-1][j] == 0):
                ins = -1
            else:
                ins = -2
            if (bt[i][j-1] == 2):
                dels = -1
            else:
                dels = -2
            arr = [x[i-1][j] + ins, x[i-1][j-1] + t, x[i][j-1] + dels]
            x[i][j] = max(arr)
            if (arr.index(x[i][j]) == 1):
                if (t == match):
                    bt[i][j] = 1 #match
                else:
                    bt[i][j] = 3 #mm
                del_open = False
                ins_open = False
                bt[i][j] = arr.index(x[i][j])
            '''else:
                if (bt[i][j] == 2):
                    ins_open = True
                else:
                    del_open = True'''
                
            if ((x[i][j] > max_val) & (j == b-1) & (i > b)):
                max_val = x[i][j]
                max_i = i
                max_j = j
    
    i = max_i
    j = max_j
    if (x[i][j] < 30):
        return 0,0,0

    seq1 = ''
    seq2 = ''
    ins = {}
    dels = {}
    snps = []
    num_match = 0
    num_ins = 0
    num_del = 0
    num_mm = 0
    p_ins = 0 # how many prior positions were insertions
    p_del = 0 # how many prior positions were deletions


    while ((j > 0) & (i > 0)):
        if (bt[i][j] == 1): # match
            seq1 = s1[i-1] + seq1
            seq2 = s2[j-1] + seq2
            i = i - 1
            j = j - 1
            num_match = num_match + 1
            if (p_ins > 0):
                ins[j] = p_ins
                p_ins = 0
            if (p_del > 0):
                dels[j] = p_del
                p_del = 0
        elif (bt[i][j] == 3): # SNP
            snps.append(j-1)
            seq1 = s1[i-1].lower() + seq1
            seq2 = s2[j-1].lower() + seq2
            i = i - 1
            j = j - 1
            num_mm = num_mm + 1
            if (p_ins > 0):
                ins[j] = p_ins
                p_ins = 0
            if (p_del > 0):
                ins[j] = p_del
                p_del = 0
        elif (bt[i][j] == 2): # insertion
            p_ins = p_ins + 1
            seq1 = '-' + seq1
            seq2 = s2[j-1] + seq2
            j = j - 1
            num_ins = num_ins + 1
        elif (bt[i][j] == 0): # deletion
            p_del = p_del + 1
            seq1 = s1[i-1] + seq1
            seq2 = '-' + seq2
            i = i - 1
            num_del = num_del + 1
    
    return snps, ins, dels, i
     #          1         2         3         4         5         6         7         8         9
#s1 ='0123456789012345678 01234567890123456789012345678901234567890123456789 12345678901234567890'
s1 = 'ATGGTAGTCGCTACGCCTACGCTACGTAACTGATCGTAGCTAGCTGACGATTTTGGTGATGCTAGCTGATAGTGATGGAGTCGACTTCTGTC' #ref
s2 =                     'GCTACGTAACTGATCGTAGCTAGCTGACGATTGTGATGCTAGCTGATAGT'  #read

print(local_align(s1, s2))

