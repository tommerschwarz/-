import numpy as np
import read
import genome as gm
import util as u
import copy

def align_trivial(reads, genome): # takes a list with two reads, and a ref genome
    reference = genome.getSequence()[0] # for now takes only first chromosome
    read1_fw = reads[0].sequence
    read1_rv = u.rev(read1_fw)
    read2_fw = reads[1].sequence
    read2_rv = u.rev(read2_fw)
    all_alng = {0:[],1:[]}           # stores all valid alignments for each end of the read
    best_aln = reads                 # contains alignments to report
    d = 0                            # max number of mismatches allowed when aligning

   # iterate only once over regions of the genome of the size of the reads
    for i in range( len(reference) - len(read1_fw) + 1 ):
        ref = reference[i:(i+len(read1_fw))]

      # check end_1
        alignment = copy.deepcopy(reads[0])
      # forward
        diff = HammingDistance( read1_fw, ref, d)
        if diff != None:
            alignment.aligned = True
            alignment.position = i
            alignment.strand = "+"
            alignment.mismatch = diff if len(diff) > 0 else ["."]
            all_alng[0].append(alignment)
      # reverse
        else:
            diff = HammingDistance( read1_rv, ref, d)
            if diff != None:
                alignment.aligned = True
                alignment.position = i
                alignment.strand = "-"
                alignment.mismatch = diff if len(diff) > 0 else ["."]
                all_alng[0].append(alignment)

      # check end_2
        alignment = copy.deepcopy(reads[1])
      # forward
        diff = HammingDistance( read2_fw, ref, d)
        if diff != None:
            alignment.aligned = True
            alignment.position = i
            alignment.strand = "+"
            alignment.mismatch = diff if len(diff) > 0 else ["."]
            all_alng[1].append(alignment)
      # reverse
        else:
            diff = HammingDistance( read2_rv, ref, d)
            if diff != None:
                alignment.aligned = True
                alignment.position = i
                alignment.strand = "-"
                alignment.mismatch = diff if len(diff) > 0 else ["."]
                all_alng[1].append(alignment)

    # check if there were any valid alignments for end_1
    if len(all_alng[0]) > 0:
        # from all valid alignments, get the best for end_1
        best_aln[0] = all_alng[0][0]
        for j in range(len(all_alng[0])):
            if len(all_alng[0][j].mismatch) <= len(best_aln[0].mismatch):
                best_aln[0] = all_alng[0][j]
    # check if there were any valid alignments for end_2
    if len(all_alng[1]) > 0:
        # from all valid alignments, get the best for end_2
        best_aln[1] = all_alng[1][0]
        for j in range(len(all_alng[1])):
            if len(all_alng[1][j].mismatch) <= len(best_aln[1].mismatch):
                best_aln[1] = all_alng[1][j]

    # add to best alignments the information from the other end
    best_aln[0].pair_position = best_aln[1].position
    best_aln[0].pair_strand   = best_aln[1].strand

    best_aln[1].pair_position = best_aln[0].position
    best_aln[1].pair_strand   = best_aln[0].strand

    return best_aln

def HammingDistance( read, reference , max_mis):
    positions = []
    for i in range(len(read)):
        if len(positions) > max_mis:
            return None
        if read[i] != reference[i]:
            positions.append(i)
    return positions

def align_seeded(read, genome):
    """aligns a read to a reference genome
    using a list of lists with the starting
    position of each kmer length 8"""
    pass

def align_bwt(reads, index, count):
    """aligns a read to a reference genome
    following the BWT index"""

    read1_fw = reads[0].sequence
    read1_rv = u.rev(read1_fw)
    read2_fw = reads[1].sequence
    read2_rv = u.rev(read2_fw)
    all_alng = {0:[],1:[]}           # stores all valid alignments for each end of the read
    best_aln = reads                 # contains alignments to report

  # check end_1
    alignment = copy.deepcopy(reads[0])
  # forward
    aln1 = perfect_match_bwt(read1_fw, index, count)
    if aln1 != ".":
        alignment.aligned = True # perfect alignment found
        alignment.position = aln1
        alignment.strand = "+"
        all_alng[0].append(alignment)
  # reverse
    else:
        aln2 = perfect_match_bwt(read1_rv, index, count)
        if aln2 != ".":
            alignment.aligned = True
            alignment.position = aln2
            alignment.strand = "-"
            all_alng[0].append(alignment)

  # check end_2
    alignment = copy.deepcopy(reads[1])
  # forward
    aln3 = perfect_match_bwt(read2_fw, index, count)
    if aln3 != ".":
        alignment.aligned = True
        alignment.position = aln3
        alignment.strand = "+"
        all_alng[1].append(alignment)
  # reverse
    else:
        aln4 = perfect_match_bwt(read2_rv, index, count)
        if aln4 != ".":
            alignment.aligned = True
            alignment.position = aln4
            alignment.strand = "-"
            all_alng[1].append(alignment)

    # save first alignment if there was any valid for both ends
    perf_val = 0 # 0 if neither is perfect, 1 if first is perfect, 2, if second is perfect, 3 if both are perfect
    if len(all_alng[0]) > 0: # first read perfect
        perf_val = perf_val + 1
        best_aln[0] = all_alng[0][0]
    if len(all_alng[1]) > 0:
        perf_val = perf_val + 2
        best_aln[1] = all_alng[1][0]
    

    if (perf_val == 1): # perfect alignment in the first read
        local_align(ref, read2)
        return
    elif (perf_val == 2): # perfect alignment in the second read
        local_align(ref, read1)
        return

        

# TODO insert local alignment here... use best_aln[0 or 1]

    local_align(ref, read)

    # add to best alignments the information from the other end
    best_aln[0].pair_position = best_aln[1].position
    best_aln[0].pair_strand   = best_aln[1].strand

    best_aln[1].pair_position = best_aln[0].position
    best_aln[1].pair_strand   = best_aln[0].strand

    return best_aln

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

def perfect_match_bwt(read, index, count):
    # initialize pointers up and down 
    #      with the last position of the read
    i = len(read) - 1
    last = read[i]
    nt = ["A","C","G","T"]
    cnt = 0
    for n in nt:
        if n != last:
            cnt += count[n]      # advance in the index to the pos of last nt
        else:
            up = cnt + 1         # pointer to the first (in 1st col)
            dn = cnt + count[n]  # pointer to the last  (in 1st col)
            break

    # search over the whole read from back to front:
    while i > 0 and up <= dn:
        i -= 1

        # save the next to last nt
        last = read[i]

        # search within range until match
        while index[up][0] != last:
            up += 1
        while index[dn][0] != last:
            dn -= 1

        # move pointers to new position in index
        cnt = 0
        for n in nt:
            if n != last:
                cnt += count[n]
            elif i > 0:                 # dont move if already end # might be causing problems...
                up = cnt + index[up][1]
                dn = cnt + index[dn][1]
                break

    if dn-up == 0 and i == 0:
        return index[up][2]        # starting position on genome
    else:
        return "."                 # no perfect alignment found
