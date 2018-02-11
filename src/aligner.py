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
    d = 1                            # max number of mismatches allowed when aligning

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

    nt = ["A","C","G","T"]

  # check end_1
    alignment = copy.deepcopy(reads[0])
  # forward

    # initialize up and down in the last position of the read
    i = len(read1_fw) - 1
    last = read1_fw[i]
    cnt = 0
    for n in nt:
        if n != last:
            cnt += count[n]      # advance in the index to the pos of last nt
        else:
            up = cnt + 1         # pointer to the first (in 1st col)
            dn = cnt + count[n]  # pointer to the last  (in 1st col)
            break
    print(up,index[up],index[dn],dn,last)

    # search over the whole read from back to front:
    while i > 0: #and up < dn:
        i -= 1

        # save the next to last nt
        last = read1_fw[i]

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
            else:
                up = cnt + index[up][1]
                dn = cnt + index[dn][1]
                break
        print(up,index[up],index[dn],dn,last)

    print(index[up],up, index[dn],dn,last,read1_fw,count)


""" uncomment here

    seq = index[0][0]
    pos = index[0][1]
    while seq[0] != "$":         # until the first character
        cpos = 0
        for n in nt:
            if n != seq[0]:
                cpos += count[n] # advance first column index per nt
            else:
                break
        cpos += pos              # advance to the position of the nt
        seq = index[cpos][0] + seq
        if index[cpos][0] != "$":
            pos = index[cpos][1]
    return seq[1:]               # omit initial '$'
"""

"""
        alignment.aligned = True
        alignment.position = i
        alignment.strand = "+"
        alignment.mismatch = diff if len(diff) > 0 else ["."]
        all_alng[0].append(alignment)

   leave for later...
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
    up = 0   # pointer to the first (in 1st col)
    dn = 0   # pointer to the last  (in 1st col)
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
"""

""" uncomment here

    # check if there were any valid alignments for end_1
    if len(all_alng[0]) > 0:
        # from all valid alignments, get the best for end_1
        best_aln[0] = all_alng[0][0]
        for j in range(len(all_alng[0])):
            if len(all_alng[0][j].mismatch) <= len(best_aln[0].mismatch):
                best_aln[0] = all_alng[0][j]
"""

"""
# report the best alignments for end 2 and update the pair end information

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
"""

#    return best_aln


























