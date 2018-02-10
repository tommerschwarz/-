"""Script to generate a pileup from aligned reads and a reference"""
import sys
import util    as u

def pileup(aligns, ref_genome):
    """take a file with aligned reads and extract variants."""
    print("Generating Pileup from file: ")
    print("   " + aligns)

    # open aligned reads file, and pileup file
    reads  = open(aligns, 'r')
    pileup = open(aligns + ".plp", "w")

    # set variables
    genome = ref_genome.getSequence()[0]
    pile   = [0] * len(genome)
    sym    = {"A":0,"C":1,"G":2,"T":3}
    let    = {0:"A",1:"C",2:"G",3:"T"}
    donor  = [list(pile),list(pile),list(pile),list(pile)]
    
    # write with the expected format
    pileup.write(">" + list(ref_genome.genome.keys())[0] +"\n")
    pileup.write(">SNP\n")

    # collect the information of ALL THE READS (modify this!)
    for read in reads:
        read = read.strip().split('\t')
        seq,init,pos = read[0],int(read[3]),read[2]
        if pos == "-":
            seq = u.rev(seq)                  # reverse if inverted
        for i in range(len(seq)):
            donor[ sym[seq[i]] ][init+i] += 1   # PILEUP!
    
    # call variants:
    depth = 0
    dom = 0
    for j in range(len(genome)):     # for each position in the reference...
        don = genome[j]
        for nt in range(4):
            depth += donor[nt][j]    # collect sequencing depth for position j
            if dom < donor[nt][j]:
                dom = donor[nt][j]   # save the max depth for any nucleotide in pos j
                don = let[nt]        #  and the identity of that nucleotide
        if don != genome[j]:         # if donor and reference diverge...
            if depth > 2:                                # QC: and donor has >2 reads
                if float(dom)/float(depth) > 0.6:        # QC: and freq(nt) >0.6
#                     print(genome[j],donor[0][j],donor[1][j],donor[2][j],donor[3][j],don,j)
                    pileup.write(genome[j] +","+ don +","+ str(j) +"\n")
        depth = 0
        dom = 0
    pileup.close()
    return 1
