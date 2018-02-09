"""Script to generate a pileup from aligned reads and a reference"""
import sys
import numpy   as np
import util    as u
import read    as rd
import genome  as gm
import variants_discoverer as vd

ref_fa  = sys.argv[1]
aligns  = sys.argv[2]
ref_seq = vd.load_reference(ref_fa)

def pileup(aligns, ref_genome):
    """take a file with aligned reads and extract variants."""
    print("Generating Pileup from file: ")
    print("   " + aligns)
    reads  = open(aligns, 'r')
    pileup = open(aligns + ".plp", "w")
    genome = ref_genome.getSequence()[0]
    pile = [0] * len(genome)
    sym = {"A":0,"C":1,"G":2,"T":3}
    let = {0:"A",1:"C",2:"G",3:"T"}
    donor = [list(pile),list(pile),list(pile),list(pile)]
    pileup.write(">" + list(ref_genome.genome.keys())[0] +"\n")
    pileup.write(">SNP\n")

    for read in reads:
        read = read.strip().split('\t')
        seq,init,pos = read[0],int(read[3]),read[2]
        if pos == "-":
            seq = u.rev(seq)
        for i in range(len(seq)):
            donor[sym[seq[i]]][init+i] += 1
    
    depth = 0
    dom = 0
    for j in range(len(genome)):
        don = genome[j]
        for nt in range(4):
            depth += donor[nt][j]
            if dom < donor[nt][j]:
                dom = donor[nt][j]
                don = let[nt]
        if don != genome[j]:
            if depth > 2:
                if float(dom)/float(depth) > 0.6:
#                     print(genome[j],donor[0][j],donor[1][j],donor[2][j],donor[3][j],don,j)
                    pileup.write(genome[j] +","+ don +","+ str(j) +"\n")
        depth = 0
        dom = 0
    pileup.close()
    return 1

pileup(aligns, ref_seq)
