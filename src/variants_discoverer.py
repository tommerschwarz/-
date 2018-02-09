"""Script to find variants from a donor in NGS paired end reads."""
import sys
import numpy   as np
import util    as u
import read    as rd
import genome  as gm
import aligner as aln

def main():
    """Main entry point for the script."""
    """
    for inp in range(1,len(sys.argv)):
        print(sys.argv[inp])
        if sys.argv[inp].startswith("-"):
            pass
        else:
        	break
            alns = generate_pileup()
    """
    ref_fa = sys.argv[1]
    reads  = sys.argv[2]
    ref_seq = load_reference(ref_fa)
#    align_reads(reads, ref_seq)
#    ref_seq.index_kmers()

# TODO: add option for unpaired reads
def align_reads(read_fn, ref_genome):
    """take a file with reads and align read by read."""
    print("Aligning reads from file: ")
    print("   " + read_fn)
    fastq   = open(read_fn, 'r')
    mapped  = open(read_fn + ".aln", "w")
    read_id = 0
    for read in fastq:
        if read_id > 0:                    # Skip the first line
            read = read.strip()    
            reads_raw = read.split(',')    # pairs are separated by comma
            read_1 = rd.read(str(read_id)+"_1", reads_raw[0])
            read_2 = rd.read(str(read_id)+"_2", reads_raw[1])
            aln_reads = aln.align_trivial([read_1,read_2], ref_genome)
            for alignment in aln_reads:
                mapped.write( alignment.summary() +"\n")
        if read_id % 10 == 0:
            print("{} reads aligned".format(read_id))
        read_id += 1
    mapped.close()
    return 1

def load_reference(ref_fn):
    fasta = open(ref_fn, 'r')
    seqfa = []
    seqid = []
    seq = ""
    fline = False
    for line in fasta:
        if line.startswith(">"):
            seqid.append(line.strip()[1:])
            if fline:
                seqfa.append(seq)
                seq = ""
            fline = True
        else:
            seq += line.strip()
    seqfa.append(seq)
    for i in range(len(seqid)):
        genome = gm.genome.add_sequence(seqid[i],seqfa[i])
    print("Reference genome loaded.")
    print("         total length: {}".format(genome.length))
    return genome

def print_variants():
    """format and print the list of variatns"""
    pass

if __name__ == '__main__':
    sys.exit(main())
