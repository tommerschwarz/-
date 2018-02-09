import util as u

class genome():

    length = 0
    kmers = None
    genome = {}

    def __init__(self, ids, seqs):
        self.genome[ids] = seqs
        self.length = len(seqs)
        self.index_kmers()

    def add_sequence(ids, seqs):
        self.genome[ids] = seqs
        self.length += len(seqs)

    def index_kmers(self):
        """creates the index of the reference"""
        ref = self.genome[list(self.genome.keys())[0]]
        k = 6        # index positionos of kmers
        index = [None] * (4**k)
        for i in range(len(ref)-k+1):
            kmer_num = u.PatternToNumber( ref[i:i+k] )
            if(index[kmer_num] == None):
                index[kmer_num] = [i]
            else:
                index[kmer_num].append(i)
        self.kmers = index

    def getSequence(self):
    	chrs = list(self.genome.keys())
    	sequence = []
    	for chrom in chrs:
    		sequence.append(self.genome[chrom])
    	return sequence

    def index_BWT(self):
        sequence = self.getSequence()
        sequence = sequence + "$"
        index_m = []
        for i in range(len(sequence)''):
            index_m.append(sequence[i:] + sequence [:i])
            print(index_m[i])

