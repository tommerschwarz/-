import numpy as np

class read():

    read_id  = None
    read_len = None
    sequence = None

    aligned  = False
    position = "."
    strand   = "."
    mismatch = ["."]
# TODO: insert del
# TODO: insert ins

    pair_position = "."
    pair_strand   = "."

    def __init__(self, rdid, seq):
        self.sequence = seq
        self.read_id = rdid
        self.read_len = len(seq)

    def summary(self):
        """returns a string with a summary of the alignment of the read"""
        return '{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(self.sequence,self.read_id,
        	self.strand,self.position,",".join(map(str, self.mismatch)),
        	self.pair_strand,self.pair_position)

