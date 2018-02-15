class read():

    read_id  = None
    read_len = None
    sequence = None

    aligned   = False
    position  =  "."
    strand    =  "."
    mismatch  = ["."]
    deletion  = ["."]
    insertion = ["."]

    pair_position = "."
    pair_strand   = "."

    def __init__(self, rdid, seq):
        self.sequence = seq
        self.read_id = rdid
        self.read_len = len(seq)

    def summary(self):
        """returns a string with a summary of the alignment of the read"""
    # for mismatch
        if self.mismatch != ["."] and len(self.mismatch)>0:
            mism_p = []
            mism_a = []
            mismatches = self.mismatch.keys()
            for msm in mismatches:
                mism_p.append(           str(msm) )
                mism_a.append( self.mismatch[msm] )
        else:
            mism_p = ["."]
            mism_a = ["."]
    # for deletions
        if self.deletion != ["."] and len(self.deletion)>0:
            dels_p = []
            dels_a = []
            deletions = self.deletion.keys()
            for dtn in deletions:
                dels_p.append(           str(dtn) )
                dels_a.append( self.deletion[dtn] )
        else:
            dels_p = ["."]
            dels_a = ["."]
    # for insertions
        if self.insertion != ["."] and len(self.insertion)>0:
            ins_p = []
            ins_a = []
            insertions = self.insertion.keys()
            for itn in insertions:
                ins_p.append(            str(itn) )
                ins_a.append( self.insertion[itn] )
        else:
            ins_p = ["."]
            ins_a = ["."]

        if self.mismatch == []:
            self.mismatch = ["."]

        return '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
            self.sequence,self.read_id,
            self.strand,self.position,
            ",".join(mism_p),",".join(mism_a),
            ",".join(dels_p),",".join(dels_a),
            ",".join( ins_p),",".join( ins_a),
            self.pair_strand,self.pair_position)
