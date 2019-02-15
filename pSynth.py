from collections import defaultdict

class DNA:
    codons = []
    DNA_double_helix = [{'',''}] # array of one item dict for each base pairs to allow for duplicates e.g. {{'A':'T'}, {'C':'G'}, {'A':'T'}}
    DNA_coding_strand = [] # 5' to 3'
    DNA_template_strand = [] # 3' to 5'
    mRNA_strand = []
    tRNA_strand = []
    length = 0
    DNANucleotideBasePairings = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    RNANucleotideComplementsTODNA = {'A':'U', 'T':'A', 'C':'G', 'G':'C'}

    # Hashmap of the 20 Amino acids, their single-letter data-base codes (SLC), and their corresponding DNA codons
    amino acids = {}

    def __init__(self, sequence):
        self.DNA_double_helix = sequence;

    # STEP 1 of TRANSCRIPTION
    def transcription_initiation(self):
        # RNA polymerase binds to a sequence of DNA called the promoter, found near the beginning of a gene.
        # Each gene (or group of co-transcribed genes, in bacteria) has its own promoter.
        # Once bound, RNA polymerase separates the DNA strands, providing the single-stranded template needed for transcription.
        for basePairs in self.DNA_double_helix: # for each {u,v} in [{w,x}, {y, z}]
            for aminoAcid in basePairs.keys(): # for each u in {u,v}
                self.DNA_coding_strand.append(aminoAcid)
            for aminoAcid in basePairs.values(): # for each v in {u,v}
                self.DNA_template_strand.append(aminoAcid)

    # STEP 2 of TRANSCRIPTION
    def transcription_elongation(self):
        # One strand of DNA, the template strand, acts as a template for RNA polymerase.
        # As it "reads" this template one base at a time, the polymerase builds an RNA molecule out of complementary nucleotides, making a chain that grows from 5' to 3'.
        # The RNA transcript carries the same information as the non-template (coding) strand of DNA, but it contains the base uracil (U) instead of thymine (T).
        for nucleotideBase in self.DNA_template_strand:
            self.mRNA_strand.append(self.RNANucleotideComplementsTODNA[nucleotideBase])

    # STEP 3 of TRANSCRIPTION
    def transcription_termination(self):
        #Sequences called terminators signal that the RNA transcript is complete.
        #Once they are transcribed, they cause the transcript to be released from the RNA polymerase.
        if len(self.mRNA_strand) == len(self.DNA_template_strand):
            pass
        else:
            print("DNA => mRNA transcription not complete")

    def transcribe(self):
        print("Original DNA:\n"+str(self.DNA_double_helix)+"\n")
        self.transcription_initiation()
        print("DNA_coding_strand:   "+str(self.DNA_coding_strand))
        print("DNA_template_strand: "+str(self.DNA_template_strand))
        self.transcription_elongation()
        print("mRNA:                "+str(self.mRNA_strand))
        self.transcription_termination()
        return self.mRNA_strand

    def translation_initiation():
        # the ribosome gets together with the mRNA and the first tRNA so translation can begin.
        pass

    def translation_elongation():
        # amino acids are brought to the ribosome by tRNAs and linked together to form a chain.
        pass

    def translation_termination():
        # the finished polypeptide is released to go and do its job in the cell.
        pass

    def translate(self):
        pass

def main():
    dna_test_seq = [{'A':'T'}, {'C':'G'}, {'T':'A'}, {'G':'C'}, {'T':'A'}, {'C':'G'}, {'A':'T'}]
    dna_test = DNA(dna_test_seq)
    dna_test.transcribe()

if __name__ == "__main__": main()
