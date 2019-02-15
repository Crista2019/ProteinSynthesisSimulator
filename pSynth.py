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
    RNANucleotideComplements = {'A':'U', 'T':'A', 'U':'A', 'C':'G', 'G':'C'}
    polypeptide_coded_from_sequence = []

    # Hashmap of the 20 Amino acids, their single-letter data-base codes (SLC), and their corresponding RNA codons
    I = ["AUU", "AUC", "AUA"]
    L = ["CUU", "CUC", "CUA", "CUG", "UUA", "UUG"]
    V = ["GUU", "GUC", "GUA", "GUG"]
    F = ["UUU", "UUC"]
    M = ["AUG"]
    C = ["UGU", "UGC"]
    A = ["GCU", "GCC", "GCA", "GCG"]
    G = ["GGU", "GGC", "GGA", "GGG"]
    P = ["CCU", "CCC", "CCA", "CCG"]
    T = ["ACU", "ACC", "ACA", "ACG"]
    S = ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"]
    Y = ["UAU", "UAC"]
    W = ["UGG"]
    Q = ["CAA", "CAG"]
    N = ["AAU", "AAC"]
    H = ["CAU", "CAC"]
    E = ["GAA", "GAG"]
    D = ["GAU", "GAC"]
    K = ["AAA", "AAG"]
    R = ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"]
    stop = ["UAA", "UAG", "UGA"]
    amino_acids = {"Isoleucine":I, "Leucine":L, "Valine":V, "Phenylalanine":F, "Methionine":M, "Cysteine":C, "Alanine":A, "Glycine":G, "Proline":P, "Threonine":T, "Serine":S, "Tyrosine":Y, "Tryptophan":W, "Glutamine":Q, "Asparagine":N, "Histidine":H, "Glutamic acid":E, "Aspartic acid":D, "Lysine":K, "Arginine":R, "Stop Codon":stop}

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
            self.mRNA_strand.append(self.RNANucleotideComplements[nucleotideBase])

    # STEP 3 of TRANSCRIPTION
    def transcription_termination(self):
        #Sequences called terminators signal that the RNA transcript is complete.
        #Once they are transcribed, they cause the transcript to be released from the RNA polymerase.
        if len(self.mRNA_strand) == len(self.DNA_template_strand):
            pass
        else:
            print("DNA => mRNA transcription not complete")

    def parsemRNAtotRNA(self, mRNA):
        tRNA = []
        for nucleotideBase in mRNA:
            tRNA.append(self.RNANucleotideComplements[nucleotideBase])
        return tRNA

    def transcribe(self):
        print("Original DNA:\n"+str(self.DNA_double_helix)+"\n")
        self.transcription_initiation()
        print("DNA_coding_strand:\n"+str(self.DNA_coding_strand)+"\n")
        print("DNA_template_strand:\n"+str(self.DNA_template_strand)+"\n")
        self.transcription_elongation()
        print("mRNA:\n"+str(self.mRNA_strand)+"\n")
        self.transcription_termination()
        print("tRNA:\n"+str(self.parsemRNAtotRNA(self.mRNA_strand))+"\n")
        return self.mRNA_strand

    def parseToCodons(self, mRNA):
        # converts mRNA into an array of 3 base codons
        codonIndex = 0
        eachCodon = ""
        for nucleotideBase in mRNA:
            if codonIndex < 3:
                eachCodon+=str(nucleotideBase)
                codonIndex+=1
            if codonIndex == 3:
                self.codons.append(eachCodon)
                codonIndex = 0
                eachCodon = ""
                # any partial nucleotide bases are not included (only full codons of 3 elements get parsed)
        print("String of codons:\n"+str(self.codons))

    def translation_initiation(self):
        # the ribosome gets together with the mRNA and the first tRNA so translation can begin.
        # An "initiator" tRNA carries the first amino acid in the protein, which is almost always methionine (Met)
        codonStartIndex = 0
        self.parseToCodons(self.mRNA_strand)
        for codon in self.codons:
            if codon == "AUG":
                codonStartIndex = self.codons.index(codon)
        print("(The start codon Methionine is codon #"+str(codonStartIndex+1)+" in the sequence)\n")
        return codonStartIndex

    def translation_elongation(self, codonStartIndex):
        # amino acids are brought to the ribosome by tRNAs and linked together to form a chain.
        self.codons = self.codons[codonStartIndex:] #splicing codon sequence to start at Met
        codonStopIndex = len(self.codons)
        for codon in self.codons:
            for amino_acid_name in self.amino_acids.keys():
                for each_amino_acid_codon in self.amino_acids[amino_acid_name]:
                    if codon == each_amino_acid_codon:
                        self.polypeptide_coded_from_sequence.append(amino_acid_name)
                        if codon in self.stop:
                            codonStopIndex = self.codons.index(codon)
        codonStopIndex +=1 # to allow the stop codon to be shown in the sequence
        return codonStopIndex

    def translation_termination(self, codonStopIndex):
        # the finished polypeptide is released to go and do its job in the cell.
        self.polypeptide_coded_from_sequence = self.polypeptide_coded_from_sequence[: codonStopIndex]
        print("polypeptide_coded_from_sequence:\n"+str(self.polypeptide_coded_from_sequence))
        if self.polypeptide_coded_from_sequence[len(self.polypeptide_coded_from_sequence)-1] == "Stop Codon":
            pass
        else:
            print("Error: sequence does not end at stop codon")

    def translate(self):
        codonStartIndex = self.translation_initiation()
        codonStopIndex = self.translation_elongation(codonStartIndex)
        self.translation_termination(codonStopIndex)

def main():
    dna_test_seq = [{'A':'T'}, {'C':'G'}, {'T':'A'}, {'A':'T'}, {'T':'A'}, {'G':'C'}, {'C':'G'}, {'C':'G'}, {'T':'A'}, {'C':'G'}, {'A':'T'}, {'T':'A'}, {'C':'G'}, {'G':'C'}, {'A':'T'}, {'A':'T'}, {'A':'T'}, {'G':'C'}, {'C':'G'}, {'G':'C'}, {'A':'T'}, {'T':'A'}, {'C':'G'}, {'G':'C'}, {'G':'C'}, {'G':'C'}, {'T':'A'}, {'G':'C'}, {'T':'A'}, {'A':'T'}, {'T':'A'}, {'A':'T'}, {'A':'T'}, {'G':'C'}, {'G':'C'}, {'G':'C'}]
    dna_test = DNA(dna_test_seq)
    dna_test.transcribe()
    dna_test.translate()

if __name__ == "__main__": main()
