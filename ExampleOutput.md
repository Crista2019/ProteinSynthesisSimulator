## Example input in main:
```
dna_test_seq1 = [{'A':'T'}, {'C':'G'}, {'T':'A'}, {'A':'T'}, {'T':'A'}, {'G':'C'}, {'C':'G'}, {'C':'G'}, {'T':'A'}, {'C':'G'}, {'A':'T'}, {'T':'A'}, {'C':'G'}, {'G':'C'}, {'A':'T'}, {'A':'T'}, {'A':'T'}, {'G':'C'}, {'C':'G'}, {'G':'C'}, {'A':'T'}, {'T':'A'}, {'C':'G'}, {'G':'C'}, {'G':'C'}, {'G':'C'}, {'T':'A'}, {'G':'C'}, {'T':'A'}, {'A':'T'}, {'T':'A'}, {'A':'T'}, {'A':'T'}, {'G':'C'}, {'G':'C'}, {'G':'C'}]
dna_test1 = DNA(dna_test_seq1)
dna_test1.synthesize_polypeptides()
```

## Output printed to the console upon execution:
```
Original DNA:
[{'A': 'T'}, {'C': 'G'}, {'T': 'A'}, {'A': 'T'}, {'T': 'A'}, {'G': 'C'}, {'C': 'G'}, {'C': 'G'}, {'T': 'A'}, {'C': 'G'}, {'A': 'T'}, {'T': 'A'}, {'C': 'G'}, {'G': 'C'}, {'A': 'T'}, {'A': 'T'}, {'A': 'T'}, {'G': 'C'}, {'C': 'G'}, {'G': 'C'}, {'A': 'T'}, {'T': 'A'}, {'C': 'G'}, {'G': 'C'}, {'G':'C'}, {'G': 'C'}, {'T': 'A'}, {'G': 'C'}, {'T': 'A'}, {'A': 'T'}, {'T': 'A'}, {'A': 'T'}, {'A': 'T'}, {'G': 'C'}, {'G': 'C'}, {'G': 'C'}]

DNA_coding_strand:
['A', 'C', 'T', 'A', 'T', 'G', 'C', 'C', 'T', 'C', 'A', 'T', 'C', 'G', 'A', 'A', 'A', 'G', 'C', 'G', 'A', 'T', 'C', 'G', 'G', 'G', 'T', 'G', 'T', 'A', 'T', 'A', 'A', 'G', 'G', 'G']

DNA_template_strand:
['T', 'G', 'A', 'T', 'A', 'C', 'G', 'G', 'A', 'G', 'T', 'A', 'G', 'C', 'T', 'T', 'T', 'C', 'G', 'C', 'T', 'A', 'G', 'C', 'C', 'C', 'A', 'C', 'A', 'T', 'A', 'T', 'T', 'C', 'C', 'C']

mRNA:
['A', 'C', 'U', 'A', 'U', 'G', 'C', 'C', 'U', 'C', 'A', 'U', 'C', 'G', 'A', 'A', 'A', 'G', 'C', 'G', 'A', 'U', 'C', 'G', 'G', 'G', 'U', 'G', 'U', 'A', 'U', 'A', 'A', 'G', 'G', 'G']

tRNA:
['U', 'G', 'A', 'U', 'A', 'C', 'G', 'G', 'A', 'G', 'U', 'A', 'G', 'C', 'U', 'U', 'U', 'C', 'G', 'C', 'U', 'A', 'G', 'C', 'C', 'C', 'A', 'C', 'A', 'U', 'A', 'U', 'U', 'C', 'C', 'C']

String of codons:
['ACU', 'AUG', 'CCU', 'CAU', 'CGA', 'AAG', 'CGA', 'UCG', 'GGU', 'GUA', 'UAA', 'GGG']
(The start codon Methionine is codon #2 in the sequence)

polypeptide_coded_from_sequence:
['Methionine', 'Proline', 'Histidine', 'Arginine', 'Lysine', 'Arginine', 'Serine', 'Glycine', 'Valine', 'Stop Codon']
```
