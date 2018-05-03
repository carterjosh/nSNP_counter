bases = ['T', 'C', 'A', 'G']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))

#pncA sequence excluding start and stop codons
seq = 'ATGCGGGCGTTGATCATCGTCGACGTGCAGAACGACTTCTGCGAGGGTGGCTCGCTGGCGGTAACCGGTGGCGCCGCGCTGGCCCGCGCCATCAGCGACTACCTGGCCGAAGCGGCGGACTACCATCACGTCGTGGCAACCAAGGACTTCCACATCGACCCGGGTGACCACTTCTCCGGCACACCGGACTATTCCTCGTCGTGGCCACCGCATTGCGTCAGCGGTACTCCCGGCGCGGACTTCCATCCCAGTCTGGACACGTCGGCAATCGAGGCGGTGTTCTACAAGGGTGCCTACACCGGAGCGTACAGCGGCTTCGAAGGAGTCGACGAGAACGGCACGCCACTGCTGAATTGGCTGCGGCAACGCGGCGTCGATGAGGTCGATGTGGTCGGTATTGCCACCGATCATTGTGTGCGCCAGACGGCCGAGGACGCGGTACGCAATGGCTTGGCCACCAGGGTGCTGGTGGACCTGACAGCGGGTGTGTCGGCCGATACCACCGTCGCCGCGCTGGAGGAGATGCGCACCGCCAGCGTCGAGTTGGTTTGCAGCTCC'

seq_codons = [[seq[x*3:x*3+3]] for x in range(0,int(len(seq)/3))]
stop_list = [0 for i in range(0,len(seq_codons))]
codon_list = [0 for i in range(0,len(seq_codons))]
unique_list = [0 for i in range(0,len(seq_codons))]

print('The input sequence is', int(len(seq)/3), 'AAs long.')

codon_count=0
for codon in seq_codons:
    AAs = []
    AAs.append(codon_table[codon[0]])
    for pos in range(0,len(codon[0])):
        SNPs = set('ACGT') - set(codon[0][pos])
        for snp in SNPs:
            test_codon = list(codon[0])
            test_codon[pos] = snp
            test_codon = ''.join(test_codon)
            if codon_table[test_codon] == '*':
                stop_list[codon_count] = stop_list[codon_count] + 1
            elif codon_table[test_codon] != codon_table[codon[0]]:
                codon_list[codon_count] = codon_list[codon_count] + 1
            if codon_table[test_codon] not in AAs and codon_table[test_codon] != '*':
                unique_list[codon_count] = unique_list[codon_count] + 1
                AAs.append(codon_table[test_codon])
    codon_count+=1
print('The total number of possible stop SNPs is', sum(stop_list))
print('The total number of possible nSNPs is', sum(codon_list))
print('The total number of possible unique nsAAs is', sum(unique_list))
print(unique_list)