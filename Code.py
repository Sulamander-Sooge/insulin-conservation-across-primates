import pandas as pd
from Bio import SeqIO, Phylo
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from collections import Counter

records = SeqIO.parse(
    'd:/Sumanyu/Relearning Python/ZZ-Mini Projects/2) Bioinformatics project on mammalian insulin sequences/BLAST Sequences/Z_combined_fasta_file.fasta', 'fasta')
for record in records:
    print(f'{record.id}: {len(record.seq)}')

alignment = AlignIO.read(
    'd:/Sumanyu/Relearning Python/ZZ-Mini Projects/2) Bioinformatics project on mammalian insulin sequences/BLAST Sequences/ZZ_clustalo MSA.fa', 'fasta')
dist_calc = DistanceCalculator('blosum62')
dist_matrix = dist_calc.get_distance(alignment)
constructor = DistanceTreeConstructor()
tree = constructor.nj(dist_matrix)
print(dist_matrix)
Phylo.draw_ascii(tree)
Phylo.write(tree, 'd:/Sumanyu/Relearning Python/ZZ-Mini Projects/2) Bioinformatics project on mammalian insulin sequences/tree.nwk', 'newick')

names = dist_matrix.names
n = len(names)

full_matrix = [[0.0] * n for _ in range(n)]

for i in range(n):
    for j in range(i + 1):
        full_matrix[i][j] = dist_matrix.matrix[i][j]
        full_matrix[j][i] = dist_matrix.matrix[i][j]

df_matrix = pd.DataFrame(full_matrix, index=names, columns=names)
df_matrix.to_csv(
    'd:/Sumanyu/Relearning Python/ZZ-Mini Projects/2) Bioinformatics project on mammalian insulin sequences/Distance Matrix.csv', index=False)


def analyse_conservation(alignment):
    results = []

    for position in range(alignment.get_alignment_length()):
        column = alignment[:, position]

        residues = [r for r in column if r != '-']

        counts = Counter(residues)
        residue, count = counts.most_common(1)[0]

        conservation = (count / len(residues)) * 100

        fully_conserved = []
        partially_conserved = []
        variable = []

        if conservation == 100:
            status = 'Fully Conserved'
            fully_conserved.append(position+1)
        elif conservation >= 75:
            status = 'Partially Conserved'
            partially_conserved.append(position+1)
        else:
            status = 'Variable'
            variable.append(position+1)

        results.append((position + 1, residue, conservation, status))

    return results, fully_conserved, partially_conserved, variable


results, fully_conserved, partially_conserved, variable = analyse_conservation(
    alignment)

for position, residue, conservation, status in results:
    print(f'Position: {position}')
    print(f"Consensus Residue: {residue}")
    print(f"Conservation: {conservation:.1f}%")
    print(f'Status: {status}')

samples = []

for position, residue, conservation, status in results:
    samples.append({
        'Position': position,
        'Consensus Residue': residue,
        'Conservation (%)': conservation,
        'Status': status
    })

df_samples = pd.DataFrame(samples)
# df_samples.to_csv('d:/Sumanyu/Relearning Python/ZZ-Mini Projects/2) Bioinformatics project on mammalian insulin sequences/Conservation Analysis.csv', index=False)
