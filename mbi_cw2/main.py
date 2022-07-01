from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable

if __name__ == '__main__':
    table = CodonTable.unambiguous_dna_by_name["Standard"]
    records = list(SeqIO.parse("hymenolepis_diminuta.PRJEB507.WBPS10.protein.fa", "fasta"))

    protein = Seq(records[1].seq)
    print(f"Białko: \n {protein}")
    dna = Seq("")
    for i in range(len(protein)):
        dna += table.back_table[f"{protein[i]}"]
    print(f"DNA: \n {dna}")
