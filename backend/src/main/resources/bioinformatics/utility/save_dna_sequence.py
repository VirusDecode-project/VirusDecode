from Bio import Entrez, SeqIO

# NCBI에 이메일 주소 등록
Entrez.email = "example@gmail.com"

class NCBISequenceFetcher:
    def fetch_nucleotide_sequence(self, nucleotide_id):
        handle = Entrez.efetch(db="nucleotide", id=nucleotide_id, rettype="gb", retmode="text")
        seq_record = SeqIO.read(handle, "genbank")
        handle.close()
        return seq_record

    def fetch_fasta_sequence(self, nucleotide_id, file_path):
        handle = Entrez.efetch(db="nucleotide", id=nucleotide_id, rettype="fasta", retmode="text")
        with open(file_path, 'w') as fasta_file:
            fasta_file.write(handle.read())
        handle.close()
        print(f"FASTA sequence saved to {file_path}")


# 사용 예시
nucleotide_id_list = ["OL672836.1", "MW642250.1", "OM958567.1"]
# nucleotide_id = "PP346415.1"  # 예: coronavirus 2

# 서열 데이터 가져오기
fetcher = NCBISequenceFetcher()

for nucleotide_id in nucleotide_id_list:
    sequence_record = fetcher.fetch_nucleotide_sequence(nucleotide_id)

    # FASTA 파일 저장
    fetcher.fetch_fasta_sequence(nucleotide_id, f"{nucleotide_id}.fasta")
