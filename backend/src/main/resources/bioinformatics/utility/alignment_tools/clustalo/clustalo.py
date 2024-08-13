# 1. $ cd clustal-omega-<version>
# 2. $ ./configure
# 3. $ make
# 4. $ sudo make install
# 5. $ make clean


# 0.   $ brew install argtable
# 2-1. $ ./configure CFLAGS='-I/opt/homebrew/include' LDFLAGS='-L/opt/homebrew/lib'


from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import SeqIO

# 입력 및 출력 파일 경로 설정
input_file = "input.fasta"
output_file = "output.fasta"

# Clustal Omega 명령어 생성
clustalomega_cline = ClustalOmegaCommandline(infile=input_file, outfile=output_file, verbose=True, auto=True)

# Clustal Omega 실행
stdout, stderr = clustalomega_cline()

# 결과 확인
print(stdout)
print(stderr)

# 정렬된 서열 출력
alignments = SeqIO.parse(output_file, "fasta")
for alignment in alignments:
    print(f">{alignment.id}\n{alignment.seq}")