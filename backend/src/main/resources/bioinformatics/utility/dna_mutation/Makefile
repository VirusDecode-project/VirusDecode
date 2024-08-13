# 변수 설정
PYTHON=python3
SCRIPT=level2.py
FASTA1=NC_045512.fasta
FASTA2=OX014251.fasta
MUSCLE=muscle

# 기본 목표
all: run

# Python 스크립트 실행
run: $(SCRIPT) $(FASTA1) $(FASTA2)
	$(PYTHON) $(SCRIPT)

# 정리 목표
clean:
	rm -f combined.fasta aligned.fasta mutations.html mutation.txt

.PHONY: all run clean
