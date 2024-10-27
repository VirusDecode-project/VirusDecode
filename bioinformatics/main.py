import sys
from analysis import Metadata, Alignment, MRNADesign, ThreeDViewer

if __name__ == "__main__":
    option = int(sys.argv[1]) if len(sys.argv) > 1 else None

    if option == 1:
        reference_id = sys.argv[2]
        metadata = Metadata(reference_id)
        metadata.run()

    elif option == 2:
        reference_id, fasta_content = sys.argv[2], sys.argv[3]
        alignment = Alignment(reference_id, fasta_content)
        alignment.run()

    elif option == 3:
        amino_acid_sequence = sys.argv[2]
        mRNADesign = MRNADesign(amino_acid_sequence)
        mRNADesign.run()

    elif option == 4:
        sequence = sys.argv[2]
        threeDViewer = ThreeDViewer(sequence)
        threeDViewer.run()
