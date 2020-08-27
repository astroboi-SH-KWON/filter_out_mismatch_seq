from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62

class Logics:
    def __init__(self):
        pass

    """
    by using the BLOSUM62 matrix, together with a gap open penalty of 10 and a gap extension penalty of 0.5 (using globalds)
    """
    def get_pairwise2_needle_result(self, asequence, bsequence, matrx=blosum62, gap_open_penalty=10, extension_penalty=0.5):
        alignments = pairwise2.align.globalds(asequence.upper().replace(" ", ""), bsequence.upper().replace(" ", ""),
                                              matrx, -gap_open_penalty, -extension_penalty)
        alignments_result = pairwise2.format_alignment(*alignments[0])
        align_arr = alignments_result.split("\n")
        return align_arr[0], align_arr[1], align_arr[2], alignments_result

    """
    LSPADKTNVKAA
      |..|..|
    --PEEKSAV---
    Score=16
    <BLANKLINE>
    """
    def get_pairwise2_localds_result(self, asequence, bsequence, matrx=blosum62, gap_open_penalty=10, extension_penalty=1):
        alignments = pairwise2.align.localds(asequence.upper().replace(" ", ""), bsequence.upper().replace(" ", ""),
                                              matrx, -gap_open_penalty, -extension_penalty)
        alignments_result = pairwise2.format_alignment(*alignments[0])
        align_arr = alignments_result.split("\n")
        return align_arr[0], align_arr[1], align_arr[2], alignments_result