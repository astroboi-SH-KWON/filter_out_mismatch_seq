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

    def exist_exact_seq(self, needle_result):
        needles = needle_result.strip()
        if " " in needles:
            return False
        elif "." in needles:
            return False
        return True

    def count_hyphen_bfore_strt(self, seq_w_needle, idx):
        if seq_w_needle[idx] == "-":
            return self.count_hyphen_bfore_strt(seq_w_needle, idx + 1)
        else:
            return idx

    def count_hyphen_aftr_end(self, seq_w_needle, idx):
        if seq_w_needle[idx] == "-":
            return self.count_hyphen_aftr_end(seq_w_needle, idx - 1)
        else:
            return idx + 1

    def count_hyphen_bfore_strt_aftr_end(self, seq_w_needle):
        return self.count_hyphen_bfore_strt(seq_w_needle, 0), self.count_hyphen_aftr_end(seq_w_needle, -1)

    def get_strt_end_idx(self, full_seq, trgt_seq):
        start_idx = full_seq.find(trgt_seq)
        return start_idx, start_idx + len(trgt_seq)

    def get_target_seq(self, full_seq, trgt_idx, len_trgt, flag=True):
        if flag:
            return full_seq[trgt_idx - len_trgt: trgt_idx]
        else:
            return full_seq[trgt_idx: trgt_idx + len_trgt]

