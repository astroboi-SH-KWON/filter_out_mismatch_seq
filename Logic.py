from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62

import LogicPrep
class Logics:
    def __init__(self, init, excel_arr, D0_D4_FLAG):
        self.SCAFFOLD_SEQ = init[0]
        self.FRONT_SCAF = init[1]
        self.FRONT_SCAF_WIN = init[2]
        self.FRONT_SCAF_POS = init[3]
        self.REAR_SCAF_WIN = init[4]

        self.LEN_GUIDE = init[5]
        self.LEN_UMI = init[6]
        self.TTTG = init[7]
        self.LEN_BRCD = init[8]
        self.LEN_RAND_BP = init[9]
        self.LEN_RAND_WIN = init[10]
        self.LEN_TRGT = init[11]
        self.D0_D4_FLAG = D0_D4_FLAG

        self.csv_list = excel_arr[0]
        self.index_list = excel_arr[1]
        self.guide_list = excel_arr[2]
        self.barcd_3bp_list = excel_arr[3]
        self.trgt_list = excel_arr[4]
        self.d0_seq_wo_scaf_list = excel_arr[5]

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

    def multi_filter_out_mismatch_seq(self, fastq_list):
        logic_prep = LogicPrep.LogicPreps()
        data_list = []
        err_list = []
        for idx in range(len(fastq_list)):
            ori_ngs_read = fastq_list[idx]
            estimate_scaf_ngs = ori_ngs_read[
                                self.FRONT_SCAF_POS - self.FRONT_SCAF_WIN: self.FRONT_SCAF_POS + len(self.SCAFFOLD_SEQ) + self.REAR_SCAF_WIN]
            ngs_read_needle, needle_result, ref_seq_needle, alignments_result = self.get_pairwise2_localds_result(
                estimate_scaf_ngs, self.SCAFFOLD_SEQ)

            strt_needle_idx, end_needle_idx = self.count_hyphen_bfore_strt_aftr_end(ref_seq_needle)
            real_scaf_ngs = logic_prep.get_real_scaf_seq(ngs_read_needle, strt_needle_idx, end_needle_idx)
            strt_scaf_idx, end_scaf_idx = self.get_strt_end_idx(ori_ngs_read, real_scaf_ngs)

            guide_seq = self.get_target_seq(ori_ngs_read, strt_scaf_idx, self.LEN_GUIDE)

            umi_TTTG_brcd_3bp_seq = self.get_target_seq(ori_ngs_read, end_scaf_idx,
                                                         self.LEN_UMI + len(self.TTTG) + self.LEN_BRCD + self.LEN_RAND_BP + self.LEN_RAND_WIN,
                                                         False)

            umi_seq_needle, umi_TTTG_needle_result, TTTG_needle, alignments_umi_TTTG = self.get_pairwise2_localds_result(
                umi_TTTG_brcd_3bp_seq, self.TTTG)

            strt_TTTG_idx, end_TTTG_idx = self.count_hyphen_bfore_strt_aftr_end(TTTG_needle)

            umi_seq = umi_TTTG_brcd_3bp_seq[:strt_TTTG_idx]
            seq_sftr_TTTG = umi_TTTG_brcd_3bp_seq[end_TTTG_idx:]
            brcd_3bp_seq = seq_sftr_TTTG[:self.LEN_BRCD + self.LEN_RAND_BP]

            umi_TTTG_brcd_3bp_trgt_seq = self.get_target_seq(ori_ngs_read, end_scaf_idx, self.LEN_UMI + len(
                self.TTTG) + self.LEN_BRCD + self.LEN_RAND_BP + self.LEN_RAND_WIN + self.LEN_TRGT, False)
            brcd_3bp_end_idx = umi_TTTG_brcd_3bp_trgt_seq.find(brcd_3bp_seq) + self.LEN_BRCD + self.LEN_RAND_BP
            trgt_seq = self.get_target_seq(umi_TTTG_brcd_3bp_trgt_seq, brcd_3bp_end_idx, self.LEN_TRGT, False)

            if self.exist_exact_seq(needle_result):
                try:
                    # check exist of barcode
                    brcd_3bp_idx = self.barcd_3bp_list.index(brcd_3bp_seq)
                    # check guide
                    if guide_seq == self.guide_list[brcd_3bp_idx]:
                        # Day 0 ==> check trgt_seq
                        if self.D0_D4_FLAG:
                            if trgt_seq == self.trgt_list[brcd_3bp_idx]:
                                data_list.append(
                                    ["", self.index_list[brcd_3bp_idx], '', guide_seq, real_scaf_ngs, umi_seq, brcd_3bp_seq,
                                     trgt_seq, ori_ngs_read])
                            else:
                                err_list.append(
                                    ["wrong_trgt", self.index_list[brcd_3bp_idx], '', guide_seq, real_scaf_ngs, umi_seq,
                                     brcd_3bp_seq, trgt_seq, ori_ngs_read])

                        # not Day 0
                        else:
                            data_list.append(
                                ["", self.index_list[brcd_3bp_idx], '', guide_seq, real_scaf_ngs, umi_seq, brcd_3bp_seq,
                                 trgt_seq, ori_ngs_read])
                    else:
                        try:
                            guide_idx = self.guide_list.index(guide_seq)
                            print("guide_idx : ", guide_idx)
                            err_list.append(["wrong_guide", self.index_list[brcd_3bp_idx], self.index_list[guide_idx], guide_seq,
                                             real_scaf_ngs, umi_seq, brcd_3bp_seq, trgt_seq, ori_ngs_read])
                        except Exception as err:
                            print("brcd_3bp_seq : ", brcd_3bp_seq, "guide_seq : ", guide_seq)
                            err_list.append(
                                ["no_matched_guide", self.index_list[brcd_3bp_idx], '', guide_seq, real_scaf_ngs, umi_seq,
                                 brcd_3bp_seq, trgt_seq, ori_ngs_read])
                except Exception as err:
                    err_list.append(
                        ["no_brcd", 'no_brcd', '', guide_seq, real_scaf_ngs, umi_seq, brcd_3bp_seq, trgt_seq,
                         ori_ngs_read])
            else:
                try:
                    brcd_3bp_idx = self.barcd_3bp_list.index(brcd_3bp_seq)
                    err_list.append(["wrong_scaffold", self.index_list[brcd_3bp_idx], '', guide_seq, real_scaf_ngs, umi_seq,
                                     brcd_3bp_seq, trgt_seq, ori_ngs_read])
                except Exception as err:
                    if real_scaf_ngs in self.d0_seq_wo_scaf_list:
                        print(real_scaf_ngs, "::::::::::::::::::::::::::::::::::::::")
                    err_list.append(
                        ["wrong_scaffold", 'no_brcd', '', guide_seq, real_scaf_ngs, umi_seq, brcd_3bp_seq, trgt_seq,
                         ori_ngs_read])

        return data_list, err_list