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
        self.barcd_randBP_list = excel_arr[3]
        self.trgt_list = excel_arr[4]
        self.d0_seq_wo_scaf_list = excel_arr[5]
        self.barcd_list = excel_arr[6]

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
        return ''.join([i for i in align_arr[0] if not i.isdigit()]), ''.join(
            [i for i in align_arr[1] if not i.isdigit()]), ''.join(
            [i for i in align_arr[2] if not i.isdigit()]), alignments_result

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
            # logic_prep.get_real_scaf_seq 에서 end_idx == 0 인 경우 처리해서 필요없음
            # if (idx + 1) == 0:
            #     return len(seq_w_needle)
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

    def multi_filter_out_mismatch_seq_wo_rand_seq_in_brcd(self, fastq_list):
        print("multi_filter_out_mismatch_seq_wo_rand_seq_in_brcd >>>>>>>>>>>>>>>>>>>>>>")
        logic_prep = LogicPrep.LogicPreps()
        data_list = []
        err_list = []
        for idx in range(len(fastq_list)):
            ori_ngs_read = fastq_list[idx]

            estimate_scaf_ngs = ori_ngs_read[self.FRONT_SCAF_POS - self.FRONT_SCAF_WIN: self.FRONT_SCAF_POS + len(
                self.SCAFFOLD_SEQ) + self.REAR_SCAF_WIN]

            ngs_read_needle, needle_result, ref_seq_needle, alignments_result = self.get_pairwise2_localds_result(
                estimate_scaf_ngs, self.SCAFFOLD_SEQ)
            # ref_seq_needle 시작과 끝을 ngs_read_needle에서 scaffold의 시작과 끝으로 추정
            strt_needle_idx, end_needle_idx = self.count_hyphen_bfore_strt_aftr_end(ref_seq_needle)
            real_scaf_ngs = logic_prep.get_real_scaf_seq(ngs_read_needle, strt_needle_idx, end_needle_idx)
            strt_scaf_idx, end_scaf_idx = self.get_strt_end_idx(ori_ngs_read, real_scaf_ngs)

            guide_seq = self.get_target_seq(ori_ngs_read, strt_scaf_idx, self.LEN_GUIDE)

            umi_TTTG_brcd_seq = self.get_target_seq(ori_ngs_read, end_scaf_idx,
                                                         self.LEN_UMI + len(self.TTTG) + self.LEN_BRCD + self.LEN_RAND_WIN,
                                                         False)
            try:
                est_umi_seq_needle, umi_TTTG_needle_result, TTTG_needle, alignments_umi_TTTG = self.get_pairwise2_localds_result(
                    umi_TTTG_brcd_seq, self.TTTG)
            except Exception as err:
                err_list.append(["no_TTTG", 'no_TTTG', '', guide_seq, real_scaf_ngs, umi_TTTG_brcd_seq, "", "", "", "",
                                 ori_ngs_read])
                continue

            strt_TTTG_idx, end_TTTG_idx = self.get_strt_end_idx(umi_TTTG_brcd_seq,
                                                                est_umi_seq_needle.replace(" ", "").replace("-", ""))
            umi_seq = umi_TTTG_brcd_seq[:strt_TTTG_idx]
            seq_sftr_TTTG = umi_TTTG_brcd_seq[end_TTTG_idx:]
            brcd_seq = seq_sftr_TTTG[:self.LEN_BRCD]
            rand_rand_seq = seq_sftr_TTTG[self.LEN_BRCD:self.LEN_BRCD + self.LEN_RAND_BP]

            umi_TTTG_brcd_rand_trgt_seq = self.get_target_seq(ori_ngs_read, end_scaf_idx, self.LEN_UMI + len(
                self.TTTG) + self.LEN_BRCD + self.LEN_RAND_BP + self.LEN_RAND_WIN + self.LEN_TRGT, False)
            brcd_rand_end_idx = umi_TTTG_brcd_rand_trgt_seq.find(brcd_seq) + self.LEN_BRCD + self.LEN_RAND_BP
            trgt_seq = self.get_target_seq(umi_TTTG_brcd_rand_trgt_seq, brcd_rand_end_idx, self.LEN_TRGT, False)

            if self.SCAFFOLD_SEQ == real_scaf_ngs:
                try:
                    # check exist of barcode
                    brcd_idx = self.barcd_list.index(brcd_seq)
                    # check guide
                    if guide_seq == self.guide_list[brcd_idx]:
                        # Day 0 ==> check trgt_seq
                        if self.D0_D4_FLAG:
                            if trgt_seq == self.trgt_list[brcd_idx]:
                                data_list.append(
                                    ['', '', self.index_list[brcd_idx], guide_seq, real_scaf_ngs, umi_seq, brcd_seq,
                                     rand_rand_seq, self.barcd_randBP_list[brcd_idx][-3:], trgt_seq, ori_ngs_read])
                            else:
                                err_list.append(
                                    ["wrong_trgt", self.index_list[brcd_idx], '', guide_seq, real_scaf_ngs, umi_seq,
                                     brcd_seq, rand_rand_seq, self.barcd_randBP_list[brcd_idx][-3:], trgt_seq,
                                     ori_ngs_read])

                        # not Day 0
                        else:
                            data_list.append(
                                ['', '', self.index_list[brcd_idx], guide_seq, real_scaf_ngs, umi_seq, brcd_seq,
                                 rand_rand_seq, self.barcd_randBP_list[brcd_idx][-3:], trgt_seq, ori_ngs_read])
                    else:
                        try:
                            guide_idx = self.guide_list.index(guide_seq)
                            err_list.append(
                                ["wrong_guide", self.index_list[brcd_idx], self.index_list[guide_idx], guide_seq,
                                 real_scaf_ngs, umi_seq, brcd_seq, rand_rand_seq, self.barcd_randBP_list[brcd_idx][-3:],
                                 trgt_seq, ori_ngs_read])
                        except Exception as err:
                            err_list.append(
                                ["no_matched_guide", self.index_list[brcd_idx], '', guide_seq, real_scaf_ngs, umi_seq,
                                 brcd_seq, rand_rand_seq, self.barcd_randBP_list[brcd_idx][-3:], trgt_seq, ori_ngs_read])
                except Exception as err:
                    err_list.append(
                        ["no_brcd", 'no_brcd', '', guide_seq, real_scaf_ngs, umi_seq, brcd_seq, rand_rand_seq, '',
                         trgt_seq, ori_ngs_read])
            else:
                try:
                    brcd_idx = self.barcd_list.index(brcd_seq)
                    err_list.append(
                        ["wrong_scaffold", self.index_list[brcd_idx], '', guide_seq, real_scaf_ngs, umi_seq, brcd_seq,
                         rand_rand_seq, self.barcd_randBP_list[brcd_idx][-3:], trgt_seq, ori_ngs_read])
                except Exception as err:
                    if real_scaf_ngs in self.d0_seq_wo_scaf_list:
                        print(real_scaf_ngs, "::::::::::::::::::::::::::::::::::::::")
                    err_list.append(
                        ["wrong_scaffold", 'no_brcd', '', guide_seq, real_scaf_ngs, umi_seq, brcd_seq, rand_rand_seq, '',
                         trgt_seq, ori_ngs_read])

        return data_list, err_list

    def multi_filter_out_mismatch_seq_with_brcd_rand_seq(self, fastq_list):
        print("multi_filter_out_mismatch_seq_with_brcd_rand_seq >>>>>>>>>>>>>>>>>>>>>>")
        logic_prep = LogicPrep.LogicPreps()
        data_list = []
        no_TTTG_err_list = []
        wrong_trgt_err_list = []
        wrong_guide_err_list = []
        no_matched_guide_err_list = []
        no_brcd_err_list = []
        wrong_scaffold_err_list = []
        for idx in range(len(fastq_list)):
            ori_ngs_read = fastq_list[idx]

            estimate_scaf_ngs = ori_ngs_read[self.FRONT_SCAF_POS - self.FRONT_SCAF_WIN: self.FRONT_SCAF_POS + len(
                self.SCAFFOLD_SEQ) + self.REAR_SCAF_WIN]

            ngs_read_needle, needle_result, ref_seq_needle, alignments_result = self.get_pairwise2_localds_result(
                estimate_scaf_ngs, self.SCAFFOLD_SEQ)
            # ref_seq_needle 시작과 끝을 ngs_read_needle에서 scaffold의 시작과 끝으로 추정
            strt_needle_idx, end_needle_idx = self.count_hyphen_bfore_strt_aftr_end(ref_seq_needle)
            real_scaf_ngs = logic_prep.get_real_scaf_seq(ngs_read_needle, strt_needle_idx, end_needle_idx)
            strt_scaf_idx, end_scaf_idx = self.get_strt_end_idx(ori_ngs_read, real_scaf_ngs)

            guide_seq = self.get_target_seq(ori_ngs_read, strt_scaf_idx, self.LEN_GUIDE)

            umi_TTTG_brcd_rand_seq = self.get_target_seq(ori_ngs_read, end_scaf_idx,
                                                        self.LEN_UMI + len(
                                                            self.TTTG) + self.LEN_BRCD + self.LEN_RAND_BP + self.LEN_RAND_WIN,
                                                        False)
            try:
                est_umi_seq_needle, umi_TTTG_needle_result, TTTG_needle, alignments_umi_TTTG = self.get_pairwise2_localds_result(
                    umi_TTTG_brcd_rand_seq, self.TTTG)
            except Exception as err:
                no_TTTG_err_list.append(
                    ["no_TTTG", 'no_TTTG', '', guide_seq, real_scaf_ngs, umi_TTTG_brcd_rand_seq, "", "", "", ori_ngs_read])
                continue

            strt_TTTG_idx, end_TTTG_idx = self.get_strt_end_idx(umi_TTTG_brcd_rand_seq,
                                                                est_umi_seq_needle.replace(" ", "").replace("-", ""))
            umi_seq = umi_TTTG_brcd_rand_seq[:strt_TTTG_idx]
            seq_sftr_TTTG = umi_TTTG_brcd_rand_seq[end_TTTG_idx:]
            brcd_rand_seq = seq_sftr_TTTG[:self.LEN_BRCD + self.LEN_RAND_BP]

            umi_TTTG_brcd_rand_trgt_seq = self.get_target_seq(ori_ngs_read, end_scaf_idx, self.LEN_UMI + len(
                self.TTTG) + self.LEN_BRCD + self.LEN_RAND_BP + self.LEN_RAND_WIN + self.LEN_TRGT, False)
            brcd_rand_end_idx = umi_TTTG_brcd_rand_trgt_seq.find(brcd_rand_seq) + self.LEN_BRCD + self.LEN_RAND_BP
            trgt_seq = self.get_target_seq(umi_TTTG_brcd_rand_trgt_seq, brcd_rand_end_idx, self.LEN_TRGT, False)

            if self.SCAFFOLD_SEQ == real_scaf_ngs:
                try:
                    # check exist of barcode
                    brcd_rand_idx = self.barcd_randBP_list.index(brcd_rand_seq)
                    # check guide
                    if guide_seq == self.guide_list[brcd_rand_idx]:
                        # Day 0 ==> check trgt_seq
                        if self.D0_D4_FLAG:
                            if trgt_seq == self.trgt_list[brcd_rand_idx]:
                                data_list.append(
                                    ['', '', self.index_list[brcd_rand_idx], guide_seq, real_scaf_ngs, umi_seq,
                                     brcd_rand_seq[:-self.LEN_RAND_BP], brcd_rand_seq[-self.LEN_RAND_BP:], trgt_seq,
                                     ori_ngs_read])
                            else:
                                wrong_trgt_err_list.append(
                                    ["wrong_trgt", self.index_list[brcd_rand_idx], '', guide_seq, real_scaf_ngs, umi_seq,
                                     brcd_rand_seq[:-self.LEN_RAND_BP], brcd_rand_seq[-self.LEN_RAND_BP:], trgt_seq,
                                     ori_ngs_read])

                        # not Day 0
                        else:
                            data_list.append(['', '', self.index_list[brcd_rand_idx], guide_seq, real_scaf_ngs, umi_seq,
                                              brcd_rand_seq[:-self.LEN_RAND_BP], brcd_rand_seq[-self.LEN_RAND_BP:],
                                              trgt_seq, ori_ngs_read])
                    else:
                        try:
                            guide_idx = self.guide_list.index(guide_seq)
                            wrong_guide_err_list.append(
                                ["wrong_guide", self.index_list[brcd_rand_idx], self.index_list[guide_idx], guide_seq,
                                 real_scaf_ngs, umi_seq, brcd_rand_seq[:-self.LEN_RAND_BP],
                                 brcd_rand_seq[-self.LEN_RAND_BP:], trgt_seq, ori_ngs_read])
                        except Exception as err:
                            no_matched_guide_err_list.append(
                                ["no_matched_guide", self.index_list[brcd_rand_idx], '', guide_seq, real_scaf_ngs,
                                 umi_seq, brcd_rand_seq[:-self.LEN_RAND_BP], brcd_rand_seq[-self.LEN_RAND_BP:], trgt_seq,
                                 ori_ngs_read])
                except Exception as err:
                    no_brcd_err_list.append(
                        ["no_brcd", 'no_brcd', '', guide_seq, real_scaf_ngs, umi_seq, brcd_rand_seq[:-self.LEN_RAND_BP],
                         brcd_rand_seq[-self.LEN_RAND_BP:], trgt_seq, ori_ngs_read])
            else:
                try:
                    brcd_rand_idx = self.barcd_randBP_list.index(brcd_rand_seq)
                    wrong_scaffold_err_list.append(
                        ["wrong_scaffold", self.index_list[brcd_rand_idx], '', guide_seq, real_scaf_ngs, umi_seq,
                         brcd_rand_seq[:-self.LEN_RAND_BP], brcd_rand_seq[-self.LEN_RAND_BP:], trgt_seq, ori_ngs_read])
                except Exception as err:
                    if real_scaf_ngs in self.d0_seq_wo_scaf_list:
                        print(real_scaf_ngs, "::::::::::::::::::::::::::::::::::::::")
                    wrong_scaffold_err_list.append(
                        ["wrong_scaffold", 'no_brcd', '', guide_seq, real_scaf_ngs, umi_seq,
                         brcd_rand_seq[:-self.LEN_RAND_BP],
                         brcd_rand_seq[-self.LEN_RAND_BP:], trgt_seq, ori_ngs_read])

            return data_list, [no_TTTG_err_list, wrong_trgt_err_list, wrong_guide_err_list, no_matched_guide_err_list,
                               no_brcd_err_list, wrong_scaffold_err_list]

    def multi_filter_out_mismatch_seq_by_scaffold_existence(self, fastq_list):
        print("multi_filter_out_mismatch_seq_by_scaffold_existence >>>>>>>>>>>>>>>>>>>>>>")
        logic_prep = LogicPrep.LogicPreps()
        data_list = []
        no_TTTG_err_list = []
        wrong_trgt_err_list = []
        wrong_guide_err_list = []
        no_matched_guide_err_list = []
        no_brcd_err_list = []
        wrong_scaffold_err_list = []
        for idx in range(len(fastq_list)):
            ori_ngs_read = fastq_list[idx]

            if self.SCAFFOLD_SEQ in ori_ngs_read:
                real_scaf_ngs = self.SCAFFOLD_SEQ
                strt_scaf_idx, end_scaf_idx = self.get_strt_end_idx(ori_ngs_read, real_scaf_ngs)

                guide_seq = self.get_target_seq(ori_ngs_read, strt_scaf_idx, self.LEN_GUIDE)

                umi_TTTG_brcd_rand_seq = self.get_target_seq(ori_ngs_read, end_scaf_idx, self.LEN_UMI + len(
                    self.TTTG) + self.LEN_BRCD + self.LEN_RAND_BP + self.LEN_RAND_WIN, False)
                try:
                    est_umi_seq_needle, umi_TTTG_needle_result, TTTG_needle, alignments_umi_TTTG = self.get_pairwise2_localds_result(
                        umi_TTTG_brcd_rand_seq, self.TTTG)
                except Exception as err:
                    no_TTTG_err_list.append(
                        ["no_TTTG", 'no_TTTG', '', guide_seq, real_scaf_ngs, umi_TTTG_brcd_rand_seq, "", "", "",
                         ori_ngs_read])
                    continue

                strt_TTTG_idx, end_TTTG_idx = self.get_strt_end_idx(umi_TTTG_brcd_rand_seq,
                                                                    est_umi_seq_needle.replace(" ", "").replace("-",
                                                                                                                ""))
                umi_seq = umi_TTTG_brcd_rand_seq[:strt_TTTG_idx]
                seq_sftr_TTTG = umi_TTTG_brcd_rand_seq[end_TTTG_idx:]
                brcd_rand_seq = seq_sftr_TTTG[:self.LEN_BRCD + self.LEN_RAND_BP]

                umi_TTTG_brcd_rand_trgt_seq = self.get_target_seq(ori_ngs_read, end_scaf_idx, self.LEN_UMI + len(
                    self.TTTG) + self.LEN_BRCD + self.LEN_RAND_BP + self.LEN_RAND_WIN + self.LEN_TRGT, False)
                brcd_rand_end_idx = umi_TTTG_brcd_rand_trgt_seq.find(brcd_rand_seq) + self.LEN_BRCD + self.LEN_RAND_BP
                trgt_seq = self.get_target_seq(umi_TTTG_brcd_rand_trgt_seq, brcd_rand_end_idx, self.LEN_TRGT, False)

                try:
                    # check exist of barcode
                    brcd_rand_idx = self.barcd_randBP_list.index(brcd_rand_seq)
                    # check guide
                    if guide_seq == self.guide_list[brcd_rand_idx]:
                        # Day 0 ==> check trgt_seq
                        if self.D0_D4_FLAG:
                            if trgt_seq == self.trgt_list[brcd_rand_idx]:
                                data_list.append(
                                    ['', '', self.index_list[brcd_rand_idx], guide_seq, real_scaf_ngs, umi_seq,
                                     brcd_rand_seq[:-self.LEN_RAND_BP], brcd_rand_seq[-self.LEN_RAND_BP:], trgt_seq,
                                     ori_ngs_read])
                            else:
                                wrong_trgt_err_list.append(
                                    ["wrong_trgt", self.index_list[brcd_rand_idx], '', guide_seq, real_scaf_ngs, umi_seq,
                                     brcd_rand_seq[:-self.LEN_RAND_BP], brcd_rand_seq[-self.LEN_RAND_BP:], trgt_seq,
                                     ori_ngs_read])

                        # not Day 0
                        else:
                            data_list.append(['', '', self.index_list[brcd_rand_idx], guide_seq, real_scaf_ngs, umi_seq,
                                              brcd_rand_seq[:-self.LEN_RAND_BP], brcd_rand_seq[-self.LEN_RAND_BP:],
                                              trgt_seq, ori_ngs_read])
                    else:
                        try:
                            guide_idx = self.guide_list.index(guide_seq)
                            wrong_guide_err_list.append(
                                ["wrong_guide", self.index_list[brcd_rand_idx], self.index_list[guide_idx], guide_seq,
                                 real_scaf_ngs, umi_seq, brcd_rand_seq[:-self.LEN_RAND_BP],
                                 brcd_rand_seq[-self.LEN_RAND_BP:], trgt_seq, ori_ngs_read])
                        except Exception as err:
                            no_matched_guide_err_list.append(
                                ["no_matched_guide", self.index_list[brcd_rand_idx], '', guide_seq, real_scaf_ngs,
                                 umi_seq, brcd_rand_seq[:-self.LEN_RAND_BP], brcd_rand_seq[-self.LEN_RAND_BP:], trgt_seq,
                                 ori_ngs_read])

                except Exception as err:
                    no_brcd_err_list.append(
                        ["no_brcd", 'no_brcd', '', guide_seq, real_scaf_ngs, umi_seq, brcd_rand_seq[:-self.LEN_RAND_BP],
                         brcd_rand_seq[-self.LEN_RAND_BP:], trgt_seq, ori_ngs_read])

            else:
                wrong_scaffold_err_list.append(
                    ["no_scaffold", '', '', '', '', '', '', '', '', ori_ngs_read])

        return data_list, [no_TTTG_err_list, wrong_trgt_err_list, wrong_guide_err_list, no_matched_guide_err_list,
                               no_brcd_err_list, wrong_scaffold_err_list]