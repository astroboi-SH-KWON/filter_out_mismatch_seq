import time
import os

import Util
import Logic
import LogicPrep
############### start to set env ################
WORK_DIR = "D:/000_WORK/KimNahye/20200827/WORK_DIR/"
# WORK_DIR = os.getcwd() + "/"
PROJECT_NAME = WORK_DIR.split("/")[-2]
FASTQ = "FASTQ/200302_PCR switching_hi-seq/"
INPUT = "input/"
GUIDE_BARCODE_CSV = "190509_FINAL.CSV"
D0_Lib_10fg = [1, 4]
D4_Gen_10ng = [2, 3, 5, 6, 7, 8]
D0_D4_FLAG_ARR = [True, False]
FASTQ_ARR = [D0_Lib_10fg, D4_Gen_10ng]
FASTQ_EXT = ".extendedFrags.fastq"

SCAFFOLD_SEQ = "GTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT"
FRONT_SCAF = SCAFFOLD_SEQ[:5]
FRONT_SCAF_WIN = 2
FRONT_SCAF_POS = 50
REAR_SCAF_WIN = 10

LEN_GUIDE = 19
LEN_UMI = 8
TTTG = "TTTG"
LEN_BRCD = 15
LEN_RAND_BP = 3
LEN_RAND_WIN = 3
LEN_TRGT = 24
############### end setting env #################


def main():
    util = Util.Utils()
    logic_prep = LogicPrep.LogicPreps()

    csv_list = util.read_csv_ignore_N_line(WORK_DIR + INPUT + GUIDE_BARCODE_CSV)
    index_list = logic_prep.make_1_arr_list_to_list(0, csv_list)
    guide_list = logic_prep.make_1_arr_list_to_list(2, csv_list)
    barcd_3bp_list = logic_prep.make_2_arr_list_to_list(6, 7, csv_list)
    trgt_list = logic_prep.make_1_arr_list_to_list(8, csv_list)
    d0_seq_wo_scaf_list = logic_prep.make_3_arr_list_to_list(3, 4, 5, csv_list)

    logic = Logic.Logics([], [])

    for d0_d4_idx in range(len(D0_D4_FLAG_ARR)):
        err_list = []
        data_list = []

        for fn_nm in FASTQ_ARR[d0_d4_idx]:
            fastq_list = util.read_fastq_to_list(WORK_DIR + FASTQ + str(fn_nm) + FASTQ_EXT)

            for idx in range(len(fastq_list)):
                ori_ngs_read = fastq_list[idx]
                estimate_scaf_ngs = ori_ngs_read[
                                    FRONT_SCAF_POS - FRONT_SCAF_WIN: FRONT_SCAF_POS + len(SCAFFOLD_SEQ) + REAR_SCAF_WIN]
                ngs_read_needle, needle_result, ref_seq_needle, alignments_result = logic.get_pairwise2_localds_result(
                    estimate_scaf_ngs, SCAFFOLD_SEQ)

                strt_needle_idx, end_needle_idx = logic.count_hyphen_bfore_strt_aftr_end(ref_seq_needle)
                real_scaf_ngs = logic_prep.get_real_scaf_seq(ngs_read_needle, strt_needle_idx, end_needle_idx)
                strt_scaf_idx, end_scaf_idx = logic.get_strt_end_idx(ori_ngs_read, real_scaf_ngs)

                guide_seq = logic.get_target_seq(ori_ngs_read, strt_scaf_idx, LEN_GUIDE)

                umi_TTTG_brcd_3bp_seq = logic.get_target_seq(ori_ngs_read, end_scaf_idx,
                                                             LEN_UMI + len(TTTG) + LEN_BRCD + LEN_RAND_BP + LEN_RAND_WIN,
                                                             False)

                umi_seq_needle, umi_TTTG_needle_result, TTTG_needle, alignments_umi_TTTG = logic.get_pairwise2_localds_result(
                    umi_TTTG_brcd_3bp_seq, TTTG)

                strt_TTTG_idx, end_TTTG_idx = logic.count_hyphen_bfore_strt_aftr_end(TTTG_needle)

                umi_seq = umi_TTTG_brcd_3bp_seq[:strt_TTTG_idx]
                seq_sftr_TTTG = umi_TTTG_brcd_3bp_seq[end_TTTG_idx:]
                brcd_3bp_seq = seq_sftr_TTTG[:LEN_BRCD + LEN_RAND_BP]

                umi_TTTG_brcd_3bp_trgt_seq = logic.get_target_seq(ori_ngs_read, end_scaf_idx, LEN_UMI + len(
                    TTTG) + LEN_BRCD + LEN_RAND_BP + LEN_RAND_WIN + LEN_TRGT, False)
                brcd_3bp_end_idx = umi_TTTG_brcd_3bp_trgt_seq.find(brcd_3bp_seq) + LEN_BRCD + LEN_RAND_BP
                trgt_seq = logic.get_target_seq(umi_TTTG_brcd_3bp_trgt_seq, brcd_3bp_end_idx, LEN_TRGT, False)

                if logic.exist_exact_seq(needle_result):
                    try:
                        # check exist of barcode
                        brcd_3bp_idx = barcd_3bp_list.index(brcd_3bp_seq)
                        # check guide
                        if guide_seq == guide_list[brcd_3bp_idx]:
                            # Day 0 ==> check trgt_seq
                            if D0_D4_FLAG_ARR[d0_d4_idx]:
                                if trgt_seq == trgt_list[brcd_3bp_idx]:
                                    data_list.append(
                                        ["", index_list[brcd_3bp_idx], '', guide_seq, real_scaf_ngs, umi_seq, brcd_3bp_seq,
                                         trgt_seq, ori_ngs_read])
                                else:
                                    err_list.append(
                                        ["wrong_trgt", index_list[brcd_3bp_idx], '', guide_seq, real_scaf_ngs, umi_seq,
                                         brcd_3bp_seq, trgt_seq, ori_ngs_read])

                            # not Day 0
                            else:
                                data_list.append(
                                    ["", index_list[brcd_3bp_idx], '', guide_seq, real_scaf_ngs, umi_seq, brcd_3bp_seq,
                                     trgt_seq, ori_ngs_read])
                        else:
                            try:
                                guide_idx = guide_list.index(guide_seq)
                                print("guide_idx : ", guide_idx)
                                err_list.append(["wrong_guide", index_list[brcd_3bp_idx], index_list[guide_idx], guide_seq,
                                                 real_scaf_ngs, umi_seq, brcd_3bp_seq, trgt_seq, ori_ngs_read])
                            except Exception as err:
                                print("brcd_3bp_seq : ", brcd_3bp_seq, "guide_seq : ", guide_seq)
                                err_list.append(
                                    ["no_matched_guide", index_list[brcd_3bp_idx], '', guide_seq, real_scaf_ngs, umi_seq,
                                     brcd_3bp_seq, trgt_seq, ori_ngs_read])
                    except Exception as err:
                        err_list.append(
                            ["no_brcd", 'no_brcd', '', guide_seq, real_scaf_ngs, umi_seq, brcd_3bp_seq, trgt_seq,
                             ori_ngs_read])
                else:
                    try:
                        brcd_3bp_idx = barcd_3bp_list.index(brcd_3bp_seq)
                        err_list.append(["wrong_scaffold", index_list[brcd_3bp_idx], '', guide_seq, real_scaf_ngs, umi_seq,
                                         brcd_3bp_seq, trgt_seq, ori_ngs_read])
                    except Exception as err:
                        if real_scaf_ngs in d0_seq_wo_scaf_list:
                            print(real_scaf_ngs, "::::::::::::::::::::::::::::::::::::::")
                        err_list.append(
                            ["wrong_scaffold", 'no_brcd', '', guide_seq, real_scaf_ngs, umi_seq, brcd_3bp_seq, trgt_seq,
                             ori_ngs_read])

        err_head = ['error_code', 'expected_index', 'index_from_NGS', 'guide_NGS', 'scaf_NGS', 'umi', 'barcode_3bp',
                    'target_NGS', 'full_NGS']
        result_head = err_head[1:]
        util.make_excel(WORK_DIR + "output/result_" + str(d0_d4_idx), result_head, data_list, 1)

        util.make_excel(WORK_DIR + "output/err_" + str(d0_d4_idx), err_head, data_list)



if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    main()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))