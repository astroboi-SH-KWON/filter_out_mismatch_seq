import time
import os
import multiprocessing as mp
import numpy as np

import Util
import Logic
import LogicPrep
############### start to set env ################
# WORK_DIR = "D:/000_WORK/KimNahye/20200827/WORK_DIR/"
WORK_DIR = os.getcwd() + "/"
PROJECT_NAME = WORK_DIR.split("/")[-2]
FASTQ = "FASTQ/20201019/"
# FASTQ = "FASTQ/200302_PCR switching_hi-seq/"
INPUT = "input/"
# GUIDE_BARCODE_CSV = "190509_FINAL.CSV"
GUIDE_BARCODE_CSV = "1st_LibraryA.CSV"

LibA_D0 = ['LibA_D0']
WT = ['WT']
# D0_Lib_10fg = [1, 4]
# D4_Gen_10ng = [2, 3, 5, 6, 7, 8]
D0_D4_FLAG_ARR = [True, False]  # Day 0 : True, non Day 0 : False
# D0_D4_FLAG_ARR = [True, False]
FASTQ_ARR = [LibA_D0, WT]
# FASTQ_ARR = [D0_Lib_10fg, D4_Gen_10ng]
FASTQ_N = ['LibA_D0', 'WT']
# FASTQ_N = ['D0_Lib_10fg', 'D4_Gen_10ng']
FASTQ_EXT = ".extendedFrags.fastq"

SCAFFOLD_SEQ = "GTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT"
FRONT_SCAF = SCAFFOLD_SEQ[:5]
FRONT_SCAF_WIN = 2
FRONT_SCAF_POS = 50
REAR_SCAF_WIN = 10

LEN_GUIDE = 19
LEN_UMI = 8
TTTG = "TTTG"
TTG = "TTTG"
LEN_BRCD = 15
# LEN_RAND_BP = 3
LEN_RAND_BP = 9
POS_SLICE_RAND_BP = 3
LEN_RAND_WIN = 3
LEN_PAM = 6
LEN_TRGT = 24 + LEN_PAM

INIT = [SCAFFOLD_SEQ, FRONT_SCAF, FRONT_SCAF_WIN, FRONT_SCAF_POS, REAR_SCAF_WIN, LEN_GUIDE, LEN_UMI, TTTG, LEN_BRCD, LEN_RAND_BP, LEN_RAND_WIN, LEN_TRGT, TTG, POS_SLICE_RAND_BP]

TOTAL_CPU = mp.cpu_count()
MULTI_CNT = int(TOTAL_CPU*0.8)
############### end setting env #################

def multi_processing():
    util = Util.Utils()
    logic_prep = LogicPrep.LogicPreps()

    excel_arr = []
    csv_list = [[x.upper() for x in tmp_arr] for tmp_arr in
                util.read_csv_ignore_N_line(WORK_DIR + INPUT + GUIDE_BARCODE_CSV)]
    # csv_list
    excel_arr.append(csv_list)
    # index_list
    excel_arr.append(logic_prep.make_1_arr_list_to_list(0, csv_list))
    # guide_list
    excel_arr.append(logic_prep.make_1_arr_list_to_list(2, csv_list))
    # barcd_randBP_list
    excel_arr.append(logic_prep.make_2_arr_list_to_list(6, 7, csv_list))
    # trgt_list
    excel_arr.append(logic_prep.make_1_arr_list_to_list(8, csv_list))
    # d0_seq_wo_scaf_list
    excel_arr.append(logic_prep.make_3_arr_list_to_list(3, 4, 5, csv_list))
    # barcd_list
    excel_arr.append(logic_prep.make_1_arr_list_to_list(6, csv_list))

    for d0_d4_idx in range(len(D0_D4_FLAG_ARR)):
        logic = Logic.Logics(INIT, excel_arr, D0_D4_FLAG_ARR[d0_d4_idx])

        for fn_nm in FASTQ_ARR[d0_d4_idx]:
            fastq_list = util.read_fastq_to_list(WORK_DIR + FASTQ + str(fn_nm) + FASTQ_EXT)

            splited_fastq_list = np.array_split(fastq_list, MULTI_CNT)

            print("total cpu_count : " + str(TOTAL_CPU))
            print("will use : " + str(MULTI_CNT))
            pool = mp.Pool(processes=MULTI_CNT)

            # pool_list = pool.map(logic.filter_out_mismatch_seq_with_brcd_rand_seq, splited_fastq_list)
            pool_list = pool.map(logic.filter_out_mismatch_seq_by_scaffold_existence1, splited_fastq_list)

            data_list, err_list = util.merge_multi_err_list(pool_list)
            pool.close()

            head = ['error_code', 'expected_index', 'index_from_NGS', 'guide_NGS', 'scaf_NGS', 'umi', 'barcode',
                    'rand_' + str(LEN_RAND_BP) + '_bp', 'target_NGS', 'full_NGS']
            util.make_excel(WORK_DIR + "output/" + str(fn_nm) + "_result_" + FASTQ_N[d0_d4_idx], head, data_list, 2)

            for err_arr in err_list:
                if len(err_arr) == 0:
                    continue

                try:
                    util.make_excel(
                        WORK_DIR + "output/" + str(fn_nm) + "_err_" + err_arr[0][0] + "_" + FASTQ_N[d0_d4_idx], head,
                        err_arr)
                except Exception as err:
                    util.make_tsv(
                        WORK_DIR + "output/" + str(fn_nm) + "_err_" + err_arr[0][0] + "_" + FASTQ_N[d0_d4_idx], head,
                        err_arr)

def multi_processing_plan_B():
    util = Util.Utils()
    logic_prep = LogicPrep.LogicPreps()

    excel_arr = []
    csv_list = [[x.upper() for x in tmp_arr] for tmp_arr in
                util.read_csv_ignore_N_line(WORK_DIR + INPUT + GUIDE_BARCODE_CSV)]
    # csv_list
    excel_arr.append(csv_list)
    # index_list
    excel_arr.append(logic_prep.make_1_arr_list_to_list(0, csv_list))
    # guide_list
    excel_arr.append(logic_prep.make_1_arr_list_to_list(2, csv_list))
    # barcd_randBP_list
    excel_arr.append(logic_prep.make_2_arr_list_to_list_after_slice(6, 7, POS_SLICE_RAND_BP, csv_list))
    # trgt_list : 20201026 trgt + PAM
    excel_arr.append(logic_prep.make_2_arr_list_to_list(8, 9, csv_list))
    # d0_seq_wo_scaf_list
    excel_arr.append(logic_prep.make_3_arr_list_to_list(3, 4, 5, csv_list))
    # barcd_list
    excel_arr.append(logic_prep.make_1_arr_list_to_list(6, csv_list))
    # randBP_list
    excel_arr.append(logic_prep.make_1_arr_list_to_list(7, csv_list))

    for d0_d4_idx in range(len(D0_D4_FLAG_ARR)):
        logic = Logic.Logics(INIT, excel_arr, D0_D4_FLAG_ARR[d0_d4_idx])

        for fn_nm in FASTQ_ARR[d0_d4_idx]:
            fastq_list = util.read_fastq_to_list(WORK_DIR + FASTQ + str(fn_nm) + FASTQ_EXT)

            splited_fastq_list = np.array_split(fastq_list, MULTI_CNT)

            print("total cpu_count : " + str(TOTAL_CPU))
            print("will use : " + str(MULTI_CNT))
            pool = mp.Pool(processes=MULTI_CNT)

            # pool_list = pool.map(logic.filter_out_mismatch_seq_by_scaffold_existence2, splited_fastq_list)
            pool_list = pool.map(logic.filter_out_mismatch_seq_by_full_brcd_seq, splited_fastq_list)

            data_list, err_list = util.merge_multi_err_list(pool_list)
            pool.close()

            head = ['error_code', 'expected_index', 'index_from_NGS', 'guide_NGS', 'scaf_NGS',
                    'umi_with_extra_T_or_not', TTG + '_barcode', 'rand_' + str(LEN_RAND_BP) + '_bp',
                    'target_NGS + PAM(' + str(LEN_PAM) + ')', 'full_NGS']
            util.make_excel(WORK_DIR + "output/" + str(fn_nm) + "_result_" + FASTQ_N[d0_d4_idx], head, data_list, 2)

            for err_arr in err_list:
                if len(err_arr) == 0:
                    continue

                try:
                    util.make_excel(
                        WORK_DIR + "output/" + str(fn_nm) + "_err_" + err_arr[0][0] + "_" + FASTQ_N[d0_d4_idx], head,
                        err_arr)
                except Exception as err:
                    util.make_tsv(
                        WORK_DIR + "output/" + str(fn_nm) + "_err_" + err_arr[0][0] + "_" + FASTQ_N[d0_d4_idx], head,
                        err_arr)


if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    # multi_processing()
    multi_processing_plan_B()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))