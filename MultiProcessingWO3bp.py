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
FASTQ = "FASTQ/200302_PCR switching_hi-seq/"
INPUT = "input/"
GUIDE_BARCODE_CSV = "190509_FINAL.CSV"
D0_Lib_10fg = [1, 4]
D4_Gen_10ng = [2, 3, 5, 6, 7, 8]
D0_D4_FLAG_ARR = [True, False]
FASTQ_ARR = [D0_Lib_10fg, D4_Gen_10ng]
FASTQ_N = ['D0_Lib_10fg', 'D4_Gen_10ng']
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

INIT = [SCAFFOLD_SEQ, FRONT_SCAF, FRONT_SCAF_WIN, FRONT_SCAF_POS, REAR_SCAF_WIN, LEN_GUIDE, LEN_UMI, TTTG, LEN_BRCD, LEN_RAND_BP, LEN_RAND_WIN, LEN_TRGT]

TOTAL_CPU = mp.cpu_count()
MULTI_CNT = int(TOTAL_CPU*0.8)
############### end setting env #################

def multi_processing_wo_randbp_in_brcd():
    util = Util.Utils()
    logic_prep = LogicPrep.LogicPreps()

    excel_arr = []
    csv_list = util.read_csv_ignore_N_line(WORK_DIR + INPUT + GUIDE_BARCODE_CSV)
    excel_arr.append(csv_list)
    excel_arr.append(logic_prep.make_1_arr_list_to_list(0, csv_list))
    excel_arr.append(logic_prep.make_1_arr_list_to_list(2, csv_list))
    excel_arr.append(logic_prep.make_2_arr_list_to_list(6, 7, csv_list))
    excel_arr.append(logic_prep.make_1_arr_list_to_list(8, csv_list))
    excel_arr.append(logic_prep.make_3_arr_list_to_list(3, 4, 5, csv_list))
    excel_arr.append(logic_prep.make_1_arr_list_to_list(6, csv_list))

    for d0_d4_idx in range(len(D0_D4_FLAG_ARR)):
        logic = Logic.Logics(INIT, excel_arr, D0_D4_FLAG_ARR[d0_d4_idx])

        for fn_nm in FASTQ_ARR[d0_d4_idx]:
            fastq_list = util.read_fastq_to_list(WORK_DIR + FASTQ + str(fn_nm) + FASTQ_EXT)

            splited_fastq_list = np.array_split(fastq_list, MULTI_CNT)

            print("total cpu_count : " + str(TOTAL_CPU))
            print("will use : " + str(MULTI_CNT))
            pool = mp.Pool(processes=MULTI_CNT)

            pool_list = pool.map(logic.filter_out_mismatch_seq_with_brcd_rand_seq, splited_fastq_list)

            data_list, err_list = util.merge_multi_list(pool_list)
            pool.close()

            head = ['error_code', 'expected_index', 'index_from_NGS', 'guide_NGS', 'scaf_NGS', 'umi', 'barcode',
                    'rand_3bp_ngs', 'rand_3bp_frm_index', 'target_NGS', 'full_NGS']
            util.make_excel(WORK_DIR + "output/wo3bp/" + str(fn_nm) + "_result_wo3bp_" + FASTQ_N[d0_d4_idx], head, data_list, 2)

            sorted_err_list = logic_prep.sort_list_by_ele(err_list, 0)
            util.make_tsv(WORK_DIR + "output/wo3bp/" + str(fn_nm) + "_err_wo3bp_" + FASTQ_N[d0_d4_idx], head, sorted_err_list)
            pool.close()


if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    multi_processing_wo_randbp_in_brcd()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))