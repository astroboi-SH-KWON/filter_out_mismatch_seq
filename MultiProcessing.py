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
FASTQ = "FASTQ/201026_PAM_variant_LibraryA/"
INPUT = "input/"
# GUIDE_BARCODE_CSV = "190509_FINAL.CSV"
GUIDE_BARCODE_CSV = "1st_LibraryA.CSV"

D0_LIST = ['LibA_D0_1'
           , 'LibA_D0_2'
           ]
NON_D0_LIST = ['NRCH_previous'
                , 'NRCH_re_Rep1'
                , 'NRCH_re_Rep2'
                , 'NRRH_Rep1'
                , 'NRRH_Rep2'
                , 'NRTH_Rep1'
                , 'NRTH_Rep2'
                , 'sc++_Rep1'
                , 'sc++_Rep2'
                , 'SpCas9_NG_Rep1'
                , 'SpCas9_NG_Rep2'
                , 'SpCas9_Rep1'
                , 'SpCas9_Rep2'
                , 'SpG_Rep1'
                , 'SpG_Rep2'
                , 'SpRY_Rep1'
                , 'SpRY_Rep2'
                , 'VRQR_Rep1'
                , 'VRQR_Rep2'
                , 'WT_with_G_RY_Rep1'
                , 'WT_with_G_RY_Rep2'
                , 'xCas9_Rep1'
                , 'xCas9_Rep2'
               ]
# D0_D4_FLAG_ARR = [True, False]  # Day 0 : True, non Day 0 : False
D0_D4_FLAG_ARR = [False]  # Day 0 : True, non Day 0 : False 20201028 request to run D0 as non_D0
FASTQ_ARR = [D0_LIST, NON_D0_LIST]
# FASTQ_N = ['LibA_D0', 'LibA_non_D0']
FASTQ_N = ['LibA_D0_without_filtering_trgt', 'LibA_non_D0']
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
            fastq_list.clear()

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
            init_split_file = {'big_file_path': WORK_DIR + FASTQ + str(fn_nm) + FASTQ_EXT
                                , 'num_row': 4000000
                                , 'splited_files_dir': WORK_DIR + FASTQ + str(fn_nm) + "/"
                                , 'output_file_nm': str(fn_nm)
                                , 'output_file_ext': FASTQ_EXT
                               }
            util.split_big_file_by_row(init_split_file)
            print("end to split : ", str(fn_nm))

            splited_fastq_arr = util.get_files_from_dir(WORK_DIR + FASTQ + str(fn_nm) + '/' + '*' + FASTQ_EXT)

            for fq_idx in range(len(splited_fastq_arr)):
                path_arr = splited_fastq_arr[fq_idx].split('/')
                print(path_arr[-2:])
                fastq_list = util.read_fastq_to_list(splited_fastq_arr[fq_idx])

                splited_fastq_list = np.array_split(fastq_list, MULTI_CNT)
                fastq_list.clear()

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
                os.makedirs(WORK_DIR + FASTQ + str(fn_nm) + "/output/", exist_ok=True)
                os.makedirs(WORK_DIR + FASTQ + str(fn_nm) + "/input/", exist_ok=True)
                fl_nm = path_arr[-1].replace(FASTQ_EXT, "")
                util.make_excel(
                    WORK_DIR + FASTQ + str(fn_nm) + "/output/" + fl_nm + "_result_" + FASTQ_N[d0_d4_idx], head,
                    data_list, 2)
                util.make_tsv(
                    WORK_DIR + FASTQ + str(fn_nm) + "/input/" + fl_nm + "_result_" + FASTQ_N[d0_d4_idx], head,
                    data_list, 2, ",")

                for err_arr in err_list:
                    if len(err_arr) == 0:
                        continue

                    try:
                        util.make_excel(
                            WORK_DIR + FASTQ + str(fn_nm) + "/output/" + fl_nm + "_err_" + err_arr[0][0] + "_" +
                            FASTQ_N[d0_d4_idx], head, err_arr)
                    except Exception as err:
                        util.make_tsv(
                            WORK_DIR + FASTQ + str(fn_nm) + "/output/" + fl_nm + "_err_" + err_arr[0][0] + "_" +
                            FASTQ_N[d0_d4_idx], head, err_arr)


def split_file_for_plan_B():
    util = Util.Utils()
    for d0_d4_idx in range(len(D0_D4_FLAG_ARR)):
        for fn_nm in FASTQ_ARR[d0_d4_idx]:
            print("start to split : ", str(fn_nm))
            init_split_file = {'big_file_path': WORK_DIR + FASTQ + str(fn_nm) + FASTQ_EXT
                                , 'num_row': 4000000
                                , 'splited_files_dir': WORK_DIR + FASTQ + str(fn_nm) + "/"
                                , 'output_file_nm': str(fn_nm)
                                , 'output_file_ext': FASTQ_EXT
                               }
            util.split_big_file_by_row(init_split_file)
            print("end to split : ", str(fn_nm))


if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    # multi_processing()

    # split_file_for_plan_B()
    multi_processing_plan_B()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))