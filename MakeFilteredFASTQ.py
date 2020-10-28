import time
import os
import platform

import Util
import LogicPrep
############### start to set env ################
SYSTEM_NM = platform.system()
if SYSTEM_NM == 'Linux':
    # REAL
    WORK_DIR = os.getcwd() + "/"
else:
    # DEV
    WORK_DIR = "D:/000_WORK/KimNahye/20201019/WORK_DIR/"
PROJECT_NAME = WORK_DIR.split("/")[-2]

FASTQ = "FASTQ/201026_PAM_variant_LibraryA/"

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
D0_D4_FLAG_ARR = [True, False]  # Day 0 : True, non Day 0 : False
FASTQ_ARR = [D0_LIST, NON_D0_LIST]
FASTQ_N = ['LibA_D0', 'LibA_non_D0']
FASTQ_EXT = ".extendedFrags.fastq"

############### end setting env #################

def make_filtered_FASTQ():
    logic_prep = LogicPrep.LogicPreps()
    util = Util.Utils()
    for d0_d4_idx in range(len(D0_D4_FLAG_ARR)):
        for fn_nm in FASTQ_ARR[d0_d4_idx]:
            splited_fastq_dir = WORK_DIR + FASTQ + str(fn_nm) + "/"

            splited_fastq_arr = util.get_files_from_dir(splited_fastq_dir + '*' + FASTQ_EXT)

            merged_fastq = FASTQ + "merged_" + str(fn_nm) + FASTQ_EXT
            try:
                os.remove(merged_fastq)
            except Exception as err:
                print(err, "os.remove(merged_fastq)", merged_fastq)

            for fq_idx in range(len(splited_fastq_arr)):
                fastq_path = splited_fastq_arr[fq_idx]
                path_arr = fastq_path.split('/')
                print(path_arr[-2:])
                fl_nm = path_arr[-1].replace(FASTQ_EXT, "_result_" + FASTQ_N[d0_d4_idx] + ".txt")

                filtered_result_path = splited_fastq_dir + "input/" + fl_nm

                full_result_list = util.read_csv_ignore_N_line(filtered_result_path)
                ngs_read_list = logic_prep.make_1_arr_list_to_list(7, full_result_list)
                full_result_list.clear()
                ngs_read_set = set(ngs_read_list)
                ngs_read_list.clear()

                with open(fastq_path, 'r') as fastq_f:
                    with open(merged_fastq, 'a') as result_f:
                        while True:
                            fastq_id = fastq_f.readline()
                            if fastq_id == '':
                                break
                            fastq_seq = fastq_f.readline()
                            fastq_flag = fastq_f.readline()
                            fastq_qual = fastq_f.readline()
                            if fastq_seq[:-1] in ngs_read_set:
                                result_f.write(fastq_id + fastq_seq + fastq_flag + fastq_qual)

            print('Done :', merged_fastq, '👍\n\n')


if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>🤞")
    make_filtered_FASTQ()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))