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


FASTQ = "FASTQ/20201019/"
OUTPUT = "output/"

LibA_D0 = ['LibA_D0']
WT = ['WT']
D0_D4_FLAG_ARR = [True, False]  # Day 0 : True, non Day 0 : False
FASTQ_ARR = [LibA_D0, WT]
FASTQ_N = ['LibA_D0', 'WT']
FASTQ_EXT = ".extendedFrags.fastq"

############### end setting env #################

def make_filtered_FASTQ():
    logic_prep = LogicPrep.LogicPreps()
    util = Util.Utils()
    for d0_d4_idx in range(len(D0_D4_FLAG_ARR)):
        for fn_nm in FASTQ_ARR[d0_d4_idx]:
            fastq_path = WORK_DIR + FASTQ + str(fn_nm) + FASTQ_EXT
            filtered_result_path = WORK_DIR + "input/" + str(fn_nm) + "_result_" + FASTQ_N[d0_d4_idx] + ".CSV"

            full_result_list = util.read_csv_ignore_N_line(filtered_result_path)
            ngs_read_list = logic_prep.make_1_arr_list_to_list(7, full_result_list)
            full_result_list.clear()
            ngs_read_set = set(ngs_read_list)
            ngs_read_list.clear()

            with open(fastq_path, 'r') as fastq_f:
                with open(fastq_path + '.result', 'w') as result_f:
                    while True:
                        fastq_id = fastq_f.readline()
                        if fastq_id == '':
                            break
                        fastq_seq = fastq_f.readline()
                        fastq_flag = fastq_f.readline()
                        fastq_qual = fastq_f.readline()
                        if fastq_seq[:-1] in ngs_read_set:
                            result_f.write(fastq_id + fastq_seq + fastq_flag + fastq_qual)


if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    make_filtered_FASTQ()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))