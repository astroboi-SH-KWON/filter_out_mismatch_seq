import time

import Util
import Logic
import LogicPrep
############### start to set env ################
WORK_DIR = "D:/000_WORK/KimNahye/20200827/WORK_DIR/"
# WORK_DIR = os.getcwd() + "/"
PROJECT_NAME = WORK_DIR.split("/")[-2]
FASTQ = "FASTQ/"
INPUT = "input/"
GUIDE_BARCODE_CSV = "190509_FINAL.CSV"

SCAFFOLD_SEQ = "GTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT"
FRONT_SCAF = SCAFFOLD_SEQ[:5]
FRONT_SCAF_WIN = 2
FRONT_SCAF_POS = 50

############### end setting env #################


def main():
    util = Util.Utils()
    csv_list = util.csv_to_list_ignr_nLine_header(WORK_DIR + INPUT + GUIDE_BARCODE_CSV)
    guide_list = [tmp_arr[2] for tmp_arr in csv_list]
    barcd_list = [tmp_arr[6] + tmp_arr[7] for tmp_arr in csv_list]
    trgt_list = [tmp_arr[8] for tmp_arr in csv_list]
    d0_seq_wo_scaf_list = [tmp_arr[3] + tmp_arr[4] + tmp_arr[5] for tmp_arr in csv_list]























if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    main()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))