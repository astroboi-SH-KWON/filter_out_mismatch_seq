import time

import Util
import Logic
import LogicPrep
############### start to set env ################
WORK_DIR = "D:/000_WORK/JangHyeWon_KimMinYung/20200703/WORK_DIR/"
# WORK_DIR = os.getcwd() + "/"
PROJECT_NAME = WORK_DIR.split("/")[-2]
INPUT = "input/"
GUIDE_BARCODE_EXCEL = "190509_FINAL.xlsx"
############### end setting env #################





















if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")

    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))