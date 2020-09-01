import glob
from Bio import SeqIO
import openpyxl

import Logic

class Utils:
    def __init__(self):
        self.ext_txt = ".txt"
        self.ext_dat = ".dat"
        self.ext_xlsx = ".xlsx"

    """
    get file lists in target dir by target ext
    :param
        path : target dir + "*." + target ext
    :return
        ['target dir/file_name.target ext', 'target dir/file_name.target ext' ...]
    """
    def get_files_from_dir(self, path):
        return glob.glob(path)

    def read_csv_ignore_N_line(self, path, deli_str=",", n_line=1):
        result_list = []
        with open(path, "r") as f:
            for ignr_line in range(n_line):
                header = f.readline()
                print(header)
            while True:
                tmp_line = f.readline().replace("\n", "")
                if tmp_line == '':
                    break

                result_list.append(tmp_line.split(deli_str))
        return result_list

    def read_fastq_to_list(self, path):
        temp = list(SeqIO.parse(path, "fastq"))
        return [str(temp[k].seq) for k in range(len(temp))]

    def make_row(self, sheet, row, data_arr, col=1):
        for idx in range(len(data_arr)):
            sheet.cell(row=row, column=(col + idx), value=data_arr[idx])

    def make_excel(self, path, header, data_list, strt_idx=0):
        workbook = openpyxl.Workbook()
        sheet = workbook.active

        row = 1
        self.make_row(sheet, row, header)

        for data_arr in data_list:
            row += 1
            self.make_row(sheet, row, data_arr[strt_idx:])

        workbook.save(filename=path + self.ext_xlsx)

    def merge_multi_list(self, pool_list):
        data_list = []
        err_list = []
        for tuple_val in pool_list:
            print("str(len(tuple_val[0])) : ", str(len(tuple_val[0])))
            print("str(len(tuple_val[1])) : ", str(len(tuple_val[1])))
            data_list.extend(tuple_val[0])
            err_list.extend(tuple_val[1])
        return data_list, err_list