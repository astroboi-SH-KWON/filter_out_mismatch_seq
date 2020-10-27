import glob
from Bio import SeqIO
import openpyxl
import os

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
        self.make_row(sheet, row, header[strt_idx:])

        for data_arr in data_list:
            row += 1
            self.make_row(sheet, row, data_arr[strt_idx:])

        workbook.save(filename=path + self.ext_xlsx)

    def merge_multi_list(self, pool_list):
        data_list = []
        err_list = []
        for tuple_val in pool_list:
            data_list.extend(tuple_val[0])
            err_list.extend(tuple_val[1])
        return data_list, err_list

    def make_tsv(self, path, header, data_list, strt_idx=0, deli='\t'):
        with open(path + self.ext_txt, 'w') as f:
            tmp_head = ''
            for head in header[strt_idx:]:
                tmp_head += (head + deli)
            f.write(tmp_head[:-1] + "\n")

            for data_arr in data_list:
                tmp_row = ''
                for row_val in data_arr[strt_idx:]:
                    tmp_row += (row_val + deli)
                f.write(tmp_row[:-1] + "\n")

    def merge_multi_err_list(self, pool_list):
        data_list = []
        err_list = []

        # make init of err_list by the size of pool_list[0][1]
        for err_idx in range(len(pool_list[0][1])):
            err_list.append([])

        for tuple_val in pool_list:
            data_list.extend(tuple_val[0])
            for err_idx in range(len(tuple_val[1])):
                err_list[err_idx].extend(tuple_val[1][err_idx])
        return data_list, err_list

    def split_big_file_by_row(self, init):
        big_file_path = init['big_file_path']
        num_row = init['num_row']
        splited_files_dir = init['splited_files_dir']
        output_file_nm = init['output_file_nm']
        output_file_ext = init['output_file_ext']

        os.makedirs(splited_files_dir, exist_ok=True)

        with open(big_file_path) as fin:
            fout = open('{}{}_{}{}'.format(splited_files_dir, output_file_nm, '0', output_file_ext), "w")
            for i, line in enumerate(fin):
                fout.write(line)
                if (i + 1) % num_row == 0:
                    fout.close()
                    fout = open('{}{}_{}{}'.format(splited_files_dir, output_file_nm, str(i // num_row + 1), output_file_ext), "w")

            fout.close()


