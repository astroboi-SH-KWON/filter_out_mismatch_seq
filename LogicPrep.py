
class LogicPreps:

    def __init__(self):
        self.ext_fa = ".fa"
        self.ext_dat = ".dat"
        self.ext_gtf = ".gtf"

    def make_1_arr_list_to_list(self, idx_1, arr_list):
        return [tmp_arr[idx_1] for tmp_arr in arr_list]

    def make_2_arr_list_to_list(self, idx_1, idx_2, arr_list):
        return [tmp_arr[idx_1] + tmp_arr[idx_2] for tmp_arr in arr_list]

    def make_3_arr_list_to_list(self, idx_1, idx_2, idx_3, arr_list):
        return [tmp_arr[idx_1] + tmp_arr[idx_2] + tmp_arr[idx_3] for tmp_arr in arr_list]