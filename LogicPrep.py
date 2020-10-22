
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

    def get_real_scaf_seq(self, ngs_read_needle, strt_idx, end_idx):
        if end_idx == 0:
            return ngs_read_needle[strt_idx:].replace("-", "").strip()
        else:
            return ngs_read_needle[strt_idx: end_idx].replace("-", "").strip()

    def sort_list_by_ele(self, data_list, ele_idx, up_down_flag=True):
        result_list = []
        for tmp_arr in sorted(data_list, key=lambda tmp_arr: tmp_arr[ele_idx], reverse=up_down_flag):
            result_list.append(tmp_arr)
        return result_list

    def make_2_arr_list_to_list_after_slice(self, idx_1, idx_2, pos_slice, arr_list):
        return [tmp_arr[idx_1] + tmp_arr[idx_2][:pos_slice] for tmp_arr in arr_list]
