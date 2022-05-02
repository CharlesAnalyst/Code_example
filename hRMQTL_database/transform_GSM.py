import os
import xlrd
import pandas as pd
from multiprocessing import Pool

def get_GSM_list():
    excel = "/home/galaxy/Table_S1.xlsx"
    wb = xlrd.open_workbook(excel)
    sheet = wb.sheet_by_index(0)
    ip_list, input_list = sheet.col_values(4), sheet.col_values(5)
    ip_GSM = [str(x).split(",") for x in ip_list if str(x).strip("'").startswith("GSM")]
    input_GSM = [str(x).split(",") for x in input_list if str(x).strip("'").startswith("GSM")]
    #
    def merge_list(in_list):
        sigle_GSM = []
        for x in in_list:
            sigle_GSM += x
        return sigle_GSM
    #
    ip_result = merge_list(ip_GSM)
    input_result = merge_list(input_GSM)
    return ip_result, input_result

def query_SRR_according_GSM(GSM):
    command="esearch -db sra -query %s |efetch -format docsum |xtract -pattern DocumentSummary -element Run@acc" % GSM
    #  return str(os.popen(command).read())
    print("\t%s\t%s" % (GSM, str(os.popen(command).read())))
    # fastq-dump -A $SRR -O sra_fastq

ip_result, input_result = get_GSM_list()
ip_srr_list = ["%s\t%s\n" % (x, query_SRR_according_GSM(x)) for x in ip_result]
input_srr_list = ["%s\t%s\n" % (x, query_SRR_according_GSM(x)) for x in input_result]
result_dir = "/home/galaxy/project/m6AQTL/data/RIP_seq"
with open(os.path.join(result_dir, "ip_id.txt"), 'w') as fw_ip:
    fw_ip.writelines(ip_srr_list)
with open(os.path.join(result_dir, "input_id.txt"), 'w') as fw_input:
    fw_input.writelines(input_srr_list)

#####
in_file = "/home/galaxy/project/m6AQTL/data/RIP_seq/m6avar_GSM.txt"
out_file = "/home/galaxy/project/m6AQTL/data/RIP_seq/m6avar_SRA.txt"
in_list = pd.read_table(in_file, sep="\t", header=None).iloc[:, 0]
pool = Pool()
for x in in_list:
	print("\t%s\t%s" % (x, query_SRR_according_GSM(x).strip()))
	pool.apply_async(query_SRR_according_GSM, (x,))
pool.close()
pool.join()