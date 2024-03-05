import os, argparse, time, click

# test
# -m /8_reads_asv_map_convert.txt -t 0.005 -o 240119_o
# -m /8_reads_asv_map_convert.txt -t 0.005 -o 240207
# no F.G., Mexico
# 240207: 增加 asv blast结果筛选 用 9_asv_all-ex_py.txt, thomas 的 py脚本生成的文件 /9_asv_all-ex_py.txt
# -m /8_reads_asv_map_convert.txt -t 0.005 -o 240210 -b /9_asv_all-ex_py.txt -l family -n Staphylinidae

parser = argparse.ArgumentParser(description='Find Top ASV/Sample/OTU')
parser.add_argument('-m', '--matrix', dest='matrix', required=True, type=str)
parser.add_argument('-b', '--blast', dest='blast', required=False, type=str)
parser.add_argument('-l', '--level', dest='level', required=False, type=str)
parser.add_argument('-n', '--name', dest='level_name', required=False, type=str)
parser.add_argument('-t', '--threshold', dest='threshold', required=True, type=float)
parser.add_argument('-o', '--outname', dest='outname', required=True, type=str)

args = parser.parse_args()
file = args.matrix
blast = args.blast
b_level = args.level
b_name = args.level_name
threshold = args.threshold
out = args.outname + "_asvSample.txt"


# out2 = args.outname + "_otu_asv.txt"


def readFile(file):
    dic = {}  # otu_sample_asv_read

    with open(file) as f:
        for line in f:
            line_t = line.strip().split('\t')
            if 'asv' not in line:
                otu, sa, asv, reads, quality, picName = line_t[2], line_t[1], line_t[0], line_t[3], line_t[6], line_t[7]
                sa_pic[sa] = picName
                if 'MSLM017' not in sa:  # no F.G., Mexico
                    if float(quality) >= threshold:
                        if otu not in dic:
                            dic[otu] = {sa: {asv: {reads: quality}}}
                        elif sa not in dic[otu]:
                            dic[otu][sa] = {asv: {reads: quality}}
                        elif asv not in dic[otu][sa]:
                            first_value = list(list(dic[otu][sa].values())[0].keys())[0]
                            if int(first_value) <= int(reads):
                                if int(first_value) == int(reads):
                                    dic[otu][sa][asv] = {reads: quality}
                                else:
                                    dic[otu][sa] = {asv: {reads: quality}}
    return dic


def writeFile(info):
    print("Print duplicate groups:\n")
    with open(out, 'w') as f_w:
        f_w.write('OTU\tSample\tASV\tReads\tSectorQuality\n')
        for outer_key, inner_dict in info.items():
            # print(f"Outer Key: {outer_key}")  # outer_key: otuxx; inner_dict: {'BIMB_MSLM006_A01': {'202211_uniq18085': {'19': '0.00996852'}}, 'BIMB_MSLM006_C01': {'202211_uniq2967': {'158': '0.091224018'}}, ...}
            for inner_key, value in inner_dict.items():  # inner_key: BIMB_MSLM006_A01; value: {'202211_uniq18085': {'19': '0.00996852'}}
                if len(value.items()) != 1:
                    print(f"{outer_key}\t{inner_key}\t{value}")
                for i in value:  # i: 202211_uniq18085
                    # print(list(value[i].keys())[0]) # 19
                    # print(list(value[i].values())[0]) # 0.00996852
                    if asv_blast_want:
                        if i in asv_blast_want:
                            f_w.write(
                                f"{outer_key}\t{inner_key}\t{i}\t{(list(value[i].keys())[0])}\t{list(value[i].values())[0]}\n")
                    else:
                        f_w.write(
                            f"{outer_key}\t{inner_key}\t{i}\t{(list(value[i].keys())[0])}\t{list(value[i].values())[0]}\n")


# def writeOTU(info):
#     with open(out2, 'w') as f_w2:
#         # f_w2.write('OTU\tASV\n')
#
#         for outer_key, inner_dict in info.items():
#             print(f"Outer Key: {outer_key}")
#             for inner_key, value in inner_dict.items():
#                 print(f"  Inner Key: {inner_key}, Value: {value}")
#                 print()


def readBlast(blastFile, taxa_level, want_taxa):
    targetCol = 0
    with open(blastFile) as b:
        for line in b:
            line_t = line.strip().split('\t')
            if 'taxid' in line:
                for colidx in range(1, len(line_t) - 1):
                    if line_t[colidx] == taxa_level:
                        # print(colidx)
                        targetCol = colidx
            else:
                asv = line_t[0].split(';')[0]
                if len(line_t) > targetCol:
                    target = line_t[targetCol]
                    if target == want_taxa:
                        asv_blast_want.append(asv)
                        # print(asv)


if __name__ == '__main__':
    sa_pic = {}
    asv_blast_want = []
    write2 = {}

    info = readFile(file)
    if blast:
        readBlast(blast, b_level, b_name)
        # asv_blast_want.sort()
    writeFile(info)
    # writeOTU(info)
    print("Done")
    print()
