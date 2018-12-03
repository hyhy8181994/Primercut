import argparse
import operator
from itertools import product
import os

# parsing the arguments and warning message
ap = argparse.ArgumentParser(description="Primer trimming tools")
ap.add_argument("-in_file",required=True,
	help="The name of the input FASTA file")
ap.add_argument("-out_file",required=True,
	help="The name of the trimmed reads file")
ap.add_argument("-unk_file",required=True,
	help="The name of the file containing unprocessed reads")
ap.add_argument("-n_mismatch",type=int ,required=True,
	help="The tolerance for mismatches in the forward or reverse primer sequence")
ap.add_argument("-min_len",type=int, required=True,
	help="The minimum length a sequence must be in order for it to be processed")
ap.add_argument("-forward", required=True,
	help="Forward primer")
ap.add_argument("-reverse", required=True,
	help="Reverse primer")

args = vars(ap.parse_args())
in_file=args["in_file"]
out_file=args["out_file"]
unk=args["unk_file"]
mismatch=args["n_mismatch"]
min_len=args["min_len"]
primer1=args["forward"]
primer2=args["reverse"]


with open(in_file,"r") as file:

    header = list()
    linelist = list()
    for line in file.readlines():

        if ">" in line:

           header.append(line.strip("\n"))
        else:
           linelist.append(line.strip("\n"))
           linelist = list(filter(None,linelist))

    dic_lines = dict(zip(header, linelist))
# construct dictionary of sequences indexed by headers
log_file = "logfile.log"


ambiguous_bases = {
    "A": ["A"],
    "C": ["C"],
    "T": ["T"],
    "G": ["G"],
    "N": ["A","T","C","G"],
    "M": ["A","C"]
}

base_translation = {
    "A": ["T"],
    "T": ["A"],
    "G": ["C"],
    "C": ["G"],
    "N": ["A","T","C","G"],
    "M": ["T","G"]
}


# invert primer2
re_primer2=primer2[::-1]
#translate the ambigious base
primer1_tr = list(map("".join, product(*map(ambiguous_bases.get,primer1))))
re_primer2_tr = list(map("".join, product(*map(base_translation.get,re_primer2))))
# caluate hamming distance(str1 has to be primer)
def hamming(str1, str2):
    b = list()
    count = 1
    if len(str1) > 1:
        for i in str1:
            assert len(i) == len(str2)
            ne = operator.ne
            a = sum(map(ne, ''.join(i), str2))
            b = b + [a]

            if len(b) > 1:
                dis = min(b)
            if count != len(str1):
                count += 1
                continue
            else:
                return dis
    else:
        assert len(''.join(str1)) == len(str2)
        ne = operator.ne
        return sum(map(ne, ''.join(str1), str2))

#create output file

output=open(out_file,"w+")
unfile=open(unk,"w+")
log = open(log_file, "w+")
log.write("Log File\r\n")
tri_count = 0
raw_length = []
tri_length = []
log1 = open("log.temp.log", "w+")
log1.write("Log")
# create a loop to scan the head and end of sequences
for key, val in dic_lines.items():
    cc1 = 0
    cc2 = 0
    len1 = len(primer1)
    len2 = len(re_primer2)
    dd1 = cc1 + len1
    dd2 = cc2 - len2
    rest = list()
    r_dis = list()
    r_dis2 = list()
    r_cor1 = list()
    r_cor2 = list()
    re_count = 0
    raw_length = raw_length + [len(val)]

    if len(val) < min_len:
        unfile.write("%s\r\n" % (key))
        unfile.write("%s\r\n" % (val))
        log1 = open("log.temp.log", "a+")
        log1.write("%s\r\n" % (key))
        log1.write("The sequence is less than %d\r\n" % (min_len)) # add %n
        continue
    else:
        for i in range(50):  # increase range if primer isn't found
            p1 = val[cc1:dd1]
            dis = hamming(primer1_tr, p1)
            if dis > mismatch:
                while dd1 == len1 + 49:
                    print(key)
                    unfile.write("%s\r\n" % (key))
                    unfile.write("%s\r\n" % (val))
                    log1 = open("log.temp.log", "a+")
                    log1.write("%s\r\n" % (key))
                    log1.write("The forward primer is not detected in this sequece\r\n")
                    break
                cc1 += 1
                dd1 += 1
                continue
            else:
                cor1 = (cc1, dd1)
                r_dis = r_dis + list(str(dis))
                r_cor1 = r_cor1 + list(str(cor1))
                if len(r_dis) > 3:  # need further inspection
                    dis_dir = dict(zip(r_dis, r_cor1))
                    min_l = min(r_dis)
                    cor1 = dis_dir[min_l]
                else:
                    p11 = val[dd1::]

            if re_count == 0:
                p2 = p11[dd2::]
                dis2 = hamming(re_primer2_tr, p2)

                if dis2 > mismatch:
                    re_count += 1
                    cc2 -= 1
                    dd2 -= 1
                else:
                    dis2_1 = dis2
                    cor2_1 = dd2
                    r_dis2 = r_dis2 + [dis2_1]
                    r_cor2 = r_cor2 + [cor2_1]

            if re_count > 0:
                p2 = p11[dd2:cc2]
                dis2 = hamming(re_primer2_tr, p2)

                if dis2 > mismatch:
                    cc2 -= 1
                    dd2 -= 1
                    re_count += 1
                    while re_count > 50:  # also remember to change this
                        unfile.write("%s\r\n" % (key))
                        unfile.write("%s\r\n" % (val))
                        log1 = open("log.temp.log", "a+")
                        log1.write("%s\r\n" % (key))
                        log1.write("The reverse primer is not detected in this sequece\r\n")
                        break
                    continue
                else:
                    cor2 = dd2
                    if cor2 is None:
                        pass
                    else:
                        dis2_2 = dis2
                        cor2_2 = dd2
                        r_dis2 = r_dis2 + [dis2_2]
                        r_cor2 = r_cor2 + [cor2_2]
            dis_dir2 = dict(zip(r_dis2, r_cor2))

            if len(r_dis2) >= 2:
                # if have time add staff if more than 1 match is dectected
                pass
                # min_l2 = min(r_dis2)
                # cor2 = r_dis2[min_l2]
            else:
                cor_reomve = dis_dir2[r_dis2[0]]
                if type(cor_reomve) is int:
                    tri_count += 1
                    p22 = p11[0:cor_reomve]
                    print("Trimming %s" % (key))
                    output.writelines("%s\r\n" % (key))
                    output.writelines("%s\r\n" % (p22))
                    tri_length = tri_length + [len(p22)]
                    break
                else:
                    print("There is something wrong!")

            cc1 += 1
            dd1 += 1
            cc2 -= 1
            dd2 -= 1
            re_count += 1

pro_reads = len(dic_lines.items())
raw_av = sum(raw_length) / float(len(raw_length))
tri_av = sum(tri_length) / float(len(tri_length))
log1 = open("log.temp.log","r")
temp_reads = log1.readlines()
log = open(log_file,"w")

log.write("Processed Reads: %d\r\n" % (pro_reads))

log.write("Trimmed Reads: %d\r\n" % (tri_count))

log.write("Average Read Length (Raw): %d\r\n" % (raw_av))

log.write("Average Read Length (Trimmed): %d\r\n" % (tri_av))

log.write("Unprovcessed Read Reasons:\r\n")

for line in temp_reads:
    log.write(line)



os.remove("log.temp.log")
log.close()
output.close()
unfile.close()