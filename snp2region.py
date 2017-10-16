# usage: py snp2region.py 200 5e-10
# 200 means 200kbp
# 5e-10 means the p value threshold

import pdb, progressbar, sys
from itertools import islice
nslice = 20 

if not len(sys.argv) == 3:
	print("err: please specify the range, e.g. py snp2region 200 5e-10")
rangename = "lociRange_" + sys.argv[1] + "kbp_"+sys.argv[2]+".txt"
mapname  = "MIGen_QC.map"
resultname = "MIGen_range_snps_" + sys.argv[1] + "kbp_"+sys.argv[2]+".txt"
resultname_rs = "MIGen_range_snps_rs_" + sys.argv[1] + "kbp_"+sys.argv[2]+".txt"

# rangefile: lociRange.txt: CHR STARTBP ENDBP

def isinrange(rangechr, snpchr, rangestart, rangeend, snpbp):
    if snpchr < rangechr:
        return -1
    elif snpchr == rangechr:
        if snpbp < rangestart:
            return -1
        elif rangestart <= snpbp <= rangeend:
            return 0
        elif snpbp > rangeend:
            return 1
        else:
            print("err1")
    elif snpchr > rangechr:
        return 1
    else:
        print("err2")

rangefile = open(rangename, mode = "r")
mapfile  = open(mapname, mode = "r")
outfile = open(resultname, mode = "w")
outfile2 = open(resultname_rs, mode = 'w')

rangelist = [] # [[CHR startBP endBP BP1 BP2 ...]]
rangelist_rs = [] # [[CHR startBP endBP rs1 rs2 ...]]
snplist  = [] # [[CHR rsSNP BP]]

while True:
    next_n_lines = list(islice(rangefile, nslice))
    if not next_n_lines:
        break

    for line in next_n_lines:
        if len(line) == 0:
            continue
        rangelist.append(line.split())
        rangelist_rs.append(line.split())

rangefile.close()
print("Finished reading the range file: " + rangename + '\n')

while True:
    next_n_lines = list(islice(mapfile, nslice))
    if not next_n_lines:
        break

    for line in next_n_lines:
        if len(line) == 0:
            continue
        snplist.append([line.split()[i] for i in [0, 1, 3]])        

mapfile.close()
print("Finished reading .map file" + mapname + '\n')

# i is the index for rangelist (e.g. 5000+)
# j is the index for snplist  (e.g. 650000+) 
i = 0
j = 0
numsnp = len(snplist)
numrange = len(rangelist)


## Notes: Cannot use the logic used in snp2gene.py
## Because now the regions are not ordered 

bar = progressbar.ProgressBar(redirect_stdout = True, max_value = len(snplist))
for j in range(len(snplist)):
    bar.update(j)
    snpchr = int(snplist[j][0])
    snpbp  = int(snplist[j][2])

    for i in range(len(rangelist)):
        bpchr   = int(rangelist[i][0])
        bpstart = int(rangelist[i][1])
        bpend   = int(rangelist[i][2])
        
        if isinrange(bpchr, snpchr, bpstart, bpend, snpbp) == 0:
            rangelist[i].append(snplist[j][2])
            rangelist_rs[i].append(snplist[j][1])
            break
bar.finish()


print("Finished the snp -> range mapping\n")

for region in rangelist:
    outfile.write('\t'.join(region[0:3]) + '\t' + str(len(region)-3)+'\t' + '\t'.join(region[3:]) + '\n')

for region in rangelist_rs:
    outfile2.write('\t'.join(region[0:3]) + '\t' + str(len(region)-3)+'\t' + '\t'.join(region[3:]) + '\n')

outfile.close()
outfile2.close()
print("Finished writing the result file: " + resultname + '\n')
print("Finished writing the result file: " + resultname_rs + '\n')


























