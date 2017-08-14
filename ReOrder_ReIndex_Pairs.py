###############################################################################
## Function:
## This file reads in original BOOST result file InteractionRecords.txt
## and the plink map file ESRD2_..._Ch1-22.map
## and generates the reformated BOOST results: results_ESRD2_BOOST.txt
################################################################################
from itertools import islice
from scipy.stats import chi2

ori_file = open(r"MIGen_InteractionRecords.txt",mode = 'r')
result_file = open(r"results_MIGen_BOOST.txt", mode = 'w+')
ori_file_map = open(r"MIGen_QC.map", mode = 'r')

index = 0

SNPs = []
BPs = []
CHRs = []

for line in ori_file_map:
    statistics = line.split()
    CHRs.append(statistics[0])
    SNPs.append(statistics[1])
    BPs.append(statistics[3])

# The indices of SNPs in the three (first, last, remain) result files start from 0
all_list = []
while True:
    next_n_lines = list(islice(ori_file,20))
    if not next_n_lines:
        break
    res_n_lines=[]    
    for line in next_n_lines:
        if len(line)==0:
            continue
        sample = line.split()

        # IND1  IND2   InteractionBOOST, SNP1, SNP2, BP1, BP2, InteractionPVal, CHR1, CHR2
        # For GBOOST results: use float(sample[5])
        # For GBOOST2 results: use float(sample[3])
        chi4GBOOST = sample[5]
        pGBOOST = 1 - chi2.cdf(float(chi4GBOOST), 4) 
        all_list.append([sample[1], sample[2], chi4GBOOST,
                         SNPs[int(sample[1])], SNPs[int(sample[2])],
                         BPs[int(sample[1])], BPs[int(sample[2])], str(pGBOOST), CHRs[int(sample[1])], CHRs[int(sample[2])]])

all_list_sorted = sorted(all_list, key=lambda tup: tup[2], reverse = True)

for ele in all_list_sorted:
    result_file.write(' '.join(ele)+'\n')

ori_file.close()
ori_file_map.close()
result_file.close()
