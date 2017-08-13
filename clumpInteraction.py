############################################################################################
## Function: 
## This script read in a BOOST result file, do clumping and then output the following 6 files
## out_file1 = open('loci.txt', mode = 'w')
## out_file2 = open('loci_rs.txt', mode = 'w')
## out_file3 = open('lociPairs.txt', mode = 'w')
## out_file4 = open('lociPairs_rs.txt', mode = 'w')
## out_file5 = open('lociPairs_rs_reformated.txt', mode = 'w')
## out_file6 = open('lociRange.txt', mode = 'w')
############################################################################################

from itertools import islice
from scipy.stats import chi2
import gzip, datetime

############################################################################################
## Need to specify the following parameters before use!
############################################################################################

# Thresholds: 5E-8/39.69705, 5E-9/44.52118, 5E-10/49.3228, 5E-11/54.10619,
# 5E-12/58.8745, 5E-13/63.62993, 5E-14/68.37653, BC/67.88988
# BC (Bonferroni Correction) p-value threshold: 1.074603E-12/62.05125 (p=305054)

p_thr = 5e-8
chisq_thr = chi2.ppf(1-p_thr, 4)
region_size = 500 # unit: kbp
buffer_size = 25  # unit: kbp (upward_buffer = downward_buffer = buffer_size)
path = "D:\\Dropbox\\Research_Projects\\Region_Interaction_cc\\CUHK_Data_Processing\\Data\\"
fName = "results_ESRD2_BOOST.txt"

############################################################################################
print '-----------------------------------------------------------'
print "Chi-square threshold (df = 4): "+str(chisq_thr)
print "Corresponding p-value threshold: "+str(p_thr)
print "Region size: "+str(region_size)+" kbp"
print "Buffer size: "+str(buffer_size)+" kbp"


############################################################################################
# loci: [lowerBound, upperBound, [snp1,snp2,...(snps inside this clump)]]
# loci_rs: [snps_rs]
# lociPairs: [{loci1, loci2}, [interInfo1, interInfo2, ...]]
# interInfo: BP1 BP2 CHISQ   SNP1(rs)    SNP2(rs)

def snpClumping(loci, loci_rs, lociPairs, interInfo):
    # lociIndex start from 0, consistent with the index in the python list "loci"
    # In current func, lociIndex will have two components
	lociIndex = []  

	for indsnp in range(2): # SNP1 and SNP2 for this snp-snp interaction
		snp = interInfo[indsnp]
		snp_rs = interInfo[indsnp+3]
		snpRange = [ele for ele in loci if ele[0] <= snp <= ele[1]]
		if len(snpRange) == 0: # This snp does not exist in any identified loci
			lociIndex.append(len(loci))
            # Include the nearby region_size (e.g. 500) kbp region
			loci.append([snp - region_size*1000/2, snp + region_size*1000/2, [snp]])
			loci_rs.append([snp_rs])
		else:
			lociIndex.append(loci.index(snpRange[0]))
			if snp not in loci[lociIndex[-1]][2]: # if this snp is not recorded before
				loci_rs[lociIndex[-1]].append(snp_rs)
				loci[lociIndex[-1]][2].append(snp)
				loci[lociIndex[-1]][0] = min(loci[lociIndex[-1]][0], snp - buffer_size*1000)
				loci[lociIndex[-1]][1] = max(loci[lociIndex[-1]][1], snp + buffer_size*1000)

	if {lociIndex[0], lociIndex[1]} in [ele[0] for ele in lociPairs]:
		pairIndex = [ele[0] for ele in lociPairs].index({lociIndex[0], lociIndex[1]})
		if interInfo not in lociPairs[pairIndex][1]: 
			lociPairs[pairIndex][1].append(interInfo)
	else:
		lociPairs.append([{lociIndex[0], lociIndex[1]}, [interInfo]])


############################################################################################
## This part is used for clumping loci
# loci: [lowerBound, upperBound, [snp1,snp2,...(snps inside this clump)]]
# lociPairs: [{loci1, loci2}, [snpPairInfo1, snpPairInfo2, ...]]

## Format of results_CHD_BOOST.txt
# IND1	IND2	CHISQ	SNP1	SNP2	BP1	BP2
# snp1: IND1(order in .map file), SNP1, BP1 (Need to check the corresponding genome assembly)
# snp2: IND2 (order in .map file), SNP2, BP2 (Need to check the corresponding genome assembly)
startTime = datetime.datetime.now()
print '-----------------------------------------------------------'
print 'Clumping... '
print 'file: '+path+fName
fileInter = open(path + fName, mode='r')
all_list = []

while True:
	next_n_lines = list(islice(fileInter, 20))
	if not next_n_lines:
		break

	res_n_lines = []

	for line in next_n_lines:
		if len(line) == 0:
			continue
		sample = line.split()
		# SNP1(bp)	SNP2(bp)	InteractionBOOST	InteractionPLINK	   singlelocusAssoc1	singlelocusAssoc2
        # IND1	IND2	CHISQ	SNP1	SNP2	BP1	BP2
        # all_list/part_list: BP1 BP2 CHISQ   SNP1(rs)    SNP2(rs)
		all_list.append(
			[int(sample[5]), int(sample[6]), float(sample[2]), sample[3], sample[4]])

fileInter.close()
# Only select snp pairs that pass the chisq_thr
part_list = [ele for ele in all_list if ele[2] >= chisq_thr]
loci = []
loci_rs = []
lociPairs = []

for i in range(len(part_list)):
	# Need a function to judge whether a new snp belongs to an existing loci
	snpClumping(loci, loci_rs, lociPairs, part_list[i])

endTime = datetime.datetime.now()
print 'Time for clumping the snps is: ' + str((endTime - startTime).seconds) + 's'

############################################################################################
## Writing the results (sample and snpFinalIndex)
startTime = datetime.datetime.now()
print '-----------------------------------------------------------'
print 'Writing the results...'
out_file1 = open('loci.txt', mode = 'w')
out_file2 = open('loci_rs.txt', mode = 'w')
out_file3 = open('lociPairs.txt', mode = 'w')
out_file4 = open('lociPairs_rs.txt', mode = 'w')
out_file5 = open('lociPairs_rs_reformated.txt', mode = 'w')
out_file6 = open('lociRange.txt', mode = 'w')

# Format of 'loci.txt': each line contains all the snps of a locus
for ele in loci:
	out_file1.write(' '.join([str(x) for x in ele[2]])+'\n')
	out_file6.write(str(ele[0])+' '+str(ele[1])+'\n')

for ele in loci_rs:
	out_file2.write(' '.join([str(x) for x in ele])+'\n')

# Format of 'lociPairs.txt': lociInd1 lociInd2 snpPair1_1 snpPair1_2 snpPair2_1 snpPair2_2 snpPair3_1 snpPair3_2 ...
for ele in lociPairs:
	# ele[0]: loci indexes
	if len(ele[0])==1: 
		out_file3.write(str(list(ele[0])[0])+' '+str(list(ele[0])[0])+' ')
		out_file4.write(str(list(ele[0])[0])+' '+str(list(ele[0])[0])+' ')
	else: 
		out_file3.write(str(list(ele[0])[0])+' '+str(list(ele[0])[1])+' ')
		out_file4.write(str(list(ele[0])[0])+' '+str(list(ele[0])[1])+' ')

	# ele[1]: interInfo: BP1 BP2 CHISQ   SNP1(rs)    SNP2(rs)
	for snpPair in ele[1]:
		out_file3.write(' '.join([str(x) for x in snpPair[0:2]])+' ')
		out_file4.write(' '.join([str(x) for x in snpPair[3:]])+' ')
	out_file3.write('\n')
	out_file4.write('\n')


# Write the results of reformated locipair_rs
for ele in lociPairs:
	pairs_rs_unique = [[],[]]
	if len(ele[0])==1: 
		ind0 = list(ele[0])[0]
		ind1 = list(ele[0])[0]
	else: 
		ind0 = list(ele[0])[0]
		ind1 = list(ele[0])[1]

	for snpPair in ele[1]:
		for rssnp in snpPair[3:]:
			if (rssnp in loci_rs[ind0]) and (rssnp not in pairs_rs_unique[0]):
				pairs_rs_unique[0].append(rssnp)
			elif (rssnp in loci_rs[ind1]) and (rssnp not in pairs_rs_unique[1]):
				pairs_rs_unique[1].append(rssnp)

	out_file5.write(str(ind0)+' '+str(ind1)+' | ')
	out_file5.write(' '.join(pairs_rs_unique[0]) + ' | ')
	out_file5.write(' '.join(pairs_rs_unique[1]) + '\n')


out_file1.close()
out_file2.close()
out_file3.close()
out_file4.close()
out_file5.close()
out_file6.close()

endTime = datetime.datetime.now()
print 'Time for writing the results is: ' + str((endTime - startTime).seconds) + 's'
