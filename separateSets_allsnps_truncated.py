################################################################################
## Function:
## This file reads in range_snps.txt and lociPairs.txt file
## and then output locipair.set and locipair.snps files
################################################################################
import datetime, sys, os, progressbar, random

if not len(sys.argv) == 3:
	print("err: please specify the range and p_thr, e.g. py separateSets_allsnps_truncated.py 200 5e-10")

region_size = int(sys.argv[1])
p_thr = sys.argv[2]
max_snps 	= 31 # sqrt(1000)

pathout = 'all_sets_'+str(region_size) + 'kbp_truncated_'+p_thr
fname = 'MIGen_range_snps_rs_'+str(region_size)+'kbp_'+p_thr+'.txt'
pairname = 'lociPairs_'+str(region_size) + 'kbp_'+p_thr+'.txt'

if os.path.exists(pathout):
	os.system('rm -rf '+pathout+'/*.set')
	os.system('rm -rf '+pathout+'/*.snps')
else:
	os.system('mkdir '+pathout)

def truncateSet(all_snpset, snpset, max_snps):
	if len(all_snpset) <= max_snps:
		for snp in all_snpset:
			snpset.append(snp)
	else:
		randindices = random.sample(range(len(all_snpset)), max_snps)
		for k in randindices:
			snpset.append(all_snpset[k])



print("----------------------------------------------------------")
print("Output path: "+pathout)
print("Locipairs file: "+ pairname)
print("range_snps file: " + fname)

print("----------------------------------------------------------")
print("Separating the sets...")
st = datetime.datetime.now()

rangelist = []
with open(fname, mode = 'r') as f:
	for line in f:
		rangelist.append(line.split()[4:])


with open(pairname, mode = 'r') as f:
	lines = f.readlines()
	bar = progressbar.ProgressBar(max_value = len(lines))
	i = 0
	for line in lines:
		bar.update(i)
		pairind1 = int(line.split()[0])
		pairind2 = int(line.split()[1])

		all_snpset1 = rangelist[pairind1]
		all_snpset2 = rangelist[pairind2]

		snpset1 = []
		snpset2 = []

		truncateSet(all_snpset1, snpset1, max_snps)
		truncateSet(all_snpset2, snpset2, max_snps)

		fileSet = open(pathout+'/locipair'+str(i)+'.set', mode = 'w')
		fileSNPs = open(pathout+'/locipair'+str(i)+'.snps', mode = 'w')

		fileSet.write('SET_1'+'\n')
		for snp1 in snpset1:
		    fileSet.write(snp1+'\n')
		    fileSNPs.write(snp1+'\n')

		fileSet.write('END'+'\n\n')
		fileSet.write('SET_2'+'\n')

		for snp2 in snpset2:
		    fileSet.write(snp2+'\n')
		    fileSNPs.write(snp2+'\n')

		fileSet.write('END'+'\n')

		fileSet.close()
		fileSNPs.close()
		i += 1
	bar.finish()

et = datetime.datetime.now()
print('Time for separating the sets: ' + str((et - st).seconds) + 's')




