################################################################################
## Function:
## This file reads in one lociParis_rs_reformated.txt file
## and then output locipair.set and locipair.snps files
################################################################################
import datetime 
# from progress.bar import Bar

#  pheno = 'CHD'
#  pathout = "/disks/raid5a/szhang/rrintcc/Data_BOOST_Results/"+pheno+"/all_sets/"
#  fname = '/disks/raid5a/szhang/rrintcc/Data_BOOST_Results/'+pheno+'/lociPairs_rs_reformated.txt'

pathout = 'all_sets'
fname = 'lociPairs_rs_reformated.txt'

print "----------------------------------------------------------"
print "Output path: "+pathout
print "Locipairs file: "+ fname
# print "----------------------------------------------------------"
# print "Getting data size.."
# fpair = open(fname, mode = 'r')
# numpair = 0
# for line in fpair:
#     numpair = numpair+1
# fpair.close()
# steplen = numpair/100

print "----------------------------------------------------------"
print "Separating the sets..."
st = datetime.datetime.now()

fpair = open(fname, mode = 'r')
i = 0 # index (starts from 1) of loci pairs 
# bar = Bar('Processing', fill='=', suffix = '%(percent)d%%')
for line in fpair:
    i = i+1
    # if i%steplen == 0: bar.next()
    tmp = line.split('|')
    snpset1 = tmp[1].split()
    snpset2 = tmp[2].split()

    fileSet = open(pathout+'locipair'+str(i)+'.set', mode = 'w')
    fileSNPs = open(pathout+'locipair'+str(i)+'.snps', mode = 'w')

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

fpair.close()

# bar.finish()
et = datetime.datetime.now()
print 'Time for separating the sets: ' + str((et - st).seconds) + 's'




