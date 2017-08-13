##############################################################################
import os, datetime
from progress.bar import Bar

###############################################################################
silent = False
barFlag = False
fname = 'CHD_HKDRGWA_6445CC_Clean_Ch1-22'
lociName = 'loci.txt'
locipairsName = 'lociPairs.txt'
buffer = 0
outputname = 'genepair_result_bt.txt'
covFile = None
filePath = os.getcwd()
corrMethod = 'svd' # {'chol' or 'svd' or 'eigen'}
loci = []
locipairs = []
bimData = []
bimbp = []
bimsnpID = []
bimchr = []
numSingleInter = 0

# Clear output
if(os.path.isfile(outputname)):
    os.system('rm '+outputname)

print '**********************************************'
print 'Loading the data...'
startTime = datetime.datetime.now()

# Loading the loci information
print 'file: '+lociName
fLoci = open(lociName, mode = 'r')
for line in fLoci:
    loci.append(line.split())
fLoci.close()

# Loading the locipairs information
print 'file: '+locipairsName
fLocipairs = open(locipairsName, mode = 'r')
for line in fLocipairs:
    locipairs.append(line.split())
fLocipairs.close()

# Loading .bim file for the relationship between bp and snpID(e.g. rs)
print 'file: '+fname+'.bim'
fBim = open(fname+'.bim', mode = 'r')
for line in fBim:
    bimData.append(line.split())
bimchr = [ele[0] for ele in bimData]
bimsnpID = [ele[1] for ele in bimData]
bimbp = [ele[3] for ele in bimData]
fBim.close()

endTime = datetime.datetime.now()
print 'Time for loading the data: '+str((endTime-startTime).seconds)+'s'

# Do analysis for each loci pair
if barFlag:
    print '**********************************************'
    print 'Generating .set and .snps files for all loci pairs with >1 snp pair...'
    startTime = datetime.datetime.now() 
    stepLen = len(locipairs)/100
    bar = Bar('Processing', fill='=', suffix = '%(percent)d%%')

indTmp = 0
for locipairTmp in locipairs[0:2]: # Add [0:1] to only check the first loci pair
    if barFlag:
        if indTmp%stepLen == 0: bar.next()
    indTmp += 1

    # Currently skip those locipair that only contains one snp pair inside
    if len(locipairTmp)/2-1 == 1:
        numSingleInter += 1
        continue

    if not silent:
        print '**********************************************'
        print 'Interaction between loci '+locipairTmp[0]+' and loci '+locipairTmp[1]+'...'

    # Get the loci indices for current loci pair
    lociInd1 = int(locipairTmp[0])
    lociInd2 = int(locipairTmp[1])

    # Get the snps that are involved in each of current loci pair
    snpset1 = []
    snpset2 = []
    for snp in locipairTmp[2:]:
        if (snp in loci[lociInd1]) & (snp not in snpset1):
            snpset1.append(snp)
        elif (snp in loci[lociInd2]) & (snp not in snpset2):
            snpset2.append(snp)
        elif (snp not in loci[lociInd1]) & (snp not in loci[lociInd2]):
            print "Error 1 occurs!"

    if not silent:
        print 'Total number of snps involved: '+str(len(snpset1)+len(snpset2))
        print 'Number of snps in loci'+locipairTmp[0]+': '+str(len(snpset1))
        print 'Number of snps in loci'+locipairTmp[1]+': '+str(len(snpset2))
        print 'Total number of snp-snp interactions after BOOST thresholding: '+str(len(locipairTmp)/2-1)

    # Transform bp to SNP(e.g. rs)
    for i in range(len(snpset1)):
        snpset1[i] = bimsnpID[bimbp.index(snpset1[i])]
    for i in range(len(snpset2)):
        snpset2[i] = bimsnpID[bimbp.index(snpset2[i])]

    # Write the two sets and corresponding snps into temp.set and temp.snps
    fileSet = open('temp.set', mode = 'w')
    fileSNPs = open('temp.snps', mode = 'w')

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

    # Extracting temp.bed, temp.bim, temp.fam for current temp.snps
    print '----------------------------------------------'
    print 'Extracting temp.bed/.bim/.fam files from current temp.snps...'
    sTime2 = datetime.datetime.now()
    os.system('plink1.9.exe --bfile '+fname+' --extract temp.snps --out temp --make-bed --silent --noweb')
    eTime2 = datetime.datetime.now()
    print 'Time for extracting temp.bed/.bim/.fam: '+str((eTime2-sTime2).seconds)+'s'

    # Calculating ld contrast test statistics and corr mat
    print '----------------------------------------------'
    print 'Calculating ld contrast test statistics and corr mat...'
    sTime2 = datetime.datetime.now()
    os.system("plink_Sen.exe --bfile temp --fast-epistasis --set-test --set temp.set --epi1 1 --epi2 1 --out " + fname + "_ld --silent --noweb")
    eTime2 = datetime.datetime.now()
    print 'Time for calculating ld contrast test statistis and corr mat: '+str((eTime2-sTime2).seconds)+'s'

    # Combining p values using R
    print '----------------------------------------------'
    print 'Combining p values using R...'
    sTime2 = datetime.datetime.now()
    os.system('Rscript.exe gene_based_inter_bt.R '+fname+'_ld.epi.cc '+fname+'_ld.epi.cc.summary.corr '+'loci'+locipairTmp[0]+' '+'loci'+locipairTmp[1]+' '+outputname+' '+filePath+' '+corrMethod)
    eTime2 = datetime.datetime.now()
    print 'Time for combining p values using R: '+str((eTime2-sTime2).seconds)+'s'

    # Arrange and clean temp files
    os.system('mv '+fname+'_ld.epi.cc loci'+locipairTmp[0]+'_loci'+locipairTmp[1]+'.cc')
    os.system('mv '+fname+'_ld.epi.cc.summary.corr loci'+locipairTmp[0]+'_loci'+locipairTmp[1]+'.cc.corr')
    # os.system('rm '+fname+'_ld.epi.cc.* temp.* *.snps *.set ')

    stop = 1

if barFlag:
    bar.finish()
    endTime = datetime.datetime.now()
    print 'Time forgenrating all the temp.set and temp.snps files: '+str((endTime-startTime).seconds)+'s'

print '----------------------------------------------'
print 'The number of locipairs that are analyzed: '+str(indTmp-numSingleInter)
print 'The number of locipairs that only contain one snp pair inside: '+str(numSingleInter)+' (skipped)'
stop = 1