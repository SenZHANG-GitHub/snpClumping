import os, datetime

# fPath = 'F:\\Data\\ToUse\\'
fname = 'Bowden_study_subjects_PLINK_QC'
gene_interloc = 'gene_test_list_example2.txt'
buffer = 0
outputname = 'genepair_result_bt.txt'
covFile = None
filePath = os.getcwd()
corrMethod = 'svd'  # {'chol' or 'svd' or 'eigen'}

if(os.path.isfile(outputname)):
    os.system("rm "+outputname)

fileGenes = open(gene_interloc, mode = 'r')
genes = []
for line in fileGenes:
    genes.append(line.split())
fileGenes.close()

fileBim = open(fname+'.bim', mode = 'r')
bimData = []
for line in fileBim:
    bimData.append(line.split())
fileBim.close()

fileFam = open(fname+'.fam', mode = 'r')
famData = []
for line in fileFam:
    famData.append(line.split())
fileFam.close()

# nested loop for get two gene information for gene_gene interaction calculate
for i in range(len(genes)):
    gene1 = genes[i][3]
    for j in range(i+1, len(genes)):
        gene2 = genes[j][3]
        print '---------------------------------------------'
        print 'interaction bettwen '+gene1+' and '+gene2
        gene1stbp = int(genes[i][1])
        gene1endbp = int(genes[i][2])
        gene1chr = genes[i][0]
        gene2stbp = int(genes[j][1])
        gene2endbp = int(genes[j][2])
        gene2chr = genes[j][0]
        # buffers
        buffstart1 = gene1stbp - buffer
        buffstart2 = gene2stbp - buffer
        buffstop1 = gene1endbp + buffer
        buffstop2 = gene2endbp + buffer

        gene1Bim = [ele for ele in bimData if ele[0]==gene1chr]
        gene2Bim = [ele for ele in bimData if ele[0]==gene2chr]
        firstbp1 = int(gene1Bim[0][3])
        lastbp1 = int(gene1Bim[-1][3])
        firstbp2 = int(gene2Bim[0][3])
        lastbp2 = int(gene2Bim[-1][3])

        # make sure regions don't overshoot the boundaries
        if buffstart1 < firstbp1: buffstart1 = firstbp1
        if buffstart2 < firstbp2: buffstart2 = firstbp2
        if buffstop1 > lastbp1: buffstop1 = lastbp1
        if buffstop2 > lastbp2: buffstop2 = lastbp2

        # only do this if both genes are within the region of our data
        # do this by checking bp, need to make sure of the genome assembly version (gene v.s. data)
        if (buffstart1 < lastbp1) & (buffstart2 < lastbp2) & (buffstop1 > firstbp1) & (buffstop2 > firstbp2):
            fileSet = open(fname+'.set', mode = 'w')
            fileSNPs = open(fname+'_temp.snps', mode = 'w')

            gene1SetBim = [ele for ele in gene1Bim if buffstart1 <= int(ele[3]) <= buffstop1]
            gene2SetBim = [ele for ele in gene2Bim if buffstart2 <= int(ele[3]) <= buffstop2]

            fileSet.write('SET_1'+'\n')
            for ele in gene1SetBim:
                fileSet.write(ele[1]+'\n')
                fileSNPs.write(ele[1]+'\n')

            fileSet.write('END'+'\n\n')
            fileSet.write('SET_2'+'\n')

            for ele in gene2SetBim:
                fileSet.write(ele[1]+'\n')
                fileSNPs.write(ele[1]+'\n')
            fileSet.write('END'+'\n')

            fileSet.close()
            fileSNPs.close()

            numsnps = len(gene1SetBim)+len(gene2SetBim)
            numInd = len(famData)

            if (len(gene1SetBim) > 0) & (len(gene2SetBim) > 0):
                print '---------------------------------------------'
                print 'Separating the snps within these two genes...'
                startTime = datetime.datetime.now()
                os.system('xwas.exe --bfile '+fname+' --extract '+
                          fname+'_temp.snps --out temp_bt --make-bed --silent --noweb')
                endTime = datetime.datetime.now()
                print 'Time for separating the snps: '+str((endTime-startTime).seconds)+'s'

                print '---------------------------------------------'
                print 'Calculating the interaction test stats and corr mat...'
                startTime = datetime.datetime.now()
                os.system("plink_Sen.exe --bfile temp_bt --fast-epistasis --set-test --set " +
                          fname+".set --epi1 1 --epi2 1 --out " + fname + "_ld --silent --noweb")
                endTime = datetime.datetime.now()
                print 'Time for calculating the interaction test stats and corr mat: '+str((endTime-startTime).seconds)+'s'

                print '---------------------------------------------'
                print 'Calculating gene-based interaction p-values using R...'
                startTime = datetime.datetime.now()
                os.system('Rscript.exe gene_based_inter_bt.R '+fname+'_ld.epi.cc '+fname+'_ld.epi.cc.summary.corr '+
                           gene1+' '+gene2+' '+outputname+' '+filePath+' '+corrMethod)
                endTime = datetime.datetime.now()
                print 'Time for calculating gene-based interaction p-values: '+str((endTime-startTime).seconds)+'s'

                os.system('mv '+fname+'_ld.epi.cc '+gene1+'_'+gene2+'.cc')
                os.system('mv '+fname+'_ld.epi.cc.summary.corr '+gene1+'_'+gene2+'.cc.corr')
                os.system('rm '+fname+'_ld.epi.cc.* temp_bt.* *.snps *.set ')


stop = 1


