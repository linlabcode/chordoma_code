#!/usr/bin/python
#pipeline_template.py

'''
The MIT License (MIT)

Copyright (c) 2016 Charles Lin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
'''

#generic pipeline template for human data


#==========================================================================
#=============================DEPENDENCIES=================================
#==========================================================================


import sys
sys.path.append('/ark/home/cl512/src/pipeline/')

import pipeline_dfci
import utils
import string
import os
from collections import defaultdict
#==========================================================================
#============================PARAMETERS====================================
#==========================================================================



projectName = 'chordoma_fix'
dataFile = '/grail/projects/%s/CHORDOMA_PRIMARY_TABLE.txt' % (projectName)
genome ='hg19'
annotFile = '/ark/home/cl512/src/pipeline/annotation/%s_refseq.ucsc' % (genome)

#project folders
projectFolder = '/grail/projects/%s/' % (projectName) #PATH TO YOUR PROJECT FOLDER

#standard folder names
gffFolder ='%sgff/' % (projectFolder)
macsFolder = '%smacsFolder/' % (projectFolder)
macsEnrichedFolder = '%smacsEnriched/' % (projectFolder)
mappedEnrichedFolder = '%smappedEnriched/' % (projectFolder)
mappedFolder = '%smappedFolder/' % (projectFolder)
wiggleFolder = '%swiggles/' % (projectFolder)
metaFolder = '%smeta/' % (projectFolder)

#making folders
folderList = [gffFolder,macsFolder,macsEnrichedFolder,mappedEnrichedFolder,mappedFolder,wiggleFolder,metaFolder]

for folder in folderList:
    pipeline_dfci.formatFolder(folder,True)


primary_list = ['PRIMARY_CHOR_269A7_H3K27AC','PRIMARY_CHOR_A269A1_H3K27AC','PRIMARY_CHOR_01192016_H3K27AC','PRIMARY_CHOR_12132013_H3K27AC','PRIMARY_CHOR_142A2_H3K27AC','PRIMARY_CHOR_194A1_H3K27AC','PRIMARY_CHOR_243A2_H3K27AC','PRIMARY_CHOR_4616_H3K27AC']



all_data_file = '%sCHORDOMA_H3K27AC_ALL_TABLE.txt' % (projectFolder)
lines_list = ['MUGCHOR1_H3K27AC','UCH1_H3K27AC','UCH2_H3K27AC','UMCHOR1_H3K27AC','JHC7_H3K27AC']



#==========================================================================
#========================FORMATTING SAMPLE TABLE===========================
#==========================================================================

##THIS SECTION CREATES A DATA TABLE FROM A WHITEHEAD ANNOTATION SPREADSHEET

##give full path
##sampleTableFile = 'YOUR_WIGTC_ANNOTATION.xls' #<- the .xls file in the seq data folder provided by WI

#dirpath = ''  <- provide full path of folder containing raw seq files
##e.g. /ark/home/jr246/raw/130925_..../QualityScore/

##bamPath <- where we store our bams.  Must have write access if you want to call bowtie
##e.g. /ark/home/jr246/bam/
#bamPath = '/ark/home/jr246/bam/'

#pipeline_dfci.makePipelineTable(sampleTableFile,dirPath,bamPath,dataFile)

#dataDict = pipeline_dfci.loadDataTable(dataFile)

#namesList = dataDict.keys()

#print(namesList)

#==========================================================================
#=======================LOADING DATA ANNOTATION============================
#==========================================================================

##THIS SECTION LOADS A DATA TABLE.  MUST BE UNCOMMENTED FOR REST OF CODE TO WORK


#LOADING THE DATA TABLE
dataDict = pipeline_dfci.loadDataTable(dataFile)
print(dataDict.keys())

pipeline_dfci.summary(dataFile)


#checking background assignment
namesList = dataDict.keys()
print(namesList)
for name in namesList:
    if name.count('H3K27AC') == 1:
        print(name)
        print(dataDict[name]['background'])



#==========================================================================
#==========================CALLING BOWTIE==================================
#==========================================================================

##THIS SECTION CALLS BOWTIE ON RAW READ FILES TO GENERATE SORTED AND INDEXED BAMS IN THE BAM FOLDER


#namesList = []  <- fill this in if you want to only map a subset of the data. otherwise leave blank

##SET LAUNCH TO False to debug
#pipeline_dfci.makeBowtieBashJobs(dataFile,namesList,launch=True)



#==========================================================================
#============================MERGING BAMS==================================
#==========================================================================


# #merge all the bams

# mergeCmd = 'samtools merge /grail/projects/chordoma_fix/bam/chordoma_H3K27AC_merged.bam '

# namesList = [name for name in dataDict.keys() if name.count('H3K27AC') == 1]
# for name in namesList:
#     mergeCmd += dataDict[name]['bam']
#     mergeCmd += ' '

# print(mergeCmd)

# os.system(mergeCmd)

#==========================================================================
#=============================CALL MACS====================================
#==========================================================================

##THIS SECTION CALLS THE MACS ERROR MODEL


# namesList = [name for name in dataDict.keys() if name.count('H3K27AC')]

# pipeline_dfci.callMacs(dataFile,macsFolder,namesList,overwrite=False,pvalue='1e-9',useBackground=True)



#line_data_file = '%sCHORDOMA_TABLE.txt' % (projectFolder)

#namesList = ['UCH1_T']
#pipeline_dfci.callMacs(line_data_file,macsFolder,namesList,overwrite=False,pvalue='1e-9',useBackground=True)






#==========================================================================
#=======================FORMAT MACS OUTPUT=================================
#==========================================================================

##THIS SECTION FORMATS THE OUTPUT FROM MACS, CREATES THE MACSENRICHED FOLDER AND MOVES WIGGLES TO THE DESTINATION

#pipeline_dfci.formatMacsOutput(dataFile,macsFolder,macsEnrichedFolder,wiggleFolder,wigLink='/ark/wiggles/',useBackground=True)



#pipeline_dfci.formatMacsOutput(line_data_file,macsFolder,macsEnrichedFolder,wiggleFolder,wigLink='/ark/wiggles/',useBackground=True)

#==========================================================================
#===========================MAKE GENE GFFS=================================
#==========================================================================

#pipeline_dfci.makeGeneGFFs(annotFile,gffFolder,species='HG19')


#==========================================================================
#=========================MAP ENRICHED GFF=================================
#==========================================================================

# setName = 'PRIMARY_CHORDOMA'
# gffList = ['%sHG19_TSS_ALL_-1000_+1000.gff' % (gffFolder)]
# cellTypeList = ['PRIMARY']
# enrichedFolder = macsEnrichedFolder

# pipeline_dfci.mapEnrichedToGFF(dataFile,setName,gffList,cellTypeList,enrichedFolder,mappedEnrichedFolder,macs=True,namesList=primary_list,useBackground=True)


# setName = 'LINES_CHORDOMA'
# gffList = ['%sHG19_TSS_ALL_-1000_+1000.gff' % (gffFolder)]

# cellTypeList = ['MUGCHOR1','UCH1','UCH2','UMCHOR1']
# enrichedFolder = macsEnrichedFolder

# pipeline_dfci.mapEnrichedToGFF(all_data_file,setName,gffList,cellTypeList,enrichedFolder,mappedEnrichedFolder,macs=True,namesList=lines_list,useBackground=True)

#==========================================================================
#==========================MAKE GENE LISTS=================================
#==========================================================================

# mappedEnrichedFile = '%sHG19_TSS_ALL_-1000_+1000/HG19_TSS_ALL_-1000_+1000_PRIMARY_CHORDOMA.txt' % (mappedEnrichedFolder)

# setList = [[x] for x in primary_list]
# print(setList)
# output = '%sgffListFiles/HG19_PRIMARY_CHORDOMA_ACTIVE.txt' % (projectFolder)
# pipeline_dfci.makeGFFListFile(mappedEnrichedFile,setList,output,annotFile)


# mappedEnrichedFile = '%sHG19_TSS_ALL_-1000_+1000/HG19_TSS_ALL_-1000_+1000_LINES_CHORDOMA.txt' % (mappedEnrichedFolder)

# setList = [[x] for x in lines_list]
# print(setList)
# output = '%sgffListFiles/HG19_LINES_CHORDOMA_ACTIVE.txt' % (projectFolder)
# pipeline_dfci.makeGFFListFile(mappedEnrichedFile,setList,output,annotFile)



#==========================================================================
#============================SUMMARIZE DATA================================
#==========================================================================


def summarizeData(dataFile,output ='',namesList= []):

    dataDict=pipeline_dfci.loadDataTable(dataFile)

    if len(namesList) == 0:
        namesList = dataDict.keys()

    if len(output) == 0:
        output = string.replace(dataFile,'.txt','_SUMMARY.txt')

    print('WRITING OUTPUT TO %s' % (output))
    readTable = [['NAME','TOTAL_READS','MAPPED_READS','PEAKS']]

    for name in namesList:
        print('GETTING DATA SUMMARY FOR %s' % (name))

        uniqueID = dataDict[name]['uniqueID']

        mappedReads = round(float(pipeline_dfci.getTONYInfo(uniqueID,'67'))/1000000,2)
        totalRaw = pipeline_dfci.getTONYInfo(uniqueID,'68')
        totalRaw = int(totalRaw.split('::')[0])
        totalReads = round(float(totalRaw)/1000000,2)
        #mappedReads = 0
        #totalReads = 0

        #getting the spot score
        #spotFile = '%sspot/%s_%s/%s_hg19.sorted.spot.out' % (projectFolder,uniqueID,name,uniqueID)
        #spotFile = '%sspot/%s_%s/%s_hg19.sorted.spot.out' % (projectFolder,uniqueID,name,uniqueID)
        #spotTable = utils.parseTable(spotFile,'\t')
        #spotScore = spotTable[1][0].split(' ')[-1]

        #get the peak count
        peakCollection = utils.importBoundRegion('%s%s' % (macsEnrichedFolder,dataDict[name]['enrichedMacs']),name)
        peakCount = len(peakCollection)



        newLine = [name,totalReads,mappedReads,peakCount]
        print(newLine)
        readTable.append(newLine)



    utils.unParseTable(readTable,output,'\t')    

output = '%stables/PRIMARY_CHORDOMA_H3K27AC_SUMMARY.txt' % (projectFolder)
# summarizeData(dataFile,output,namesList= [name for name in dataDict.keys() if name.count('H3K27AC') == 1])

#==========================================================================
#======================TEST MAPPING W/ BATCH===============================
#==========================================================================
# gffList = ['%sHG19_CHORDOMA_FIGURE_GENES.gff' % (gffFolder)]
# namesList = ['PRIMARY_CHOR_A269A1_H3K27AC']

# pipeline_dfci.mapBamsBatch(dataFile,gffList,mappedFolder,False,namesList,extension=200)


#==========================================================================
#==============================CALLING ROSE================================
#==========================================================================

# #CALLING ROSE ON H3K27AC DATA
# namesList = [name for name in dataDict.keys() if name.upper().count("H3K27AC") == 1]

# parentFolder = '%srose' % (projectFolder)
# parentFolder = utils.formatFolder(parentFolder,True)

# maskFile ='/raider/index/hg19/Masks/hg19_encode_blacklist.bed'
# bashFileName = '%srose/primary_chordoma_h3k27ac_rose.sh' %(projectFolder)

# #all datasets
# #2500 tss exclusion
# #auto stitching
# pipeline_dfci.callRose2(dataFile,macsEnrichedFolder,parentFolder,namesList,[],'',2500,'',bashFileName,maskFile)





#==========================================================================
#===========================CALLING META ROSE==============================
#==========================================================================

#  #for primary
# analysisName = 'primary_chordoma_h3K27ac_final'

# primary_list = ['PRIMARY_CHOR_A269A1_H3K27AC',
#                 'PRIMARY_CHOR_01192016_H3K27AC',
#                 'PRIMARY_CHOR_12132013_H3K27AC',
#                 'PRIMARY_CHOR_142A2_H3K27AC',
#                 'PRIMARY_CHOR_194A1_H3K27AC',
#                 'PRIMARY_CHOR_243A2_H3K27AC',
#                 'PRIMARY_CHOR_4616_H3K27AC',
#                 '160831_CHOR_H3K27Ac']

# bamFileList = [dataDict[name]['bam'] for name in primary_list]
# bamString = string.join(bamFileList,',')

# controlBams = [dataDict[name]['background'] for name in primary_list]

# controlFileList = [dataDict[name]['bam'] for name in controlBams]
# controlBamString = string.join(controlFileList,',')

# bedFileList = [macsEnrichedFolder + dataDict[name]['enrichedMacs'] for name in primary_list]
# bedString = string.join(bedFileList,',')

# roseFolder = '%smeta_rose/' % (projectFolder)
# roseFolder = utils.formatFolder(roseFolder,True)

# outputFolder = '%s%s/' % (roseFolder,analysisName)
# bashFileName = '%s%s_meta_rose.sh' % (roseFolder,analysisName)

# bashFile = open(bashFileName,'w')
# bashFile.write('#!/usr/bin/bash\n\n')
# bashFile.write('cd /ark/home/cl512/pipeline/\n')


# maskFile ='/raider/index/hg19/Masks/hg19_encode_blacklist.bed'
# metaRoseCmd = 'python /ark/home/cl512/pipeline/ROSE2_META.py -g hg19 -i %s -r %s -c %s -o %s -n %s -t 2500 --mask %s' % (bedString,bamString,controlBamString,outputFolder,analysisName,maskFile)

# bashFile.write(metaRoseCmd + '\n')




#==========================================================================
#===========================PLOTTING BY GFF================================
#==========================================================================

# figureGFF = [['chr7','CAV','CAV',115827544,116446736,'','+','','CAV'],
#              ['chr6','T','T',166508099,166598920,'','+','','T'],
#              ['chr6','T_AMPLICON','T_AMPLICON',165328878,167941884,'','+','','T_AMPLICON'],
#              ['chr6','CTGF','CTGF',132262895,132281762,'','-','','CTGF'],
#              ['chr7','TWIST1','TWIST1',19141898,19166920,'','-','','TWIST1'],
#              ['chr1','PDE4B','PDE4B',66615950,66984899,'','+','','PDE4B'],
#              ]

# figureGFFFile = '%sHG19_PRIMARY_CHORDOMA_FIGURE_GENES.gff' % (gffFolder)
# utils.unParseTable(figureGFF,figureGFFFile,'\t')
# outputFolder = '%sgenePlot' % (projectFolder)


# #making meta plots of k27ac
# plotList = [name for name in dataDict.keys() if name.count('H3K27AC') == 1]
# plotName = 'PRIMARY_FIGURE_GENES_H3K27AC_META'
# groupString = string.join(['H3K27AC']*len(plotList),',')
# print(groupString)

# pipeline_dfci.callBatchPlot(dataFile,figureGFFFile,plotName,outputFolder,namesList=plotList,uniform=True,bed ='',plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString)



# #individual k27ac


# plotList = ['PRIMARY_CHOR_269A7_H3K27AC','PRIMARY_CHOR_A269A1_H3K27AC','PRIMARY_CHOR_01192016_H3K27AC','PRIMARY_CHOR_12132013_H3K27AC','PRIMARY_CHOR_142A2_H3K27AC','PRIMARY_CHOR_194A1_H3K27AC','PRIMARY_CHOR_243A2_H3K27AC','PRIMARY_CHOR_4616_H3K27AC','160831_CHOR_H3K27Ac','160831_NAT_H3K27Ac','SF8894_H3K27Ac']
# plotName = 'PRIMARY_FIGURE_GENES_H3K27AC_MULTIPLE'


# pipeline_dfci.callBatchPlot(dataFile,figureGFFFile,plotName,outputFolder,namesList=plotList,uniform=True,bed ='',plotType= 'MULTIPLE',extension=200,multiPage = False,debug=False,nameString = '')


# plotList = ['PRIMARY_CHOR_269A7_H3K27AC','PRIMARY_CHOR_A269A1_H3K27AC']
# plotName = 'MATCHED_FIGURE_GENES_H3K27AC_MULTIPLE'


# pipeline_dfci.callBatchPlot(dataFile,figureGFFFile,plotName,outputFolder,namesList=plotList,uniform=True,bed ='',plotType= 'MULTIPLE',extension=200,multiPage = False,debug=False,nameString = '')




# #lines k27ac
# line_data_file = '%sCHORDOMA_TABLE.txt' % (projectFolder)

# line_dataDict = pipeline_dfci.loadDataTable(line_data_file)

# plotList = [name for name in line_dataDict.keys() if name.count('H3K27AC') == 1]

# plotName = 'LINES_FIGURE_GENES_H3K27AC_MULTIPLE'
# pipeline_dfci.callBatchPlot(line_data_file,figureGFFFile,plotName,outputFolder,namesList=plotList,uniform=True,bed ='',plotType= 'MULTIPLE',extension=200,multiPage = False,debug=False,nameString = '')


# plotName = 'LINES_FIGURE_GENES_H3K27AC_MULTIPLE_RELATIVE'
# pipeline_dfci.callBatchPlot(line_data_file,figureGFFFile,plotName,outputFolder,namesList=plotList,uniform=False,bed ='',plotType= 'MULTIPLE',extension=200,multiPage = False,debug=False,nameString = '')


# plotList = ['UCH1_T']

# plotName = 'T_FIGURE_GENES'
# pipeline_dfci.callBatchPlot(line_data_file,figureGFFFile,plotName,outputFolder,namesList=plotList,uniform=True,bed ='',plotType= 'MULTIPLE',extension=200,multiPage = False,debug=False,nameString = '')



#==========================================================================
#==========================CALLING CLUSTER=================================
#==========================================================================


#across all primary datasets
#analysis name
analysisName = 'primary_chordoma_h3k27ac_final'

#dataFile

primary_data_file = '%sCHORDOMA_PRIMARY_TABLE.txt' % (projectFolder)
dataDict = pipeline_dfci.loadDataTable(dataFile)


#nameslist
namesList = ['PRIMARY_CHOR_A269A1_H3K27AC',
             'PRIMARY_CHOR_01192016_H3K27AC',
             'PRIMARY_CHOR_12132013_H3K27AC',
             'PRIMARY_CHOR_142A2_H3K27AC',
             'PRIMARY_CHOR_194A1_H3K27AC',
             'PRIMARY_CHOR_243A2_H3K27AC',
             'PRIMARY_CHOR_4616_H3K27AC',
             '160831_CHOR_H3K27Ac',
             ]



namesString = string.join(namesList,',')

#set up the output
utils.formatFolder('%sclustering/' % (projectFolder),True)
outputFolder = '%sclustering/%s/' % (projectFolder,analysisName)
outputFolder = utils.formatFolder(outputFolder,True)

#get the rose folder
roseFolder ='%srose/' % (projectFolder)
roseFolder = utils.formatFolder('%srose/' % (projectFolder),True)

#get the mask file
maskFile ='/raider/index/hg19/Masks/hg19_encode_blacklist.bed'

#set up the bash file
bashFileName = '%s%s_clustering.sh' % (outputFolder,analysisName)
bashFile = open(bashFileName,'w')

bashFile.write('#!/usr/bin/bash\n')

#for now change into pipelinedir just to be safe
bashFile.write('cd /ark/home/cl512/pipeline/\n')


clusterCmd = 'python /ark/home/cl512/pipeline/clusterEnhancer.py -d %s -i %s -r %s -o %s -n %s -e super -t 2500 --mask %s' % (primary_data_file,namesString,roseFolder,outputFolder,analysisName,maskFile)

bashFile.write(clusterCmd + '\n')
print(clusterCmd)

bashFile.close()









#=======================================================

#across all  datasets
#analysis name
analysisName = 'all_chordoma_h3k27ac_final'

#dataFile

all_data_file = '%sCHORDOMA_H3K27AC_ALL_TABLE.txt' % (projectFolder)
dataDict = pipeline_dfci.loadDataTable(all_data_file)
pipeline_dfci.summary(all_data_file)


#nameslist


namesList = ['PRIMARY_CHOR_A269A1_H3K27AC',
             'PRIMARY_CHOR_01192016_H3K27AC',
             'PRIMARY_CHOR_12132013_H3K27AC',
             'PRIMARY_CHOR_142A2_H3K27AC',
             'PRIMARY_CHOR_194A1_H3K27AC',
             'PRIMARY_CHOR_243A2_H3K27AC',
             'PRIMARY_CHOR_4616_H3K27AC',
             '160831_CHOR_H3K27Ac',
             ]

namesList = namesList + lines_list + ['160831_NAT_H3K27Ac']
print(namesList)
namesString = string.join(namesList,',')

#set up the output
utils.formatFolder('%sclustering/' % (projectFolder),True)
outputFolder = '%sclustering/%s/' % (projectFolder,analysisName)
outputFolder = utils.formatFolder(outputFolder,True)

#get the rose folder
roseFolder ='%srose/' % (projectFolder)
roseFolder = utils.formatFolder('%srose/' % (projectFolder),True)

#get the mask file
maskFile ='/raider/index/hg19/Masks/hg19_encode_blacklist.bed'

#set up the bash file
bashFileName = '%s%s_clustering.sh' % (outputFolder,analysisName)
bashFile = open(bashFileName,'w')

bashFile.write('#!/usr/bin/bash\n')

#for now change into pipelinedir just to be safe
bashFile.write('cd /ark/home/cl512/pipeline/\n')


clusterCmd = 'python /ark/home/cl512/pipeline/clusterEnhancer.py -d %s -i %s -r %s -o %s -n %s -e super -t 2500 --mask %s' % (all_data_file,namesString,roseFolder,outputFolder,analysisName,maskFile)

bashFile.write(clusterCmd + '\n')
print(clusterCmd)

bashFile.close()


#==========================================================================
#==========================MAKING STITCHED GFF=============================
#==========================================================================

#making stitched gff for primaries

# primary_data_dict = pipeline_dfci.loadDataTable(dataFile)

# all_loci = []
# for name in primary_list:
#     print(name)
#     regions = utils.importBoundRegion('%s%s' % (macsEnrichedFolder,primary_data_dict[name]['enrichedMacs']),name)
#     all_loci+= regions.getLoci()

# all_collection = utils.LocusCollection(all_loci,50)

# stitched_collection = all_collection.stitchCollection()

# gff = utils.locusCollectionToGFF(stitched_collection)

# gff_path = '%sHG19_PRIMARY_CHORDOMA_STITCHED_ALL_-0_+0.gff' % (gffFolder)

# utils.unParseTable(gff,gff_path,'\t')



# lines_data_dict = pipeline_dfci.loadDataTable(all_data_file)

# print(lines_data_dict.keys())

# all_loci = []
# for name in lines_list:
#     print(name)
#     regions = utils.importBoundRegion('%s%s' % (macsEnrichedFolder,lines_data_dict[name]['enrichedMacs']),name)
#     all_loci+= regions.getLoci()

# all_collection = utils.LocusCollection(all_loci,50)

# stitched_collection = all_collection.stitchCollection()

# gff = utils.locusCollectionToGFF(stitched_collection)

# gff_path = '%sHG19_LINES_CHORDOMA_STITCHED_ALL_-0_+0.gff' % (gffFolder)

# utils.unParseTable(gff,gff_path,'\t')


#==========================================================================
#======================RUNNING ENHANCER PROMOTER===========================
#==========================================================================

# #for all primary

# parentFolder = utils.formatFolder('%senhancerPromoter/' % (projectFolder),True)
# analysis_name = 'primary_chordoma'

# parentFolder = utils.formatFolder('%senhancerPromoter/%s/' % (projectFolder,analysis_name),True)

# bash_path = '%s%s_enhancerPromoter.sh' % (parentFolder,analysis_name)

# bash_file = open(bash_path,'w')

# bash_file.write('#!/usr/bin/bash\n\n')
# bash_file.write('cd /ark/home/cl512/pipeline/\n\n')
# namesList = list(primary_list)
# dataDict = pipeline_dfci.loadDataTable(dataFile)
# gff_path = '%sHG19_PRIMARY_CHORDOMA_STITCHED_ALL_-0_+0.gff' % (gffFolder) 
# activity_path = '%sgffListFiles/HG19_PRIMARY_CHORDOMA_ACTIVE.txt' % (projectFolder)
# for name in namesList:

#     bam_path = dataDict[name]['bam']
#     background_name = dataDict[name]['background']
#     background_bam_path = dataDict[background_name]['bam']
#     output_folder = utils.formatFolder('%s%s/' % (parentFolder,name),True)

#     cmd = '# python enhancerPromoter.py -i %s -g HG19 -a %s --name %s -o %s -b %s -c %s' % (gff_path,activity_path,name,output_folder,bam_path,background_bam_path)
#     bash_file.write(cmd+'\n')


# bash_file.close()





# #for all lines

# parentFolder = utils.formatFolder('%senhancerPromoter/' % (projectFolder),True)
# analysis_name = 'lines_chordoma'

# parentFolder = utils.formatFolder('%senhancerPromoter/%s/' % (projectFolder,analysis_name),True)

# bash_path = '%s%s_enhancerPromoter.sh' % (parentFolder,analysis_name)

# bash_file = open(bash_path,'w')

# bash_file.write('#!/usr/bin/bash\n\n')
# bash_file.write('cd /ark/home/cl512/pipeline/\n\n')
# namesList = list(lines_list)
# dataDict = pipeline_dfci.loadDataTable(all_data_file)
# gff_path = '%sHG19_LINES_CHORDOMA_STITCHED_ALL_-0_+0.gff' % (gffFolder) 
# activity_path = '%sgffListFiles/HG19_LINES_CHORDOMA_ACTIVE.txt' % (projectFolder)
# for name in namesList:

#     bam_path = dataDict[name]['bam']
#     background_name = dataDict[name]['background']
#     background_bam_path = dataDict[background_name]['bam']
#     output_folder = utils.formatFolder('%s%s/' % (parentFolder,name),True)

#     cmd = '# python enhancerPromoter.py -i %s -g HG19 -a %s --name %s -o %s -b %s -c %s' % (gff_path,activity_path,name,output_folder,bam_path,background_bam_path)
#     bash_file.write(cmd+'\n')


# bash_file.close()



# #for primary merge
# # #for all primary

# parentFolder = utils.formatFolder('%senhancerPromoter/' % (projectFolder),True)
# analysis_name = 'primary_chordoma_merge'

# parentFolder = utils.formatFolder('%senhancerPromoter/%s/' % (projectFolder,analysis_name),True)

# bash_path = '%s%s_enhancerPromoter.sh' % (parentFolder,analysis_name)

# bash_file = open(bash_path,'w')

# bash_file.write('#!/usr/bin/bash\n\n')
# bash_file.write('cd /ark/home/cl512/pipeline/\n\n')
# namesList = list(primary_list)
# namesList.remove('PRIMARY_CHOR_269A7_H3K27AC')
# print(namesList)
# dataDict = pipeline_dfci.loadDataTable(dataFile)
# gff_path = '%sHG19_PRIMARY_CHORDOMA_STITCHED_ALL_-0_+0.gff' % (gffFolder) 
# activity_path = '%sgffListFiles/HG19_PRIMARY_CHORDOMA_ACTIVE.txt' % (projectFolder)

# bamList= []
# controlList = []
# for name in namesList:

#     bam_path = dataDict[name]['bam']
#     background_name = dataDict[name]['background']
#     background_bam_path = dataDict[background_name]['bam']
#     bamList.append(bam_path)
#     controlList.append(background_bam_path)

# bam_string = ' '.join(bamList)
# control_string = ' '.join(controlList)


# cmd = 'python enhancerPromoter.py -i %s -g HG19 -a %s --name %s -o %s -b %s -c %s' % (gff_path,activity_path,analysis_name,parentFolder,bam_string,control_string)
# bash_file.write(cmd+'\n')


# bash_file.close()






# #for a merge of all lines

# parentFolder = utils.formatFolder('%senhancerPromoter/' % (projectFolder),True)
# analysis_name = 'lines_chordoma_merge'

# parentFolder = utils.formatFolder('%senhancerPromoter/%s/' % (projectFolder,analysis_name),True)

# bash_path = '%s%s_enhancerPromoter.sh' % (parentFolder,analysis_name)

# bash_file = open(bash_path,'w')

# bash_file.write('#!/usr/bin/bash\n\n')
# bash_file.write('cd /ark/home/cl512/pipeline/\n\n')
# namesList = list(lines_list)

# print(namesList)
# dataDict = pipeline_dfci.loadDataTable(all_data_file)
# gff_path = '%sHG19_LINES_CHORDOMA_STITCHED_ALL_-0_+0.gff' % (gffFolder) 
# activity_path = '%sgffListFiles/HG19_LINES_CHORDOMA_ACTIVE.txt' % (projectFolder)

# bamList= []
# controlList = []
# for name in namesList:

#     bam_path = dataDict[name]['bam']
#     background_name = dataDict[name]['background']
#     background_bam_path = dataDict[background_name]['bam']
#     bamList.append(bam_path)
#     controlList.append(background_bam_path)

# bam_string = ' '.join(bamList)
# control_string = ' '.join(controlList)


# cmd = 'python enhancerPromoter.py -i %s -g HG19 -a %s --name %s -o %s -b %s -c %s' % (gff_path,activity_path,analysis_name,parentFolder,bam_string,control_string)
# bash_file.write(cmd+'\n')


# bash_file.close()


#==========================================================================
#========================PEAK TABLE BINNING================================
#==========================================================================



def binPeakTable(peak_table_path,activity_path,binSize = 1000000,output = ''):

    '''
    calculates the promoter/enahncer AUC signal
    across bins
    sets the output to the same path unless otherwise specified
    '''


    if len(output) == 0:
        output = string.replace(peak_table,'.txt','bin_table.txt')

        
    binSize = int(binSize)
    
    stepSize = binSize/2

    activityTable = utils.parseTable(activity_path,'\t')
    startDict = utils.makeStartDict(annotFile)
    tssLoci = []

    print('making tss collection for active genes')
    for line in activityTable:
        tssLoci.append(utils.makeTSSLocus(line[1],startDict,0,0))

    tssCollection = utils.LocusCollection(tssLoci,50)
    

    promoterDict = {}
    enhancerDict = {}
    tssDict = {}
    #hard wired for hg19
    chrom_path = '/ark/home/cl512/pipeline/annotation/hg19.chrom.sizes'

    chrom_table = utils.parseTable(chrom_path,'\t')

    chromDict = {}
    for line in chrom_table:
        chromDict[line[0]] = int(line[1])

    chromList = ['chr'+str(i) for i in range(1,23)] + ['chrX','chrY'] #set the hg19 chroms
    #need to seed the dict
    for chrom in chromList:
        promoterDict[chrom] = defaultdict(float)
        enhancerDict[chrom] = defaultdict(float)
        tssDict[chrom] =defaultdict(int) # dict to count active promoters
    #now as we iterate through the peak table

    peak_table = utils.parseTable(peak_table_path,'\t')
    print('filling in enhancer dict')
    for line in peak_table[1:]:

        chrom = line[1]
        
        signal = float(line[9])*int(line[4])

        #for approximation use the center coordinate to assign bin
        #every region should be in 2 bins
        center = (int(line[2]) + int(line[3]))/2

        first_bin = center/stepSize

        if center % stepSize < stepSize:
            second_bin = first_bin - 1
        else:
            second_bin = first_bin + 1

        if int(line[5]) == 1:
            promoterDict[chrom][first_bin] +=signal
            promoterDict[chrom][second_bin] +=signal
        else:
            enhancerDict[chrom][first_bin] +=signal
            enhancerDict[chrom][second_bin] +=signal
        

    #now load up the new peak table
    outTable = [['BIN','CHROM','START','STOP','TSS_COUNT','PROMOTER','ENHANCER']]
    print('making out table')
    for chrom in chromList:
        print(chrom)
        chromLength = chromDict[chrom]

        for i in range(chromLength/stepSize):
            bin_start = i*stepSize + 1
            bin_stop =  i*stepSize + binSize
            bin_locus = utils.Locus(chrom,bin_start,bin_stop,'.')
            overlapTSSCount = len(tssCollection.getOverlap(bin_locus,'both'))

            bin_id = '%s_%s' % (chrom,str(i+1))

            promoterSignal = promoterDict[chrom][i]
            enhancerSignal = enhancerDict[chrom][i]
            
            newLine = [bin_id,chrom,bin_start,bin_stop,overlapTSSCount,promoterSignal,enhancerSignal]
            outTable.append(newLine)


    utils.unParseTable(outTable,output,'\t')
    return outTable

# peak_file = '/grail/projects/chordoma_fix/enhancerPromoter/primary_chordoma_merge/primary_chordoma_merge/primary_chordoma_merge_PEAK_TABLE.txt'
# out_path = '%stables/primary_chordoma_merge_bins_2m.txt' % (projectFolder)

# activity_path = '%sgffListFiles/HG19_PRIMARY_CHORDOMA_ACTIVE.txt' %(projectFolder)


# outTable = binPeakTable(peak_file,activity_path,binSize = 2000000,output = out_path)


# peak_file = '/grail/projects/chordoma_fix/enhancerPromoter/lines_chordoma_merge/lines_chordoma_merge/lines_chordoma_merge_PEAK_TABLE.txt'

# out_path = '%stables/lines_chordoma_merge_bins_2m.txt' % (projectFolder)

# activity_path = '%sgffListFiles/HG19_LINES_CHORDOMA_ACTIVE.txt' %(projectFolder)


# outTable = binPeakTable(peak_file,activity_path,binSize = 2000000,output = out_path)




#==========================================================================
#================================GEO UPLOAD================================
#==========================================================================



# def makeGEOTable(dataFile,wiggleFolder,macsFolder,namesList,geoName,outputFolder =''):

#     '''
#     makes a geo table and a bash script to format everything
#     '''
#     dataDict = pipeline_dfci.loadDataTable(dataFile)

#     #first make a reverse wce dict

#     backgroundDict = {}
#     for name in namesList:

#         background = dataDict[name]['background']
#         backgroundDict[background] = name


#     outputFolder = pipeline_dfci.formatFolder(outputFolder,True)
#     bashFileName = '%s%s_bash.sh' % (outputFolder,geoName)

#     bashFile = open(bashFileName,'w')




#     geoTable = [['SAMPLE_NAME','TITLE','CELL_TYPE','PROCESSED_FILE','RAW_FILE','BARCODE']]


#     namesList.sort()

#     for name in namesList:

#         sampleName = dataDict[name]['uniqueID']
#         title = name
#         cell_type = name.split('_')[0]
#         processed_file = "%s.wig.gz" % (name)
#         raw_file = "%s.fastq.gz" % (name)
        
#         fastqFile = dataDict[name]['fastq']
#         uniqueID = dataDict[name]['uniqueID']
#         try:
#             barcode = pipeline_dfci.getTONYInfo(uniqueID,38)
#         except IndexError:
#             barcode = ''
#         newLine = [sampleName,title,cell_type,processed_file,raw_file,barcode]
#         geoTable.append(newLine)
        

#     utils.unParseTable(geoTable,"%s%s_meta.xls" % (outputFolder,geoName),'\t')

#     #now make the folder to hold everything and the relevant bash script
#     if len(outputFolder) == 0:
#         outputFolder ='./%s/' % (geoName)
        
#     else:
#         outputFolder = outputFolder + geoName + '/'
    
#     pipeline_dfci.formatFolder(outputFolder,True)

#     wiggleFolder = pipeline_dfci.formatFolder(wiggleFolder,False)
#     macsFolder = pipeline_dfci.formatFolder(macsFolder,False)


#     #now make the bash file
#     bashFile.write('#!/usr/bin/bash\n')
#     bashFile.write('cd %s\n' %(outputFolder))
#     bashFile.write('\n')

#     #write the untar command
#     for name in namesList:
#         fastqFile = dataDict[name]['fastq']
#         if len(fastqFile) == 0:
#             print "NO FASTQ FILE FOR %s" % (name)
#             continue
#         tarCmd = 'cp %s %s.fastq.gz\n' % (fastqFile,name)
#         bashFile.write(tarCmd)

#     bashFile.write('\n\n\n')
#     #write the wiggle cp command
#     for name in namesList:
#         if name.count('WCE') == 1:
#             refName = backgroundDict[name]
#             controlWiggleFile = '%s%s/%s_MACS_wiggle/control/%s_control_afterfiting_all.wig.gz' % (macsFolder,refName,refName,refName)
#             wigCmd = "cp '%s' %s.wig.gz\n" % (controlWiggleFile,name)
#             #wigCmd = "cp '%swceWiggles/%s_control_afterfiting_all.wig.gz' %s.wig.gz\n" % (wiggleFolder,refName,name)
#         else:
#             wigCmd = "cp '%s%s_treat_afterfiting_all.wig.gz' %s.wig.gz\n" % (wiggleFolder,name,name)
#         bashFile.write(wigCmd)

#     #write the md5sums for the wiggles
#     bashFile.write('\n\n\n')
#     bashFile.write("echo '' > md5sum.txt\n")
#     for name in namesList:
#         md5Cmd = 'md5sum %s.wig.gz >> md5sum.txt\n' % (name)
#         bashFile.write(md5Cmd)

#     #write md5sums for the fastqs    
#     for name in namesList:
#         md5Cmd = 'md5sum %s.fastq.gz >> md5sum.txt\n' % (name)
#         bashFile.write(md5Cmd)

#     #the big tar command
#     tarCmd = '#tar -cvzf %s.tar.gz %s\n' % (geoName,outputFolder)
#     bashFile.write(tarCmd)
#     bashFile.close()



# namesList = ['MUGCHOR','UCH2_1','UCH2_2']
# geoName = 'chordoma'
# #makeGEOTable(dataFile,wiggleFolder,macsEnrichedFolder,namesList,geoName,outputFolder ='%sgeo/' % (projectFolder))

#==========================================================================
#=============================CALLING CRC==================================
#==========================================================================



# enhancerFile = '%srose/UCH2_H3K27AC_ROSE/UCH2_H3K27AC_peaks_AllEnhancers.table.txt' % (projectFolder)

# bamFile = dataDict['UCH2_H3K27AC']['bam']
# backgroundBamFile = dataDict['UCH2_WCE']['bam']

# genome = 'HG19'

# outputFolder = '%snetwork/' % (projectFolder)
# outputFolder = utils.formatFolder(outputFolder,True)

# analysisName = 'UCH2'

# bashFileName = '%s%s_crc2.sh' % (outputFolder,analysisName)

# subPeakFile = '%smacsEnriched/UCH2_1_peaks.bed' % (projectFolder)

# crcCmd = 'python CRC2.py -s %s -e %s -b %s -g HG19 -o %s -n %s --promoter TRUE' % (subPeakFile,enhancerFile,bamFile,outputFolder,analysisName)



# bashFile = open(bashFileName,'w')

# bashFile.write('#!/usr/bin/bash\n')

# bashFile.write('cd /ark/home/af661/src/coreTFnetwork/\n')


# bashFile.write(crcCmd)

# bashFile.close()




#==========================================================================
#=========================CALLING META CRC=================================
#==========================================================================


#for meta 
#cd /ark/home/af661/src/medullo/
#python CRC2_fimo.py -e /grail/projects/medullo_final/meta_rose/shh/shh_AllEnhancers.table.targetGenes.txt -b /grail/projects/medullo_final/bam/medullo_shh.bam -g hg19 -o /grail/projects/medullo_final/network/meta_fimo/shh/ -n shh -x 50 -l 200 --promoter True &

#==========================================================================
#===========================FILTERING EDGE LISTS===========================
#==========================================================================


# enrichedCliqueFile = '%snetwork/UCH2_ENRICHED_CLIQUE_FACTORS.txt' %(projectFolder)


# degreeFile = '%snetwork/UCH2_DEGREE_TABLE.txt' % (projectFolder)


# edgeFile = '%snetwork/UCH2_EDGE_LIST.txt' % (projectFolder)
# newEdgeFile = '%snetwork/UCH2_EDGE_LIST_FILTERED.txt' % (projectFolder)


# enrichedCliqueTable = utils.parseTable(enrichedCliqueFile,'\t')

# degreeTable = utils.parseTable(degreeFile,'\t')

# #crcTFs = [line[0] for line in enrichedCliqueTable if float(line[1]) > 0.2]


# crcTFs = [line[0] for line in degreeTable[1:] if int(line[3]) > 40]

# print(crcTFs)





# edgeTable = utils.parseTable(edgeFile,'\t')
# newEdgeTable = [['FROM','TO']]

# for line in edgeTable[1:]:
#     if crcTFs.count(line[0]) > 0 and crcTFs.count(line[1]) > 0:
#         newEdgeTable.append(line)

# utils.unParseTable(newEdgeTable,newEdgeFile,'\t')


#==========================================================================
#====================ADDITIONAL PIPELINE ANALYSIS==========================
#==========================================================================


