import sys
import numpy as np

#
#	mergeBed.py:
#		- Combines multiple bed files, creating an output with non-overlapping regions
#		- Preserves all metadata present in columns 4+
#		- When regions overlap, precedence is given to the smaller region
#
#	Written by: Zach Stephens (04/19/2016)
#

if len(sys.argv) == 1:
	print '\npython mergeBed.py outFile.bed inFile1.bed [infile2.bed...]\n'
	exit(1)

outFile = sys.argv[1]
inFiles = sys.argv[2:]

# lexicographic sorting is annoying
chr_ind = {'chr1':1,   'chr2':2,   'chr3':3,   'chr4':4,   'chr5':5,   'chr6':6,   'chr7':7,   'chr8':8,   'chr9':9,   'chr10':10,
           'chr11':11, 'chr12':12, 'chr13':13, 'chr14':14, 'chr15':15, 'chr16':16, 'chr17':17, 'chr18':18, 'chr19':19, 'chr20':20,
           'chr21':21, 'chr22':22, 'chr23':23, 'chr24':24, 'chr25':25, 'chr26':26, 'chr27':27, 'chr28':28, 'chr29':29, 'chr30':30,
           'chr31':31, 'chr32':32, 'chr33':33, 'chr34':34, 'chr35':35, 'chr36':36, 'chr37':37, 'chr38':38, 'chr39':39, 'chr40':40,
           'chrX':41,  'chrY':42,  'chrM':43,  'chrC':44}

max_pos_found    = {}
currentDataInd   = 2
big_dict_of_data = {}
ALL_REGIONS      = {}

#
#	READ IN REGION DATA FROM BED FILES
#
for fn in inFiles:
	f = open(fn,'r')
	for line in f:
		splt  = line.strip().split('\t')
		myChr = splt[0]
		coord = (int(splt[1]),int(splt[2]))

		myLen = coord[1] - coord[0]
		if myLen <= 0:
			print 'Skipping a line because region length was non-positive for some reason...'
			continue

		if myChr not in max_pos_found:
			max_pos_found[myChr]    = coord[1]
			big_dict_of_data[myChr] = {'':1}
			ALL_REGIONS[myChr]      = []

		if coord[1] > max_pos_found[myChr]:
			max_pos_found[myChr] = coord[1]

		if len(splt) >= 4:
			myDat = '\t'.join(splt[3:])
			if myDat in big_dict_of_data[myChr]:
				myDataInd = big_dict_of_data[myChr][myDat]
			else:
				big_dict_of_data[myChr][myDat] = currentDataInd
				myDataInd = currentDataInd
				currentDataInd += 1
		else:
			myDataInd = 1

		ALL_REGIONS[myChr].append((myLen,coord,myDataInd))
	f.close()

#
#	CREATE SOME BIG DATASTRUCTURES THAT WE'LL NEED
#
reverse_data_dict = {}
for k in sorted(big_dict_of_data.keys()):
	reverse_data_dict[k] = {}
	for k2 in big_dict_of_data[k].keys():
		reverse_data_dict[k][big_dict_of_data[k][k2]] = k2

for k in sorted(ALL_REGIONS.keys()):
	ALL_REGIONS[k].sort(reverse=True)

currentSortInd = len(chr_ind)+1
sortedChrList  = []
for k in ALL_REGIONS.keys():
	if k.lower() in chr_ind:
		sortedChrList.append([chr_ind[k.lower()],k])
	else:
		sortedChrList.append([currentSortInd,k])
		currentSortInd += 1
sortedChrList = [n[1] for n in sorted(sortedChrList)]

#
#	MERGE ALL BED DATA, ONE REF AT A TIME. SO I HEARD YOU LIKE GIANT ARRAYS???
#
f = open(outFile,'w')
for k in sortedChrList:
	# 0 means no region here, 1 means no data, otherwise value is index of data in reverse_data_dict
	datInds = np.zeros(max_pos_found[k]+100,dtype='<u4')

	for r in ALL_REGIONS[k]:
		for i in xrange(r[1][0],r[1][1]):
			datInds[i] = r[2]

	regStart = 0
	prevVal  = 0
	for i in xrange(len(datInds)):
		if datInds[i] > 0:
			myBedVal = datInds[i]
			if myBedVal != prevVal:
				if prevVal > 0:
					#bedRegion = (k,regStart,i,prevVal)
					#print bedRegion
					f.write(k + '\t' + str(regStart) + '\t' + str(i) + '\t' + reverse_data_dict[k][prevVal] + '\n')
				regStart = i
				prevVal  = myBedVal
		else:
			if prevVal > 0:
				#bedRegion = (k,regStart,i,prevVal)
				#print bedRegion
				f.write(k + '\t' + str(regStart) + '\t' + str(i) + '\t' + reverse_data_dict[k][prevVal] + '\n')
				prevVal = 0
f.close()

