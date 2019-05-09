#!/usr/bin/env python
import os, sys, string
import boctopus_makeplot

## 1. to obtain data from the ALL dataset file.
## 2. to obtain predinction outcomes of svms from prf files
## 3. to obtain predictions from the hmm model
## 4. to plot 1,2,3


def readdataset(pDict, filename, dataHome, outHome):
	f     = open(dataHome + filename, "r")
	lines = f.readlines()
	f.close()

	for data in lines:
		data = data.strip()
		data = string.split(data)
		#print data[:7]
		protein = data[0]
		if protein in pDict.keys():
			pDict[protein].append(data[:7])

	for pr in pDict:
		print pr, len(pDict[pr])
	return


def readprfFiles(pDict, pid, line):
	f    = open(line, "r")
        data = f.readlines()
        f.close()

	start = 0
	print "pid", pid

	labels = []
	profiles = []

        for prefs in data:
	        prefs = prefs.strip()
		prefs = string.split(prefs)

		if len(prefs) == 0:
			continue

		if prefs[0] == "ALPHABET:":
	                start = 1
			for label in prefs[1:-4]:
				labels.append(label)

		if start == 1:
			##print prefs[0]
			if prefs[0] == "COL":
				#print pid, prefs[2:-4]
				profiles = prefs[2:-4]
				profiles.append(prefs[-1])
				pDict[pid].append(profiles)

	#print labels	
	return labels


def getProteinDict(data):
	pDict = {}

	for pName in data:
		pName     = pName.strip()
		pName     = string.split(pName)
		prfName   = pName[0].split("/")[-1]
		proteinID = prfName#.split(".")[0]
		proteinID = proteinID#.split("_")[0]

		print proteinID
		pDict[proteinID] = []

	return pDict


def readPredTopo(proDict_predTopo, proteinid, filename, dataHome, outHome):
        f = open(outHome + proteinid + "_stateName.txt", "r")
        lines = f.readlines()
        f.close()

	data = []
	for line in lines:
		line = line.strip()
		if line.startswith(">"):
			continue
		#line = string.split(line)
		proDict_predTopo[proteinid].append(line)

	return
	

def start(fileList, path_seqfile):
        Home = path_seqfile + "/"
        svmoutHome = path_seqfile + "/svmoutput/"
        outHome    = path_seqfile + "/output/"
        tmpHome    = path_seqfile + "/tmp/"
        dataHome = Home

	f     = open(fileList, "r")
	lines = f.readlines()
	f.close()

	proDict_dataset  = getProteinDict(lines)
        proDict_profiles = getProteinDict(lines)
	proDict_predTopo = getProteinDict(lines)
	
        #readdataset(proDict_dataset, dataFile)

	for pr in proDict_profiles:
		print "In dict", pr

	for line in lines:
		line = line.strip()

                pName = string.split(line)
                prfName   = pName[0].split("/")[-1]
                proteinID = prfName#.split(".")[0]
                proteinID = proteinID#.split("_")[0]
	
		labels = readprfFiles(proDict_profiles, proteinID, line)
		print proteinID, labels

		readPredTopo(proDict_predTopo, proteinID, line, dataHome, outHome)

	print "Generating plots"
	boctopus_makeplot.start(proDict_dataset, proDict_profiles, proDict_predTopo, labels, proteinID, path_seqfile)

	return


##fileList = sys.argv[1]
#fileList = "prFile.txt"
#start(fileList)
