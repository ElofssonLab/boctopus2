import os, sys, string
import re


strandEncoding = {'1':'D', '0':'X', '-1':'U'}
stateEncoding = {'1':'N', '0':'M'}

reverseAATable = { 'ALA' : 'A',   'CYS' : 'C',   'ASP' : 'D',   'GLU' : 'E',
                   'PHE' : 'F',   'GLY' : 'G',   'HIS' : 'H',   'ILE' : 'I',
                   'LYS' : 'K',   'MET' : 'M',   'ASN' : 'N',   'PRO' : 'P',
                   'GLN' : 'Q',   'ARG' : 'R',   'SER' : 'S',   'THR' : 'T',
                   'VAL' : 'V',   'TRP' : 'W',   'TYR' : 'Y',   'LEU' : 'L'}


def writeHeader(filepath, pName, segment, predictions, prob, probKey, fastaSeq):
        filename = pName + "." + segment + ".prf.txt"
        f = open(filepath + filename, "w")
        f.write("Sequence: " + pName +".fa" + "\n")
        f.write("NR of aligned sequences: 100" + "\n")
        f.write("Length of query sequence: " + str(len(prob)) + "\n")
        f.write("\n\n")
        f.write("START 1" + "\n")
        f.write("ALPHABET:\t")

        #print probKey, len(probKey[0])
        for i in range(0, len(probKey[1])):
                f.write(str(probKey[1][i])+"\t")
        f.write("-\t<SPACE>\t<LABEL>\t<QUERY>\n")

        print pName

        k = 0
        for i in range(0, len(fastaSeq)):
                f.write("COL\t" + str(i+1) + ":\t")
                for j in range(0, len(prob[k])):
                        f.write(str("%.2f"%(float(prob[k][j])))+"\t")
                f.write("0.00\t0.00\t.\t" + fastaSeq[i][0]+"\n")
                k += 1

        f.write("END 1")
        f.close()

        return filename
	#COL    1:       0.091   0.355   0.918   0.182   0.00   0.00   .       L


def readFile(filepath, filename):
	print "inside readfile", filename

	f     = open(filepath + filename, "r")
	lines = f.readlines()
	f.close()

	scount = 0
	ecount = 0
	classPred = 0
	probFlag  = 0
	probKey   = []

	data  = []
	result= []
	prob  = []
	for line in lines:
		line = line.strip()
		line = string.split(line)

		if len(line) == 2:
			flag = line[1]
			flag = flag[1:-1]
			if flag == "Start":
				scount += 1
			if flag == "End":
				ecount += 1
			if scount == 2:
				classPred = 1

			if flag == "Prob":
				probFlag = 1
				#print probFlag
	
		if classPred == 1 and ecount == 1 and len(line) == 1:
			data.append(line[0])

		if probFlag == 1:
			if len(line) == 2:
				#print probFlag, line, len(line)
				probKey.append(line)
			if len(line) == 3:
				prob.append(line[1:])

	for i in range(1, len(data)):
		if i%2 != 0:
			result.append(data[i])

	return result, prob, probKey


def getList(filepath, filename, fastaSeq, segment):
	proteinName = filename.split("_")[0]

	predictions, prob, probKey = readFile(filepath, filename)

	pred2prfFile = writeHeader(filepath, proteinName, segment, predictions, prob, probKey, fastaSeq)

        print filename, len(predictions), len(fastaSeq)

	return pred2prfFile


def start(filepath, filename, vSeq, segment):
	pred2prfFile = getList(filepath, filename, vSeq, segment)         ## returns predictions and probability of class
	
	return pred2prfFile
