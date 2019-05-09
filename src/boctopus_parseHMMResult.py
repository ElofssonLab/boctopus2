import os, sys, string


reverseAATable = { 'ALA' : 'A',   'CYS' : 'C',   'ASP' : 'D',   'GLU' : 'E',
                   'PHE' : 'F',   'GLY' : 'G',   'HIS' : 'H',   'ILE' : 'I',
                   'LYS' : 'K',   'MET' : 'M',   'ASN' : 'N',   'PRO' : 'P',
                   'GLN' : 'Q',   'ARG' : 'R',   'SER' : 'S',   'THR' : 'T',
                   'VAL' : 'V',   'TRP' : 'W',   'TYR' : 'Y',   'LEU' : 'L'}


def getObservations(oList, Home, outHome, dataHome):
        f = open(dataHome + dataFile, "r")
        lines = f.readlines()
        f.close()

        for line in lines:
                line = line.strip()
                line = string.split(line)
		#print line
		index = 4
                if line[0] in oList.keys():
                        if line[index] == "U":
                                #pList[line[0]].append("0")
                                oList[line[0]].append("M")
                        if line[index] == "D":
                                #pList[line[0]].append("0")
                                oList[line[0]].append("M")
                        if line[index] == "X":
                                #pList[line[0]].append("0")
                                oList[line[0]].append("M")
                        if line[index] == "i":
                                #pList[line[0]].append("1")
	                        oList[line[0]].append("i")
	                if line[index] == "o":
                                #pList[line[0]].append("0")			
                                oList[line[0]].append("o")

        return


def getObsSeq(Home, outHome, dataHome):
        f = open(dataHome + dataFile, "r")
        lines = f.readlines()
        f.close()

	fastaList = {}
	pList = {}
	oList = {}

        for line in lines:
                line = line.strip()
                line = string.split(line)
		proteinName = line[0]

		pList[proteinName] = []
		oList[proteinName] = []
		fastaList[proteinName] = [] 

        for line in lines:
                line = line.strip()
                line = string.split(line)
                fastaList[line[0]].append(reverseAATable[line[1]])

        return oList, pList, fastaList


def parsefile(filename, proteinName, Home, outHome, dataHome):#{{{
	f     = open(filename, "r")
	lines = f.readlines()
	f.close()
	
	flag = 0
	path = 0
	pred = []
	statePath = []
	pName = proteinName

	statefile = outHome + proteinName + "_stateName.txt"
	f = open(statefile, "w")

	scorefile = outHome + proteinName + "_score.txt"
        fscore = open(scorefile, "w")

	for line in lines:
                checkline = line.strip()
		if line.startswith("<lab"):
			flag = 1
		if line.startswith("</lab"):
			flag = 0

		if line.startswith("<path"):
			path = 1
		if line.startswith("</path"):
			path = 0
		if line.startswith("<pure_seq_name_a>"):
			line = line.strip()
			pName = (line.split(">")[1]).split(".")[0].split("_")[0]
			print "protein Name", pName

		if flag:
			line = line.strip()
			line = line.split(">")
			if len(line) > 1:
				data = line[1]
				if len(data) > 1:
					data = data.split("<")[0]
                                        if data in ['o', 'i', 'I', 'O', 'P', 'X', "M"]:
						#print data
						pred.append(data)

                if checkline.startswith("<normalizedLogLikelihood>"):
                        fscore.write("normalizedLogLikelihood " + checkline.split(">")[1].split("<")[0] + "\n")
                if checkline.startswith("<logodds>"):
                        fscore.write("logodds " + checkline.split(">")[1].split("<")[0]+ "\n")
		
		if path:		
                        line = line.strip()
                        line = line.split(">")
                        if len(line) > 1:
                                data = line[1]
				if len(data) > 1:
	                                data = data.split("<")[0]
        	                        #print data
                	                statePath.append(data)
		
        fscore.close()

	mypredicted_topology = "".join(pred)
	
	mypredicted_topology_newlabels = mypredicted_topology.replace("i","p").replace("o", "L")
	
	print "predicted topology", mypredicted_topology_newlabels
	print len(pred), len(statePath)

	f.write(">" + pName + "\n")
	for i in range(0, len(mypredicted_topology_newlabels)):
		f.write(mypredicted_topology_newlabels[i])
	f.close()

	os.system("cp " + statefile + " " + outHome + "/" + pName + "_topologies.txt")

	return statePath

#}}}
def getACC(oList, pList, filename, proteinName):#{{{
        totalCorrect = 0
        totalRes     = 0

	iaso = 0 
	oasi = 0 
	dasi = 0
	daso = 0
	oasd = 0
	iasd = 0
	dasd = 0
	
        uasi = 0
        uaso = 0
        oasu = 0
        iasu = 0
        uasu = 0

	iasx = 0
	xasi = 0
	dasx = 0
	xasd = 0
	oasx = 0
	xaso = 0
	uasx = 0
	xasu = 0

        for protein in oList:
                correct = 0
		if len(oList[protein]) == 0 or  len(pList[protein]) == 0:
			continue
		
                print protein, len(oList[protein]), len(pList[protein])
		print "".join(oList[protein])
		print "".join(pList[protein])

                for i in range(0, len(oList[protein])):
			obs  = oList[protein][i]
			pred = pList[protein][i]
			#print i, obs, pred
			
                        if obs == pred:
                                correct += 1
                                totalCorrect += 1
				if obs == "U":
					uasu += 1
				if obs == "D":
					dasd += 1
			if obs == "i" and pred == "o":
				iaso += 1
                        if obs == "o" and pred == "i":
				oasi += 1

                        if obs == "D" and pred == "i":
                                dasi += 1
                        if obs == "D" and pred == "o":
                                daso += 1
                        if obs == "i" and pred == "D":
                                iasd += 1
                        if obs == "o" and pred == "D":
                                oasd += 1

                        if obs == "U" and pred == "i":
				uasi += 1
                        if obs == "U" and pred == "o":
				uaso += 1
                        if obs == "i" and pred == "U":
				iasu += 1
                        if obs == "o" and pred == "U":
				oasu += 1

                        if obs == "i" and pred == "X":
                                iasx += 1
                        if obs == "X" and pred == "i":
                                xasi += 1
                        if obs == "o" and pred == "X":
                                oasx += 1
                        if obs == "X" and pred == "o":
                                xaso += 1
                        if obs == "U" and pred == "X":
                                uasx += 1
                        if obs == "X" and pred == "U":
                                xasu += 1
                        if obs == "D" and pred == "X":
                                dasx += 1
                        if obs == "X" and pred == "D":
                                xasd += 1

		print "iaso, oasi, dasi, daso, oasd, iasd, uasi, uaso, oasu, iasu: MasM"
		print iaso, oasi, dasi, daso, oasd, iasd, uasi, uaso, oasu, iasu, ":", dasd+uasu
		print "iasx, xasi, dasx, xasd, oasx, xaso, uasx, xasu"
                print iasx, xasi, dasx, xasd, oasx, xaso, uasx, xasu

                #print "Result", protein, correct, len(oList[protein]), "%.2f"%float((100.0*correct)/len(oList[protein]))
                totalRes += len(oList[protein])

        print "----"
        print "Result", proteinName, totalCorrect, totalRes, "%.2f"%float((100.0*totalCorrect)/totalRes)
	print

        return
#}}}

def getProteinName(inputFile):#{{{
	f = open(inputFile, "r")
	lines = f.readlines()
	f.close()

	for line in lines:
		line        = line.strip()
		prFile      = line.split("/")[-1]
		proteinName = prFile#.split(".")[0].split("_")[0]

	return proteinName
 #}}}

def start(filename, inputFile, path_seqfile):
        Home = path_seqfile + "/"
        svmoutHome = path_seqfile + "/svmoutput/"
        outHome    = path_seqfile + "/output/"
        tmpHome    = path_seqfile + "/tmp/"
        dataHome = Home
	
	proteinName = getProteinName(inputFile)
	print "proteinName ", proteinName

	statepath = parsefile(filename, proteinName, Home, outHome, dataHome)

	f = open(outHome + proteinName + "_statepath.txt", "w")
	for state in statepath:
		f.write(str(state) + " ")
	f.close()

	return

#start(sys.argv[0], sys.argv[1])

