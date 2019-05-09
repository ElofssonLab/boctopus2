## standard leave one out on the selected parameters.

import os, sys, string, random
import boctopus_writepred2prf

rundir = os.path.realpath(os.path.dirname(__file__))



threshold = 0.3


def write_Rscript():#{{{
	return#}}}


def write_profilefile(pName, iData, oData, inloopData, outloopData, fastaSeq, start, end, Home, svmoutPath, tmpHome):#{{{
        svmfile = pName + "_ioIOS.prf.txt"

        #print fastaSeq
        f = open(svmoutPath+svmfile, "w")
        f.write("Sequence: " + pName +".fa" + "\n")
        f.write("NR of aligned sequences: 100" + "\n")
        f.write("Length of query sequence: " + str(len(fastaSeq)) + "\n")
        f.write("\n\n")
        f.write("START 1" + "\n")
        f.write("ALPHABET:\t")

        f.write(str(iData[0])+"\t"+ str(oData[0])+"\t" + str(inloopData[0])+"\t"+str(outloopData[0])+"\tS\t")
        f.write("-\t<SPACE>\t<LABEL>\t<QUERY>\n")

        for i in range(1, len(iData)):
                f.write("COL\t" + str(i) + ":\t")
                if i >= start and i <= end:
                        f.write(str(iData[i]) +"\t"+ str(oData[i])+"\t"+ str(inloopData[i]) +"\t"+ str(outloopData[i])+"\t0.0")
                else:
                        f.write(str(0.00) +"\t"+ str(0.00)+"\t"+ str(1.0) +"\t"+ str(0.0)+"\t0.0")
                f.write("\t0.00\t0.00\t.\t" + fastaSeq[i-1]+"\n")

        f.write("END 1")
        f.close()

        return
#}}}

def getrunning_average(pname, data1, data2, pflag):
        mydata = []
        for i in range(1, len(data1)):
                mydata.append(float(data1[i]) + float(data2[i]))

        ws = 100

        x = []
        y = []

        for i in range(0, len(mydata)):
                start = i- ws/2
                end   = i + ws/2

                if start < 0:
                        start = 0

                if end > len(mydata):
                        end = len(mydata)

                if len(mydata[start:i]) > 0:
                        avg1 = sum(mydata[start:i])/len(mydata[start:i])
                else:
                        avg1 = 0.0

                if len(mydata[i:end]) > 0:
                        avg2 = sum(mydata[i:end])/len(mydata[i:end])
                else:
                        avg2 = 0.0

                if pflag:
                        print start, i, end,  avg1+avg2

                x.append(i)
                y.append(avg1+avg2)

        return x, y


def find_barrelregio(avgx, avgy):
        flag = 0
        start= 0
        startflag = 0
        end  = len(avgx)

        for i in range(0, len(avgx)):
                #print avgx[i], avgy[i]
                if avgy[i] > threshold and startflag == 0:
                        startflag  = 1
                        start = i
                if avgy[i] > threshold:
                        flag = 1
                if flag and avgy[i] < threshold:
                        end = i
                        #print "end", end
                        flag = 0
                if flag and avgy[i] > threshold:
                        end = len(avgx)
                        #print "**", flag

        return start, end


def readprfFiles(filename, pDict, pidcid):
	f    = open(filename, "r")
        data = f.readlines()
        f.close()
        
	start = 0

	labels   = []
	profiles = []

        for prefs in data:
	        prefs = prefs.strip()
		prefs = string.split(prefs)
                
		if len(prefs) == 0:
			continue

		if prefs[0] == "ALPHABET:":
	                start = 1
			for label in prefs[1:-3]:
				labels.append(label)

		if start == 1:
			if prefs[0] == "COL":
				#print pidcid, prefs, prefs[2:-3], prefs[-1]
				profiles = prefs[2:-3]
				profiles.append(prefs[-1])
				pDict[pidcid].append(profiles)
    
	return labels


def get_labeledTestFiles(testX, windowSize, testfilename, Home, svmoutPath, tmpHome):#{{{
	windowSize = int(windowSize)
	factor     = 10.0

        print "Remember how you are dividing the trainX parameter", factor, "windowSize", windowSize

        f = open(tmpHome + testfilename, "w")
        for i in range(0,len(testX)):
                for j1 in range(1, windowSize+1):
                        zero, one, flagNA = 0, 0, 0
                        for j2 in range(20*(j1-1)+0,20*(j1-1)+20):
                                if testX[i][j2] != "NA":
                                        f.write(str(float(testX[i][j2])/factor)+" ")
                                else:
                                        f.write("0 ")
                                        flagNA  = 1
                f.write("\n")

        f.close()

	return
#}}}

def actualPred(proteinDict, ws_cytosolic, ws_extracellular, ws_lipidfacing, ws_porefacing, outHome, prefix, seq, Home, svmoutPath, tmpHome):#{{{
	proteins = proteinDict.keys()

	pDict_in = {}
	pDict_ot = {}
	pDict_top= {}
	pDict_inloop = {}
	pDict_otloop = {}

	for i in range(0, len(proteins)):
		pTest = proteins[i]

                testProtein = []
                testProtein.append(pTest)

		print "SVM-stage", testProtein

		testX  = []

		for protein in testProtein:
			for data in proteinDict[protein]:
				testX.append(data[1:])
			if len(proteinDict[protein]) == len(testX):
				print "Running SVMs for", protein, len(proteinDict[protein]), len(testX)
			else:
				print "boctopus_predict length error ", protein, len(proteinDict[protein]), len(testX)
				sys.exit()
		
		pDict_in[pTest] = []
	    	pDict_ot[pTest] = []
    		pDict_otloop[pTest] = []
    		pDict_inloop[pTest] = []

		#print protein, pTest

                testfilename = pTest + "_cyto.dat"
               	get_labeledTestFiles(testX, ws_cytosolic, testfilename, Home, svmoutPath, tmpHome)
               	#os.system("cp " + tmpHome + testfilename + " " + tmpHome + "test.dat")
               	outputPred = pTest + "_" + str(ws_cytosolic) + "_cytoloop_svmradialResult.txt"
               	rcommand1 = prefix + "R --vanilla --slave < " + "%s/traninedsvmRadial_cyto.R --args "%(rundir) + tmpHome+testfilename + " > " + svmoutPath + outputPred
		print rcommand1
		os.system(rcommand1)
		print "reading prf file Innerloop"
		Ifile =boctopus_writepred2prf.start(svmoutPath, outputPred, seq, "I")
                readprfFiles(svmoutPath+Ifile, pDict_inloop, pTest)
		print

	
               	testfilename = pTest + "_extra.dat"
               	get_labeledTestFiles(testX, ws_extracellular, testfilename,  Home, svmoutPath, tmpHome)
               	#os.system("cp " + tmpHome + testfilename + " " + tmpHome + "test.dat")
		outputPred = pTest + "_" + str(ws_extracellular) + "_extraloop_svmradialResult.txt"
               	rcommand1 = prefix + "R --vanilla --slave < " + "%s/traninedsvmRadial_extra.R --args "%(rundir) + tmpHome+testfilename + " > " + svmoutPath + outputPred
		print rcommand1
		os.system(rcommand1)
		print "reading prf file Outloop"
		Ofile = boctopus_writepred2prf.start(svmoutPath, outputPred, seq, "O")
                readprfFiles(svmoutPath+Ofile, pDict_otloop, pTest)
		print


                testfilename = pTest + "_lipid.dat"
                get_labeledTestFiles(testX, ws_lipidfacing, testfilename,  Home, svmoutPath, tmpHome)
                #os.system("cp " + tmpHome + testfilename + " " + tmpHome + "test.dat")
                outputPred = pTest + "_" + str(ws_lipidfacing) + "_lipidfacing_svmradialResult.txt"
                rcommand1 = prefix + "R --vanilla --slave < " + "%s/traninedsvmRadial_lipid.R --args "%(rundir) + tmpHome+testfilename + " > " + svmoutPath + outputPred
		print rcommand1
                os.system(rcommand1)
		print "reading prf file lipid-facing"
                ofile = boctopus_writepred2prf.start(svmoutPath, outputPred, seq, "o")
                readprfFiles(svmoutPath+ofile, pDict_ot, pTest)
		print
		

                testfilename = pTest + "_pore.dat"
                get_labeledTestFiles(testX, ws_porefacing, testfilename,  Home, svmoutPath, tmpHome)
                #os.system("cp " + tmpHome + testfilename + " " + tmpHome + "test.dat")
                outputPred = pTest + "_" + str(ws_porefacing) + "_porefacing_svmradialResult.txt"
                rcommand1 = prefix + "R --vanilla --slave < " + "%s/traninedsvmRadial_pore.R --args "%(rundir) + tmpHome+testfilename + " > " + svmoutPath + outputPred
		print rcommand1
                os.system(rcommand1)
		print "reading prf file pore-facing"
                ifile = boctopus_writepred2prf.start(svmoutPath, outputPred, seq, "i")
                readprfFiles(svmoutPath+ifile, pDict_in, pTest)
		print


		print "---------------"
		print ofile, len(pDict_ot[pTest])
		print ifile, len(pDict_in[pTest])
		print Ifile, len(pDict_inloop[pTest])
		print Ofile, len(pDict_otloop[pTest])
		print "---------------"

		pflag = 0
		print
		for pname in pDict_in:
			#print "pname", pname, len(pDict_in[pname]), len(pDict_ot[pname])

		        iData = ["i"]
		        oData = ["o"]
		        inloopData  = ["I"]
		        outloopData = ["O"]

	        for i in range(0, len(pDict_in[pname])):
	            iData.append(float(pDict_in[pname][i][1]))
        	    oData.append(float(pDict_ot[pname][i][1]))
	            inloopData.append(float(pDict_inloop[pname][i][0]))
        	    outloopData.append(float(pDict_otloop[pname][i][1]))
            
	        avgx, avgy = getrunning_average(pTest, iData, oData, pflag)
		start, end = find_barrelregio(avgx, avgy)    

    		print "Barrel region: start, end, protein_length - ", start, end, len(seq)

    		write_profilefile(pname, iData, oData, inloopData, outloopData, seq, start, end, Home, svmoutPath, tmpHome)

	return
	#}}}

def readPSSM21Files_3class(filename, proteinName):
	print filename

        f     = open(filename, "r")
        lines = f.readlines()
        f.close()

	proteinDict = {}
	proteinDict[proteinName] = []

        for line in lines:
                line = line.strip()
                line = string.split(line)
		
	        proteinDict[proteinName].append(line)
			
        return proteinDict


def startPredicting(outHome, ws_cytosolic, ws_extracellular, ws_lipidfacing, ws_porefacing, filename, proteinName, prefix, validaaSeq, path_seqfile):#{{{
        Home = path_seqfile + "/"
        svmoutHome = path_seqfile + "/svmoutput/"
        outHome    = path_seqfile + "/output/"
        tmpHome    = path_seqfile + "/tmp/"
        svmoutPath = svmoutHome
        print "startPredicting in boctopus_predictsvm.py"

        proteinDict = readPSSM21Files_3class(filename, proteinName)

        #for protein in proteinDict.keys():
        #        print protein, len(proteinDict[protein])
        actualPred(proteinDict, ws_cytosolic, ws_extracellular, ws_lipidfacing, ws_porefacing, outHome, prefix, validaaSeq, Home, svmoutPath, tmpHome)

        return
#}}}
