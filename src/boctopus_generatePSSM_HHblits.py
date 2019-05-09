#!/usr/bin/env python
import os, sys, string
import filter_gappyseqs
import boctopus_neighborJoin_profile_pssm

#use hhblits to get homologoues sequences in psi format
#/afs/pdc.kth.se/home/s/shayat/shayat9/hhsuite/hhsuite-2.0.13-linux-x86_64/bin/hhblits -i 1kmpA_seq.txt -d /afs/pdc.kth.se/home/s/shayat/shayat9/hhsuite/database/nr20_12Aug11 -opsi pssm.psi

#convert psi to clustal format
#/afs/pdc.kth.se/home/s/shayat/shayat9/hhsuite/hhsearch_1.5.0/./reformat.pl pssm.psi pssm.clu

#use clustal format with a fake db to get pssm
#/afs/pdc.kth.se/home/a/arnee/MODULES/seqtools/blast-2.2.18/bin/blastpgp -i 1kmpA_seq.txt -B pssm.clu -j 1 -d ../../BLAST_NR_70/example.fasta -Q aa.pssm



def fillPssmData(filename, dataHolder):#{{{
        f     = open(filename, "r")
        lines = f.readlines()
        f.close()

        for line in lines:
                line = line.strip()
                line = string.split(line)

                dataHolder.append(line[1:])

        return
#}}}

def makeNeighborPssmData(pssmName, filename, aaHolder, neighborbasePath):#{{{
        f     = open(filename, "r")
        lines = f.readlines()
        f.close()

        f       = open(neighborbasePath + "tmp/" + pssmName + ".temp.filtered", "w")
        counter = 0

        for line in lines:
                line = line.strip()
                line = string.split(line)

                if len(line) == 44:
                        aaHolder.append(line[1:22])
                        f.write(str(counter))
                        for data in line[2:22]:
                                f.write(" " + data)
                        f.write("\n")
                        counter += 1

        f.close()

        boctopus_neighborJoin_profile_pssm.myneighbourMain(pssmName, pssmName + ".temp.filtered", neighborbasePath)

        return
#}}}

def getPSSM(pname, seq, Home, outHome, pssmHome, tmpHome, neighborbasePath):#{{{
        pssmData     = []
        pssmAAHolder = []
        pssm41_path  = ""

        pssm41_path = neighborbasePath + "pssm/" + pname + "_41.filtered"
 
	filename = pssmHome + pname + ".pssm.filtered"

	makeNeighborPssmData(pname, filename, pssmAAHolder, neighborbasePath)
        
	fillPssmData(pssm41_path, pssmData)

        return pssmData, pssmAAHolder, pssm41_path
#}}}

def parsefirstline(pname, reformatpath, Home, outHome, pssmHome, tmpHome, neighborbasePath):#{{{
        os.system("cp " + Home + pname + ".psi.filtered " + Home + "temp.psi")

        command = reformatpath + " %s/temp.psi %s/temp.clu"%(Home, Home)
        print command
        os.system(command)

        os.system("mv %s/temp.clu "%(Home) + pssmHome+pname + ".clu.filtered")

	if 0:
		f = open(tmpHome + pname + ".psi", "r")
		lines = f.readlines()
		f.close()
	
		emptystr = ""

		print pname
		for i in range(0, 21-len(pname)):
			emptystr += " "

		f = open(tmpHome + "temp" + ".psi", "w")
		for line in lines:
			line = line.strip()
			#print line
			## there is an extra - at the begining of the seqwuence. remove that.
			## limit the number of characters to compare in the pname
			if line.startswith(pname[:10]):
				line = string.split(line)
				#print line
				#print line[0] + emptystr+ line[1][1:]
				f.write(line[0] + emptystr + line[1][1:] + "\n")
			else:
				#print line
				f.write(line + "\n")

		f.close()

		command = hhsearchpath + " " + tmpHome + "temp.psi " + tmpHome + pname + ".clu"
		print command
		os.system(command)	

	return
#}}}

def formatcluFile(pname, seqs, seqFile, blastpgppath, fakedbpath, Home, outHome, pssmHome, tmpHome, neighborbasePath):
        tempEntry = []
        count = 0
        for aa in seqs:
                tempEntry.append([aa, count])
                count += 1

        if not os.path.exists(pssmHome + pname + ".pssm.filtered"):
                f     = open(pssmHome + pname + ".clu.filtered", "r")
                lines = f.readlines()
                f.close()

                flag = 0
                f = open(Home + "temp" + ".clu.filtered", "w")
                for line in lines:
                        line = line.strip()
                        if line.startswith("CLUSTAL"):
                                continue
                        if len(line) > 0:
                                flag = 1
                        if flag:
                                f.write(line + "\n")

                f.close()

                command = blastpgppath + " -i " + outHome+"/"+seqFile + " -B " + Home + "temp" + ".clu.filtered -j 1 -d " + fakedbpath + " -Q " + pssmHome+pname + ".pssm.filtered"
                print command
                os.system(command)

        if not os.path.exists(pssmHome + pname + ".filtered.pssmvals"):
                getpssmvalues(pname, seqs, tempEntry, Home, outHome, pssmHome, tmpHome, neighborbasePath)

        return


def getpssmvalues(pname, seqs, tempEntry, Home, outHome, pssmHome, tmpHome, neighborbasePath):
        f     = open(pssmHome+pname + ".pssm.filtered", "r")
        lines = f.readlines()
        f.close()

        result, tempIndex = [], 0
        for i in range(0,len(lines)):
                temp = string.split(lines[i])
                if len(temp) > 40:
                        if temp[1] != "X":
                                #print ">> ", temp, tempEntry, tempIndex
                                result.append([tempEntry[tempIndex][1]]+temp[2:22])
                        else:
                                tempValue = [0]*20
                                tempValue[aaIndex[tempEntry[tempIndex][0]]] = 1
                                result.append([tempEntry[tempIndex][1]]+tempValue)
                        tempIndex += 1

        if len(seqs) != len(result):
                print "Something wrong. len(seqs) != len(result)", len(seqs), len(result)
                print tempEntry
                print result
                sys.exit()
        else:
                print "Writing PSSM", pname
                f = open(pssmHome+pname + ".filtered.pssmvals", "w")
                for i in range(0, len(result)):
                        f.write(str(result[i][0]) + " ")
                        for j in range(1, 21):
                                f.write(str(result[i][j]) + " ")
                        f.write("\n")
                f.close()

        return


def _getpssmvalues(pname, seqs, tempEntry, Home, outHome, pssmHome, tmpHome, neighborbasePath):
        f     = open(pssmHome + pname + ".pssm", "r")
        lines = f.readlines()
        f.close()

        result, tempIndex = [], 0
        for i in range(0,len(lines)):
                temp = string.split(lines[i])
                if len(temp) > 40:
                        if temp[1] != "X":
                                #print ">> ", temp, tempEntry, tempIndex
                                result.append([tempEntry[tempIndex][1]]+temp[2:22])
                        else:
                                tempValue = [0]*20
                                tempValue[aaIndex[tempEntry[tempIndex][0]]] = 1
                                result.append([tempEntry[tempIndex][1]]+tempValue)
                        tempIndex += 1

        if len(seqs) != len(result):
                print "Something wrong. len(seqs) != len(result)", len(seqs), len(result)
                print tempEntry
                print result
                sys.exit()
        else:
                print "Writing PSSM", pname
                f = open(pssmHome + pname + ".pssmvals", "w")
                for i in range(0, len(result)):
                        f.write(str(result[i][0]) + " ")
                        for j in range(1, 21):
                                f.write(str(result[i][j]) + " ")
                        f.write("\n")
                f.close()

        return


def start(entry, mySequence, fakedbpath, dbpath, blastpgppath, hhsearchpath, hhblitspath, path_seqfile):#{{{
	#print "fakedbpath", fakedbpath
	#print "dbpath", dbpath
	#print "blastpgppath", blastpgppath
	#print "hhsearchpath", hhsearchpath
	#print "hhblitspath", hhblitspath

        Home = path_seqfile + "/"
        svmoutHome = path_seqfile + "/svmoutput/"
        outHome    = path_seqfile + "/output/"
        tmpHome    = path_seqfile + "/tmp/"
        pssmHome = path_seqfile + "/pssm/"
        neighborbasePath = Home

        seqFile = entry + ".seq"
        f = open(outHome + seqFile, "w")
        f.write(">" + entry + "\n")
        f.write(mySequence + "\n")
        f.close()

	if not os.path.exists(Home + entry + ".pssm"):	
		#command = hhblitspath + " -i " + outHome + seqFile + " -d " + dbpath + " -opsi " + tmpHome + entry + ".psi" 
		print
		command = hhblitspath + "/bin/hhblits -i " + outHome + seqFile + " -d " + dbpath + " -opsi " + Home + entry+".psi -oa3m " + Home + entry+".a3m -all -n 4 -maxfilt 100000 -realignoldhits -realign_max 5000 -Z 5000 -B 5000 -cpu 4 -e 1E-3 > " + tmpHome + entry + "_hhblits.stdout"
                print command
		os.system(command)

		filter_gappyseqs.start(Home + entry + ".psi")
		print
		
		if not os.path.exists(Home + entry + ".clu.filtered"):
			parsefirstline(entry, hhsearchpath, Home, outHome, pssmHome, tmpHome, neighborbasePath)

		formatcluFile(entry, mySequence, seqFile, blastpgppath, fakedbpath, Home, outHome, pssmHome, tmpHome, neighborbasePath)

		pssmData, pssmAAHolder, pssm21_path = getPSSM(entry, mySequence, Home, outHome, pssmHome, tmpHome, neighborbasePath)

	return pssmData, pssmAAHolder, pssm21_path
#}}}

#entry = "1P4T"
#mySequence = "EGASGFYVQADAAHAKASSSLGSAKGFSPRISAGYRINDLRFAVDYTRYKNYKAPSTDFKLYSIGASAIYDFDTQSPVKPYLGARLSLNRASVDLGGSDSFSQTSIGLGVLTGVSYAVTPNVDLDAGYRYNYIGKVNTVKNVRSGELSAGVRVKF"

entry       = sys.argv[1]
mySequence  = sys.argv[2]
fakedbpath  = sys.argv[3]
dbpath      = sys.argv[4]
blastpgppath= sys.argv[5]
hhsearchpath=sys.argv[6]
hhblitspath = sys.argv[7]

#print entry
#print mySequence
#sys.exit()

#start(entry, mySequence, fakedbpath, dbpath, blastpgppath, hhsearchpath, hhblitspath)

