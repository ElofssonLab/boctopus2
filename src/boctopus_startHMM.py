import os, sys, string, time, random

import boctopus_neighborJoin_profile_pssm
import boctopus_generatePSSM_HHblits
import boctopus_predictsvm
import boctopus_callHMM

# argument parsing
numArgv = len(sys.argv)

if numArgv != 15:
    print >> sys.stderr, "Syntax error for %s!"%(sys.argv[0])
    sys.exit(1)

fastainputfile = sys.argv[1]
blastpath      = sys.argv[2]

modHome        = sys.argv[3]
hmmfilename    = sys.argv[4]
ws_cytosolic   = sys.argv[5]
ws_extracellular = sys.argv[6]
ws_lipidfacing   = sys.argv[7]
ws_porefacing    = sys.argv[8]
rpath            = sys.argv[9]

fakedbpath   = sys.argv[10]
dbpath       = sys.argv[11]
blastpgppath = sys.argv[12]
hhsearchpath =  sys.argv[13]
hhblitspath  =  sys.argv[14]

path_seqfile = os.path.dirname(os.path.realpath(fastainputfile))



outHome  = path_seqfile + "/output/"
svmoutPath = path_seqfile +"/svmoutput/"
tmpHome  = path_seqfile + "/tmp/"
pssmHome = path_seqfile + "/pssm/"

if not os.path.exists(outHome):
    os.makedirs(outHome)
if not os.path.exists(svmoutPath):
    os.makedirs(svmoutPath)
if not os.path.exists(tmpHome):
    os.makedirs(tmpHome)
if not os.path.exists(pssmHome):
    os.makedirs(pssmHome)

neighborbasePath = path_seqfile
neighborbaseHome = path_seqfile

aaTable = { 'A' : 'ALA',   'C' : 'CYS',   'D' : 'ASP',   'E' : 'GLU',
            'F' : 'PHE',   'G' : 'GLY',   'H' : 'HIS',   'I' : 'ILE',
            'K' : 'LYS',   'M' : 'MET',   'N' : 'ASN',   'P' : 'PRO',
            'Q' : 'GLN',   'R' : 'ARG',   'S' : 'SER',   'T' : 'THR',
            'V' : 'VAL',   'W' : 'TRP',   'Y' : 'TYR',   'L' : 'LEU'}


def makeNeighborPssmData(pssmName, filename, aaHolder, ivalue, cddflag):#{{{
    f     = open(filename, "r")
    lines = f.readlines()
    f.close()

    f       = open(neighborbasePath + "tmp/" + pssmName + ".temp", "w")
    counter = 0

    iterations  = 3

    checklength = 22

    if iterations > 1:        
        checklength = 44
    if cddflag == 0:
        checklength = 44
 
    for line in lines:
        line = line.strip()
        line = string.split(line)

        if len(line) == checklength:
            aaHolder.append(line[1:22])
            f.write(str(counter))
            for data in line[2:22]:
                    f.write(" " + data)
            f.write("\n")
            counter += 1

        f.close()

        boctopus_neighborJoin_profile_pssm.myneighbourMain(pssmName, neighborbaseHome)

    return
#}}}

def fillPssmData(filename, dataHolder):#{{{
    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    for line in lines:
        line = line.strip()
        line = string.split(line)

        dataHolder.append(line[1:])

    return
#}}}

def generatePSSM(name, seq, blastpath, evalue, ivalue, cddflag, fakedbpath, dbpath, blastpgppath, hhsearchpath, hhblitspath):#{{{
    ## hhblits with gaps removed
    print "Calling boctopus_generatePSSM_HHblits"

    return boctopus_generatePSSM_HHblits.start(name, seq, fakedbpath, dbpath, blastpgppath, hhsearchpath, hhblitspath, path_seqfile)
#}}}

def getPSSM(proteinName, seq, blastpath, evalue, ivalue, cddflag, fakedbpath, dbpath, blastpgppath, hhsearchpath, hhblitspath):#{{{
    pssmData     = []
    pssmAAHolder = []
    pssm21_path  = ""
        
    pssmmat_path = path_seqfile + "/tmp/" + proteinName + "_temp.mat"

    pssmData, pssmAAHolder, pssm21_path = generatePSSM(proteinName, seq, blastpath, evalue, ivalue, cddflag, fakedbpath, dbpath, blastpgppath, hhsearchpath, hhblitspath)

    print "PSSM generated at", pssm21_path
    
    chkPath = ""

    return pssmData, pssmAAHolder, pssm21_path, chkPath
#}}}

def readFasta(fastaname, seq, blastpath, evalue, ivalue, cddflag, fakedbpath, dbpath, blastpgppath, hhsearchpath, hhblitspath):#{{{
    pssmData, pssmAAHolder, pssm21_path, chkPath = getPSSM(fastaname, seq, blastpath, evalue, ivalue, cddflag, fakedbpath, dbpath, blastpgppath, hhsearchpath, hhblitspath)

    dictData = []

    if pssm21_path == "":
        return dictData, pssm21_path

    for i in range(0, len(seq)):
        if pssmAAHolder[i][0] in aaTable.keys():
            if pssmAAHolder[i][0] == seq[i]:
                    dictData.append([seq[i], str(i+1), pssmData[i]])
                    ## try to keep the pssm as the last value
            else:
                print ">> PSSM error EXITING ", aaTable[pssmAAHolder[i][0]], seq[i], fastaname
                return -1, -1
        else:
            #print "invalid AA", pssmData[i][:20], pssmAAHolder[i][0], seq[i], "EXCLUDED FROM DATASET"
            dictData.append([seq[i], str(i+1), pssmData[i]])
            ## try to keep the pssm as the last value

    return dictData, pssm21_path, chkPath
#}}}

def getvalidSeq(seqs):#{{{
    temp = ""
    for i in range(0, len(seqs)):
        if seqs[i] in aaTable.keys():
            temp+= seqs[i]

    return temp
#}}}

def getsequences(lines):#{{{
    tmpfasta   = []    ## name obtained from the user
    tmpSeq     = []
    fastaindex = []    ## name given by us

    sequence = ""
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            tmpfasta.append(line)
            #tempname = line[1:] + "_" + str(random.randrange(1000)) + "-"+str(time.time())
            #tempname = tempname[:15]
            tempname = "query"

            fastaindex.append(str(tempname))
            if len(sequence) > 0:
                tmpSeq.append(sequence)
                sequence = ""
        else:
            sequence += line

    tmpSeq.append(sequence)

    return tmpfasta, tmpSeq, fastaindex
    #}}}

def readList(inputfile, blastpath, modHome, hmmfilename, ws_cytosolic, ws_extracellular, ws_lipidfacing, ws_porefacing, rpath, fakedbpath, dbpath, blastpgppath, hhsearchpath, hhblitspath):#{{{
    cddflag = 0
    evalue  = 1
    ivalue  = 3

    print
    print "Input parameters" 
    print inputfile
    print rpath
    print fakedbpath
    print dbpath
    print blastpgppath
    print hhsearchpath
    print hhblitspath
    
    f     = open(inputfile, "r")
    lines = f.readlines()
    f.close()
    
    print
    print "Reading list to process"
    fastaList, seqList, fastaName = getsequences(lines)
    print "Files to process", fastaName

    prfSuffix = "_ioIOS.prf.txt"

    for i in range(0, len(fastaList)):
        validaaSeq  = getvalidSeq(seqList[i])
        proteinName = fastaName[i]
        pName       = proteinName #.split("_")[0] + proteinName.split("_")[1]

        print
        print "Generating PSSM ..."
        data, pssm21_path, chkPath = readFasta(fastaName[i], validaaSeq, blastpath, evalue, ivalue, cddflag, fakedbpath, dbpath, blastpgppath, hhsearchpath, hhblitspath)

        print outHome, ws_cytosolic, ws_extracellular, ws_lipidfacing, ws_porefacing, pssm21_path, proteinName, rpath, validaaSeq

        print "Predicting using SVMs"
        boctopus_predictsvm.startPredicting(outHome, ws_cytosolic, ws_extracellular, ws_lipidfacing, ws_porefacing, pssm21_path, proteinName, rpath, validaaSeq, path_seqfile)

        print "Predicting using HMM"
        print hmmfilename, pName, modHome, prfSuffix
# from path_seqfile, one can locate svmoutHome, outHome, tmpHome
        boctopus_callHMM.start(hmmfilename, pName, modHome, prfSuffix, path_seqfile)

    return 1

#}}}



readList(fastainputfile, blastpath, modHome, hmmfilename, ws_cytosolic, ws_extracellular, ws_lipidfacing, ws_porefacing, rpath, fakedbpath, dbpath, blastpgppath, hhsearchpath, hhblitspath)
