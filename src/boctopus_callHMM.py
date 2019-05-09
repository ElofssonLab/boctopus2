import os, sys, string


import boctopus_parseHMMResult
import boctopus_getplotdata



def readSinglehmg(hmmList, fileList, modHome, pName, Home, svmoutHome, outHome, tmpHome):
        #os.system("ls " + Home + hmmList + " > " + tmpHome + "minimal.txt")
        os.system("ls " +  hmmList + " > " + tmpHome + "minimal.txt")
	os.system("ls " + svmoutHome + fileList + " > " + tmpHome + "prFile.txt")

        command = modHome + "./modhmms_octopus -m " + tmpHome + "minimal.txt -s " + tmpHome + "prFile.txt -f prf -o " + "resultDir -r " +  modHome + "replacement_letter_multi.rpl -L -M DP -v --nopostout -u --nolabels --viterbi --path --labeloddsout --labelllout > " + outHome + pName + "_outfile.xml"
        os.system(command)
	
        boctopus_parseHMMResult.start(outHome + pName + "_outfile.xml", tmpHome+"prFile.txt", Home)
	boctopus_getplotdata.start(tmpHome+"prFile.txt", Home)

	return


def start(hmmList, pName, modHome, prfSuffix, path_seqfile):
        Home = path_seqfile + "/"
        svmoutHome = path_seqfile + "/svmoutput/"
        outHome    = path_seqfile + "/output/"
        tmpHome    = path_seqfile + "/tmp/"
        print hmmList, pName, modHome, prfSuffix

        fileList = pName + prfSuffix

        readSinglehmg(hmmList, fileList, modHome, pName, Home, svmoutHome, outHome, tmpHome )

        return


#start("1", "simple_preseq_iom_mo.hmg", "2pilA")
