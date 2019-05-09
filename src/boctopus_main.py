#!/usr/bin/env python
# ChangeLog 2015-05-21
#   Nanjiang added simple argument parsing and error handelling
#   too many hardcoded places need to be changed to make it really standalone
#   and do not write anything within the folder where the scripts are located.
#   I will just make a quick and dirty solution.

import os, sys, string
import yaml

sys.path.append("/usr/local/lib/python2.7/dist-packages")

import tempfile
import shutil
import myfunc

## PATHs to change begin


TMPPATH = "/tmp"
if os.path.exists("/scratch"):
    TMPPATH = "/scratch"

rundir = os.path.realpath(os.path.dirname(__file__))

usage = """
USAGE: %s SEQFILE OUTPATH
"""%(sys.argv[0])

## List of fasta files to run
# argument parsing
numArgv = len(sys.argv)
if numArgv != 3:
    print >> sys.stderr, "Syntax error for %s!"%(sys.argv[0])
    print >> sys.stderr, usage
    sys.exit(1)


infile    = sys.argv[1]
outpath   = sys.argv[2]

if not os.path.exists(infile) or os.path.getsize(infile)<1:
    print >> sys.stderr, "%s: infile '%s' does not exist or empty. Exit."%(sys.argv[0], infile)
    sys.exit(1)
if not os.path.exists(outpath):
    try:
        os.makedirs(outpath)
    except OSError:
        print >> sys.stderr, "%s: failed to create outpath '%s'. Exit."%(sys.argv[0], outpath)

infile = os.path.realpath(infile)
outpath = os.path.realpath(outpath)

##################### Load configuration ######################
with open("config.yml", 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

HHBLITS_PATH = cfg['HHBLITS_PATH']
HHLIB = HHBLITS_PATH + "/lib/hh"
reformatpath = "%s/scripts/reformat.pl"%(HHLIB)
os.environ["HHLIB"] = HHLIB

HHBLITS_DB_PATH = cfg['HHBLITS_DB_PATH']

blastpath = cfg['blastpath']
blastpgppath = "%s/blastpgp"%(blastpath)

rpath = cfg['rpath']

modHome   = "%s/modhmm/"%(rundir)
##################  End of configuration ####################


try:
    tmpdir = tempfile.mkdtemp(prefix="%s/boctopus2_"%(TMPPATH))
except OSError:
    print >> sys.stderr, "%s: Failed to create tmpdir"%(sys.argv[0])
    sys.exit(1)

## required by blastpgp to convert alignment formats
fakedbpath = "%s/dummy_db/1a0s.fas"%(rundir)

## PATHs to change end

## ------------------------------------------------------
## Do not change code below this line.
## ------------------------------------------------------


def start_boctopus(infile, blastpath, modHome, hmmfilename, ws_cytosolic, ws_extracellular, ws_lipidfacing, ws_porefacing, \
    fakedbpath, dbpath, blastpgppath, hhsearchpath, hhblitspath, rpath):
    print "boctopus2 will start with ", infile

#     f = open(infile, "r")#{{{ DELETED
#     lines = f.readlines()
#     f.close()
# 
#     pname   = []
#     seqname = []
#     tempseq = ""
#     for line in lines:
#         line = line.strip()
# 
#         if line.startswith(">"):
#             pname.append(line[1:])
#             if len(tempseq) > 0:
#                 seqname.append(tempseq)
#             tempseq = ""
#         else:
#             tempseq += line
# 
#     if len(tempseq) > 0:
#         seqname.append(tempseq)
# 
#     print pname
#     print seqname
# 
#     if len(pname) != len(seqname):
#         print "number of pnames and seqs not the same."
#     else:#}}}

    # rewrite sequence reading part
    (seqidlist, seqannolist, seqlist) = myfunc.ReadFasta(infile)
    if len(seqidlist) <= 0:
        print >> sys.stderr, "No valid sequences read from file '%s'"%(infile)
        return 1

    #for i in range(0, len(pname)):
    for i in xrange(len(seqidlist)):
        seqid = seqidlist[i]
        seq = seqlist[i]
        seqanno = seqannolist[i]
        print "processing ", i , seqanno

        subtmpdir = "%s/seq_%d"%(tmpdir, i)
        if os.path.exists(subtmpdir):
            shutil.rmtree(subtmpdir)
        os.makedirs(subtmpdir)

        singleseqfile = "%s/query.fa"%(subtmpdir)
        myfunc.WriteFile(">%s\n%s\n"%(seqanno, seq), singleseqfile, mode="w", isFlush=True)

        if not os.path.exists(singleseqfile):
            print >> sys.stderr, "Failed to write to singleseqfile %s"%(singleseqfile)
            continue

        command = "python "+ "%s/boctopus_startHMM.py "%(rundir) + singleseqfile + " " + blastpath + " " + modHome + " " + hmmfilename + " " + ws_cytosolic + " " + ws_extracellular + " " + ws_lipidfacing + " " + ws_porefacing + " " + rpath+ " " +fakedbpath+\
" " + dbpath+ " " + blastpgppath+ " " + hhsearchpath + " " + hhblitspath
        print command
        os.system(command)
        outpath_this_seq = "%s/seq_%d"%(outpath, i)
        if not os.path.exists(outpath_this_seq):
            os.makedirs(outpath_this_seq)
        filepair_to_copy = [
                ("%s/query.fa"%subtmpdir, "%s/query.fa"%outpath_this_seq),
                ("%s/output/query_ioIOS.prf.txt_svm_topo.png"%subtmpdir, "%s/query.predict.png"%(outpath_this_seq)),
                ("%s/output/query_topologies.txt"%(subtmpdir), "%s/query_topologies.txt"%outpath_this_seq),
                ("%s/svmoutput/query_ioIOS.prf.txt"%subtmpdir, "%s/profile.txt"%outpath_this_seq),
                ("%s/pssm/query.filtered.pssmvals"%subtmpdir, "%s/pssm.txt"%(outpath_this_seq))

        ]
        for tup in filepair_to_copy:
            shutil.move(tup[0], tup[1])


    return


sys.stderr.write("""
****************************************************************************
BOCTOPUS_v2: improved topology prediction of transmembrane beta-barrel proteins
****************************************************************************

If you use BOCTOPUS_v2, please cite the following paper:
Hayat, Sikander and Elofsson, Arne, (to submit)

Related -
Hayat, Sikander and Elofsson, Arne, Bioinformatics, 28(4), 516--522, 2012
"BOCTOPUS: improved topology prediction of transmembrane beta-barrel proteins"

-----------------------------------------------------------------------------


""")



## HMM file used
hmmfilename = "%s/st911_Noss.hmg"%(rundir)

## Do not change window size. Else you have to re-train SVMs.
#ws_cytosolic     = "31"
#ws_extracellular = "31"
#ws_porefacing    = "13"
#ws_lipidfacing   = "11"
ws_cytosolic     = "23"
ws_extracellular = "31"
ws_porefacing    = "13"
ws_lipidfacing   = "9"


currdir = os.getcwd()
os.chdir(rundir)

start_boctopus(infile, blastpath, modHome, hmmfilename, ws_cytosolic, ws_extracellular, ws_lipidfacing, ws_porefacing, \
        fakedbpath, HHBLITS_DB_PATH, blastpgppath, reformatpath, HHBLITS_PATH, rpath)

os.chdir(currdir)

if os.path.exists(tmpdir):
    shutil.rmtree(tmpdir)
