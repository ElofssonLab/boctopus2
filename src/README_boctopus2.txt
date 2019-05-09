boctopus_V2 NOTES_v1


Dependencies:
	a. hhblits executibles (version hhsuite-2.0.16) and sequence database (version uniprot20_2013_03) ftp://toolkit.genzentrum.lmu.de/pub/HH-suite/
	b. blastpgp (version blast-2.2.26) ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ 
	c. modhmm (executible provided in this folder modhmm/)

	SET THESE PATHS IN boctopus_main.py


Usage:

    USAGE: ./boctopus_main.py SEQFILE OUTPATH

Examples:

$ unzip boctopus2_newset_hhblits.zip

$ cd boctopus2_newset_hhblits/test

$ ../boctopus_main.py LPTD.fa out1


Output description:

	a. Predicted topology "[pname]_ioIOS.prf.txt_stateName.txt"
	I=inner loop
	O=outer loop
	i=pore-facing residue in a membrane strand
	o=lipid-facing residue in a membrane strand

	b. Image of predicted topology "[pname]_ioIOS.prf.txt_svm_topo.png"

	c. HMM path travesered "[pname]_ioIOS.prf.txt_statepath.txt"

	d. HHblits output "[pname].hhr"



FILES:
SVM: objects trained on the dataset with window sizes estimated using a cross-validation test.
svmtrain42_innerloop.Rdata
svmtrain42_lipidfacing.Rdata
svmtrain42_outerloop.Rdata
svmtrain42_porefacing.Rdata

R scripts (for SVM-stage predictions):
traninedsvmRadial_cyto.R  traninedsvmRadial_extra.R  traninedsvmRadial_lipid.R  traninedsvmRadial_pore.R


HMM: architecture and parameters st911_Noss.hmg


Python scripts:
boctopus_main.py         		## set paths and variables, calls boctopus_startHMM.py
boctopus_startHMM.py			## calls programs to generate sequence alignment, filter out gappy sequences, svm prediction, hmm prediction
boctopus_generatePSSM_HHblits.py	## generate sequence alignment
filter_gappyseqs.py			## filter out gappy sequences
boctopus_neighborJoin_profile_pssm.py	## create input file with PSSM values
boctopus_predictsvm.py			## svm prediction
boctopus_writepred2prf.py		## convert svm output to a profile that is used as input to HMM
boctopus_callHMM.py               	## call modhmm with input profile
boctopus_parseHMMResult.py		## parse HMM output and generate text output
boctopus_getplotdata.py  		## gather output data for plotting results
boctopus_makeplot.py                   	## plot output


ISSUES:
1. The header of input fasta file can have spaces. handle it correctly.

