import os, sys, string

gap_threshold = 60.0

def start(filename):
	f     = open(filename, "r")
	lines = f.readlines()
	f.close()
	
	count = 0
	f = open(filename + ".filtered", "w")
	for line in lines:
		line  = line.strip()
		aline = string.split(line)
		frac_gap = (aline[1].count("-")*100.0)/len(aline[1])

		if frac_gap < gap_threshold:
			#print line
			f.write(line + "\n")
			count += 1

	print filename, "discared", len(lines)-count, len(lines)
	f.close()

	return


#filename = sys.argv[1]
#start(filename)
