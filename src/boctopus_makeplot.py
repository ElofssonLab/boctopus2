#!/usr/bin/env python
import os, sys, string, math

sys.path.append("/usr/local/lib/python2.7/dist-packages")


import matplotlib
from matplotlib import *
matplotlib.use('Agg')
from pylab import *

#print "directory for matplotlib", os.environ.get('MPLCONFIGDIR')

color  = ['red', 'blue', 'orange', 'green', 'cyan']
marker = ['o', 'v', '.', 's', 'x', 's', 'd', '<', '1']

majorLocator   = MultipleLocator(25)
majorFormatter = FormatStrFormatter('%d')
minorLocator   = MultipleLocator(5)
minorLocator_y = MultipleLocator(0.1)

dataEncoding    = {'I':1.38, 'O':1.42, 'D': 1.4, 'U': 1.4, 'M':1.4, 'X':1.4}
dataEncodingObs = {'i':1.28, 'o':1.32, 'D': 1.3, 'U': 1.3, 'M':1.3, 'X':1.3}
colorEncoding   = {'i':'green', 'o':'blue', 'X':'gray', 'D': 'red', 'U': 'orange', 'M': 'magenta'}
markerEncoding  = {'i':'o', 'o':'v', 'X':'.' , 'D': 's' , 'U': 'x', 'M': 'x'}

label_disc_dict = {'i':"pore-facing", 'o':"lipid-facing", 'I':"inner-loop", 'O':"outer-loop", 'M':"Membrane", "S":"Psipred(E)"}


def appendObsTopo(proDict_dataset):
	x     = []
	count = 0
	y     = []

	for data in proDict_dataset:
		status = data[4]
		#print status, dataEncoding[status]		
		x.append(count)
		count += 1
		y.append(dataEncoding[status])

	return y, x


def appendPredTopo(proDict_predTopo):
        x     = []
        count = 0
        y     = []

        for data in proDict_predTopo:
		for i in range(0, len(data)):
	                status = data[i]
			#print status, dataEncoding[status]
                	x.append(count)
	                count += 1
        	        y.append(dataEncoding[status])

        return y, x


def appendPredTopo_top(proDict_predTopo):	## topcons style output figure
        x     = []
        count = 0
        y     = []

        colorcode  = []
        linewidths = []

        lines  = []
        linesx = []

        templinem = []
        templinei = []
        templineo = []

        templinemx = []
        templineix = []
        templineox = []

        for data in proDict_predTopo:
                for i in range(0, len(data)):
                        status = data[i]

                        if status not in ["I", "O"]:#print status, dataEncoding[status]         
                                templinem.append(dataEncoding["I"])
                                templinemx.append(count)
                        else:
                                if len(templinem) > 0:
                                        lines.append(templinem)
                                        linesx.append(templinemx)
                                        colorcode.append("gray")
                                        linewidths.append(0.06)
                                templinem  = []
                                templinemx = []

                        if status in ["I"]:#print status, dataEncoding[status]         
                                templinei.append(dataEncoding[status])
                                templineix.append(count)
                        else:
                                if len(templinei) > 0:
                                        lines.append(templinei)
                                        linesx.append(templineix)
                                        colorcode.append("orange")
                                        linewidths.append(0.02)
                                templinei  = []
                                templineix = []

                        if status in ["O"]:#print status, dataEncoding[status]         
                                templineo.append(dataEncoding[status])
                                templineox.append(count)
                        else:
                                if len(templineo) > 0:
                                        lines.append(templineo)
                                        linesx.append(templineox)
                                        colorcode.append("green")
                                        linewidths.append(0.02)
                                templineo  = []
                                templineox = []

                        #x.append(count)
                        count += 1
                        #y.append(dataEncoding[status])

                if len(templinei) > 0:
                        lines.append(templinei)
                        linesx.append(templineix)
                        colorcode.append("orange")
                        linewidths.append(0.02)

                if len(templineo) > 0:
                        lines.append(templineo)
                        linesx.append(templineox)
                        colorcode.append("green")
                        linewidths.append(0.02)

                if len(templinem) > 0:
                        #print templinem
                        lines.append(templinem)
                        linesx.append(templinemx)
                        colorcode.append("black")
                        linewidths.append(0.06)

        return lines, linesx, colorcode, linewidths


def appendProfile(labels, proDict_profiles):
	data = zeros([len(labels), len(proDict_profiles)], float)
        y    = zeros(len(proDict_profiles), int)

        for i in range(0, len(data)):
        	for j in range(0, len(data[0])):
                        y[j] = j

       	for i in range(0, len(proDict_profiles[0])-1):
       		for j in range(0, len(proDict_profiles)):
                	data[i][j] = float(proDict_profiles[j][i])

	return data, y


def start(proDict_dataset, proDict_profiles, proDict_predTopo, labels, prName, path_seqfile): #{{{
        Home = path_seqfile + "/"
        os.environ['MPLCONFIGDIR'] = path_seqfile
        svmoutHome = path_seqfile + "/svmoutput/"
        outHome    = path_seqfile + "/output/"
        tmpHome    = path_seqfile + "/tmp/"

	for protein in proDict_profiles:
                if len(proDict_predTopo[protein][0]) != len(proDict_profiles[protein]):
                       print "makeplot: Length error detected. Possible cause in getData", len(proDict_profiles[protein]), len(proDict_predTopo[protein])
                       sys.exit()
		else:
			#print protein, labels
			data, y = appendProfile(labels, proDict_profiles[protein])

			#fig = figure()

			figxlength = len(proDict_profiles[protein])/20
			fig = figure(figsize=(figxlength, 5))
			ax  = fig.add_subplot(111)

			plt.tight_layout()

			for i in range(0, len(data)-1):
				plotValues(ax, y, data[i], color[i], 1.0)#, marker[i])

                        if len(labels) == 5:
                                legend((label_disc_dict[labels[0]], label_disc_dict[labels[1]], label_disc_dict[labels[2]], label_disc_dict[labels[3]]), prop={'size':12}, bbox_to_anchor=(0., 0.95, 1., .102), loc=3, ncol=4, mode="expand", borderaxespad=0.)
                        if len(labels) == 6:
                                legend((label_disc_dict[labels[0]], label_disc_dict[labels[1]], label_disc_dict[labels[2]], label_disc_dict[labels[3]], label_disc_dict[labels[4]]), bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, mode="expand", borderaxespad=0.)



			y_all, x_all, colorcode, lwds = appendPredTopo_top(proDict_predTopo[protein])

                        for i in range(0, len(y_all)):
                                p = Rectangle([x_all[i][0], y_all[i][0]], x_all[i][-1]-x_all[i][0]+1, lwds[i], fill=True, linewidth=0, edgecolor = None, facecolor=colorcode[i])
                                ax.add_patch(p)

                        ax.text(x_all[i][-1]/3.0, 1.55, '%s' % ("Predicted topology (HMM-output)"))
                        ax.text(x_all[i][-1]/3.0, 1.06, '%s' % ("per-residue prediction (SVM-output)"))

                        x1,x2,y1,y2 = ax.axis()
                        ax.axis((x1, x2, 0, 1.8))

			ax.xaxis.set_major_locator(majorLocator)

			locs, xlabels = plt.xticks()
			plt.xticks(rotation=45)

			#for the minor ticks, use no labels; default NullFormatter
			ax.xaxis.set_minor_locator(minorLocator)
                        ax.yaxis.set_minor_locator(minorLocator_y)

			ax.grid(True)
                        ax.set_ylim(-0.1)

		        ylabel('Value')
		        xlabel('Residue Number (' + 'PROTEIN NAME:' + protein + ")")

			hmmimage = outHome + protein + '_svm_topo.png'
			fig.savefig(hmmimage, dpi=300)
		
			print outHome + protein + '_svm_topo.png Finished' 

	return
#}}}

def plotValues(ax, x1, y1, c, mylw):#, m):
	#ax.plot(x1, y1, color=c, linewidth=mylw)#, marker=m)#, x1, y1, 'green', "k")
	ax.bar(x1, y1, color=c, alpha=0.5, edgecolor = "none", linewidth=0)

        return


def showPlot():
	ylabel('Value')
        xlabel('Residue Number')
        #title('boctopus predictions')
        show()

	return

	
def test(x,y):
	x1=arange(0,20,0.01)
	y1=2*sin(2*pi*(x1-1/4))

	x2= x1-3
	y2= y1-3
	aa = figure(1)

	plot(x1, y1, '.', x2, y2, '-')

	plot(x1, y1, '.')
	xlabel('x-axis')
	ylabel('y-axis')
	title(r'$boctopus predictions$')

	show()

	return 

#test()
