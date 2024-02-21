#!/usr/bin/env python3
################################################################################
## pylabplot_pmpar.py: Make pmpar plots, pylab-style
## Adam Deller, 26 Aug 2011
################################################################################
from matplotlib.ticker import MaxNLocator
import matplotlib
matplotlib.use('Agg')
import os, sys, math, pylab
from optparse import OptionParser

################################################################################
## Functions
################################################################################

# Fake an ellipse using an N-sided polygon
def Ellipse(xxx_todo_changeme, xxx_todo_changeme1, resolution=40, orientation=0, **kwargs):
    (x,y) = xxx_todo_changeme
    (rx, ry) = xxx_todo_changeme1
    xs = []
    ys = []
    for i in range(resolution):
        theta = 2*math.pi/resolution*i + orientation
        xs.append(x + rx * math.cos(theta))
        ys.append(y + ry * math.sin(theta))
    return pylab.Polygon(list(zip(xs, ys)), **kwargs)

# Run pmpar, optionally with proper motion subtracted, get results
def pmpar(pmparexec, filename, subtractpm=False, fixlines=[]):
    os.system("rm -f junkpmparout.jjj")
    os.system("%s %s > junkpmparout.jjj" % (pmparexec, filename))
    junklines = open("junkpmparout.jjj").readlines()
    fitlines = []
    for line in junklines:
        if "pi" in line or "RA" in line or "Dec" in line or "mu_a" in line or "mu_d" in line or "poch" in line:
            fitlines.append(line.split('+')[0] + '\n')
    os.system("rm -f junkpmparout.jjj") 
    runfilename = filename
    if len(fixlines) > 0:
        runfilename = "junkpmparin.jjj"
        os.system("rm -f " + runfilename)
        output = open(runfilename, "w")
        for line in fixlines:
            output.write(line)
        for line in open(filename).readlines():
            if not "pi" in line and not "RA" in line and not "Dec" in line and not "mu_a" in line and not "mu_d" in line and not "poch" in line:
                output.write(line)
        output.close()
    runline = "%s %s > /dev/null" % (pmparexec, runfilename)
    if subtractpm:
        runline += " -om"
    os.system(runline)
    os.system("rm -f junkpmparin.jjj")
    tin = open("pmpar_t")
    tlines = tin.readlines()
    tin.close()
    ein = open("pmpar_e")
    elines = ein.readlines()
    ein.close()
    ttimes = []
    tras =   []
    tdecs =  []
    etimes = []
    eras =   []
    edecs =  []
    eraerr = []
    edecerr= []
    pras =   []
    pdecs =  []
    for line in tlines:
        splitline = line.split()
        ttimes.append(float(splitline[0]))
        tras.append(float(splitline[1]))
        tdecs.append(float(splitline[2]))
    for line in elines:
        splitline = line.split()
        etimes.append(float(splitline[0]))
        eras.append(float(splitline[1]))
        edecs.append(float(splitline[3]))
        eraerr.append(float(splitline[2]))
        edecerr.append(float(splitline[4]))
        pras.append(float(splitline[5]))
        pdecs.append(float(splitline[6]))
    return ttimes, tras, tdecs, etimes, eras, edecs, eraerr, edecerr, pras, pdecs, fitlines

################################################################################
## Main code
################################################################################
if __name__ == "__main__":
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option("-f", "--filename", dest="filename", default="",
                      help="The pmpar input file")
    parser.add_option("--file2", dest="file2", default="",
                      help="A second pmpar input file to overplot")
    parser.add_option("--plottype", dest="plottype", default="png",
                      help="The file type for the plot (def. png)")
    parser.add_option("--pmparexec", dest="pmparexec", default="pmpar",
                      help="Path to the pmpar executable")
    parser.add_option("--docolour", dest="docolour", default=False,
                      action="store_true", help="Make colour, not greyscale")
    parser.add_option("--target", dest="target", default="", help="Target name")
    parser.add_option("--nopmsubtract", dest="nopmsubtract", default=False,
                      action="store_true", help="Don't subtract PM for " + \
                      "vs. time plots")
    parser.add_option("--fontsize", dest="fontsize", default="16",
                      help="Font size for the plot")
    parser.add_option("--dotitle", dest="dotitle", default=False,
                      action="store_true", help="Print a title")
    parser.add_option("--ralim1", dest="ralim1", default="",
                     help="Limits for the ra offset in RA/dec plot (mas)")
    parser.add_option("--ralim2", dest="ralim2", default="",
                     help="Limits for the ra offset in RA vs time plot (mas)")
    parser.add_option("--declim1", dest="declim1", default="",
                     help="Limits for the dec offset in RA/dec plot (mas)")
    parser.add_option("--legendlocation", dest="legendlocation", default="",
                      help="Where to put the legend in each plot")
    parser.add_option("-b", "--bootstrapplot", dest="bootstrapplot",default="",
                      help="Plot the contents of this bootstrap file")
    parser.add_option("--fitlabel", dest="fitlabel", default="Best astrometric fit",
                      help="Label for the fit line on the plot legend")
    parser.add_option("--measlabel", dest="measlabel", default="Measured positions",
                      help="Label for the measured points line on the plot legend")
    (options, junk) = parser.parse_args()
    pmparexec       = options.pmparexec
    filename        = options.filename
    file2           = options.file2
    plottype        = options.plottype
    ralim1          = options.ralim1
    ralim2          = options.ralim2
    declim1         = options.declim1
    docolour        = options.docolour
    dotitle         = options.dotitle
    nopmsubtract    = options.nopmsubtract
    font_size       = float(options.fontsize)
    target          = options.target
    legendlocation  = options.legendlocation
    fitlabel        = options.fitlabel
    measlabel       = options.measlabel
    bootstrapplot   = options.bootstrapplot
    if len(junk) > 0:
        parser.error("You can only supply the allowed options")
    if filename=="":
        parser.error("You must supply a filename to run on with -f or --filename")
    if legendlocation == "":
        legendlocation="best"
    fitlinestyle = 'k--'
    fitpointsstyle = 'kx'
    bothfitstyle = 'kx--'
    fitpointmarkersize = 5
    if docolour:
        fitlinestyle = 'b--'
        fitpointsstyle = 'bx'
        bothfitstyle = 'bx--'
    errorlinecolour ='k'
    error2linecolour='g'
    if docolour:
        errorlinecolour = 'r'
    measurementlinestyle = 'ko'
    measurement2linestyle = 'go'
    if docolour:
        measurementlinestyle = 'ko'
    
    ### FIRST DO RA VS DEC
    ttimes, tras, tdecs, etimes, eras, edecs, eraerr, edecerr, pras, pdecs, fitlines = pmpar(pmparexec, filename)
    if not file2 == "":
        ttimes2, tras2, tdecs2, etimes2, eras2, edecs2, eraerr2, edecerr2, pras2, pdecs2, fitlines2 = pmpar(pmparexec, file2, False, fitlines)
    plotfile = "radec." + plottype
    if not target == "":
        plotfile = "radec." + target + "." + plottype
    junkforlegend = pylab.plot([],[],linestyle='None', marker='o',markeredgecolor=errorlinecolour, mfc="None", markersize=10, markeredgewidth=1.5)
    if not file2 == "":
        junkforlegend2 = pylab.plot([],[],linestyle='None', marker='o',markeredgecolor=error2linecolour, mfc="None", markersize=10, markeredgewidth=1.5)
    pylab.clf()
    pylab.rc("axes", linewidth=2.0)
    pylab.rc("lines", markeredgewidth=2.0)  # this also grows the points/lines in the plot..
    ax  = pylab.subplot(111)
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(font_size*0.9)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(font_size*0.9)
    if dotitle:
        pylab.title("Motion of " + target)
    #pylab.ylim(-beamhwhm*2,beamhwhm*2)
    #pylab.xlim(-beamhwhm*2,beamhwhm*2)
    if ralim1 != "":
        print("Setting ralims")
        splitralim = ralim1.split(',')
        pylab.xlim(float(splitralim[0]), float(splitralim[1]))
    if declim1 != "":
        print("Setting declims")
        splitdeclim = declim1.split(',')
        pylab.ylim(float(splitdeclim[0]), float(splitdeclim[1]))
    pylab.xlabel("Right Ascension offset (mas)", fontsize=font_size)
    pylab.ylabel("Declination offset (mas)", fontsize=font_size)
    fitline = pylab.plot(tras,tdecs,fitlinestyle)
    if ralim1 != "":
        print("Setting ralims")
        splitralim = ralim1.split(',')
        pylab.xlim(float(splitralim[0]), float(splitralim[1]))
    #ax.xaxis.set_major_locator(MaxNLocator(8, prune='lower'))
    #ax.yaxis.set_major_locator(MaxNLocator(8))
    fitpoints = pylab.plot(pras,pdecs,fitpointsstyle,markersize=fitpointmarkersize)
    for ra, dec, raerr, decerr in zip(eras, edecs, eraerr, edecerr):
        errorellipse = Ellipse((ra,dec), (raerr, decerr), facecolor='none', edgecolor=errorlinecolour, linewidth=1.5)
        ax.add_patch(errorellipse)
    if not file2 == "":
        for ra, dec, raerr, decerr in zip(eras2, edecs2, eraerr2, edecerr2):
            errorellipse = Ellipse((ra,dec), (raerr, decerr), facecolor='none', edgecolor=error2linecolour, linewidth=1.5)
            ax.add_patch(errorellipse)
    #ellipseforlegend = pylab.Circle((0,0), 1, facecolor='r', edgecolor=errorlinecolour, linewidth=2)
    #if file2 == "":
    #    pylab.legend((fitline,junkforlegend),("Best astrometric fit", "Error ellipses for observations"), numpoints=1,loc=legendlocation)
    #else:
    #    pylab.legend((fitline,junkforlegend,junkforlegend2),("Best astrometric fit", "Error ellipses for observations","Secondary error ellipses"), numpoints=1,loc=legendlocation)
    # Save the file
    pylab.savefig(plotfile, dpi=300)
    pylab.clf()
    
    ### NOW DO THE RA VS TIME WITH PM SUBTRACTED
    if not nopmsubtract:
        ttimes, tras, tdecs, etimes, eras, edecs, eraerr, edecerr, pras, pdecs, fitlines = pmpar(pmparexec, filename, True)
        if not file2 == "":
            ttimes2, tras2, tdecs2, etimes2, eras2, edecs2, eraerr2, edecerr2, pras2, pdecs2, fitlines2 = pmpar(pmparexec, file2, True, fitlines)
    if nopmsubtract:
        plotfile = "ratime." + plottype
        if not target == "":
            plotfile = "ratime." + target + "." + plottype
    else:
        plotfile = "ratime_nopm." + plottype
        if not target == "":
            plotfile = "ratime_nopm." + target + "." + plottype
    pylab.rc("axes", linewidth=2.0)
    pylab.rc("lines", markeredgewidth=2.0)  # this also grows the points/lines in the plot..
    ax  = pylab.subplot(111)
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(font_size*0.9)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(font_size*0.9)
    if dotitle:
        pylab.title("Parallax signature of %s in right ascension" % target)
    pylab.xlabel("Time (MJD)", fontsize=font_size)
    pylab.ylabel("Right ascension offset (mas)", fontsize=font_size)
    fitline = pylab.plot(ttimes,tras,fitlinestyle)
    fitpoints = pylab.plot(etimes,pras,fitpointsstyle, markersize=fitpointmarkersize)
    junk = pylab.plot([],[],bothfitstyle,markersize=fitpointmarkersize,label=fitlabel)
    if len(ralim2) > 1:
        splitralim = ralim2.split(',')
        pylab.ylim(float(splitralim[0]), float(splitralim[1]))
    #ax.xaxis.set_major_locator(MaxNLocator(6, prune='lower'))
    #ax.yaxis.set_major_locator(MaxNLocator(6))
    measurementline = pylab.errorbar(etimes, eras, yerr=eraerr, fmt=measurementlinestyle, markersize=3, label=measlabel)
    if not file2 == "":
        measurementline2 = pylab.errorbar(etimes2, eras2, yerr=eraerr2, fmt=measurement2linestyle, markersize=3, mfc=error2linecolour, mec=error2linecolour, label="Measured positions (alternate)")
    pylab.legend(numpoints=1, loc=legendlocation)
    #pylab.legend(numpoints=1, loc=legendlocation, prop={'size':(font_size*4)/5})
    pylab.savefig(plotfile)
    pylab.clf()
    
    ### NOW DO THE DEC VS TIME WITH PM SUBTRACTED
    if nopmsubtract:
        plotfile = "dectime." + plottype
        if not target == "":
            plotfile = "dectime." + target + "." + plottype
    else:
        plotfile = "dectime_nopm." + plottype
        if not target == "":
            plotfile = "dectime_nopm." + target + "." + plottype
    pylab.rc("axes", linewidth=2.0)
    pylab.rc("lines", markeredgewidth=2.0)  # this also grows the points/lines in the plot..
    ax  = pylab.subplot(111)
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(font_size*0.9)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(font_size*0.9)
    if dotitle:
        pylab.title("Parallax signature of %s in declination" % target)
    pylab.xlabel("Time (MJD)", fontsize=font_size)
    pylab.ylabel("Declination offset (mas)", fontsize=font_size)
    fitline = pylab.plot(ttimes,tdecs,fitlinestyle,label=fitlabel)
    fitpoints = pylab.plot(etimes,pdecs,fitpointsstyle,markersize=fitpointmarkersize)
    measurementline = pylab.errorbar(etimes, edecs, fmt=measurementlinestyle, yerr=edecerr, markersize=3, label=measlabel)
    if not file2 == "":
        measurementline = pylab.errorbar(etimes2, edecs2, fmt=measurement2linestyle, yerr=edecerr2, markersize=3, label="Measured positions (alternate)")
    #ax.xaxis.set_major_locator(MaxNLocator(6, prune='lower'))
    #ax.yaxis.set_major_locator(MaxNLocator(6))
    #pylab.legend(numpoints=1, loc=legendlocation, prop={'size':(font_size*4)/5})
    pylab.legend(numpoints=1, loc=legendlocation)
    pylab.savefig(plotfile)
    pylab.clf()
    
