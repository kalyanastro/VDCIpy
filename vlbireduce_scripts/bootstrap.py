#!/usr/bin/python3
"""
Adam Deller
Usage:
> astrometric_bootstrap.py -p J0332+5434 -n 100000 -f J0332+5434.1.modified.pmpar.in 
                    -s /home/kalyan/Desktop/projects/BD174/BD174A/results/BD174_01
for help
> astrometric_bootstrap.py -h
"""


import os, sys, re, subprocess, random, numpy, math
import astro_utils
from optparse import OptionParser
import matplotlib.pyplot as pyplot
from pylab_pmparplot import pmpar


def smooth(x, window_len=11, window='hanning'):
    """smooth the data using a window with requested size.
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len<3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window must be one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s=numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='valid')
    return y[(int(window_len/2-1)):-int((window_len/2))]    # Ashish added int()


def get_1sigma_confidence(a, numbins, valname, binvals=[], docolour=True, fullfitvals=[], useefac=False):
    sorteda = sorted(a)
    print("Median for", valname, "is ", sorteda[int(len(a)/2)])   # Ashish added int()
    numvals = len(a)
    offset = int(0.68*numvals)
    print("Offset is " + str(offset))
    numtrials = len(a) - offset
    minwidth = 9e99
    for i in range(numtrials):
        width = math.fabs(sorteda[i+offset] - sorteda[i])
        if width < minwidth:
            minwidth = width
            bestindex = i
    meanval = (sorteda[bestindex] + sorteda[bestindex+offset])/2.0
    xlimlow = meanval - 2.5*minwidth
    xlimhigh = meanval + 2.5*minwidth
    if len(binvals) == 0:
        binvals = []
        binincrement = float(5*minwidth)/numbins
        for i in range(numbins+1):
            binvals.append(xlimlow + i*binincrement)
    # (pdf, binedges) = numpy.histogram(a, bins=binvals, normed=True)
    (pdf, binedges) = numpy.histogram(a, bins=binvals, density=True)   # Ashish added it
    smoothedpdf = smooth(pdf,5,'hanning')
    #if valname == "inclination":
    #    print binedges, binvals
    binwidth = binedges[1] - binedges[0]
    xvals = []
    yvals = []
    print(len(pdf), numbins)
    for i in range(numbins):
        xvals.append(binedges[i] + binwidth/2)
        yvals.append(pdf[i])
    numpy.savetxt("pdf_" + valname + ".txt", (xvals, yvals))
    mode = 0
    smoothedmode = 0
    pdfmax = 0
    smoothedpdfmax = 0
    pdfsum = 0
    count = 0
    for p, s in zip(pdf, smoothedpdf):
        pdfsum += p
        if p > pdfmax:
            pdfmax = p
            mode = binedges[count] + binwidth/2.0
        if s > smoothedpdfmax:
            smoothedpdfmax = s
            smoothedmode = binedges[count] + binwidth/2.0
        count += 1

    modelinestyle = 'dotted'
    modelinecolour = 'k'
    confidencelinestyle='dashdot'
    confidencelinecolour='k'
    bootstraplinestyle='solid'
    bootstraplinecolour = 'k'
    leastsquareslinestyle = 'solid'
    leastsquareslinecolour = 'silver'
    if docolour:
        modelinestyle = 'solid'
        modelinecolour = 'r'
        confidencelinestyle='solid'
        confidencelinecolour='g'
        bootstraplinestyle='solid'
        bootstraplinecolour = 'k'
        leastsquareslinestyle = 'solid'
        leastsquareslinecolour = 'b'

    pyplot.clf()
    pyplot.plot(xvals, yvals, color=bootstraplinecolour, linestyle=bootstraplinestyle, linewidth=1.75)
    pyplot.gca().yaxis.set_visible(False)
    #pyplot.ylabel("Probability density")
    if valname == "parallax":
        # pyplot.xlabel("Parallax (mas)")
        pyplot.xlabel("Parallax [mas]", fontsize=16) # Ashish added fontsize
    elif valname == "pm_ra":
        pyplot.xlabel("Proper motion (R.A., [mas/yr])", fontsize=16)
    elif valname == "pm_dec":
        pyplot.xlabel("Proper motion (Decl., [mas/yr])", fontsize=16)
    if len(fullfitvals) > 0:
        tempmu = float(fullfitvals[2])
        tempsigma = float(fullfitvals[4])
        yvals = []
        for x in xvals:
            try:
                tempmodifier = 1/(math.sqrt(2*math.pi)*tempsigma)
                tempexponent = ((x - tempmu)**2)/(2*tempsigma*tempsigma)
            except ZeroDivisionError:
                tempmodifier = 0
                tempexponent = 0
            yvals.append(tempmodifier*math.exp(-tempexponent))
        pyplot.plot(xvals, yvals, color=leastsquareslinecolour, linestyle=leastsquareslinestyle, linewidth=1.75)
    if len(binvals) == 0:
        pyplot.xlim(binvals[0], binvals[-1])
    pyplot.axvline(mode, color=modelinecolour, linestyle=modelinestyle, linewidth=1.75)
    pyplot.axvline(sorteda[bestindex], color=confidencelinecolour, linestyle=confidencelinestyle, linewidth=1.75)
    pyplot.axvline(sorteda[bestindex+offset], color=confidencelinecolour, linestyle=confidencelinestyle, linewidth=1.75)
    efacstring = ""
    if useefac:
        efacstring = ".efac"
    pyplot.ylabel("Unnormalised pdf")
    pyplot.tight_layout() # Ashish added it
    pyplot.savefig('bootstrap_%s%s.pdf' % (valname, efacstring))
    pyplot.clf()

    # Now the smoothed version
    yvals = []
    for i in range(numbins):
        yvals.append(smoothedpdf[i])
    pyplot.plot(xvals, yvals, color=bootstraplinecolour, linestyle=bootstraplinestyle)
    pyplot.gca().yaxis.set_visible(False)
    #pyplot.ylabel("Probability density")
    if valname == "parallax":
        # pyplot.xlabel("Parallax (mas)")
        pyplot.xlabel("Parallax [mas]", fontsize=16) # Ashish added fontsize
    elif valname == "pm_ra":
        pyplot.xlabel("Proper motion (R.A., [mas/yr])", fontsize=16)
    elif valname == "pm_dec":
        pyplot.xlabel("Proper motion (Decl., [mas/yr])", fontsize=16)
    if len(fullfitvals) > 0:
        tempmu = float(fullfitvals[2])
        tempsigma = float(fullfitvals[4])
        yvals = []
        for x in xvals:
            try:
                tempmodifier = 1/(math.sqrt(2*math.pi)*tempsigma)
                tempexponent = ((x - tempmu)**2)/(2*tempsigma*tempsigma)
            except ZeroDivisionError:
                tempmodifier = 0
                tempexponent = 0
            yvals.append(tempmodifier*math.exp(-tempexponent))
        pyplot.plot(xvals, yvals, color=leastsquareslinecolour, linestyle=leastsquareslinestyle)
    if not len(binvals) == 0:
        pyplot.xlim(binvals[0], binvals[-1])
    
    pyplot.axvline(smoothedmode, color=modelinecolour, linestyle=modelinestyle)
    pyplot.axvline(sorteda[bestindex], color=confidencelinecolour, linestyle=confidencelinestyle)
    pyplot.axvline(sorteda[bestindex+offset], color=confidencelinecolour, linestyle=confidencelinestyle)
    pyplot.tight_layout() # Ashish added it
    pyplot.savefig('smoothedbootstrap_%s.pdf' % (valname))
    bracket = (sorteda[bestindex], sorteda[bestindex+offset], smoothedmode)
    return bracket


usage = "usage: %prog [options]"
pmpar_exe = 'pmpar '
parser = OptionParser(usage)

parser.add_option("-p", "--pulsar", dest="pulsar", default="local",
                  help="pulsar name eg 1023+0038", metavar="PULSAR")
parser.add_option("-n", "--numtrials", dest="numtrials", default=4000,
                  help="Number of bootstrap trials")
parser.add_option("-f", "--filename", dest="filename", default="",
                  help ="Filename (default <pulsar>.jmfit.pmpar.in)")
parser.add_option("-s", "--soldir", dest="soldir", default="",
                  help="The directory where the pmpar input file can be found")
parser.add_option("--nopmsubtract", dest="nopmsubtract", default=False,
                  action="store_true", help="In plots, don't subtract proper motion")
parser.add_option("--useefac", dest="useefac", default=False,
                  action="store_true", help="Use scale factor rather than quadrature addition for lsq comparison in plot (True/False)")
parser.add_option("-o", "--orbitfitfile", dest="orbitfitfile", default="",
                  help="Fit for orbital motion as well using this pulsar config file")
parser.add_option("--orbitalbins", dest="orbitalbins", default="-1",
                  help="Override number of bins for Omega and inc.")
parser.add_option("-m", "--minepochs", dest="minepochs", default=6,
                  help="Minimum number of different dates for a valid trial")
parser.add_option("-r", "--readfromfile", dest="readfromfile", default="",
                  help="Read existing results from this file rather than rerunning")
parser.add_option("--parallaxbinvals", dest="parallaxbinvals", default="",
                  help="binning details for parallax in form minval,maxval,numbins")
(options, args) = parser.parse_args()

if len(args) > 0:
     parser.error("You can only supply the permitted options!")
pulsar = options.pulsar
soldir = options.soldir
if soldir == "": soldir = os.getcwd()

scriptdir = os.path.dirname(os.path.abspath(__file__)) # Ashish added it
os.chdir(soldir) # Ashish added it


numtrials = int(options.numtrials)
minepochs = int(options.minepochs)
# print(minepochs)
orbitalbins = int(options.orbitalbins)
subtractpm = not options.nopmsubtract
readfromfile = False
if not options.readfromfile == "":
    readfromfile = True
    bootstrapfile = options.readfromfile
    if not bootstrapfile[0] == '/':
        bootstrapfile = soldir + '/' + bootstrapfile
pmparfile = soldir + '/' + pulsar + ".jmfit.pmpar.in"
if options.filename != "":
    pmparfile = soldir + '/' + options.filename
doorbitalfit = False
if options.orbitfitfile != "":
    doorbitalfit = True
    
if not os.path.exists(pmparfile):
    print("Could not find file " + pmparfile + " - aborting!!!")
    sys.exit(1)

# Parse in pmpar file
pmparin = open(pmparfile, 'r')
pmparlines = pmparin.readlines()
pmparin.close()

# Get some info we need
headerendlinenum = 0
haveepoch = False
havedm    = False
numkeywords = 0
while re.search('=', pmparlines[headerendlinenum]) or \
      pmparlines[headerendlinenum].rstrip() == "" or \
      pmparlines[headerendlinenum].lstrip()[0] == '#':
    if re.search('DM', pmparlines[headerendlinenum]) or \
       re.search('dm', pmparlines[headerendlinenum]) or \
       re.search('Dm', pmparlines[headerendlinenum]):
        havedm = True
    if re.search('epoch', pmparlines[headerendlinenum]) or \
       re.search('Epoch', pmparlines[headerendlinenum]) or \
       re.search('EPOCH', pmparlines[headerendlinenum]):
        haveepoch = True
    if re.search('=', pmparlines[headerendlinenum]):
        if not "pi" in pmparlines[headerendlinenum]:
            numkeywords = numkeywords + 1
    headerendlinenum = headerendlinenum + 1
if not haveepoch:
    numkeywords = numkeywords+1
if not havedm:
    numkeywords = numkeywords+1

# Set up the arrays for positions, results etc
astrometric_points = []
dates = []
numpoints = 0
for line in pmparlines[headerendlinenum:]:
    if not line.lstrip() == "" and not line.lstrip()[0] == '#':
        astrometric_points.append(line)
        date = line.split()[0]
        if not date in dates:
            dates.append(date)
        numpoints = numpoints + 1
ref_ra = []
ref_dec = []
# print(astrometric_points)
if len(dates) <= minepochs: 
    print("Requested a minimum of %d epochs, but input file only contains %d" % (minepochs, len(dates)))
    sys.exit()

if not readfromfile:
    parallax = numpy.zeros((numtrials), numpy.float64)
    pm_dec = numpy.zeros((numtrials), numpy.float64)
    pm_ra = numpy.zeros((numtrials), numpy.float64)
    ref_ra_hour = numpy.zeros((numtrials), numpy.float64)
    ref_dec_deg = numpy.zeros((numtrials), numpy.float64)
    Omega = numpy.zeros((numtrials), numpy.float64)
    inclination = numpy.zeros((numtrials), numpy.float64)
    foldedinclination = numpy.zeros((numtrials), numpy.float64)
    xdot = numpy.zeros((numtrials), numpy.float64)
    #Run numtrials, saving the parallax, proper motion and position from each
    i = 0
    while i < numtrials:
        if numtrials > 10 and i%(numtrials/10) == 0:
            print(str(i/(numtrials/10)) + "0% done")
        tempout = open('bootstrapjunk.txt', 'w')
        tempout.writelines(pmparlines[0:headerendlinenum])
        if not havedm:
            tempout.write("#DM = ??\n\n")
        if not haveepoch:
            tempout.write("Epoch = 2007.0\n\n")
        uniquedates = []
        # print(numpoints)
        for j in range(numpoints):
            r = random.randint(0, numpoints-1)
            entrysplit = astrometric_points[r].split()
            found = False
            for d in uniquedates:
                if entrysplit[0] == d:
                    found = True
                    break
            if not found:
                uniquedates.append(entrysplit[0])
            tempout.write(astrometric_points[r] + '\n')
        tempout.close()
        # print(uniquedates)
        if len(uniquedates) >= minepochs:
            #Run pmpar, get the values
            if doorbitalfit: # Need to first search over orbital parameters
                os.system(scriptdir + '/' + "orbitalfit.py " + options.orbitfitfile + " bootstrapjunk.txt --noplot")
                os.system("mv best.txt bootstrapjunk.txt")
                orbitalvals = open("bestorbitalfit.txt").readlines()
                splitorbitalvals = orbitalvals[0].split()
                Omega[i] = float(splitorbitalvals[0])
                inclination[i] = float(splitorbitalvals[1])
                if inclination[i] > 90:
                    foldedinclination[i] = 180 - inclination[i]
                else:
                    foldedinclination[i] = inclination[i]
                if len(splitorbitalvals) > 2:
                    xdot[i] = float(splitorbitalvals[2])
            output = subprocess.getoutput(pmpar_exe + 'bootstrapjunk.txt').split('\n')
            startline = 1
            ref_ra.append(output[startline].split()[2])
            ref_dec.append(output[startline+1].split()[2])
            pm_ra[i] = float(output[startline+4].split()[2])
            pm_dec[i] = float(output[startline+5].split()[2])
            parallax[i] = float(output[startline+8].split()[2])
            i = i+1
        else:
            print("Skipping because only " + str(len(uniquedates)) + \
                  " dates were chosen...")
    
    #Make RA and dec numeric, and print to a file
    resultfile = open(soldir + '/' + pulsar + '.bootstrapresults', 'w')
    for i in range(numtrials):
        rasplit = ref_ra[i].split(':')
        decsplit = ref_dec[i].split(':')

        ref_ra_hour[i] = (float(rasplit[0]) + float(rasplit[1])/60.0 +
                          float(rasplit[2])/3600.0)
        if float(decsplit[0]) < 0:
            ref_dec_deg[i] = (float(decsplit[0]) - float(decsplit[1])/60.0 -
                              float(decsplit[2])/3600.0)
        else:
            ref_dec_deg[i] = (float(decsplit[0]) + float(decsplit[1])/60.0 +
                              float(decsplit[2])/3600.0)
        resultfile.write("%.15f " % ref_ra_hour[i])
        resultfile.write(str(ref_dec_deg[i]) + ' ')
        resultfile.write(str(pm_ra[i]) + ' ')
        resultfile.write(str(pm_dec[i]) + ' ')
        resultfile.write(str(parallax[i]) + ' ')
        if doorbitalfit:
            resultfile.write(str(Omega[i]) + ' ')
            resultfile.write(str(inclination[i]) + ' ')
            resultfile.write(str(xdot[i]) + ' ')
        resultfile.write('\n')
    resultfile.close()
else:
    resultfile = open(bootstrapfile)
    lines = resultfile.readlines()
    resultfile.close()
    numtrials = len(lines)
    if numtrials > 0 and len(lines[0].split()) > 5:
        doorbitalfit = True
    parallax = numpy.zeros((numtrials), numpy.float64)
    pm_dec = numpy.zeros((numtrials), numpy.float64)
    pm_ra = numpy.zeros((numtrials), numpy.float64)
    ref_ra_hour = numpy.zeros((numtrials), numpy.float64)
    ref_dec_deg = numpy.zeros((numtrials), numpy.float64)
    Omega = numpy.zeros((numtrials), numpy.float64)
    inclination = numpy.zeros((numtrials), numpy.float64)
    foldedinclination = numpy.zeros((numtrials), numpy.float64)
    xdot = numpy.zeros((numtrials), numpy.float64)
    values = []
    names = []
    values.append(parallax)
    names.append("Parallax")
    values.append(pm_ra)
    names.append("PM_RA")
    values.append(pm_dec)
    names.append("PM_Dec")
    values.append(ref_ra_hour)
    names.append("Pos_RA")
    values.append(ref_dec_deg)
    names.append("Pos_Dec")
    if doorbitalfit:
        values.append(Omega)
        names.append("Omega")
        values.append(inclination)
        names.append("inclination")
        values.append(foldedinclination)
        names.append("foldedinclination")
        values.append(xdot)
        names.append("xdot")
    count = 0
    for line in lines:
        splitline = line.split()
        ref_ra_hour[count] = splitline[0]
        ref_dec_deg[count] = splitline[1]
        pm_ra[count] = splitline[2]
        pm_dec[count] = splitline[3]
        parallax[count] = splitline[4]
        if doorbitalfit:
            Omega[count] = splitline[5]
            inclination[count] = splitline[6]
            if inclination[count] > 90:
                foldedinclination[count] = 180 - inclination[count]
            else:
                foldedinclination[count] = inclination[count]
            xdot[count] = splitline[7]
        count += 1
    for name, value in zip(names,values):
        # sort into order
        value = numpy.sort(value)
        median = value[int(numtrials/2)]   # Ashish added int()
        discardedbottom = 0
        bestbottom = -1
        minspan = 99999.9
        while discardedbottom < int(numtrials/3):    # Ashish added int()
            span = value[discardedbottom + int((2*numtrials)/3)] - value[discardedbottom]  # Ashish added int()
            if span < minspan:
                minspan = span
                bestbottom = discardedbottom
            discardedbottom += 1
        bestplus = value[bestbottom + int((2*numtrials)/3)] - median   # Ashish added int()
        bestminus = median - value[bestbottom]
        if not "Pos" in name:
            print("For %s, best span was %.3f, and the range is %.3f + %.3f - %.3f" % (name, minspan, median, bestplus, bestminus))
        else:
            hh = int(median)
            mm = int((median-hh)*60.0)
            ss = median*3600.0 - (hh*3600.0 + mm*60.0)
            bestplus *= 3600.0
            bestminus *= 3600.0
            print("For %s, best span was %.9f, and the range is %02d:%02d:%09.6f + %.6f - %.6f" % (name, minspan, hh, mm, ss, bestplus, bestminus))
    #sys.exit()

# Try running pmpar2pmparsystematics also
os.system("rm -f full.fit")
os.system(scriptdir + '/' + "pmpar2pmpar_systematics.py -f " + pmparfile + " > full.fit")
os.system("mv junk.txt " + pulsar + ".equad.pmparin")
fullfitlines = open("full.fit").readlines()
if options.useefac:
    rchisq = -1
    for line in fullfitlines:
        if "Reduced Chi^2 = " in line:
            rchisq = float(line.split()[-1])    
            break
    if rchisq < 0:
        print("Something went wrong when trying to use efac")
        sys.exit()
    os.system(scriptdir + '/' + "pmpar2pmpar_systematics.py -f " + pmparfile + " --scalefactor=" + str(math.sqrt(rchisq)) + "> full.fit")
    os.system("mv junk.txt " + pulsar + ".efac.pmparin")
fullfitlines = open("full.fit").readlines()
fullpxvals = []
fullpmravals = []
fullpmdecvals = []
for line in fullfitlines:
    if "pi" in line:
        fullpxvals = line.split()
    if "mu_a" in line:
        fullpmravals = line.split()
    if "mu_d" in line:
        fullpmdecvals = line.split()

#Do a histogram-based approach on the parallax
docolour = False
numbins = int(numtrials/200)
if numbins > 200: numbins = 200
if numbins < 20:
    numbins = 20
if options.parallaxbinvals != "":
    parallaxbinvals = numpy.linspace(float(options.parallaxbinvals.split(',')[0]), float(options.parallaxbinvals.split(',')[1]), int(options.parallaxbinvals.split(',')[2])+1)
    parallax_bracket = get_1sigma_confidence(parallax, len(parallaxbinvals)-1, "parallax", parallaxbinvals, docolour, fullpxvals, options.useefac)
else:
    parallax_bracket = get_1sigma_confidence(parallax, numbins, "parallax", [], docolour, fullpxvals, options.useefac)
pmra_bracket     = get_1sigma_confidence(pm_ra,    numbins, "pm_ra", [], docolour, fullpmravals, options.useefac)
pmdec_bracket    = get_1sigma_confidence(pm_dec,   numbins, "pm_dec", [], docolour, fullpmdecvals, options.useefac)
if doorbitalfit:
    usebins = numbins
    if orbitalbins > 0:
        usebins = orbitalbins
    omegabinvals = []
    incbinvals = []
    foldedincbinvals = []
    omegainc = 360.0/usebins
    incinc = 180.0/usebins
    foldedincinc = 90.0/usebins
    for i in range(usebins+1):
        omegabinvals.append(i*omegainc)
        incbinvals.append(i*incinc)
        foldedincbinvals.append(i*foldedincinc)
    Omega_bracket = get_1sigma_confidence(Omega, usebins, "Omega", omegabinvals, docolour)
    inclination_bracket = get_1sigma_confidence(inclination, usebins, "inclination", incbinvals, docolour)
    foldedinclination_bracket = get_1sigma_confidence(foldedinclination, usebins, "foldedinclination", foldedincbinvals, docolour, [])
    xdot_bracket = get_1sigma_confidence(xdot, usebins, "xdot", [], docolour)

#Work out the RA/Dec values
mean_ra_hour = numpy.mean(ref_ra_hour)
mean_dec_deg = numpy.mean(ref_dec_deg)
negdecstring = ""
if mean_dec_deg < 0:
    negdecstring = "-"
    mean_dec_deg = -mean_dec_deg
print("Mean ra (hours) is " + str(mean_ra_hour))
meanrahh = int(mean_ra_hour)
meanramm = int((mean_ra_hour - float(meanrahh))*60.0)
meanrass = (mean_ra_hour - float(meanrahh) - float(meanramm)/60.0)*3600.0
#meanrass = 3600*mean_ra_hour -(3600*meanrahh + 60*meanramm)
meandecdd = int(mean_dec_deg)
meandecmm = int((mean_dec_deg - float(meandecdd))*60)
#meandecmm = int(60*mean_dec_deg - 60*meandecdd)
meandecss = 3600*mean_dec_deg - (3600*meandecdd + 60*meandecmm)
if negdecstring == "-":
     meandecdd = -meandecdd

#Run once more with full file
output = subprocess.getoutput(pmpar_exe + pmparfile).split('\n')

#Print comparison to screen
#for i in range(startline):
#    print output[i] + "\n"
print("PMPAR (FULL DATASET) VALUES")
print("===========================")
os.system(pmpar_exe + pmparfile)
#print output[startline] + "\n" + output[startline+1] + "\n" + output[startline+4]
#print output[startline+5] + "\n" + output[startline+8] + "\n" +  output[startline+16] + "\n"
print("     BOOTSTRAP VALUES     ")
print("===========================")
print(("RA    = %02d:%02d:%09.7f +- %09.7f" % (meanrahh, meanramm, \
       meanrass, (numpy.std(ref_ra_hour))*3600.0)))
print(("Dec   = %02d:%02d:%09.7f +- %09.7f" % (meandecdd, meandecmm, \
       meandecss, (numpy.std(ref_dec_deg))*3600.0)))
#print "RA    = " + str(meanrahh) + ":" + str(meanramm) + ":" + \
#      str(meanrass) + " +/- " + ("%09.8f" % ((numpy.std(ref_ra_hour))*3600.0))
#print "Dec   = " + negdecstring + str(meandecdd) + ":" + str(meandecmm) + \
#      ":" + str(meandecss) + " +/- " + ("%09.8f" % ((numpy.std(ref_dec_deg))*3600.0))
print("mu_a  = " + str(numpy.mean(pm_ra)) + " +/- " + \
      ("%09.8f" % numpy.std(pm_ra)))
print("mu_d  = " + str(numpy.mean(pm_dec)) + " +/- " + \
      ("%09.8f" % numpy.std(pm_dec)))
print("pi    = " + str(numpy.mean(parallax)) + " +/- " + \
      ("%09.8f" % numpy.std(parallax)))
dist         = 1000.0/numpy.mean(parallax)
plus1sigdist = 1000.0/(numpy.mean(parallax) - numpy.std(parallax))
min1sigdist  = 1000.0/(numpy.mean(parallax) + numpy.std(parallax))
print("dist  = %6.1f + %5.1f - %5.1f pc" % (dist, plus1sigdist-dist, dist-min1sigdist))
if doorbitalfit:
    print("Omega = " + str(numpy.mean(Omega)) + " +/- " + \
      ("%09.8f" % numpy.std(Omega)))
    print("incl. = " + str(numpy.mean(inclination)) + " +/- " + \
      ("%09.8f" % numpy.std(inclination)))
    print("folded incl. = " + str(numpy.mean(foldedinclination)) + " +/- " + \
      ("%09.8f" % numpy.std(foldedinclination)))
    print("xdot. = " + str(numpy.mean(xdot)) + " +/- " + \
      ("%.2E" % numpy.std(xdot)))

print("===========================\n")

print("  BOOTSTRAP VALUES (HIST)  ")
print("===========================")
#centre_mura = (numpy.median(pm_ra) + numpy.mean(pm_ra))/2.0
#centre_mudec = (numpy.median(pm_dec) + numpy.mean(pm_dec))/2.0
#centre_parallax = (numpy.median(parallax) + numpy.mean(parallax))/2.0
centre_mura = pmra_bracket[2]
centre_mudec = pmdec_bracket[2]
centre_parallax = parallax_bracket[2]
#print inclination_bracket
#print sorted(inclination)
if doorbitalfit:
    centre_Omega = Omega_bracket[2]
    centre_inc = inclination_bracket[2]
    centre_foldedinc = foldedinclination_bracket[2]
    centre_xdot = xdot_bracket[2]
    if centre_Omega < Omega_bracket[0]:
        centre_Omega = Omega_bracket[0]
    if centre_Omega > Omega_bracket[1]:
        centre_Omega = Omega_bracket[1]
    if centre_inc < inclination_bracket[0]:
        centre_inc = inclination_bracket[0]
    if centre_inc > inclination_bracket[1]:
        centre_inc = inclination_bracket[1]
    if centre_foldedinc < foldedinclination_bracket[0]:
        centre_foldedinc = foldedinclination_bracket[0]
    if centre_foldedinc > foldedinclination_bracket[1]:
        centre_foldedinc = foldedinclination_bracket[1]
    if centre_xdot < xdot_bracket[0]:
        centre_xdot = xdot_bracket[0]
    if centre_xdot > xdot_bracket[1]:
        centre_xdot = xdot_bracket[1]

print("mu_a  = %f + %09.8f - %09.8f" % \
      (centre_mura, pmra_bracket[1] - centre_mura,
       centre_mura - pmra_bracket[0]))
print("mu_d  = %f + %09.8f - %09.8f" % \
      (centre_mudec, pmdec_bracket[1] - centre_mudec,
       centre_mudec - pmdec_bracket[0]))
print("pi    = %f + %09.8f - %09.8f" % \
      (centre_parallax, parallax_bracket[1] - centre_parallax,
       centre_parallax - parallax_bracket[0]))
if doorbitalfit:
    print("Omega  = %f + %09.8f - %09.8f" % \
          (centre_Omega, Omega_bracket[1] - centre_Omega,
           centre_Omega - Omega_bracket[0]))
    print("inc   = %f + %09.8f - %09.8f" % \
          (centre_inc, inclination_bracket[1] - centre_inc,
           centre_inc - inclination_bracket[0]))
    print("folded inc = %f + %09.8f - %09.8f" % \
          (centre_foldedinc, foldedinclination_bracket[1] - centre_foldedinc,
           centre_foldedinc - foldedinclination_bracket[0]))
    print("xdot   = %.2E + %.2E - %.2E" % \
          (centre_xdot, xdot_bracket[1] - centre_xdot,
           centre_xdot - xdot_bracket[0]))
    orbitparamlines = open(options.orbitfitfile).readlines()
    timingxdotval = ""
    timingxdotsig = ""
    for line in orbitparamlines:
        if "xdot:" in line and not line[0] == "#":
            timingxdotval = line.split(':')[-1]
        if "xdot_sigma:" in line and not line[0] == "#":
            timingxdotsig = line.split(':')[-1]
    if not timingxdotval == "" and not timingxdotsig == "":
        print("Timing xdot value is %s +/- %s" % (timingxdotval, timingxdotsig))

# Plot the results
if doorbitalfit:
    print("Not doing fancy plot for binary systems yet")
    sys.exit()

docolour = False
maxplottrials = 1000
font_size = 16
plottype = 'pdf'
legendlocation="best"
fitlinestyle = 'k-'
fitpointsstyle = 'kx'
bothfitstyle = 'kx-'
fitpointmarkersize = 5
if docolour: 
    fitlinestyle = 'b-'
    fitpointsstyle = 'bx'
    bothfitstyle = 'bx-'
errorlinecolour ='k'
if docolour:
    errorlinecolour = 'r'
measurementlinestyle = 'ko'
if docolour:
    measurementlinestyle = 'ko'
pyplot.clf()
pyplot.rc("axes", linewidth=2.0)
pyplot.rc("lines", markeredgewidth=2.0)  # this also grows the points/lines in the plot..
pyplot.figure(1)
ax  = pyplot.subplot(111)
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(font_size*0.9)
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(font_size*0.9)
pyplot.xlabel("Time [MJD]", fontsize=font_size)
pyplot.ylabel("R.A. offset [mas]", fontsize=font_size)
plotfile1 = "ratime.bootstrap." + plottype
if not pulsar == "":
    plotfile1 = "ratime.bootstrap." + pulsar + "." + plottype
pyplot.figure(2)
ax  = pyplot.subplot(111)
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(font_size*0.9)
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(font_size*0.9)
pyplot.xlabel("Time [MJD]", fontsize=font_size)
pyplot.ylabel("Decl. offset [mas]", fontsize=font_size)
plotfile2 = "dectime.bootstrap." + plottype
if not pulsar == "":
    plotfile2 = "dectime.bootstrap." + pulsar + "." + plottype

junkfile = 'bootstrapjunk.txt'
todo = maxplottrials
alphaval = 10.0/maxplottrials
meanrarad = numpy.mean(ref_ra_hour)*math.pi/12
meandecrad = numpy.mean(ref_dec_deg)*math.pi/180
rastr, decstr = astro_utils.posradians2string(meanrarad, meandecrad)
maxrapx = []
maxdecpx = []
if numtrials < maxplottrials:
    todo = numtrials
    alphaval = 10.0/numtrials
for  i in range(todo):
    if i%(todo/10) == 0:
        print(str(i/(todo/10)) + "0% done")

    tempout = open(junkfile, 'w')
    #rastr, decstr = astro_utils.posradians2string(ref_ra_hour[i]*math.pi/12, ref_dec_deg[i]*math.pi/180)
    tempout.writelines(pmparlines[0:headerendlinenum])
    tempout.write("pi = %f\n" % parallax[i])
    tempout.write("mu_a = %f\n" % pm_ra[i])
    tempout.write("mu_d = %f\n" % pm_dec[i])
    tempout.write("ra = %s\n" % rastr)
    tempout.write("dec = %s\n" % decstr)
    tempout.write(astrometric_points[0] + '\n')
    tempout.write(astrometric_points[-1] + '\n')
    tempout.close()
    ttimes, tras, tdecs, etimes, eras, edecs, eraerr, edecerr, pras, pdecs, fitlines = pmpar(pmpar_exe, junkfile, subtractpm)
    tras += ((ref_ra_hour[i]*math.pi/12 - meanrarad)*180*3.6e6*math.cos(meandecrad))/math.pi
    tdecs += (ref_dec_deg[i] - meandecrad*180/math.pi)*3.6e6
    pyplot.figure(1)
    fitline = pyplot.plot(ttimes,tras,fitlinestyle,alpha=alphaval)
    pyplot.figure(2)
    fitline = pyplot.plot(ttimes,tdecs,fitlinestyle,alpha=alphaval)
    maxrapx.append(max(tras))
    maxdecpx.append(max(tdecs))
    #fitpoints = pyplot.plot(etimes,pras,fitpointsstyle, markersize=fitpointmarkersize)
# Plot with a label once for the legend
pyplot.figure(1)
junk = pyplot.plot([],[],bothfitstyle,markersize=fitpointmarkersize,label="Astrometric fits")
pyplot.figure(2)
junk = pyplot.plot([],[],bothfitstyle,markersize=fitpointmarkersize,label="Astrometric fits")

# Now the actual measurements
tempout = open(junkfile, 'w')
#rastr, decstr = astro_utils.posradians2string(ref_ra_hour[i]*math.pi/12, ref_dec_deg[i]*math.pi/180)
tempout.writelines(pmparlines[0:headerendlinenum])
tempout.write("pi = %f\n" % numpy.mean(parallax))
tempout.write("mu_a = %f\n" % numpy.mean(pm_ra))
tempout.write("mu_d = %f\n" % numpy.mean(pm_dec))
rastr, decstr = astro_utils.posradians2string(numpy.mean(ref_ra_hour)*math.pi/12, numpy.mean(ref_dec_deg)*math.pi/180)
tempout.write("ra = %s\n" % rastr)
tempout.write("dec = %s\n" % decstr)
for i in range(len(astrometric_points)):
    tempout.write(astrometric_points[i] + '\n')
tempout.close()
ttimes, tras, tdecs, etimes, eras, edecs, eraerr, edecerr, pras, pdecs, fitlines = pmpar(pmpar_exe, junkfile, subtractpm)

# Calculate the yrange for RA and dec, leave some space at the top for legend, round up to nearest 0.2 mas
maxrapx.sort()
maxdecpx.sort()

ylimra = 2*maxrapx[int(len(maxrapx)*0.5)] + 0.1
ylimdec = 2*maxdecpx[int(len(maxdecpx)*0.5)] + 0.1
if ylimra < max(eras):
    ylimra = max(eras)*1.2
if ylimdec < max(edecs):
    ylimdec = max(edecs)*1.2

ylimra = float(int(5*ylimra+1))/5.0
ylimdec = float(int(5*ylimdec+1))/5.0

pyplot.figure(1)
pyplot.ylim(-ylimra, ylimra)
measurementline = pyplot.errorbar(etimes, eras, yerr=eraerr, fmt=measurementlinestyle, markersize=3, capsize=3, markeredgewidth=1.5, label="Measured positions")
pyplot.legend(numpoints=1, loc=legendlocation)
pyplot.tight_layout()
pyplot.savefig(plotfile1)
pyplot.figure(2)
pyplot.ylim(-ylimdec, ylimdec)
measurementline = pyplot.errorbar(etimes, edecs, yerr=edecerr, fmt=measurementlinestyle, markersize=3, capsize=3, markeredgewidth=1.5, label="Measured positions")
pyplot.legend(numpoints=1, loc=legendlocation)
pyplot.tight_layout()
pyplot.savefig(plotfile2)

#end======================================
##
###