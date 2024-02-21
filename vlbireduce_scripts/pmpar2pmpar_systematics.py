#!/usr/bin/python3
# Adam Deller


import os, sys, re, subprocess, math
import numpy
from optparse import OptionParser

## GLOBAL VARIABLES #####################################
astrometric_points = []
numpoints = 0
numepochs = 0
lastmjd = 0
mjds = []
ra = []
dec = []
ra_err = []
dec_err = []
wmean_ra = []
wmean_dec = []
wmean_ra_errs = []
wmean_dec_errs = []
ra_systematic_err = []
dec_systematic_err = []
leadingras = []
leadingdecs = []
epoch_points = []

## SUBROUTINE ###########################################
def run_syserr_pmpar_file(filename, syserrra, syserrdec):
    tempout = open(filename, 'w')
    tempout.writelines(pmparlines[0:headerendlinenum])
    if not haveepoch:
        tempout.write("epoch = 2007.0\n\n")
    for i in range(numepochs):
        tempout.write(mjds[i] + " " + leadingras[i] + str(wmean_ra[i]) + \
        " " + str(math.sqrt(wmean_ra_errs[i]*wmean_ra_errs[i] + \
        ra_systematic_err[i]*ra_systematic_err[i] + syserrra*syserrra)) + \
        " " + leadingdecs[i] + str(wmean_dec[i]) + " " + \
        str(math.sqrt(wmean_dec_errs[i]*wmean_dec_errs[i] + \
        dec_systematic_err[i]*dec_systematic_err[i] + \
        syserrdec*syserrdec)) + "\n")
    tempout.close()
    output = subprocess.getoutput(pmpar_exe + filename).split('\n')
    return output

## MAIN CODE ###########################################
usage = "usage: %prog -p <pulsar name> [options]"
pmpar_exe = 'pmpar '
parser = OptionParser(usage)

parser.add_option("-d", "--directory", dest="directory",
                  default="/export/home/marathon2/data/1023+0038/solutions/",
                  help="pulsar name eg 0437-4715")
parser.add_option("-p", "--pulsar", dest="pulsar", default="1023+0038",
                  help="pulsar name eg 0437-4715", metavar="PULSAR")
parser.add_option("-f", "--filename", dest="filename", default="",
                  help ="Filename (default <pulsar>.jmfit.pmpar.in)")
parser.add_option("--nointraepoch", dest="nointraepoch", default=False,
                  action="store_true", help ="Kill the intra-epoch estimate")
parser.add_option("--scalefactor", dest="scalefactor", default="1.0",
                  help="Scale all errors by this factor")
parser.add_option("--sysfloor", dest="sysfloor", default="0.0",
                  help="Add this many microarcsecond as a floor error")
parser.add_option("--removelineartrend", dest="removelineartrend", default="",
                  help ="Remove linear trend between bands, F for fit, " + \
                        "otherwise give raslope[,decslope] in radians")

(options, args) = parser.parse_args()
if len(args) > 0:
     parser.error("You can only supply the permitted options!")

pulsar        = options.pulsar
directory     = options.directory
no_intraepoch = options.nointraepoch
fitlinear     = False
removelinear  = False
sysfloor      = float(options.sysfloor)
scalefactor   = float(options.scalefactor)
raslope       = 0.0
decslope      = 0.0
if options.removelineartrend != "":
    removelinear = True
    if options.removelineartrend == "F":
        fitlinear = True
    else:
        trendsplit = options.removelineartrend.split(',')
        if len(trendsplit) > 2:
            parser.error("--removelineartrend=F or --removelineartrend=raslope[,decslope]")
        raslope = float(trendsplit[0])
        if len(trendsplit) == 2:
            decslope = float(trendsplit[1])
print(sysfloor)
pmparfile = directory + pulsar + ".jmfit.pmpar.in"
if options.filename != "":
    pmparfile = options.filename
if not os.path.exists(pmparfile):
    print("Could not find file " + pmparfile + " - aborting!!!")
    sys.exit(1)

# Parse in pmpar file
pmparin = open(pmparfile, 'r')
pmparlines = pmparin.readlines()
pmparin.close()

headerendlinenum = 0
haveepoch = False
numkeywords = 0
while re.search('=', pmparlines[headerendlinenum]) or \
      pmparlines[headerendlinenum].rstrip() == "" or \
      pmparlines[headerendlinenum].lstrip()[0] == '#':
    if (re.search('epoch', pmparlines[headerendlinenum]) or \
       re.search('Epoch', pmparlines[headerendlinenum]) or \
       re.search('Epoch', pmparlines[headerendlinenum])) and \
       (not pmparlines[headerendlinenum].lstrip()[0] == '#') :
        haveepoch = True
    if re.search('=', pmparlines[headerendlinenum]):
        numkeywords = numkeywords + 1
    print("Skipping line " + pmparlines[headerendlinenum][0:-1])
    headerendlinenum = headerendlinenum + 1
if not haveepoch:
    numkeywords = numkeywords+1
startline = numkeywords

for line in pmparlines[headerendlinenum:]:
    if not line.lstrip() == "" and not line.lstrip()[0] == '#':
        astrometric_points.append(line)
        ra.append([])
        dec.append([])
        ra_err.append([])
        dec_err.append([])
        numpoints = numpoints + 1
        splitline = line.split()
        mjd = float(splitline[0])
        if abs(mjd-lastmjd) > 1.0 or lastmjd < 50000:
            numepochs = numepochs + 1
            lastmjd = mjd
            mjds.append(splitline[0])
            wmean_ra.append(0.0)
            wmean_dec.append(0.0)
            wmean_ra_errs.append(0.0)
            wmean_dec_errs.append(0.0)
            ra_systematic_err.append(0.0)
            dec_systematic_err.append(0.0)
            epoch_points.append(0)
            leadingras.append(splitline[1].split(':')[0] + ":" + \
                              splitline[1].split(':')[1] + ":")
            leadingdecs.append(splitline[3].split(':')[0] + ":" + \
                               splitline[3].split(':')[1] + ":")
            decradians = (math.fabs(float(splitline[3].split(':')[0])) + \
                          float(splitline[3].split(':')[1])/60.0)*\
                          math.pi/180.0
        ra[numepochs-1].append(float(splitline[1].split(':')[2]))
        dec[numepochs-1].append(float(splitline[3].split(':')[2]))
        ra_err[numepochs-1].append(float(splitline[2])*scalefactor)
        dec_err[numepochs-1].append(float(splitline[4])*scalefactor)
        epoch_points[numepochs-1] = epoch_points[numepochs-1] + 1

#Work out the systematic contribution in ra and dec for each epoch
averagesyserror = 0.0
average_wmean_error_ra = 0.0
average_wmean_error_dec = 0.0
averagesyserror_ra = 0.0
averagesyserror_dec = 0.0
for i in range(numepochs):
    sigmax_sum = 0.0
    sigmay_sum = 0.0

    if removelinear:
        if fitlinear:
            # Make a fit to any linear trend in RA and dec across the bands
            ramodel = numpy.polyfit(list(range(epoch_points[i])), ra[i], 1)
            decmodel = numpy.polyfit(list(range(epoch_points[i])), dec[i], 1)
            rapredictor = numpy.poly1d(ramodel)
            decpredictor = numpy.poly1d(decmodel)
            raslope = rapredictor[1]
            decslope = decpredictor[1]
            print("RA predictor: " + str(ramodel))
            print("Dec predictor: " + str(decmodel))
        for j in range(epoch_points[i]):
            ra[i][j] += (epoch_points[i]-j)*raslope
            dec[i][j] += (epoch_points[i]-j)*decslope

    #First get the weighted mean of ra and dec for this epoch
    for j in range(epoch_points[i]):
        wmean_ra[i] = wmean_ra[i] + ra[i][j]/(ra_err[i][j]*ra_err[i][j])
        wmean_dec[i] = wmean_dec[i] + dec[i][j]/(dec_err[i][j]*dec_err[i][j])
        sigmax_sum = sigmax_sum + 1.0/(ra_err[i][j]*ra_err[i][j])
        sigmay_sum = sigmay_sum + 1.0/(dec_err[i][j]*dec_err[i][j])
    wmean_ra[i] = wmean_ra[i] / sigmax_sum
    wmean_dec[i] = wmean_dec[i] / sigmay_sum
    wmean_ra_errs[i] = 1.0/math.sqrt(sigmax_sum)
    wmean_dec_errs[i] = 1.0/math.sqrt(sigmay_sum)
    print("weighted mean position for epoch " + str(i) + " is " + \
          str(wmean_ra[i]) + " " + str(wmean_dec[i]))
    print("weighted mean errors for epoch " + str(i) + " are " + \
          str(wmean_ra_errs[i]) + " " + str(wmean_dec_errs[i]))
    average_wmean_error_ra = average_wmean_error_ra + \
                             15.0*math.cos(decradians)*\
                             wmean_ra_errs[i]/numepochs
    average_wmean_error_dec = average_wmean_error_dec + \
                              wmean_dec_errs[i]/numepochs
   #Then iterate towards systematic  std dev in ra and dec for this epoch
    average_ra_stddevs = 0.0
    average_dec_stddevs = 0.0
    average_ra_dist = 0.0
    average_dec_dist = 0.0
    for j in range(epoch_points[i]):
        average_ra_stddevs = average_ra_stddevs + \
                             math.fabs(ra[i][j] - \
                             wmean_ra[i])/(wmean_ra_errs[i]*epoch_points[i])
        average_dec_stddevs = average_dec_stddevs + math.fabs(dec[i][j] - \
                              wmean_dec[i])/(wmean_dec_errs[i]*epoch_points[i])
        average_ra_dist = average_ra_dist + math.fabs(ra[i][j] - \
                          wmean_ra[i])/(epoch_points[i])
        average_dec_dist = average_dec_dist + math.fabs(dec[i][j] - \
                           wmean_dec[i])/(epoch_points[i])

    if average_ra_stddevs > 1.0:
        while math.fabs(average_ra_stddevs - 1.0) > 0.0000001:
            ra_systematic_err[i] = ra_systematic_err[i] + \
                                   (average_ra_stddevs - 1.0)*\
                                   (0.5*average_ra_dist)
            average_ra_stddevs = 0.0
            for j in range(epoch_points[i]):
                average_ra_stddevs = average_ra_stddevs + math.fabs(ra[i][j] \
                                     - wmean_ra[i])/\
                                     (math.sqrt(wmean_ra_errs[i]*\
                                     wmean_ra_errs[i] + \
                                     ra_systematic_err[i]*\
                                     ra_systematic_err[i])*epoch_points[i])
    if average_dec_stddevs > 1.0:
        while math.fabs(average_dec_stddevs - 1.0) > 0.0000001:
            dec_systematic_err[i] = dec_systematic_err[i] + \
                                    (average_dec_stddevs - 1.0)*\
                                    (0.5*average_dec_dist)
            average_dec_stddevs = 0.0
            for j in range(epoch_points[i]):
                average_dec_stddevs = average_dec_stddevs + \
                                      math.fabs(dec[i][j] - \
                                      wmean_dec[i])/\
                                      (math.sqrt(wmean_dec_errs[i]*\
                                      wmean_dec_errs[i] + \
                                      dec_systematic_err[i]*\
                                      dec_systematic_err[i])*epoch_points[i])
    print("systematic errors for epoch " + str(i) + " are " + \
          str(15.0*math.cos(decradians)*ra_systematic_err[i]) + " " + \
          str(dec_systematic_err[i]))
    averagesyserror = averagesyserror + \
                      math.sqrt(225.0*math.cos(decradians)*\
                                math.cos(decradians)*ra_systematic_err[i]*\
                                ra_systematic_err[i] + dec_systematic_err[i]*\
                                dec_systematic_err[i])/numepochs
    averagesyserror_ra = averagesyserror_ra + 15.0*math.cos(decradians)*\
                         ra_systematic_err[i]/numepochs
    averagesyserror_dec = averagesyserror_dec + dec_systematic_err[i]/numepochs
    if no_intraepoch:
        print("Resetting intra-epoch errors to zero")
        ra_systematic_err[i] = 0
        dec_systematic_err[i] = 0

print("Total average RA intra-epoch systematic error is " + \
      str(averagesyserror_ra))
print("Total average Dec intra-epoch systematic error is " + \
      str(averagesyserror_dec))
print("Total average RA weighted error is " + str(average_wmean_error_ra))
print("Total average Dec weighted error is " + str(average_wmean_error_dec))
print("Decradians was " + str(decradians))

#Write out a file for the original
tempout = open('junk.txt', 'w')
tempout.writelines(pmparlines[0:headerendlinenum])
if not haveepoch:
    tempout.write("epoch = 55000.0\n\n")
tempout.writelines(pmparlines[headerendlinenum:])

#Write out a file for the original plus estimated systematic 
#errors at each epoch
tempout = open('full_indbands.txt', 'w')
tempout.writelines(pmparlines[0:headerendlinenum])
if not haveepoch:
    tempout.write("epoch = 55000.0\n\n")
for i in range(numepochs):
    for j in range(epoch_points[i]):
        tempout.write(mjds[i] + " " + leadingras[i] + str(ra[i][j]) + " " + \
                      str(math.sqrt(ra_err[i][j]*ra_err[i][j] + \
                      ra_systematic_err[i]*ra_systematic_err[i])) + " " + \
                      leadingdecs[i] + str(dec[i][j]) + " " + \
                      str(math.sqrt(dec_err[i][j]*dec_err[i][j] + \
                      dec_systematic_err[i]*dec_systematic_err[i])) + "\n")

#Write out a file for the weighted mean epoch points and errors
tempout = open('wmean.txt', 'w')
tempout.writelines(pmparlines[0:headerendlinenum])
if not haveepoch:
    tempout.write("epoch = 55000.0\n\n")
for i in range(numepochs):
    tempout.write(mjds[i] + " " + leadingras[i] + str(wmean_ra[i]) + " " + \
                  str(wmean_ra_errs[i]) + " " + leadingdecs[i] + \
                  str(wmean_dec[i]) + " " + str(wmean_dec_errs[i]) + "\n")
tempout.close()

# If specified, do a wmean+sysfloor file
if sysfloor > 0.0:
    sysfloor /= 1e6
    sysfloorraerr = sysfloor/(15.0*math.cos(decradians))
    print("sysfloorerr is %f, sysfloorraerr is %f" % (sysfloor, sysfloorraerr))
    tempout = open('meanplusfloor.txt', 'w')
    tempout.writelines(pmparlines[0:headerendlinenum])
    if not haveepoch:
        tempout.write("epoch = 55000.0\n\n")
    for i in range(numepochs):
        raerr = math.sqrt(wmean_ra_errs[i]*wmean_ra_errs[i] + sysfloorraerr*sysfloorraerr)
        decerr = math.sqrt(wmean_dec_errs[i]*wmean_dec_errs[i] + sysfloor*sysfloor)
        tempout.write(mjds[i] + " " + leadingras[i] + str(wmean_ra[i]) + " " + \
                      str(raerr) + " " + leadingdecs[i] + \
                      str(wmean_dec[i]) + " " + str(decerr) + "\n")
    tempout.close()


#Write out a file for the same plus the intra-epoch systematic 
#errors as estimated above
tempout = open('full.txt', 'w')
tempout.writelines(pmparlines[0:headerendlinenum])
if not haveepoch:
    tempout.write("epoch = 55000.0\n\n")
for i in range(numepochs):
    tempout.write(mjds[i] + " " + leadingras[i] + str(wmean_ra[i]) + " " + \
                  str(math.sqrt(wmean_ra_errs[i]*wmean_ra_errs[i] + \
                  ra_systematic_err[i]*ra_systematic_err[i])) + " " + \
                  leadingdecs[i] + str(wmean_dec[i]) + " " + \
                  str(math.sqrt(wmean_dec_errs[i]*wmean_dec_errs[i] + \
                  dec_systematic_err[i]*dec_systematic_err[i])) + "\n")
tempout.close()

#Print comparison to screen
startline = 1 #set startline to appropriate value for pmpar output

#Original first
output = subprocess.getoutput(pmpar_exe + "junk.txt").split('\n')
print("===========================")
print("PMPAR (FULL DATASET) VALUES")
print("===========================")
for line in output:
    if not line == "" and not line[0] == "#":
        print(line)

#Then with the systematic errors added to the band estimates
output = subprocess.getoutput(pmpar_exe + " full_indbands.txt").split('\n')
print("========================================================")
print("PMPAR (INDIVIDUAL BANDS + INTRA EPOCH SYSTEMATIC) VALUES")
print("========================================================")
for line in output:
    if not line == "" and not line[0] == "#":
        print(line)

#Then weighted mean - shoud be the same as the original dataset
output = subprocess.getoutput(pmpar_exe + " wmean.txt").split('\n')
print("====================================")
print("PMPAR (WEIGHTED MEAN DATASET) VALUES")
print("====================================")
for line in output:
    if not line == "" and not line[0] == "#":
        print(line)

#Then with systematic estimate
output = subprocess.getoutput(pmpar_exe + " full.txt").split('\n')
print("=====================================================")
print("PMPAR (WEIGHTED MEAN + INTRA EPOCH SYSTEMATIC) VALUES")
print("=====================================================")
for line in output:
    if not line == "" and not line[0] == "#":
        print(line)

#Now iteratively calculate a single estimate of inter-epoch systematic error, 
#add and run until a chi-squared of one is obtained
rchisq = float(output[-2].split()[3])
if rchisq < 1.0:
     print("Reduced chi--squared is already < 1.0!! No inter-epoch " + \
           "systematic will be added")
     sys.exit(0)

intersys = 0.0
if averagesyserror == 0:
    averagesyserror = 0.001
while math.fabs(rchisq - 1.0) > 0.0001:
    intersys = intersys + (rchisq - 1.0)*averagesyserror*0.05
    intersysra = intersys/(15.0*math.cos(decradians))
    tempout = open('junk.txt', 'w')
    tempout.writelines(pmparlines[0:headerendlinenum])
    if not haveepoch:
        tempout.write("epoch = 55000.0\n\n")
    for i in range(numepochs):
        tempout.write(mjds[i] + " " + leadingras[i] + str(wmean_ra[i]) + \
                      " " + str(math.sqrt(wmean_ra_errs[i]*wmean_ra_errs[i] + \
                      ra_systematic_err[i]*ra_systematic_err[i] + \
                      intersysra*intersysra)) + " " + leadingdecs[i] + \
                      str(wmean_dec[i]) + " " + \
                      str(math.sqrt(wmean_dec_errs[i]*wmean_dec_errs[i] + \
                      dec_systematic_err[i]*dec_systematic_err[i] + \
                      intersys*intersys)) + "\n")
    tempout.close()
    output = subprocess.getoutput(pmpar_exe + " junk.txt").split('\n')
    rchisq = float(output[-2].split()[3])
    print("rchisq was " + str(rchisq))

#Print estimate with inter-epoch systematic error apportioned evenly between RA and dec
print("Estimated inter-epoch systematic error is " + str(intersys*1000.0) + \
      " mas for both RA and dec, thus total inter-epoch systematic error is " + str(math.sqrt(2.0)*intersys*1000.0) + \
      " mas")
print("============================================================================================")
print("PMPAR (WEIGHTED MEAN + INTRA EPOCH SYSTEMATIC + ESTIMATED INTER EPOCH (EQ RA AND DEC) VALUES")
print("============================================================================================")
for line in output:
    if not line == "" and not line[0] == "#":
        print(line)

os.system("mv -f junk.txt wmeanplussys.txt")

#Work out minimal sysra^2 + sysdec^2 by allowing each to float but keeping rchi^2 at 1
deltaangle = 20.0*math.pi/180.0
angle = 45.0*math.pi/180.0
amp = math.sqrt(2.0)*intersys
while deltaangle > 0.1*math.pi/180.0 and angle > 0.0 and angle < math.pi/2.0:
    orgamp = amp
    print("Trying around angle " + str(int(angle*180.0/math.pi + 0.5)))
    #try increasing angle
    cur_angle = angle + deltaangle
    if cur_angle > math.pi/2.0:
        cur_angle = math.pi/2.0
    cur_rasys = amp*math.cos(cur_angle)/(15.0*math.cos(decradians))
    cur_decsys = amp*math.sin(cur_angle)
    output = run_syserr_pmpar_file('junk.txt',cur_rasys,cur_decsys)
    rchisq = float(output[-2].split()[3])
    print("First guess at angle " + str(180.0*cur_angle/math.pi) + " degrees has rchisquared=" + str(rchisq))
    if rchisq < 1.0:
        while math.fabs(rchisq - 1.0) > 0.0001:
            amp = 0.5*amp + 0.5*amp*math.sqrt(rchisq)
            cur_rasys = amp*math.cos(cur_angle)/(15.0*math.cos(decradians))
            cur_decsys = amp*math.sin(cur_angle)
            output = run_syserr_pmpar_file('junk.txt',cur_rasys,cur_decsys)
            rchisq = float(output[-2].split()[3])
    posamp = amp
    #Try decreasing angle
    amp = orgamp
    cur_angle = angle - deltaangle
    if cur_angle < 0.0:
        cur_angle = 0.0
    cur_rasys = amp*math.cos(cur_angle)/(15.0*math.cos(decradians))
    cur_decsys = amp*math.sin(cur_angle)
    output = run_syserr_pmpar_file('junk.txt',cur_rasys,cur_decsys)
    rchisq = float(output[-2].split()[3])
    print("First guess at angle " + str(180.0*cur_angle/math.pi) + " degrees has rchisquared=" + str(rchisq))
    if rchisq < 1.0:
        while math.fabs(rchisq - 1.0) > 0.0001:
            amp = 0.5*amp + 0.5*amp*math.sqrt(rchisq)
            cur_rasys = amp*math.cos(cur_angle)/(15.0*math.cos(decradians))
            cur_decsys = amp*math.sin(cur_angle)
            output = run_syserr_pmpar_file('junk.txt',cur_rasys,cur_decsys)
            rchisq = float(output[-2].split()[3])
    negamp = amp
    if posamp < negamp:
        angle = angle+deltaangle
        amp = posamp
    elif negamp < posamp:
        angle = angle-deltaangle
        amp = negamp
    else:
        deltaangle = deltaangle*2.0/3.0
        amp = orgamp

while math.fabs(rchisq - 1.0) > 0.0001:
    amp = 0.5*amp + 0.5*amp*math.sqrt(rchisq)
    cur_rasys = amp*math.cos(cur_angle)/(15.0*math.cos(decradians))
    cur_decsys = amp*math.sin(cur_angle)
    output = run_syserr_pmpar_file('junk.txt',cur_rasys,cur_decsys)
    rchisq = float(output[-2].split()[3])

#Print final best estimate with inter-epoch systematic error variable 
#between RA and dec
print("Estimated inter-epoch systematic error is " + \
      str(cur_rasys*1000.0*15.0*math.cos(decradians)) + " for RA and " + \
      str(cur_decsys*1000.0) + " mas for dec, thus total inter-epoch " + \
      "systematic error is " + str(math.sqrt(cur_rasys*cur_rasys*15.0*\
      math.cos(decradians)*15.0*math.cos(decradians) + \
      cur_decsys*cur_decsys)*1000.0) + " mas")
output = run_syserr_pmpar_file('junk.txt',cur_rasys,cur_decsys)
print("================================================================" + \
      "============================")
print("PMPAR (WEIGHTED MEAN + INTRA EPOCH SYSTEMATIC + ESTIMATED INTER " + \
      "EPOCH (EQ RA AND DEC) VALUES")
print("================================================================" + \
      "============================")
for line in output:
    if not line == "" and not line[0] == "#":
        print(line)
