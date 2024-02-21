"""
Ashish Kalyan, Hao Ding, and Adam Deller

The script assembles the AIPS tasks and other functions like primary beam correction
and pulsar scintillation correction (taken from Hao's GitHub repository),
https://github.com/ashishkalyan/psrvlbireduce/tree/vlbireduce_python3
AIPSImage(imgdata.name, beamklass, imgdata.disk, imgdata.seq)
"""

from scipy.special import jn # Bessel function
import math, os, sys, ephem, datetime, subprocess, yaml
import AIPSTV
from AIPSTask import AIPSTask
from AIPSData import AIPSUVData, AIPSImage
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
import vlbi_calibration_support as vcs
import astro_utils as ast


yamldir = os.path.dirname(os.path.abspath(__file__)) + "/examples"
expconfig  = yaml.load(open(yamldir + '/expconfig.yaml'), yaml.SafeLoader)
aipsver = expconfig['aipsver']
print("You are using AIPS version: ", aipsver)


vlbadiameter = 25.47415721 # metres
gbtdiameter = 100 # metres
vlbabeams = {}
# Squint(R-L) Az, Squint(R-L) EL, RCP Az beamwidth, RCP EL beamwidth, LCP Az beamwidth, LCP El beamwidth, reffreq
vlbabeams['SC'] = [-1.53, -0.56, 29.58, 29.57, 29.61, 29.36, 1438.0e6]
vlbabeams['HN'] = [-1.35, -0.67, 29.17, 29.21, 29.26, 29.09, 1438.0e6]
vlbabeams['NL'] = [-1.56, -0.62, 29.44, 29.55, 29.51, 29.47, 1438.0e6]
vlbabeams['FD'] = [-1.56, -0.58, 29.27, 29.24, 29.41, 29.27, 1438.0e6]
vlbabeams['LA'] = [-1.62, -0.67, 29.12, 29.03, 29.20, 29.00, 1438.0e6]
vlbabeams['PT'] = [-1.59, -0.63, 29.27, 29.18, 29.28, 29.14, 1438.0e6]
vlbabeams['KP'] = [-1.61, -0.59, 29.30, 29.32, 29.19, 29.04, 1438.0e6]
vlbabeams['OV'] = [-1.76, -0.95, 29.57, 29.64, 29.67, 29.62, 1438.0e6]
vlbabeams['BR'] = [-1.44, -0.50, 29.19, 29.36, 29.55, 29.44, 1438.0e6]
vlbabeams['MK'] = [-1.56, -0.60, 29.20, 29.26, 29.24, 29.25, 1438.0e6]
vlbabeams['GB'] = [0, 0, 7.3, 7.3, 7.3, 7.3, 1438.0e6] # Total guess, based on scaling from VLBA

vlbalocations = {}
# Latitude, longitude, elevation(m)
vlbalocations['SC'] = ['17:45:23.68', '-64:35:01.07',  16  ]
vlbalocations['HN'] = ['42:56:00.99', '-71:59:11.69',  296 ]
vlbalocations['NL'] = ['41:46:17.13', '-91:34:26.88',  222 ]
vlbalocations['FD'] = ['30:38:06.11', '-103:56:41.34', 1606]
vlbalocations['LA'] = ['35:46:30.45', '-106:14:44.15', 1962]
vlbalocations['PT'] = ['34:18:03.61', '-108:07:09.06', 2365]
vlbalocations['KP'] = ['31:57:22.70', '-111:36:44.72', 1902]
vlbalocations['OV'] = ['37:13:53.95', '-118:16:37.37', 1196]
vlbalocations['BR'] = ['48:07:52.42', '-119:40:59.80', 250 ]
vlbalocations['MK'] = ['19:48:04.97', '-155:27:19.81', 3763]
vlbalocations['GB'] = ['38:26:16.16', '-100:09:51.20', 844]


def delAIPSCat (indata):
    """
    ZAP: Delete AIPSCat entry
    """
    if indata.exists():
        indata.zap()
        print(f"VLBI PIPE: delAIPSCat: {indata} deleted")
    else:
        print(f"VLBI PIPE: delAIPSCat: {indata} does not found in AIPSCat")
            
            
def deletetable(uvdataset, tabletype, tableversion):
    """
    ZAP: Delete an extension table of input uvdata
    """
    try:
        uvdataset.table(tabletype, tableversion).zap()
        print(f'VLBI PIPE: deletetable: {tabletype}{tableversion} table deleted')
        return True
    except IOError:
        print(f"VLBI PIPE: deletetable: Apparently {tabletype} version {tableversion} didn't exist...")
        return False
    

def fitld_vlba (fits_file, outdata, clint=0.166667, wtthreshhold=0.05):
    """
    FITLD: Load only VLBA uvdata
    """
    if os.path.exists(fits_file):
        fitld = AIPSTask('fitld', version=aipsver)
        fitld.datain     = fits_file
        fitld.outdata    = outdata
        fitld.douvcomp   = 1
        fitld.digicor    = 1
        fitld.clint      = clint
        fitld.wtthresh   = wtthreshhold
        # fitld.refdate    = refdate
        fitld.antname[1] = 'VLBA'
        if outdata.exists():
            if vcs.yesno(f"Delete existing UV dataset {outdata.name} ? (No will abort pipeline)"):
                outdata.zap()
            else: sys.exit()
        fitld.inp()
        fitld()
    else: print(f"VLBI PIPE: fitld_vlba: {fits_file} file does not found")


def fitld_image (fits_file, outdata):
    """
    FITLD: Load image
    """
    if os.path.exists(fits_file):
        fitld = AIPSTask('fitld', version=aipsver)
        fitld.datain   = fits_file
        fitld.outdata  = outdata
        fitld.douvcomp = 1
        fitld.digicor  = -1
        if outdata.exists():
            if vcs.yesno(f"Delete existing image file {outdata.name}? (No will abort pipeline)"):
                outdata.zap()
            else: sys.exit()
        fitld.inp()
        fitld()
    else: print(f"VLBI PIPE: fitld_image: {fits_file} file not found")
    
    
def splat_uvdata(indata, outdata, doband, bif, eif, solint):
    """
    SPLAT: Split uvdata based on IFs (BIF, EIF) without frequency averaging
    """
    if outdata.exists():
        print(f'VLBI PIPE: splat_uvdata: {outdata} is already in AIPSCat, I am going to abort the task')
        sys.exit()
    if indata.exists():
        splat = AIPSTask('splat', version=aipsver)
        splat.indata  = indata
        splat.outdata = outdata
        splat.bif     = bif
        splat.eif     = eif
        splat.solint  = solint
        splat.docalib = 1
        splat.doband  = doband
        splat.bpver   = 0
        splat.aparm[1:] = [0,0,0,0,1,0]
        splat.inp()
        splat()
    else: 
        print(f"VLBI PIPE: splat_uvdata: {indata} does not exist")
        sys.exit()


def split_uvdata(indata, sources, gainuse, outclass, outseq):
    """
    SPLIT: Split calibrated uvfits file to single source file by averaging channels
    (excise first and last two channels)
    Return: split message
    """
    outdata = AIPSUVData(sources[0], outclass, 1, outseq)
    if outdata.exists():
        if vcs.yesno(f"Delete existing UV dataset {outdata.name}? (No will abort pipeline)"):
            outdata.zap()
        else: sys.exit()
    split = AIPSTask('split', version = aipsver)
    split.indata   = indata
    split.outclass = outclass
    split.outseq   = outseq
    split.sources[1:] = sources
    try:
        nchan = indata.table('FQ', 1)[0].total_bandwidth[0]/(indata.table('FQ', 1)[0].ch_width[0])
    except (AttributeError, TypeError):
        nchan = indata.table('FQ', 1)[0].total_bandwidth/(indata.table('FQ', 1)[0].ch_width)
    if nchan > 8:
        split.bchan = 2
        split.echan = nchan-2
    split.docalib   = 1
    split.gainuse   = gainuse
    split.doband    = 1
    split.bpver     = 0
    split.aparm[1:] = [2,0,0,0,0,0]
    split.douvcomp  = -1
    split.inp()
    split.go()
    lines = split.message()
    return lines
    
def splittoseq(indata, sources, gainuse, outclass, outseq):
    """
    SPLIT: Split calibrated uvfits file to single source file by averaging channels
    set outseq and if outdata is already exit then it will be deleted
    Return: split message
    """
    splitdata = AIPSUVData(sources[0], outclass, 1,1)
    if splitdata.exists():
        splitdata.zap()
    split_uvdata(indata, sources, gainuse, outclass, outseq)


####### SOURCE POSITION CORRECTION #############################################
def shift_source(uvdataset, source, rashift, decshift, clversion):
    """
    Functionality
    -------------
    Shift the source position  by rashift and decshift. 
    It shifts RAEPO/DECEPO and RAAPP/DECAPP values of SU table.

    Input parameters
    ----------------
    rashift : float 
        in mas.
    decshift : float
        in mas.
    """
    clcor = AIPSTask('clcor', version = aipsver)
    clcor.indata = uvdataset
    clcor.sources[1] = source
    clcor.gainver = clversion
    clcor.gainuse = clversion
    clcor.opcode = 'ANTC'
    clcor.clcorprm[1:] = [0]
    clcor.clcorprm[5] = rashift/1000.0
    clcor.clcorprm[6] = decshift/1000.0
    clcor.clcorprm[8] = 0.0
    clcor.clcorprm[9] = 0.0
    clcor.infile = ""
    clcor.inp()
    clcor()



#====================== Data examination and plotting =========================
aips_tv = AIPSTV.AIPSTV()
def run_AIPSTV():
    if aips_tv.exists() == False:
        aips_tv.start()


def lwpla(indata, outputfile):
    """
    LWPLA: Save plot to disk
    """
    if os.path.exists(outputfile):
        os.system("rm -f " + outputfile)
    lwpla = AIPSTask('lwpla', version = aipsver)
    lwpla.indata   = indata
    lwpla.plver    = 1
    lwpla.invers   = 9999
    lwpla.outfile  = outputfile
    lwpla.lpen     = 5
    lwpla.dparm[5] = 1
    lwpla.dparm[5] = 0
    lwpla.docolor  = 1
    # lwpla.plcolors[10][1:] = [0,0,0]  # overall white background
    lwpla.plcolors[10][1:] = [1,1,1]
    lwpla()
    for table in indata.tables:
        if table[1] == 'AIPS PL':
            deletetable(indata, 'PL', table[0])


def uvcoverage(indata, outputfile, sources):
    """
    UVPLT: Plot uvcoverage
    """
    uvplt = AIPSTask('uvplt', version=aipsver)
    uvplt.indata = indata
    uvplt.sources[1:] = sources
    uvplt.docalib   = -1
    uvplt.bpver     = 0
    uvplt.dotv      = -1
    uvplt.do3color  = 1
    uvplt.bparm[1:] = [6,7,2]
    uvplt.go()
    lwpla(indata, outputfile)


def snplt(indata, outputfile, inext, invers, optype, bif=0, eif=0, timerang=[0]):
    """
    SNPLT: Plot AIPS tables (SN/CL/TY)
    """
    snplt = AIPSTask('snplt', version=aipsver)
    snplt.indata = indata
    # snplt.sources[1:] = sources
    snplt.inext  = inext
    snplt.invers = invers
    snplt.optype = optype
    snplt.opcode = 'ALSI'
    snplt.bif    = bif
    snplt.eif    = eif
    snplt.timerang[1:] = timerang
    snplt.dotv   = -1 
    snplt.do3color = 1
    # snplt.inp()
    snplt.go()
    lwpla(indata, outputfile)


def plotbp(indata, outputfile):
    """
    POSSM: Plot BP table
    """
    possm = AIPSTask('possm', version=aipsver)
    possm.indata    = indata
    possm.aparm[1:] = [0,1,0,0,-180,180,0,0,1]
    possm.aparm[8]  = 2
    possm.codetype  = 'A&P'
    possm.dotv      = -1
    possm.nplot     = 2
    possm.go()
    lwpla(indata, outputfile)


def possm(indata, outputfile, calibrate, gainuse, doband, sources=[''], flagver=0, nplot=1, 
          bif=0, eif=0, bchan=0, echan=0, solint=-1, timerang=[0]):
    """
    POSSM: Plot visibility amplitude and phase as a function of channels (frequency)
    """
    possm = AIPSTask('possm', version=aipsver)
    possm.indata = indata
    possm.sources[1:]  = sources
    possm.solint    = solint
    possm.bchan     = bchan
    possm.echan     = echan
    possm.flagver   = flagver
    possm.aparm[1:] = [0,1,0,0,-180,180,0,0,1]
    if calibrate == 1:
        possm.docalib = 1
        possm.gainuse = gainuse
        if doband == 1:
            possm.doband = 1
            possm.bpver  = 0
    possm.aparm[8]  = 0
    possm.codetype  = 'A&P'
    possm.timerang[1:] = timerang
    possm.bif    = bif
    possm.eif    = eif
    possm.dotv   = -1
    possm.nplot  = nplot
    possm.go()
    lwpla(indata, outputfile)


def uvplt_ampuvdist(indata, outputfile, doband, gainuse, sources, stokes='I', solint=0, nchavg=0):
    """
    UVPLT: Plot visibility amplitude as a function of uvdistance
    """
    uvplt = AIPSTask('uvplt', version=aipsver)
    uvplt.indata = indata
    uvplt.sources[1:] = sources
    uvplt.docalib   = 1
    uvplt.gainuse   = gainuse
    # uvplt.flagver   = flagver
    uvplt.doband    = doband
    uvplt.nchav     = nchavg
    uvplt.solint    = solint
    uvplt.stokes    = stokes
    uvplt.bpver     = 0
    uvplt.dotv      = -1
    uvplt.do3color  = 1
    uvplt.bparm[1:] = [0]
    uvplt.go()
    lwpla(indata, outputfile)


def uvplt_amptime(indata, outputfile, doband, gainuse, sources, solint=0, nchavg=0):
    """
    UVPLT: Plot visibility amplitude as a function of time
    """
    uvplt = AIPSTask('uvplt', version=aipsver)
    uvplt.indata = indata
    uvplt.sources[1:] = sources
    uvplt.docalib   = 1
    uvplt.gainuse   = gainuse
    # uvplt.flagver   = flagver
    uvplt.doband    = doband
    uvplt.nchav     = nchavg
    # uvplt.solint    = solint
    uvplt.bpver     = 0
    uvplt.dotv      = -1
    uvplt.do3color  = 1
    uvplt.bparm[1:] = [11,0]
    uvplt.go()
    lwpla(indata, outputfile)


def vplot_amptime(indata, outputfile, sources, gainuse, doband, solint, bif, eif):
    """
    VPLOT: Plot visibility amplitude v/s time for different baselines
    """
    vplot = AIPSTask('vplot', version=aipsver)
    vplot.indata = indata
    vplot.sources[1:] = sources
    vplot.bif = bif
    vplot.eif = eif 
    vplot.avgchan = 1
    vplot.docalib = 1
    vplot.gainuse = gainuse
    vplot.doband  = doband #Not working
    vplot.solint  = solint 
    vplot.bparm[1:] = [0,-1] 
    vplot.dotv   = -1
    vplot.crowded = 1
    vplot.do3color = 1
    vplot.nplot = 4
    vplot.go()
    lwpla(indata, outputfile)


def vbrfi(indata, outputfile, sources, solint, flagver, bif, eif):
    """
    VPLOT: Plot visibility amplitude v/s time for different baselines
    """
    if os.path.exists(outputfile):
        os.system("rm -f " + outputfile)
    vbrfi = AIPSTask('vbrfi', version=aipsver)
    vbrfi.indata = indata
    vbrfi.sources[1:] = sources
    vbrfi.bif = bif
    vbrfi.eif = eif 
    vbrfi.flagver = flagver
    vbrfi.solint  = solint  
    vbrfi.doplot = 1
    vbrfi.dotv   = -1
    vbrfi.outfile = outputfile.replace('.ps', '.txt')
    vbrfi.go()
    lwpla = AIPSTask('lwpla', version = aipsver)
    lwpla.indata = indata
    lwpla.plver = 1
    lwpla.invers = 9999
    lwpla.outfile = outputfile
    lwpla.lpen = 5
    lwpla.dparm[5] = 0
    lwpla.docolor = 1
    lwpla.plcolors[1][1:] = [0,0,0]
    lwpla.plcolors[2][1:] = [0.06275,1,0]
    lwpla.plcolors[3][1:] = [1,0.06275,1]
    lwpla.plcolors[4][1:] = [0,1,1]
    lwpla.plcolors[9][1:] = [0,0,0] #[1,1,0]
    lwpla.plcolors[10][1:] = [1,1,1]
    lwpla.inp()
    lwpla()
    for table in indata.tables:
        if table[1] == 'AIPS PL':
            deletetable(indata, 'PL', table[0])
    

def plotan(indata, outputfile):
    """
    PRTAN: Plot antennas location on UV-plane
    """
    prtan = AIPSTask('prtan', version = aipsver)
    prtan.indata = indata
    prtan.inver  = 1
    prtan.doplot  = 1
    prtan.dotv = -1
    prtan.outprint = outputfile
    prtan.docrt = -1
    # prtan.inp()
    prtan()
    lwpla(indata, outputfile)


def fringerateflag(indata, suppressionfactor, flagver=1):
    """
    UVFLG: Fringe rate based flagging
    Return: uvflg message
    """
    uvflg = AIPSTask('uvflg', version = aipsver)
    uvflg.indata    = indata
    uvflg.intext    = ''
    uvflg.outfgver  = flagver
    uvflg.aparm[1:] = [0]
    uvflg.aparm[3]  = suppressionfactor
    uvflg.opcode    = 'FLAG'
    uvflg.reason    = 'FRINGERATE'
    uvflg.inp()
    uvflg()
    lines = uvflg.message()
    return lines


##### User-specified additional flags ##########################################
def userflag(uvdataset, outflagver, filename):
    uvflg = AIPSTask('uvflg', version = aipsver)
    uvflg.indata = uvdataset
    uvflg.intext = filename
    uvflg.outfgver = outflagver
    # uvflg.opcode = 'FLAG'      # Ashish added it, already added in flag file
    # uvflg.reason = 'USERFLAG'  # Ashish added it
    uvflg.inp()
    uvflg()


def listr(indata, outprint, optype):
    """
    LISTR: Observation scan information
    """
    listr = AIPSTask('listr', version = aipsver)
    listr.indata   = indata  
    listr.outprint = outprint
    listr.optype   = optype
    listr.docalib  = 1
    listr()


def antflag(indata, antenna, flagver=1, sources=''):
    """
    UVFLG: Flag antennas
    Return: uvflg message
    """
    uvflg = AIPSTask('uvflg', version = aipsver)
    uvflg.indata      = indata
    uvflg.outfgver    = flagver
    uvflg.antenna[1:] = antenna
    uvflg.aparm[1:]   = [0]
    uvflg.reason      = 'TSYS'
    uvflg.sources[1]  = sources
    uvflg.opcode      = 'FLAG'
    uvflg.inp()
    uvflg()
    lines = uvflg.message()
    return lines
    
    
def snedt(indata, snver):
    """
    SNEDT: Edit SN table
    """
    snedt = AIPSTask('snedt', version = aipsver)
    snedt.indata  = indata
    snedt.inext   = 'sn'
    snedt.invers  = snver
    snedt.dodelay = 1
    snedt.freqid  = -1
    snedt.dotwo   = 1
    run_AIPSTV()
    snedt()


##### Find a timerange when the main calibrator is up ##########################
def get_ampcal_timer(listrfile, srcname, scanno):
    timerang = []
    listrjunk = open(listrfile)
    listrmsgs = listrjunk.readlines()
    listrjunk.close()
    
    atline = 0
    while len(listrmsgs[atline].split()) < 5 or \
              not listrmsgs[atline].split()[0] == "Scan":
        atline = atline + 1
    for i in range(scanno):
        atline = atline + 1
        while atline < len(listrmsgs):
            if len(listrmsgs[atline].split()) <= 1:
                atline += 1
            else:
                if listrmsgs[atline].split()[1] == srcname:
                     print("Found " + srcname + " for scan " + \
                           str(i+1) + "/" + str(scanno))
                     break
                atline += 1
    if atline >= len(listrmsgs):
        print("Couldn't find scan " + str(scanno) + " on source " + srcname)
        print("Sorry, but I'll have to abort the pipeline")
        sys.exit()
    print(listrmsgs[atline])
    starttime = listrmsgs[atline][41:51]
    endtime = listrmsgs[atline][56:66]
    print("Starttime of ampcal scan is " + starttime)
    print("Endtime of ampcal scan is " + endtime)

    startday  = int(starttime[0])
    starthour = int(starttime[2:4])
    startmin  = int(starttime[5:7])
    startsec  = int(starttime[8:10])
    endday    = int(endtime[0])
    endhour   = int(endtime[2:4])
    endmin    = int(endtime[5:7])
    endsec    = int(endtime[8:10])
    #Clip some, if we think we can get away with it
    timediff  = (endday-startday)*86400 + (endhour-starthour)*3600 + (endmin-startmin)*60 + \
                            endsec - startsec
    if timediff >= 70:
        endsec   -= 10
        timediff -= 10
        if endsec < 0:
            endsec += 60
            endmin -= 1
            if endmin < 0:
                endmin += 60
                endhour -= 1
                if endhour < 0:
                    endhour += 24
                    endday -= 1
        tosubtract = 0
        if timediff > 600:
            tosubtract = 240
        elif timediff > 360:
            tosubtract = 120
        elif timediff > 180:
            tosubtract = 60
        elif timediff > 120:
            tosubtract = 30
        elif timediff > 90:
            tosubtract = 20
        elif timediff > 70:
            tosubtract = 10
        startday  = endday
        starthour = endhour
        startmin  = endmin
        startsec  = endsec - (timediff - tosubtract)
        while startsec < 0:
            startmin -= 1
            startsec += 60
        while startmin < 0:
            startmin += 60
            starthour -= 1
        if starthour < 0:
            starthour += 24
            startday -= 1
    timerang.append(startday)
    timerang.append(starthour)
    timerang.append(startmin)
    timerang.append(startsec)
    timerang.append(endday)
    timerang.append(endhour)
    timerang.append(endmin)
    timerang.append(endsec)
    return timerang


##### Do the pulse cal corrections #############################################
def pccor(uvdataset, listrfile, srcname, snversion, ampcalscan, refant):
    pccor = AIPSTask('pccor', version = aipsver)
    pccor.indata = uvdataset
    pccor.timerang[1:] = get_ampcal_timer(listrfile, srcname, ampcalscan)
    print ("Running PCCOR over timerange:")
    print (pccor.timerang[1:])
    pccor.snver = snversion
    pccor.inver = 0
    pccor.refant = refant
    pccor.cutoff = 0
    pccor.delcorr = 0
    pccor.inp()
    pccor()


def tecor_iono(indata, ionodir, vexfile, ionextype, petrov=False):
    """
    TECOR: Ionospheric correction using TEC model
    ionextypes = ["jplg","esag","codg","upcg","igsg","c1pg","c2pg","u2pg","e1pg"]
    """
    sdoy, syear, edoy, eyear = vcs.get_obs_start_end_doy_year(vexfile)
    numfile = edoy - sdoy

    infile = ionodir + '/' + ionextype + '%3.3d0.%si' % (sdoy,str(syear)[2:4]) 
    tecor  = AIPSTask('tecor', version=aipsver)
    tecor.indata  = indata
    tecor.infile  = infile
    tecor.nfiles  = int(numfile) + 1
    tecor.gainver = 0
    tecor.gainuse = 0
    if petrov:
        tecor.aparm[1:] = [1, 0.2, 1, 0.85, 56.7, 0.9782]
    else:
        tecor.aparm[1:] = [1,0.2,0,0,0,0]
    tecor.inp()
    tecor.go()


def eops_clcor(indata, eopspath):
    """
    CLCOR: EOPS correction
    """
    clcor = AIPSTask('clcor', version=aipsver)
    clcor.indata = indata     
    clcor.opcode = 'EOPS'
    clcor.clcorprm[1] = 1
    clcor.clcorprm[2] = 5
    clcor.infile = eopspath
    clcor.inp()
    clcor()


def pang_clcor(indata):
    """
    CLCOR: Parallactic angle correction
    """
    clcor = AIPSTask('clcor', version=aipsver)
    clcor.indata = indata
    clcor.opcode = 'PANG'
    clcor.clcorprm[1] = 1
    clcor.inp()
    clcor()


def accor(indata, solint):
    """
    ACCOR: Auto-correlation correction
    """
    accor=AIPSTask('accor', version=aipsver)
    accor.indata = indata
    accor.solint = solint
    accor.inp()
    accor()
    

def apcal(indata):
    """
    APCAL: Amplitude calibration
    """
    apcal = AIPSTask('apcal', version=aipsver)
    apcal.indata    = indata
    apcal.aparm[1]  = 1
    apcal.dofit[1:] = [-1]
    apcal.inp()
    apcal()


def snsmo(indata, refant, invers, smotype, timemins, ampdev, phasedev, delaydev, passthrough, 
          ratesmoothmin, ratedev=0):
    """
    SNSMO: Smoothing and clipping of the SN table
    """
    snsmo = AIPSTask('snsmo', version=aipsver)
    snsmo.indata    = indata
    snsmo.samptype  = 'MWF'
    snsmo.smotype   = smotype
    snsmo.refant    = refant
    snsmo.invers    = invers
    snsmo.bparm[1:] = [0]
    snsmo.doblank   = -1
    snsmo.dobtween  = 1 # by default 1
    if passthrough==1:
        snsmo.doblank = 1
        snsmo.dobtween = 1
        snsmo.bparm[1] = timemins/60.0
        snsmo.bparm[2] = timemins/60.0
        snsmo.bparm[3] = timemins/60.0
        snsmo.bparm[4] = timemins/60.0
        snsmo.bparm[5] = timemins/60.0

    snsmo.cparm[1:] = [0]
    if ampdev > 0.0:
        snsmo.cparm[1] = timemins/60
        snsmo.cparm[6] = ampdev
    if phasedev > 0.0:
        snsmo.cparm[2] = timemins/60
        snsmo.cparm[7] = phasedev
    if delaydev > 0.0:
        snsmo.cparm[4] = timemins/60
        snsmo.cparm[9] = delaydev
    if ratedev > 0.0:
        snsmo.cparm[3] = timemins/60
        snsmo.cparm[8] = ratedev

    if ratesmoothmin > 0:
        snsmo.samptype = 'BOX'
        snsmo.smotype  = 'VLRI'
        snsmo.bparm[3] = ratesmoothmin/60.0
    snsmo.inp()
    snsmo()
    
    
def fring_fitting(indata, refant, in2data, calsour, doband, zerorate, solint, dostokesi, snr, 
                  delaywin=400,ratewin=30,inttimesecs=2):
    """
    FRING: Fringe-fitting
    """
    fring = AIPSTask('fring', version=aipsver)
    fring.indata      = indata  
    fring.calsour[1:] = calsour
    fring.doband      = doband
    fring.bpver       = 0
    fring.docalib     = 1
    fring.refant      = refant
    fring.solint      = solint 
    fring.smooth[1]   = 5         # Hanning smoothing, width 4 channels
    fring.in2data     = in2data
    try:
        nchan = indata.table('FQ', 1)[0].total_bandwidth[0]/indata.table('FQ', 1)[0].ch_width[0]
    except (AttributeError, TypeError):
        nchan = indata.table('FQ', 1)[0].total_bandwidth/indata.table('FQ', 1)[0].ch_width
    if nchan > 8:
        fring.bchan = 2
        fring.echan = nchan - 1
    else:
        fring.bchan = 1
        fring.echan = nchan

    fring.aparm[1:] = [3,0,0,0,0,2,snr,0,1,0]
    if dostokesi > 0:
        fring.aparm[3] = 1
 
    fring.dparm[1:] = [0,delaywin,ratewin,inttimesecs,0,0,0,0]
    fring.dparm[8] = zerorate 
    fring.inp() 
    fring()
    lines = fring.message()
    return lines


def bpass(indata, refant, in2data, calsour):
    """
    BPASS: Bandpass Calibration
    """
    bpass = AIPSTask('bpass', version=aipsver)
    bpass.indata  = indata
    bpass.docalib = 1
    bpass.gainuse = 0
    bpass.solint  = 0
    bpass.refant  = refant
    bpass.soltype = 'L1R'
    bpass.calsour[1:] = calsour
    bpass.doband  = -1
    bpass.bpver   = 0
    bpass.bpassprm[5]  = 1
    bpass.bpassprm[9]  = 1
    bpass.bpassprm[10] = 1
    bpass.in2data = in2data
    bpass.inp()
    bpass.go()


def clcal_applySN(indata, refant, sources, snver, interpol='2PT', opcode='CALI'):
    """
    CLCAL: Apply SN table
    """
    clcal = AIPSTask('clcal', version=aipsver)
    clcal.indata      = indata 
    clcal.sources[1:] = sources
    clcal.interpol    = interpol
    clcal.refant      = refant
    clcal.snver       = snver
    clcal.dobtween    = -1
    clcal.doblank     = 1
    clcal.opcode      = opcode
    clcal.inp()
    clcal()

 
def uvsub_divide (indata, in2data):
    """
    UVSUB: Divide data by model
    """
    uvsubdata = AIPSUVData(indata.name, "uvsub", 1, 1)
    if uvsubdata.exists():
        print(f"VLBI PIPE: uvsub_divide: {uvsubdata} already exists")
        if vcs.yesno(f"Delete existing UV dataset {indata.name}? (No will abort pipeline)"):
            uvsubdata.zap()
        else: sys.exit()
    uvsub = AIPSTask('uvsub', version=aipsver)
    uvsub.indata  = indata
    uvsub.in2data = in2data
    if not in2data.exists():
        print(f"VLBI PIPE: uvsub_divide: {in2data} model data is not present, I am going to abort the task")
        sys.exit()
    uvsub.outdata = uvsubdata
    uvsub.nmaps   = 1
    uvsub.channel = 0
    uvsub.flux    = 0
    uvsub.cmethod = 'DFT'
    uvsub.cmodel  = 'COMP'
    uvsub.opcode  = 'DIV'
    uvsub.inp()
    uvsub()
    return uvsubdata


def fittp (indata, dataout):
    """
    FITTP: Store fits file to disk
    """
    if os.path.exists(dataout):
        print(f"VLBI PIPE: fittp: {dataout.split('/')[-1]} file already exist")
        if vcs.yesno(f"Remove disk file {dataout.split('/')[-1]} to proceed (no will abort pipeline)? "):
            os.remove(dataout)
            print(f"VLBI PIPE: fittp: {dataout.split('/')[-1]} removed")
        else:
            sys.exit()
    fittp = AIPSTask('fittp')
    fittp.indata  = indata
    fittp.dataout = dataout
    fittp.inp()
    fittp.go()


####### GIVES THE GREAT CIRCLE ANGULAR DISTANCE BETWEEN TWO POINTS #############
def haversine(az1, el1, az2, el2):
    """
    Determine the great circle angular distance between two points in spherical polar coordinates
    This should give the distance from the pointing centre (az1, el1) of a 2nd point (az2, el2)
    See the Ipad (Ashish) notes for more details of this formula
    """
    a = math.sin((el1 - el2)/2)**2 + (math.sin((az1 - az2)/2)**2 * math.cos(el1)*math.cos(el2))
    return 2*math.atan2(math.sqrt(a), math.sqrt(1 - a))


####### GIVES THE INITIAL BEARING FROM NORTH OF POINT 1 FROM POINT 2 ###########
# Bearing is the angle measured clockwise from North.
def azbearing(az1, el1, az2, el2):
    """
    Initial bearing from North of point 1 from point 2
    For an az/el mounted telescope, North is always 'up' so this should give the bearing within 
    the primary beam of (az2, el2) from a pointing centre (az1, el1).
    The exact form is such that if the two points share the same RA, and the decl(2) > decl(1),
    the angle returned is the Parallactic Angle
    NB not tested for southern hemisphere antennas.
    """
    return -math.atan2(math.sin(az2 - az1)*math.cos(el1), 
                    (math.cos(el2)*math.sin(el1)) - (math.sin(el2)*math.cos(el1)*math.cos(az2 - az1)))
        
        
####### GIVES THE SEPARATION BETWEEN TWO POINTINGS IN AZ AND EL ################
def get_beam_position(pointing_az, pointing_al, target_az, target_al):
    """
    return cartesian offset (N and E) of pointing from target
    inputs in radians
    return position of target w.r.t. pointing in cartesian coordinates
    az, al both in radians
    """
    # work out distance from pointing to target (doesn't really need recalculating every time but this saves 
    # passing it to this function)
    r = haversine(target_az, target_al, pointing_az, pointing_al)
    # Work out angle from North in the primary beam
    theta = azbearing(target_az, target_al, pointing_az, pointing_al)
    return r*math.sin(theta), r*math.cos(theta)


####### THEORETICAL RESPONSE OF A FILLED DISH  #################################
def airyresponse(theta, D, lam):
    u = math.sin(theta)
    if u==0:
        return 1.0
    else:
        return (2*jn(1,(math.pi*u*D/lam))/(math.pi*u*D/lam))**2
        
        
####### GIVES THE ANTENNA RESPONSE FOR A GIVEN POINTING ########################
def ant_response(now, antcode, cephem, tephem, wavelength, verbose):
    """
    Calculates the beam response of an antenna given its
    pointing centre, the target coordinates, bmaj in Az and El, and the
    amount and direction of beam squint
    """
    # Get the antenna and beam location, and the wavelength and scale compared to reference wavelength
    ant  = vlbalocations[antcode]
    beam = vlbabeams[antcode]
    f = 299792458/wavelength
    beamcorr = beam[-1]/f

    # get date and location of antenna and calculate the separation between pointing centre and target
    date = ephem.Date(now)
    obs  = ephem.Observer()
    obs.lat, obs.long, obs.elevation = ant[0], ant[1], int(ant[2])
    obs.date = date

    # pyephem calculates Az/El using an atmospheric model, set pressure to zero to switch this off
    # obs.pressure=0.0
    # this is the "international standard atmosphere"
    # see http://de.wikipedia.org/wiki/Barometrische_H%C3%B6henformel#Internationale_H.C3.B6henformel
    
    obs.pressure = 1013.25*(1 - 0.0065*obs.elevation/288.15)**5.255
    cephem.compute(obs)
    tephem.compute(obs)

    # calculate how much the beams differ in rad
    # FIXME: This is inverted, but seems to work??
    # division factor is not clear to me
    DAz = -beam[0]*beamcorr*(math.pi/180)/2/60
    DEl = -beam[1]*beamcorr*(math.pi/180)/2/60

    # Get the separation between pointing center and target
    sep = get_beam_position(float(cephem.az), float(cephem.alt), float(tephem.az), float(tephem.alt))

    # Add squint offsets, assuming everything is in a plane as angles are small
    Az_LCP = sep[0] + DAz
    Az_RCP = sep[0] - DAz
    El_LCP = sep[1] + DEl
    El_RCP = sep[1] - DEl

    # calculate the total offsets, in radians
    sep_LCP = math.sqrt(Az_LCP**2 + El_LCP**2)
    sep_RCP = math.sqrt(Az_RCP**2 + El_RCP**2)

    # get the primary beam responses as a function of that separation
    if antcode == "GB":
        LCP = airyresponse(sep_LCP, gbtdiameter, wavelength)
        RCP = airyresponse(sep_RCP, gbtdiameter, wavelength)
    else:
        LCP = airyresponse(sep_LCP, vlbadiameter, wavelength)
        RCP = airyresponse(sep_RCP, vlbadiameter, wavelength)

    # make lots of noise if desired
    if verbose:
        print("\nTime              : %s" % date)
        print("Station             : %s Long: %s Lat: %s Elev: %im Press: %.1f" % (antcode, obs.long, obs.lat, obs.elevation, obs.pressure))
        print("Pointing centre     : ra=%s  dec=%s  az=%s  el=%s" % (cephem.ra, cephem.dec, cephem.az, cephem.alt))
        print("LCP pointing centre : az=%s  el=%s" % (Az_LCP, El_LCP))
        print("RCP pointing centre : az=%s  el=%s" % (Az_RCP, El_RCP))
        print("Target              : ra=%s  dec=%s  az=%s  el=%s" % (tephem.ra, tephem.dec, tephem.az, tephem.alt))
        print("Scaling factor      : %.4f" % beamcorr)
        print("Beam offsets        : %.4f'  %.4f'" % (beam[0], beam[1]))
        print("Scaled beam offsets : %.4f'  %.4f'" % (beam[0]*beamcorr, beam[1]*beamcorr))
        print("Beam sizes          : %.4f'  %.4f'  %.4f'  %.4f'" % (beam[2], beam[3], beam[4], beam[5]))
        print("Scaled beam sizes   : %.4f'  %.4f'  %.4f'  %.4f'" % (beam[2]*beamcorr, beam[3]*beamcorr, beam[4]*beamcorr, beam[5]*beamcorr))
        print("LCP-Target          : %.4f'" % (sep_LCP*180*60/math.pi))
        print("RCP-Target          : %.4f'" % (sep_RCP*180*60/math.pi))
        print("LCP response        : %.5f" % LCP)
        print("RCP response        : %.5f" % RCP)
    return (RCP, LCP)


def correct_primarybeam(uvdata, examplesnversion, vexfile, fieldsourcenames, iscal, phasecentrenum, 
                issearch=True, isonepointing=False, onlygettimes=False, skipmissingsources=False):
    """
    Make an amplitude correction for primary beam effects
    """
    scanlist       = vcs.getvexscaninfo(vexfile) # Ashish added this line
    numantennas    = len(uvdata.antennas)
    wizuvdata      = WizAIPSUVData(uvdata)
    examplesntable = wizuvdata.table('SN', examplesnversion)
    num_if         = examplesntable.keywords['NO_IF']
    num_pol        = examplesntable.keywords['NO_POL']
    num_snterms    = num_if*num_pol

    antable   = uvdata.table('AN', 1)
    annumbers = []
    annames   = []
    maxannumber = 0
    for row in antable:
        annumbers.append(row.nosta)
        annames.append(row.anname.strip())
        if row.nosta > maxannumber:
            maxannumber = row.nosta

    freqs   = []
    lambdas = []
    equivtimeonsource = {}
    fqtable = uvdata.table('FQ', 1)
    for row in fqtable:
        try:
            for iffreq, ifbw, ifsideband in zip(row.if_freq, row.total_bandwidth, row.sideband):
                freqentry = float(iffreq) + float(uvdata.header.crval[2]) + float(ifbw)*float(ifsideband)/2.0
                freqs.append(freqentry)
        except (AttributeError, TypeError):
            freqentry = float(row.if_freq) + float(uvdata.header.crval[2]) + float(row.total_bandwidth)*float(row.sideband)/2.0
            freqs.append(freqentry)
    for f in freqs:
        lambdas.append(299792458.0/f)

    if not onlygettimes:
        newsntable = wizuvdata.attach_table('SN', examplesnversion+1, no_term=num_snterms)
        newsntable.keywords['NO_IF']   = num_if
        newsntable.keywords['NO_POL']  = num_pol
        newsntable.keywords['NO_TERM'] = num_snterms
        newsntable.keywords['NO_ANT']  = maxannumber
    row = examplesntable[0]
    gainlength = num_if
    newdelay1 = []
    newdelay2 = []
    newrate1  = []
    newrate2  = []
    newimag1  = []
    newimag2  = []
    for i in range(gainlength):
        newdelay1.append(0.0)
        newdelay2.append(0.0)
        newrate1.append(0.0)
        newrate2.append(0.0)
        newimag1.append(0.0)
        newimag2.append(0.0)
    row.delay_1 = newdelay1
    row.delay_2 = newdelay2
    row.rate_1 = newrate1
    row.rate_2 = newrate2
    row.imag1 = newimag1
    row.imag2 = newimag2
    row.i_far_rot = 0
    row.time_interval = 0.0
    row.source_id = 1

    aipsstartmjd = vcs.getAIPSMJD(uvdata.header.date_obs)
    sutable = uvdata.table('SU', 1)
    numsusources = 0
    for srow in sutable:
        numsusources += 1
    for scan in scanlist:
        meancorr = 0.0
        meancorrcount = 0
        sutablerow = None
        scansource = scan.source
        print("Looking for match to vexname ", scansource.name)
        
        if iscal:
            localname = scan.source.name
            if numsusources == 1: #Probably a bd161 experiment, must be this source
                for srow in sutable:
                    sutablerow = srow
                    break
            if (localname[-2] == '-' or "PT" in localname[-4:]) and numsusources > 1:
                continue #Target source, not important
            for srow in sutable:
                if srow.source.strip() == localname:
                    sutablerow = srow
                    break
            if sutablerow == None:
                print("Couldn't find source " + localname + " in SU table!")
                if skipmissingsources:
                    continue
                else:
                    sys.exit()
        else:
            if issearch:
                if not scansource.name[-2] == '-' and not "-PT" in scansource.name[-4:-1] and not isonepointing: #A calibrator
                    print ("Skipping source " + scansource.name)
                    continue
                ssname = scansource.name[:-2]
                if "-PT" in scansource.name[-4:-1]:
                    ssname = scansource.name[:-4]
                tcount = 0
                # Figure out what source this is, calculate PB attenuation
                if len(fieldsourcenames) > 1:
                    for s in fieldsourcenames:
                        if ssname == s[0][0] or (len(ssname) == 5 and ssname in s[0][0]):
                            break
                        elif ssname[-1] == "A" and ssname[:-1] in s[0][0] and s[0][0][-1] == "A":
                            break
                        tcount += 1
                    if tcount == len(fieldsourcenames):
                        print ("Couldn't find source " + scansource.name + " (was looking for " + ssname + ")")
                        if skipmissingsources:
                            continue
                        else:
                            sys.exit()
                else:
                    print ("There is only one target, ssname is " + ssname)
                fcount = int(scansource.name[-1]) - 1
                print (fcount)
                
                if phasecentrenum >= len(fieldsourcenames[tcount][fcount]):
                    print ("No phase centre %d for pointing %s" % (phasecentrenum, scansource.name))
                    continue
                localname = fieldsourcenames[tcount][fcount][phasecentrenum]
            else:
                try:
                    localname = fieldsourcenames[scan.source.name]
                except KeyError:
                    print ("No inbeam for source " + scan.source.name + " for this pass")
                    continue
            for srow in sutable:
                print (srow.source.strip())
                if srow.source.strip() == localname:
                    sutablerow = srow
                    break
            if sutablerow == None:
                print ("Couldn't find source " + localname + " in SU table!")
                if issearch:
                    sys.exit()
                else:
                    continue
                
        # create centre and target objects in pyephem
        crastr, cdecstr = ast.posradians2string(scansource.ra, scansource.dec)
        cephem = ephem.readdb("%s,f|M|F7,%s, %s,2.02,2000" % (scansource.name, crastr, cdecstr))
        trastr, tdecstr = ast.posradians2string(math.pi*sutablerow.raepo/180.0, math.pi*sutablerow.decepo/180.0)
        tephem = ephem.readdb("%s,f|M|F7,%s, %s,2.02,2000" % (scansource.name, trastr, tdecstr))

        #Do three points - 5 sec after start of scan, one in middle of scan, one 5 sec from end of scan
        scanpoints = []
        scanpoints.append(scan.startmjd - aipsstartmjd + 5.0/86400.0)
        scanpoints.append((scan.startmjd+scan.stopmjd)/2.0 - aipsstartmjd)
        scanpoints.append(scan.stopmjd - aipsstartmjd - 5.0/86400.0)
        
        for scanpoint in scanpoints:
            year, month, day, hour, minute, second = ast.mjd2ymdhms(scanpoint + aipsstartmjd)
            now = datetime.datetime(year, month, day, hour, minute, int(second))
            row.time = scanpoint
            row.source_id = sutablerow.id__no
            # Put in an entry for every antenna
            for i in range(numantennas):
                row.antenna_no = annumbers[i]
                newreal1 = []
                newreal2 = []
                newweight1 = []
                newweight2 = []
                for j in range(gainlength):
                    rcpcorr, lcpcorr = ant_response(now, annames[i], cephem, tephem, lambdas[j], False)
                    meancorr += rcpcorr + lcpcorr
                    meancorrcount += 2
                    newreal1.append(1.0/math.sqrt(rcpcorr))
                    newreal2.append(1.0/math.sqrt(lcpcorr))
                    newweight1.append(rcpcorr)
                    newweight2.append(lcpcorr)
                row.real1 = newreal1
                row.real2 = newreal2
                row.weight_1 = newweight1
                row.weight_2 = newweight2
                print ("Doing time %f for antenna %d" % (scanpoint, i+1))
                if not onlygettimes:
                    newsntable.append(row)
        try:
            test = equivtimeonsource[localname]
        except KeyError:
            equivtimeonsource[localname] = 0.0
        scanlen = (scan.stopmjd - scan.startmjd)*86400.0
        meancorr /= meancorrcount
        equivtimeonsource[localname] += scanlen*meancorr*meancorr
        print ("Adding %f seconds (weighted to %f) to %s onsourcetime, it is now %f, meancorr was %f" % \
              (scanlen, scanlen*meancorr*meancorr, localname, equivtimeonsource[localname], meancorr))
    if not onlygettimes:
        newsntable.close()
    print("VLBI PIPE: correct_primarybeam: Writing SN table %d for subarray  1" %(examplesnversion+1))
    print("VLBI PIPE: correct_primarybeam: Appears to have ended successfully")
    return equivtimeonsource


def calib_selfcal(indata, in2data, refant, calsour, solint, dostokesi, doamp, avgif, calibrate, 
                  snr, weightit=0, normalize=0):
    """
    CALIB: Self-Calibration (Amplitude and Phase)
    if you want to use a unit Jy point source then set in2data=None, and 
    if you want to vary flux 
    Return: CALIB message
    """
    calib = AIPSTask('calib', version=aipsver)
    calib.indata = indata
    calib.calsour[1:] = calsour    
    if in2data == None:
        calib.smodel[1] = 1
        calib.smodel[2:] = [0]
    elif isinstance(in2data, float) or isinstance(in2data, int):
        print("Calibrating with a delta function of amp " + str(in2data))
        calib.smodel[1] = in2data
        calib.smodel[2:] = [0]
    calib.cmethod   = 'DFT'
    calib.aparm[1:] = [3,0,0,0,0,1,snr,0,0,0]
    calib.solint    = solint
    calib.refant    = refant
    calib.soltype   = 'L1R'
    calib.in2data   = in2data
    calib.weightit  = weightit
    calib.doapply   = -1
    # if data set has no CL and BP 
    if calibrate > 0:
        calib.docalib = 1
        calib.doband  = 1
        calib.bpver   = 0
    else:
        calib.docalib = 0
        calib.doband  = 0
        calib.bpver   = 0
    if dostokesi > 0:
        calib.aparm[3] = 1
    if avgif > 0:
        calib.aparm[5] = 1
    if doamp > 0:
        calib.solmode = 'A&P'
        if normalize > 0:
            calib.normaliz = 1
    else:
        calib.solmode = 'P'
        if normalize:
            print("VLBI PIPE: calib_selfcal: You can't normalise phase-only solutions")
    calib.inp()
    calib()
    lines = calib.message()
    return lines


def tacop(indata, outdata, inext, inver):
    """
    TACOP: Copy an AIPS table from an AIPSCat to another
    """
    tacop = AIPSTask('tacop')
    tacop.indata  = indata
    tacop.outdata = outdata
    tacop.inext   = inext
    tacop.inver   = inver
    tacop.ncount  = 1
    tacop.inp()
    tacop()

 
def morif (indata, outname, npiece):
    """
    MORIF: Split IFs into more IFs
    """
    outdata = AIPSUVData(outname, 'MORIF4',1,1)
    if outdata.exists():
        if vcs.yesno(f"Delete existing UV dataset {outdata} ? (No will abort pipeline)"):
            outdata.zap()
        else: sys.exit()
    morif = AIPSTask('morif')
    morif.indata  = indata
    morif.npiece  = npiece
    morif.outdata = outdata
    morif.inp()
    user_ip = str(input("Do you want to split IF (MORIF [y/n])?:"))
    if user_ip.upper()=='Y':
        morif()
        
        
def writetable(uvdataset, tabletype, tableversion, filename):
    """
    TBOUT: Write an AIPS table to disk
    """
    tbout = AIPSTask('tbout', version = aipsver)
    tbout.indata  = uvdataset
    tbout.inext   = tabletype
    tbout.invers  = tableversion
    if os.path.exists(filename):
        if vcs.yesno(f"Delete existing disk SN table '{filename.split('/')[-1]}' ? (No will abort pipeline)"):
            os.remove(filename)
        else: sys.exit()
    tbout.outtext = filename
    # tbout.inp()
    tbout()


def loadtable(uvdataset, tablename, tableversion):
    """
    Load a table to AIPS
    """
    tables = uvdataset.tables
    tabletype = tablename.split('.')[-1].upper()
    for table in tables:
        if table[1][-2:] == tabletype and table[0] >= tableversion:
            if vcs.yesno(f"Can i delete {table[1][-2:]} table version {str(table[0])}? (No will abort pipeline)"):
                uvdataset.table(tabletype, table[0]).zap()
            else:
                sys.exit()
    tbin = AIPSTask('tbin', version = aipsver)
    tbin.outdata = uvdataset
    tbin.intext = tablename
    tbin()


def getimagestats(image, cutofffrac):
    """
    Use imean to get stats from an image
    cutofffrac = 1 (full image), 0.7 (70 % of image)
    Return: [peak, rms, peakx, peaky, rastr, decstr]
    """
    imean = AIPSTask('imean', version = aipsver)
    imean.indata = image
    xguardpix = int(image.header.naxis[0]*cutofffrac)
    yguardpix = int(image.header.naxis[1]*cutofffrac)
    print("Avoiding a guard band of %dx%d pixels around edge of image..." % \
          (xguardpix, yguardpix))
    imean.blc[1:]  = [xguardpix, yguardpix]
    imean.trc[1:]  = [image.header.naxis[0]-xguardpix,image.header.naxis[1]-yguardpix]
    imean.dohist   = -1
    imean.doinvers = -1
    imean.outtext  = ''
    imean.nboxes   = 0

    try:
        imean()
    except RuntimeError:
        print("IMEAN failed - this image must be bogus!")
        return [0.0, 1.0, image.header.naxis[0]/2, image.header.naxis[1]/2, "", ""]
    
    imeanlines = imean.message()
    numlines = len(imeanlines)
    atline = 0
    while atline < numlines and imeanlines[atline].split()[-1] != "histogram":
        atline += 1
    if atline >= numlines:
        atline = 0
        noise1 = -1
    else:
        val = imeanlines[atline].split()[4]
        if val == "****":
            val = imeanlines[atline].split()[3]
        noise1 = float(val)
    while atline < numlines and imeanlines[atline].split()[-1] != "pixels":
        atline += 1
    val = imeanlines[atline].split()[4]
    if val == "JY/BEAM":
        val = imeanlines[atline].split()[3]
    noise2 = float(val)
    while atline < numlines and imeanlines[atline].split()[1] != "Skypos:":
        atline += 1
    while atline < numlines and not "Maximum" in imeanlines[atline]:
        atline += 1
    val1 = imeanlines[atline].split()[2]
    val2 = imeanlines[atline].split()[4]
    val3 = imeanlines[atline].split()[5]
    if val1 == "at":
        val1 = (imeanlines[atline].split()[1]).split('=')[-1]
        val2 = imeanlines[atline].split()[3]
        val3 = imeanlines[atline].split()[4]
    peak  = float(val1)
    peakx = int(val2)
    peaky = int(val3)
    while atline < numlines and not "Skypos" in imeanlines[atline]:
        atline += 1
    splitline = imeanlines[atline].split()
    rastr = splitline[3] + ":" + splitline[4] + ":" + splitline[5]
    decstr = splitline[7] + ":" + splitline[8] + ":" + splitline[9]
    if noise1 < 0:
        return [peak, noise2, peakx, peaky, rastr, decstr]
    return [peak, noise1, peakx, peaky, rastr, decstr]


def wizCorrectScint(parentuvdata, examplesnversion, outputsnversion, splituvdata, solmins, outputsntable):
    """
    Pulsar scintillation correction
    """
    # Get an example SN table, prefill all the entries that can be prefilled
    print("VLBI PIPE: wizCorrectScint: Running scintillation correction task")
    wizuvdata      = WizAIPSUVData(parentuvdata)
    examplesntable = wizuvdata.table('SN', examplesnversion)
    num_if         = examplesntable.keywords['NO_IF']
    num_pol        = examplesntable.keywords['NO_POL']
    num_snterms    = num_if*num_pol
    numantennas    = len(parentuvdata.antennas)
    examplerow     = examplesntable[0]
    annumbers      = []
    maxanno        = 0
    antable = parentuvdata.table('AN', 1)
    for row in antable:
        if row.nosta > maxanno:
            maxanno = row.nosta
    for i in range(maxanno):
        annumbers.append(i+1)
    tables = parentuvdata.tables
    for table in tables:
        if table[1][-2:] == 'SN' and table[0] >= outputsnversion:
            if vcs.yesno("Can i delete " + table[1][-2:] + " table version " + \
                     str(table[0]) + "? (No will abort pipeline)"):
                parentuvdata.table('SN', table[0]).zap()
            else:
                sys.exit()
    wizuvdata  = WizAIPSUVData(parentuvdata)
    newsntable = wizuvdata.attach_table('SN', outputsnversion, no_term=num_snterms)
    newsntable.keywords['NO_IF']   = num_if
    newsntable.keywords['NO_POL']  = num_pol
    newsntable.keywords['NO_TERM'] = num_snterms
    newsntable.keywords['NO_ANT']  = maxanno
    newimag1  = []
    newimag2  = []
    newdelay1 = []
    newrate1  = []
    newdelay2 = []
    newrate2  = []
    for j in range(num_if):
        newimag1.append(0.0)
        newimag2.append(0.0)
        newdelay1.append(0.0)
        newdelay2.append(0.0)
        newrate1.append(0.0)
        newrate2.append(0.0)
    examplerow.imag1   = newimag1
    examplerow.imag2   = newimag2
    examplerow.delay_1 = newdelay1
    examplerow.delay_2 = newdelay2
    examplerow.rate_1  = newrate1
    examplerow.rate_2  = newrate2
    examplerow.time_interval = solmins/1400.0  

    # Run SPLAT to generate the dataset that will be used
    outputuvdata = AIPSUVData(splituvdata.name,"scntt",1,1)
    if outputuvdata.exists():
        outputuvdata.zap()
    splat =  AIPSTask('splat', version = aipsver)
    splat.indata  = splituvdata
    splat.outdata = outputuvdata
    splat.solint  = solmins 
    splat()

    # Loop through the data, summing all amplitudes of all baselines
    wizuvdata2 = WizAIPSUVData(outputuvdata)
    ampsum = 0
    count  = 0
    for visibility in wizuvdata2:
        for i in range(num_if):
            for j in range(num_pol):
                if visibility.visibility[i][0][j][2] > 0.0:
                    amp = math.sqrt(visibility.visibility[i][0][j][0]**2 + visibility.visibility[i][0][j][1]**2)
                    ampsum += amp
                    count += 1
    avamp = ampsum/count
    # print('avamp %d' %avamp) #*************
    lasttime = -999
    ifampsum = []
    ifcount  = []
    for i in range(num_if):
        ifampsum.append(0)
        ifcount.append(0)
    for visibility in wizuvdata2:
        if visibility.time > lasttime + examplerow.time_interval/2.0:
            if lasttime > 0:
                examplerow.time = lasttime
                newreal1 = []
                newreal2 = []
                newweight1 = []
                newweight2 = []
                for i in range(num_if):
                    correction = 1.0
                    corrweight = 1.0
                    if ifcount[i] > 0:
                        ifampsum[i] /= ifcount[i]
                        correction = math.sqrt(avamp/ifampsum[i])
                        corrweight = 1.0/(correction*correction)
                    newreal1.append(correction)
                    newreal2.append(correction)
                    newweight1.append(corrweight)
                    newweight2.append(corrweight)
                examplerow.real1 = newreal1
                examplerow.real2 = newreal2
                examplerow.weight_1 = newweight1
                examplerow.weight_2 = newweight2
                for annum in annumbers:
                    examplerow.antenna_no = annum
                    newsntable.append(examplerow)
            lasttime = visibility.time
            for i in range(num_if):
                ifampsum[i] = 0
                ifcount[i]  = 0
        for i in range(num_if):
            for j in range(num_pol):
                if visibility.visibility[i][0][j][2] > 0.0:
                    amp = math.sqrt(visibility.visibility[i][0][j][0]**2 + visibility.visibility[i][0][j][1]**2)
                    ifampsum[i] += amp
                    ifcount[i] += 1
    newsntable.close()
    if os.path.exists(outputsntable):
        os.system("rm -f " + outputsntable)
    writetable(parentuvdata, 'SN', outputsnversion, outputsntable)
    print("VLBI PIPE: wizCorrectScint: Scintillation corrections are written in SN%d" %outputsnversion)
    print("VLBI PIPE: wizCorrectScint: Pulsar scintillation correction is done")


def imagejmfit(modelfile, file_dir, psrmodel, pixwindow=20):
    """
    JMFIT: Gaussian fitting to image plane to determine the source parameters
    First, it load image to AIPS
    """
    outdata   = AIPSImage(psrmodel.name, 'JMFIT', 1, 1)
    if outdata.exists():
        outdata.zap()

    fitld_image(modelfile, psrmodel)

    imean = AIPSTask("imean", version = aipsver)
    imean.indata = psrmodel
    imean.inp()
    imean.go()

    lines = imean.message()
    xpix = 128
    ypix = 128
    for line in lines:
        if "Maximum" in line:
            splitline = line.split()
            xpix = float(splitline[4])
            ypix = float(splitline[5])
    print(xpix, ypix)
    blc = [None,0,0]
    trc = [None,0,0]
    blc[1:] = [xpix-pixwindow/2, ypix-pixwindow/2]
    trc[1:] = [xpix+pixwindow/2, ypix+pixwindow/2]
    print(blc)
    print(trc)
    if blc[1] < 0: 
        blc[1] = 1
    if blc[2] < 0: 
        blc[2] = 1
    if trc[1] > int(psrmodel.header.naxis[0]): 
        trc[1] = int(psrmodel.header.naxis[0])-1
    if trc[2] > int(psrmodel.header.naxis[1]): 
        trc[2] = int(psrmodel.header.naxis[1])-1
    
    jmfit = AIPSTask('jmfit', version=aipsver)
    jmfit.indata  = psrmodel
    jmfit.outdata = outdata
    jmfit.niter   = 4000
    jmfit.ngauss  = 1
    jmfit.fwidth[1][1]  = 0
    jmfit.fwidth[1][2]  = 0
    jmfit.fwidth[1][3]  = 0
    jmfit.domax[1:]     = [1]
    jmfit.dopos[1][1]   = 1
    jmfit.dopos[1][2]   = 1
    jmfit.dowidth[1][1] = 1
    jmfit.dowidth[1][2] = 1
    jmfit.dowidth[1][3] = 1
    jmfit.blc = blc
    jmfit.trc = trc
    jmfit.inp()
    jmfit()
    jmfitmassage = jmfit.message()
    file = modelfile.replace('.fits', '.jmfit')
    with open(file, 'w') as f:
        for lines in jmfitmassage:
            f.write(lines) 
            f.write('\n')  


def dbcon(inputuvdatas, outputuvdata):
    """
    DBCON: Concatenate AIPSCat enteries
    """
    errormsg = "ANTENNA tables don't match! Can't DBCON with doarray=1.  " + \
               "This is set so that there are no subarrays, and difmap " + \
               "can still process the data.  Write a manual dbcon to get " + \
               "around this if subarrays are ok."
    ocount = 1
    junko = AIPSUVData('JUNKO','JUNKO',1,1)
    if junko.exists():
        junko.zap()
    dbcon = AIPSTask('dbcon', version = aipsver)
    dbcon.reweight[1:] = [0]
    dbcon.dopos[1:][1:] = [0]
    dbcon.doarray = 1
    dbcon.fqtol = -1
    missingan1 = []
    inputuvdatas.sort(key=lambda a: len(a.antennas), reverse=True)
    for i in range(len(inputuvdatas)-1):
        for j in range(i+1,len(inputuvdatas)):
            antable1 = WizAIPSUVData(inputuvdatas[i]).table('AN', 1)
            antable2 = WizAIPSUVData(inputuvdatas[j]).table('AN', 1)
            if len(inputuvdatas[i].antennas) < len(inputuvdatas[j].antennas):
                print(errormsg)
                sys.exit()
            for a in inputuvdatas[j].antennas:
                index1 = -1
                index2 = -1
                for row in antable1:
                    if row.anname.strip() == a.strip():
                        index1 = row.nosta
                        break
                for row in antable2:
                    if row.anname.strip() == a.strip():
                        index2 = row.nosta
                        break
                if index1 < 0:
                    missingan1.append(a)
                if index1 != index2 and index1 >= 0 and index2 >= 0:
                    print(a)
                    print(index1, index2)
                    print(inputuvdatas[i].antennas)
                    for row in antable1:
                        print(row)
                    print(inputuvdatas[j].antennas)
                    for row in antable2:
                        print(row)
                    print(errormsg)
                    sys.exit()
                if index1 >= 0 and index2 >= 0:
                    for row1 in antable1:
                        if row1.anname.strip() == a.strip():
                            staxof = row1.staxof
                            x = row1.stabxyz[0]
                            y = row1.stabxyz[1]
                            z = row1.stabxyz[2]
                            break
                    for row in antable2:
                        if row.anname.strip() == a.strip():
                            print("UPDATING " + a)
                            staxdiff = row.staxof - staxof
                            xdiff = row.stabxyz[0] - x
                            ydiff = row.stabxyz[1] - y
                            zdiff = row.stabxyz[2] - z
                            if math.fabs(staxdiff) > 0.5 or math.fabs(xdiff) > 0.5 or \
                               math.fabs(ydiff) > 0.5 or math.fabs(zdiff) > 0.5:
                                print(row)
                                print(row1)
                                print(errormsg)
                                sys.exit()
                            row.staxof = staxof
                            row.stabxyz[0] = x
                            row.stabxyz[1] = y
                            row.stabxyz[2] = z
                            row.update()
            print("VLBI PIPE: dbcon: Missing antenna list,", missingan1)         # Ashish has added this line******
            for a in missingan1:
                for row in antable2:
                    if row.anname.strip() == a.strip():
                        antable1.append(row)
            antable1.close()
            antable2.close()
            missingan1.clear()         # Ashish has added this line******
    # Merge the two antenna tables manually
    dbcon.indata  = inputuvdatas[0]
    dbcon.in2data = inputuvdatas[1]
    print(inputuvdatas[0].antennas)
    print(inputuvdatas[1].antennas)
    dbcon.outdata = junko
    dbcon.inp()
    dbcon()
    for uvdata in inputuvdatas[2:]:
        ocount += 1
        dbcon.indata  = junko
        dbcon.in2data = uvdata
        print(junko.antennas)
        print(uvdata.antennas)
        junko = AIPSUVData('JUNKO','JUNKO',1,ocount)
        if junko.exists():
            junko.zap()
        dbcon.outdata = junko
        dbcon.inp()
        dbcon()
    print(junko.antennas)

    if outputuvdata.exists():
        if vcs.yesno(f"Delete existing UV dataset {outputuvdata.name} ? (No will abort pipeline)"):
            outputuvdata.zap()
        else: sys.exit()
    uvsrt = AIPSTask('UVSRT', version = aipsver)
    uvsrt.indata = junko
    uvsrt.outdata = outputuvdata
    uvsrt.sort = 'TB'
    uvsrt.inp()
    uvsrt()
    for i in range(ocount):
        junko = AIPSUVData('JUNKO','JUNKO',1,i+1)
        junko.zap()


##### Run quack for the specified time range #########
def quack(uvdata, sources, type, seconds):
    quack = AIPSTask('quack', version = aipsver)
    quack.indata = uvdata
    quack.sources[1:] = sources
    quack.timerang[1:] = [0,0,0,0,0,0,0,0]
    quack.antennas[1:] = [0]
    # for i in range(len(antennas)):
    #     quack.antennas[i+1] = antennas[i]
    quack.flagver = 1
    quack.opcode = type #E.g. BEG OR END OR ENDB
    quack.reason = 'QUACK'
    quack.aparm[1] = 0
    quack.aparm[2] = seconds/60.0
    quack.aparm[3] = 0.5
    quack.inp()
    quack()


def difmap_maptarget(uvfile, pixsize, dogaussian, mapsize=1024, uvweightstr="0,-1", 
         uvaverstr='20,True', uvtaperstr='0.99,1000'):

    uvfilesplit = uvfile.split('/')
    os.chdir(os.path.dirname(uvfile))
    difmap = subprocess.Popen('difmap', stdin=subprocess.PIPE, text=True)
    
    difmap.stdin.write("float pkflux\n")
    difmap.stdin.write("float peakx\n")
    difmap.stdin.write("float peaky\n")
    difmap.stdin.write("float finepeakx\n")
    difmap.stdin.write("float finepeaky\n")
    difmap.stdin.write("float rmsflux\n")
    difmap.stdin.write("integer ilevs\n")
    difmap.stdin.write("float lowlev\n")

    difmap.stdin.write("obs " + uvfile + "\n")
    difmap.stdin.write("select i\n")
    difmap.stdin.write("uvaver " + uvaverstr + "\n")

    psfile = uvfilesplit[-1].replace('aips.fits', 'radplot.png')
    os.system("rm -f " + psfile)
    difmap.stdin.write("device %s/PNG\n" % psfile)
    difmap.stdin.write("radplot\n")

    difmap.stdin.write("mapsize " + str(mapsize) + "," + str(pixsize) + "\n")
    difmap.stdin.write("uvweight " + uvweightstr + "\n")
    difmap.stdin.write("uvtaper " + uvtaperstr + "\n")
    difmap.stdin.write("peakx = peak(x,abs)\n")
    difmap.stdin.write("peaky = peak(y,abs)\n")
    difmap.stdin.write("unshift\n")
    difmap.stdin.write("shift -peakx,-peaky\n")

    difmap.stdin.write("pkflux = peak(flux,abs)\n")
    difmap.stdin.write("finepeakx = peak(x,abs)\n")
    difmap.stdin.write("finepeaky = peak(y,abs)\n")
    difmap.stdin.write("clrmod true\n")
    if dogaussian:
        difmap.stdin.write("addcmp pkflux, true, finepeakx, finepeaky, true, 0.5, " + \
                        "true, 1, true, 0, true\n")
    else:
        difmap.stdin.write("addcmp pkflux, true, finepeakx, finepeaky, true, 0, " + \
                        "false, 1, false, 0, false, 0, 0, 0\n")
        
    difmap.stdin.write("modelfit 50\n")
    difmap.stdin.write("restore\n")
    difmap.stdin.write("selfcal\n")
    # difmap.stdin.write("modelfit 40\n")
    difmap.stdin.write("gscale\n")
    difmap.stdin.write("selfcal true,true,20\n")
    difmap.stdin.write("selfcal\n")
    difmap.stdin.write("selfcal true,true,10\n")

    difmap.stdin.write("rmsflux = imstat(rms)\n")
    difmap.stdin.write("pkflux = peak(flux,abs)\n")
    difmap.stdin.write("ilevs = pkflux/rmsflux\n")
    difmap.stdin.write("lowlev = 300.0/ilevs\n")
    difmap.stdin.write("loglevs lowlev\n")

    psfile = uvfilesplit[-1].replace('aips.fits', 'difmap.ps')
    os.system("rm -f " + psfile)
    difmap.stdin.write("device %s/CPS\n" % psfile)
    difmap.stdin.write("mappl cln\n")
    imagename = uvfilesplit[-1].replace('aips', 'difmap')
    difmap.stdin.write("wmap " + imagename + "\n")

    difmap.stdin.write("exit\n\n")


def vplot_difmap(uvfile, numifs, nplot=4, uvaverstr='0,True'):
    uvfilesplit = uvfile.split('/')
    os.chdir(os.path.dirname(uvfile))
    difmap = subprocess.Popen('difmap', stdin=subprocess.PIPE, text=True)
    difmap.stdin.write("obs " + uvfile + "\n")
    difmap.stdin.write("select i\n")
    difmap.stdin.write("uvaver " + uvaverstr + "\n")
    for i in range(numifs):
        psfile = uvfilesplit[-1].replace('aips.fits', f"vplot_if{i+1}.ps")
        os.system("rm -f " + psfile)
        difmap.stdin.write("device %s/PS\n" % psfile)
        difmap.stdin.write(f"vplot {nplot},,{i+1}\n")
    difmap.stdin.write("exit\n\n")

#end===================================================
##
###
