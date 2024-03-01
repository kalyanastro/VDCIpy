"""
Ashish Kalyan, Hao Ding, and Adam Deller
"""

import math, sys, re, os
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import astro_utils as ast


class VexScan:
    def __init__(self, scanname, starttime, stoptime, source, modename, vlaapdet):
        self.scanname  = scanname
        self.starttime = starttime
        self.stoptime  = stoptime
        self.source    = source
        self.modename  = modename
        self.startmjd  = getVexMJD(starttime)
        self.stopmjd   = getVexMJD(stoptime)
        self.vlaautophasedetermine = vlaapdet

    def incStopTime(self, incsecs):
        refmjd = int(self.stopmjd)
        doy = self.getStartDOY()
        self.stopmjd += incsecs/86400.0
        fracmjd = self.stopmjd - refmjd
        if fracmjd > 1:
            doy += 1
            fracmjd -= 1
        hh = int(fracmjd*24.0)
        mm = int(fracmjd*1440.0 - 60*hh)
        ss = int(fracmjd*86400.0 - (3600*hh + 60*mm))
        self.stoptime = "%s%03dd%02dh%02dm%02ds" % (self.stoptime[:5], doy, hh, mm, ss)

    def setStopTime(self, stoptime):
        self.stoptime = stoptime
        self.stopmjd = getVexMJD(stoptime)

    def getStartDOY(self):
        return int(self.starttime[5:8])

    def getStartHour(self):
        return int(self.starttime[9:11])

    def getStartMinute(self):
        return int(self.starttime[12:14])

    def getStartSecond(self):
        return int(self.starttime[15:17])

    def getStopDOY(self):
        return int(self.stoptime[5:8])

    def getStopHour(self):
        return int(self.stoptime[9:11])

    def getStopMinute(self):
        return int(self.stoptime[12:14])

    def getStopSecond(self):
        return int(self.stoptime[15:17])
    
    
class Source:
    def __init__(self, name, ra, dec):
        self.name = name
        self.ra   = ra
        self.dec  = dec

    def getSepArcmin(self, rarad, decrad):
        sinsqdecdiff = math.sin((decrad - self.dec)/2.0)
        sinsqdecdiff = sinsqdecdiff*sinsqdecdiff
        sinsqradiff  = math.sin((rarad - self.ra)/2.0)
        sinsqradiff  = sinsqradiff*sinsqradiff
        return 180.0*60.0*2*math.asin(math.sqrt(sinsqdecdiff + math.cos(decrad)*math.cos(self.dec)*sinsqradiff))/math.pi


def getVexMJD(vexdate):
    """
    CONVERT A VEX HEADER DATE STRING TO MJD
    """
    year = int(vexdate[0:4])
    doy  = int(vexdate[5:8])
    hh   = int(vexdate[9:11])
    mm   = int(vexdate[12:14])
    ss   = int(vexdate[15:17])
    jan1mjd = ast.ymd2mjd(year, 1, 1)
    return jan1mjd + (doy - 1) + float(hh)/24.0 + float(mm)/1440.0 + float(ss)/86400.0


def getrefdate(vexfile):
    """
    PARSE A VEX FILE AND GET THE REFERENCE DATE
    """
    refdate = ""
    vexlines = open(vexfile).readlines()
    for line in vexlines:
        if "exper_nominal_start" in line:
            splitline = line.split('=')
            vexmjd = getVexMJD(splitline[-1][:-1])
            year, month, day, hour, minute, second = ast.mjd2ymdhms(vexmjd)
            refdate = "%04d%02d%02d" % (year, month, day)
            break
    return refdate


############### Checks if user wants to do something or not ###############
alwaysuse = False
def yesno(prompt):
    if alwaysuse:
        return True
    ans = ""
    while ans == "":
        ans = input(prompt)
    return (ans[0]=='y' or ans[0]=='Y')


def getvexscaninfo(vexfile):
    """
    PARSE A VEX FILE INTO A LIST OF SCANS
    """
    scanlist = []
    srclist  = []
    vexin    = open(vexfile)
    vexlines = vexin.readlines()
    vexin.close()

    vlaautophasedetermine = False
    atline = 0
    numlines = len(vexlines)
    while atline < numlines and not '$SOURCE;' in vexlines[atline]:
        atline += 1
    if atline == numlines:
        print ("Couldn't find source block in vex file")
        sys.exit()

    atline += 1
    while atline < numlines and not '$' in vexlines[atline]:
        splitline = vexlines[atline].split()
        if len(splitline) > 1 and splitline[0] == 'def':
            atline += 1
            while atline < numlines and not 'source_name' in vexlines[atline]:
                atline += 1
            if atline == numlines:
                print ("Got lost in vex source block")
                sys.exit()
            srcname = vexlines[atline].split()[-1][:-1]
            atline += 1
            while atline < numlines and not 'J2000' in vexlines[atline]:
                atline += 1
            if atline == numlines:
                print ("Got lost in vex source block")
                sys.exit()
            splitline = vexlines[atline].split()
            ra  = splitline[2][:-1]
            dec = splitline[5][:-1]
            hh = float(ra[:2])
            mm = float(ra[3:5])
            ss = float(ra[6:-1])
            rarad = (hh/24.0 + mm/1440.0 + ss/86400.0)*2.0*math.pi
            if dec[0] == '-':
                multiplier = -1.0
                dd = float(dec[1:3])
                mm = float(dec[4:6])
                ss = float(dec[7:-1])
            else:
                multiplier = 1.0
                dd = float(dec[:2])
                mm = float(dec[3:5])
                ss = float(dec[6:-1])
            decrad = multiplier*(dd + mm/60.0 + ss/3600.0)*math.pi/180.0
            srclist.append(Source(srcname, rarad, decrad))
        atline += 1
    print ("Found " + str(len(srclist)) + " sources")

    while atline < numlines and not '$SCHED;' in vexlines[atline]:
        atline += 1
    if atline == numlines:
        print ("Couldn't find scan block in vex file")
        sys.exit()

    lastsec = 0
    while atline < numlines:
        while atline < numlines and (not 'start' in vexlines[atline] or '*' in vexlines[atline]):
            if 'sec:' in vexlines[atline]:
                lastsec = int(vexlines[atline].split('sec:')[1].strip())
            if "VLA:AUTOPHASE_DETERMINE" in vexlines[atline]:
                vlaautophasedetermine = True
            if "VLA:AUTOPHASE_APPLY" in vexlines[atline] or "VLA:AUTOPHASE_OFF" in vexlines[atline]:
                vlaautophasedetermine = False
            atline += 1
        if atline == numlines:
            break
        scanname = "unknown"
        if "scan" in vexlines[atline-1]:
            scanname = vexlines[atline-1].split()[1][:-1]
        splitline = vexlines[atline].split()
        starttime = splitline[0].split('=')[1][:-1]
        if len(scanlist) > 0:
            scanlist[-1].incStopTime(lastsec)
        stoptime = starttime
        scansourcename = splitline[2].split('=')[1][:-1]
        scanmodename = splitline[1].split('=')[1][:-1]
        source = None
        for src in srclist:
            if src.name == scansourcename:
                source = src
                break
        if source == None:
            print ("Couldn't find source " + scansourcename + " in source list!")
            sys.exit()
        scanlist.append(VexScan(scanname, starttime, stoptime, source, scanmodename, vlaautophasedetermine))
        atline += 1
    if len(scanlist) > 0:
        scanlist[-1].incStopTime(lastsec)
    if "bd179" in vexfile:
        for source in srclist:
            if source.name[-2:] == "AP":
                source.name = source.name + "T"
    print ("Found " + str(len(scanlist)) + " scans in vex file " + vexfile)
    return scanlist


def getAIPSMJD(aipsdate):
    """
    CONVERT AN AIPS HEADER DATE STRING TO MJD
    """
    year  = int(aipsdate[0:4])
    month = int(aipsdate[5:7])
    day   = int(aipsdate[8:10])
    return ast.ymd2mjd(year, month, day)
    
    
def obsyear(vexfile):
    """
    Input: vexfile
    Output: tuple (observation epoch, decimal year, modified julian date)
    """
    vexin    = open(vexfile)
    vexlines = vexin.readlines()
    vexin.close()
    atline = 0
    while 'exper_name' not in vexlines[atline]:
        atline += 1
    epoch = re.split(' +', vexlines[atline].strip())[-1][:-1]
    while 'year' not in vexlines[atline]:
        atline += 1
    if 'year' in vexlines[atline]:
        dateline = re.split(' +', vexlines[atline].strip())
        year = int(dateline[3][:-1])
        doy  = dateline[4]
        if int(year)%4 == 0:
            decimalyear = year + int(doy)/366
        else:
            decimalyear = year + int(doy)/365
    while 'MJD' not in vexlines[atline]:
        atline += 1
    mjd = vexlines[atline].split(':')[-1].strip()
    return epoch, decimalyear, mjd


def get_obs_start_end_doy_year(vexfile):
    """
    Get the observation start and stop day of year (doy) and year
    Return: int: start_doy, start_year, end_doy, end_year
    """
    vexin = open(vexfile)
    vexlines = vexin.readlines()
    vexin.close()
    startfound = False
    stopfound = False
    for line in vexlines:
        if "exper_nominal_start" in line:
            splitline = line.split('=')
            syear     = int(splitline[1][:4])
            sdoy      = int(splitline[1][5:8])
            startfound = True
        elif "exper_nominal_stop" in line:
            splitline = line.split('=')
            eyear     = int(splitline[1][:4])
            edoy      = int(splitline[1][5:8])
            stopfound = True
        if startfound and stopfound:
            break
    return sdoy, syear, edoy, eyear


def on_src_time (source, listrfile):
    """
    Determine the on-source time from the listr file
    Return: on-source time [min]
    """
    listrr = open(listrfile, 'r')
    lines  = listrr.readlines()
    listrr.close()

    atline = 0
    while "Scan summary" not in lines[atline]:
        atline += 1
    atline+2

    starttimes = []
    stoptimes  = []
    while "Source summary" not in lines[atline]:
        if source in lines[atline]:
            starttimes.append(re.split(' +', lines[atline].strip())[6])
            stoptimes.append(re.split(' +', lines[atline].strip())[8])
        atline += 1

    startmins  = np.array(timestring2min(starttimes))
    stopmins   = np.array(timestring2min(stoptimes))
    onsrctimes = stopmins - startmins
    tot = 0
    for t in onsrctimes:
        tot += t
    return tot


def src_coord_from_listrfile(listrfile, source):
    """
    Read LISTR file and get a particular source coordinates
    """
    file = open(listrfile, 'r')
    lines = file.readlines()
    numlines = len(lines)
    atline = 0
    while 'Source summary' not in lines[atline]:
        atline += 1
    
    while source not in lines[atline] and atline < numlines-1:
        atline += 1

    if source in lines[atline]:
        ra = re.split('\s+', lines[atline].strip())[5]
        dec = re.split('\s+', lines[atline].strip())[6]
        return ra, dec
    else: 
        print("Source not found in LISTR file")


def jmfitout_reader(jmfitfile):
    """
    Input: JMFIT message file
    Output: str: ra, raerr, dec, decerr
    """
    jmfitin = open(jmfitfile, 'r')
    lines   = jmfitin.readlines()
    jmfitin.close()
    atline = 0
    while 'Solution from JMFIT' not in lines[atline]:
        atline += 1
    while 'RA' not in lines[atline]:
        atline += 1
    if 'RA' in lines[atline]:
        raline = re.split(' +', lines[atline].strip())
        ra = raline[2] + ':' + raline[3] + ':' + f"{float(raline[4]):.8f}"
        raerr = f"{float(raline[6]):.8f}"
    atline += 1  
    if 'DEC' in lines[atline]:
        decline = re.split(' +', lines[atline].strip())
        dec = decline[2] + ':' + decline[3] + ':' + f"{float(decline[4]):.8f}"
        decerr = f"{float(decline[6]):.8f}"
    return ra, raerr, dec, decerr


def su_src_reader(sufile):
    """
    Input: SU extension table, which is stored using TBOUT task
    Output: list: sources
    """
    sutable = open(sufile)
    lines = sutable.readlines()
    sutable.close()

    atline = 0
    while 'NO_IF' not in lines[atline]:
        atline += 1

    ifs = int(lines[atline].strip()[-1])
    while 'BEGIN*PASS' not in lines[atline]:
        atline += 1

    atline += 1
    srcs = []
    while 'END*PASS' not in lines[atline]:
        src = re.split('\s+', lines[atline].strip())[2][1:]
        srcs.append(src)
        atline += ifs
    return srcs


def src_coord_from_sutable(sutable, source):
    """
    It reads the SU table and then get the coordinates of the given source
    RAEPO/DECEPO, source coordinates at equinox given in epoch
    RAAPP/DECAPP, source coordinates at the start of the observation
    RAOBS/DECOBS, antenna pointing center
    Return: float: raepo, decepo, raapp, decapp, raobs, decobs [deg]
    """
    print("Reading ... ", sutable.split('/')[-1])
    f = open(sutable)
    lines = f.readlines()
    f.close()
    atline = 0

    endpass = False
    for line in lines:
        if source in line:
            sourcenum = lines[atline].strip().split(' ')[0]
            break
        if '***END*PASS***' in line:
            endpass = True
            break
        atline += 1
    if endpass:
        print("VLBI PIPE: Source not found, I'm aborting the pipeline")
        sys.exit()

    while '***BEGIN*PASS***' not in lines[atline]:
        atline += 1

    endpass = False
    for i in range(atline, len(lines)):
        splitline = re.split(r'\s+', lines[i].strip())
        if sourcenum in splitline[0]:
            if splitline[3] != "''":
                raepo = float(splitline[3].replace('D', 'E'))
                decepo = float(splitline[4].replace('D', 'E'))
            
        if '***END*PASS***' in lines[i]:
            endpass = True
            break
        atline += 1
        
    endpass = False
    for i in range(atline+1, len(lines)):
        splitline = re.split(r'\s+', lines[i].strip())
        if sourcenum in splitline[0]:
            if splitline[3] != "''":
                raapp = float(splitline[1].replace('D', 'E'))
                decapp = float(splitline[2].replace('D', 'E'))
                raobs = float(splitline[3].replace('D', 'E'))
                decobs = float(splitline[4].replace('D', 'E'))
            
        if '***END*PASS***' in lines[i]:
            endpass = True
            break
        atline += 1

    return raepo, decepo, raapp, decapp, raobs, decobs


def psrpi_reader(psrpi_file, source):
    """
    It reads the 'psrpi_info.txt' file and get the source info
    Return: str: mjd, rastring, decstring, pmra, pmdec
    """
    f = open(psrpi_file)
    lines = f.readlines()
    lines = [line for line in lines if line.strip() and not line.strip().startswith('#')]
    f.close()
    atline = 0
    numlines = len(lines)
    mjd = lines[0].split('=')[1].strip()

    srcline = re.split('\s+', lines[atline].strip())
    while source not in lines[atline]:
        atline += 1
        if atline == numlines:
            print("VLBI PIPE: source doesn't found")
            sys.exit()
    srcline   = re.split('\s+', lines[atline].strip())
    rastring  = srcline[1]
    decstring = srcline[2]
    pmra  = srcline[3]
    pmdec = srcline[4]
    return mjd, rastring, decstring, pmra, pmdec


def mas2deg(a):
    return (a * 10**(-3)) / (3600)


def position_offset_due2_proper_motion(modelvexfile, vexfile1, psrpi_file, source):
    """
    It computes a shift in the pulsar's position due to the proper motion
    Ashish updated it on 18DEC23, earlier it was not working for t<0 means not back tracing
    the pulsar's position. Now, I have added minus if t<0 that will push source towards the other direction.
    Return: float: rashift, decshift [deg]
    """

    mjd, rastring, decstring, pmra, pmdec = psrpi_reader(psrpi_file, source)
    epoch0, decimal_year0, mjd_obs0 = obsyear(modelvexfile)
    epoch1, decimal_year1, mjd_obs1 = obsyear(vexfile1)

    t = (float(mjd_obs1) - float(mjd_obs0)) / 365.25  # [years]
    # t = (float(mjd_obs1) - float(mjd_obs0)) / 365  # [years]
    rashift  = mas2deg(float(pmra)) * t
    decshift = mas2deg(float(pmdec)) * t
    return rashift, decshift


def shift_image_center(fitsfile, outputfits, rashift, decshift):
    """
    It adds rashift and decshift (deg) to image center
    """
    hdul = fits.open(fitsfile)

    # update fits header
    print(f"FITS header [deg]: ", (hdul[0].header['CRVAL1'], hdul[0].header['CRVAL2']))
    hdul[0].header['CRVAL1'] += rashift
    hdul[0].header['CRVAL2'] += decshift
    print("VLBI PIPE: updated fits header [deg]", (hdul[0].header['CRVAL1'], hdul[0].header['CRVAL2']))
    hdul[1].header['XTENSION'] = "BINTABLE"
    hdul.writeto(outputfits, overwrite=True)
    hdul.close()


def pmparin_reader(inpfile):
    """
    Extract source position information from PMPAR input file
    Return: string: decimalyear, ras [hms], raerrs [sec], decs [dms], 
                    decerrs [arcsec], epochs, mjd
    """
    inpfile = open(inpfile, 'r')
    lines   = inpfile.readlines()
    inpfile.close()
    
    numlines = len(lines)
    atline   = 0
    decimalyear, ras, decs, raerrs, decerrs, epochs, mjd = [],[],[],[],[],[],[]

    while '# Source' not in lines[atline]:
        atline += 1
    for atline in range(atline+1, numlines):
        if lines[atline] != ' ' and lines[atline].lstrip()[0] != '#':
            line  = re.split('\s+', lines[atline].strip())
            decimalyear.append(line[0])
            # print(line[0])
            ras.append(line[1])     # hms
            raerrs.append(line[2])  # sec
            decs.append(line[3])    # dms
            decerrs.append(line[4]) # arcsec
            epochs.append(line[-2][1:])
            mjd.append(line[-1][1:])
    return decimalyear, ras, raerrs, decs, decerrs, epochs, mjd


def pmparout_reader(pmparout):
    """
    The function reads the pmpar output
    Return: string: RA [hms], error_RA [sec], Dec [dms], error_Dec [arcsec], 
                    mu_a [mas/year], error_mu_a, mu_d [mas/year], error_mu_d, pi [mas],
                    error_pi, rchsq, mjd
    """
    rchsq = 0
    f = open(pmparout, 'r')
    lines = f.readlines()
    f.close()
    for line in lines:
        if 'epoch' in line:
            epoch = line.split('=')[1].strip()
        if 'Reduced' in line:
            rchsq = float(line.split('=')[1].strip())
        for estimate in ['mu_a', 'mu_d', 'pi']:
            if estimate in line:
                # print(line)
                value = line.split('=')[-1].split('+')[0].strip()
                exec("%s='%s'" % (estimate.strip(), value), globals())
                error = line.split('+-')[1].strip().split(' ')[0]
                exec("error_%s = %s" % (estimate.strip(), error), globals())
                exec("%s = float(%s)" % (estimate.strip(), estimate.strip()), globals())
        if 'RA' in line:
            RA = line.split('=')[-1].split('+')[0].strip()
            error_RA = float(line.split('+-')[1].strip().split(' ')[0])
        if 'Dec  ' in line:
            Dec = line.split('=')[-1].split('+')[0].strip()
            error_Dec = float(line.split('+-')[1].strip().split(' ')[0])
    return RA, error_RA, Dec, error_Dec, mu_a, error_mu_a, mu_d, error_mu_d, pi, error_pi, rchsq, float(epoch)


def pmpare_reader(pmpare):
    """
    It reads pmpar_e file 
    Return: string: etimes, eras, edecs, eraerr, edecerr, pras, pdecs
    """
    print(f"Reading ... {pmpare}")
    f = open(pmpare)
    elines = f.readlines()
    f.close()
    etimes, eras, edecs, eraerr, edecerr, pras, pdecs = [],[],[],[],[],[],[]
    for line in elines:
        splitline = line.split()
        etimes.append(float(splitline[0]))
        eras.append(float(splitline[1]))
        edecs.append(float(splitline[3]))
        eraerr.append(float(splitline[2]))
        edecerr.append(float(splitline[4]))
        pras.append(float(splitline[5]))
        pdecs.append(float(splitline[6]))
    print("done")
    return etimes, eras, edecs, eraerr, edecerr, pras, pdecs


def pmpart_reader(pmpart):
    """
    It reads pmpar_t file 
    Return: string: ttimes, tras, tdecs
    """
    print(f"Reading ... {pmpart}")
    f = open(pmpart)
    elines = f.readlines()
    f.close()
    ttimes, tras, tdecs = [],[],[]
    for line in elines:
        splitline = line.split()
        ttimes.append(float(splitline[0]))
        tras.append(float(splitline[1]))
        tdecs.append(float(splitline[2]))
    print("done")
    return ttimes, tras, tdecs


def run_pmpar(pmparfile, subtractpm=False):
    """
    Run pmpar with and without proper motion removal
    """
    pmparsplit = pmparfile.split('/')
    os.chdir(os.path.dirname(pmparfile))
    if subtractpm:
        os.system("pmpar " + pmparsplit[-1] + " -om")
    else:
        # os.system("pmpar " + pmparsplit[-1] + " -c > " + pmparsplit[-1].replace(".in", ".full.out"))
        os.system("pmpar " + pmparsplit[-1] + " > " + pmparsplit[-1].replace(".in", ".out"))


def replace_word_in_file(file_name, old_word, new_word):
    """
    Replace a given word of a text file with the new word
    """
    with open(file_name, 'r') as file:
        file_contents = file.read()
        
    file_contents = file_contents.replace(old_word, new_word)
    
    with open(file_name, 'w') as file:
        file.write(file_contents)


def pmparin_copy_updated_ra_dec_columns(source_file, dest_file, new_ra_entries, new_dec_entries):
    """
    Make a copy of pmparin file with updated RA and DEC columns enteries
    """
    with open(source_file, 'r') as f:
        lines = f.readlines()
        atline = 0
        numlines = len(lines)
        header = [] 

        while '# Source' not in lines[atline]:
            header.append(lines[atline])
            atline += 1
        header.append(lines[atline])

        data_lines = []
        for atline in range(atline+1, numlines):
            if lines[atline] != ' ' and lines[atline].lstrip()[0] != '#':
                data_lines.append(lines[atline])

    updated_data = []
    for i, line in enumerate(data_lines):
        columns = re.split(r'\s+', line)
        psr_coord = SkyCoord(new_ra_entries[i], new_dec_entries[i], unit='deg')
        ra_hms = psr_coord.ra.to_string(unit="hour", sep=":", pad=True)
        dec_dms = psr_coord.dec.to_string(unit="degree", sep=":", pad=True, alwayssign=True)
        columns[1] = ra_hms
        columns[3] = dec_dms
        updated_data.append('  '.join(columns))
    
    # Create a copy with updated content
    with open(dest_file, 'w') as f:
        for line in header:
            f.write(line)
        for line in updated_data:
            f.write(line + '\n')


def joint_astrometric_fit(pmparin_joint, pmparin1, pmparin2):
    """
    It make a joint pmpar file from two pmpar files by adjusting the reference position offset
    it will store the resultant output in the joint_fit directory.
    joint.pmpar file will have two data section: top will be pmparin1 section and
    bottom will be pmparin2 section 
    """

    print("pmparin1: ", pmparin1)
    
    pmparoutf1 = pmparin1.replace("pmpar.in", "pmpar.out")
    pmparoutf2 = pmparin2.replace("pmpar.in", "pmpar.out")
    if not os.path.exists(pmparoutf1):
        run_pmpar(pmparin1)
    if not os.path.exists(pmparoutf2):
        run_pmpar(pmparin2)
    
    pmparout1 = pmparout_reader(pmparoutf1)
    pmparout2 = pmparout_reader(pmparoutf2)
    print("MJD: ", pmparout1[-1])
    if not pmparout1[-1] == pmparout2[-1]:
        print("Both pmpar.in should have same ref epoch")
        sys.exit()
    coord1 = ast.radec_format_conversion(pmparout1[0], pmparout1[2])
    coord2 = ast.radec_format_conversion(pmparout2[0], pmparout2[2])
    s1 = SkyCoord(coord1[0], coord1[1], frame='icrs')
    s2 = SkyCoord(coord2[0], coord2[1], frame='icrs')
    print("Reference positions difference [mas]: ", round(s1.separation(s2).to(u.arcsec).value*1000, 4))

    # Calculate the RA and DEC differences
    ra_diff  = (s1.ra.deg - s2.ra.deg)    #deg
    dec_diff = (s1.dec.deg - s2.dec.deg)  #deg
    decimalyear, ras, raerrs, decs, decerrs, epochs, mjd = pmparin_reader(pmparin1)
    ra_update, dec_update = [],[]
    for i in range(len(ras)):
        radeg, decdeg = ast.radec2deg(ras[i], decs[i])
        ra_update.append(radeg - ra_diff)
        dec_update.append(decdeg - dec_diff)
    pmparin_copy_updated_ra_dec_columns(pmparin1, pmparin_joint, ra_update, dec_update)

    # append pmparin2 to the joint.pmpar file
    f1 = open(pmparin_joint, 'a')
    f2 = open(pmparin2)
    lines = f2.readlines()
    f2.close()
    atline = 0
    while '# Source' not in lines[atline]:
        atline += 1
    atline +=1
    for i in range(atline, len(lines)):
        if i==atline:
            f1.write('\n' + lines[i])
        else: f1.write(lines[i])
    f1.close()
    return pmparin_joint


def timestring2min(timestringlist):
    """
    Timestringlist element should be in 'day/hh:mm:ss' form.
    Return: minlist
    """
    minlist = []
    for element in timestringlist:
        timestring = element.split('/')
        hh = int(timestring[1].split(':')[0])
        mm = int(timestring[1].split(':')[1])
        ss = int(timestring[1].split(':')[2])
        minlist.append(hh*60 + mm + ss/60)
    return minlist
#end=================================================
##
###