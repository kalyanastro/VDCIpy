#!/usr/bin/python3
##########################################################################################
## To compute the systematic error contribution in astrometry
## Hao Ding, Ashish Kalyan
##########################################################################################

import os, yaml, math, re, argparse, sys
from argparse import RawTextHelpFormatter
import numpy as np
np.set_printoptions(suppress=True)
import vlbi_calibration_support as vcs
import astro_utils as ast


def ellipse2xy(LA, SA, PA):
    """
    LA ->full long axis; 
    SA ->full short axis; 
    PA ->position angle
    """
    b = 0.5 * LA
    a = 0.5 * SA
    d = PA * math.pi / 180
    # for y_max:
    lamda = (a**2 - b**2) * math.sin(d) * math.cos(d) / ((b*math.cos(d))**2 + (a*math.sin(d))**2) #lamda=x/y
    y_max = ((lamda * math.cos(d) + math.sin(d)) / a)**2 + ((math.cos(d) - lamda * math.sin(d)) / b)**2
    y_max = y_max**(-0.5)
    # for x_max:
    mu = (a**2 - b**2) * math.sin(d) * math.cos(d) / ((a*math.cos(d))**2 + (b*math.sin(d))**2)
    x_max = ((math.cos(d) + mu*math.sin(d)) / a)**2 + ((mu*math.cos(d) - math.sin(d)) / b)**2
    x_max = x_max**(-0.5)
    return x_max, y_max

def deprojectbeam2xy(LA, SA, PA):
    [x_max, y_max] = ellipse2xy(LA, SA, PA)
    xyratio = x_max / y_max
    x = (math.pi * xyratio * LA/2 * SA/2)**0.5
    y = (math.pi * LA/2 * SA/2/xyratio)**0.5
    return x, y

def dms2deg(dec):  
    dec1 = dec.split(':')
    decd = float(dec1[0])   # Dec of source [hour]
    decm = float(dec1[1])   # Dec of source [min]
    decs = float(dec1[2])   # Dec of source [s]
    decs1 = abs(decd)*60*60 + decm*60 + decs  # in [arcs]
    decdeg = decs1/3600
    if decd < 0:    
        decdeg = -(decs1/3600)
    return decdeg

def mas2ms(ErrArray_mas, Dec):
    d = dms2deg(Dec)*math.pi/180
    ErrArray_ms = ErrArray_mas/15/math.cos(d)
    return ErrArray_ms

def ms2mas(ErrArray_ms, Dec):
    d = dms2deg(Dec)*math.pi/180
    d = np.array(d)
    ErrArray_ms = np.array(ErrArray_ms)
    ErrArray_mas = ErrArray_ms * 15 * np.cos(d)
    return ErrArray_mas

def separation(RA1, Dec1, RA2, Dec2): 
    """
    Functionality
    -------------
    calculate angular separation (in arcmin) given RAs/Decs in hms/dd:mm:ss.ssss format;
    """
    RA0 = np.array([dms2deg(str(RA1)), dms2deg(str(RA2))]) # multipy by 15 done below
    Dec = np.array([dms2deg(str(Dec1)), dms2deg(str(Dec2))])
    Dec_rad = Dec*math.pi/180.
    RA = RA0*15 #hour to deg
    RA_rad = RA*math.pi/180
    cos_sep = np.cos(Dec_rad[0])*np.cos(Dec_rad[1])*np.cos(RA_rad[0]-RA_rad[1]) + \
                        np.sin(Dec_rad[0])*np.sin(Dec_rad[1])
    sep = math.acos(cos_sep)
    sep = sep*180/math.pi*60 # rad to arcmin
    return sep


class expno_sysErr:
    paraA = 0.001
    paraB = 0.6
    def __init__(self, obscode, targetname, RA_target, Dec_target, RA_prIBC, Dec_prIBC, sumfile, prIBCstatsfile):
        self.prIBCstatsfile = prIBCstatsfile
        self.RA_target  = RA_target
        self.Dec_target = Dec_target
        self.targetname = targetname
        self.RA_prIBC   = RA_prIBC
        self.Dec_prIBC  = Dec_prIBC
        self.sumfile    = sumfile
        self.obscode    = obscode


    def sysErr(self):
        """
        Output parameter
        ----------------
        sysErrRA  : float [mas]
        sysErrDec : float [mas]
        sysErrRA_ms : float [ms]
        """
        [beamRA, beamDec] = self.beam_in_RA_Dec()
        delta_sys1  = self.delta_sys()    
        sysErrRA    = beamRA*delta_sys1  # [mas]
        sysErrDec   = beamDec*delta_sys1 # [mas]
        sysErrRA_ms = mas2ms(sysErrRA, self.Dec_target)
        return sysErrRA, sysErrDec, sysErrRA_ms
    

    def delta_sys(self):
        """
        Input parameters
        ----------------
        RA_target : str
            RA of the target in HH:MM:SS.SSSS
        Dec_prIBC : str
            Dec of the primary inbeam calibrator in dd:mm:ss.ssss
        """
        RA_target  = self.RA_target
        Dec_target = self.Dec_target
        RA_prIBC   = self.RA_prIBC
        Dec_prIBC  = self.Dec_prIBC

        sep = separation(RA_target, Dec_target, RA_prIBC, Dec_prIBC) #unit: arcmin
        print(f"\nSeparation between pulsar and IBC [arcmin]: {sep}")

        Av_cscEl1 = self.Av_cscEl()
        [junk1, junk2, junk3, SNprIBC] = self.targetbeam()
        delta_sys = self.paraA * sep * Av_cscEl1 + self.paraB / SNprIBC
        return delta_sys
    

    def beam_in_RA_Dec(self):
        [beamPA, beamSA, beamLA, junk1] = self.targetbeam()
        #full-width deprojection on RA/Dec from beam
        [beamRA, beamDec] = deprojectbeam2xy(beamLA,beamSA,beamPA)
        return beamRA, beamDec


    def targetbeam(self): #also get SNprIBCs, using prIBC statsfile,
        """
        Note
        ----
        1. statsfile for divided IBC fitsfile is not used, as the image S/N normally increases 
            after dividing the model.
        """
        print(f"\nGetting beam and SNR info from {self.prIBCstatsfile.split('/')[-1]}")
        lines = open(self.prIBCstatsfile).readlines()[-10:]
        for line in lines:
            if 'S/N' in line:
                SNprIBC = line.split(':')[-1].strip()
                SNprIBC = float(SNprIBC)
            if 'beam' in line:
                line = line.split('beam')[-1].strip().split(' ')
                beamPA   = float(line[-2])
                beamsize = line[0]
                beamLA   = float(beamsize.split('x')[-1])
                beamSA   = float(beamsize.split('x')[0])
        return beamPA, beamSA, beamLA, SNprIBC
    
    
    def Av_cscEl(self):
        targetname = self.targetname
        obscode = self.obscode
        TelAv_csc_Els = np.array([])
        
        if not targetname[-1].isdigit():
            targetname = targetname[:-1]
        keyword = targetname
        if obscode == "BD174":
            keyword = targetname[0:5] + 'PT'

        print(f"Getting Av_cscEl info from {self.sumfile}")
        f = open(self.sumfile)
        lines = f.readlines()
        f.close()
        
        for line in lines:
            if keyword in line:
                if ':' in line.split(keyword)[0]:
                    try:
                        elevations = line.split('-')[-1].strip().split('    ')
                        elevations = list(map(float, elevations))
                    except ValueError:
                        try:
                            elevations = line.split('-')[-1].strip().split('   ')
                            elevations = list(map(float, elevations))
                        except ValueError:
                            continue
                    elevations = np.asarray(elevations)
                    elevations_rad = elevations*math.pi/180
                    csc_elevations = (np.sin(elevations_rad))**(-1)
                    TelAv_csc_El = np.average(csc_elevations) #telescope averaged csc(el) for each scan
                    TelAv_csc_Els = np.append(TelAv_csc_Els, TelAv_csc_El)
        Av_cscEl = np.average(TelAv_csc_Els) #averaged csc(el) over telescopes and scans
        return Av_cscEl
    

def update_raerr_decerr_columns_in_copy(source_file, dest_file, ra_syserr, dec_syserr):
    """
    update raerr and decerr columns enteries and write a new pmparin file
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
        columns[2] = ra_syserr[i]
        columns[4] = dec_syserr[i]
        updated_data.append('  '.join(columns))
    
    # Create a copy with updated content
    with open(dest_file, 'w') as f:
        for line in header:
            f.write(line)
        for line in updated_data:
            f.write(line + '\n')


def jmfit_stats_reader(jmfitstatsfile):
    """
    It reads the stats file and get the ra and dec of the image center
    Return: rastring, decstring
    """
    statsfile = open(jmfitstatsfile, 'r')
    lines = statsfile.readlines()
    statsfile.close()
    for line in lines:
        if 'Centre RA' in line:
            rastring = line.split(' ')[-1].strip()
        if 'Centre Dec' in line:
            decstring = line.split(' ')[-1].strip()
    return rastring, decstring


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make a PMPAR input file from preliminary file', 
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument('-pf', '--pmparfile', metavar='', help='Preliminary pmpar file')
    args = parser.parse_args()

    pmparfile  = args.pmparfile
    yamldir    = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "/examples"
    expconfig  = yaml.load(open(yamldir + '/' + 'expconfig.yaml'), yaml.SafeLoader)
    obscode    = expconfig['obscode'].upper()
    projectdir = expconfig['projectdir']
    targetname = expconfig['target']
    ibcname    = expconfig['primaryinbeam']
    sumfiledir = projectdir + f"support_files/sumfiles"
    target  = expconfig['target']
    band = str(expconfig['band'])
    bif = str(expconfig['bif'])
    eif = str(expconfig['eif'])
    suffix = band
    if expconfig['splitifs']:
        suffix =  bif + eif
    original_stdout = sys.stdout

    decimalyear, ras, raerrs, decs, decerrs, epochs, mjd = vcs.pmparin_reader(pmparfile)  # target

    # write a log file in result directory
    logfilepath = os.path.dirname(pmparfile) + f"/{target}.sysErr.log"
    logfile = open(logfilepath, 'w')
    sys.stdout = ast.Logger(logfile)

    ibcs = expconfig['ibcs']
    prIBCkey = 'ibc' + str(ibcs.index(expconfig['primaryinbeam']) + 1)
    ibcfitsfile = expconfig[prIBCkey + 'modelfile']
    ibcstatsfile = ibcfitsfile.replace(".fits", ".jmfit.stats")
    print(ibcstatsfile)
    if not os.path.exists(ibcstatsfile):
        os.system(f"jmfitfromfile.py {ibcfitsfile} {ibcstatsfile.split('/')[-1]}")
    raibcstring, decibcstring = jmfit_stats_reader(ibcstatsfile)
    print("IBC coordinates: ", [raibcstring, decibcstring])

    delta_ra, delta_dec = [], []
    for i, epoch in enumerate(epochs):
        sumfile = sumfiledir + f"/{obscode.lower()}{epoch.lower()}.sum"
        if os.path.exists(sumfile):
            pass
        else:
            print(".sum file not found")
            sys.exists() 
        ibcstatsfile = projectdir + f"{epoch}/source_models/{epoch}{prIBCkey}{suffix}_{ibcname}_model/" + \
                                        f"{epoch}{prIBCkey}{suffix}_{ibcname}_difmap.jmfit.stats"
        if os.path.exists(ibcstatsfile):
            pass
        else:
            print("ibcstats file not found")
            sys.exists() 
       
        print(ibcstatsfile)
        syserr_class = expno_sysErr(obscode, targetname, ras[i], decs[i], raibcstring, decibcstring,
                                sumfile, ibcstatsfile)
        sysErrRA, sysErrDec, sysErrRA_ms = syserr_class.sysErr()

        delta_ra.append(f"{np.round(np.sqrt(float(raerrs[i])**2 + (sysErrRA_ms/1000)**2), 8):.8f}")
        delta_dec.append(f"{np.round(np.sqrt(float(decerrs[i])**2 + (sysErrDec/1000)**2), 8):.8f}")

    dest_file = pmparfile.replace("pmpar.in.preliminary", "pmpar.in")
    update_raerr_decerr_columns_in_copy(pmparfile, dest_file, delta_ra, delta_dec)
    print("sysErr.py: Appears to have ended successfully")
    logfile.close()
    sys.stdout = original_stdout
# end=======================
##
###