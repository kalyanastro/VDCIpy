"""
Ashish Kalyan and Adam Deller
"""

import math, sys
import numpy as np
from astropy.time import Time
from astropy.constants import c, pc
from astropy import units as u

light_speed = c.value
parsec = pc.value


class Logger:
    """
    used to make the log file
    """
    def __init__(self, logfile):
        self.terminal = sys.stdout
        self.logfile = logfile
    def write(self, message):
        self.terminal.write(message)       
        self.logfile.write(message)
    def flush(self):
        self.terminal.flush()
        self.logfile.flush()
    def isatty(self):
        return False


def mas2rad(x):
    return x * (1e-3) * (1 / 3600) * (np.pi / 180)


def radec_format_conversion(ra, dec):
    """
    Convert ra/dec (:) to hms/dms
    """
    rasplit = ra.split(':')
    decsplit = dec.split(':')
    rahms = rasplit[0] + 'h' + rasplit[1] + 'm' +rasplit[2] + 's'
    decdms = decsplit[0] + 'd' + decsplit[1] + 'm' + decsplit[2] + 's'
    return rahms, decdms


def decimal_year_to_mjd(decimal_year):
    """
    Convert decimal year to a string representation of a Gregorian date
    Return: mjd
    """
    base_year = int(decimal_year)
    remainder = decimal_year - base_year
    start_of_year = Time(f'{base_year}-01-01T00:00:00', format='isot', scale='utc')
    end_of_year = Time(f'{base_year+1}-01-01T00:00:00', format='isot', scale='utc')
    delta_time = (end_of_year - start_of_year).jd
    result_date = start_of_year + remainder * delta_time
    
    # Convert Gregorian date to Julian Date
    jd = result_date.jd
    mjd = jd - 2400000.5
    return mjd


def posdiff(targetrarad, targetdecrad, calrarad, caldecrad):
    sinsqdecdiff = math.sin((targetdecrad-caldecrad)/2.0)
    sinsqdecdiff = sinsqdecdiff*sinsqdecdiff
    sinsqradiff  = math.sin((targetrarad-calrarad)/2.0)
    sinsqradiff  = sinsqradiff*sinsqradiff

    return 2*math.asin(math.sqrt(sinsqdecdiff +
                       math.cos(targetdecrad)*math.cos(caldecrad)*sinsqradiff))


def posradians2string(rarad, decrad):
    """
    Convert RA, DEC (in rad) to string form
    Return: tuple (rastring (:), decstring (:))
    """
    rah  = rarad * 12/math.pi
    rhh  = int(rah)
    rmm  = int(60*(rah - rhh))
    rss  = 3600*rah - (3600*rhh + 60*rmm)
    decd = decrad * 180/math.pi
    decformat = "+%02d:%02d:%010.7f"
    if decd < 0:
        decd = -decd
        decformat = '-' + decformat[1:]
    ddd = int(decd)
    dmm = int(60*(decd - ddd))
    dss = 3600*decd - (3600*ddd + 60*dmm)
    rastring  = "%02d:%02d:%011.8f" % (rhh,rmm,rss)
    decstring = decformat % (ddd, dmm, dss)
    return rastring, decstring


def radec2rad(rastr, decstr):
    """
    Convert RA (h:m:s) and Dec (deg:m:s) to radians
    Return: array([ra_rad, dec_rad])
    """
    # position of source (RA & Dec)
    ra1  =  rastr.split(':')
    rah  =  float(ra1[0])                  # RA of source [hour]
    ram  =  float(ra1[1])                  # RA of source [min]
    ras  =  float(ra1[2])                  # RA of source [s]
    ras1 =  rah*60*60 + ram*60 + ras       # in [s]
    rarad =  np.radians((ras1/3600)*15)    # RA [rad]

    dec1   =  decstr.split(':')
    decd   =  float(dec1[0])                    # Dec of source [hour]
    decm   =  float(dec1[1])                    # Dec of source [min]
    decs   =  float(dec1[2])                    # Dec of source [s]
    decs1  =  abs(decd)*60*60 + decm*60 + decs  # in [arcs]
    decrad =  np.radians(decs1/3600)
    if decd < 0:    
        decrad =  -np.radians(decs1/3600)       # Dec [rad]
    return np.array([rarad, decrad])


def dms2deg(dec):
    """
    Convert dms to decimal degress
    """
    dec1   =  dec.split(':')
    decd   =  float(dec1[0])                    # Dec of source [hour]
    decm   =  float(dec1[1])                    # Dec of source [min]
    decs   =  float(dec1[2])                    # Dec of source [s]
    decs1  =  abs(decd)*60*60 + decm*60 + decs  # in [arcs]
    decdeg =  decs1/3600
    if decd < 0:    
        decdeg =  -decs1/3600       # Dec [rad]
    return decdeg


def radec2deg(rastr, decstr):
    """
    Convert RA (h:m:s) and Dec (deg:m:s) to decimal degree
    Return: array([radeg, decdeg])
    """
    # position of source (RA & Dec)
    ra1   =  rastr.split(':')
    rah   =  float(ra1[0])             # RA of source [hour]
    ram   =  float(ra1[1])             # RA of source [min]
    ras   =  float(ra1[2])             # RA of source [s]
    ras1  =  rah*60*60 + ram*60 + ras  # in [s]
    radeg =  (ras1/3600)*15            # RA [rad]

    dec1   =  decstr.split(':')
    decd   =  float(dec1[0])                    # Dec of source [hour]
    decm   =  float(dec1[1])                    # Dec of source [min]
    decs   =  float(dec1[2])                    # Dec of source [s]
    decs1  =  abs(decd)*60*60 + decm*60 + decs  # in [arcs]
    decdeg =  decs1/3600
    if decd < 0:    
        decdeg =  -decs1/3600       # Dec [rad]
    return np.array([radeg, decdeg])


def mjd2ymdhms(mjd):
    """
    Convert modified Julian day to year, month, day, hour, minutes and seconds
    Return: tuple (year, month, day, hour, minute, second)
    """
    imjd = int(mjd)
    fmjd = mjd - imjd

    j  = imjd + 32044 + 2400001
    g  = int(j/146097)
    dg = j%146097
    c  = int(((int(dg/36524) + 1)*3)/4)
    dc = dg - c*36524
    b  = int(dc/1461)
    db = dc%1461
    a  = int(((int(db/365) + 1)*3)/4)
    da = db - a*365
    y  = g*400 + c*100 + b*4 + a
    m  = int((da*5 + 308)/153) - 2
    d  = da - int(((m + 4)*153)/5) + 122

    year  = y - 4800 + int((m + 2)/12)
    month = (m + 2)%12 + 1
    day   = d + 1
    hour  = int(24*fmjd)
    minute = int(1440*fmjd - 60*hour)
    second = 86400*fmjd - (3600*hour + 60*minute)
    return year, month, day, hour, minute, second


def ymd2mjd(year, month, day):
    """
    Convert year, month and day to modified julian day
    Return: MJD
    """
    return year*367 - int(7*(year + int((month + 9)/12))/4) + int(275*month/9) + day - 678987


def resol_pixelsize(freq, baseline):
    """
    Computes the resolution of synthesis array and pixel size
    https://science.nrao.edu/facilities/vlba/docs/manuals/oss/ang-res
    Input: frequency [GHz], baseline [km]
    Return: resolution [mas], pixel size [mas]
    """

    waveln = 100*(light_speed/(freq*1e9))   # cm
    reso   = np.round(2063*(waveln/baseline), 2) # in arcsecond
    print('Resolution/beam size [mas]: ', reso)

    # pixel size atleast should be 1/3 of the beam
    pix_size = np.round(reso/6, 2)
    print('Pixel/cell size (six pixels in synthesized beam) [mas]: ', pix_size)
    return reso, pix_size


def data_rate(ifs, bw, npol, nbit_sampling, t_acc):
    """
    Input: IFS: number of IFs, 
           BW: bandwidth of observation [MHz]
           NPOL: Number of polarization
           NBIT_SAMPLING: n-bit sampling
           T_ACC: integration time [s]
    Return: data rate [Mbps], product of all the inputs
    """
    return ifs * bw * npol * nbit_sampling * t_acc


def sensitivity(sefd, neeta, nant, bw, tint):
    """
    Compute the image RMS
    Input: System Equivalent Flux Density (SEFD [Jy]), system efficieny factor (neeta)
           nuber of antennas (nant), bandwidth (bw [Hz]), on-source time (tint [sec])
    """
    return sefd/(neeta*np.sqrt(nant*(nant - 1)*bw*tint))


def parallax2dist(parallax, upper_error, lower_error):
    """
    Compute distance from the parallax value
    if parallax [mas] > distance [kpc]
    Return: distance, uperror, loerror
    """
    distance = 1 / parallax

    # Compute Upper/lower Error
    upper_distance = 1 / (parallax - upper_error)
    lower_distance = 1 / (parallax + lower_error)
    uperror = upper_distance - distance
    loerror = distance - lower_distance
    return distance, uperror, loerror

def transverse_velocity(pm_alpha, pm_delta, distance):
    """
    Determine the transverse valocity from the proper motion and distance
    pm_alpha/delta [mas/yr]
    distance [pc]
    return: transverse velocity [km/s]
    """
    yrtosec = 1 * u.year.to(u.s)

    proper_motion = np.sqrt(pm_alpha**2 + pm_delta**2)
    mu_rad = mas2rad(proper_motion)       # [radians/yr]

    # Compute transverse velocity
    v_t = (mu_rad * distance * parsec / 1000) / yrtosec  # [km/s]
    return v_t

#end=============================
##
###