#!/usr/bin/python3

"""
Ashish Kalyan

It downloads the ionex (igsg only), EOPS, .vex and .sum files for a given observation.
It takes inputs (epochs, obscode, projectdir, band) from the expconfig file 
so do the needful changes over there.
"""

import os, yaml, glob, requests, datetime, wget
from astropy.io import fits
import vlbi_calibration_support as vcs


def ionexDownload(vexfile, outputdir, group_name='igs'):
    """
    Download the ionex file of the observation
    Input: vexfile: vexfile for which ionex file is needed
           outputdir: ionex will be placed here
    Return: ionoex file list
    """
    sdoy, syear, edoy, eyear = vcs.get_obs_start_end_doy_year(vexfile)
    print(f"VLBI PIPE: Start and end day of year of the observation {sdoy, edoy}")
    
    ionopath1 = outputdir + f"/{group_name}g{sdoy:03}0.{str(syear)[2:4]}i"
    ionopath2 = outputdir + f"/{group_name}g{edoy:03}0.{str(eyear)[2:4]}i"
    # print(ionopath1, ionopath2)
    ionolist = []
    ionolist.append(ionopath1.split('/')[-1])
    if not os.path.exists(ionopath1):
        iono_url = "curl -u anonymous:daip@nrao.edu --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/gps/" +\
                    f"products/ionex/{syear}/{sdoy}/{group_name}g{sdoy}0.{str(syear)[2:4]}.Z " +\
                        f"> {ionopath1}.Z"
        os.system(iono_url)
        os.chdir(outputdir)
        os.system('uncompress *.Z')
        print(f"VLBI PIPE: '{ionopath1.split('/')[-1]}' file has been downloaded")

    if sdoy != edoy:
        ionolist.append(ionopath2.split('/')[-1])
        if not os.path.exists(ionopath2):
            iono_url = "curl -u anonymous:daip@nrao.edu --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/gps/" +\
                        f"products/ionex/{eyear}/{edoy}/{group_name}g{edoy}0.{str(eyear)[2:4]}.Z " +\
                            f"> {ionopath2}.Z"
            os.system(iono_url)
            os.chdir(outputdir)
            os.system('uncompress *.Z')
            print(f"VLBI PIPE: '{ionopath2.split('/')[-1]}' file has been downloaded")
    return ionolist


def eopsdownload(outputdir):
    """
    Download latest EOPS file
    """
    eopspath = outputdir + '/usno_finals.erp' 
    if not os.path.exists(eopspath):
        eop_link = "curl -u anonymous:daip@nrao.edu --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/" +\
                    f"vlbi/gsfc/ancillary/solve_apriori/usno_finals.erp > {eopspath}"
        os.system(eop_link) 
        print(f"VLBI PIPE: '{eopspath.split('/')[-1]}' file has been downloaded")
    else:
        print(f"VLBI PIPE: '{eopspath.split('/')[-1]}' file already available")
    return eopspath


def check_link_exists(url):
    try:
        response = requests.head(url)
        # If the response status code is 200, return True
        return response.status_code == 200
    except requests.RequestException:
        return False
    

def vex_sum_download(idifile, obscode, epoch, outputdir, suffix='.vex'):
    """
    Download .vex or .sum file of the observation
    """
    vexpath = outputdir + "/" + obscode.lower() + epoch.lower() + suffix
    if not os.path.exists(vexpath):
        idi = fits.open(idifile)
        obsdate = idi[0].header["DATE-OBS"]
        datesplit = obsdate.split('-')
        pydate = datetime.datetime(int(datesplit[0]), int(datesplit[1]), int(datesplit[2]))
        mm = pydate.strftime("%B")[:3]
        yy = str(pydate.year)[2:]
        
        vexurl = f"http://www.vlba.nrao.edu/astro/VOBS/astronomy/{mm.lower()}{yy}/{obscode}{epoch}/"+\
                                                f"MARK5C/{obscode}{epoch}{suffix}"

        if not check_link_exists(vexurl):
            vexurl = f"http://www.vlba.nrao.edu/astro/VOBS/astronomy/{mm.lower()}{yy}/{obscode}{epoch}/"+\
                                        f"{obscode}{epoch}{suffix}"
        print(vexurl)
        
        wget.download(vexurl, out=outputdir)
        print(f"\nVLBI PIPE: '{vexpath.split('/')[-1]}' file has been downloaded")
    else:
        print(f"VLBI PIPE: '{vexpath.split('/')[-1]}' file already available")
    return vexpath


if __name__ == "__main__":
    yamldir    = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "/examples"
    expconfig  = yaml.load(open(yamldir + '/expconfig.yaml'), yaml.SafeLoader)
    epochs     = expconfig['epochs']
    obscode    = expconfig['obscode']
    idikey     = 'ibc1'
    projectdir = expconfig["projectdir"]
    band       = expconfig['band']

    for epoch in epochs:
        idiname = obscode + epoch + idikey + str(band) + '.idifits'
        idipath = projectdir + f"{epoch}/{epoch}_data/{idiname}"
        support_filesdir = projectdir + 'support_files'
        if not os.path.exists(support_filesdir):
            os.mkdir(support_filesdir)
        vexdir = support_filesdir + '/vexfiles'
        if not os.path.exists(vexdir):
            os.mkdir(vexdir)
        sumdir = support_filesdir + '/sumfiles'
        if not os.path.exists(sumdir):
            os.mkdir(sumdir)
        ionodir = support_filesdir + '/ionofiles'
        if not os.path.exists(ionodir):
            os.mkdir(ionodir)

        vexfile = vex_sum_download(idipath, obscode.lower(), epoch.lower(), vexdir, '.vex')
        sumfile = vex_sum_download(idipath, obscode.lower(), epoch.lower(), sumdir, '.sum')
        ionexDownload(vexfile, ionodir, group_name='igs')
        if glob.glob(support_filesdir + "/*.erp") == []:
            eopsdownload(support_filesdir)

#end===================================
##
###
