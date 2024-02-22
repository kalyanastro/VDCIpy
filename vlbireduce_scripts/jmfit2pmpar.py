#!/usr/bin/ParselTongue
"""
Ashish Kalyan
This script writes an input file for PMPAR software
"""

import os, yaml, sys, re, argparse, shutil
from argparse import RawTextHelpFormatter
import AIPS
from AIPSData import AIPSImage
import vlbi_calibration_support as vcs
import vlbi_calibration_tasks as vct
import astro_utils as ast


class jmfit2pmpar:
    def __init__(self, obscode, epoch, projectdir, idikey, suffix, target, dm):
        epochdir  = projectdir + epoch
        modeldir  = epochdir + f"/source_models/{epoch}{idikey}{suffix}_{target}_model"
        resultdir = projectdir + 'results/' + obscode + '_' + suffix
        pmparfile = resultdir + '/' + target + '.' + suffix + '.pmpar.in.preliminary'
        
        self.pmparfile  = pmparfile
        self.projectdir = projectdir
        self.epochdir   = epochdir
        self.modeldir   = modeldir
        self.resultdir  = resultdir
        self.idikey  = idikey
        self.suffix  = suffix
        self.target  = target
        self.obscode = obscode
        self.dm = dm


    def jmfitting(self, epoch):
        """
        Load ...difmap.fits file to AIPS and determine the target position by JMFITting
        In case of inverse referencing, it will take the modeldivided data 
        (if model_divided_fits_for_pmpar: True)
        Return: ra, raerr, dec, decerr
        """
        modeldir = self.modeldir
        idikey   = self.idikey
        suffix   = self.suffix
        files    = os.listdir(modeldir)
        modelfilefound = False
        for file in files:
            if '_difmap.fits' in file:
                if idikey == 'psrg':
                    modelfile = modeldir + '/' + file
                    modelfilefound = True
                    break
                elif 'ibc' in idikey:
                    if expconfig['model_divided_fits_for_pmpar']:
                        if 'modeldiv_difmap.fits' in file:
                            modelfile = modeldir + '/' + file
                            modelfilefound = True
                            break
                    else:
                        if 'modeldiv_difmap.fits' not in file:
                            modelfile = modeldir + '/' + file
                            modelfilefound = True
                            break
        if not modelfilefound:
            print(f"{epoch} model file does not found")
            sys.exit()

        inname = epoch + idikey + suffix
        srcmodel = AIPSImage(inname, 'icln', 1, 1)
        jmfitfile = modelfile.replace('.fits', '.jmfit')
        if not os.path.exists(jmfitfile):
            vct.imagejmfit(modelfile, modeldir, srcmodel, pixwindow=20)
            vct.delAIPSCat(srcmodel)  
        ra, raerr, dec, decerr = vcs.jmfitout_reader(jmfitfile)
        return ra, raerr, dec, decerr
    

    def pmpar_fixed_section(self, refepoch):
        """
        Write the header part of the PMPAR input file
        """
        pmparfile = self.pmparfile
        target = self.target
        dm = self.dm
        pmparin = open(pmparfile, 'a')
        
        pmparin.write(f'# This is a pmpar input file for the {target}.\n')
        if dm != None:
            pmparin.write(f'dm    = {dm}\n')
        pmparin.write(f'epoch = {refepoch}.0\n')
        print("Reference Epoch: ", refepoch)
        pmparin.write('# Source positions at different epochs [date, ra, ra_error, dec, dec_error,'+
                       ' #epoch #mjd]\n')
        pmparin.close()
        

    def pmpar_data_section(self, epoch):
        """
        Write the data section of the PMPAR input file
        """
        projectdir = self.projectdir
        pmparfile = self.pmparfile
        obscode = self.obscode
        pmparin = open(pmparfile, 'a')
        vexfile = projectdir + '/support_files/vexfiles/' + obscode + epoch + '.vex'
        if os.path.exists(vexfile):
            pass
        else:
            print(f"{epoch} vex file does not found")
        obsepoch, decimalyear, mjd = vcs.obsyear(vexfile)
        ra, raerr, dec, decerr = self.jmfitting(epoch)
        pmparin.write(f"{decimalyear:.5f}  {ra}  {float(raerr):.8f}  {dec}  {float(decerr):.8f}"+ 
                       f"  #{epoch}  #{mjd}\n")
        pmparin.close()


# write the pulsar pmpar.in file utilizing the IBC pmpar.in file (in the case of inverse referencing)
class ibc2target_pmpar(jmfit2pmpar):
    """
    convert IBC pmpar file to target pmpar file, if inverse referencing is true
    """
    def __init__(self, obscode, epochs, projectdir, idikey, suffix, target, dm, prIBC):
        super().__init__(obscode, epochs[0], projectdir, idikey, suffix, target,dm)
        self.epochs = epochs
        self.prIBC = prIBC


    def write_target_motion_shift_to_textfile(self, vexfile0, vexfile1, psrpi_file):
        """
        It writes the target shift due to proper motion into a text file named 'target_motion_shift.txt'
        and save in the resultdir
        'psrpi_file' have the pulsar proper motion information
        """
        epochs = self.epochs
        target = self.prIBC
        vexname = vexfile1.split('/')[-1]
        motionfile = open(self.resultdir + '/target_motion_shift.txt', 'w')

        motionfile.write('epoch  rashift (deg)  decshift (deg)\n')
        for i in range(0, len(epochs)):
            epochvexfile = os.path.dirname(vexfile1) + '/' + vexname.replace(vexname[5:7], epochs[i])
            rashift, decshift = vcs.position_offset_due2_proper_motion(vexfile0, epochvexfile, psrpi_file, target)
            motionfile.write(f"{epochs[i]}  {rashift}  {decshift}\n")
        motionfile.close()


    def target_pmpar_from_ibc(self, raepo_ibc, decepo_ibc, raepo_psr, decepo_psr):
        """
        It takes ibc pmpar file as an input and convert to the target pmpar file 
        First, it subtracts ibc coordinates (infered from the global model of IBC) from the 
        pmpar position list and the resultant difference is further subtracted from the
        target position (pulsar model position plus the proper motion)
        Return: mjdlist, ra_psr, dec_psr
        """
        epochs = self.epochs
        ibcpmpar = self.pmparfile
        resultdir = self.resultdir

        # read pulsar motion shift file and append ra and dec shift to list
        f = open(resultdir + '/target_motion_shift.txt', 'r')
        lines = f.readlines()
        lines = [line for line in lines if line.strip() and not line.strip().startswith('#')]
        f.close()
        ra_motion_shift = []
        dec_motion_shift = []
        for i in range(1, len(lines)):
            linesplit = re.split(' +', lines[i].strip())
            ra_motion_shift.append(linesplit[1])
            dec_motion_shift.append(linesplit[2])
            if i == 1:
                print(f"\nVLBI PIPE: Shift due to proper motion for epoch: (rashift, decshift [deg])")
            print(f"{linesplit[0]}: ({linesplit[1]}, {linesplit[2]})")

        # read pmpar file and get ra and dec
        decimalyear, ras, raerrs, decs, decerrs, epochs, mjd = vcs.pmparin_reader(ibcpmpar)

        # subtract ibc coordinate from ra, dec list and subtract the resultant from the psr coordinate
        ra_psr, dec_psr = [],[]
        for i in range(len(ras)):
            ra, dec = ast.radec2deg(ras[i], decs[i]) #deg
            radiff_ibc = ra - raepo_ibc
            decdiff_ibc = dec - decepo_ibc
            if i==0:
                print("\nVLBI PIPE: Measured IBC position offset from the IBC fixed position for epoch:"+
                                " (radiff, decdiff [deg])")
            print(f"{epochs[i]}: ({radiff_ibc}, {decdiff_ibc})")

            ra_psr.append(raepo_psr + float(ra_motion_shift[i]) - radiff_ibc)
            dec_psr.append(decepo_psr + float(dec_motion_shift[i]) - decdiff_ibc)
        return mjd, ra_psr, dec_psr
    

def jmfit_stats_reader(jmfitstatsfile):
    """
    It reads the stats file and get the ra and dec of the image center
    """
    statsfile = open(jmfitstatsfile, 'r')
    lines = statsfile.readlines()
    statsfile.close()
    for line in lines:
        if 'Centre RA' in line:
            rastring = line.split(' ')[-1].strip()
        if 'Centre Dec' in line:
            decstring = line.split(' ')[-1].strip()
    radeg, decdeg = ast.radec2deg(rastring, decstring)
    return radeg, decdeg
    

def position_of_psr_ibc_model(ibckey):
    """
    It calls the jmfit_stats_reader and get the coordinates of the pulsar and ibc
    by supplying the pulsar model and global model of IBC
    Return: psrra, psrdec, ibcra, ibcdec [deg]
    """
    # pulsar.jmfit.stats
    pulsarmodel = expconfig['targetmodelfile']
    psrjmfitstatsfile = pulsarmodel.split('/')[-1].replace(".fits", ".jmfit.stats")
    if not os.path.exists(psrjmfitstatsfile):
        os.system(f"jmfitfromfile.py {pulsarmodel} {psrjmfitstatsfile}")

    psrjmfitstatsfile = os.path.dirname(pulsarmodel) + '/' + pulsarmodel.split('/')[-1].replace(".fits", ".jmfit.stats")
    psrra, psrdec = jmfit_stats_reader(psrjmfitstatsfile)
    
    # ibc.jmfit.stats
    ibcmodel = expconfig[ibckey + 'modelfile']
    ibcjmfitstatsfile = ibcmodel.split('/')[-1].replace(".fits", ".jmfit.stats")
    if not os.path.exists(ibcjmfitstatsfile):
        os.system(f"jmfitfromfile.py {ibcmodel} {ibcjmfitstatsfile}")

    ibcjmfitstatsfile = os.path.dirname(ibcmodel) + '/' + ibcmodel.split('/')[-1].replace(".fits", ".jmfit.stats")
    ibcra, ibcdec = jmfit_stats_reader(ibcjmfitstatsfile)

    return psrra, psrdec, ibcra, ibcdec


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make a PMPAR input file', 
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument('-re', '--refepoch', metavar='', help='Reference epoch')
    parser.add_argument('-i',  '--ibckey',    metavar='', help='ibckey, ibc1/2/3')
    args = parser.parse_args()

    yamldir   = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "/examples"
    expconfig = yaml.load(open(yamldir + '/expconfig.yaml'), yaml.SafeLoader)
    obscode   = expconfig['obscode']
    epochs    = expconfig['epochs']
    projectdir = expconfig['projectdir']
    band       = str(expconfig['band'])
    AIPS.userno = expconfig['userno']
    bif = str(expconfig['bif'])
    eif = str(expconfig['eif'])
    suffix = band
    if expconfig['splitifs']:
        suffix =  bif + eif
    try: 
        dm = expconfig['dm']
    except KeyError: 
        dm = None
    original_stdout = sys.stdout

    refepoch = args.refepoch
    target  = expconfig['target']
    idikey = 'psrg'

    resultdir = projectdir + f"results/{obscode}_{suffix}"
    if os.path.exists(resultdir):
        if vcs.yesno(f"Delete existing /results/{obscode}_{suffix} directory ? (No will abort pipeline)"):
            shutil.rmtree(resultdir)
        else: sys.exit()

    if not os.path.exists(resultdir):
        if not os.path.exists(projectdir + 'results'):
            os.mkdir(projectdir + 'results')
        os.mkdir(resultdir)

    # write a log file in result directory
    logfilepath = resultdir + f"/jmfit2pmpar.{target}.log"
    logfile = open(logfilepath, 'w')
    sys.stdout = ast.Logger(logfile)

    print(f"Inverse Referencing: {expconfig['inverse_referencing']}")
    if expconfig['inverse_referencing']:
        ibckey = args.ibckey
        ibcs   = expconfig['ibcs']
        prIBC  = ibcs[int(ibckey[-1])-1]
        target, prIBC = prIBC, target
        idikey = ibckey
        
        print("VLBI PIPE: IBC and PSR coordinates taken from corresponding stats file (RA, DEC (deg))")
        psrra, psrdec, ibcra, ibcdec = position_of_psr_ibc_model(ibckey)
        print(f"PSR: {(psrra, psrdec)}")
        print(f"IBC: {(ibcra, ibcdec)}")
    
    for index, epoch in enumerate(epochs):
        pmpar_inp = jmfit2pmpar(obscode, epoch, projectdir, idikey, suffix, target, dm)
        if index == 0:
            pmpar_inp.pmpar_fixed_section(refepoch)
        pmpar_inp.pmpar_data_section(epoch)

    if expconfig['inverse_referencing']:
        psrpi_file = yamldir + '/psrpm_info.txt'
        vexfile0 = expconfig['modelepochvexfile']
        vexfile1 = projectdir + 'support_files/vexfiles/' + obscode + epochs[0] + '.vex'
        print("\nvexfile1 will be updated during the run and will be iterate over epochs")
        print("vexfile0: ", vexfile0)
        print("vexfile1: ", vexfile1)

        ibc2psr = ibc2target_pmpar(obscode, epochs, projectdir, idikey, suffix, target, dm, prIBC)
        target_file = ibc2psr.pmparfile
        dest_file = target_file.replace(target, prIBC)
        
        ibc2psr.write_target_motion_shift_to_textfile(vexfile0, vexfile1, psrpi_file)

        mjdlist, ras_psr, decs_psr = ibc2psr.target_pmpar_from_ibc(ibcra, ibcdec, psrra, psrdec)
        
        vcs.pmparin_copy_updated_ra_dec_columns(target_file, dest_file, ras_psr, decs_psr)
        vcs.replace_word_in_file(dest_file, target, prIBC)
        
    print("jmfit2pmpar.py: Appears to have ended successfully")
    logfile.close()
    sys.stdout = original_stdout
#end=================================================================
##
###