#!/usr/bin/ParselTongue

"""
Ashish Kalyan
The "vlbireduce.py" is the principal compiler of the pipeline and has separate utility for a 
		a source calibration. Here is details of all the routes

        1) FFC, calibrates the fringe finder calibrator, and running it over all epoch you can
		   make a global model (by concatenating all epochs' calibrated data).
		2) PRC, calibrates the phase reference calibrator, and running over all epochs you can make a 
		   global model of PRC
		3) IPR (Inverse Phase referencing), calibrates the target source by utlizing the solutions
		   of calibrators (obtained with first two steps)
		4) SELFCAL_IBC1/2/3, self-calibrate the inbeam calibrator (IBC1/2/3).
        5) IBC1TOIBC2/3, apply IBC1 solutions to the IBC2/3
		6) IPR2IBC1/2/3, apply the inverse phase referencing solutions to the IBC12/3,
		
		NOTE: First, make a global model for all the calibrators, using the required methods mentioned 
			  above, and then do the target calibration.
"""

import AIPS
from AIPSData import AIPSUVData
import os, yaml, argparse, sys
from argparse import RawTextHelpFormatter
import calibration_methods as cm
import vlbi_calibration_tasks as vct


parser = argparse.ArgumentParser(description='Calibrate VLBI data', formatter_class=RawTextHelpFormatter)
parser.add_argument('source_cal', metavar='', help='Source Calibration Method (mandatory),\n' 
                            'available options are,\n'
                            'FFC:               Fringe Finder Calibrator usage,\n'
                            'PRC:               Phase Reference Calibrator usage,\n'
                            'IBC1TOIBC2/3:      Apply IBC1 solutions to IBC2/3\n'
                            'SELFCAL_IBC1/2/3:  IBC1/2/3 Self-Calibration,\n'
                            'IPR:               Inverse Phase Referencing,\n'
                            'IPR2IBC1/2/3:      Apply Inverse Phase Referencing Solutions to IBC1/2/3,\n'
                            'IBPR:              In-Beam Phase Referencing,\n'
                            'IBPR2TARGET:       Apply inbeam calibrator solutions to target,\n'
                            )
args   = parser.parse_args()

yamldir    = os.path.dirname(os.path.abspath(__file__))
expconfig  = yaml.load(open(yamldir + '/expconfig.yaml'), yaml.SafeLoader)
epochs     = expconfig['epochs']
obscode    = expconfig['obscode']
projectdir = expconfig['projectdir']
band = str(expconfig['band'])
bif  = expconfig['bif']
eif  = expconfig['eif']
refant  = expconfig['refant']
AIPS.userno = expconfig['userno']
original_stdout = sys.stdout

# contains the pulsar proper motion information, used to build pulsar preliminary model
psrpi_file = yamldir + '/psrpm_info.txt'


def build_epoch_variables(expconfig, obscode, epoch, idifilename):
    """
    Prepare variables for pipeline
    output: tabledir: idifile aips tables directory, which is a subdirectory od epochtabledir
    Return: uvdata, epochdir, vexfile, tabledir, supportfiledir
    """
    catname  = idifilename[5:12]     # aips catalogue name A1IBC10 
    epochdir = projectdir + epoch
    datadir  = epochdir + '/' + epoch + '_data'
    supportfiledir = projectdir + 'support_files'
    if not os.path.exists(supportfiledir):
        print("VLBI PIPE: Supportfiledir not found, I'm aborting the pipeline ...\n"+
                        "           use 'support_file_download.py', to create it" )
        sys.exit()
    vexfile = supportfiledir + '/vexfiles/' + obscode + epoch + ".vex"        
    bif = expconfig['bif']
    eif = expconfig['eif']

    # load data to AIPS
    uvdata = AIPSUVData(catname, 'uvdata', 1,1)
    if not uvdata.exists():
        vct.fitld_vlba(datadir + '/' + idifilename, uvdata)

    # sntables will be stored in the ibc1 folder
    tabledir = epochdir + "/sntables"
    if not os.path.exists(tabledir):
        os.mkdir(tabledir)
    
    if expconfig['splitifs']:
        if bif == 0:
            bif = 1
        tabledir = tabledir + f"/{epoch}{bif}{eif}_tables"
        if not os.path.exists(tabledir):
            os.mkdir(tabledir)
        uvsplaton = uvdata.name[:-1] + str(bif) + str(eif)
        uvsplatdata = AIPSUVData(uvsplaton, 'splat', 1,1)
        if not uvsplatdata.exists():
            doband = 0
            vct.splat_uvdata(uvdata, uvsplatdata, doband, bif, eif, solint=0)
        uvdata = uvsplatdata
    return uvdata, epochdir, vexfile, tabledir, supportfiledir


if __name__ == "__main__":
    def initiate_srccal(epoch, idikey, logdir):
        """
        Return: tuple: srccal, supportfiledir, tabledir
        """
        idifilename = obscode + epoch + idikey + band + '.idifits'

        uvdata, epochdir, vexfile, tabledir, supportfiledir = build_epoch_variables(expconfig, obscode, epoch, idifilename)
        
        srccal = cm.srccal(uvdata, expconfig, epochdir, obscode, vexfile, tabledir, logdir)
        return srccal, supportfiledir, tabledir


    # FRING (with non zero rate) and BPASS on FFC
    # Idea is that first run fring with non zeroing rate and get good solution and then apply 
    # them on data and run bpass to get the bandpass solutions
    # and later delete SN and CL tables obtained with fring on FFC
    if args.source_cal.upper() =='FFC': 
        for epoch in epochs:
            idikey = 'ibc1'
            logfile, logdir = cm.open_logger(epoch, idikey, expconfig)
            
            srccal,supportfiledir,tabledir = initiate_srccal(epoch, idikey, logdir)
            gainuse, snver = srccal.ffc_cal_non_zeroing_rate(supportfiledir)

            logfile.close()
            sys.stdout = original_stdout


    # FRING (with zeroing rate) on FFC, and FRING and CALIB on PRC
    elif args.source_cal.upper() =='PRC':
        for epoch in epochs:
            idikey = 'ibc1'
            logfile, logdir = cm.open_logger(epoch, idikey, expconfig)
            
            srccal, supportfiledir, tabledir = initiate_srccal(epoch, idikey, logdir)
            srccal.ffc_cal_zeroing_rate()
            gainuse, snver = srccal.prc_cal()

            logfile.close()
            sys.stdout = original_stdout

    
    # Apply IBC1 solutions (FFC and PRC solutions) to the IBC2/3
    elif args.source_cal.upper() in ('IBC1TOIBC2', 'IBC1TOIBC3'):
        ibcs = expconfig['ibcs']
        for epoch in epochs:
            if 'IBC2' in args.source_cal.upper():
                idikey = 'ibc2'
            else: 
                idikey = 'ibc3'
            logfile, logdir = cm.open_logger(epoch, idikey, expconfig)
            srccal, supportfiledir, tabledir = initiate_srccal(epoch, idikey, logdir)
            # gainuse, snver = srccal.ibc1_solutions2IBC2_3(supportfiledir, ibcs[int(idikey[-1])-1])
            gainuse, snver = srccal.apply_ibc1_solutions2others(supportfiledir, ibcs[int(idikey[-1])-1])
            logfile.close()
            sys.stdout = original_stdout

            
    # Self-calibration on IBC
    elif "SELFCAL_IBC" in args.source_cal.upper():
        idikey = args.source_cal.lower().split('_')[-1]
        for epoch in epochs:
            logfile, logdir = cm.open_logger(epoch, idikey, expconfig)
            srccal, supportfiledir, tabledir = initiate_srccal(epoch, idikey, logdir)
            ibcmodel = cm.load_model2aips(expconfig[idikey+'modelfile'], expconfig[idikey+'modelon'])
            gainuse, snver = srccal.ibc_selfcal(idikey, ibcmodel)
            logfile.close()
            sys.stdout = original_stdout


    # Apply IBC1 solutions (FFC and PRC solutions) to the pulsar and run CALIB on pulsar
    elif args.source_cal.upper() =='IPR':
        for epoch in epochs:
            idikey = 'psrg'
            srckey = 'target'
            source = expconfig[srckey]
            logfile, logdir = cm.open_logger(epoch, idikey, expconfig)

            srccal, supportfiledir, tabledir = initiate_srccal(epoch, idikey, logdir)
            vexfile = srccal.vexfile
            gainuse, snver = srccal.check_gainuse_snver()
            if snver == 0:
                # gainuse, snver = srccal.target_or_ibc2_cal(supportfiledir, expconfig['target'])
                gainuse, snver = srccal.apply_ibc1_solutions2others(supportfiledir, expconfig['target'])
            else:
                print("VLBI PIPE: I am assuming that the initial calibration is done, " +
                    "directly moving to the pulsar self-cal")
            srccal.pulsar_selfcal_advanceed(epoch, vexfile, psrpi_file)
            
            logfile.close()
            sys.stdout = original_stdout


    # Apply IPR solutions to IBC1/2/3
    elif args.source_cal.upper() in ('IPR2IBC1', 'IPR2IBC2', 'IPR2IBC3'):
        ibcs = expconfig['ibcs']
        for epoch in epochs:
            ibcnum = args.source_cal[-1]
            idikey = 'ibc' + ibcnum
            source = ibcs[int(ibcnum) - 1]
            
            logfile, logdir = cm.open_logger(epoch, idikey, expconfig)
            srccal, supportfiledir, tabledir = initiate_srccal(epoch, idikey, logdir)
            gainuse, snver = srccal.target_utils()
            logfile.close()
            sys.stdout = original_stdout


    # Apply IBPR solutions (self-calibration solutions on IBC) to pulsar
    elif args.source_cal.upper() == 'IBPR2TARGET':
        idikey = 'psrg'
        ibcs   = expconfig['ibcs']
        prIBCkey = 'ibc' + str(ibcs.index(expconfig['primaryinbeam']) + 1)
        for epoch in epochs:
            logfile, logdir = cm.open_logger(epoch, idikey, expconfig)
            srccal, supportfiledir, tabledir = initiate_srccal(epoch, idikey, logdir)
            srccal.target_cal_inbeam_referencing(supportfiledir, prIBCkey)
            logfile.close()
            sys.stdout = original_stdout
#end================================================
##
###          