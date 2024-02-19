# Ashish Kalyan


import os, sys, shutil
from AIPSData import AIPSUVData, AIPSImage
import vlbi_calibration_tasks as vct
import vlbi_calibration_support as vcs
import astro_utils as ast


vlbaant = {1:"BR", 2:"FD", 3:"HN", 4:"KP", 5:"LA", 6:"MK", 7:"NL", 8:"OV", 9:"PT", 10:"SC"}
    

def open_logger(epoch, idikey, expconfig):
    """
    Write logs
    Return: logfile, logdir
    """
    projectdir = expconfig['projectdir']
    band = str(expconfig['band'])
    bif  = expconfig['bif']
    eif  = expconfig['eif']
    logdir = projectdir + epoch + '/aips_logs'
    if not os.path.exists(logdir):
        os.mkdir(logdir)

    logfilepath = logdir + '/' + epoch + idikey + band + ".log"
    if expconfig['splitifs']:
        logfilepath = logdir + '/' + epoch + idikey + str(bif) + str(eif) + ".log"
    try:
        logfilemode = expconfig['logfilemode']
    except KeyError:
        logfilemode = 'a'
        
    logfile = open(logfilepath, logfilemode)
    sys.stdout = ast.Logger(logfile)
    return logfile, logdir


def aips2disk(uvdata, source, epochdir, catclass='split'):
    """
    split and store AIPSCat to disk
    """
    splitdata = AIPSUVData(source, catclass, 1, 1)
    vct.split_uvdata(uvdata, [source], 0, catclass, 1)
    modelfolder = epochdir + '/source_models/'
    srcmodelfolder = modelfolder + uvdata.name + '_' + source + '_model'
    fitsdataout = srcmodelfolder + '/' + uvdata.name + '_' + source + '_aips.fits'
    if os.path.exists(srcmodelfolder):
        print('VLBI PIPE: AIPS2DISK: %s folder already exists' %srcmodelfolder.split('/')[-1])
        if vcs.yesno("Delete existing UV dataset folder ? (No will abort pipeline)"):
            shutil.rmtree(srcmodelfolder)
    if not os.path.exists(fitsdataout):
        if not os.path.exists(srcmodelfolder):
            if not os.path.exists(modelfolder):       
                os.mkdir(modelfolder)
            os.mkdir(srcmodelfolder)
        vct.fittp(splitdata, fitsdataout)
        vct.delAIPSCat(splitdata)
    

def snedt(uvdata, snver):
    """
    Edit SN table
    Return: snver (one higher to the input)
    """
    print(f"VLBI PIPE: Do a manual eidting of SN{snver}, if required")
    user_ip = str(input('Do you want to run the SNEDT task (y/n)?:'))
    if user_ip.upper()=='Y':
        vct.snedt(uvdata, snver)
        snver += 1
    return snver


def load_model2aips(modelfile, modelon):
    """
    Load calibrator model to AIPSCat
    Input: modelfile: Clean Componenets fits file
           modelon: AIPSCat name (outname)
    Return: modeldata
    """
    if modelon not in ('', None):
        modeldata = AIPSImage(modelon, 'icln',1,1)
        if not modeldata.exists():
            vct.fitld_image(modelfile, modeldata)
            if not modeldata.exists():
                sys.exit()
    else:
        if vcs.yesno("model file doesn't supplied in config file, do you want to "+
                     "procced without it? (No will abort pipeline)"):
            modeldata = AIPSImage('', '',1,1)
        else: sys.exit()
    return modeldata


class calMethods:
    # class variables
    def __init__(self, uvdata, expconfig, tabledir, obscode):
        # instance variables
        self.uvdata    = uvdata
        self.tabledir  = tabledir
        self.expconfig = expconfig
        self.obscode   = obscode


    def check_gainuse_snver(self):
        """
        Determines the SN and CL version
        Return: gainuse, snver
        """
        uvdata  = self.uvdata
        gainuse = 0
        snver   = 0
        for table in uvdata.tables:
            if table[1] == 'AIPS SN':
                snver += 1
            elif table[1] == 'AIPS CL':
                gainuse += 1
        return gainuse, snver


    def data_inspect(self):
        """
        Get dignostic plots, lister info, and shift the source (if shift is defined in expconfig file)
		Return: source list
        """
        expconfig = self.expconfig
        tabledir  = self.tabledir
        uvdata    = self.uvdata
        listrfile = tabledir + '/' + uvdata.name + '_listr.txt'
        sufile    = tabledir + '/' + uvdata.name + '_su_tbout.txt'
        prtanplot = tabledir + '/' + uvdata.name + '_prtan.ps'
        
        if not os.path.exists(listrfile):
            vct.listr(uvdata, listrfile, optype='SCAN')
        if not os.path.exists(sufile):
            vct.writetable(uvdata, 'SU', 1, sufile)
        srclist = vcs.su_src_reader(sufile)

        # source shifting [optional]
        try:
            srcshift = expconfig['shift']
        except KeyError:
            srcshift = None
        if srcshift != None:
            source = expconfig['shift'][0]
            if source in srclist:
                rashift = expconfig['shift'][1]
                decshift = expconfig['shift'][2]
                print(f"Shifting {source} by {rashift} in RA [mas] and {decshift} in Decl. [mas]")
                vct.shift_source(uvdata, source, rashift, decshift, 0)
                sufile_updated = tabledir + '/' + uvdata.name + '_su_source_shifted_tbout.txt'
                vct.writetable(uvdata, 'SU', 1, sufile_updated)

        if 'ibc1' in uvdata.name:
            if not os.path.exists(prtanplot):
                vct.plotan(uvdata, prtanplot)
        self.sufile = sufile
        return srclist
    

    def flagging(self, uvdata, epoch):
        """
        It does the fringe rate based flagging and user defined flagging.
        """
        expconfig = self.expconfig
        projectdir = expconfig['projectdir']
        obscode = expconfig['obscode']
        vct.fringerateflag(uvdata, suppressionfactor=7, flagver=1)
        flagfile = projectdir + f"flags/{obscode}{epoch}.flag"
        if os.path.exists(flagfile):
            vct.userflag(uvdata, 1, flagfile)


    def reduce_without_calibrator(self, supportfiledir, vexfile):
        """
        It does the fringe rate based flagging and then run TECOR, CLCOR (EOPS and PANG)
		Return: gainuse
        """
        uvdata  = self.uvdata
        expconfig = self.expconfig
        gainuse = 1
        epoch = uvdata.name[0:2]

        self.flagging(uvdata, epoch)
        vct.tecor_iono(uvdata, supportfiledir + '/ionofiles', vexfile, 'igsg', expconfig['petrov'])
        gainuse += 1

        vct.eops_clcor(uvdata, supportfiledir + '/usno_finals.erp')
        gainuse += 1
        vct.pang_clcor(uvdata)
        gainuse += 1
        return gainuse


    def amplitude_calibration(self):
        """
        It runs ACCOR and APCAL and store the smoothed SN tables to disk
        Return: gainuse, snver
        """
        tabledir  = self.tabledir
        uvdata    = self.uvdata
        expconfig = self.expconfig
        refant    = expconfig['refant']
        flagver   = 1
        gainuse, snver = self.check_gainuse_snver()

        # accor
        accorsn = tabledir + "/accor.sn"
        if not os.path.exists(accorsn):
            vct.accor(uvdata, expconfig['accorsolmins'])
            snver += 1
            vct.snsmo(uvdata, refant, snver, 'BOTH', timemins=20, ampdev=0.03, phasedev=0, 
                        delaydev=0, passthrough=1, ratesmoothmin=0, ratedev=0)
            snver += 1
            vct.snplt(uvdata, tabledir + '/accor.ps', 'SN', snver, 'AMP')
            if not expconfig['skipsnedt']:
                print("VLBI PIPE: Examine ACCOR SN table (amplitude) and flag the outliers")
                snver = snedt(uvdata, snver)
            vct.writetable(uvdata, 'SN', snver, accorsn)
        else:
            snver += 1
            vct.loadtable(uvdata, accorsn, snver)
        vct.clcal_applySN(uvdata, refant, sources=[''], snver=snver, interpol='2PT', opcode='CALI')
        gainuse += 1

        # apcal
        apcalsn = tabledir + "/apcal.sn"
        if not os.path.exists(apcalsn):
            vct.apcal(uvdata)
            snver += 1
            vct.snsmo(uvdata, refant, snver, 'BOTH', timemins=20, ampdev=1.5, phasedev=0, 
                        delaydev=0, passthrough=1, ratesmoothmin=0, ratedev=0)
            snver += 1
            vct.snplt(uvdata, tabledir + '/apcal.ps', 'SN', snver, 'AMP')

            # flag bad data points and edit APCAL SN table if required
            print("VLBI PIPE: Examine APCAL SN table and if its values are >50 or zero for"+ 
                "some IFs and polirazations then manually flag those values")
            snver = snedt(uvdata, snver)

            print("VLBI PIPE: Examine APCAL SN table and if its values are >50 or zero in"+
                "all IFs and polirazations then flag the antenna")
            user_ip = str(input('Do you want to flag any antenna (UVFLG [y/n])?:'))
            if user_ip.upper()=='Y':
                antennas = []
                antennas = [int(item) for item in input("List of antennas that need to be"+ 
                                        "flagged (don't use comma (eg 1 2 3)): ").split()]
                print("Flag Ant List: ", antennas)
                vct.antflag(uvdata, antennas, flagver, sources='')

                # add antenna flag information to flag file
                flagfile = expconfig['projectdir'] + f"flags/{expconfig['obscode']}{uvdata.name[0:2]}.flag"
                if not os.path.exists(flagfile):
                    file = open(flagfile, 'w')
                    file.write(f"opcode = 'FLAG'\n")
                    file.write(f"dtimrang = 1  timeoff = 0\n")
                    file.close()
                file = open(flagfile, 'a')
                for ant in antennas:
                    file.write(f"ant_name='{vlbaant[ant]}' timerang=000,00,00,00, 900,00,00,00 bif=0 eif=0  sources='' Reason='TSYS' /\n")
                file.close()

            vct.writetable(uvdata, 'SN', snver, apcalsn)
        else:
            snver += 1
            vct.loadtable(uvdata, tabledir + '/apcal.sn', snver)
        vct.clcal_applySN(uvdata, refant, sources=[''], snver=snver, interpol='2PT', opcode='CALI')
        gainuse += 1 
        return gainuse, snver


    def primary_beamcorr(self, vexfile, source):
        """
        Primary beam correction
        NOTE: source is required only in case of pulsar/inbean2/3 calibrtation
		Return: gainuse, snver
        """
        uvdata    = self.uvdata
        expconfig = self.expconfig
        refant    = expconfig['refant']
        gainuse, snver = self.check_gainuse_snver()
        print(f"VLBI PIPE: Executing 'primary_beamcorr' method with SN{snver} and CL{gainuse}")
        
        if '174' in expconfig['obscode']:
            targetpt = expconfig['target'][0:5] + 'PT'
        else:
            targetpt = expconfig['target'] + 'PT'
        
        fieldsourcenames = {}
        if 'ibc1' in uvdata.name:
            fieldsourcenames[expconfig['ffc']] = expconfig['ffc']
            fieldsourcenames[expconfig['prc']] = expconfig['prc']
            fieldsourcenames[targetpt] = expconfig['ibcs'][0]
            print(f"VLBI PIPE: FieldSources {fieldsourcenames}")
        
            vct.correct_primarybeam(uvdata, snver, vexfile, fieldsourcenames, iscal=False, phasecentrenum=0, 
                        issearch=False, isonepointing=False, onlygettimes=False, skipmissingsources=False)
            snver += 1
            vct.clcal_applySN(uvdata, refant, sources=[''], snver=snver, interpol='2PT', opcode='CALI')
            gainuse += 1
        else:
            fieldsourcenames = {targetpt: source}
            print(f"VLBI PIPE: FieldSources {fieldsourcenames}")
            
            vct.correct_primarybeam(uvdata, snver, vexfile, fieldsourcenames, iscal=False, phasecentrenum=0, 
                        issearch=False, isonepointing=False, onlygettimes=False, skipmissingsources=False)
            snver += 1
            vct.clcal_applySN(uvdata, refant, sources=[''], snver=snver, interpol='SELN', opcode='CALI')
            gainuse += 1
        return gainuse, snver
    

    def pccor_run(self, uvdata, listrfile, srcname, ampcalscan):
        """
        Run PCCOR if PC table is present, and apply solution to uvdata
        Also write the sntable to disk
        """
        expconfig = self.expconfig
        refant    = expconfig['refant']
        tabledir  = self.tabledir
        gainuse, snver = self.check_gainuse_snver()
        pctable = False
        for table in uvdata.tables:
            if table[1] == 'AIPS PC':
                pctable = True
                break
        if pctable:
            pccorsn = tabledir + "/pccor.sn"
            if not os.path.exists(pccorsn):
                snver += 1
                vct.pccor(uvdata, listrfile, srcname, snver, ampcalscan, refant)
                vct.snplt(uvdata, tabledir + '/pccor.ps', 'SN', snver, 'PHAS')
                vct.writetable(uvdata, 'SN', snver, pccorsn)
            else:
                snver += 1
                vct.loadtable(uvdata, pccorsn, snver)
            vct.clcal_applySN(uvdata, refant, sources=[''], snver=snver, interpol='2PT', opcode='CALI')
            gainuse += 1
        else: print('VLBI PIPE: PC table is not found, skipping PCCOR task')
        return snver, gainuse
        

    def ffc_fring_bpass_non_zeroing_rate(self, ffcmodel):
        """
        It runs FRING and BPASS on FFS, keeping rate to nonzero
		Return: gainuse, snver
        """
        tabledir  = self.tabledir
        uvdata    = self.uvdata
        expconfig = self.expconfig
        refant    = expconfig['refant']
        gainuse, snver = self.check_gainuse_snver()
        print(f"VLBI PIPE: I am using 'FFC_FRING_BPASS_NON_ZEROING_RATE' method "+ 
                        f"with SN{snver} and CL{gainuse}")
        
        vct.fring_fitting(uvdata, refant, ffcmodel, [expconfig['ffc']], -1, 0, expconfig['fringffcsolmin1'], 
                          dostokesi=0, snr=6, delaywin=400, ratewin=30, inttimesecs=2)
        snver += 1
        vct.snsmo(uvdata, refant, snver, smotype='VLRI', timemins=20, ampdev=0, phasedev=0, delaydev=10, 
                                passthrough=0, ratesmoothmin=3, ratedev=10)
        snver += 1

        vct.snplt(uvdata, tabledir + '/ffc.fring0.ps', 'SN', snver, 'DELA')

        if not expconfig['skipsnedt']:
            print("VLBI PIPE: Examine FFC FRING SN table (delay) and flag the outliers")
            snver = snedt(uvdata, snver)

        vct.clcal_applySN(uvdata, refant, [expconfig['ffc']], snver, 'AMBG', opcode='CALP')
        gainuse += 1

        vct.bpass(uvdata, refant, ffcmodel, [expconfig['ffc']])
        vct.plotbp(uvdata, tabledir + '/bpass.ps')
        vct.writetable(uvdata, 'BP', 1, tabledir + "/bpass.bp")
        
        return gainuse, snver


    def ffc_fring_zeroing_rate(self, ffcmodel):
        """
        It runs FRING on FFC keeping rate to zero
        NOTE: if you have ran "ffc_fring_bpass_non_zeroing_rate" then make sure that 
        you have deleted the last CL and last two SN
		Return: gainuse, snver
        """
        tabledir   = self.tabledir
        uvdata     = self.uvdata
        expconfig  = self.expconfig
        refant     = expconfig['refant']
        gainuse, snver = self.check_gainuse_snver()
        print(f"VLBI PIPE: I am using 'FFC_FRING_ZEROING_RATE' method with SN{snver} and CL{gainuse}")

        vct.fring_fitting(uvdata, refant, ffcmodel, [expconfig['ffc']], 1, 1, expconfig['fringffcsolmin2'], 
                          dostokesi=0, snr=6, delaywin=400, ratewin=30, inttimesecs=2)
        snver += 1
        vct.snsmo(uvdata, refant, snver, smotype='DELA', timemins=20, ampdev=0, phasedev=0, delaydev=10, 
                                passthrough=0, ratesmoothmin=0, ratedev=0)
        snver += 1
        
        vct.snplt(uvdata, tabledir + '/ffc.fring.ps', 'SN', snver, 'DELA')

        if not expconfig['skipsnedt']:
            print("VLBI PIPE: Examine FFC FRING SN table (zeroing rate, delay) and flag the outliers")
            snver = snedt(uvdata, snver)

        vct.writetable(uvdata, 'SN', snver, tabledir + "/ffc.fring.sn")
        vct.clcal_applySN(uvdata, refant, sources=[''], snver=snver, interpol='2PT', opcode='CALP')
        gainuse += 1
        return gainuse, snver


    def prc_fring(self, prcmodel):
        """
        It runs FRING on the PRC
		Return: gainuse, snver
        """
        uvdata    = self.uvdata
        tabledir  = self.tabledir
        expconfig = self.expconfig
        refant    = expconfig['refant']
        gainuse, snver = self.check_gainuse_snver()
        print(f"VLBI PIPE: I am using 'PRC_FRING' method with SN{snver} and CL{gainuse}")

        vct.fring_fitting(uvdata, refant, prcmodel, [expconfig['prc']], 1, 0, solint=3, dostokesi=0, 
                          snr=5, delaywin=400,ratewin=30,inttimesecs=2)
        snver += 1
        vct.snsmo(uvdata, refant, snver, 'VLRI', timemins=20, ampdev=0, phasedev=0, delaydev=10, 
                  passthrough=0, ratesmoothmin=3, ratedev=10)
        snver += 1

        vct.snplt(uvdata, tabledir + '/prc.fring.ps', 'SN', snver, 'DELA')

        # flag bad data points and edit prc_fring SN table if required
        if not expconfig['skipsnedt']:
            print("VLBI PIPE: Examine PRC FRING SN table (delay) and flag the outliers")
            snver = snedt(uvdata, snver)

        vct.writetable(uvdata, 'SN', snver, tabledir + "/prc.fring.sn")
        vct.clcal_applySN(uvdata, refant, [''], snver, 'AMBG', 'CALI')
        gainuse += 1
        return gainuse, snver
        
        
    def prc_selfcal(self, prcmodel):
        """
        It runs phase and amplitude self-calibration on PRC.
		Return: gainuse, snver
        """
        uvdata    = self.uvdata
        tabledir  = self.tabledir
        expconfig = self.expconfig
        refant    = expconfig['refant']
        gainuse, snver = self.check_gainuse_snver()
        print(f"VLBI PIPE: I am using 'PRC_SELFCAL' method with SN{snver} and CL{gainuse}")
        
        # phase self-calibration
        prccalibdata = AIPSUVData(expconfig['prc'], 'calib', 1,1)
        vct.splittoseq(uvdata, [expconfig['prc']], gainuse, 'calib', prccalibdata.seq)
        prccalibdata.table('NX', 1).zap()

        vct.calib_selfcal(prccalibdata, prcmodel, refant, [expconfig['prc']], solint=1.5, dostokesi=0, 
                                  doamp=0, avgif=0, calibrate=0, snr=5, weightit=0, normalize=0)
        vct.snsmo(prccalibdata, refant, invers=1, smotype='AMP', timemins=300, ampdev=0.2, phasedev=0,
                   delaydev=0, passthrough=0, ratesmoothmin=0, ratedev=0)

        # Copy phase self-cal SN table to caldata and apply the solutions
        vct.tacop(prccalibdata, uvdata, 'SN', inver=2)
        snver += 1

        vct.snplt(uvdata, tabledir + '/prc.calib.p.ps', 'SN', snver, 'PHAS')

        if not expconfig['skipsnedt']:
            print("VLBI PIPE: Examine PRC PHASE CALIB SN table (phase) and flag the outliers")
            snver = snedt(uvdata, snver)

        vct.writetable(uvdata, 'SN', snver, tabledir + "/prc.calib.p.sn")
        
        vct.clcal_applySN(uvdata, refant, [''], snver, '2PT', opcode='CALI')
        gainuse += 1
        vct.delAIPSCat (prccalibdata)

        # amplitude and phase self-calibration
        vct.splittoseq(uvdata, [expconfig['prc']], gainuse, 'calib', prccalibdata.seq)
        prccalibdata.table('NX', 1).zap()

        vct.calib_selfcal(prccalibdata, prcmodel, refant, [expconfig['prc']], solint=20, dostokesi=0, 
                                  doamp=1, avgif=0, calibrate=0, snr=6, weightit=0, normalize=0)
        vct.snsmo(prccalibdata, refant, invers=1, smotype='AMP', timemins=300, ampdev=0.2, phasedev=0,
                   delaydev=0, passthrough=0, ratesmoothmin=0, ratedev=0)

        # Copy amp & phase self-cal SN table to caldata and apply the solutions
        vct.tacop(prccalibdata, uvdata, 'SN', inver=2)
        snver += 1

        vct.snplt(uvdata, tabledir + '/prc.calib.ap.ps', 'SN', snver, 'AMP')

        if not expconfig['skipsnedt']:
            print("VLBI PIPE: Examine PRC AP CALIB SN table (amplitude) and flag the outliers")
            snver = snedt(uvdata, snver)
        
        vct.writetable(uvdata, 'SN', snver, tabledir + "/prc.calib.ap.sn")
        
        vct.clcal_applySN(uvdata, refant, [''], snver, '2PT', opcode='CALI')
        gainuse += 1
        vct.delAIPSCat (prccalibdata)
        return gainuse, snver


    def ibc_selfcal(self, srckey, ibcmodel):
        """
        It runs phase and amplitude self-calibration on IBC.
		Return: gainuse, snver
        """
        uvdata    = self.uvdata
        tabledir  = self.tabledir
        expconfig = self.expconfig
        refant    = expconfig['refant']
        pribc     = expconfig['ibcs'][int(srckey[-1]) - 1]
        gainuse, snver = self.check_gainuse_snver()
        print(f"VLBI PIPE: I am using 'IBC_SELFCAL' method with SN{snver} and CL{gainuse}")
        
        splitibcdata = AIPSUVData(pribc, 'temp', 1, 1)
        vct.splittoseq(uvdata, [pribc], gainuse, 'temp', splitibcdata.seq)
        splitibcdata.table('NX', 1).zap()

        # Normalize split data by dividing model
        uvsubibcdata = AIPSUVData(pribc, 'uvsub', 1, 1)
        vct.uvsub_divide(splitibcdata, ibcmodel)
        vct.delAIPSCat(splitibcdata)
        
        # phase self-calibration by averaging Stokes and IFs
        vct.calib_selfcal(uvsubibcdata, ibcmodel, refant, [pribc], solint=2, dostokesi=1, 
                                  doamp=0, avgif=1, calibrate=0, snr=5, weightit=0, normalize=0)
        
        # phase self-calibration by averaging only Stokes 
        vct.calib_selfcal(uvsubibcdata, ibcmodel, refant, [pribc], solint=20, dostokesi=1, 
                                  doamp=0, avgif=0, calibrate=0, snr=5, weightit=0, normalize=0)
        
        # Copy phase self-cal SN tables to caldata and apply the solutions
        vct.tacop(uvsubibcdata, uvdata, 'SN', inver=1)
        snver += 1

        vct.snplt(uvdata, tabledir + f"/{srckey}.calib.p1.ps", 'SN', snver, 'PHAS')

        if not expconfig['skipsnedt']:
            print("VLBI PIPE: Examine PRC FRING SN table (phase) and flag the outliers")
            snver = snedt(uvdata, snver)
        
        vct.writetable(uvdata, 'SN', snver, tabledir + f"/{srckey}.calib.p1.sn")
        
        vct.clcal_applySN(uvdata, refant, sources=[''], snver=snver, interpol='2PT', opcode='CALP')
        gainuse += 1

        vct.tacop(uvsubibcdata, uvdata, 'SN', inver=2)
        snver += 1

        vct.snplt(uvdata, tabledir + f"/{srckey}.calib.p2.ps", 'SN', snver, 'PHAS') 

        if not expconfig['skipsnedt']:
            print("VLBI PIPE: Examine PRC FRING SN table (phase) and flag the outliers")
            snver = snedt(uvdata, snver)
        
        vct.writetable(uvdata, 'SN', snver, tabledir + f"/{srckey}.calib.p2.sn")
             
        vct.clcal_applySN(uvdata, refant, sources=[''], snver=snver, interpol='2PT', opcode='CALP')
        gainuse += 1
        vct.delAIPSCat(uvsubibcdata)

        # amplitude-and-phase self-calibration
        ibc1apcalibdata = AIPSUVData(pribc, 'calib', 1, 1)
        vct.splittoseq(uvdata, [pribc], gainuse, 'calib', ibc1apcalibdata.seq)
        ibc1apcalibdata.table('NX', 1).zap()

        vct.calib_selfcal(ibc1apcalibdata, ibcmodel, refant, [pribc], solint=10, dostokesi=1,
                        doamp=1, avgif=0, calibrate=0, snr=7, weightit=0, normalize=0) 
        vct.snsmo(ibc1apcalibdata, refant, invers=1, smotype='AMP', timemins=300, ampdev=0.5, phasedev=0,
                        delaydev=0, passthrough=0, ratesmoothmin=0, ratedev=0)
        vct.tacop(ibc1apcalibdata, uvdata, 'SN', inver=2)
        snver += 1

        vct.snplt(uvdata, tabledir + f"/{srckey}.calib.ap.ps", 'SN', snver, 'AMP')

        if not expconfig['skipsnedt']:
            print("VLBI PIPE: Examine PRC FRING SN table (amplitude) and flag the outliers")
            snver = snedt(uvdata, snver)
        
        vct.writetable(uvdata, 'SN', snver, tabledir + f"/{srckey}.calib.ap.sn")
        
        vct.clcal_applySN(uvdata, refant, sources=[''], snver=snver, interpol='2PT', opcode='CALP')
        gainuse += 1
        vct.delAIPSCat(ibc1apcalibdata)
        return gainuse, snver    


    def apply_scintillation_correction(self):
        """
        Apply the scintillation corrrection to pulsar
        Return: tuple: (gainuse, snver)
        """
        uvdata    = self.uvdata
        expconfig = self.expconfig
        tabledir  = self.tabledir
        refant    = expconfig['refant']
        target    = expconfig['target']
        gainuse, snver = self.check_gainuse_snver()

        splituvdata = AIPSUVData(target, 'scint', 1,1)
        vct.splittoseq(uvdata, [target], gainuse, 'scint', splituvdata.seq)
        splituvdata.table('NX', 1).zap()

        outputfile = tabledir + '/' + uvdata.name + '_' + target + '_ScintCorr.txt'
        vct.wizCorrectScint(uvdata, 1, snver+1, splituvdata, expconfig['scntsolmins'], outputfile)
        snver += 1
        vct.snplt(uvdata, tabledir + '/psr.scint.ps', 'SN', snver, 'AMP')
        vct.clcal_applySN(uvdata, refant, [''], snver, '2PT', opcode='CALP')
        gainuse += 1
        scntuvdata = AIPSUVData(target, 'scntt', 1,1)
        vct.delAIPSCat(splituvdata)
        vct.delAIPSCat(scntuvdata)
        return gainuse, snver 
    

    def pulsar_selfcal(self, srcmodel):
        """
        It runs phase and amplitude self-calibration with pulsar
		Return: gainuse, snver
        """
        uvdata    = self.uvdata
        expconfig = self.expconfig
        tabledir  = self.tabledir
        refant    = expconfig['refant']
        srckey    = 'target'
        source    = expconfig[srckey]
        gainuse, snver = self.check_gainuse_snver()
        print(f"VLBI PIPE: I am using 'PULSAR_SELFCAL' method with SN{snver} and CL{gainuse}")
    
        psrcalibdata = AIPSUVData(source, 'calib', 1, 1)
        vct.splittoseq(uvdata, [source], gainuse, 'calib', psrcalibdata.seq)
        psrcalibdata.table('NX', 1).zap()
        
        # phase self-calibration by averaging Stokes and IFs
        vct.calib_selfcal(psrcalibdata, srcmodel, refant, [source], solint=2, dostokesi=1, 
                                doamp=0, avgif=1, calibrate=0, snr=5, weightit=0, normalize=0)
        # Copy phase self-cal SN tables to caldata and apply the solutions
        vct.tacop(psrcalibdata, uvdata, 'SN', inver=1)
        snver += 1

        vct.snplt(uvdata, tabledir + f"/{srckey}.calib.p1.ps", 'SN', snver, 'PHAS')

        if not expconfig['skipsnedt']:
            print("VLBI PIPE: Examine PRC FRING SN table (phase) and flag the outliers")
            snver = snedt(uvdata, snver)
        
        vct.writetable(uvdata, 'SN', snver, tabledir + f"/{srckey}.calib.p1.sn")
        
        vct.clcal_applySN(uvdata, refant, sources=[''], snver=snver, interpol='2PT', opcode='CALP')
        gainuse += 1
        vct.delAIPSCat(psrcalibdata)

        vct.splittoseq(uvdata, [source], gainuse, 'calib', psrcalibdata.seq)
        psrcalibdata.table('NX', 1).zap()

        # phase self-calibration by averaging only Stokes 
        vct.calib_selfcal(psrcalibdata, srcmodel, refant, [source], solint=20, dostokesi=1, 
                                doamp=0, avgif=0, calibrate=0, snr=5, weightit=0, normalize=0)
        # Copy phase self-cal SN tables to caldata and apply the solutions
        vct.tacop(psrcalibdata, uvdata, 'SN', inver=1)
        snver += 1

        vct.snplt(uvdata, tabledir + f"/{srckey}.calib.p2.ps", 'SN', snver, 'PHAS')

        if not expconfig['skipsnedt']:
            print("VLBI PIPE: Examine PRC FRING SN table (phase) and flag the outliers")
            snver = snedt(uvdata, snver)
        
        vct.writetable(uvdata, 'SN', snver, tabledir + f"/{srckey}.calib.p2.sn") 
           
        vct.clcal_applySN(uvdata, refant, sources=[''], snver=snver, interpol='2PT', opcode='CALP')
        gainuse += 1
        vct.delAIPSCat(psrcalibdata)
        return gainuse, snver
# =================================================


class srccal(calMethods):
    def __init__(self, uvdata, expconfig, epochdir, obscode, vexfile, tabledir, logdir):
        super().__init__(uvdata, expconfig, tabledir, obscode)   # Call the constructor of vlbireduce
        self.uvdata   = uvdata
        self.vexfile  = vexfile
        self.epochdir = epochdir
        self.tabledir = tabledir
        self.logdir   = logdir
        self.obscode  = obscode
        

    def ffc_cal_non_zeroing_rate(self, supportfiledir):
        """
        It clubs the 'data_inspect', 'reduce_without_calibrator', 'amplitude_calibration', 
        'primary_beamcorr', 'ffc_nonzero_rate.' It is designed for FFC calibration and by 
        running over all epoch you can make a global model of FFC.
        Return: tuple (gainuse, snver)
        """
        uvdata    = self.uvdata
        tabledir  = self.tabledir
        vexfile   = self.vexfile
        expconfig = self.expconfig
        srclist   = self.data_inspect()
        gainuse   = self.reduce_without_calibrator(supportfiledir, vexfile)

        self.amplitude_calibration()
        gainuse, snver = self.primary_beamcorr(vexfile, srclist[0])

        listrfile = tabledir + '/' + uvdata.name + '_listr.txt'
        self.pccor_run(uvdata, listrfile, expconfig['ffc'], expconfig['ampcalscan'])

        ffcmodel = load_model2aips(expconfig['ffcmodelfile'], expconfig['ffcmodelon'])
        gainuse, snver = self.ffc_fring_bpass_non_zeroing_rate(ffcmodel)

        if ffcmodel.name == '':
            aips2disk(uvdata, expconfig['ffc'], self.epochdir, 'split')

        # if expconfig['save_calibrated_data']:
        #     aips2disk(uvdata, expconfig['ffc'], self.epochdir, 'split')

        # outputfile = tabledir + '/' + uvdata.name + '_FFC_CL' + str(gainuse) + '_bp_possm.ps'
        # vct.possm(uvdata, outputfile, 1, gainuse, 1, [expconfig['ffc']], nplot=2)
        return gainuse, snver
    

    def ffc_cal_zeroing_rate(self):
        """
        It assums that "ffc_cal_non_zeroing_rate" has been already execuated 
        so it first delete the last two SN and last CL and then runs FRING with FFC,
        keeping rate to zero
        """
        uvdata    = self.uvdata
        expconfig = self.expconfig
        gainuse, snver = self.check_gainuse_snver()

        vct.deletetable(uvdata, 'CL', gainuse)
        vct.deletetable(uvdata, 'SN', snver) 
        vct.deletetable(uvdata, 'SN', snver-1)

        ffcmodel = load_model2aips(expconfig['ffcmodelfile'], expconfig['ffcmodelon'])
        gainuse, snver = self.ffc_fring_zeroing_rate(ffcmodel)


    def prc_cal(self):
        """
        It clubs the 'prc_cal' and 'prc_selfcal' and ultimately do the fringe-fitting and 
        self-calibration (only if model exists) with PRC.
		Return: gainuse, snver
        """
        expconfig = self.expconfig
        prcmodel  = load_model2aips(expconfig['prcmodelfile'], expconfig['prcmodelon'])
        gainuse, snver = self.prc_fring(prcmodel)

        if prcmodel.name == '':
            aips2disk(self.uvdata, expconfig['prc'], self.epochdir, 'split')

        if prcmodel.exists():
            gainuse, snver = self.prc_selfcal(prcmodel)

            # if expconfig['save_calibrated_data']:
            #     aips2disk(self.uvdata, expconfig['prc'], self.epochdir, 'split')
        return gainuse, snver
    

    def ibc1_utils(self):
        """
        It does copy FRING and BPASS table obtained with FFC and FRING and 
        CALIB solution with PRC and apply to uvdata
        Additionally apply the PCCOR solution if pccor.sn table is present
        Return: gainuse, snver
        """
        uvdata    = self.uvdata
        tabledir  = self.tabledir
        expconfig = self.expconfig
        refant    = expconfig['refant']
        gainuse, snver = self.check_gainuse_snver()

        print(f"VLBI PIPE: I am using 'IBC1_UTILS' method with SN{snver} and CL{gainuse}")

        # apply pccor solution if pccor.sn table exists
        if os.path.exists(tabledir + '/pccor.sn'):
            snver += 1
            vct.loadtable(uvdata, tabledir + '/pccor.sn', snver)
            vct.clcal_applySN(uvdata, refant, [''], snver, '2PT', opcode='CALI')
            gainuse += 1
        
        snver += 1
        vct.loadtable(uvdata, tabledir + '/ffc.fring.sn', snver)
        vct.clcal_applySN(uvdata, refant, [''], snver, '2PT', opcode='CALP')
        gainuse += 1

        vct.loadtable(uvdata, tabledir + '/bpass.bp', 1)

        snver += 1
        vct.loadtable(uvdata, tabledir + '/prc.fring.sn', snver)
        vct.clcal_applySN(uvdata, refant, [''], snver, 'AMBG', 'CALI')
        gainuse += 1

        snver += 1
        vct.loadtable(uvdata, tabledir + '/prc.calib.p.sn', snver)
        vct.clcal_applySN(uvdata, refant, [''], snver, '2PT', 'CALI')
        gainuse += 1
        
        snver += 1
        vct.loadtable(uvdata, tabledir + '/prc.calib.ap.sn', snver)
        vct.clcal_applySN(uvdata, refant, [''], snver, '2PT', 'CALI')
        gainuse += 1
        return gainuse, snver


    # def ibc1_solutions2IBC2_3(self, supportfiledir, source):
    def apply_ibc1_solutions2others(self, supportfiledir, source):
        """
        Apply the IBC1 solutions to the other IBCs / pulsar
        """
        vexfile = self.vexfile
        srclist = self.data_inspect()
        
        self.reduce_without_calibrator(supportfiledir, vexfile)
        self.amplitude_calibration()
        self.primary_beamcorr(vexfile, source)
        gainuse, snver = self.ibc1_utils()
        return gainuse, snver
    
    
    def target_cal_inbeam_referencing(self, supportfiledir, prIBCkey):
        """
        In-beam phase referencing
        It copies the required SN tables from the inbeam1 data set and apply them to pulsar.
        Return: tuple (gainuse, snver)
        """
        uvdata    = self.uvdata
        expconfig = self.expconfig
        refant    = expconfig['refant']
        tabledir  = self.tabledir
        
        # gainuse, snver = self.target_or_ibc2_cal(supportfiledir, expconfig['target'])
        gainuse, snver = self.apply_ibc1_solutions2others(supportfiledir, expconfig['target'])

        snver += 1
        vct.loadtable(uvdata, tabledir + f"/{prIBCkey}.calib.p1.sn", snver)
        vct.clcal_applySN(uvdata, refant, [''], snver, '2PT', 'CALP')
        gainuse += 1

        snver += 1
        vct.loadtable(uvdata, tabledir + f"/{prIBCkey}.calib.p2.sn", snver)
        vct.clcal_applySN(uvdata, refant, [''], snver, '2PT', 'CALP')
        gainuse += 1

        snver += 1
        vct.loadtable(uvdata, tabledir + f"/{prIBCkey}.calib.ap.sn", snver)
        vct.clcal_applySN(uvdata, refant, [''], snver, '2PT', 'CALP')
        gainuse += 1
        
        gainuse, snver = self.apply_scintillation_correction()
        return gainuse, snver
    

    def target_utils(self):
        """
        Apply inverse phase-referencing solution to the IBC
        Return: tuple: gainuse, snver
        """
        uvdata    = self.uvdata
        tabledir  = self.tabledir
        expconfig = self.expconfig
        refant    = expconfig['refant']
        gainuse, snver = self.check_gainuse_snver()

        print(f"VLBI PIPE: I am using 'TARGET_UTILS' method with SN{snver} and CL{gainuse}")
        
        snver += 1
        vct.loadtable(uvdata, tabledir + '/target.calib.p1.sn', snver)
        vct.clcal_applySN(uvdata, refant, [''], snver, '2PT', opcode='CALP')
        gainuse += 1

        snver += 1
        vct.loadtable(uvdata, tabledir + '/target.calib.p2.sn', snver)
        vct.clcal_applySN(uvdata, refant, [''], snver, '2PT', 'CALP')
        gainuse += 1
        return gainuse, snver
    

    def pulsar_selfcal_advanceed(self, epoch, vexfile, psrpm_file):
        """
        It runs the pulsar selfcal function but if different epoch model is supplied then
        first, it add a shift to the image centre based on the pulsar proper motion taken from psrpm_info file
        Ashish updated it on 18DEC23, earlier it was taking the first epoch from epochs so the model
        should belong to the same experiment but now you can supply pulsar model that may belong to other experiment.
        eg: I have used BD174 pulsar model to calibrate the BD152 epochs. 
        """
        expconfig  = self.expconfig
        projectdir = expconfig['projectdir']
        targetmodelfile = expconfig['targetmodelfile']
        targetmodelon = expconfig['targetmodelon']
        epochdir = self.epochdir
        uvdata = self.uvdata
        source = expconfig['target']
        modelepochvexfile = expconfig['modelepochvexfile']

        # I will use the pulsar model spacified in expconfig to calibrate all epochs
        # but with a shifted position of source so the delay should be minimum
        if targetmodelfile != None and targetmodelfile != '':
            # It will store the model in the targetmodel directory
            targetmodeldir = projectdir + 'targetmodel/'
            if not os.path.exists(targetmodeldir):
                os.mkdir(targetmodeldir)
            oldmodelfits  = targetmodelfile.split('/')[-1]
            destfile = targetmodeldir + oldmodelfits
            if not os.path.exists(destfile):
                shutil.copy(expconfig['targetmodelfile'], destfile)
            
            # add a shift to image centre if fits epoch and current reduction epoch are not same
            if epoch != oldmodelfits[0:2]:
                rashiftdeg, decshiftdeg = vcs.position_offset_due2_proper_motion(modelepochvexfile, vexfile, psrpm_file, source)
                print("VLBI PIPE: target position shift (mas)", (rashiftdeg*36e5, decshiftdeg*36e5))

                newmodelfits = targetmodeldir + oldmodelfits.replace(oldmodelfits[0:2], epoch)
                vcs.shift_image_center(targetmodelfile, newmodelfits, rashiftdeg, decshiftdeg)
                targetmodelon = targetmodelon.replace(targetmodelon[0:2], epoch)
                targetmodelfile = newmodelfits

            targetmodel = load_model2aips(targetmodelfile, targetmodelon)
            if epoch != oldmodelfits[0:2]: 
                os.remove(targetmodelfile)
            if targetmodel.exists():
                self.pulsar_selfcal(targetmodel)
            vct.delAIPSCat(targetmodel)
        else:
            aips2disk(uvdata, expconfig['target'], epochdir, 'split')

#end================================
##
###