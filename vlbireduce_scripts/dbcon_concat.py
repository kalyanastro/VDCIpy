#!/usr/bin/ParselTongue

"""
Ashish Kalyan & Adam Deller

This script is build on top of the dbcon function. 
Ashish faced an issue while dbconing 29 file, "obit error: Too many open files". 
To resolve it, this script is written, it does the same job but concatenate in chunks.

It first split the desired source by applying calibration solutions and then perform the 
concatenation job.

output 
    missing an list shows the antennas those are missing in concatenating file.

    usages: 
> ParselTongue dbcon_concat.py -h
"""

from AIPS import AIPS
from AIPSData import AIPSUVData
import vlbi_calibration_tasks as vct
import astro_utils as ast
import yaml, os, argparse, sys
from argparse import RawTextHelpFormatter


yamldir     = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "/examples"
expconfig   = yaml.load(open(yamldir + '/expconfig.yaml'), yaml.SafeLoader)
aipsver     = expconfig['aipsver']
AIPS.userno = expconfig['userno']
epochs  = expconfig['epochs']
projdir = expconfig['projectdir']
band    = str(expconfig['band'])
bif = str(expconfig['bif'])
eif = str(expconfig['eif'])
original_stdout = sys.stdout
splitclass = 'temp'
suffix = band
if expconfig['splitifs']:
    suffix = bif + eif


parser = argparse.ArgumentParser(description='Concatnate AIPS data sets', formatter_class=RawTextHelpFormatter)
parser.add_argument('-s',    '--source',   metavar='', help='Source that you want to concatenate')
parser.add_argument('-ic',   '--inclass',  metavar='', help='AIPSCat (multisource file) class')
parser.add_argument('-ikey', '--idikey',   metavar='', help='idikey, IBC1/2/3 or PSRG/U')
args = parser.parse_args()

if args.source  != None: source  = args.source.upper()
if args.inclass != None: inclass = args.inclass
if args.idikey  != None: idikey  = args.idikey


gblmodeldir  = projdir + '/globalmodels'
sourcegbldir = gblmodeldir + '/gbl' + idikey + suffix + '_' + source
if not os.path.exists(sourcegbldir):
    if not os.path.exists(gblmodeldir):
        os.mkdir(gblmodeldir)
    os.mkdir(sourcegbldir)
diskoutname = sourcegbldir + '/gbl' + idikey + suffix + '_' + source + '_aips.fits'


# write a log file in globalmodel directory
logfilepath = sourcegbldir + f"/dbcon.{source}.log"
logfile = open(logfilepath, 'w')
sys.stdout = ast.Logger(logfile)


# split the source that has to concatenate
for epoch in epochs:
    inname = epoch + idikey + suffix 
    uvdata = AIPSUVData(inname, inclass, 1,1)
    vct.split_uvdata(uvdata, [source], 0, splitclass, outseq=0)


numcat = len(epochs)
chunksize = 10
if numcat < chunksize:
    batches = 1
    chunksize = numcat
else:
    remander = numcat%chunksize
    batches  = int(numcat/chunksize)
    if remander != 0:
        batches = int((numcat/chunksize)) + 1
     
if batches > chunksize:
    print(f"VLBI PIPE: Bactches should be less then {chunksize} otherwise edit the script,"+
                    "you may increase the chunksize")
print(f"VLBI PIPE: batches and chunksize {batches}, {chunksize}")

inputuvdatas = []
for i in range(numcat):
    inputuvdatas.append(AIPSUVData(source, splitclass, 1, i+1))

inputuvdatas.sort(key=lambda a: len(a.antennas), reverse=True)
print(*inputuvdatas, sep = "\n") #*******    


tempdboc    = 'DBCT'
uvdatachunk, uvdatabatch = [], []
for i in range(batches):
    for j in range(chunksize):
        k = i*chunksize + j
        if k > len(inputuvdatas)-1:
            break
        uvdatachunk.append(inputuvdatas[k])
    # print(uvdatachunk)
    if len(uvdatachunk) > 1:   
        dboutdata  = AIPSUVData(source, tempdboc, 1, i+1)
        uvdatabatch.append(dboutdata)
        vct.dbcon(uvdatachunk, dboutdata)
    # print(uvdatabatch)
    elif len(uvdatachunk) == 1:
        uvdatabatch.append(uvdatachunk[0])
    uvdatachunk.clear()

print(*uvdatabatch, sep = "\n")
if len(uvdatabatch) == 1:
    vct.fittp(uvdatabatch[0], diskoutname)
else:
    finaldboc = 'DBC' + suffix
    dboutdata = AIPSUVData(source, finaldboc, 1, 1)
    vct.dbcon(uvdatabatch, dboutdata)
    vct.fittp(dboutdata, diskoutname)


for indata in uvdatabatch:
    if indata.klass == tempdboc:
        vct.delAIPSCat(indata)
        

for i in range(numcat):
    splitdata = AIPSUVData(source, splitclass, 1, i+1)
    vct.delAIPSCat(splitdata)

logfile.close()
sys.stdout = original_stdout

#end====================================
##
###