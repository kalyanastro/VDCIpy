#!/usr/bin/ParselTongue
"""
Ashish Kalyan
This script does the imaginig of calibrated data.
"""

import os, yaml, argparse
from argparse import RawTextHelpFormatter
import vlbi_calibration_tasks as vct


yamldir    = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "/examples"
expconfig  = yaml.load(open(yamldir + '/expconfig.yaml'), yaml.SafeLoader)
epochs     = expconfig['epochs']
obscode    = expconfig['obscode']
projectdir = expconfig['projectdir']
band = str(expconfig['band'])
bif  = str(expconfig['bif'])
eif  = str(expconfig['eif'])
suffix = band
pixsize = expconfig['pixsize']
mapsize = expconfig['mapsize']
if expconfig['splitifs']:
    suffix = bif + eif

parser = argparse.ArgumentParser(description='DIFMAP Imaging', formatter_class=RawTextHelpFormatter)
parser.add_argument('-s',  '--source',     metavar='', help='Source')
parser.add_argument('-i',  '--idikey',     metavar='', help='idikey, ibc1/2/3 or psrg/u')
parser.add_argument('-md', '--modeldiv',   metavar='', help='if Y, means, it is model-divided data')
parser.add_argument('-dg', '--dogaussian', metavar='', help='if T, use gaussian model otherwise point source model')
args = parser.parse_args()

if args.idikey  != None: idikey = args.idikey
if args.source  != None: source = args.source.upper()
if args.dogaussian.upper() == 'T':
    dogaussian = True
else: dogaussian = False

for epoch in epochs:
    foldername = epoch + idikey + suffix + '_' + source + '_model'
    if args.modeldiv.upper() == 'Y':
        fitsname = epoch + idikey + suffix + '_' + source + '_modeldiv_aips.fits'
    else:
        fitsname = epoch + idikey + suffix + '_' + source + '_aips.fits'
    uvfile = projectdir + '/' + epoch + '/source_models/' + foldername + '/' + fitsname
    vct.difmap_maptarget(uvfile, pixsize, dogaussian, mapsize)

#end ==========================
##
###