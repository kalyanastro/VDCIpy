#!/usr/bin/python3

"""
Ashish Kalyan

The script renames the idifits file and make their name compatible with the pipeline.
"idi_rename" function looks for '_0_'/'_1_' and '_un/gated_'/'_inbeam1/2/3_' keys
in idifits file name and then rename accordingly and transfer all data files in 
'<epoch>_data' directory.
eg. output name: <obscode><epoch><idikey><band>.idifits == BD174A1IBC10.idifits, 
    where band is equal to 0/1 (generally, it will be zero, except for the dual band observation like BD174.)
"""

import os, shutil, argparse

parser = argparse.ArgumentParser('usages: file_rename -oc bd152 -pd /home/Desktop/project/bd152/ -ep a1,a2')
parser.add_argument('-oc', '--obscode',  metavar='',  help='OBSCODE')
parser.add_argument('-pd', '--projdir',  metavar='',  help='project directory')
parser.add_argument('-ep', '--epochs',   metavar='',  help='epochs (list eg. -ep a1,a2,a3)')
args   = parser.parse_args()

obscode = args.obscode
projdir = args.projdir
epochs  = args.epochs.split(',')


def rename_file(old_name, new_name):
    """
    Rename directory/file
    Return: string: new directory/file
    """
    try:
        os.rename(old_name, new_name)
    except FileNotFoundError:
        print(f"The file '{old_name}' does not exist.")
    return new_name


def transfer_files2folder(dirpath, folder_name):
    """
    It transfers all files in the "dirpath" directory to a "folder_name" folder,
    which will be in the same directory
    Return: string: Target folder path
    """
    try:
        folder_path = os.path.join(dirpath, folder_name)
        if not os.path.exists(folder_path):
            os.mkdir(folder_path)

        for item in os.listdir(dirpath):
            item_path = os.path.join(dirpath, item)
            
            # Check if the item is a file and it's not the target folder itself
            if os.path.isfile(item_path) and item != folder_name:
                shutil.move(item_path, folder_path)
        return folder_path

    except FileNotFoundError:
        print(f"The directory '{dirpath}' does not exist.")


def idi_rename(datadir):
    """
    It renames the idifiles stored in the "datadir" directroy
    """
    idilist = os.listdir(datadir)
    print(f"VLBI PIPE: idi_rename: {len(idilist)} idifits file/s found in {datadir}")
    for filename in idilist:
        oldfile = datadir + '/' + filename
        if '_0_' in filename and '_gated' in filename:
            newfile = datadir + '/' + obscode + epoch + 'psrg0.idifits'
            rename_file(oldfile, newfile)

        elif '_0_' in filename and '_ungated' in filename:
            newfile = datadir + '/' + obscode + epoch + 'psru0.idifits'
            rename_file(oldfile, newfile)
            
        elif '_0_' in filename and '_inbeam1' in filename:
            newfile = datadir + '/' + obscode + epoch + 'ibc10.idifits'
            rename_file(oldfile, newfile)
            
        elif '_0_' in filename and 'inbeam2' in filename:
            newfile = datadir + '/' + obscode + epoch + 'ibc20.idifits'
            rename_file(oldfile, newfile)
            
        elif '_0_' in filename and 'inbeam3' in filename:
            newfile = datadir + '/' + obscode + epoch + 'ibc30.idifits'
            rename_file(oldfile, newfile)

        elif '_1_' in filename and '_gated' in filename:
            newfile = datadir + '/' + obscode + epoch + 'psrg1.idifits'
            rename_file(oldfile, newfile)
        
        elif '_1_' in filename and '_ungated' in filename:
            newfile = datadir + '/' + obscode + epoch + 'psru1.idifits'
            rename_file(oldfile, newfile)
        
        elif '_1_' in filename and '_inbeam1' in filename:
            newfile = datadir + '/' + obscode + epoch + 'ibc11.idifits'
            rename_file(oldfile, newfile)
        
        elif '_1_' in filename and '_inbeam2' in filename:
            newfile = datadir + '/' + obscode + epoch + 'ibc21.idifits'
            rename_file(oldfile, newfile)
        
        elif '_1_' in filename and '_inbeam3' in filename:
            newfile = datadir + '/' + obscode + epoch + 'ibc31.idifits'
            rename_file(oldfile, newfile)


if __name__ == '__main__':
    for epoch in epochs:
        olddir = projdir + '/' + obscode.upper() + epoch.upper()
        epochdir = projdir + '/' + epoch.lower()
        rename_file(olddir, epochdir)

        new_folder = epoch.lower() + '_data'
        datadir = transfer_files2folder(epochdir, new_folder)

        idi_rename(datadir)

#end===================================
##
###
