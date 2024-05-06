# read in csv with encelia manual scores
# move data into folders

import pandas as pd
import shutil
import os

os.chdir('/Users/singhal/Dropbox (Personal)/research_projects-ACTIVE/Encelia/analysis/black_yellow/machine_learning')

outdir = '/Users/singhal/Desktop/machine_learning/'

d = pd.read_csv('Encelia.csv')

# only keeping brown, yellow & no flowers
coltypes = ["B", "Y", "N"]

for coltype in coltypes:
    subdir = os.path.join(outdir, coltype)
    if not os.path.isdir(subdir):
        os.mkdir(subdir)

    dd = d[d.flowers == coltype]

    print(dd.head())

    
    for i, row in dd.iterrows():
        # old file name
        origfile = os.path.join('inat_photos', row["file"])
        # new file name
        newfile = os.path.join(outdir, coltype, row["file"])

        shutil.copyfile(origfile, newfile)
