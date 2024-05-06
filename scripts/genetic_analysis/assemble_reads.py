# import pandas as pd
import subprocess
import os
import re

dir = '/Volumes/heloderma4/sonal/encelia/'

file = os.path.join(dir, "encelia_samples_v3.csv")
f = open(file, 'r')
head = f.next()

for line in f:
	d = re.split(',', line)
	sample = d[0]

	call = "/Volumes/MacHD3.5/sonal/bin/VelvetOptimiser/VelvetOptimiser.pl -k 'tbp' -s 21 -e 91 -x 10 -t 16 -d %s -f '-shortPaired -fastq.gz -separate %strim_reads/%s_R1.final.fq.gz %strim_reads/%s_R2.final.fq.gz -short -fastq.gz %strim_reads/%s_unpaired.final.fq.gz'\n" % (outdir, dir, sample, dir, sample, dir, sample)
	subprocess.call(call, shell=True)
	outfile = os.path.join(outdir, 'contigs.fa')
	subprocess.call("cp %s %s" % (outfile, new), shell=True)

