from lightning_model import EnceliaModel
from data_processing import getAllImagesDataset
import torch
import torch.utils.data
import torch.nn.functional as F
import csv
import re
import glob
from pathlib import Path
import os.path
from argparse import ArgumentParser

argp = ArgumentParser(
    description='Uses a cross-validation ensemble to analyze images.'
)
argp.add_argument(
    '-c', '--cv_dir', type=str, required=True, 
    help='The path to a cross-validation output directory.'
)
argp.add_argument(
    '-i', '--images', type=str, required=True,
    help='The path to a collection of images.'
)
argp.add_argument(
    '-o', '--output', type=str, required=False, default='',
    help='The path of an output CSV file.'
)

args = argp.parse_args()

# Get the cross-validation model checkpoints.
ckpts = []
cv_dir = args.cv_dir
ckpts = glob.glob(cv_dir + "/*/*/*/checkpoints/*")
print(ckpts)

all_images = getAllImagesDataset(args.images)

if args.output != '':
    writer = csv.DictWriter(
        open(args.output, 'w'),
        ['file', 'prediction', '0', '1', '2']
    )
    writer.writeheader()
else:
    writer = None

rowout = {}
models = []

with torch.no_grad():
    for i, ckpt in enumerate(ckpts):
        print(f'Loading best model from fold {i}...')
        model = EnceliaModel.load_from_checkpoint(ckpt, lr=0.1)
        model.eval()
        models.append(model)

    i = 0
    correct_cnt = 0
    for img, label in all_images:
        imgfile = all_images.samples[i][0]

        outputs = []
        for model in models:
            output = model(torch.unsqueeze(img, 0))
            output = F.softmax(output, 1)
            outputs.append(torch.squeeze(output))

        outputs = torch.stack(outputs)
        #print(outputs)
        model_avg = torch.mean(outputs, 0)
        #print(model_avg)
        p_label = torch.max(model_avg, 0).indices
        print(i, int(p_label), imgfile)
        #print(p_labels, all_images.samples[i])

        if writer is not None:
            rowout['file'] = os.path.basename(imgfile)
            rowout['prediction'] = int(p_label)
            rowout['0'] = float(model_avg[0])
            rowout['1'] = float(model_avg[1])
            rowout['2'] = float(model_avg[2])
            writer.writerow(rowout)

        i = i + 1
print(i, 'total images processed.')
