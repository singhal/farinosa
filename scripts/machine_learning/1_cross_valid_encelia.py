#Salamander color morph machine learning 

import os
from lightning_model import EnceliaModel
from data_processing import getDatasets, getDataLoaders, getAllImagesDataset
import numpy as np
import pytorch_lightning as pl
from pytorch_lightning import loggers as pl_loggers
from pytorch_lightning.callbacks import ModelCheckpoint
import sys
import os.path
import pathlib
from sklearn.model_selection import KFold
import argparse

parser = argparse.ArgumentParser(
                    prog='Runs machine laerning',
                    description='Runs machine laerning')
# folder with data
parser.add_argument('-i', '--img') 
# output folder
parser.add_argument('-o', '--outdir') 
# learning rate
parser.add_argument('-l', '--learn') 

args = parser.parse_args()
learning_rate = float(args.learn)
image_dir = args.img

if not os.path.isdir(args.outdir):
    os.mkdir(args.outdir)

kf = KFold(n_splits=4, shuffle=True) 
all_images = getAllImagesDataset(image_dir)
num_images = len(all_images)
indices = list(range(num_images))

loop_count = 0
for train_idx, valid_idx in kf.split(indices): 
    # print('train_idx: %s, valid_idx: %s' % (train_idx, valid_idx))

    fold_folder = os.path.join(args.outdir, 'fold_%s' % loop_count)
    os.mkdir(fold_folder)
    print(f'Cross-validation fold {loop_count}; saving results to {fold_folder}.')
    
    #rng = np.random.default_rng(seed=12345)

    print("* getting datasets")
    train_data, val_data = getDatasets(
        image_dir,
        train_idx=train_idx, valid_idx=valid_idx
    )
    print("* loading datasets")
    trainloader, valloader = getDataLoaders(train_data, val_data)

    print("* building model")
    model = EnceliaModel(lr=learning_rate)
    print("* logging tensorboard")
    tb_logger = pl_loggers.TensorBoardLogger(
        fold_folder, f'encelia_exp_cross_val-{model.lr}'
    )

    print("* checkpointing model")
    checkpoint_callback = ModelCheckpoint(
        dirpath=fold_folder,
        save_top_k=1,
        verbose=True,
        monitor='valid_loss',
        mode='min',
        filename=f'weight-{model.lr}'
    )

    print("* training model")
    #trainer = pl.Trainer(logger=tb_logger)
    trainer = pl.Trainer(logger=tb_logger, max_epochs=100)

    print("* validating model")
    trainer.fit(model, trainloader, valloader)

    loop_count += 1
