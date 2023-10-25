# Reconstruction of fluorescence map
from datetime import datetime

import numpy as np
import torch
from PIL import Image
from matplotlib import pyplot as plt
from pytorch_msssim import ssim
from scipy.io import loadmat, savemat
from torch import nn
from torch.optim.lr_scheduler import CosineAnnealingWarmRestarts
from torch.utils.data import DataLoader, Dataset
from torchvision.transforms import Resize, ToTensor, Compose, Normalize
from torchvision.utils import make_grid

from utils import CompressiveImagingReal,gradient_loss
from model import Siren, INR
from skimage import filters
import wandb

image_size = 128
filename = 'data4rec_real_2021_07_mixed_particles_NMC_LCO_LNO'
element = 'Ni'

loss_re, recons_re = [],[]
# load measurement matrix
A = loadmat(filename + '.mat')['Patterns']
A = torch.from_numpy(A).to(torch.float32)
A = torch.transpose(A, 0,1)
# load the XRF measurement
compressed = loadmat(filename + '.mat')['XRF_amount_'+element]
compressed = torch.from_numpy(compressed).to(torch.float32).T
alpha = 8 # may be tuned for better performance, there is a scalling inconsistence (to be explored)
compressed = torch.squeeze(compressed)*alpha
# load the TXM image
TXM = loadmat(filename + '.mat')['TXM']
TXM = Image.fromarray(TXM)
transform = Compose([
    Resize(image_size),
    ToTensor(),
    # Normalize(torch.Tensor([0.5]), torch.Tensor([0.5]))
])
TXM = transform(TXM)
TXM = torch.transpose(TXM,1,2)
data = CompressiveImagingReal(image_size, TXM)
dataloader = DataLoader(data, batch_size=1, pin_memory=True, num_workers=0)
model_input, TXM = next(iter(dataloader))
model_input, TXM, compressed, A = model_input.cuda(),TXM.cuda(),compressed.cuda(), A.cuda()

total_steps = 20000
steps_til_summary = 1000
num_restarts = 1
exit_window = 50

psnr_list = []
taskname = "rec_" + filename + "_hidden_feature_1024_layer_5_lr_5e-6_lambda_1e-5"
wandb.init(project="RealExp",name=taskname)

for j in range(num_restarts):

    img_siren = Siren(in_features=2, out_features=1, hidden_features=1024,  #### related to noise level
                      hidden_layers=5, outermost_linear=True)
    img_siren = img_siren.cuda()

    optim = torch.optim.Adam(lr=5e-6, params=img_siren.parameters()) # important parameter
    # optim = torch.optim.AdamW(img_siren.parameters(), lr=sampling_ratio * 1e-5)

    loss_iter = []
    recons_iter = []

    for i in range(total_steps):

        model_output, coords = img_siren(model_input)
        # model_output = torch.clip(model_output, 0, torch.inf)
        model_output_vec = torch.flatten(torch.transpose(model_output,1,2))
        compressed_out = A @ model_output_vec

        y_loss = ((compressed_out - compressed) ** 2).mean()
        loss = y_loss + 1e-5*gradient_loss(model_output.view(image_size, image_size),
                                             TXM.view(image_size, image_size)) # important parameter

        if i >= total_steps - exit_window:
            recons_iter.append(model_output.cpu().detach().numpy())
            loss_iter.append(loss.cpu().detach().numpy())

        if not i % steps_til_summary:
            print("Step %d, Total loss %0.6f" % (i, loss))
            wandb.log({'loss': loss})
            imglist = [model_output.view(image_size, image_size).t(),TXM.view(image_size, image_size).t()]
            image_array = torch.cat(imglist, 1)
            images = wandb.Image(
                image_array,
                caption="Recon, TXM"
            )
            wandb.log({"examples": images})
            #np.save('process/'+element+'_itr_'+str(i)+'.npy',model_output.view(image_size, image_size).t().cpu().detach().numpy())

        optim.zero_grad()
        loss.backward()
        optim.step()

    idx_itr = np.argmin(loss_iter, axis=0)
    recons_re.append(recons_iter[idx_itr][0])
    loss_re.append(y_loss.data.cpu().numpy())

idx_re = np.argmin(loss_re, axis=0)
x_hat = np.reshape(recons_re[idx_re],(image_size,image_size))

mdic = {"xhat": x_hat, "TXM": TXM.cpu().detach().numpy(), "y": compressed.cpu().detach().numpy(), "A": A.cpu().detach().numpy()}

savemat(taskname+".mat", mdic)


