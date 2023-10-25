import datetime
import random

import torch
from scipy.io import loadmat
from torch import nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset
import os
# Scientific computing
import numpy as np
import scipy as sp
import scipy.linalg as lin
import scipy.ndimage as ndim
from scipy import io
from scipy.sparse.linalg import svds
from scipy import signal

import torch

# Plotting
import cv2
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image
from torchvision.transforms import Resize, Compose, ToTensor, Normalize
import numpy as np
import skimage
import matplotlib.pyplot as plt

import time
from torch.autograd import Variable

from model import Siren

def get_mgrid(sidelen, dim=2):
    '''Generates a flattened grid of (x,y,...) coordinates in a range of -1 to 1.
    sidelen: int
    dim: int'''
    tensors = tuple(dim * [torch.linspace(-1, 1, steps=sidelen)])
    mgrid = torch.stack(torch.meshgrid(*tensors), dim=-1)
    mgrid = mgrid.reshape(-1, dim)
    return mgrid

class CompressiveImagingReal(Dataset):
    def __init__(self, sidelength,TXM):
        super().__init__()

        self.coords = get_mgrid(sidelength, 2)

        self.TXM = TXM

    def __len__(self):
        return 1

    def __getitem__(self, idx):

        if idx > 0: raise IndexError

        return self.coords, self.TXM

# def laplace(y, x):
#     grad = gradient(y, x)
#     return divergence(grad, x)

def set_dtype(CUDA):
    if CUDA: # if cuda is available
        return torch.cuda.FloatTensor
    else:
        return torch.FloatTensor
def divergence(y, x):
    div = 0.
    for i in range(y.shape[-1]):
        div += torch.autograd.grad(y[..., i], x, torch.ones_like(y[..., i]), create_graph=True)[0][..., i:i+1]
    return div

def gradient_loss_fn(gen_frames, gt_frames, gt_dx, gt_dy, alpha=1):

    # gradient
    gen_dx, gen_dy = cal_gradient(gen_frames)
    #
    # gt_dx, gt_dy = cal_gradient(gt_frames)

    grad_diff_x = torch.abs(gt_dx - gen_dx)
    grad_diff_y = torch.abs(gt_dy - gen_dy)

    # condense into one tensor and avg
    return torch.mean(grad_diff_x ** alpha) + torch.mean(grad_diff_y ** alpha)

def gradient_loss(gen_frames, gt_frames, alpha=1):

    def gradient(x):
        # idea from tf.image.image_gradients(image)
        # https://github.com/tensorflow/tensorflow/blob/r2.1/tensorflow/python/ops/image_ops_impl.py#L3441-L3512
        # x: (b,c,h,w), float32 or float64
        # dx, dy: (b,c,h,w)

        h_x = x.size()[-2]
        w_x = x.size()[-1]
        # gradient step=1
        left = x
        right = F.pad(x, [0, 1,0,0])[:, 1:]
        top = x
        bottom = F.pad(x, [0,0,0, 1])[1:, :]

        # dx, dy = torch.abs(right - left), torch.abs(bottom - top)
        dx, dy = right - left, bottom - top
        # dx will always have zeros in the last column, right-left
        # dy will always have zeros in the last row,    bottom-top
        dx[:, -1] = 0
        dy[-1, :] = 0

        return dx, dy

    # gradient
    gen_dx, gen_dy = gradient(gen_frames)
    gt_dx, gt_dy = gradient(gt_frames)
    #
    grad_diff_x = torch.abs(gt_dx - gen_dx)
    grad_diff_y = torch.abs(gt_dy - gen_dy)

    # condense into one tensor and avg
    return torch.mean(grad_diff_x ** alpha) + torch.mean(grad_diff_y ** alpha)

def cal_gradient(x):
    # idea from tf.image.image_gradients(image)
    # https://github.com/tensorflow/tensorflow/blob/r2.1/tensorflow/python/ops/image_ops_impl.py#L3441-L3512
    # x: (b,c,h,w), float32 or float64
    # dx, dy: (b,c,h,w)

    h_x = x.size()[-2]
    w_x = x.size()[-1]
    # gradient step=1
    left = x
    right = F.pad(x, [0, 1,0,0])[:, 1:]
    top = x
    bottom = F.pad(x, [0,0,0, 1])[1:, :]

    # dx, dy = torch.abs(right - left), torch.abs(bottom - top)
    dx, dy = right - left, bottom - top
    # dx will always have zeros in the last column, right-left
    # dy will always have zeros in the last row,    bottom-top
    dx[:, -1] = 0
    dy[-1, :] = 0

    return dx, dy

class CompressiveImaging(Dataset):
    def __init__(self, sidelength,A,noiselevel,Gx, Gy, img):
        super().__init__()

        self.Gx = Gx
        self.Gy = Gy

        self.ground_truth = img
        self.pixels = img.permute(1, 2, 0).view(-1, 1)
        self.coords = get_mgrid(sidelength, 2)

        self.compressed = A@self.pixels

        # eta_sig = noiselevel  # set value to induce noise
        # eta = np.random.normal(0, eta_sig * (1.0 / self.compressed[0].size()[0]), self.compressed[0].size()[0])
        # eta = torch.from_numpy(eta).float()
        # eta = eta.unsqueeze(1)
        # self.compressed = self.compressed + eta

    def __len__(self):
        return 1

    def __getitem__(self, idx):

        np.random.seed(datetime.datetime.now().second + datetime.datetime.now().microsecond)
        if idx > 0: raise IndexError

        return self.coords, self.ground_truth, self.compressed, self.Gx, self.Gy

def load_data(filename, sampling_ratio, image_size, noise_level=0):
    # load measurement matrix
    A = loadmat(filename+'_size_' + str(image_size) + "_sampling_" + str(sampling_ratio) + '.mat')['h']
    A = torch.from_numpy(A).to(torch.float32).T
    img = loadmat(filename+'_size_' + str(image_size) + "_sampling_" + str(sampling_ratio) + '.mat')['z']
    img = Image.fromarray(img)
    transform = Compose([
        Resize(image_size),
        ToTensor(),
        # Normalize(torch.Tensor([0.5]), torch.Tensor([0.5]))
    ])
    img = transform(img)
    Gx, Gy = cal_gradient(img[0,:,:])

    data = CompressiveImaging(image_size, A, noise_level, Gx, Gy, img)

    dataloader = DataLoader(data, batch_size=1, pin_memory=True, num_workers=0)
    model_input, ground_truth, compressed, Gx, Gy = next(iter(dataloader))
    model_input, ground_truth, compressed, Gx, Gy,A = model_input.cuda(), ground_truth.cuda(), compressed.cuda(), Gx.cuda(), Gy.cuda(), A.cuda()
    return model_input, ground_truth, compressed, Gx, Gy,A

def normalize(x, fullnormalize=False):
    '''
        Normalize input to lie between 0, 1.

        Inputs:
            x: Input signal
            fullnormalize: If True, normalize such that minimum is 0 and
                maximum is 1. Else, normalize such that maximum is 1 alone.

        Outputs:
            xnormalized: Normalized x.
    '''

    if x.sum() == 0:
        return x

    xmax = x.max()

    if fullnormalize:
        xmin = x.min()
    else:
        xmin = 0

    xnormalized = (x - xmin) / (xmax - xmin)

    return xnormalized


def rsnr(x, xhat):
    '''
        Compute reconstruction SNR for a given signal and its reconstruction.

        Inputs:
            x: Ground truth signal (ndarray)
            xhat: Approximation of x

        Outputs:
            rsnr_val: RSNR = 20log10(||x||/||x-xhat||)
    '''
    xn = lin.norm(x.reshape(-1))
    en = lin.norm((x - xhat).reshape(-1))
    rsnr_val = 20 * np.log10(xn / en)

    return rsnr_val


def psnr(x, xhat):
    ''' Compute Peak Signal to Noise Ratio in dB

        Inputs:
            x: Ground truth signal
            xhat: Reconstructed signal

        Outputs:
            snrval: PSNR in dB
    '''
    err = x - xhat
    denom = np.mean(pow(err, 2))

    snrval = 10 * np.log10(np.max(x) / denom)

    return snrval


def measure(x, noise_snr=40, tau=100):
    ''' Realistic sensor measurement with readout and photon noise

        Inputs:
            noise_snr: Readout noise in electron count
            tau: Integration time. Poisson noise is created for x*tau.
                (Default is 100)

        Outputs:
            x_meas: x with added noise
    '''
    x_meas = np.copy(x)

    noise = np.random.randn(x_meas.size).reshape(x_meas.shape) * noise_snr

    # First add photon noise, provided it is not infinity
    if tau != float('Inf'):
        x_meas = x_meas * tau

        x_meas[x > 0] = np.random.poisson(x_meas[x > 0])
        x_meas[x <= 0] = -np.random.poisson(-x_meas[x <= 0])

        x_meas = (x_meas + noise) / tau

    else:
        x_meas = x_meas + noise

    return x_meas


def build_montage(images):
    '''
        Build a montage out of images
    '''
    nimg, H, W = images.shape

    nrows = int(np.ceil(np.sqrt(nimg)))
    ncols = int(np.ceil(nimg / nrows))

    montage_im = np.zeros((H * nrows, W * ncols), dtype=np.float32)

    cnt = 0
    for r in range(nrows):
        for c in range(ncols):
            h1 = r * H
            h2 = (r + 1) * H
            w1 = c * W
            w2 = (c + 1) * W

            if cnt == nimg:
                break

            montage_im[h1:h2, w1:w2] = normalize(images[cnt, ...], True)
            cnt += 1

    return montage_im


def count_parameters(model):
    return sum(p.numel() for p in model.parameters() if p.requires_grad)


def get_coords(H, W, T=None):
    '''
        Get 2D/3D coordinates
    '''
    if T is None:
        X, Y = np.meshgrid(np.linspace(-1, 1, W), np.linspace(-1, 1, H))
        coords = np.hstack((X.reshape(-1, 1), Y.reshape(-1, 1)))
    else:
        X, Y, Z = np.meshgrid(np.linspace(-1, 1, W),
                              np.linspace(-1, 1, H),
                              np.linspace(-1, 1, T))
        coords = np.hstack((X.reshape(-1, 1),
                            Y.reshape(-1, 1),
                            Z.reshape(-1, 1)))

    return torch.tensor(coords.astype(np.float32))


def resize(cube, scale):
    '''
        Resize a multi-channel image

        Inputs:
            cube: (H, W, nchan) image stack
            scale: Scaling
    '''
    H, W, nchan = cube.shape

    im0_lr = cv2.resize(cube[..., 0], None, fx=scale, fy=scale)
    Hl, Wl = im0_lr.shape

    cube_lr = np.zeros((Hl, Wl, nchan), dtype=cube.dtype)

    for idx in range(nchan):
        cube_lr[..., idx] = cv2.resize(cube[..., idx], None,
                                       fx=scale, fy=scale,
                                       interpolation=cv2.INTER_AREA)
    return cube_lr



@torch.no_grad()
def get_layer_outputs(model, coords, imsize,
                      nfilters_vis=16,
                      get_imag=False):
    '''
        get activation images after each layer

        Inputs:
            model: INR model
            coords: 2D coordinates
            imsize: Size of the image
            nfilters_vis: Number of filters to visualize
            get_imag: If True, get imaginary component of the outputs

        Outputs:
            atoms_montages: A list of 2d grid of outputs
    '''
    H, W = imsize

    if model.pos_encode:
        coords = model.positional_encoding(coords)

    atom_montages = []

    for idx in range(len(model.net) - 1):
        layer_output = model.net[idx](coords)
        layer_images = layer_output.reshape(1, H, W, -1)[0]

        if nfilters_vis is not 'all':
            layer_images = layer_images[..., :nfilters_vis]

        if get_imag:
            atoms = layer_images.detach().cpu().numpy().imag
        else:
            atoms = layer_images.detach().cpu().numpy().real

        atoms_min = atoms.min(0, keepdims=True).min(1, keepdims=True)
        atoms_max = atoms.max(0, keepdims=True).max(1, keepdims=True)

        signs = (abs(atoms_min) > abs(atoms_max))
        atoms = (1 - 2 * signs) * atoms

        # Arrange them by variance
        atoms_std = atoms.std((0, 1))
        std_indices = np.argsort(atoms_std)

        atoms = atoms[..., std_indices]

        atoms_min = atoms.min(0, keepdims=True).min(1, keepdims=True)
        atoms_max = atoms.max(0, keepdims=True).max(1, keepdims=True)

        atoms = (atoms - atoms_min) / np.maximum(1e-14, atoms_max - atoms_min)

        atoms[:, [0, -1], :] = 1
        atoms[[0, -1], :, :] = 1

        atoms_montage = build_montage(np.transpose(atoms, [2, 0, 1]))

        atom_montages.append(atoms_montage)
        coords = layer_output

    return atom_montages