import numpy as np
import torch
from torch import nn
import pdb
import math

import numpy as np

import torch
from torch import nn

class SineLayer(nn.Module):
    # See paper sec. 3.2, final paragraph, and supplement Sec. 1.5 for discussion of omega_0.

    # If is_first=True, omega_0 is a frequency factor which simply multiplies the activations before the
    # nonlinearity. Different signals may require different omega_0 in the first layer - this is a
    # hyperparameter.

    # If is_first=False, then the weights will be divided by omega_0 so as to keep the magnitude of
    # activations constant, but boost gradients to the weight matrix (see supplement Sec. 1.5)

    def __init__(self, in_features, out_features, bias=True,
                 is_first=False, omega_0=30):
        super().__init__()
        self.omega_0 = omega_0
        self.is_first = is_first

        self.in_features = in_features
        self.linear = nn.Linear(in_features, out_features, bias=bias)


        self.init_weights()

    def init_weights(self):
        with torch.no_grad():
            if self.is_first:
                self.linear.weight.uniform_(-1 / self.in_features,
                                            1 / self.in_features)
            else:
                self.linear.weight.uniform_(-np.sqrt(6 / self.in_features) / self.omega_0,
                                            np.sqrt(6 / self.in_features) / self.omega_0)

    def forward(self, input):
        return torch.sin(self.omega_0 * self.linear(input))

    def forward_with_intermediate(self, input):
        # For visualization of activation distributions
        intermediate = self.omega_0 * self.linear(input)
        return torch.sin(intermediate), intermediate


class Siren(nn.Module):
    def __init__(self, in_features, hidden_features, hidden_layers, out_features, outermost_linear=False,
                 first_omega_0=30, hidden_omega_0=30.):
        super().__init__()

        self.net = []
        self.net.append(SineLayer(in_features, hidden_features,
                                  is_first=True, omega_0=first_omega_0))

        # self.net.append(nn.Dropout(0.1))
        for i in range(hidden_layers):
            self.net.append(SineLayer(hidden_features, hidden_features,
                                      is_first=False, omega_0=hidden_omega_0))


        if outermost_linear:
            final_linear = nn.Linear(hidden_features, out_features)

            with torch.no_grad():
                final_linear.weight.uniform_(-np.sqrt(6 / hidden_features) / hidden_omega_0,
                                             np.sqrt(6 / hidden_features) / hidden_omega_0)

            # self.net.append(final_linear)
            # self.net.append(nn.ReLU())

            self.net.append(final_linear)

        else:
            self.net.append(SineLayer(hidden_features, out_features,
                                      is_first=False, omega_0=hidden_omega_0))

        self.net = nn.Sequential(*self.net)



    def forward(self, coords):
        coords = coords.clone().detach().requires_grad_(True) # allows to take derivative w.r.t. input
        output = self.net(coords)
        return output, coords


class ReLULayer(nn.Module):
    '''
        Drop in replacement for SineLayer but with ReLU non linearity
    '''

    def __init__(self, in_features, out_features, bias=True,
                 is_first=False, omega_0=30, scale=10.0):
        '''
            is_first, and omega_0 are not used.
        '''
        super().__init__()
        self.in_features = in_features
        self.omega_0 = omega_0
        self.is_first = is_first
        self.linear = nn.Linear(in_features, out_features, bias=bias)

    def forward(self, input):
        return nn.functional.relu(self.linear(input))


class PosEncoding(nn.Module):
    '''Module to add positional encoding as in NeRF [Mildenhall et al. 2020].'''

    def __init__(self, in_features, sidelength=None, fn_samples=None, use_nyquist=True):
        super().__init__()

        self.in_features = in_features

        if self.in_features == 3:
            self.num_frequencies = 10
        elif self.in_features == 2:
            assert sidelength is not None
            if isinstance(sidelength, int):
                sidelength = (sidelength, sidelength)
            self.num_frequencies = 4
            if use_nyquist:
                self.num_frequencies = self.get_num_frequencies_nyquist(min(sidelength[0], sidelength[1]))
        elif self.in_features == 1:
            # assert fn_samples is not None
            fn_samples = sidelength
            self.num_frequencies = 4
            if use_nyquist:
                self.num_frequencies = self.get_num_frequencies_nyquist(fn_samples)
        else:
            self.num_frequencies = 4

        self.out_dim = in_features + 2 * in_features * self.num_frequencies

    def get_num_frequencies_nyquist(self, samples):
        nyquist_rate = 1 / (2 * (2 * 1 / samples))
        return int(math.floor(math.log(nyquist_rate, 2)))

    def forward(self, coords):
        coords = coords.view(coords.shape[0], -1, self.in_features)

        coords_pos_enc = coords
        for i in range(self.num_frequencies):
            for j in range(self.in_features):
                c = coords[..., j]

                sin = torch.unsqueeze(torch.sin((2 ** i) * np.pi * c), -1)
                cos = torch.unsqueeze(torch.cos((2 ** i) * np.pi * c), -1)

                coords_pos_enc = torch.cat((coords_pos_enc, sin, cos), axis=-1)

        return coords_pos_enc.reshape(coords.shape[0], -1, self.out_dim)


class INR(nn.Module):
    def __init__(self, in_features,
                 hidden_features, hidden_layers,
                 out_features, outermost_linear=True,
                 first_omega_0=30, hidden_omega_0=30., scale=10.0,
                 pos_encode=True, sidelength=512, fn_samples=None,
                 use_nyquist=True):
        super().__init__()
        self.pos_encode = pos_encode

        self.complex = False
        self.nonlin = ReLULayer

        if pos_encode:
            self.positional_encoding = PosEncoding(in_features=in_features,
                                                   sidelength=sidelength,
                                                   fn_samples=fn_samples,
                                                   use_nyquist=use_nyquist)
            in_features = self.positional_encoding.out_dim

        self.net = []
        self.net.append(self.nonlin(in_features, hidden_features,
                                    is_first=True, omega_0=first_omega_0,
                                    scale=scale))

        #

        for i in range(hidden_layers):
            self.net.append(self.nonlin(hidden_features, hidden_features,
                                        is_first=False, omega_0=hidden_omega_0,
                                        scale=scale))
        self.net.append(nn.Dropout(0.1))

        if outermost_linear:
            if self.complex:
                dtype = torch.cfloat
            else:
                dtype = torch.float
            final_linear = nn.Linear(hidden_features,
                                     out_features,
                                     dtype=dtype)

            self.net.append(final_linear)
        else:
            self.net.append(self.nonlin(hidden_features, out_features,
                                        is_first=False, omega_0=hidden_omega_0,
                                        scale=scale))

        self.net = nn.Sequential(*self.net)

    def forward(self, coords):
        if self.pos_encode:
            coords = self.positional_encoding(coords)

        coords = coords.clone().detach().requires_grad_(True)
        output = self.net(coords)

        return output,coords