# Class for block anaysis

import os
import sys
import numpy as np
import math
from numpy import log2, zeros, mean, var, sum, arange, array, cumsum,  floor
import pyblock
import json

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

def load_json(file:str, flag:str='r'):

    return json.load(open(file,flag))

def block(x):
    # preliminaries
    d = log2(len(x))
    if (d - floor(d) != 0):
        x = x[:2**int(floor(d))]
    d = int(floor(d))
    n = 2**d
    s, gamma = zeros(d), zeros(d)
    mu = mean(x)
    # estimate the auto-covariance and variances
    # for each blocking transformation
    for i in arange(0, d):
        n = len(x)
        # estimate autocovariance of x
        gamma[i] = (n)**(-1)*sum((x[0:(n-1)]-mu)*(x[1:n]-mu))
        # estimate variance of x
        s[i] = var(x)
        # perform blocking transformation
        x = 0.5*(x[0::2] + x[1::2])

    # generate the test observator M_k from the theorem
    M = (cumsum(((gamma/s)**2*2**arange(1, d+1)[::-1])[::-1]))[::-1]

    # we need a list of magic numbers
    q = array([6.634897,  9.210340,  11.344867, 13.276704, 15.086272,
            16.811894, 18.475307, 20.090235, 21.665994, 23.209251,
            24.724970, 26.216967, 27.688250, 29.141238, 30.577914,
            31.999927, 33.408664, 34.805306, 36.190869, 37.566235,
            38.932173, 40.289360, 41.638398, 42.979820, 44.314105,
            45.641683, 46.962942, 48.278236, 49.587884, 50.892181])

    # use magic to determine when we should have stopped blocking
    for k in arange(0, d):
        if(M[k] < q[k]):
            break
    if (k >= d-1):
        print("Warning: Use more data")

    return (s[k]/2**(d-k))

def chunkIt(a:int, num:int):
    avg = a / float(num)
    out = []
    last = 0.0
    while last < a-1:
        out.append([int(last), int(last+avg)])
        last += avg
    return out

def get_blockerror(Data:np.array):
    data = Data
    average = np.average(data)
    be = block(data)**.5

    return average, be


def get_blockerror_pyblock(Data:np.array):
    average = np.average(Data)
    if (average != 0) and (average != 1):
        reblock_data = pyblock.blocking.reblock(Data)
        opt = pyblock.blocking.find_optimal_block(len(Data), reblock_data)[0]
        be = reblock_data[opt][4]
    else:
        be = 0

    return average, float(be)


def get_blockerror_pyblock_nanskip_rw_(Data:np.array,weights:np.array=None):
    if weights is not None : average = np.dot(Data,weights)
    else : average = np.average(Data)

    if (average != 0) and (average != 1):
        reblock_data = pyblock.blocking.reblock(Data,weights=weights)
        opt = pyblock.blocking.find_optimal_block(len(Data), reblock_data)[0]
        if(math.isnan(opt)):
            be_max = 0
            for i in range(0, len(reblock_data)):
                be = reblock_data[i][4]
                if(be > be_max):
                    be_max = be
        else:
            be = reblock_data[opt][4]
    else:
        be = 0

    return average, float(be)

def get_blockerror_pyblock_max(Data:np.array):
    average=np.average(Data)
    if (average!=0) and (average!=1):
        reblock_data = pyblock.blocking.reblock(Data)
        be_max=0
        for i in range(0,len(reblock_data)): 
            be=reblock_data[i][4]
            if(be > be_max):
                be_max=be
    else:
        be=0

    return average,float(be)

def get_blockerrors(Data:np.array, bound_frac:float):
    n_data = len(Data[0])
    block_errors = []
    ave = []
    for i in range(0, n_data):
        data = Data[:, i]
        average = np.average(data)
        be = block(data)**.5
        ave.append(np.average(data))
        block_errors.append(be)
    ave_bf = np.asarray(ave)/bound_frac
    be_bf = np.asarray(block_errors)/bound_frac

    return ave_bf, be_bf


def get_blockerrors_pyblock(Data:np.array, bound_frac:float):
    n_data = len(Data[0])
    block_errors = []
    ave = []
    for i in range(0, n_data):
        data = Data[:, i]
        average = np.average(data)
        if (average != 0) and (average != 1):
            reblock_data = pyblock.blocking.reblock(data)
            opt = pyblock.blocking.find_optimal_block(
                len(data), reblock_data)[0]
            opt_block = reblock_data[opt]
            be = opt_block[4]
        else:
            be = 0
        ave.append(average)
        block_errors.append(be)

    ave_bf = np.asarray(ave)/bound_frac
    be_bf = np.asarray(block_errors)/bound_frac

    return ave_bf, be_bf


def get_blockerrors_pyblock_nanskip_rw_(Data:np.array, bound_frac:float, Weights:np.array=[]):
    n_data = len(Data[0])
    block_errors = []
    ave = []
    for i in range(0, n_data):
        data = Data[:, i]
        weights=Weights
        if len(weights)>0 : average = np.dot(data,weights)
        else : average = np.average(Data)
        
        if (average != 0) and (average != 1):
            reblock_data = pyblock.blocking.reblock(data,weights=weights)
            opt = pyblock.blocking.find_optimal_block(len(data), reblock_data)[0]
#          opt_block = reblock_data[opt]
#          be = opt_block[4]
            if (math.isnan(opt)):
                be_max=0
                for i in range(0, len(reblock_data)):
                    be = reblock_data[i][4]
                    if (be > be_max) :
                        be_max=be
            else:
                be = reblock_data[opt][4]
            
        else:
            be = 0
        ave.append(average)
        block_errors.append(be)

    ave_bf = np.asarray(ave)/bound_frac
    be_bf = np.asarray(block_errors)/bound_frac
    return ave_bf, be_bf

def compute_temperatures(temp_range:tuple, nreps:int):
    import numpy as np
    from math import exp, log
    tlow, thigh = temp_range
    temps = []
    for i in range(nreps):
        temps.append(np.round(tlow*exp((i)*log(thigh/tlow)/(nreps-1)),3))

    return np.array(temps)

def sequence_ticks_1(sequence):

    aa_dict = {'ALA' : ["A", "Alanine"], 'ARG' : ["R", "Arginine"], 'ASN' : ["N", "Asparagine"], 'ASP' : ["D", "Aspartic-acid"], 'CYS' : ["C", "Cysteine"], 
            'GLU' : ["E", "Glutamic-acid"], 'GLN' : ["Q", "Glutamine"], 'GLY' : ["G", "Glycine"], 'HIS' : ["H", "Histidine"], 'ILE' : ["I", "Isoleucine"],
            'LEU' : ["L", "Leucine"], 'LYS' : ["K", "Lysine"], 'MET' : ["M", "Methionine"], 'PHE' : ["F", "Phenylalanine"], 'PRO' : ["P", "Proline"], 
            'SER' : ["S", "Serine"], 'THR' : ["T", "Threonine"], 'TRP' : ["W", "Tryptophan"], 'TYR' : ["Y", "Tyrosine"], 'VAL' : ["V", "Valine"],
            'CALA' : ["A", "Alanine"], 'NASP' : ["D", "Aspartic-acid"] }

    import numpy as np

    def split_temp(inp:np.ndarray):
        import re

        exp=r"([a-z]+)([0-9]+)"
        out=[]

        for i in range(len(inp)):
            match = re.match(exp, inp[i], re.I)

            if match:
                items = match.groups()

            out.append(list(items))

        return np.array(out)

    a=split_temp(np.array(sequence))

    b=np.zeros(len(a),dtype=object)
    c=np.zeros(len(a),dtype=object)

    for i in range(0,len(a)):

        for k in aa_dict.keys():

            if k in a[i]:
                
                a[i][0] = aa_dict[k][0]

        b[i]=''.join(a[i])
        c[i]=a[i][0]

    return np.array(b) , np.array(c)

def make_dir(dir_name:str):
    import os
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)