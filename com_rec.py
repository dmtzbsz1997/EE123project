# Import functions and libraries
import numpy as np
import matplotlib.pyplot as plt
from imageio import imread
import time
#import pywt
import PIL
# from numpy.lib.stride_tricks import as_strided
import gzip
import os
import shutil
from scipy import interpolate
import scipy.misc
import imageio


import numpy.ctypeslib as npct
from ctypes import c_int
from ctypes import c_float


import queue as Queue
import time
import sys

from numpy import pi
from numpy import sin
from numpy import zeros
from numpy import r_
from numpy import ones
from scipy import signal
from scipy import integrate
import threading

from numpy import mean
from numpy import power
from numpy.fft import fft
from numpy.fft import fftshift
from numpy.fft import ifft
from numpy.fft import ifftshift
import bitarray
from  scipy.io.wavfile import read as wavread
import newax25 as ax25

import multiprocessing

from math import gcd
import sounddevice as sd
import RPi.GPIO as GPIO
from functools import reduce
from numpy import ones,zeros, pi, cos, exp, sign

import numpy as np

from PIL import Image
import time
import glob
import os
import bokeh.plotting as bk
from bokeh.io import push_notebook
from bokeh.resources import INLINE
from bokeh.models import GlyphRenderer

bk.output_notebook(INLINE)

import re

numbers = re.compile(r'(\d+)')



array_1d_int = npct.ndpointer(dtype=np.int, ndim=1, flags='CONTIGUOUS')

libcd = npct.load_library("./libpll", ".")
libcd.pll.restype = c_int
libcd.pll.argtypes= [array_1d_int, c_int, array_1d_int,array_1d_int,  array_1d_int,array_1d_int, c_int, c_float]


import numpy.ctypeslib as npct
from ctypes import c_int
from ctypes import c_float

array_1d_int = npct.ndpointer(dtype=np.int, ndim=1, flags='CONTIGUOUS')

libcd = npct.load_library("./libpll", ".")
libcd.pll.restype = c_int
libcd.pll.argtypes= [array_1d_int, c_int, array_1d_int,array_1d_int,  array_1d_int,array_1d_int, c_int, c_float]

# function to compute least common multipler
def lcm(numbers):
    return reduce(lambda x, y: (x*y)//gcd(x,y), numbers, 1)

class TNCaprs:
    
    def __init__(self, fs = 48000.0, Abuffer = 512, Nchunks=10,baud=1200,center=3500):
        
        #  Implementation of an afsk1200 TNC. 
        #
        #  The TNC processes a `Abuffer` long buffers, till `Nchunks` number of buffers are collected into a large one.
        #  This is because python is able to more efficiently process larger buffers than smaller ones.
        #  Then, the resulting large buffer is demodulated, sampled and packets extracted.
        #
        # Inputs:
        #    fs  - sampling rate
        #   TBW  -  TBW of the demodulator filters
        #   Abuffer - Input audio buffers from Pyaudio
        #   Nchunks - Number of audio buffers to collect before processing
        #   plla    - agressivness parameter of the PLL
        
        
        ## compute sizes based on inputs
        self.center_f = center
        self.space_f = self.center_f + 500
        self.mark_f = self.center_f - 500
        self.baud = baud 
        self.TBW = 2.0   # TBW for the demod filters
        self.N = (int(fs/self.baud*self.TBW)//2)*2+1   # length of the mark-space filters for demod
        self.fs = fs     # sampling rate   
        self.BW = 1200      # BW of filter based on TBW
        self.Abuffer = Abuffer             # size of audio buffer
        self.Nchunks = Nchunks             # number of audio buffers to collect
        self.Nbuffer = Abuffer*Nchunks+(self.N*3-3)         # length of the large buffer for processing
         
        self.Ns = 1.0*fs/self.baud # samples per symbol
        ## state variables for the modulator
        self.prev_ph = 0  # previous phase to maintain continuous phase when recalling the function
                         
        ##  Generate Filters for the demodulator
        self.h_lp = signal.firwin(self.N,self.BW/self.fs*1.0,window='hanning')
        self.h_lpp = signal.firwin(self.N,self.BW*2*1.2/self.fs,window='hanning')
       
        self.h_space = self.h_lp*exp(1j*2*pi*(self.space_f)*r_[-self.N/2:self.N/2]/self.fs)
        self.h_mark = self.h_lp*exp(1j*2*pi*(self.mark_f)*r_[-self.N/2:self.N/2]/self.fs)
        self.h_bp = (signal.firwin(self.N,self.BW/self.fs*2.2,window='hanning'))*exp(1j*2*pi*self.center_f*r_[-self.N/2:self.N/2]/self.fs)
        

        

        ## PLL state variables  -- so conntinuity between buffers is preserved
        self.dpll = np.round(2.0**32 / self.Ns).astype(np.int32)    # PLL step
        self.pll =  0                # PLL counter
        self.ppll = -self.dpll       # PLL counter previous value -- to detect overflow
        self.plla = 0.74             # PLL agressivness (small more agressive)
        

        ## state variable to NRZI2NRZ
        self.NRZIprevBit = bool(1)  
        
        ## State variables for findPackets
        self.state='search'   # state variable:  'search' or 'pkt'
        self.pktcounter = 0   # counts the length of a packet
        self.packet = bitarray.bitarray([0,1,1,1,1,1,1,0])   # current packet being collected
        self.bitpointer = 0   # poiter to advance the search beyond what was already searched in the previous buffer

        ## State variables for processBuffer
        self.buff = zeros(self.Nbuffer)   # large overlapp-save buffer
        self.chunk_count = 0              # chunk counter
        self.oldbits = bitarray.bitarray([0,0,0,0,0,0,0])    # bits from end of prev buffer to be copied to beginning of new
        self.Npackets = 0                 # packet counter
        
        
    
    
    def NRZ2NRZI(self,NRZ, prevBit = True):
        NRZI = NRZ.copy() 
        for n in range(0,len(NRZ)):
            if NRZ[n] :
                NRZI[n] = prevBit
            else:
                NRZI[n] = not(prevBit)
            prevBit = NRZI[n]
        return NRZI
    



    def NRZI2NRZ(self, NRZI):  
        NRZ = NRZI.copy() 
    
        for n in range(0,len(NRZI)):
            NRZ[n] = NRZI[n] == self.NRZIprevBit
            self.NRZIprevBit = NRZI[n]
    
        return NRZ
    
    def KISS2bits(self,KISS):
        # function that takes a KISS frame sent via TCP/IP and converts it to an APRSpacket bit stream.
        
        bits = bitarray.bitarray(endian="little")
        bits.frombytes(KISS)
        fcs = ax25.FCS()
        for bit in bits:
            fcs.update_bit(bit)
            
        bits.frombytes(fcs.digest())
        return bitarray.bitarray('01111110') + ax25.bit_stuff(bits) + bitarray.bitarray('01111110') 
     
    def bits2KISS(self,bits):
        # function that takes a bitstream of an APRS-packet, removes flags and FCS and unstuffs the bits
        bitsu = ax25.bit_unstuff(bits[8:-8])
        return  bitsu[:-16].tobytes() 
    
    
    def modulate(self,bits):
    # the function will take a bitarray of bits and will output an AFSK1200 modulated signal of them, sampled at fs Hz
    #  Inputs:
    #         bits  - bitarray of bits
    #         fs    - sampling rate
    # Outputs:
    #         sig    -  returns afsk1200 modulated signal 
        # For you to complete
        fss = lcm((self.baud,self.fs))
        deci = fss//self.fs
        Nb = fss//self.baud
        print("Nb:",Nb)
        nb = len(bits)
        NRZ = ones((nb,Nb))
        for n in range(0,nb):
            if bits[n]:
                NRZ[n,:]=-NRZ[n,:]
    
#         freq = 1700 + 500*NRZ.ravel()
        freq = self.center_f + 500*NRZ.ravel()
        ph = 2.0*pi*integrate.cumtrapz(freq)/fss
        sig = cos(ph[::deci])
        return sig
        

    
    def modulatPacket(self, callsign, digi, dest, info, preflags=80, postflags=80 ):
        
        # given callsign, digipath, dest, info, number of pre-flags and post-flags the function contructs
        # an appropriate aprs packet, then converts them to NRZI and calls `modulate` to afsk 1200 modulate the packet. 
        
        packet = ax25.UI(destination=dest,source=callsign, info=info, digipeaters=digi.split(b','),)
        prefix = bitarray.bitarray(np.tile([0,1,1,1,1,1,1,0],(preflags,)).tolist())
        suffix = bitarray.bitarray(np.tile([0,1,1,1,1,1,1,0],(postflags,)).tolist())
        sig = self.modulate(self.NRZ2NRZI(prefix + packet.unparse()+suffix))

        return sig
    
    

    def demod(self, buff):
        #Demodulates a buffer and returns valid NRZ
        # Similar to afsk1200_demod,  for you to complete

        
        buff = np.convolve(buff.copy(),self.h_bp,'valid')
        mark = abs(np.convolve(buff,self.h_mark,mode='valid'))
        space = abs(np.convolve(buff,self.h_space,mode='valid'))
        NRZ = mark-space
        NRZ = np.convolve(NRZ,self.h_lpp,'valid')
        return NRZ



    def FastPLL(self,NRZa):
        recbits = np.zeros(len(NRZa)//(self.fs//self.baud)*2,dtype=np.int32)
        pll = np.zeros(1,dtype = np.int32)
        pll[0] = self.pll
        ppll = np.zeros(1,dtype = np.int32)
        ppll[0] = self.ppll
        
        #print("pll = ",pll,"   ppll=",ppll)
        
        
        NRZb = (NRZa > 0).astype(np.int32)
        tot = libcd.pll(NRZb,len(NRZb),recbits,recbits,pll,ppll,self.dpll,self.plla)
        
        self.ppll = ppll.copy()
        self.pll = pll.copy()
        
        #print("post: pll = ",pll,"   ppll=",ppll)
        
        return bitarray.bitarray(recbits[:tot].tolist())
    
    def PLL(self, NRZa):
       #print("running PLL")
        idx = zeros(len(NRZa)//int(self.Ns)*2)   # allocate space to save indexes        
        c = 0
        
        for n in range(1,len(NRZa)):
            if (self.pll < 0) and (self.ppll >0):
                idx[c] = n
                c = c+1
        
            if (NRZa[n] >= 0) !=  (NRZa[n-1] >=0):
                self.pll = np.int32(self.pll*self.plla)
            
        
            self.ppll = self.pll
            self.pll = np.int32(self.pll+ self.dpll)
    
        return idx[:c].astype(np.int32) 
    
   

    def findPackets(self,bits):
        # function take a bitarray and looks for AX.25 packets in it. 
        # It implements a 2-state machine of searching for flag or collecting packets
        flg = bitarray.bitarray([0,1,1,1,1,1,1,0])
        packets = []
        n = self.bitpointer
        
        # Loop over bits
        while (n < len(bits)-7) :
            # default state is searching for packets
            if self.state is 'search':
                # look for 1111110, because can't be sure if the first zero is decoded
                # well if the packet is not padded.
                if bits[n:n+7] == flg[1:]:
                    # flag detected, so switch state to collecting bits in a packet
                    # start by copying the flag to the packet
                    # start counter to count the number of bits in the packet
                    self.state = 'pkt'
                    self.packet=flg.copy()
                    self.pktcounter = 8
                    # Advance to the end of the flag
                    n = n + 7
                else:
                    # flag was not found, advance by 1
                    n = n + 1            
        
            # state is to collect packet data. 
            elif self.state is 'pkt':
                # Check if we reached a flag by comparing with 0111111
                # 6 times ones is not allowed in a packet, hence it must be a flag (if there's no error)
                if bits[n:n+7] == flg[:7]:
                    # Flag detected, check if packet is longer than some minimum
                    if self.pktcounter > 200:
                        #print('packet found!')
                        # End of packet reached! append packet to list and switch to searching state
                        # We don't advance pointer since this our packet might have been
                        # flase detection and this flag could be the beginning of a real packet
                        self.state = 'search'
                        self.packet.extend(flg)
                        packets.append(self.packet.copy())
                    else:
                        # packet is too short! false alarm. Keep searching 
                        # We don't advance pointer since this this flag could be the beginning of a real packet
                        self.state = 'search'
                # No flag, so collect the bit and add to the packet
                else:
                    # check if packet is too long... if so, must be false alarm
                    if self.pktcounter < 2680:
                        # Not a false alarm, collect the bit and advance pointer        
                        self.packet.append(bits[n])
                        self.pktcounter = self.pktcounter + 1
                        n = n + 1
                    else:  #runaway packet
                        #runaway packet, switch state to searching, and advance pointer
                        self.state = 'search'
                        n = n + 1
        
        self.bitpointer = n-(len(bits)-7) 
        return packets

    
    # function to generate a checksum for validating packets
    def genfcs(self,bits):
        # Generates a checksum from packet bits
        fcs = ax25.FCS()
        for bit in bits:
            fcs.update_bit(bit)
    
        digest = bitarray.bitarray(endian="little")
        digest.frombytes(fcs.digest())

        return digest




    # function to parse packet bits to information
    def decodeAX25(self,bits, deepsearch=False):
        ax = ax25.AX25()
        ax.info = "bad packet"
    
    
        bitsu = ax25.bit_unstuff(bits[8:-8])
    
        
        #foundPacket = False
        #if (self.genfcs(bitsu[:-16]).tobytes() == bitsu[-16:].tobytes()):
        #        foundPacket = True
        #elif deepsearch: 
        #    tbits = bits[8:-8]
        #    for n in range(0,len(tbits)):
        #        tbits[n] = not tbits[n]
        #        if (self.genfcs(bitsu[:-16]).tobytes() == bitsu[-16:].tobytes()):
        #            foundPacket = True
        #            print("Success deep search")
        #            break
        #        tbits[n] = not tbits[n]
        # 
        #if foundPacket == False:
        #    return ax
        
        if (self.genfcs(bitsu[:-16]).tobytes() == bitsu[-16:].tobytes()) == False:
            #print("failed fcs")
            return ax
                  
    
        bytes = bitsu.tobytes()
        ax.destination = ax.callsign_decode(bitsu[:56]).decode('ascii')
        source = ax.callsign_decode(bitsu[56:112]).decode('ascii')
    
        if source[-1].isdigit() and source[-1]!="0":
            ax.source = "".join((source[:-1],'-',source[-1]))
        else:
            ax.source = source[:-1]
    
        digilen=0    
    
        if bytes[14]=='\x03' and bytes[15]=='\xf0':
            digilen = 0
        else:
            for n in range(14,len(bytes)-1):
                if bytes[n] & 1:
                    digilen = (n-14)+1
                    break

        #    if digilen > 56:
        #        return ax
        ax.digipeaters =  ax.callsign_decode(bitsu[112:112+digilen*8]).decode('ascii')
        ax.info = bitsu[112+digilen*8+16:-16].tobytes()
    
    
        return ax

    def processBuffer(self, buff_in):
        
        # function processes an audio buffer. It collect several small into a large one
        # Then it demodulates and finds packets.
        #
        # The function operates as overlapp and save
        # The function returns packets when they become available. Otherwise, returns empty list
        
        N = self.N
        NN = (N*3 -3 )
        
        
        Nchunks = self.Nchunks
        Abuffer = self.Abuffer
        fs = self.fs
        Ns = self.Ns
        
        validPackets=[]
        packets=[]
        NRZI=[]
        idx = []
        bits = []
        
        # Fill in buffer at the right place
        self.buff[NN+self.chunk_count*Abuffer:NN+(self.chunk_count+1)*Abuffer] = buff_in.copy()
        self.chunk_count = self.chunk_count + 1
        
        
        # number of chunk reached -- process large buffer
        if self.chunk_count == Nchunks:
            # Demodulate to get NRZI
            NRZI = self.demod(self.buff)
            # compute sampling points, using PLL
            #idx = self.PLL(NRZI)
            # Sample and make a decision based on threshold
            #bits = bitarray.bitarray((NRZI[idx]>0).tolist())
            
            bits = self.FastPLL(NRZI)
            #bits = self.PLL(NRZI)
            # In case that buffer is too small raise an error -- must have at least 7 bits worth
            if len(bits) < 7:
                raise ValueError('number of bits too small for buffer')
            
            # concatenate end of previous buffer to current one
            bits = self.oldbits + self.NRZI2NRZ(bits)
            
            # store end of bit buffer to next buffer
            self.oldbits = bits[-7:].copy()
            
            # look for packets
            packets = self.findPackets(bits)
            
            # Copy end of sample buffer to the beginning of the next (overlapp and save)
            self.buff[:NN] = self.buff[-NN:].copy()
            
            # reset chunk counter
            self.chunk_count = 0
            
            # checksum test for all detected packets
            for n in range(0,len(packets)):
                if len(packets[n]) > 200: 
                    try:
                        ax = self.decodeAX25(packets[n])
                    except:
                        ax = ax25.AX25()
                        ax.info = "bad packet"
                    if ax.info != 'bad packet':
                        validPackets.append(packets[n])
                        
            
        return validPackets


def file_detection(img,m):
    #k compression rate of PCA : Frame/k
    #m compression rate of Huffman coding
    #return d downsample rate
#     Frame, Height, Width, Color = img.shape
#     k = int(Frame/3)
#     d = int(np.sqrt(Frame*Height*Width*Color/(256*35)/m/3)) 
#     return k,d
    Frame, Height, Width, Color = img.shape
    PSNR = np.zeros((Frame, 10))
    d_final=3
    for k in range(2,Frame+1):
        for d in range(d_final,11):
            print(k,d)
            data = compress(img,downrate=d, rank=k)
            img_rec = recon(data)
            np.save('simple1.npy',data)
            f = gzip.GzipFile('simple1.gz','w')
            np.save(f, np.load('simple1.npy'))
            f.close()
            if os.path.getsize('simple1.gz') < 9000:
                PSNR[k-1,d-1] = psnr(img, img_rec, maxVal=255)
                d_final=d
                break
                
    k,d = np.where(PSNR == np.max(PSNR))
    k += 1
    d += 1
    return int(np.array(k)),int(np.array(d))



def compress(img, downrate = 10, rank = 20):
    
    ########### uint8 ###########
    # rank     determine from detection function
    #downrate  upsample rate, also determine from detection function
    #dim     original dimension (flame,heigh,width,color)
    #un,uv   negative number index
    
    ########### float ###########
    #s       sigular value
    #mu,mv   max of u and v for normalization
    
    def reshape1(img):
        shape = img.shape
        img = img.reshape((shape[0],np.prod(shape[1:])))
        return img
    def downsample(img,n):
        img = img[::,::n,::n]
        return img,img.shape
    odim = img.shape[1:3]
    img, dim = downsample(img,downrate)
    img = reshape1(img)
    u,s,v = np.linalg.svd(img,full_matrices=False)
    ud = u[:,:rank]
    vd = v[:rank,:]
#     un = np.nonzero(ud<0)
#     vn = np.nonzero(vd<0)
    umin = np.min(ud)
    vmin = np.min(vd)
    ud = ud - umin
    vd = vd - vmin

    
#     umax = np.max(ud)
#     vmax = np.max(vd)
    
#     ud = abs(ud)
#     vd = abs(vd)
    s = s[:rank]
    ku = 255/np.max(ud)
    kv = 255/np.max(vd)
    d = [ud*ku,vd*kv]
#     midu = np.unravel_index(np.argmax(d[0], axis=None), d[0].shape)
#     midv = np.unravel_index(np.argmax(d[1], axis=None), d[1].shape)
#     d[0][un] = (d[0][un] + 1)
#     d[1][vn] = (d[1][vn] + 1)
#     d[0][midu] = 255
#     d[1][midv] = 255
    d = [i.astype('uint8') for i in d]
#     d[1] = averager(d[1],3)
    d[1] = averager(d[1], 5, 20)
#     print(d)

#     nid = [un,vn] 
    return [s,[ku,kv],[umin,vmin]], d, [dim,odim]


def recon(data, uprate = 5):
    #mu,mv   max of u and v for normalization
    #un,uv   negative number index
    #s       sigular value
    #dim     original dimension (flame,heigh,width,color)
    #uprate  upsample rate
    s = data[0][0]
    ku = data[0][1][0]
    kv = data[0][1][1]   
    um = data[0][2][0]
    vm = data[0][2][1]
    d = data[1]
    dim = data[2][0]
    odim = data[2][1]



    def upsample(img,odim):
#         lanczos:16.54804770571443, 23.220852341612467
#           bilinear  16.670716307372775, 23.548733868332896, 17.037282941006527
#         bicubic 16.60428051658223, 23.332048513937874
#        cubic 16.60428051658223
#      nearest 16.10182141434048
        shape = img.shape
        if shape[3] == 4:
            return  [scipy.misc.imresize(i, odim, interp = 'bicubic', mode = 'RGBA') for i in img]
        elif shape[3] == 3:
            return  [scipy.misc.imresize(i, odim, interp = 'bicubic', mode = 'RGB') for i in img]
#     [np.array(Image.fromarray(i, mode = 'RGB').resize(odim, PIL.Image.BICUBIC)) for i in img]
#         return img.repeat(n, axis=1).repeat(n, axis=2)
    def reshape2(img,dim):
        img = img.reshape((dim))
        return img      
    d = [d[0]/ku,d[1]/kv]
#     d[0][un] = -1*(d[0][un])
#     d[1][vn] = -1*(d[1][vn])
    d[0] = d[0] + um
    d[1] = d[1] + vm
    x = np.matmul(np.matmul(d[0],np.diag(s)),d[1]) 
    x = x/np.max(x)*255
    x = reshape2(x,dim)
#     print(x.shape)
    im = upsample(x,odim)
#     print(im.shape)
#     im = np.array(abs(im)).astype(int)
    return np.array(im) #im

def recon_display(filename,img, imgor):
    path = filename +'_out/'
    os.system("rm -rf {:s}".format(path))
    os.system("mkdir {:s}".format(path))
    for i in range(len(img)):
        imageio.imwrite(path+'img{}.tiff'.format(i+1),img[i])
    print('PSNR: ',psnr(img, imgor))
    Tiff_play(path= path+'*.tiff', display_size=400, frame_rate=15)
        
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts


# Tiff stack player

def Tiff_play(path, display_size, frame_rate):
    image_files = sorted(glob.glob(path), key=numericalSort)
    Nframe = len(image_files)

    im = Image.open(image_files[0])
    xdim, ydim= im.size
    display_array = np.zeros((Nframe,ydim,xdim,4),dtype='uint8')
    
    # load image stack
    for i in range(0,Nframe):
        im = Image.open(image_files[i])
        im = im.convert("RGBA")
        imarray = np.array(im)
        display_array[i] = np.flipud(imarray)
        
    # Play video

    wait_time = 1/frame_rate
    normalized_size = display_size
    max_size = np.maximum(xdim,ydim)
    width = (xdim/max_size * normalized_size).astype('int')
    height = (ydim/max_size * normalized_size).astype('int')
    
    counter = 0
    first_round = True
    try:
        while True:
            if counter == 0 and first_round:
                p = bk.figure(x_range=(0,xdim), y_range=(0,ydim), plot_height = height, plot_width = width)
                p.image_rgba(image=[display_array[counter]], x=0, y=0, dw=xdim, dh=ydim, name='video')
                bk.show(p, notebook_handle=True)
                counter += 1
                first_round = False
            else:
                renderer = p.select(dict(name='video', type=GlyphRenderer))
                source = renderer[0].data_source
                source.data['image'] = [display_array[counter]]
                push_notebook()
                if counter == Nframe-1:
                    counter = 0
                else:
                    counter += 1
            time.sleep(wait_time)
            
    except KeyboardInterrupt:
        pass
            
# Image_stack loader

def Tiff_load(path):
    image_files = sorted(glob.glob(path), key=numericalSort)
    Nframe = len(image_files)
    
    im = Image.open(image_files[0])
    xdim, ydim= im.size
    image_stack = np.zeros((Nframe,ydim,xdim,3),dtype='uint8')
    
    for i in range(0,Nframe):
        im = Image.open(image_files[i])
        image_stack[i] = np.array(im)
    
    return image_stack

# Image_stack loader with ffmpeg

def imageStack_load(filename):
    path        = filename[:filename.find('.')]+'/'
    os.system("rm -rf {:s}".format(path))
    os.system("mkdir {:s}".format(path))
    os.system("ffmpeg -i {:s} {:s}frame_%2d.tiff".format(filename, path))
    image_files = sorted(glob.glob(path+"*.tiff"), key=numericalSort)
    Nframe = len(image_files)
    
    im = Image.open(image_files[0])
    xdim, ydim= im.size
    image_stack = np.zeros((Nframe,ydim,xdim,3),dtype='uint8')
    
    for i in range(0,Nframe):
        im = Image.open(image_files[i])
        image_stack[i] = np.array(im)
    
    return image_stack

# Save gif with ffmpeg

def GIF_save(path, framerate):
    os.system("ffmpeg -r {:d} -i {:s}frame_%2d.tiff -compression_level 0 -plays 0 -f apng {:s}animation.png".format(framerate, path,path))

# Compute video PSNR

def psnr(ref, meas, maxVal=255):
    assert np.shape(ref) == np.shape(meas), "Test video must match measured vidoe dimensions"


    dif = (ref.astype(float)-meas.astype(float)).ravel()
    mse = np.linalg.norm(dif)**2/np.prod(np.shape(ref))
    psnr = 10*np.log10(maxVal**2.0/mse)
    return psnr
    

     
def averager(np_array, num_bins, threshold, ignore_first = False):
    #np_array is the large PCA array
    #num_bins is the number of bins for the histogram
    #threshold is how far away from the averag you would be before you are sent to it
    for vector in np_array[1:]:
        hist, bin_edges = np.histogram(vector, bins = num_bins, range = (0, 255));
        indexes = {}
        for idx, bin_edge in enumerate(bin_edges):
            indexes[idx] = []
        for v_idx, value in enumerate(vector):
            for idx, bin_edge in enumerate(bin_edges):
                if idx < num_bins:
                    if value >= bin_edge and value <= bin_edges[idx + 1]:
                        indexes[idx] += [v_idx]
        for key in indexes:
            new_vec = vector[indexes[key]]
            if len(new_vec) > 0:
                avg = np.average(new_vec)
                for idx in indexes[key]:
                    if vector[idx] > avg - threshold and vector[idx] < avg + threshold:
                        vector[idx] = avg
#         plt.plot(vector)
#         plt.show()
            
        
    return np_array
            
def zip_up(data, filename):
    f = gzip.GzipFile(filename + '.gz', "w")
    np.save(f, data)
    f.close()

def unzip(pathname):
    f = gzip.GzipFile(pathname, "r")
    data = np.load(pathname)
    f.close()
    return data
