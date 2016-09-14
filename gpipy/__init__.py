# Licensed under a 3-clause BSD style license - see LICENSE.rst

""" Gemini Planet Imager Python data reduction tools (GPIPy)

Or should this be GPI Python Infrastructure (GPIPI?) perhaps.

"""


#import data
from .data import read, IFSData, IFSSpectralCube,IFSPolarimetryPair,IFSStokesCube,DataCollection
from .pipeline import *


# Only try to import the GUI classes if wx is present.
try:
    import wx as _wx
    _HAVE_WX = True
except:
    _HAVE_WX = False

if _HAVE_WX:
    import guis
