#!/usr/bin/env python
#import os, sys
import numpy as np
import matplotlib.pyplot as pl
import matplotlib
import os
import astropy
try:
    from IPython.core.debugger import Tracer; stop = Tracer()
except:
    pass

try:
    import astropy.io.fits as fits
except:
    import pyfits as fits

import logging
import glob
import copy
_log = logging.getLogger('gpidata')

__version__ = '0.1'


__doc__ ="""
GPI Data File interface


Various Pythonic, object interfaces to GPI data files

Simple Usage: 
    >>> dat = gpidata.read('somefilename.fits')

The type of object that is returned will vary depending on what kind of file you are loading; there are subclasses for
2D and 3D GPI files, polarimetry data, wavelength calibrations. What you can do with a gpidata object depends on the subclass. 
At a minimum you can access *the most commonly used* keywords as attributes: 

    >>> dat.filter
    'H'
    etc.


There is also a DataCollection class that lets you manipulate and query large groups of FITS files at once:

    >>> dc = gpidata.DataCollection('S20130909S*.fits')
    >>> dc.itimes
    [numpy ndarray of all the itimes for all those files]


The DataCollection.ParseToRecipes()  function replicates much of 
the functionality of the GPI Data Pipeline's Data Parser. It will generate 
recipes to reduce your data. 

"""


def read(filename, loadData=True, verbose=False, **kwargs):
    """ Open a GPI data file and return it as an object.
    
    This is a factory function that returns an object 
    of the most specific appropriate type (e.g. IFSSpectralCube, IFSPolarimetryCube etc)
    depending on the header keywords found therein.

    Examples
    ---------
    >>> df = gpidata.read('/path/to/myfile.fits')
    >>> df
    <instance of gpidata.IFSSpectralCube>

    
    Parameters
    ----------
    loadData : bool
        if True (default), upon opening loads into memory the data associated with this file.
        If False, only the headers are loaded.

    verbose : bool
        Print verbose output about file type determination. 
    

    """
    head = fits.getheader(filename,ext=0)

    try:
        exthead = fits.getheader(filename,ext=1)
    except:
        exthead=head # if we only have one header (in early devel. GPI data files) just use it for both sets of keywords

    try:
        if exthead['NAXIS'] == 3 and head['FILETYPE'] == 'Wavelength Solution Cal File':
            if verbose: _log.info("File %s appears to be a Wavelength Solution Cal File" % filename)
            return IFSwavecal(filename, loadData=loadData, **kwargs)
    except: pass

    try: 
        if exthead['NAXIS'] == 3 and 'PRISM' in head['DISPERSR']:
            if verbose: _log.info("File %s appears to be a Spectral Data Cube" % filename)
            return IFSSpectralCube(filename, loadData=loadData, **kwargs)
    except: pass

    if exthead['NAXIS'] == 3 and 'WOLL' in head['DISPERSR']:
        if verbose: _log.info("File %s appears to be a Polarimetric Data Cube" % filename)
        return IFSPolarimetryCube(filename, loadData=loadData, **kwargs)
    elif exthead['NAXIS'] == 2 and exthead['XTENSION'] == 'IMAGE':
        if verbose: _log.info("File %s appears to be a 2D IFS Image" % filename)
        return IFSData2D(filename, loadData=loadData, **kwargs)
    elif exthead['NAXIS'] == 2 and exthead['XTENSION'] == 'BINTABLE' and head['FILETYPE']=='1D Spectrum':
        if verbose: _log.info("File %s appears to be a 2D IFS Image" % filename)
        return IFSSpectrum1D(filename, loadData=loadData, **kwargs)
    else:
        if verbose: _log.info("Cannot determine file type for %s; returning as generic data" % filename)
        return GPI_generic_data(filename, loadData=loadData, **kwargs)


#-------------------------------------------------------------------------------

class GPI_generic_data(object):
    """ Object class for GPI IFS data

    Provides a simple interface to manipulate files

    Attributes:
        .data
        .var
        .dq
        .priheader
        .extheader

        .filter
        .prism
        .occulter
        .apodizer


    Parameters
    ------------
    loadData : bool
        If True, load data of the file. If False, just load the headers.
        Default is True. Note that if you select False then later attempt
        to access the data attribute, the data will be loaded at that time.
        Thus you can lazily access data while minimizing memory overhead. 


    """
    def __init__(self, filename, loadData=True):
        self.name = os.path.basename(filename)
        self.abspath = os.path.abspath(filename)
        "Filename or other description"
        if loadData:
            self._loadData()
        else:
            # just load headers
            ph = fits.PrimaryHDU(header=fits.getheader(filename,ext=0))
            try:
                ih = fits.ImageHDU(header=fits.getheader(filename,ext=1))
            except:
                ih = ph # if we only have one header (in early devel. GPI data files) just use it for both sets of keywords
            self._HDUlist=fits.HDUList([ph,ih])
            self._data_loaded=False


        if len(self._HDUlist) > 1:
            self._sci_ext = 'SCI'
        else:
            self._sci_ext = 0
 
    def _loadData(self):
        """ Load the data in a file """
        self._HDUlist= fits.open(self.abspath)
        self._data_loaded=True


    def __repr__(self):
        return "<%s.%s : %s>" % (self.__class__.__module__, self.__class__.__name__, self.name) 
    ### Data properties ###
    @property
    def data(self):
        "Data in the file"
        if self._data_loaded is False: self._loadData()
        return self._HDUlist[self._sci_ext].data
    @data.setter
    def data(self, value):
        self._HDUlist[self._sci_ext].data = value

    @property
    def var(self):
        return self._HDUlist['VAR'].data
    @property
    def dq(self):
        return self._HDUlist['DQ'].data


    @property
    def shape(self):
        return self.data.shape
    @shape.setter
    def shape(self):
        self.data.shape = value
 
    @property
    def priheader(self):
        return self._HDUlist[0].header
    @property
    def extheader(self):
        return self._HDUlist[self._sci_ext].header
    def Copy(self):
        "Return a copy of the data file as a different object."
        return copy.deepcopy(self)
    def __getitem__(self, key):
        """ Cleverly allow access to either data or keywords depending
        on whether the key is a string or not.
        """
        if isinstance(key, basestring):
            try:
                return self.priheader[key]
            except:
                return self.extheader[key]
        else:
            return self.data[key]
    ### Header properties ###

    def _extKeyword(self, keyword, docstr=''):
        #def actualfn(self, keyword):
            docstr
            try: 
                return self.extheader[keyword]
            except:
                return None
        #return actualfn
    def _priKeyword(self, keyword, docstr=""):
        #def actualfn(self, keyword):
            docstr
            try: 
                return self.priheader[keyword]
            except:
                return None
        #return actualfn


    @property
    def filter(self):
        val = self._priKeyword('IFSFILT', 'Filter name')
        try:
            if val.startswith('IFSFILT_'): val = val.split('_')[1]
        except: pass
        return val
    @property
    def itime(self):
        val = self._extKeyword('ITIME', 'Integration time per coadd')
        if val is None: val = self._priKeyword('ITIME')
        if val is None: val = self._priKeyword('TRUITIME')
        return val
    @property
    def prism(self):
        val = self._priKeyword('DISPERSR', 'Dispersing prism')
        try:
            if val.startswith('DISP_'): val = val.split('_')[1]
        except: pass
        return val
 
    disperser=prism
    @property
    def object(self):
        return self._priKeyword('OBJECT', '')
 
    @property
    def lyot(self):
        val = self._priKeyword('LYOTMASK', 'Lyot Mask in IFS')
        try:
            if val.startswith('LYOT_'): val = val.split('_')[1]
        except: pass
        return val
 
    @property
    def occulter(self):
        return self._priKeyword('OCCULTER', 'Occulter FPM Mask in Cal')
    @property
    def apodizer(self):
        val = self._priKeyword('APODIZER', 'Apodizer PPM Mask in AO/OMSS')
        try:
            if val.startswith('APOD_'): val = val.split('_')[1]
        except: pass
        return val
 
    @property
    def filetype(self):
        return self._priKeyword('FILETYPE', 'Type of GPI Data File')
    @property
    def dateobs(self):
        return self._priKeyword('DATE-OBS', 'UT Date of Observation Start')
    @property
    def datetime(self):
        """ return observation start time as astropy.time object"""
        import astropy.time
        return astropy.time.Time( self._priKeyword('DATE-OBS')+" "+self._priKeyword('UTSTART') , format='iso', scale='utc')

    @property
    def version(self):
        return self._priKeyword('DRPVER', 'DRP version used to produce this file')
    @property
    def ncombined(self):
        val = self._priKeyword('DRPNFILE', '# of input files used to create this file.')
        return val if val is not None else 1
    @property
    def readmode(self):
        
        val = self._extKeyword('SAMPMODE')
        readmodes_table = {None: 'Unknown', 1:'single', 2: 'CDS', 3:'MCDS', 4:'UTR'}
        try:
            longval = readmodes_table[val]
        except:
            longval = str(val)
        if longval[0] == 'M' or longval[0]=='U': longval = longval + '-'+str(self._extKeyword('READS'))
        return longval


    def verifyKeywords(self, verbose=False, return_info=False):
        import gpi_pipeline
        import os
        import xlrd

        print "Now verifying "+self.name
        c = gpi_pipeline.GPIConfig()
        fn = os.path.join(c['GPI_DRP_CONFIG_DIR'], 'keywordconfig.txt')
        keywordlist = [t.strip().split('\t') for t in open(fn).readlines()]

        xlf = xlrd.open_workbook( os.path.join(c['GPI_DRP_CONFIG_DIR'], 'GPI_FITS_Header_Keywords.xls'))
        xls = xlf.sheet_by_index(0)
        goodkeywordct=0
        missingkeywordct=0
        optionalkeywordct=0
        value_ok =0
        value_bad=0
        #for parts in keywordlist:
        results = dict()
        results['missingkeywords'] = []
        results['optionalkeywords'] = []
        results['invalidtypes'] = []
        for i in range(xls.nrows):
            keywordname = str(xls.cell_value(i,0)) # cast to str to de-unicode the value
            if keywordname == '': continue
            if keywordname != keywordname.upper(): continue # probably not a keyword
            if keywordname.startswith("#"): continue
            mandatory_or_optional = xls.cell_value(i,2)
            mandatory=(mandatory_or_optional == 'Mandatory')
            extension = xls.cell_value(i,3)
            keywordsrc = xls.cell_value(i,4)
            if verbose: print "Verifying ", (keywordname, extension, mandatory)

            if extension=='Extension' or extension=='Ext.':
                present = keywordname in self._HDUlist[1].header
            elif extension=='PHU': 
                present = keywordname in self._HDUlist[0].header
            elif extension=='Both' or extension=='PHU & Ext.': 
                present = ((keywordname in self._HDUlist[0].header)  and
                         (keywordname in self._HDUlist[1].header))

            if present:
                valid=True # unless we learn otherwise...
                # validate the type
                value = self[keywordname]

                required_type_name=xls.cell_value(i,6)
                req_types = {u'string': str, u'int':int,
                        u'float':float, u'char':str, u'bool':bool,
                        u'double':float, u'real':float, u'integer':int}
                req_type = req_types[required_type_name]
                if not isinstance(value, req_type):
                    print ("ERROR: Incorrect data type for %s: found %s but expected %s\t\t(source=%s)" % (keywordname, type(value).__name__, required_type_name, keywordsrc))
                    results['invalidtypes'].append(keywordname)
                    valid=False

                allowed_values_str = str(xls.cell_value(i,5))
                if '|' in allowed_values_str: # this is an enumerated type
                    allowed_values = [str(i).strip() for i in value.split('|')]
                    if value not in allowed_values:
                        valid=False
                        print ("ERROR: Incorrect data value for %s: found %s but that is not an allowed value.\t\t(source=%s)" % (keywordname, value, keywordsrc))
                        value_bad+=1
                    else:
                        if verbose: print "Verified value OK for %s" % keywordname
                        value_ok+=1

                
                



            if mandatory:
                if present:
                    goodkeywordct+=1
                else:

                    print "ERROR: mandatory keyword %s is missing from the %s header\t\t(source=%s)" % (keywordname,extension, keywordsrc)

                    missingkeywordct+=1
                    results['missingkeywords'].append(keywordname)
            else:
                if present:
                    optionalkeywordct+=1
                    results['optionalkeywords'].append(keywordname)

        print "Verified %d mandatory keywords and %d optional keywords present. %d mandatory keywords missing." % (goodkeywordct, optionalkeywordct, missingkeywordct)
        print "Verified value OK for %d keywords with defined acceptable values. Invalid values found for %d keywords." % (value_ok, value_bad)
        if return_info: return results

    def extractRecipe(self):
        """ Given a pipeline-reduced image which has a Recipe embedded in its
        header history, extract a copy of that recipe for re-use

        Returns the recipe as a string. 
        TODO: also add option for writing to disk as ASCII file?
        """
        if fits.__name__ == 'astropy.io.fits':
            import distutils.version as ver

            try: 
                if not ver.LooseVersion(fits.__version__) > ver.LooseVersion('0.2') : 
                    raise ImportError("extractRecipe requires astropy >= 0.2 (or pyfits >= 3.1)")
            except: pass 

        else:
            if float(fits.__version__ ) < 3.1: raise ImportError("extractRecipe requires pyfits >= version 3.1 (or astropy >= 0.2)")

        comments =  self.priheader['COMMENT'] # being able to retrieve multiple COMMENT values in this way is what sets the
                                              # above requirements for version numbers.
        comments = [c for c in comments if c.startswith('DRF')]
        output = []
        name=''
        for line in comments:
            if line.startswith('DRFN'): name=line[4:]
            elif line.startswith('DRF '): output.append(line[4:])
            elif line.startswith('DRFC<'): output.append(line[4:])
            elif line.startswith('DRFC'): output[-1] = output[-1]+ line[4:] # append continued lines

        output = "\n".join(output)
        return output


class IFSData(GPI_generic_data):
    """ Object class for GPI IFS data

    Provides a simple interface to manipulate files

    """
    def __init__(self, filename, **kwargs):
        GPI_generic_data.__init__(self, filename, **kwargs)
        self.units = 'ADU/coadd'


    def set_units(self, units):
        "Change the units to either ADU per coadd or per sec"
        if self.units == units: return

        if self.units == 'ADU/coadd' and units == 'ADU/s':
            self._HDUlist[self._sci_ext].data /= self.itime
            self.units = 'ADU/s'
        elif self.units == 'ADU/s' and units == 'ADU/coadd':
            self._HDUlist[self._sci_ext].data *= self.itime
            self.units = 'ADU/coadd'



    ### Arithmatic operations ###

    def __add__(self, other):
        new = self.Copy()
        if np.isscalar(other):
            new.name += " + "+str(other)
            new.HDUlist[self._sci_ext].data += other 
        elif isinstance(other, IFSData) : 
            new.name = self.name +" + "+other.name
            new.data = self.data - other.data
        return new

    def __sub__(self, other):
        new = self.Copy()
        if np.isscalar(other):
            new.name += " - "+str(other)
            new.data -= other 
        elif isinstance(other, IFSData) : 
            new.name = self.name +" - "+other.name
            new.data = self.data - other.data
        return new

    def __mul__(self, other):
        new = self.Copy()
        if np.isscalar(other):
            new.name += " * "+str(other)
            new.data *= other 
        elif isinstance(other, IFSData) : 
            new.name = self.name +" * "+other.name
            new.data = self.data * other.data
        return new
    def __div__(self, other):
        new = self.Copy()
        if np.isscalar(other):
            new.name += " / "+str(other)
            new.data /= other 
        elif isinstance(other, IFSData) : 
            new.name = self.name +" / "+other.name
            new.data = self.data / other.data
        return new

    def writeto(self, *args, **kwargs):
        self._HDUlist.writeto(*args, **kwargs)


#--------------------------------------------------------------------------------

class IFSSpectralCube(IFSData):
    """ 3D spectrally dispersed IFS data from GPI

    """

    ### wavelength data
    @property
    def wavelengths(self):
        """ Wavelengths for the slices in this datacube, in microns """
        # in accordance with Gemini standards for WCS, try CDi_j first.
        # only if that's missing would we fall back to CDELTi...
        try:
            cd3 = self['CD3_3']
        except:
            cd3 = self['CDELT3']

        crpix3 = self['CRPIX3']
        if crpix3 == 0:
            _log.info('found crpix=3; assuming this is an older, non-WCS-compliant file')
            _log.info('interpreting as if CRPIX3=1 and continuing...')
            crpix3=1


        return (np.arange(self['NAXIS3']) - (crpix3-1))*cd3 + self['CRVAL3']

    def ExtractSpectrum(self, center=None, full=None, radius=10, method='mean'):
        """ Extract a 1D spectrum from a datacube. 

        Properties
        ------------
        center : 2 element iterable
            Center of location for spectral extraction, in (y,x) Pythonic order
        radius : int
            Radius for spectral extraction
        full : bool
            if True, return mean spectrum for entire cube instead of just a part
        method : string
            either 'mean', 'median' or 'sum', for how to compute it
        """
        from scipy.stats.stats import nanmean, nanmedian
        from numpy import nansum
        if center == None and full== None: full=True

        if method.lower() == 'mean':
            combinefunc = nanmean
        elif method.lower() == 'median':
            combinefunc = nanmedian
        elif method.lower() == 'sum':
            combinefunc = nansum
        else: raise ValueError("Unknown method: "+method+"; must be one of mean, median, or sum")


        if full:
            data_realshape  = self.data.shape # save this aside for a moment...
            # take advantage of array-relabeling to net us easily do this using nanmean
            self.data.shape = (data_realshape[0], data_realshape[1]*data_realshape[2])
            myspec = combinefunc(self.data, axis=1)
            self.data.shape=data_realshape # put back what we saved
        else:
            y,x = np.indices(self.data.shape[1:])
            y-= center[0]
            x-= center[1]
            r = np.sqrt(y**2+x**2)
            wg = np.where(r < radius)

            myspec = np.zeros(self.data.shape[0])
            for i in range(self.data.shape[0]):
                myspec[i] = combinefunc(self.data[i][wg])
     
            #raise NotImplementedError('not yet...')
            
        return myspec


    def PlotSpectrum(self, *args, **kwargs):
        """
        Plot a 1D spectrum extracted from the cube.
        See ExtractSpectrum for optional arguments
        """

        spec = self.ExtractSpectrum(*args, **kwargs)
        pl.plot(self.wavelengths,spec)

    def Fit1Peak(self, *args, **kwargs):
        """ Fit a single Gaussian function to a 1D spectrum extracted from the cube.
        See ExtractSpectrum for optional arguments
        """

        import scipy.optimize

        spec = self.ExtractSpectrum(*args, method='sum', **kwargs)

        spec_uncert = np.sqrt(spec.clip(0)) # FIXME add in gain here...
        spec_uncert[spec_uncert == 0] = np.median(spec_uncert) # this is a quick hack
        waves = self.wavelengths
        pl.clf()
        pl.plot(waves, spec)


        def gaussian(lam, cen, sig, height, offset):
            return np.exp(-((lam-cen)/sig)**2) *height + offset
        initial_guess = [(spec*waves).sum()/spec.sum(), 0.1, spec.max(), 0]

        popt, pcov = scipy.optimize.curve_fit(gaussian, self.wavelengths, spec, p0=initial_guess, sigma=spec_uncert)

        pl.plot(self.wavelengths, gaussian(self.wavelengths, *popt), "r--")
        pl.text(popt[0]+0.01, popt[2]*0.95, "Center: %.3f $\mu m \pm $ %.3f \nFWHM: %.3f $\mu m \pm $ %.3f" % 
                (popt[0], np.sqrt(pcov[0,0]), popt[1], np.sqrt(pcov[1,1])))
        #stop()

        pl.suptitle(self.name)
        return popt

    def FitNPeaks(self, N, *args, **kwargs):
        """ Fit multiple Gaussian functions to a 1D spectrum extracted from the cube.
        All must have same FWHM.
        See ExtractSpectrum for optional arguments
        """

        import scipy.optimize

        spec = self.ExtractSpectrum(*args, method='sum', **kwargs)

        spec_uncert = np.sqrt(spec.clip(0)) # FIXME add in gain here...
        spec_uncert[spec_uncert == 0] = np.median(spec_uncert) # this is a quick hack
        waves = self.wavelengths
        pl.clf()
        pl.plot(waves, spec)


        def gaussian(lam, cen, sig, height, offset):
            """ one gaussian plus constant offset """
            return np.exp(-((lam-cen)/sig)**2) *height + offset

        def ngaussians(lam, cens, heights, sig, offset):
            """ several gaussians, all with same FWHM, plus constant offset """
            result = np.zeros_like(lam)+offset
            for cen,height in zip(cens,heights):
                result += gaussian(lam, cen, sig, height,0)
            return result


        initial_guess = [(spec*waves).sum()/spec.sum(), 0.1, spec.max(), 0]

        popt, pcov = scipy.optimize.curve_fit(gaussian, self.wavelengths, spec, p0=initial_guess, sigma=spec_uncert)

        pl.plot(self.wavelengths, gaussian(self.wavelengths, *popt), "r--")
        pl.text(popt[0]+0.01, popt[2]*0.95, "Center: %.3f $\mu m \pm $ %.3f \nFWHM: %.3f $\mu m \pm $ %.3f" % 
                (popt[0], np.sqrt(pcov[0,0]), popt[1], np.sqrt(pcov[1,1])))
        #stop()

        pl.suptitle(self.name)
        return popt


    def PlotArcLampSpectra(self, filter=None):

        import gpi_pipeline
        import astropy.io.ascii
        gc = gpi_pipeline.GPIConfig()
        waves = self.wavelengths
        if filter is None: filter=self.filter

        for name, color in zip( ('Xe', 'Ar'), ('blue', 'magenta')):
            arcfile_dst = os.path.join(gc.get_directory('GPI_DST_DIR'), name+"ArcLampG.txt")
            table = astropy.io.ascii.read(arcfile_dst)
            table['col1'].name = 'wavelength'
            table['col2'].name = 'intensity'
            table['wavelength'] /= 1e4 # convert angstroms to microns

            #pl.plot(table['wavelength'], table['intensity'], "--", color=color)
            for wave, flux in table:
                pl.plot( [wave, wave], [0,flux],  color=color, linewidth=2)
                #stop()
    
            # rebin array onto observed spectral channels, taking into account sums over
            # each wavelength bin.
            wdelta = waves[1]-waves[0]
            edge_locations = np.concatenate( (waves-wdelta/2, [waves[-1]+wdelta/2]) ) 
            #edge_locations = np.searchsorted( table['wavelength'], wave_bin_edges)
            rebinned_flux = np.zeros_like(waves)
            for i in range(len(waves)):
                wbin = np.where( (table['wavelength'] > edge_locations[i] ) & 
                                 (table['wavelength'] < edge_locations[i+1] ) )
                rebinned_flux[i] = table['intensity'][wbin].sum()
            pl.plot( waves, rebinned_flux, ":", drawstyle='steps-mid', color=color)

        ax = pl.gca()

        # now overplot the filter transmission profile
        ax2 = ax.twinx()
        filterfile= fits.open( os.path.join(gc.get_directory('GPI_DRP_CONFIG_DIR'), "filters","GPI-filter-%s.fits" % filter) )
        ax2.plot( filterfile[1].data['WAVELENGTH'].flat, filterfile[1].data['TRANSMISSION'].flat, "k:")
        ax2.set_ybound(0,1)
        ax2.set_yticklabels([])
        ax.set_xlabel("Wavelength [$\mu m$]")

        if filter=='Y': ax.set_xbound(0.9, 1.2)
        elif filter=='J': ax.set_xbound(1.05, 1.4)
        elif filter=='H': 
            ax.set_xbound(1.4, 1.9)
            ax.set_ybound(0,5e4)
        elif filter=='K1': 
            ax.set_xbound(1.8, 2.3)
            ax.set_ybound(0,2e4)
        elif filter=='K2': 
            ax.set_xbound(2.0, 2.5)
            ax.set_ybound(0,2e4)
 

        # now read in the information the DRP uses for fitting wavelengths
        linelist_file = os.path.join(gc.get_directory('GPI_DRP_CONFIG_DIR'), "lampemissionlines.txt")
        linelist = open(linelist_file).readlines()
        # restrict to the lines we care about for the current filter:
        lines = [l.strip() for l in linelist if l.startswith(filter)]
        yrange = pl.gca().get_ybound()
        for line in lines:
            lineparts = line.split() # will be filter, Xe/Ar, then a set of wavelengths
            filtername, source, wavelens = lineparts[0], lineparts[1], lineparts[2:]
            if source == "Xe": texty = 0.9*yrange[1]
            elif source == "Ar": texty = 0.8*yrange[1]
            else: continue

            for wavelen in wavelens:
                pl.text(float(wavelen), texty, source)

 
        pl.title(filter)


        pl.draw()
        pl.show()





    def ShowSlice(self, index=0, wavelength=None, scale='log', vmin=0, vmax=None, verbose=False):

        if wavelength is not None:
            # find closest wavelength that matches
            raise NotImplementedError('')

        if vmax is None: 
            vmax = np.nanmax(self.data[index])
            if verbose: print "Max value of data: "+str(vmax)

        if scale == 'linear':
            norm=matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        else:
            if vmin == 0: vmin=vmax/1e6
            norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)


        pl.imshow( self.data[index], norm=norm) 
        pl.title(self.name+", slice "+str(index)+": %.3f $\mu m$" % self.wavelengths[index])


#--------------------------------------------------------------------------------

class IFSPolarimetryCube(IFSData):
    """ 3D polarimetry cube IFS data from GPI

    """

    ### wavelength data
    @property
    def stokes(self):
        pass
 
#---------------------------------------------------------------

class IFSData2D(IFSData):
    """ 2D (i.e. dispersed detector format) IFS data for GPI

    """
    def Destripe(self, method='simple'):
        """ Remove horizontal banding using H2RG reference pixels 

        Parameters
        -----------
        method : string
            'simple' or 'linear'

        """
        if self.data.shape != (2048, 2048):
            raise ValueError("Input image doesn't appear to be 2048 x 2048.")

        # fold data into 32




    def DisplayZooms(self, format=(3,3), size=64, vmin=0, vmax=None, scale='log', verbose=False, _right_side=False,
            title=True):
        """ Display zoomed in sub regions of the image 

        Parameters
        -----------
        format : tuple
            How many zoomed in images to display. Default is 3 x 3
        size : int
            How many pixels to show in each zoomed in image

        """
        if self.data.shape != (2048, 2048):
            raise ValueError("Input image doesn't appear to be 2048 x 2048.")

        npix = 2048 - size - 1
        pos_x = np.linspace(size/2, 2047-size/2, format[0])
        pos_y = np.linspace(size/2, 2047-size/2, format[1])

        if vmax is None: 
            vmax = np.nanmax(self.data)
            if verbose: print "Max value of data: "+str(vmax)

        pl.clf()

        if scale == 'linear':
            norm=matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        else:
            if vmin == 0: vmin=vmax/1e6
            norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)


        #fig = pl.figure(num=1, figsize=(7.1,6.9))
        if _right_side: 
            # subplots order is Y, X
            fig, axarr0 = pl.subplots( format[0], format[1]*2,  figsize=(17, 7.65), num=2)
            axarr = axarr0[ :, format[1]:]
        else:
            fig, axarr = pl.subplots( format[0], format[1],  figsize=(7.98,7.65), num=2)
            pl.subplots_adjust(wspace=0.01, hspace=0.03)


        for ix in range(format[0]):
            for iy in range(format[1]):
                subim = self.data[ pos_y[iy]-size/2:pos_y[iy]+size/2, pos_x[ix]-size/2:pos_x[ix]+size/2]
                extent = [pos_x[ix]-size/2, pos_x[ix]+size/2, pos_y[iy]-size/2, pos_y[iy]+size/2]

                axarr[format[1]-1-iy,ix].imshow(subim, norm=norm, extent=extent)

                axarr[format[1]-1-iy,ix].xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=4, integer=True))
                axarr[format[1]-1-iy,ix].yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=4, integer=True))

        for ix in range(1, format[0]):
            for iy in range(format[1]):
                axarr[format[1]-1-iy,ix].yaxis.set_ticklabels("")
        for iy in range(1, format[0]):
            for ix in range(format[1]):
                axarr[format[1]-1-iy,ix].xaxis.set_ticklabels("")

        if title:
            pl.suptitle("%s\n %s, %s" % (self.name, self.prism.capitalize(), self.filter), fontsize=16)
        return axarr


    def DisplayBoth(self, format=(3,3), vmax=None, vmin=0, scale='linear', verbose=False, **kwargs):
        """ Display a full field image at left, zooms at right """
        if self.data.shape != (2048, 2048):
            raise ValueError("Input image doesn't appear to be 2048 x 2048.")

        if vmax is None: 
            vmax = np.nanmax(self.data)
            if verbose: print "Max value of data: "+str(vmax)

        if scale == 'linear':
            norm=matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        else:
            if vmin == 0: vmin=vmax/1e6
            norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)

        axarr = self.DisplayZooms(format=format, vmax=vmax, vmin=vmin, scale=scale, verbose=verbose, _right_side=True, **kwargs)

        # get top and bottom coords of that array
        ymax = axarr[0,0,].get_position().corners()[3,1]
        xmin = axarr[0,0,].get_position().corners()[0,0]
        ymin = axarr[format[0]-1,format[1]-1].get_position().corners()[0,1]
        xmax = axarr[format[0]-1,format[1]-1].get_position().corners()[3,0]
        #print xmin, xmax, ymin, ymax
        #pl.subplot(121)

        #bigax = pl.axes( (0.1, ymin, xmax-xmin, ymax))
        bigax = pl.subplot(121)

        bigax.imshow(self.data, norm=norm)

#---------------------------------------------------------------

class IFSSpectrum1D(IFSData):
    """ 1D spectrum, extracted from a datacube for GPI

    """

    @property
    def wavelength(self):
        return self.data['WAVELENGTH'].ravel()

    @property
    def flux(self):
        return self.data['FLUX'].ravel()

    @property
    def uncert(self):
        return self.data['FLUX_UNCERT'].ravel()

    def plot(self, format='', clear=False, **kwargs):
        """ Plot spectrum """
        if clear: pl.clf()
        pl.plot(self.wavelength.ravel(), self.flux.ravel(), format, label=self.name, **kwargs)

    def DisplayZooms(self, format=(3,3), size=64, vmin=0, vmax=None, scale='log', verbose=False, _right_side=False,
            title=True):
        """ Display zoomed in sub regions of the image 

        Parameters
        -----------
        format : tuple
            How many zoomed in images to display. Default is 3 x 3
        size : int
            How many pixels to show in each zoomed in image

        """
        if self.data.shape != (2048, 2048):
            raise ValueError("Input image doesn't appear to be 2048 x 2048.")

        npix = 2048 - size - 1
        pos_x = np.linspace(size/2, 2047-size/2, format[0])
        pos_y = np.linspace(size/2, 2047-size/2, format[1])

        if vmax is None: 
            vmax = np.nanmax(self.data)
            if verbose: print "Max value of data: "+str(vmax)

        pl.clf()

        if scale == 'linear':
            norm=matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        else:
            if vmin == 0: vmin=vmax/1e6
            norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)


        #fig = pl.figure(num=1, figsize=(7.1,6.9))
        if _right_side: 
            # subplots order is Y, X
            fig, axarr0 = pl.subplots( format[0], format[1]*2,  figsize=(17, 7.65), num=2)
            axarr = axarr0[ :, format[1]:]
        else:
            fig, axarr = pl.subplots( format[0], format[1],  figsize=(7.98,7.65), num=2)
            pl.subplots_adjust(wspace=0.01, hspace=0.03)


        for ix in range(format[0]):
            for iy in range(format[1]):
                subim = self.data[ pos_y[iy]-size/2:pos_y[iy]+size/2, pos_x[ix]-size/2:pos_x[ix]+size/2]
                extent = [pos_x[ix]-size/2, pos_x[ix]+size/2, pos_y[iy]-size/2, pos_y[iy]+size/2]

                axarr[format[1]-1-iy,ix].imshow(subim, norm=norm, extent=extent)

                axarr[format[1]-1-iy,ix].xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=4, integer=True))
                axarr[format[1]-1-iy,ix].yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=4, integer=True))

        for ix in range(1, format[0]):
            for iy in range(format[1]):
                axarr[format[1]-1-iy,ix].yaxis.set_ticklabels("")
        for iy in range(1, format[0]):
            for ix in range(format[1]):
                axarr[format[1]-1-iy,ix].xaxis.set_ticklabels("")

        if title:
            pl.suptitle("%s\n %s, %s" % (self.name, self.prism.capitalize(), self.filter), fontsize=16)
        return axarr


    def DisplayBoth(self, format=(3,3), vmax=None, vmin=0, scale='linear', verbose=False, **kwargs):
        """ Display a full field image at left, zooms at right """
        if self.data.shape != (2048, 2048):
            raise ValueError("Input image doesn't appear to be 2048 x 2048.")

        if vmax is None: 
            vmax = np.nanmax(self.data)
            if verbose: print "Max value of data: "+str(vmax)

        if scale == 'linear':
            norm=matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        else:
            if vmin == 0: vmin=vmax/1e6
            norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)

        axarr = self.DisplayZooms(format=format, vmax=vmax, vmin=vmin, scale=scale, verbose=verbose, _right_side=True, **kwargs)

        # get top and bottom coords of that array
        ymax = axarr[0,0,].get_position().corners()[3,1]
        xmin = axarr[0,0,].get_position().corners()[0,0]
        ymin = axarr[format[0]-1,format[1]-1].get_position().corners()[0,1]
        xmax = axarr[format[0]-1,format[1]-1].get_position().corners()[3,0]
        #print xmin, xmax, ymin, ymax
        #pl.subplot(121)

        #bigax = pl.axes( (0.1, ymin, xmax-xmin, ymax))
        bigax = pl.subplot(121)

        bigax.imshow(self.data, norm=norm)


 
#---------------------------------------------------------------


class IFSwavecal(GPI_generic_data):
    """GPI IFS wavelength calibration file

       Wavelength solution File Format:
       Dispersion for each spectrum is defined as
       lambda=w3*(sqrt((x-x0)^2+(y-y0)^2))+lambda0
          Slice 1:  Y-positions (Y0) of spectra (Y:spectral direction) at wavelength=lambda0
          Slice 2:  X-positions (X0) of spectra at wavelength=lambda0
          Slice 3:  lambda0 [um]
          Slice 4:  Dispersion [um/pixel]
          Slice 5:  tilts of spectra [radians, with 0=+Y direction]
    """

    def change_ref_wavelength(self, new_ref_wavelength):
        """ Change the reference wavelength of a given wavelength solution
            I.e. the wavelength at which the X0 and Y0 positions are specified)
        """

        if new_ref_wavelength < 0.9 or new_ref_wavelength>2.5:
            raise ValueError('Invalid reference wavelength. Must be in range 0.9-2.5 microns, unit=microns')

        new_y0, new_x0 = self.get_positions_at_wavelength(new_ref_wavelength)

        #dist = (new_ref_wavelength-self.wavecal_lambda0)/self.wavecal_dispersion
        #new_x0 = self.wavecal_xpos + dist * np.sin(self.wavecal_tilt)
        #new_y0 = self.wavecal_ypos - dist * np.cos(self.wavecal_tilt)

        self.data[1] = new_x0
        self.data[0] = new_y0
        self.data[2][np.isfinite(self.data[2])] = new_ref_wavelength

        _log.info("Reference wavelength updated to {0:.3f} microns".format(new_ref_wavelength))


    def get_positions_at_wavelength(self, wavelength):
        """ return the positions Y,X for each lenslet at a given requested wavelength """


        if wavelength < 0.9 or wavelength>2.5:
            raise ValueError('Invalid requested wavelength. Must be in range 0.9-2.5 microns, specified in microns')
        dist = (wavelength-self.wavecal_lambda0)/self.wavecal_dispersion
        wave_x0 = self.wavecal_xpos + dist * np.sin(self.wavecal_tilt)
        wave_y0 = self.wavecal_ypos - dist * np.cos(self.wavecal_tilt)

        return (wave_y0, wave_x0)
 

    
    @property
    def wavecal_xpos(self):
        return self.data[1]
    @property
    def wavecal_ypos(self):
        return self.data[0]
    @property
    def wavecal_lambda0(self):
        return self.data[2]
    @property
    def wavecal_dispersion(self):
        return self.data[3]
    @property
    def wavecal_tilt(self):
        return self.data[4]

    def _show1(self, data, pos=None, vmin=None, vmax=None, title="Plot", norm=None, mask=None, **kwargs):
        """ Display one slice of wavelength solution. Utility function used in plot, PlotResolution etc """
        import gpi_utils
        if pos is not None:
            pl.subplot(2,3,pos)
        ax = gpi_utils.imshow_with_mouseover(data,  vmin=vmin, vmax=vmax, norm=norm, **kwargs)
        if mask is not None:
            gooddata = data[mask]
            meanval = gooddata.mean()
            minval = gooddata.min()
            maxval = gooddata.max()
            stdval = gooddata.std()
            #print title
            #print "Data mean=%5g +- %5g, min=%5g, max=%5g" % (meanval, stdval, minval, maxval)
            pl.text(0.5, 0.1, "Mean=%5g +- %5g,\n min=%5g, max=%5g" % (meanval, stdval, minval, maxval), transform=ax.axes.transAxes, 
                horizontalalignment='center', fontsize=10)
        
        pl.title(title)

    def _hist1(self, data, pos=None, range=None, title="Histogram", bins=20, facecolor='none', histtype='step', **kwargs):
        """ Display 1 histogram. Utility function used in PlotHistogram """
        if pos is not None:
            pl.subplot(2,3,pos)
        wg = np.isfinite(data)
        pl.hist(data[wg], range=range, bins=bins, facecolor=facecolor, histtype=histtype, **kwargs)
        pl.title(title)

 
    def plot(self):
        """ Plot the 5 components of a wavelength solution """

        mask = np.isfinite(self.data.sum(axis=0)) 

        pl.clf()
        self._show1(self.wavecal_ypos, 1, title='Y positions of spectra \n[pix]', vmin=0, vmax=2048, mask=mask)
        self._show1(self.wavecal_xpos, 2, title='X positions of spectra \n[pix]', vmin=0, vmax=2048, mask=mask)
        self._show1(self.wavecal_lambda0, 4, title='Initial wavelength \n[micron]', vmin=0.8, vmax=2.5, mask=mask)
        self._show1(self.wavecal_dispersion, 5, title='Dispersion \n[micron/pix]', vmin=0.01, vmax=0.02, mask=mask)
        self._show1(np.rad2deg(self.wavecal_tilt), 6, title='Tilt \n[deg]', vmin=-20, vmax=20, mask=mask)

        pl.gcf().suptitle(self.name, fontsize=16)

    def plotResolution(self, vmin=30, vmax=90):
        """ Plot spectral resolution over the FOV"""
        mask = np.isfinite(self.data.sum(axis=0)) 

        pl.clf()
        self._show1(self.wavecal_lambda0/self.wavecal_dispersion/2, title='Spectral Resolution per 2 pixels', vmin=vmin, vmax=vmax, mask=mask)

    def plotSeparations(self, vmin=2, vmax=5):
        if self.filter=='Y': lrange = (0.95, 1.14)
        elif self.filter=='J': lrange = (1.118, 1.345)
        elif self.filter=='H': lrange = (1.495, 1.794)
        elif self.filter=='K1': lrange = (1.899, 2.186)
        elif self.filter=='K2': lrange = (2.118, 2.384)

        #compute positon of longest end of spectrum:
        # The closest adjacent spectrum to this is going to be either one higher in Y lenslet index, 
        # or one higher in Y and one lower in X. 
        #dist0 = (lrange[0]-self.wavecal_lambda0)/ self.wavecal_dispersion
        #xend0 = self.wavecal_xpos + dist * np.sin(self.wavecal_tilt)
        #yend0 = self.wavecal_ypos + dist * np.cos(self.wavecal_tilt)

        # compute the X difference for the lenslet 1 higher in Y value
        adjacent1_ypos = np.roll(self.wavecal_ypos, -1, axis=0)  
        adjacent1_xpos = np.roll(self.wavecal_xpos, -1, axis=0)  
        # compute offset in Y to adjacent spectrum
        deltay_1 =self.wavecal_ypos - adjacent1_ypos
        # compute x position of current spectrum at that Y
        deltax_1 = deltay_1 * np.tan(self.wavecal_tilt) 
        # compute difference from next spectrum position
        xoffset_1 = np.abs(self.wavecal_xpos + deltax_1 - adjacent1_xpos)

        # repeat calc for adjacent pixel, one lower in X
        adjacent2_ypos = np.roll(np.roll(self.wavecal_ypos, -1, axis=0), 1, axis=1)  
        adjacent2_xpos = np.roll(np.roll(self.wavecal_xpos, -1, axis=0), 1, axis=1)  
        # compute offset in Y to adjacent spectrum
        deltay_2 =self.wavecal_ypos - adjacent2_ypos
        # compute wavelength at which that is the case
        deltax_2 = deltay_2 * np.tan(self.wavecal_tilt) 
        xoffset_2 = np.abs(self.wavecal_xpos + deltax_2 - adjacent2_xpos)

        xoffsets = np.zeros((2, xoffset_2.shape[0], xoffset_2.shape[1]))
        xoffsets[0] = xoffset_1
        xoffsets[1] = xoffset_2
        minxoffset = np.min(xoffsets, axis=0)
        pl.imshow(minxoffset, vmin=vmin, vmax=vmax)
        pl.title("Minimum separation between adjacent spectra for: "+self.name)

        cb = pl.colorbar()
        cb.set_label("Separation in pixels")
        cb.set_ticks([2,2.5, 3,3.5, 4,4.5, 5])
        pl.text(20, 20, "Median: %.2f" % np.median(minxoffset[np.isfinite(minxoffset)]))
        pl.text(200, 40, "Minimum: %.2f" % np.nanmin(minxoffset))
        pl.text(200, 20, "Maximum: %.2f" % np.nanmax(minxoffset))

    
    def plotComparison(self, otherfile):
        """ Compare 2 different wavecal files, with plots
        """
        import gpi_pipeline
        gc = gpi_pipeline.GPIConfig()
        if isinstance(otherfile, IFSwavecal):
            other = otherfile
        else:
            other = IFSwavecal(otherfile)
        difference = self.data - other.data

        #mask = (self.data[0] != 0)
        mask = np.isfinite(difference.sum(axis=0)) 

        pl.clf()
        #--- Create 5 plots with 2D images showing differences of the files
        self._show1(difference[0], 1, title='Y positions of spectra \n[pix]', mask=mask, vmin=-3, vmax=3)
        self._show1(difference[1], 2, title='X positions of spectra \n[pix]', mask=mask, vmin=-1, vmax=1)
        self._show1(difference[2], 4, title='Initial wavelength \n[micron]',  mask=mask, vmin=-0.2, vmax=0.2)
        self._show1(difference[3], 5, title='Dispersion \n[micron/pix]', mask=mask, vmin=-0.001, vmax=0.001)
        self._show1(np.rad2deg(difference[4]), 6, title='Tilt \n[deg]', mask=mask, vmin=-3, vmax=3)

        pl.gcf().suptitle("Diff between %s (%s) and %s (%s)" % ( self.name, self.priheader['DATE-OBS'], other.name, other.priheader['DATE-OBS']), fontsize=16)

        #--- create one more panel showing a plotted representation of the wavelength solution
        # now read in the information the DRP uses for fitting wavelengths
        linelist_file = os.path.join(gc.get_directory('GPI_DRP_CONFIG_DIR'), "lampemissionlines.txt")
        linelist = open(linelist_file).readlines()
        # restrict to the lines we care about for the current filter:
        lines = [l.strip() for l in linelist if l.startswith(self.filter)]


        ax3 = pl.subplot(2,3,3)
        for wavcal in [self, other]:
            # find center pixel
            radius = np.sqrt((wavcal[0]-1024)**2+(wavcal[1]-1024)**2)
            radius[ ~np.isfinite(wavcal[0])] = 1e6 # mask out regions outside the array. 
            wm = np.where(radius == radius.min())
            #wm = np.argmin(radius)
            x0 = wavcal[1][wm]
            y0 = wavcal[0][wm]
            lambda0 = wavcal[2][wm]
            disp0 = wavcal[3][wm]
            tilt0 = wavcal[4][wm] # radians
            

            dist = 0.3 / disp0 # FIXME should be not hard coded.
            xend = x0 + dist * np.sin(tilt0)
            yend = y0 + dist * np.cos(tilt0)

                #DISP[0,0]=X2[0,ii]  
                #DISP[1,0]=X2[1,ii]
                #d2=(lambc-wavcal(ii,jj,2))/wavcal(ii,jj,3)
                #DISP[0,1]=d2*sin(wavcal(ii,jj,4))+wavcal(ii,jj,1)
                #DISP[1,1]=-d2*cos(wavcal(ii,jj,4))+wavcal(ii,jj,0)

            _log.info("Central spectrum goes from :"+ str([x0,y0]) )
            _log.info("                      to   :"+str([xend,yend]) )
            plotline = pl.plot([x0,xend],[y0,yend], label=wavcal['DATE-OBS'])[0]
            # get current color of plot line here
            col = plotline.get_color()


            for il, line in enumerate(lines):
                lineparts = line.split() # will be filter, Xe/Ar, then a set of wavelengths
                filtername, source, wavelens = lineparts[0], lineparts[1], lineparts[2:]
                if source == "Xe": marker = 'x'
                elif source == "Ar": marker='.'
                else: continue

                for wavelen in wavelens:
                    dist = (float(wavelen)-lambda0)/disp0
                    xend = x0 + dist * np.sin(tilt0)
                    yend = y0 + dist * np.cos(tilt0)
                    plotline = pl.plot(xend,yend,marker=marker, color=col )
                    if wavcal is self:
                        pl.text(xend+2, yend, source +" "+wavelen)

                
        ax3.set_xbound(1020, 1040)
        ax3.set_ybound(1020, 1044)
        for axis in [ax3.xaxis, ax3.yaxis]:
            axis.set_major_locator(matplotlib.ticker.MaxNLocator(5))
            axis.set_minor_locator(matplotlib.ticker.MultipleLocator(1))
        ax3.grid(True, which='minor', color='black', linestyle=':', alpha=0.4)
        ax3.grid(False, which='major')
        ax3.set_aspect('equal')
    
        pl.legend(loc='lower left', frameon=False)

        stop()

        # display a 2D region with equal aspect ratio...
 
    def plotHistograms(self, erase=True, suptitle=True, **kwargs):
        """ Display histograms of the quantities of interest in a wavelength solution """
        if erase: pl.clf()
        pl.subplots_adjust(hspace=0.4)
        self._hist1(self.wavecal_ypos, 1, title='Y positions of spectra \n[pix]', range=[0,2048], **kwargs)
        self._hist1(self.wavecal_xpos, 2, title='X positions of spectra \n[pix]', range=[0,2048], **kwargs)
        self._hist1(self.wavecal_lambda0, 4, title='Initial wavelength \n[micron]', range=[0.8, 2.5], **kwargs)
        self._hist1(self.wavecal_dispersion, 5, title='Dispersion \n[micron/pix]', range=[0.01, 0.02], **kwargs)
        self._hist1(np.rad2deg(self.wavecal_tilt), 6, title='Tilt \n[deg]', range=[-20,20], **kwargs)

        if suptitle: pl.gcf().suptitle(self.name, fontsize=16)



#-------------------------------------------------------------------------------

class DataCollection(object):
    def __init__(self, filelist_or_directory_or_glob, loadData=False, simplifyKeywords=False):
        """ A collection of several GPI data files

        This acts in many ways like an OrderedDict, the contents of which are
        some number of gpidata objects. 

        >>> dc = DataCollection('/some/path/files*.fits')
        >>> len(dc)
        36
        >>> dc[4]
        <gpidata.IFSData2D at 0x104e12fd0>

        That said, it also acts a bit like a table, in that you can access various
        keywords for the whole set of files at once
        >>> dc.ifsfilters
        ['IFSFILT_H_G1213','IFSFILT_H_G1213','IFSFILT_H_G1213',...]



        Parameters
        ------------
        loadData : bool
            Defaults to False
        simplifyKeywords : bool
            if True, keyword values will be returned with the Gemini-required
            formatting cruft trimmed off. Default is False, return the actual
            value strings.

        """
        self.simplifyKeywords=simplifyKeywords
        self.loadData=loadData
        import glob
        import collections
        import os
        if isinstance(filelist_or_directory_or_glob, str) and os.path.isdir(filelist_or_directory_or_glob):
            #we have a directory
            filenames = glob.glob(os.path.join(filelist_or_directory_or_glob, "*.fits"))
        elif isinstance(filelist_or_directory_or_glob, str) and (('*' in filelist_or_directory_or_glob) or ('?' in filelist_or_directory_or_glob)):
            # we have a wildcard spec
            filenames = glob.glob(filelist_or_directory_or_glob)
        elif hasattr(filelist_or_directory_or_glob, '__iter__'):
            filenames =filelist

        _log.info("Found %d files" % len(filenames))
        # FIXME validate existence of files?
        self._contents = collections.OrderedDict( [(os.path.abspath(filename),read(filename,loadData=loadData)) for filename in filenames])

    def __repr__(self):
        return "<%s.%s : containing %d files >" % (self.__class__.__module__, self.__class__.__name__, len(self._contents)) 

    @property
    def filenames(self):
        return self._contents.keys()
 
    def __getitem__(self,key):
        """ If passed an integer N, return the Nth file
        If passed a string, try to return a file with that name
        """
        import gpi_utils
        if gpi_utils.is_integer(key):
            return self._contents[self._contents.keys()[key]]
        else:
            return self._contents[key]

    def __len__(self):
        return len(self._contents)

    def AddFile(self, filename):
        """ Add one file to the DataCollection """
        self._contents[os.path.abspath(filename)] = read(filename,loadData=self.loadData)
    def RemoveFile(self, filename):
        """ Remove one file from the DataCollection """
        self._contents.pop(os.path.abspath(filename))


    def keys(self):
        return self._contents.keys()
    def values(self):
        return self._contents.values()

    # allow iteration over the ordered dict 
    def __iter__(self):
        return self._contents.__iter__()
    def next(self):
        return self._contents.next()



    def getKeywordValues(self, keyword, simplify=None):
        """ Query all the files for some keyword and 
        return a list of the keyword values.

        Any files which lack that keyword will have a None in their place
        in the returned array.

        Optionally, can trim off Gemini-standard keyword cruft.
        """
        results = []
        for item in self._contents.values():
            try:
                value = item[keyword]
                if simplify is not None:
                    if str(value).startswith(simplify): value=value.split('_')[1]

                results.append(value)
            except:
                results.append(None)
        return np.asarray(results) # cast to ndarray so we can use where()
        #return [item[keyword] for item in self._contents.values()]

    @property 
    def ifsfilters(self):
        return self.getKeywordValues('IFSFILT', simplify='IFSFILT_' if self.simplifyKeywords else None)
    @property 
    def itimes(self):
        return self.getKeywordValues('ITIME')
    @property 
    def apodizers(self):
        return self.getKeywordValues('APODIZER', simplify='APOD_' if self.simplifyKeywords else None)
    @property 
    def lyotmasks(self):
        return self.getKeywordValues('LYOTMASK', simplify='LYOT_' if self.simplifyKeywords else None)
    @property 
    def dispersers(self):
        return self.getKeywordValues('DISPERSR', simplify='DISP_' if self.simplifyKeywords else None)
    @property 
    def obstypes(self):
        return self.getKeywordValues('OBSTYPE')
    @property 
    def obsclasses(self):
        return self.getKeywordValues('OBSCLASS')
    @property 
    def objects(self):
        return self.getKeywordValues('OBJECT')
    @property 
    def occulters(self):
        return self.getKeywordValues('OCCULTER', simplify='FPM_' if self.simplifyKeywords else None)
    @property 
    def dates(self):
        return self.getKeywordValues('DATE-OBS')

    def getTimes(self):
        """ Returns the times for all files as astropy.time.Time objects """
        import astropy.time
        utds = [d+' '+u for u, d in zip(self.getKeywordValues('UTSTART'), self.dates)]
        return astropy.time.Time(utds, format='iso', scale='utc')


    def getObjects(self):
        """ Return a list of gpidata IFSdata objects (of appropriate sub-types) for all the files in this data collection """
        return self._contents.values()

    def Describe(self):
        """ Print to the screen a nice tabular summary of the files in the data collection
        """
        oldSimplifyKeywords=self.simplifyKeywords
        self.simplifyKeywords=True # the following code assumes simplify is True
 
        objects= self.objects
        filters= self.ifsfilters
        obstypes= self.obstypes
        obsclasses=self.obsclasses
        occulters=self.occulters
        dispersers=self.dispersers
        itimes=self.itimes
        readmodes0 = self.getKeywordValues('SAMPMODE')
        readmodes_table = {None: 'Unknown', 1:'single', 2: 'CDS', 3:'MCDS', 4:'UTR'}
        readmodes1 = [readmodes_table[r] for r in readmodes0]
        nreads = self.getKeywordValues('READS')
        readmodes = [(readmodes_table[r]+"-"+str(n) if r >2 else "CDS") for r, n in zip(readmodes0, nreads)]
        print "%-30s\t%2s\t%s\t%5s\t%s" % ("Filename", "Filt","Obstype","ITime", "Readmode")
        for i in range(len(self)):
            print "%-30s\t%2s\t%s\t%.1f\t%s" % (self.filenames[i], str(filters[i]), str(obstypes[i]), itimes[i], readmodes[i])
        self.simplifyKeywords=oldSimplifyKeywords
 

    def ValidateFiles(self):
        """ Validate that all files in the collection are valid GPI data files """
        raise NotImplementedError("not yet")

    def ParseToRecipes(self, autosave=False, verbose=False):
        """ Parse GPI data file to generate recipes
        
        Given some set of files, parse the headers and generate some 
        set of recipes that would optimally reduce them. The resulting
        recipes are not written to disk, but are instead returned as a
        list. The order of elements in this list matters.

        This is a translation to Python of the core parsing algorithm from
        IDL parsergui__define.pro

        """
        
        import re
        import gpi_pipeline
        logging.basicConfig(level=logging.INFO)
        _log.info("Now parsing %d data files " % len(self))

        #--- Zeroth stage of parsing: Validate all files have the required keywords...


        #--- First stage of parsing: clean up and tweak the header keyword values into the forms we care about
        oldSimplifyKeywords=self.simplifyKeywords
        self.simplifyKeywords=True # the following code assumes simplify is True
        # create local copies of some of the keyword arrays, for cases where we need to adjust the values
        objects= self.objects
        filters= self.ifsfilters
        obstypes= self.obstypes
        obsclasses=self.obsclasses
        occulters=self.occulters
        dispersers=self.dispersers
        itimes=self.itimes
        for i in range(len(self)):

            if filters[i] is None: continue # this file is no good.
            if obstypes[i] is None: obstypes[i] = 'SCIENCE'
            if obsclasses[i] is None: obsclasses[i] = 'Object'
            if objects[i] is None: objects[i] = '-Unknown-'


            # we want Xenon&Argon considered as the same 'lamp' object for Y,K1,K2bands (for H&J, better to do separately to keep only meas. from Xenon)
            if re.match('[YK]+[12]*',filters[i]) and (objects[i].lower() == 'xenon' or objects[i].lower() == 'argon'): objects[i] = 'Lamp'
            # we also want Flat considered for wavelength solution in Y band
            if filters[i] =='Y' and objects[i].lower().startswith('flat'): objects[i]='Lamp'
            #Mark filter as irrelevant for Dark exposures
            if obstypes[i].lower() =='dark': 
                filters[i] ='-'

            #Lots of mucking around with obsclass required
            ot = obstypes[i].lower()
            if 'object' in ot: obsclasses[i] = 'science'
            #if 'standard' in ot: 
                #if 'WOLL' in dispersers[i]: obsclasses[i] = 'POLARSTD'
                #else: obsclasses[i] = 'SPECSTD'
            #TODO check ASTROMTC keyword here??

        #--- Second stage of parsing: Generate uniq ordered lists to iterate over
        # in some cases these are just hard coded to all the avilable values.
        #uniqfilters=['Y','J','H','K1','K2']
        uniqfilters=set(filters)
        uniqobstypes=['DARK', 'ARC', 'FLAT', 'STANDARD', 'OBJECT'] # 'BIAS' also allowed but irrelevant/not used for GPI
        uniqobsclasses=set(obsclasses)
        uniqdispersers = ['PRISM','WOLLASTON','OPEN']
        uniqocculters=set(occulters)
        uniqitimes=list(set(itimes))
        uniqitimes.sort() #always sort itimes
        uniqobjects=set(objects)

        _log.info("  Found data for filters: "+ ",".join(set(filters)))
        _log.info("  Found data for itimes: "+ ",".join([str(t) for t in uniqitimes]))
        _log.info("  Found data for objects: "+ ",".join(uniqobjects))
        _log.info("  Found data for obstypes: "+ ",".join(set(obstypes)))


        # TODO: iterate obstypes in a particular order
        # Dark, Arc, Flat, Object, in that order

        #TODO there was also some notion in the IDL code of using Y band flats in wavelength solution
        # but that appears to have been commented out and now deprecated?

        #--- Third stage: Iterate over the unique combinations and generate the output Recipes from the templates
        OutputRecipes = []


        for ifilter, current_filter in enumerate(uniqfilters):
          for iobstype, current_obstype in enumerate(uniqobstypes):
            for idisperser, current_disperser in enumerate(uniqdispersers):
              for iocculter, current_occulter in enumerate(uniqocculters):
                for iobsclass, current_obsclass in enumerate(uniqobsclasses):
                  for iitime, current_itime in enumerate(uniqitimes):
                    for iobject, current_object in enumerate(uniqobjects):
                      wmatch = np.where( (filters == current_filter) & 
                                         (obstypes == current_obstype) &
                                         (dispersers == current_disperser) &
                                         (occulters == current_occulter) &
                                         (obsclasses == current_obsclass) &
                                         (itimes == current_itime) &
                                         (objects == current_object) )
                      matchct = len(wmatch[0])
                      if verbose:
                          print "For parameters: "
                          print ("    "+ "-".join([current_filter,current_obstype,current_disperser,
                              current_occulter,current_obsclass,str(current_itime),current_object]))
                          print "%d files match. " % matchct
                      if matchct == 0: continue
                      current_files = np.asarray(self.filenames)[wmatch]
                      #if matchct > 0: stop()

                      current_obstype=current_obstype.upper()
                      if current_obstype == 'DARK':
                          templatenames=['Dark']
                      elif current_obstype == 'ARC':
                          if current_disperser =='WOLLASTON':
                              templatenames=['Create Polarized Flat-field']
                          else:
                              templatenames=['Wavelength Solution']
                      elif current_obstype == 'FLAT':
                          if current_disperser =='WOLLASTON':
                              templatenames=["Calibrate Polarization Spots Locations",'Create Polarized Flat-field']
                          else:
                              tempplatenames=['Flat-field Extraction']
                      elif current_obstype == 'STANDARD':
                          if current_disperser =='WOLLASTON':
                            raise NotImplementedError('not yet')
                          else:
                            raise NotImplementedError('not yet')

                      elif current_obstype == 'OBJECT':
                          if current_disperser =='WOLLASTON':
                              tempplatenames=['Basic Polarization Sequence']
                          elif current_disperser =='PRISM':
                              #TODO may wish to add various extra calibration steps here?
                              # science data: 
                              if matchct >= 5:
                                  templatenames=['Calibrated Data-cube extraction, ADI reduction']
                              else:
                                  templatenames=['Calibrated Datacube Extraction']
                          else:
                              templatenames=["Simple Undispersed Image Extraction"]
                      else:
                          _log.warn("Not sure how to process these files... No recipe selection criteria match.")
                          _log.warn("obstype=%s, obsclass=%s, filter=%s" % (current_obstype, current_obsclass, current_filter))
                          templatenames=[]
  
                      for template in templatenames:
                          OutputRecipes.append( gpi_pipeline.RecipeFromTemplate(template, current_files) )
  

        print "After parsing, generated %d recipes to process " % len(OutputRecipes)
        #--- clean up and return results
        self.simplifyKeywords=oldSimplifyKeywords

        if autosave:
            for i, rec in enumerate(OutputRecipes):
                outputfn = 'parsed_recipe_%d.xml' % i
                rec.write(outputfn)
                _log.info("Output recipe saved to "+outputfn)

        return OutputRecipes




#-------------------------------------------------------------------------------


def logsheet(dir=".", output="./ifs_files_log.txt", pattern="ifs*.fits"):
    import re
    files_all = glob.glob(pattern)
    # ignore individual reads
    files = [f for f in files_all if re.search("_[0-9][0-9][0-9].fits", f) is None]


    for f in files:
        info = IFSData(f)
        try: 
            #stop()
            print "%-20s\t%s\t%s\t%.2f\t%s\t%s" % (info.name, info.filter, info.lyot , info.itime, info.apodizer, info._priKeyword('TARGET'))
        except:
            print "%-20s\t%s\t%s\t%s\t%s\t%s" % (info.name, info.filter, info.lyot , str(info.itime), info.apodizer, info._priKeyword('TARGET'))



__all__ = [IFSData, logsheet]

if __name__ == "__main__":
        
    logging.basicConfig(level=logging.INFO, format='%(name)-12s: %(levelname)-8s %(message)s',)

    pass

