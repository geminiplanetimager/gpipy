#!/usr/bin/env python
import os
import collections
import numpy as np
import matplotlib.pyplot as pl
import matplotlib
try:
    import astropy.io.fits as fits
except:
    import pyfits as fits
from IPython.core.debugger import Tracer; stop = Tracer()
import logging
from . import utils as gpi_utils

_log = logging.getLogger('gpi')

try:
    from lxml import etree
except ImportError:
    _log.warn("Unable to import lxml.etree, falling back to regular xml.etree instead")
    import xml.etree.cElementTree as etree




""" Python GPI Data Pipeline code

 - Python wrapper to invoke the IDL DRP
 - Pythonic equivalents to various bits of the IDL GPI DRP 

Key Classes:

 - GPIConfig : Interface to GPI Pipeline Configuration System
 - GPICalDB : Interface to GPI Pipeline Calibration Database
 - Recipe : Object interface to represent & manipulate Recipe files


 
HISTORY:

  2012-10-10 MP: started
"""

#---------------------------------------------------------------------------------
class GPIConfig(object):
    """ GPI DRP Configuration

    Settings are accessible by  using this object like a hash.
    This works for both config settings and directory names.

    >>> result = gpiconfig['memsripple']
    >>> result = gpiconfig['GPI_RAW_DATA_DIR']

    These settings are read only; one must change the pipeline configuration files
    directly or adjust environment variables to modify settings.

    """
    def __init__(self):

        self._read_config_files()
        self.paths = {}


        for path in ['RAW_DATA','REDUCED_DATA', 'DRP_TEMPLATES', 'DRP_CONFIG', 'DRP_LOG', 'RECIPE_OUTPUT' , 'CALIBRATIONS', 'DRP_QUEUE']:
            key = 'GPI_%s_DIR' % path
            self.paths[key] = self.get_directory(key)


    def _read_config_files(self):
        user_settings_file = os.path.expanduser("~")+os.sep+".gpi_pipeline_settings"
        global_settings_file = self.get_directory("GPI_DRP_CONFIG_DIR")+os.sep+"pipeline_settings.txt"

        self.settings = {}

        def _parse_one_file(filename):
            _log.debug("Reading in config from "+filename)
            for line in open(filename):
                if line.startswith('#'): continue
                try:
                    parts = line.rstrip().split()
                    # Cast to integer or float if possible.
                    if parts[1].isdigit(): parts[1] = int(parts[1])
                    else:
                        try:
                            parts[1] = float(parts[1])
                        except:
                            pass
                    self.settings[parts[0].lower()] = parts[1]
                    
                except:
                    _log.debug("Could not parse line: '%s' in file %s" % (line, filename))
                    pass


        if os.path.exists(global_settings_file):
            _parse_one_file(global_settings_file)
        if os.path.exists(user_settings_file):
            _parse_one_file(user_settings_file)
        # note that we have to parse the user file *last* so its
        # settings will override the global one, for any values
        # listed in both. 

           


    def __getitem__(self, item):
        """ Look in either config files or paths from
        environment variables etc, and return the
        requested value """
        try:
            return self.settings[item.lower()]
        except:
            try:
                return self.paths[item.upper()]
            except:
                return None
                # NOTE: don't raise an error here - th
                #raise KeyError('Does not appear to be a valid configuration string or directory name: '+item)


    def get_directory(self, dirname, return_method=False, expand=True):
        """ Get a GPI directory absolute path from a directory name string.

        Python equivalent to IDL `gpi_get_directory()`

        Parameters
        -------------
        return_method : bool
            If set, return a tuple containing both the expanded path and a
            string description of where that path was set.
            Default False
        expand : bool
            If set, return an absolute version of the fully expanded pathname
            with all shortcuts such as ~ for the home directory removed. 
            Default True.
        """
        dirname = dirname.upper()
        if not dirname.endswith('_DIR'): dirname=dirname+"_DIR"
        if not dirname.startswith('GPI_'): dirname="GPI_"+dirname
        result=None

        # Case 1: try environment variables
        res = os.getenv(dirname)
        if res is not None and res != "":
            result=res
            method='environment variable'

        # Case 2: Try config files, for the ones which don't cause infinite recursion
        if result is None:
            if (dirname != 'GPI_DRP_CONFIG_DIR') and (dirname != 'GPI_DRP_DIR'):
                try:
                #if self[dirname] is not None:
                    result = self[dirname]
                    method='configuration files'
                except:
                    pass
                    #result = self[dirname]

        # case 3: Try defaults
        if result is None:
            method = 'default value'
            if dirname == 'GPI_DRP_DIR': 
                # find where this current file is
                #raise NotImplementedError('not yet')
                result = '/Users/mperrin/Dropbox/GPI/sw/pipeline'
            elif dirname == 'GPI_DRP_TEMPLATES_DIR':   
                result = self.get_directory("GPI_DRP_DIR")+os.sep+"recipe_templates"
            elif dirname == 'GPI_DRP_CONFIG_DIR':      
                result = self.get_directory("GPI_DRP_DIR")+os.sep+"config"
            elif dirname == 'GPI_DRP_LOG_DIR':         
                result = self.get_directory("GPI_REDUCED_DATA_DIR")+os.sep+"logs"
            elif dirname == 'GPI_RECIPE_OUTPUT_DIR':   
                result = self.get_directory("GPI_REDUCED_DATA_DIR")+os.sep+"recipes"
            elif dirname == 'GPI_DRF_OUTPUT_DIR':      
                result = self.get_directory("GPI_RECIPE_OUTPUT_DIR") # back compatible alias
            elif dirname == 'GPI_CALIBRATIONS_DIR':    
                result = self.get_directory("GPI_REDUCED_DATA_DIR")+os.sep+"calibrations"
            elif dirname == 'GPI_DST_DIR': 
                result = os.path.dirname(gpi_utils.idl_path_find_file('dst'))
                method='location of dst.pro in $IDL_PATH'
            elif dirname == 'GPI_DST_DIR': 
                # grandparent directory of gpipiperun.pro should be pipeline root
                result = os.path.dirname(os.path.dirname(gpi_utils.idl_path_find_file('gpipiperun')))


            else:
                raise ValueError('Could not find default value for '+dirname+'; that is an unknown directory name.')

        if expand:
            result = os.path.abspath(os.path.expandvars(os.path.expanduser(result)))

        if return_method:
            return (result, method)
        else:
            return result

    def expand_path(self, pathname):
        """ Expand paths using GPI named directories
        This includes handling the case of direcory definitions not as 
        environment variables but part of the GPI data pipeline configurations

        """
        if '$GPI' in pathname or '${GPI' in pathname:
            parts = pathname.split(os.sep)


        else:  #no expansion needed of special GPI vars
            return os.path.abspath(os.path.expanduser(os.path.expandvars(string)))

#---------------------------------------------------------------------------------
class GPIPrimitivesConfig(object):
    """ A class for the primitives configuration file
    
    Not a lot of functionality here yet. 

    >>> pc = GPIPrimitivesConfig()
    >>> pc.ListPrimitives()
    [... list of available primitives is displayed...
    >>> pc.primitive['Accumulate Images']
    [... returns a dict of information about that primitive...]


    """
    def __init__(self):
        conf = GPIConfig()
        
        self.configdir = conf.get_directory('GPI_DRP_CONFIG')
        self.filename = os.path.join(self.configdir, 'gpi_pipeline_primitives.xml')

        self._tree = etree.parse(self.filename)
        # Define an internal helper class
        class _primitive_accessor(object):
            """ Wrapper interface class to allow accessing the named
            <module> (AKA "primitive") elements in an XML Recipe file using a Dict-like interface. 

            Don't use this on it's own, use it as the primitive element of a DRF object: 

            mydrf = DRF('some_file.xml')
            mydrf.module['Dark subtraction']

            """

            def __init__(self, parent_XML):
                self._XML= parent_XML

            def __getitem__(self, i):

                if gpi_utils.is_integer(i):
                    # return the i'th module
                    primitives = [mod for mod in self._XML._tree.getroot().iter('Primitive')]
                    if i <0 or i> len(primitives):
                        raise KeyError("Module index must be between 0 and %d for this DRF." %(len(primitives)-1) )
                    else:
                        return self._xml2dict(primitives[i])

                else:
                    # try using it as a string index like a hash

                    for mod in self._XML._tree.getroot().iter('Primitive'):
                        if mod.attrib['Name'].lower() == str(i).lower():
                            return self._xml2dict(mod)
                    raise KeyError("No primitive found with name %s" % str(i))
            def _xml2dict(self, primitive_xml):
                mydict = dict()
                for key in ['Name', 'IDLFunc', 'Comment', 'Order', 'Type','Sequence']:
                    try:
                        mydict[key] = primitive_xml.attrib[key]
                        if key == 'Order': mydict[key] = float(mydict[key]);
                    except:
                        mydict[key] = None

                arguments = []
                for arg in primitive_xml.iter('Argument'):
                    argdict = dict()
                    for key in ['Name', 'Type','Range','Default', 'Desc']:
                        longkey = 'Description' if key=='Desc' else key
                        try:
                            argdict[longkey] = arg.attrib[key]
                        except:
                            argdict[longkey] = None
                    arguments.append(argdict)
                mydict['Arguments'] = arguments
                
                return mydict

    
        self.primitive= _primitive_accessor(self)


    def ListPrimitives(self):
        """ List all available primitives """
        return [mod.attrib['Name'] for mod in self._tree.getroot().iter('Primitive')]

    def primitives(self, testing=True, deprecated=True):

        result = [Primitive(self.primitive[name]) for name in self.ListPrimitives()] # for mod in self._tree.getroot().iter('Primitive')]
        if testing == False:
            result = [p for p in result if 'TESTING' not in p.type.upper()]
        if deprecated == False:
            result = [p for p in result if 'DEPRECATED' not in p.type.upper()]
 
        return result


    def ListOrderedPrimitives(self, return_list=False, type=None, testing=False,deprecated=False):
        """ Return a list of all primitives, in order of recommended execution 

        Parameters:
        -------------
        return_list : bool
            return a list? or else just print to scrren
        type : string
            primitive type (e.g. 'SpectralScience') to restrict the list
            leave undefined for all
        testing : bool
            Include primitives of category='Testing'?
        deprecated : bool
            Include primitives which have been deprecated?
        """
        all_primitives = self.primitives(testing=testing, deprecated=deprecated)

        all_primitives=[a for a in all_primitives if not(a.idlfunc.startswith('test')) ]

        if type is not None:
            all_primitives=[a for a in all_primitives if ((type in a.type) or ('ALL' in a.type ))]

        names = np.asarray( [a.name for a in all_primitives]) 
        orders = np.asarray( [a.order for a in all_primitives] )
        sorter = np.argsort(orders)
        if return_list:
            return zip(orders[sorter], names[sorter])
        else:
            for i in sorter:
                print("{0:7.3f}    {1}".format(float(orders[i]), names[i]))

        #stop()

    def unique_types(self):
        all_primitives = self.primitives()
        types = [p.type for p in all_primitives]
        return set(types)

    def __len__(self):
        return len(self.ListPrimitives())


#---------------------------------------------------------------------------------
class Primitive(object):
    """ A class for pipeline primitives

    Allows attribute-style access to various primitive properties

    Must be initialized using a primitive info dict returned from
    the GPIPrimitivesConfig object.
    """
    def __init__(self, infodict):
        for k in infodict.keys():
            self.__dict__[k.lower()] = infodict[k]

    @property
    def full_path(self):
        return os.path.join( GPIConfig().get_directory('GPI_DRP_DIR'), "primitives", self.idlfunc+".pro")

    def __repr__(self):
        return "<gpi_pipeline.Primitive: "+self.name+" >"

    #property
    def argument_names(self):
        return [item['Name'] for item in self.arguments]

    def get_idl_header(self, trim=False):
        header = []
        started=False
        for line in open(self.full_path).readlines():
            if started: 
                if line.startswith(";-"): 
                    return header

                if trim:
                    if 'PIPELINE' in line: continue
                    elif 'INPUT' in line: continue
                    elif 'OUTPUT' in line: continue
                    elif 'NAME' in line: continue

                header.append(line.lstrip().lstrip(';').rstrip() ) # remove leading whitespace and comment symobl ;
            elif not started and line.startswith(';+'):
                started=True
                continue
        else:
            return header

    def get_idl_code(self):
        return open(self.full_path).readlines()
            



    def doc(self):
        """ Output documentation for this primitive in RST format """


        header = self.get_idl_header()

        def get_header_setting(starttext):
            if header is None:
                print("No header for "+self.name+", "+self.idlfunc)
                return "Unknown"
            for line in header:
                if starttext in line:
                    return line.split(starttext)[1]
            else:
                return "Not specified"


        # figure out the output suffix
        suffix_code = [l.strip() for l in self.get_idl_code() if (('suffix' in l) and (not l.lstrip().startswith(';')) and ('=' in l))]
        if len(suffix_code) == 0:
            suffix=None
        elif len(suffix_code) == 1:
            suffix = suffix_code[0].split('=')[1]
        else:
            suffix = "Could not be determined automatically"

        output = ""
        output += self.name+"\n"
        output += "-"*len(self.name)+"\n\n"
        output += self.comment
        output += "\n\n**Category**: "+self.type+"      **Order**: "+str(self.order)
        output += "\n\n**Inputs**: "+get_header_setting('INPUTS:') 
        output += "\n\n**Outputs**: "+get_header_setting('OUTPUTS:')
        if suffix is not None:
            output += "      **Output Suffix**: "+suffix
        output += "\n\n**Notes**:\n\n.. code-block:: idl\n\n"
        output +=  "\n".join (self.get_idl_header(trim=True) )
        if len(self.arguments) >= 1:
            output += "\n\n**Parameters**:\n\n"
            head_line = ('Name', 'Type', 'Range', 'Default', 'Description')
            rows = [(a['Name'], a['Type'], a['Range'], a['Default'], a['Description']) for a in self.arguments]
            rows.insert(0, head_line)
            output+= self._toRSTtable(rows, header=True)

        output += "\n\n**IDL Filename**: "+self.idlfunc+".pro\n\n\n"

        return output

    def _toRSTtable(self, rows, header=True, vdelim="  ", padding=1, justify='right'):
        """ Outputs a list of lists as a Restructured Text Table

        - rows - list of lists
        - header - if True the first row is treated as a table header
        - vdelim - vertical delimiter betwee columns
        - padding - padding nr. of spaces are left around the longest element in the
          column
        - justify - may be left,center,right

        HISTORY:
        ---------
        From http://code.activestate.com/recipes/267662-table-indentation/
        originally by Attila Vasarhelyi with contributions from Danny G
        """
        import string

        border="=" # character for drawing the border
        justify = {'left':string.ljust,'center':string.center, 'right':string.rjust}[justify.lower()]

        # calculate column widths (longest item in each col
        # plus "padding" nr of spaces on both sides)
        #cols = zip(*rows)      #works if all rows are same len

        # to allow for different sized row lists :
        # add a blank space for every missing cell.
        cols = map(lambda *row: [elem or ' ' for elem in row], *rows)

        colWidths = [max([len(str(item))+2*padding for item in col]) for col in cols]

        # the horizontal border needed by rst
        borderline = vdelim.join([w*border for w in colWidths])

        # outputs table in rst format
        output = ""
        output += borderline+"\n"
        for row in rows:
            output += vdelim.join([justify(str(item),width) for (item,width) in zip(row,colWidths)])
            output += "\n"
            if header: output+=borderline+"\n"; header=False
        output += borderline+"\n"
        return output


#---------------------------------------------------------------------------------
class GPICalDB(object):
    """ Pythonic interface to GPI Calibrations DB

    Currently read-only; can't add new calibration files this way.

    """

    def __init__(self):
        import atpy
        conf = GPIConfig()
        
        self.caldir = conf.get_directory('GPI_CALIBRATIONS')
        self.dbfilename=self.caldir + os.sep + "GPI_Calibs_DB.fits"


        self._colnames_nice2raw = collections.OrderedDict( (
            ('Filename', "FILENAME"), ("Type", "TYPE"),
            ('Disperser',"PRISM"), ("IFS Filter", "FILTER"),
            ("Lyot Mask", "LYOT"), ('Apodizer', 'APODIZER'),
            ('ITime', "ITIME"), ("Date", "JD"), ("Julian Date", "JD")) )


        if not os.path.exists( self.dbfilename ):
            _log.warning('Cal DB filename does not exist: '+self.dbfilename+"!")
            pass
            #calibration DB dir is empty - create a blank one
        else:
            self.table = atpy.Table( self.dbfilename)
            # we have to rstrip every string type field, for some reason they have whitespace from IDL:
            for key in self.table.columns.keys:
                if 'S' in str(self.table[key].dtype):
                    self.table[key] = np.asarray( [a.rstrip() for a in self.table[key] ])

            _log.info("Loaded Calibrations DB from "+self.dbfilename)
            _log.info("    Found "+str(self.nfiles)+ " entries")

        self._verify()

    def _verify(self):
        """ Scan each file in DB and makes sure it exists. Drop any indexed files that do not exist. """
        valids = []
        fully_qualified_filenames = []
        for path, fn in zip(self.table['PATH'], self.table['FILENAME']):
            fully_qualified_filenames.append(os.path.join( path,fn))
            valids.append(os.path.exists(os.path.join( path,fn)))
        #print valids
        if False in valids:
            wvalid = np.where(valids)
            self.table.data = self.table.data[wvalid]
            fully_qualified_filenames = fully_qualified_filenames[wvalid]
            _log.warning("    Some indexed file(s) did not exist. Only "+str(self.nfiles)+" valid entries.")
        self.full_filenames = np.asarray(fully_qualified_filenames)

    def __getitem__(self,linenum):
        """ Return a tuple summarizing the properties of the Nth item in the database """

        rawline = self.table[linenum]
        output1 = [(keystring,rawline[self._colnames_nice2raw[keystring]]) for keystring in self._colnames_nice2raw.keys()]
        output = collections.OrderedDict(output1)

        output['Date'] = gpi_utils.jd2datestr(output['Date'], longform=True)

        return output.values()


    @property
    def nfiles(self):
        """ Number of files currently indexed in the calibration database."""
        return len(self.table)
    def __len__(self):
        return len(self.table)

    @property
    def unique_caltypes(self):
        """ Set of unique calibration type strings currently in database"""
        return set(self.table['TYPE'])
    @property
    def unique_dispersers(self):
        """ Set of unique DISPERSR keyword values currently in database"""
        return set(self.table['PRISM'])
    @property
    def unique_itimes(self):
        """ Set of unique ITIME keyword values currently in database"""
        return set(self.table['ITIME'])
    @property
    def unique_apodizers(self):
        """ Set of unique APODIZER keyword values currently in database"""
        return set(self.table['APODIZER'])
    @property
    def unique_readoutmodes(self):
        """ Set of unique READOUTMODE keyword values currently in database"""
        return set(self.table['READOUTMODE'])
    @property
    def unique_jds(self):
        """ Set of unique Julian Date (JD) values currently in database"""
        return set(self.table['JD'])

    @property
    def unique_datestrs(self):
        """ Set of unique date string YYMMDD values currently in database"""
        jds = self.unique_jds
        datestrs = [gpi_utils.jd2datestr(d) for d in jds]
        return set(datestrs)

    #--------
    def rescan_directory(self):
        """ Rescan a calibration directory to detect new/changed files"""
        raise NotImplementedError('not yet')

    def _cal_info_from_header(self):
        pass

    def get_best_cal(self, type_, fits_data, **kwargs):
        """ determine the best available cal of requested type for a given file

        Parameters
        ------------
        type_ : str
            String defining cal file type, e.g. "Dark", "Wavelength Solution", etc.
        fits_data : either a pyfits.HDUList or a gpidata.IFSdata object
            The file whose header(s) we're trying to match
        """
        pass

    def get_readoutmode_string(self, fits_data):
        """ format a string describing unambiguously the IFS readout mode. 

        """
        raise NotImplementedError("Not implemented yet")


    def select_calfiles(self, caltype=None, disperser=None, ifsfilt=None, apodizer=None,
            lyotmask=None, itime=None, readoutmode=None, datestr=None,
            return_mask=False, return_indices=False, return_objects=False):
        """ Return a list of one or more calibration files that match some chosen query criteria 
        
        Provide one or more selection criteria when calling. These will be combined with logical AND such that only 
        files matchings all critera are returned. 

        Parameters
        -----------
        various parameters for the different criteria

        return_mask : bool
            If true, return a mask array you can use for subscripting the index table
            instead of returning just the filenames.
        caltype : string
            Value for calibration type to select results
        disperser : string
            Value for DISPERSR keyword to select results
        ifsfilt : string
            Value for IFSFILT keyword to select results
        apodizer : string
            Value for APODIZER keyword to select results
        lyotmask : string
            Value for LYOTMASK keyword to select results
        itime : string
            Value for ITIME keyword to select results
        readoutmode: string
            Value for readout mode string to select results
        datestr : string
            Value for date string YYMMDD keyword to select results
               
        """

        mask = np.ones(len(self.table), dtype=bool)

        # simple matches on exact column values
        if caltype is not None: mask = mask & (self.table['TYPE'] == caltype)
        if disperser is not None: mask = mask & (self.table['PRISM'] == disperser)
        if ifsfilt is not None: mask = mask & (self.table['FILTER'] == ifsfilt)
        if apodizer is not None: mask = mask & (self.table['APODIZER'] == apodizer)
        if lyotmask is not None: mask = mask & (self.table['LYOT'] == lyotmask)
        if itime is not None: mask = mask & (self.table['ITIME'] == itime)
        if readoutmode is not None: mask = mask & (self.table['READOUTMODE'] == readoutmode)

        # matches on computed or derived values
        if datestr is not None: 
            datestrs = np.asarray([gpi_utils.jd2datestr(jd) for jd in self.table['JD']])
            mask = mask & (datestrs == datestr)

        # selection by ranges:
            
        if return_indices: return np.where(mask)[0]
        elif return_mask: return mask
        elif return_objects:
            import gpidata
            return [gpidata.IFSData(filename, loadData=False) for filename in self.full_filenames[mask]]
        else: return self.full_filenames[mask]


    def _describe(self, mask, sortby='itime'):
        w = np.where(mask)

        items = self.table[mask]

        items = self._sort(items, sortby)

        for line in items :
            gd = gpi_utils.jd2gd(line['JD'])
            print("%30s\t%20s\t%8s\t%6.1f\t%04d-%02d-%02d %2d:%02d" % (line['FILENAME'], line['TYPE'], line['PRISM']+"-"+line['FILTER'], line['ITIME'],gd[0], gd[1],gd[2],gd[3],gd[4]))

    def _sort(self, items, sortby='itime'):

        if sortby.lower() == 'itime':
            # sort by itime
            sorter = np.argsort(items['ITIME'])
            # for same itime, then sort by date, too...

        elif sortby.lower() == 'date' or sortby.lower() == 'jd':
            sorter = np.argsort(items['JD'])

        return items[sorter]


    def _list(self, mask, sortby='itime'):
        w = np.where(mask)
        items = self.table[mask]
        items = self._sort(items, sortby)
        names = []
        for line in items :
            names.append(line['FILENAME'])
        return names

    def report(self, by='type'):
        """ Display on screen a nicely formatted listing of calibration database contents

        Parameters
        ------------
        by : string
            How to organize/sort the report? You can organize by 'type' of file or by 'date'. 
        """


        if by.lower() =='type':
            print(" *****  Darks ***** ")
            self._describe( self.select_calfiles(caltype='Dark File', return_mask=True))
            print("")
            print(" *****  Wave Cals ***** ")
            self._describe( self.select_calfiles(caltype='Wavelength Solution Cal File', return_mask=True))
            print("")
            print(" ***** Flats ***** ")
            self._describe( self.select_calfiles(caltype='Flat field', return_mask=True))
        elif by.lower() == 'date':
            dates = list(self.unique_datestrs)
            dates.sort()
            for d in dates:
                print(" ***** %s ***** " % d)
                self._describe( self.select_calfiles(datestr=d, return_mask=True))
 

    def _closest_itime(self, requested_itime):
        itimes = np.asarray(list(self.unique_itimes))
        wm = np.argmin( np.abs(itimes-requested_itime))
        _log.info("Closest itime is "+str(itimes[wm]))
        return itimes[wm]

    def darks_by_time(self, itime=1.45):
        """ List all dark files for a given itime"""

        return db._list(db.select_calfiles(caltype='Dark File', 
            itime=self._closest_itime(itime), return_mask=True), sortby='jd')


#---------------------------------------------------------------------------------
class Recipe(object):
    """ Class for parsing/manipulating data reduction recipe files 


    >>> rec = Recipe('somefile.xml')

    """
    def __init__(self, filename=None, fromstring=None):
        """ Read a DRF from disk """
        if filename is not None:
            self.last_read_filename=filename
            self._tree = etree.parse(filename)
        elif fromstring is not None:
            self.last_read_filename=None
            self._tree = etree.fromstring(fromstring)
        else:
            raise ValueError("Must provide either a filename or an XML string to instantiate a Recipe")
            # TODO allow blank recipes somehow? Use etree factory class to instantiate a blank recipe maybe?
 


        self.last_saved_filename=None

        # Define an internal helper class
        class _primitive_accessor(object):
            """ Wrapper interface class to allow accessing the named
            <module> (AKA "primitive") elements in an XML Recipe file using a Dict-like interface. 

            Don't use this on it's own, use it as the primitive element of a DRF object: 

            mydrf = DRF('some_file.xml')
            mydrf.module['Dark subtraction']

            """

            def __init__(self, parent_Recipe):
                self._Recipe= parent_Recipe
            def __getitem__(self, i):

                if gpi_utils.is_integer(i):
                    # return the i'th module
                    modules = [mod for mod in self._Recipe._tree.getroot().iter('module')]
                    if i <0 or i> len(modules):
                        raise KeyError("Module index must be between 0 and %d for this DRF." %(len(modules)-1) )
                    else:
                        return modules[i]


                else:
                    # try using it as a string index like a hash

                    for mod in self._Recipe._tree.getroot().iter('module'):
                        if mod.attrib['name'].lower() == str(i).lower():
                            return mod
                    raise KeyError("No module found with name %s" % str(i))


    
        self.primitive= _primitive_accessor(self)


        #self.input_path =   self._xml_dataset.attrib['InputDir']
        #self.output_path =  self._xml_dataset.attrib['OutputDir']
        #self.name=          self._xml_dataset.attrib['Name']

    def __repr__(self):
        return "<gpi_pipeline.Recipe : %s, with %d files and %d primitives>" % (self.name, len(self.dataset), len(self.primitives_names))

    # there must be a more elegant way to build this with an iterative loop creating the
    # properties, but I'm ignoring it for now:

    @property
    def _xml_dataset(self): 
        dataset = list(self._tree.getroot().iter('dataset'))
        if len(dataset)> 1:
            raise ValueError('Invalid DRF: contains more than one <dataset> element')
        return dataset[0]

    @property
    def name(self): return self._tree.getroot().attrib['name']
    @name.setter
    def name(self,value):  self._tree.getroot().attrib['name'] = value 

    @property
    def reductionType(self): return self._tree.getroot().attrib['ReductionType']
    @reductionType.setter
    def reductionType(self,value):  self._tree.getroot().attrib['ReductionType'] = value 



    @property
    def input_path(self): return self._xml_dataset.attrib['InputDir']
    @input_path.setter
    def input_path(self,value):  self._xml_dataset.attrib['InputDir'] = value 
    @property
    def output_path(self): return self._xml_dataset.attrib['OutputDir']
    @output_path.setter
    def output_path(self,value):  self._xml_dataset.attrib['OutputDir'] = value 



    def write(self,filename, paths='relative', use_vars=False, check_exist=True):
        """ Write the Recipe to some location on disk

        Parameters
        ----------
        paths : string
            either 'relative' or 'absolute' to describe what sort of filenames should
            be used in the Recipe file
        use_vars : bool
            use environment variables in path names (True) or not (False)
        check_exist : bool
            check the dataset files all exist
        """
        if check_exist: self._check_dataset_exists()


        #todo implement relative/absolute paths
        if use_vars is True:
            # set the input/output dirs to use the environment variables.
            # note the existence of os.path.expandvars for string expansions, but
            # we need to write the inverse
            raise NotImplemented('not yet')


        self._tree.write(filename)
        self.last_saved_filename=filename

    def Queue(self, filename=None):
        """ Add the Recipe to the reduction queue. 
        Will use last saved filename, if available, otherwise will use a
        default name. """


        if filename is not None: 
            queue_fn = filename
        else:
            if self.last_saved_filename is not None:
                queue_fn = os.path.basename(self.last_saved_filename)
            elif self.last_read_filename is not None:
                queue_fn = os.path.basename(self.last_read_filename)
            else:
                raise NotImplementedError('need to save first to set filename')

        queue_fn = gpi_utils.expand_path(GPIConfig().get_directory('GPI_DRP_QUEUE_DIR')+os.sep+ queue_fn )

        if not queue_fn.endswith(".waiting.xml"):
            # FIXME more careful handling here...
            queue_fn = os.path.splitext(queue_fn)[0] + ".waiting.xml"

        self._tree.write(queue_fn)
        _log.info("Queued to %s" % queue_fn)



    # modules
    @property
    def primitives_names(self):
        "get the Primitive names from the XML"
        modules = list(self._tree.getroot().iter('module'))
        if len(modules)== 0:
            raise ValueError('Invalid DRF: contains no <module> elements')
        module_names = [el.attrib['name'] for el in list(modules)]
        return module_names 
    @primitives_names.setter
    def set_primitives(self, value):
        raise NotImplemented('Cannot set primitives')

    # 



    # dataset
    @property
    def dataset(self):
        "get the FITS filenames in a dataset from the XML Recipe file"
        #try:
        if 1:
            # this is a cleaner way but ends up being case sensitive:
            #dataset_files = [el.attrib['FileName'] for el in list(self._xml_dataset) if el.attrib['FileName'] != '']
            # this way is more of a pain but ends up being case insensitive:
            # it also enforces only looking at <fits> elements...
            dataset_files=[]
            for el in list(self._xml_dataset):
                if el.tag.lower() == 'fits':
                    key0 = el.keys()[0].lower()
                    if key0 == 'filename': 
                        dataset_files.append(el.attrib[el.keys()[0]])

        #except:
            #dataset_files = None
        return dataset_files
    @dataset.setter
    def dataset(self, filelist):
        """ Set the dataset to contain these files """
        if not hasattr(filelist,'__iter__'):
            raise ValueError("You can only set a dataset to some iterable (list, tuple, etc) set of filenames")
        
        # remove all existing items
        for el in self._xml_dataset: self._xml_dataset.remove(el)

        path0=None
        #add in the new ones
        for fn in filelist:
            path, fname = os.path.split(os.path.abspath(fn))
            if path0 is None: path0 = path
            elif path != path0:
                _log.error( "ERROR: files in dataset are not all in same directory")
            el = etree.SubElement(self._xml_dataset, "fits", FileName=fname)
        if path0 is None or filelist == []: 
            _log.warning("Setting dataset to a null dataset containing no files.")
            path0 = '' # for the case of a blank dataset
        self.input_path = path0
        


    def _check_dataset_exists(self):
        """ Check that all the FITS files listed in the dataset exist """
        error_filenames=[]
        for fn in self.dataset:
            if not os.path.exists(gpi_utils.expand_path(self.input_path+os.sep+fn)) : error_filenames.append(fn)
        if error_filenames != []:
            raise IOError("File(s) %s do not exist." % ", ".join(error_filenames))

    
    def __str__(self):
        return etree.tostring(self._tree)


#---------------------------------------------------------------------------------
def RecipeFromTemplate(templatename, datafiles=None):
    """ Create a Recipe object from a template

    Parameters
    ------------
    templatename : string
        Descriptive name (**not** filename) for the template
    dataset : iterable
        Some list of filenames to apply this template to.


    """
    import glob
    C = GPIConfig()
    templatedir = C['GPI_DRP_TEMPLATES_DIR']
    templates = glob.glob(templatedir+os.sep+"*.xml")
    result = None
    for filename in templates:
        r = Recipe(filename)
        if r.name.lower() == templatename.lower():
            result = r
            break
    if result is None:
        raise KeyError("Could not find any template named '%s'." % templatename)

    if datafiles is not None:
        result.dataset = datafiles
    return result


def ListRecipeTemplates(ReductionType=None):
    """ Return a list of all available recipe templates """
    import glob
    C = GPIConfig()
    templatedir = C['GPI_DRP_TEMPLATES_DIR']
    templates = glob.glob(templatedir+os.sep+"*.xml")

    if ReductionType is None:
        # return all
        out = []
        for filename in templates:
            rr = Recipe(filename)
            if hasattr(rr, 'name'): 
                out.append(rr.name)
            else:
                out.append('-no name')
        return out
        #return [Recipe(filename).name for filename in templates]
    else:
        return [Recipe(filename).name for filename in templates if Recipe(filename).reductionType.lower() == ReductionType.lower()]


def validate_templates():
    """ Validate if templates all contain valid primitives """
    import glob
    C = GPIConfig()
    PC = GPIPrimitivesConfig()
    templatedir = C['GPI_DRP_TEMPLATES_DIR']
    templates = glob.glob(templatedir+os.sep+"*.xml")

    avail_primitives = PC.ListPrimitives()

    for template_filename in templates:
        rec = Recipe(template_filename)
        print("Validating Recipe: "+rec.name)

        for prim in rec.primitives_names:
            if prim not in avail_primitives:
                print("    In template: "+template_filename)
                print("    WARNING: primitive not found for '{0}' ".format(prim))


            print(prim.argument_names)


 

#--------------------------------------------------------------------------------
class PipelineDriver(object):
    """ Object class for driving the GPI DRP """

    def __init__(self, indir = None, outdir=None):

        if indir is None: indir=gpi_utils.expand_path(os.getenv('GPI_RAW_DATA_DIR'))
        if outdir is None: outdir=gpi_utils.expand_path(os.getenv('GPI_DRP_OUTPUT_DIR'))
            

        if not indir.endswith(os.sep):
            indir += os.sep
        self.indir = indir
        self.outdir=outdir
        self.dropcounter=1 

        self.scan_templates()

    def scan_templates(self):
        """Scan all the templates in the templates directory into memory """
        templatedir = gpi_utils.expand_path(os.getenv('GPI_DRF_TEMPLATES_DIR'))
        templatefiles = glob.glob(templatedir+os.sep+"*.xml")
        templates = {}

        for tfile in templatefiles:
            tdrf = DRF(tfile)
            _log.info("  Loaded template =%s, %s" % (tfile, tdrf.name))
            templates[tdrf.name] = tdrf

        _log.info("Loaded %d templates" % len(templates.keys()))
        self.templates=templates

    def apply_template(self, templatename, datafiles):
        """ Create a new DRF based on the specified template. """
        drf = deepcopy(self.templates[templatename])
        drf.input_path = self.indir
        drf.output_path = self.outdir
        drf.dataset = datafiles
        return drf

    def parse_directory(self, path):
        import fitshdrs

        filenam

        headers = fitshdrs.FITSHeaderArray(path+os.sep+"*.fits")

        # darks
        dark_files = []
        self.apply_template('Dark template', dark_files).write('path')

        for filter in ['Y','J','H','K1','K2']:
            for prism in ['PRISM','WOLLASTON']:
                pass

   

if __name__ == "__main__":
        
    logging.basicConfig(level=logging.INFO)

    pass


#---------------------------------------------------------------------------------


if __name__ == "__main__":
    if 0:

        import argparse
        import pidly
        parser = argparse.ArgumentParser(description='Command line args for gpi-pipeline.py, a Python wrapper for the GPI Data Pipeline (in IDL)')
        parser.add_argument('--single', help='Process one single reduction recipe (DRF), given as a filename argument after this switch, and then exit.', default='')
        parser.add_argument('--nogui', help='Do not display the "Status Console" GUI window.', action='store_const', const=True, default=False)
        parser.add_argument('--rescanCalDB', help='Rescan and re-index Calibration DB directory on startup', action='store_const', const=True, default=False)
        parser.add_argument('--verbose', help='Display more text output than usual, mostly for debugging', action='store_const', const=True, default=False)

        res = parser.parse_args()

        idl = pidly.IDL()
        try:
            idl.pro('gpipiperun', nogui=res.nogui, single=res.single, verbose=res.verbose, rescan=res.rescanCalDB)
        except IOError:
            print("IDL session exited")
    else:
        logging.basicConfig(level=logging.INFO)
        C = GPIConfig()
        db = GPICalDB()

