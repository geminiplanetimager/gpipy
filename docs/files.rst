
Working with GPI Files
#################################

Reading in GPI files
--------------------------


This is simple::

        import gpipy
        cube = gpipy.read("mygpifilename.fits")

The ``read()`` function will examine the FITS file header, and 
return an instance of the most appropriate subclass, for example
an ``IFSSpectralCube`` or ``IFSStokesCube`` instance. 

Each of these
classes has various attributes appropriate for the type of
data. For instance a 



There are subclasses for
2D and 3D GPI files, polarimetry data, wavelength calibrations. What you can do with a gpidata object depends on the subclass.
At a minimum you can access *the most commonly used* keywords as attributes:

    >>> dat.filter
    'H'
    >>> dat.apodizer
    'H_apod'
    etc.


Working with many files at once
--------------------------------

There is also a ``DataCollection`` class that lets you manipulate and query large groups of FITS files at once:

    >>> dc = gpidata.DataCollection('S20130909S*.fits')
    >>> dc.itimes
    [numpy ndarray of all the itimes for all those files]


The ``DataCollection.ParseToRecipes()`` function replicates much of
the functionality of the GPI Data Pipeline's Data Parser. It can generate
recipes to reduce your data. However, this functionality is not entirely 
kept in sync with the functionality in IDL, so its use is only recommended in
specific circumstances if you know what you're doing. 



The ``logsheet()`` function prints out a nicely formatted simple text log sheet describing all files
in a given directory. ``logsheet2()`` is a variant with different columns and formatting.

