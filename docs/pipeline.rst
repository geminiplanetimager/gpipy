
Interacting with the Pipeline, Config, and Recipes 
#####################################################


This documentation is very incomplete, but see the doc strings of the various
functions and objects mentioned below for more details. 


Pipeline Configuration Data
------------------------------

It's often useful in Python to be able to access the same configuration data that's used by the pipeline in IDL.
The ``GPIConfig()`` class lets you do this easily.


Settings are accessible by using an instance of ``GPIConfig`` as a dict.
This works for both config settings and directory names.

>>> gpiconfig = gpipy.GPIConfig()

>>> result = gpiconfig['memsripple']
>>> result = gpiconfig['GPI_RAW_DATA_DIR']

>>> dirpath = gpiconfig.get_directory_name('GPI_CALIBRATIONS_DIR')

These settings are read only; one must change the pipeline configuration files
directly or adjust environment variables to modify settings.


Similarly, the ``GPIPrimitivesConfig`` class provides access to the 
index of primitives and their parameters as stored in the ``gpi_pipeline_primitives.xml`` file.


>>> pc = GPIPrimitivesConfig()
>>> pc.ListPrimitives()
[... list of available primitives is displayed...
>>> pc.primitive['Accumulate Images']
[... returns a dict of information about that primitive...]



Calibrations Database
-------------------------


The ``GPICalDB()`` class provides a Python interface to the
contents of your pipeline calibration database. 

Note, this hasn't been updated in a while and is probably out of date...



Working with Recipes
---------------------

The ``Recipe`` class provides a Pythonic object interface to
create and edit GPI reduction recipes. 

The ``RecipeFromTemplate()`` function is a helper to instantiate a 
Recipe from a named recipe template and a list of FITS files to apply it to. 
``ListRecipeTemplates()`` will give you a list of the available templates in case you need it. 



Controlling the Pipeline
----------------------------

The ``PipelineDriver`` class provides a basic interface for controlling the pipeline in IDL. 

Its ``run_recipe`` method provides a simple way to run a single recipe in the pipeline, including
returning to Python information about the execution status (success/failure) and any output FITS files
from that recipe. 

