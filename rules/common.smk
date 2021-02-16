#Snakemake common functions across the pipeline
# 
#
# 18/11/2019
#
# Author: massimiliano [dot] Cocca [at] burlo [dot] trieste [dot] it
import pandas as pd
import ntpath
import collections
from snakemake.utils import validate
from snakemake.utils import min_version


#Define the config file
# configfile: "config.yaml"
# configfile: "annotation.yaml"
#We need to validate the config file we use, to be sure we are not missing any parameter
# validate(config, schema="../schemas/config.schema.yaml")
#
#
# We will need to fill this section with other stuff useful for the whole pipeline
#
#
#

	
