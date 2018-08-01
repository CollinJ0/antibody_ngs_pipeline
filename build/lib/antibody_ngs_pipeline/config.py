import os
import .horseshoe as scripps

scripps.print_the_splash()

PROJECT             = input('Basespace Project: ')
ABSTAR_ARGS         = scripps.abstar_params(PROJECT)
MONGO_ARGS          = scripps.mongo_params(PROJECT, ABSTAR_ARGS)