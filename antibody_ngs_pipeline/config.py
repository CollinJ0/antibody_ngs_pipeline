import os
from .horseshoe import abstar_params, mongo_params, print_the_splash, parse_arguments



ARGS                = parse_arguments()
print_the_splash()
PROJECT             = input('Basespace Project: ')
ABSTAR_ARGS         = abstar_params(PROJECT)
MONGO_ARGS          = mongo_params(PROJECT, ABSTAR_ARGS)