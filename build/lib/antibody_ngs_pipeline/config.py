import os
from .horseshoe import abstar_params, mongo_params, print_the_splash

print_the_splash()

PROJECT             = input('Basespace Project: ')
ABSTAR_ARGS         = abstar_params(PROJECT)
MONGO_ARGS          = mongo_params(PROJECT, ABSTAR_ARGS)