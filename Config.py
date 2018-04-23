import json

with open('config.json', 'r') as f:
    config = json.load(f)

global gama
gama = config['CONSTANTS']['GAMA']

global cfl_coefficient
cfl_coefficient = config['CONSTANTS']['CFL_COEFFICIENT']

global parallel
parallel = config['CONSTANTS']['PARALLEL']
