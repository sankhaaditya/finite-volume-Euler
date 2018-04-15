import json

with open('config.json', 'r') as f:
    config = json.load(f)

global gama
gama = config['GAMA']

global cfl_coefficient
cfl_coefficient = config['CFL_COEFFICIENT']
