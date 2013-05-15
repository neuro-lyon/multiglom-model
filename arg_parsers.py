# -*- coding:utf-8 -*-
from argparse import ArgumentParser

SIM_PARSER = ArgumentParser(description="Run a multi-glomerular simulation.")
SIM_PARSER.add_argument('psfile')
SIM_PARSER.add_argument('--no-plot', action='store_true')
SIM_PARSER.add_argument('--no-summary', action='store_true')
SIM_PARSER.add_argument('--full-ps', action='store_true')
SIM_PARSER.add_argument('--no-brian-output', action='store_true')

MULTISIM_PARSER = ArgumentParser(description="Run simulations in parallel.")
MULTISIM_PARSER.add_argument('nproc', type=int)
MULTISIM_PARSER.add_argument('pset_dir')
MULTISIM_PARSER.add_argument('--testing', action='store_true')
