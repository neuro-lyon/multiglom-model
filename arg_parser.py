# -*- coding:utf-8 -*-
import argparse

APARSER = argparse.ArgumentParser(description="Run a multi-glomerular simulation.")
APARSER.add_argument('psfile')
APARSER.add_argument('--no-plot', action='store_true')
APARSER.add_argument('--no-indexes', action='store_true')
APARSER.add_argument('--full-ps', action='store_true')
