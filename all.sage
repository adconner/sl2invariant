from IPython import embed
from cProfile import run

import resource
import sys
resource.setrlimit(resource.RLIMIT_STACK, [0x10000000, resource.RLIM_INFINITY])
sys.setrecursionlimit(0x100000)

load_attach_mode(load_debug=True)

from itertools import *

attach('basis.py')
attach('tensor_utils.sage')

# vim: ft=python
