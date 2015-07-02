#!/usr/bin/env python
import sys
import logging
logging.basicConfig(filename='error.log', level=logging.INFO,
    format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %H:%M:%S')

try:
    from pyWAD.base import pluginWrapper
    pluginWrapper(sys.argv[1])
except:
    logging.exception('Exception on main handler:')
    raise
