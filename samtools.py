#!/usr/bin/env python
#from collections import defaultdict
from subprocess import check_output
#import os.path as path
import multiprocessing
#import os
#import sys
#import datetime
#import hashlib
#import inspect
#from inspect import signature
#import inspect
#from pathlib import Path
import utils

def index_bam(bam,cores=multiprocessing.cpu_count()):
    check_output(["samtools index -@{cores} {bam}".format(bam=bam,cores=cores) ], shell=True)

