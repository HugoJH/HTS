#!/usr/bin/env python
#from collections import defaultdict
#from subprocess import check_output
#import os.path as path
#import multiprocessing
#import os
#import sys
#import datetime
import hashlib
#import inspect
from inspect import signature
#import inspect
from pathlib import Path

def build_args_from_dict(args_dict, is_single_dashed=True):
    arguments = ""
    dashes = "-" if is_single_dashed else "--"
    for k, v in args_dict.items():
        if type(v) is list:
            arguments += " ".join([(" " + dashes + k + " " + i) for i in v])
        else:
            arguments += " " + dashes + k + " " + v
    return arguments

def function_args_dict(function,arguments):
    sig = signature(function)
    return sig.bind(*arguments)
#def reference_path(ref):
#    if ref == "b37":
#        return "/data/resources/genome_references/b37/b37.fa.gz"

def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def exists(fname):
    def wrapper(*args,**kwargs):
        output = function_args_dict(fname, args).arguments["output"]
        if not Path(output).exists():
            fname(**function_args_dict(fname, args).arguments)
        print(output)
    return wrapper
