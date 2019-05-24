#!/usr/bin/env python
from subprocess import check_output
import multiprocessing
import utils

def index_bam(bam,cores=multiprocessing.cpu_count()):
    check_output(["samtools index -@{cores} {bam}".format(bam=bam,cores=cores) ], shell=True)

