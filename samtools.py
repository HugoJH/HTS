#!/usr/bin/env python
from subprocess import check_output
import multiprocessing
import utils

def index_bam(bam, output, cores=multiprocessing.cpu_count()):
    check_output(["samtools index -@{cores} {bam} {output}".format(bam=bam,output=output,cores=cores) ], shell=True)

