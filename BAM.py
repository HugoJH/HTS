#!/usr/bin/env python
#from collections import defaultdict
#from subprocess import check_output
#import os.path as path
#import multiprocessing
#import os
#import sys
#import datetime
#import hashlib
#import inspect
#from inspect import signature
#import inspect
#from pathlib import Path
from utils import exists, md5
from GATK import *
from bwa import bwaMem

def UBAM_to_BAM(input_ubam, output, reference, interval_list, known_sites, tmp_dir):
    file_id = md5(input_bam)
    markAdapters(input_bam, f"{tmp_dir}/MarkedAdapters_{file_id}.bam", f"{tmp_dir}/MarkedAdapters_{file_id}.metrics.txt","{tmp_dir}")

    samToFastq(f"{tmp_dir}/MarkedAdapters_{file_id}.bam", f"{tmp_dir}/SamToFastq_{file_id}.fq", "{tmp_dir}")

    bwaMem(f"{tmp_dir}/SamToFastq_{file_id}.fq", f"{tmp_dir}/bwaMem_{file_id}.bam", reference)

    mergeBamAlignment(f"{tmp_dir}/MarkedAdapters_{file_id}.bam", f"{tmp_folder}/bwaMem_{file_id}.bam", f"{tmp_dir}/mergeBamAlignment_{file_id}.bam", reference, "{tmp_dir}")

    markDuplicates(f"{tmp_dir}/mergeBamAlignment_{file_id}.bam", f"{tmp_folder}/markDuplicates_{file_id}.bam", f"{tmp_folder}/markduplicates_{file_id}.bam.txt", "{tmp_dir}")

    BaseRecalibrator(f"{tmp_dir}/markDuplicates_{file_id}.bam", f"{tmp_folder}/BaseRecalibrator_{file_id}.bam", reference, interval_list, known_sites)

    applyBQSR(f"{tmp_dir}/markDuplicates_{file_id}.bam", output, f"{tmp_folder}/BaseRecalibrator_{file_id}.bam", reference)

    try:
        validateSamFile(output, f"{tmp_dir}/ValidateSamFile_{file_id}.out", reference)
    except:
        print(f"{tmp_dir}/applyBQSR_{file_id}.bam didn't pass validation")

    os.rename(f"{tmp_dir}/applyBQSR_{file_id}.bam", output)

def BAM_to_UBAM(input_bam, output, tmp_dir):
    file_id = md5(input_bam)
    revertSam(input_bam, output, tmp_dir)

def FASTQS_to_UBAM(fastq1, fastq2, output, readgroup_name, sample_name):
    file_id = md5(input_bam)
    fastqToSam(fastq1, fastq2, output, readgroup_name, sample_name)
