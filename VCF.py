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
from BAM import *
from bwa import bwaMem
from samtools import index_bam

@exists
def BAM_to_VCF(tumor_bam, normal_bam, output, germline_resource, reference, interval_list="", pon="", ip="50"):
    small_exac_vcf = "/home/hugo/MN4SHARED/shared/vainas/Path4/small_exac_common_3_b37.vcf.gz"
    workdir = path.dirname(output)
    file_id = getSampleName(tumor_bam,"/tmp/variant.id")
    if not ios.path.isfile(tumor_bam + ".bai"):
        index_bam(tumor_bam)
    if not ios.path.isfile(normal_bam + ".bai"):
        index_bam(normal_bam)
    mutect(tumor_bam, normal_bam, f"{workdir}/Mutect2_{file_id}.vcf", germline_resource, default_af, reference, interval_list, pon, ip)
    learnReadOrientationModel(f"{workdir}/Mutect2_{file_id}.vcf.f1r2.tar.gz",f"{workdir}/{file_id}-artifact-prior.tsv")
    GetPileupSummaries(tumor_bam, f"{workdir}/GetPileupSummaries_{file_id}.table" , small_exac_vcf, interval_list, ip="50")
    GetPileupSummaries(normal_bam, f"{workdir}/GetPileupSummaries_normal_{file_id}.table" , small_exac_vcf, interval_list, ip="50")
    CalculateContamination(f"{workdir}/GetPileupSummaries_{file_id}.table", f"{workdir}/GetPileupSummaries_normal_{file_id}.table", f"{workdir}/CalculateContamination_{file_id}.table", f"{workdir}/Segments_{file_id}.tsv")
    FilterMutectCalls(f"{workdir}/Mutect2_{file_id}.vcf", f"{workdir}/FilterMutectCalls_{file_id}.vcf", reference, f"{workdir}/CalculateContamination_{file_id}.table",f"{workdir}/Segments_{file_id}.tsv")
    CollectSequencingArtifactMetrics(tumor_bam, f"{workdir}/CollectSequencingArtifactMetrics_{file_id}.artifacts", reference, extension=".txt")
    FilterByOrientationBias(f"{workdir}/FilterMutectCalls_{file_id}.vcf", f"{workdir}/FilterByOrientationBias_{file_id}.vcf", f"{workdir}/CollectSequencingArtifactMetrics_{file_id}.artifacts.pre_adapter_detail_metrics.txt")
    SelectVariants(f"{workdir}/FilterByOrientationBias_{file_id}.vcf", output, reference)

