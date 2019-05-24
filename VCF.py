#!/usr/bin/env python
from subprocess import check_output
import os.path as path
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
    if not path.isfile(tumor_bam + ".bai"):
        index_bam(tumor_bam)
    if not path.isfile(normal_bam + ".bai"):
        index_bam(normal_bam)
    mutect(tumor_bam, normal_bam, f"{workdir}/Mutect2_{file_id}.vcf", germline_resource, reference, interval_list, pon, ip)
    learnReadOrientationModel(f"{workdir}/Mutect2_{file_id}.vcf.f1r2.tar.gz",f"{workdir}/{file_id}-artifact-prior.tsv.tar.gz")
    GetPileupSummaries(tumor_bam, f"{workdir}/GetPileupSummaries_{file_id}.table" , small_exac_vcf, interval_list, ip="50")
    GetPileupSummaries(normal_bam, f"{workdir}/GetPileupSummaries_normal_{file_id}.table" , small_exac_vcf, interval_list, ip="50")
    CalculateContamination(f"{workdir}/GetPileupSummaries_{file_id}.table", f"{workdir}/GetPileupSummaries_normal_{file_id}.table", f"{workdir}/CalculateContamination_{file_id}.table", f"{workdir}/Segments_{file_id}.tsv")
    FilterMutectCalls(f"{workdir}/Mutect2_{file_id}.vcf", f"{workdir}/FilterMutectCalls_{file_id}.vcf", reference, f"{workdir}/CalculateContamination_{file_id}.table",f"{workdir}/Segments_{file_id}.tsv",f"{workdir}/{file_id}-artifact-prior.tsv.tar.gz")
    SelectVariants(f"{workdir}/FilterMutectCalls_{file_id}.vcf", output, reference)

