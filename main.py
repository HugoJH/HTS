#!/usr/bin/env python
from os.path import basename, splitext
from utils import exists, md5
from GATK import *
from BAM import *
from VCF import *
from bwa import bwaMem

if __name__ == '__main__':
    from argparse import ArgumentParser
    
    parser = ArgumentParser()
    parser.add_argument("-t", "--tumor", dest="tumor",
                                help="tumor bam file")
    parser.add_argument("-n", "--normal", dest="normal",
                                help="normal bam file")
    parser.add_argument("-o", "--output", dest="output",
                                help="output vcf file")
    parser.add_argument("-tfqs", "--tumor-fastqs", dest="tumor_fastqs",
                                help="tumor fastq files")
    parser.add_argument("-nfqs", "--normal-fastqs", dest="normal_fastqs",
                                help="normal fastq files")
    parser.add_argument("-gr", "--germline-resources", dest="germline_resources",
                                help="germline resource file")
    parser.add_argument("-tmp", "--temporary-directory", dest="results_folder", help="directory where results are stored", default="/data/tmp/")
    
    parser.add_argument("-r", "--reference", dest="reference", help="genome reference")
    
    parser.add_argument("-il", "--interval-list", dest="interval_list", help="interval list file")
    parser.add_argument("-pon", "--panel-of-normals", dest="pon", help="panel of normals file",default="/data/resources/PoN_broadinstitute/384.PoN.GDC.vcf.gz")

    parser.add_argument("-ks", "--known-sites", dest="known_sites", help="known sites files", default="/data/resources/germline_resources/1000G_phase1.snps.high_confidence.b37.vcf /data/resources/germline_resources/dbsnp_138.b37.vcf /data/resources/germline_resources/Mills_and_1000G_gold_standard.indels.b37.vcf")
    args = parser.parse_args()
    if args.tumor_fastqs:
        tumor_fastqs=args.tumor_fastqs.split()
        normal_fastqs=args.normal_fastqs.split()

    reference = args.reference
    interval_list = args.interval_list
    known_sites = args.known_sites
    tumor_bam_output = ""
    normal_bam_output = ""
    results_folder = args.results_folder
    if args.tumor_fastqs:
        tumor_bam_output = "tumor_{fastq_id}.bam".format(fastq_id=md5(tumor_fastqs[0]))
        normal_bam_output = "normal_{fastq_id}.bam".format(fastq_id=md5(normal_fastqs[0]))
        tumor_ubam_output = "tumor_ubam_{fastq_id}.bam".format(fastq_id=md5(tumor_fastqs[0]))
        normal_ubam_output = "normal_ubam_{fastq_id}.bam".format(fastq_id=md5(normal_fastqs[0]))
        FASTQS_to_UBAM(tumor_fastqs[0], tumor_fastqs[1], tumor_ubam_output, "TumorDefaultGroupName", "TumorDefaultSampleName")
        FASTQS_to_UBAM(normal_fastqs[0], normal_fastqs[1], normal_ubam_output, "NormalDefaultGroupName", "NormalDefaultSampleName")
    else:
        tumor_bam_output = results_folder + '/' +  basename(splitext(args.tumor)[0]) + ".bam.output.bam"
        normal_bam_output = results_folder + '/' + basename(splitext(args.normal)[0]) + ".bam.output.bam"    
        tumor_ubam_output = results_folder + '/' +  basename(splitext(args.tumor)[0]) + ".ubam.output.bam"
        normal_ubam_output = results_folder + '/' + basename(splitext(args.normal)[0]) + ".ubam.output.bam"    
        BAM_to_UBAM(args.tumor, tumor_ubam_output, results_folder)
        BAM_to_UBAM(args.normal, normal_ubam_output, results_folder)
    
    UBAM_to_BAM(tumor_ubam_output, tumor_bam_output, reference, interval_list, known_sites, results_folder)
    UBAM_to_BAM(normal_ubam_output, normal_bam_output, reference, interval_list, known_sites, results_folder)
    
    BAM_to_VCF(tumor_bam_output, normal_bam_output, args.output, args.germline_resources, reference, args.interval_list, args.pon)


