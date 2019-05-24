#!/usr/bin/env python
from collections import defaultdict
from subprocess import check_output
from utils  import exists, build_args_from_dict
from datetime import datetime
@exists
def fastqToSam(fastq1, fastq2, output, read_group_name, sample_name, library_name="DefaultLibraryName", platform_unit="DefaultPlatformUnit", platform="DefaultPlatform", sequencing_center="DefaultSequencingCenter", run_date=datetime.now().isoformat()):
    args = defaultdict(list)
    args["FASTQ"] = fastq1
    args["FASTQ2"] = fastq2
    args["OUTPUT"] = output
    args["READ_GROUP_NAME"] = read_group_name
    args["SAMPLE_NAME"] = sample_name
    args["LIBRARY_NAME"] = library_name
    args["PLATFORM_UNIT"] = platform_unit
    args["PLATFORM"] = platform
    args["SEQUENCING_CENTER"] = sequencing_center
    args["RUN_DATE"] = run_date
    check_output(["gatk FastqToSam " + build_args_from_dict(args)], shell=True)
    return output

@exists
def revertSam(input_bam, output, tmp_dir):
    args = defaultdict(list)
    args["I"] = input_bam
    args["O"] = output
    args["SANITIZE"] = "true"
    args["MAX_DISCARD_FRACTION"] = "0.005"
    args["ATTRIBUTE_TO_CLEAR"].append("XA")
    args["ATTRIBUTE_TO_CLEAR"].append("BD")
    args["ATTRIBUTE_TO_CLEAR"].append("XS")
    args["ATTRIBUTE_TO_CLEAR"].append("BI")
    args["SORT_ORDER"] = "queryname"
    args["RESTORE_ORIGINAL_QUALITIES"] = "true"
    args["REMOVE_DUPLICATE_INFORMATION"] = "true"
    args["REMOVE_ALIGNMENT_INFORMATION"] = "true"
    args["TMP_DIR"]= tmp_dir
    check_output(["gatk RevertSam" + build_args_from_dict(args)], shell=True)
    return output

@exists
def markAdapters(input_bam, output, metrics_file, tmp_dir):
    args = defaultdict(list)
    args["I"] = input_bam
    args["O"] = output
    args["M"] = metrics_file
    args["TMP_DIR"] = tmp_dir
    check_output(["gatk MarkIlluminaAdapters " + build_args_from_dict(args)], shell=True)
    return output

@exists
def samToFastq(input_bam, output, tmp_dir):
    args = defaultdict(list)
    args["I"] = input_bam
    args["FASTQ"] = output
    args["CLIPPING_ATTRIBUTE"] = "XT"
    args["CLIPPING_ACTION"] = "2"
    args["INTERLEAVE"] = "true"
    args["TMP_DIR"] = tmp_dir
    check_output(["gatk SamToFastq " + build_args_from_dict(args)], shell=True)
    return output

@exists
def mergeBamAlignment(input_unmapped_bam, input_aligned_bam, output, reference, tmp_dir):
    args = defaultdict(list)
    args["R"] = reference
    args["UNMAPPED_BAM"] = input_unmapped_bam
    args["ALIGNED_BAM"] = input_aligned_bam
    args["O"] = output
    args["CREATE_INDEX"] = "true"
    args["ADD_MATE_CIGAR"] = "true"
    args["CLIP_ADAPTERS"] = "true"
    args["CLIP_OVERLAPPING_READS"] = "true"
    args["INCLUDE_SECONDARY_ALIGNMENTS"] = "true"
    args["MAX_INSERTIONS_OR_DELETIONS"] = "-1"
    args["PRIMARY_ALIGNMENT_STRATEGY"] = "MostDistant"
    args["ATTRIBUTES_TO_RETAIN"] = "XS"
    args["TMP_DIR"] = tmp_dir
    check_output(["gatk MergeBamAlignment" + build_args_from_dict(args)], shell=True)
    return output


@exists
def markDuplicates(input_bam, output, metrics_file, tmp_dir):
    args = defaultdict(list)
    args["METRICS_FILE"] = metrics_file
    args["O"] = output
    args["CREATE_INDEX"] = "true"
    args["OPTICAL_DUPLICATE_PIXEL_DISTANCE"] = "100"
    args["TMP_DIR"] = tmp_dir
    args["INPUT"] = input_bam
    check_output(["gatk MarkDuplicates" + build_args_from_dict(args)], shell=True)

@exists
def BaseRecalibrator(input_bam, output, reference, interval_list="", known_sites="", interval_padding="50"):
    args = defaultdict(list)
    args["R"] = reference
    args["I"] = input_bam 
    args["O"] = output
    if interval_list:
        args["L"] = interval_list
    args["ip"] = interval_padding
    if known_sites:
        for known_site in known_sites.split():
            args["-known-sites"].append(known_site)

    check_output(["gatk BaseRecalibrator" + build_args_from_dict(args)], shell=True)
    return output

@exists
def applyBQSR(input_bam, output, recalibration_table, reference):
    args = defaultdict(list)
    args["R"] = reference
    args["I"] = input_bam 
    args["O"] = output
    args["bqsr"] = recalibration_table
    check_output(["gatk ApplyBQSR" + build_args_from_dict(args)], shell=True)
    return output

def validateSamFile(input_bam, output, reference):
    if not Path(output):
        args = defaultdict(list)
        args["R"] = reference
        args["I"] = input_bam
        args["O"] = output
        args["MODE"] = "SUMMARY"
        check_output(["gatk ValidateSamFile" + build_args_from_dict(args)], shell=True)
    else:
        with open(output,"r") as f:
            print("Validation results",f.read())

def getSampleName(input_bam, output):
    args = defaultdict(list)
    args["I"] = input_bam 
    args["O"] = output
    check_output(["gatk GetSampleName" + build_args_from_dict(args)], shell=True)
    samplename = ""
    with		open(output, 'r') as f:
        samplename = f.read()
    return samplename

@exists
def mutect(tumor_bam, normal_bam, output, germline_resource, reference, interval_list="",pon="", ip="50"):
    args = defaultdict(list)
    args["R"] = reference
    args["I"].append(normal_bam)
    args["I"].append(tumor_bam)
    args["O"] = output
    args["normal"] = getSampleName (normal_bam, "/tmp/normal.name")
    args["-f1r2-tar-gz"] = output + ".f1r2.tar.gz"
    args["bamout"] = f"{output}.bamout.bam"
    args["-germline-resource"] = germline_resource
    if pon:
        args["pon"] = pon
    if interval_list:
        args["L"] = interval_list
    args["ip"] = ip

    check_output(["gatk Mutect2" + build_args_from_dict(args)], shell=True)
    return output

@exists

def learnReadOrientationModel(input_f1r2, output):
    args = defaultdict(list)
    args["I"] = input_f1r2
    args["O"] = output

    check_output(["gatk LearnReadOrientationModel" + build_args_from_dict(args)], shell=True)
    return output

@exists
def GetPileupSummaries(input_bam, output, variant, interval_list, ip="50"):
    args = defaultdict(list)
    args["I"] = input_bam
    args["O"] = output
    args["V"] = variant
    args["L"] = interval_list
    args["ip"] = ip

    check_output(["gatk GetPileupSummaries" + build_args_from_dict(args)], shell=True)
    return output

@exists
def CalculateContamination(tumor_pileupsummaries, normal_pileupsummaries, output, segments):
    args = defaultdict(list)
    args["I"] = tumor_pileupsummaries
    args["matched"] = normal_pileupsummaries
    args["O"] = output
    args["segments"] = segments

    check_output(["gatk CalculateContamination" + build_args_from_dict(args)], shell=True)
    return output

@exists
def FilterMutectCalls(input_variants, output, reference, contamination_table, segments,ob_priors):
    args = defaultdict(list)
    args["V"] = input_variants
    args["O"] = output
    args["R"] = reference
    args["-orientation-bias-artifact-priors"] = ob_priors
    args["-contamination-table"] = contamination_table
    args["-tumor-segmentation"] = segments
    check_output(["gatk FilterMutectCalls" + build_args_from_dict(args)], shell=True)
    return output


@exists
def CollectSequencingArtifactMetrics(input_bam, output, reference, extension=".txt"):
    args = defaultdict(list)
    args["I"] = input_bam
    args["O"] = output
    args["R"] = reference
    args["EXT"] = extension
    check_output(["gatk CollectSequencingArtifactMetrics" + build_args_from_dict(args)], shell=True)
    return output

@exists
def FilterByOrientationBias(input_variants, output, metrics):
    args = defaultdict(list)
    args["V"] = input_variants
    args["O"] = output
    args["-artifact-modes"].append("C/A")
    args["-artifact-modes"].append("T/G")
    args["P"] = metrics

    check_output(["gatk FilterByOrientationBias" + build_args_from_dict(args)], shell=True)
    return output

@exists
def SelectVariants(input_variants, output, reference):
    args = defaultdict(list)
    args["V"] = input_variants
    args["O"] = output
    args["R"] = reference
    args["-exclude-filtered"] = "true"
    args["-exclude-non-variants"] = "true"

    check_output(["gatk SelectVariants" + build_args_from_dict(args)], shell=True)
    return output
