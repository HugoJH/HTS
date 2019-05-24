
#!/usr/bin/env python
from subprocess import check_output
from collections import defaultdict
from multiprocessing import cpu_count
from utils import exists, build_args_from_dict

@exists
def bwaMem(input_fastq, output, reference, threads=str(cpu_count())):
    args = defaultdict(list)
    args["t"] = threads
    check_output(["bwa mem -M -p " + build_args_from_dict(args) + " " + reference + " " + input_fastq +  " > " + output], shell=True)
