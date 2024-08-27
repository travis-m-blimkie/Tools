import os
import subprocess as sp
import pandas as pd
from glob import glob
from re import search
from os.path import exists


# Inputs to be filled out by the user before running the script. Only R1 will be used, excluding R2 from STAR and
# everything downstream. Make sure to leave the trailing "/" for directories. The 'threads' parameter needs to be
# string to work with concatenation.
my_fastq_dir = "path/to/fastq/"
my_genome_dir = "path/to/genome/"
my_gtf_file = my_genome_dir + "file.gtf"
my_fasta_file = my_genome_dir + "file.fa"
threads = "8"

# Should Samtools/CRAM be run (True/False)?
do_cram = False

# To run the script, modify the above parameters as desired then use the below command to run it in the background,
# and save the output to a log file.
# python pipeline.py &> pipeline.log &

# Here's where the script starts. No editing should be required past here! Unless you want to modify the parameters
# of, for example, STAR


# FastQC
def run_fastqc(input_dir):
    os.makedirs(os.path.dirname("FastQC/"), exist_ok=True)
    fastqc_command = "fastqc -q -t " + threads + " -o FastQC/ " + input_dir + "*.fastq.gz"
    sp.run(fastqc_command, shell=True)


# STAR
def run_star(input_dir, genome_dir):
    os.makedirs(os.path.dirname("STAR/"), exist_ok=True)
    fastq_files = glob(input_dir + "*R1.fastq.gz")

    for f in fastq_files:
        prefix = os.path.basename(f.replace("_R1.fastq.gz", ""))
        check_name = "STAR/" + prefix + "_Aligned.sortedByCoord.out.bam"

        star_command = "STAR --runThreadN " + threads + " --runMode alignReads --genomeLoad LoadAndKeep --limitBAMsortRAM " \
            "48000000000 --outSAMtype BAM SortedByCoordinate --genomeDir " + genome_dir + " --readFilesIn " \
            + f + " --readFilesCommand zcat --outFileNamePrefix STAR/" + prefix + "_"

        if exists(check_name):
            print("Skipping, STAR output already exists for", prefix)
        else:
            sp.run(star_command, shell=True)

    sp.run("STAR --genomeDir " + genome_dir + " --genomeLoad Remove", shell=True)
    sp.run("rmdir _STARtmp/", shell=True)
    sp.run("rm Aligned.out.sam Log.out Log.final.out Log.progress.out SJ.out.tab", shell=True)
    sp.run("rm STAR/*SJ.out.tab", shell=True)


# HTSeq
def run_htseq(gtf_file):
    os.makedirs(os.path.dirname("HTSeq/"), exist_ok=True)

    htseq_parallel = "find STAR/ -name '*.bam' | parallel --jobs " + threads + " 'htseq-count -s reverse -a 10 -f bam -r pos {} " \
                     + gtf_file + " > HTSeq/{/.}.count'"
    sp.run(htseq_parallel, shell=True)

    # Rename the counts files (remove "_Aligned.sortedByCoord.out")
    count_files = glob("HTSeq/*.count")
    for c in count_files:
        os.rename(c, c.replace("_Aligned.sortedByCoord.out", ""))


# MultiQC
def run_multiqc():
    multiqc_cmd = "multiqc FastQC/ STAR/ HTSeq/"
    sp.run(multiqc_cmd, shell=True)


# Samtools/CRAM
def run_samtools(fasta_file):
    bam_files = glob("STAR/*.bam")
    for b in bam_files:
        cram_name = b.replace("bam", "cram")
        cram_command = "samtools view --threads " + threads + " -C " + b + " -T " + fasta_file + " > " + cram_name

        if exists(cram_name):
            print("Skipping, CRAM file already exists for", cram_name)
        else:
            sp.run(cram_command, shell=True)


# Version information
def run_versions(run_cram):
    version_fastqc = sp.check_output("fastqc --version", shell=True).decode("utf-8").replace("\n", "").replace("FastQC v", "")

    version_star = sp.check_output("STAR --version", shell=True).decode("utf-8").replace("\n", "")

    version_htseq_raw = sp.check_output("htseq-count --help | grep version", shell=True, stderr=sp.STDOUT).decode("utf-8")
    version_htseq = search("[0-9]\\.[0-9]{1,2}\\.[0-9]{1,2}", version_htseq_raw).group(0)

    version_multiqc_raw = sp.check_output("multiqc --version", shell=True).decode("utf-8").replace("\n", "")
    version_multiqc = search("[0-9]{1,2}\\.[0-9]{1,2}", version_multiqc_raw).group(0)

    if not(run_cram):
        program_version_dict = {
            "Program": ["FastQC", "STAR", "HTSeq", "MultiQC"],
            "Version": [version_fastqc, version_star, version_htseq, version_multiqc]
        }
    elif run_cram :
        version_samtools_raw = sp.check_output("samtools --version | grep samtools", shell=True).decode("utf-8").replace("\n", "")
        version_samtools = search("[0-9]{1,2}\\.[0-9]{1,2}", version_samtools_raw).group(0)

        program_version_dict = {
            "Program": ["FastQC", "STAR", "HTSeq", "MultiQC", "Samtools"],
            "Version": [version_fastqc, version_star, version_htseq, version_multiqc, version_samtools]
        }
    
    program_version_df = pd.DataFrame(program_version_dict)
    program_version_df.to_csv("pipeline_program_versions.csv", index=False)


# Run all the functions
#run_fastqc(input_dir=my_fastq_dir)
#run_star(input_dir=my_fastq_dir, genome_dir=my_genome_dir)
#run_htseq(gtf_file=my_gtf_file)
#run_multiqc()
#if do_cram:
#    run_samtools(fasta_file=my_fasta_file)
run_versions(run_cram=do_cram)

