#!/usr/bin/env python
import sys
import re
import os
import subprocess
import glob
 
ALIGN_CMD = "/home/simon/git/soladmin/script/soladmin_align.rb -i {0} -o {1} -g {2} -a {3} -r {4}"

def fastq2bam(fqs, bam, genome, aligner="", genome_dir="", force=False):
    bams = [] 
    
    if os.path.exists(bam) and not force:
        sys.stderr.write("{0} already exists, skipping...\n".format(bam))
        return

    for fq in fqs:
        bname = fq.replace(".fq.gz", "")
        
        cmd = ALIGN_CMD.format(fq, bname, genome, aligner, genome_dir)

        sys.stderr.write("Mapping {0} to {1}\n".format(fq, genome))
        p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        stdout, stderr = p.communicate()
        
        if not os.path.exists("{0}.bam".format(bname)):
            sys.stderr.write("Alignment failed\n")
            sys.stderr.write(stderr)
            sys.exit(1)
        bams.append("{0}.bam".format(bname))
    
    if len(bams) == 1:
        os.rename(bams[0], bam)
        os.rename(bams[0] + ".bai", bam + ".bai")

    elif len(bams) > 1:
        cmd = "{0} merge -f {1} {2} && samtools index {1}"
        cmd = cmd.format(
                         "samtools",
                         bam,
                         " ".join(bams),
                         )
        ret = subprocess.call(cmd, shell=True)
        if not ret:
            for bamfile in bams:
                os.unlink(bamfile)
                os.unlink("{0}.bai".format(bamfile))
        else:
            sys.stderr.write("Merging failed!\n")
            sys.exit(1)

def sra2fastq(sra, name, outdir=".", keep_sra=False):
    try:
        FASTQ_DUMP = "fastq-dump"
        cmd = "{0} --split-files --gzip {1} -O {2}".format(
                                                  FASTQ_DUMP,
                                                  sra,
                                                  outdir,
                                                  )
    
        
        p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        stdout, stderr = p.communicate()
        if stderr:
            raise Exception(stderr)
        
        sys.stderr.write("Successfully converted {0} to fastq\n".format(sra))
        
        if not keep_sra:
            os.unlink(sra)
            
        base = os.path.splitext(os.path.basename(sra))[0]
        #os.unlink(sra)
        p = re.compile(r'(SRR.+)\.sra')
        m = p.search(sra)
        srr = m.group(1)
  
        fqs = []
        for old_fq in glob.glob(os.path.join(outdir, "*{0}*fastq.gz".format(base))):
            fastq = os.path.join(outdir, "{0}.{1}.fq.gz".format(srr, name))
            os.rename(old_fq, fastq)
            fqs.append(fastq)
    
        return fqs
    
    except Exception as e:
        sys.stderr.write("fastq-dump of {0} failed :(\n".format(sra))
        sys.stderr.write("{0}\n".format(e))
        return []

def bam2bw(bam, bw, library=None, force=False):
    """ Convert bam to UCSC bigWig 
    """
    # Don't overwrite existing bigWig if not explicitly stated
    if os.path.exists(bw) and not force:
        sys.stderr.write("{0} already exists, skipping...\n".format(bw))
        return None, None
    
    if library == "RNA-Seq":
        # Do not extend, keep duplicates
        cmd = "bam2bw -i {0} -o {1} -D".format(bam, bw)
    elif library == "ChIP-Seq":
        # Extend to estimated fragment size and normalize
        cmd = "bam2bw -i {0} -o {1} -e auto -c".format(bam, bw)
    else:
        # Same as ChIP-seq for now
        cmd = "bam2bw -i {0} -o {1} -e auto -c".format(bam, bw)
    
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout,stderr = p.communicate()
    return stdout, stderr
