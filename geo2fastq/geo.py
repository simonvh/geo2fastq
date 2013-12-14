#!/usr/bin/env python
from urllib2 import urlopen
from Bio import Entrez
import sys
import StringIO
import gzip
from ftplib import FTP
import re
import os
import subprocess as sp
import argparse
import glob
import pp
from collections import defaultdict
from trackhub import Hub, GenomesFile, Genome, TrackDb, Track
from trackhub.upload import upload_track, upload_hub
from trackhub.helpers import show_rendered_files

HUB_URLBASE = 'http://mbpcsimon.azn.nl/trackhubs'
GEOFTP_URLBASE = "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/by_series/{0}/{0}_family.soft.gz"
FTP_ROOT = "ftp-trace.ncbi.nlm.nih.gov"
VERSION = '1.0'
ALIGN_CMD = "/home/simon/git/soladmin/script/soladmin_align.rb -i {0} -o {1} -a bwa -g {2} -r /usr/share/genomes"

tax2genome = defaultdict(str)
tax2genome["8364"] = "xenTro3"
tax2genome["7955"] = "danRer7"

class Geo:
    def __init__(self):
		pass

    @classmethod
    def search(self, term):
        """ Search NCBI GEO and return a dictionary of GEO accessions
		describing series (GSE) and samples (GSM).
		"""
		Entrez.email = "s.vanheeringen@ncmls.ru.nl"      
        handle = Entrez.esearch("gds", term)
        record = Entrez.read(handle)
        results = {}
        for gid in record['IdList']:
            handle = Entrez.esummary(db="gds", id=gid)
            r = Entrez.read(handle)[0]
            if r['entryType'] == 'GSE':
                gse = "GSE{0}".format(r['GSE'])
                results[gse] = {
                                'title': r['title'].encode('utf8'),
                                'samples':{}
                                }
                for sample in r['Samples']:
                    results[gse]['samples'][sample['Accession']] = sample['Title'].encode('utf8')
        return results
    
        
    def get_sample_info(accessions):
		""" Accepts a single accession or a list of accessions and returns
		detailed sample information.
		"""
        if type("") == type(accessions):
            accessions = [accessions]
        
        for accession in accessions:
            #print URLBASE.format(accession)
            fh = urlopen(GEOFTP_URLBASE.format(accession))
            fo = StringIO.StringIO(fh.read())
            g = gzip.GzipFile(fileobj=fo)
            record = soft_read(g)
            
            for gsm, data in record['SAMPLE'].items():
                #for k,v in data.items():
                #    print k,v
                sample = {'gsm':gsm}
                sample['tax_id'] = data['Sample_taxid_ch1'][0]
                sample['sra'] = []
                sample['name'] = data['Sample_title'][0]
                sample['library'] = data['Sample_library_strategy'][0]
                sample['info'] = data['Sample_characteristics_ch1']
                for sra_link in [x for x in data['Sample_relation'] if x.startswith("SRA")]:
                    sample['sra'].append(sra_link)
                yield sample
            
    @classmethod
    def download(self, accession):
		pass

def soft_read(fh):
    soft = {}
    current = soft
    for line in fh:
        ltype = line[0]
        try:
            key, val = line[1:].strip().split(" = ")
            if line.startswith("^"):
                soft.setdefault(key, {})[val] = {}
                current =  soft[key][val]
            else:
                current.setdefault(key, []).append(val)
        except:
            pass
    return soft

def retrieve_sra(srx, outdir):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    ftp = FTP(FTP_ROOT)
    ftp.login()
    rootdir = "/sra/sra-instant/reads/ByExp/sra/SRX/{0}/{1}".format(srx[:6], srx)
    ftp.cwd(rootdir)
    dirs  = []
    ftp.retrlines('LIST', callback=lambda x: dirs.append(x.split(" ")[-1]))
    for dirname in dirs:
        ftp.cwd(os.path.join(rootdir, dirname))
        fnames = []
        ftp.retrlines('LIST', callback=lambda x: fnames.append(x.split(" ")[-1]))
        for fname in fnames:
            local_name = os.path.join(outdir,fname)
            sys.stderr.write("Downloading {0}...\n".format(local_name))
            f = open(local_name, "w")
            ftp.retrbinary(
                           "RETR {0}".format(os.path.join(rootdir, dirname, fname)),
                           f.write
                           )
            f.close()
            yield local_name


def download_sra(sra_link, outdir):
    p = re.compile(r'term=(\w+)')
    m = p.search(sra_link)
    if m:
        srx = m.group(1)
        for fname in retrieve_sra(srx, outdir):
            yield fname
                    
    else:
        sys.stderr.write("No SRA link found for {0}\n".format(gsm))

def sra2fastq(sra, name, outdir="."):
    try:
        FASTQ_DUMP = "fastq-dump"
        cmd = "{0} --split-files -A {1} -O {2}".format(
                                                  FASTQ_DUMP,
                                                  sra,
                                                  outdir,
                                                  )
    
        fqs = []
        p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        stdout, stderr = p.communicate()
        if stderr:
            raise Exception(stderr)
        
        sys.stderr.write("Successfully converted {0} to fastq\n".format(sra))
        base = os.path.splitext(os.path.basename(sra))[0]
        #os.unlink(sra)
        p = re.compile(r'(SRR.+)\.sra')
        m = p.search(sra)
        srr = m.group(1)
  
        for old_fq in glob.glob(os.path.join(outdir, "*{0}*fastq".format(base))):
            fastq = os.path.join(outdir, "{0}.{1}.fq".format(srr, name))
            os.rename(old_fq, fastq)
            cmd = "pigz {0}".format(fastq)
            p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stdout, stderr = p.communicate()
            fqs.append("{0}.gz".format(fastq))
    
        return fqs
    
    except Exception as e:
        sys.stderr.write("fastq-dump of {0} failed :(\n".format(sra))
        sys.stderr.write("{0}\n".format(e))
        return []



def bam2bw(bam, bw):
    cmd = "bam2bw -i {0} -o {1} -e auto -c".format(bam, bw)
    print cmd
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout,stderr = p.communicate()
    return stdout, stderr

def create_hub(gse, samples, upload_dir, user, host,  email):
    hub = Hub(
        hub=gse,
        short_label=gse,
        long_label="Hub for {0}".format(gse),
        email=email)

    genomes_file = GenomesFile()
    
    trackdb = TrackDb()

    local_dir = gse

    hub.remote_fn = os.path.join(upload_dir, gse, os.path.basename(hub.local_fn))
    
    all_tracks = {}
    
    for sample in samples:
        genome = sample['genome']
        all_tracks.setdefault(genome, [])

        name = re.sub('[^0-9a-zA-Z]+', '_',sample['name'])
        track = Track(
            name=name,
            url=os.path.join(HUB_URLBASE, gse, genome, "{0}.bw".format(sample['gsm'])),
            tracktype='bigWig',
            short_label=sample['gsm'],
            long_label=name,
            color='128,128,0',
            maxHeightPixels='30:30:11',
            )
        basename = os.path.basename(track.url)
        track.local_fn = os.path.join(local_dir, basename)
        track.remote_fn = os.path.join(upload_dir, gse, genome, basename)
        all_tracks[genome].append(track)
    
    for build,tracks in all_tracks.items(): 

        genome = Genome(build)
        trackdb.add_tracks(tracks)
        genome.add_trackdb(trackdb)
        genomes_file.add_genome(genome)
        hub.add_genomes_file(genomes_file)

    results = hub.render()

    for track in trackdb.tracks:
        upload_track(track=track, host=host, user=user)
    
    upload_hub(hub=hub, host=host, user=user)


if __name__ == "__main__":
    print Geo.search("Heeringen AND Veenstra")

