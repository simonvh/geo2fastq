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
from collections import defaultdict
import pp
import pickle

from geo2fastq.convert import sra2fastq


GEOFTP_URLBASE = "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/by_series/{0}/{0}_family.soft.gz"
FTP_ROOT = "ftp-trace.ncbi.nlm.nih.gov"

tax2genome = defaultdict(str)
tax2genome["8364"] = "xenTro3"
tax2genome["7955"] = "danRer7"
tax2genome["9606"] = "hg19"
tax2genome["10090"] = "mm10"
tax2genome["10116"] = "rn4"
tax2genome["8355"] = "XENLA_JGIv7b"

class Geo:
    def __init__(self, gse="", email=""):
        """ Create Geo object to search GEO and parse GEO experiment data
        
        :param gse: GEO accession or filehandle to GEO soft file
        :type accession: string or file
        """
        
        self.gse = gse
        self.email = email
        self.samples = {}
        if gse:
            try:
                callable(gse.read)
                gen = self.get_sample_info(fo=gse)
            except:
                self.gse = gse
                gen = self.get_sample_info(accession=self.gse) 
            for sample in gen:
                self.samples[sample['gsm']] = sample

    @classmethod
    def search(self, term):
        """ Search NCBI GEO and return a dictionary of GEO accessions
        describing series (GSE) and samples (GSM).
        
        :param term: search term(s)
        :type term: string
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
    
    def get_sample_info(self, accession=None, fo=None):
        """ Accepts a single accession or a list of accessions and returns
        an generator of detailed sample information.
       
        Supply either an accession or a SOFT filehandle.    
     
        :param accession: GEO accession
        :type accession: string
        :params fo: filehandle to gzipped SOFT file
        :type fo: file
        """
        
        if not fo:
            if not accession:
                sys.stderr.write("Please supply either accession or fo")
                sys.exit(1)
            # Dowload GEO SOFT file
            fh = urlopen(GEOFTP_URLBASE.format(accession))
            fo = StringIO.StringIO(fh.read())
        
        # Parse gzipped SOFT file
        g = gzip.GzipFile(fileobj=fo)
        record = self._soft_read(g)
        
        self.gse = record['SERIES'].keys()[0]
        
        for gsm, data in record['SAMPLE'].items():
            #for k,v in data.items():
            #    print k,v
            sample = {'gsm':gsm}
            sample['tax_id'] = data['Sample_taxid_ch1'][0]
            sample['genome'] = tax2genome[sample['tax_id']]
            sample['sra'] = []
            sample['name'] = data['Sample_title'][0]
            sample['library'] = data['Sample_library_strategy'][0]
            sample['info'] = data['Sample_characteristics_ch1']
            for sra_link in [x for x in data['Sample_relation'] if x.startswith("SRA")]:
                sample['sra'].append(sra_link)
            yield sample
 
    def download_srx(self, srx, outdir):
        
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
   
    def download_sra(self, sra_link, outdir="./"):
        p = re.compile(r'term=(\w+)')
        m = p.search(sra_link)
        if m:
            srx = m.group(1)
            for fname in self.download_srx(srx, outdir):
                yield fname
                    
        else:
            sys.stderr.write("No SRA link found for {0}\n".format(gsm))

           
    def download(self, gsm="", outdir="./", format="fastq"):
        job_server = pp.Server(secret="polentafries") 
        jobs = []
 
        outdir = os.path.join(outdir, self.gse)
        samples = self.samples.values()
        if gsm:
            if not self.samples.has_key[gsm]:
                raise Exeception
            samples = [self.samples[gsm]]
        
        for sample in samples:
            for fname in self._download_sample(sample, outdir=outdir):
                job = job_server.submit(
                                        sra2fastq,
                                        (fname, sample['gsm'], self.gse),
                                        (),
                                        ("os", "sys", "re", "glob", "subprocess")
                                        )
                jobs.append(job)
        
        for job in jobs:

            job()
        
    def _download_sample(self, sample, outdir="."):
        for sra_link in sample['sra']:
            for fname in self.download_sra(sra_link, outdir):
                yield fname

    def _soft_read(self, fh):
        """ Parses a filehandle of a SOFT formatted file
        """
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

if __name__ == "__main__":
    #for k,v in Geo.search("Heeringen AND Veenstra").items():
    #    print k,v
    
    #x = Geo("GSE14025")
    x = Geo(open("tests/data/GSE14025_family.soft.gz"))
    for sample in x.samples.values():
        print sample
    

