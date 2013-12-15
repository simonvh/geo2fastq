#!/usr/bin/env python
import sys
import re
import os
from trackhub import Hub, GenomesFile, Genome, TrackDb, Track
from trackhub.upload import upload_track, upload_hub
from trackhub.helpers import show_rendered_files

HUB_URLBASE = 'http://mbpcsimon.azn.nl/trackhubs'

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
