hub = """
hub gene-conversion
shortLabel gene-conversion
longLabel gene-conversion
genomesFile genomes.txt
email mvollger.edu
"""

genomes = """
genome {ref}
trackDb trackDb.txt
"""

track_db_header = """
track gene-conversion
compositeTrack off
shortLabel gene-conversion
longLabel gene-conversion
visibility hide
type bigBed 9 +
itemRgb on
maxItems 100000
"""

track = """
    track g-c-{sm}
    parent gene-conversion
    bigDataUrl gene-conversion/{sm}.bb
    shortLabel {sm} gc
    longLabel {sm} gene conversion
    type bigBed 9 +
    itemRgb on
    visibility dense
"""


track_db_interact_header = """
track interact-gene-conversion
compositeTrack off
shortLabel interact-gc
longLabel gene conversion interactions
visibility hide
priority 30
type bigInteract
maxItems 100000
"""

track_interact = """
    track interact-g-c-{sm}
    parent interact-gene-conversion
    bigDataUrl gene-conversion/{sm}.interact.bb
    shortLabel {sm} interact-gc
    longLabel {sm} gene conversion interactions
    type bigInteract
    maxHeightPixels 100:20:5
    visibility full
"""

track_super = """
track gene-conversion-by-sample
superTrack on show
shortLabel gc-by-sample
longLabel gene conversion by sample

"""

track_comp = """
    track {sm}
    parent gene-conversion-by-sample
    compositeTrack on
    shortLabel {sm}-gc
    longLabel {sm} gene conversion
    type bigWig
    visibility dense
        
        track gc-{sm}
        parent {sm}
        bigDataUrl gene-conversion/{sm}.bb
        shortLabel {sm} gc
        longLabel {sm} gene conversion
        type bigBed 9 +
        itemRgb on
        visibility dense

        track interact-{sm}
        parent {sm}
        bigDataUrl gene-conversion/{sm}.interact.bb
        shortLabel {sm} interact
        longLabel {sm} interactions
        type bigInteract
        maxHeightPixels 100:20:5
        visibility full
   
"""


all_tracks = """
track g-c-interact
bigDataUrl all_candidate_interactions.bb
shortLabel all gc interact
longLabel all gene conversion interactions
type bigInteract
visibility hide

track Donor 
bigDataUrl all_candidate_windows_donor.bw
shortLabel Donor 
longLabel Donor
type bigWig
color 211,144,0
visibility full

track Acceptor 
bigDataUrl all_candidate_windows_acceptor.bw
shortLabel Acceptor 
longLabel Acceptor
type bigWig
color 0,127,211
visibility full

track gene-conversion-windows
bigDataUrl all_candidate_windows.bb
shortLabel all g-c windows
longLabel all gene conversion windows
type bigBed 9 +
itemRgb on
visibility dense
maxItems 100000

"""


with open(snakemake.output.track, "w") as out:
    out.write(all_tracks)
    if False:
        out.write(track_db_header)
        out.write(track_db_interact_header)
        [
            out.write((track + track_interact).format(sm=sm))
            for idx, sm in enumerate(snakemake.params.samples)
        ]
    else:
        out.write(track_super)
        [out.write(track_comp.format(sm=sm)) for sm in snakemake.params.samples]

open(snakemake.output.hub, "w").write(hub)

ref = snakemake.wildcards.ref
if ref == "CHM13_V1.1":
    print("changing ref")
    ref = "t2t-chm13-v1.1"
open(snakemake.output.genomes, "w").write(genomes.format(ref=ref))
