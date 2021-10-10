track_db_header = """
track gene-conversion
compositeTrack off
shortLabel gene-conversion
longLabel gene-conversion
visibility hide
priority 30
type bigBed 9 +
itemRgb on
maxItems 100000
"""

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

track = """
    track gc {sm}
    parent gene-conversion
    bigDataUrl gene-conversion/{sm}.bb
    shortLabel {sm} gc
    longLabel {sm} gene conversion
    type bigBed 9 +
    itemRgb on
    visibility dense
"""

all_tracks = """
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
shortLabel gene conversion
longLabel gene conversion
type bigBed 9 +
itemRgb on
visibility dense
maxItems 100000

"""


with open(snakemake.output.track, "w") as out:
    out.write(all_tracks)
    out.write(track_db_header)
    [out.write(track.format(sm=sm)) for sm in snakemake.params.samples]

open(snakemake.output.hub, "w").write(hub)

ref = snakemake.wildcards.ref
if ref == "CHM13_V1.1":
    print("changing ref")
    ref = "t2t-chm13-v1.1"
open(snakemake.output.genomes, "w").write(genomes.format(ref=ref))
