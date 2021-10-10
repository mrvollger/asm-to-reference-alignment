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
    track gene-conversion_{sm}
    parent gene-conversion
    bigDataUrl gene-conversion/{sm}.bb
    shortLabel {sm} gene conversion
    longLabel {sm} gene conversion
    type bigBed 9 +
    itemRgb on
    visibility dense
"""

with open(snakemake.output.track, "w") as out:
    out.write(track_db_header)
    [out.write(track.format(sm=sm)) for sm in snakemake.params.samples]

open(snakemake.output.hub, "w").write(hub)
open(snakemake.output.genomes, "w").write(genomes.format(ref=snakemake.wildcards.ref))
