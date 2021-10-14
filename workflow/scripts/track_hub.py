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
filter.score 5:1000
filterByRange.score on
filterLimits.score 0:1000
filterLabel.score Minimum decrease in mismatches
"""

track = """
    track g-c-{sm}
    parent gene-conversion
    bigDataUrl gene-conversion/{sm}.bb
    shortLabel {sm} gc D/A
    longLabel {sm} gene conversion
    type bigBed 9 +
    itemRgb on
    priority {pri}
    visibility dense
"""


track_db_interact_header = """
track interact-gene-conversion
compositeTrack off
shortLabel interact-gc
longLabel gene conversion interactions
visibility hide
type bigInteract
maxItems 100000
filter.score 5:1000
filterByRange.score on
filterLimits.score 0:1000
"""

track_interact = """
    track interact-g-c-{sm}
    parent gene-conversion
    bigDataUrl gene-conversion/{sm}.interact.bb
    shortLabel {sm} gc
    longLabel {sm} gene conversion interactions
    type bigInteract
    maxHeightPixels 100:30:5
    priority {pri2}
    visibility full
"""

track_super = """
track gene-conversion-by-sample
superTrack on show
shortLabel gc-by-sample
longLabel gene conversion by sample
filter.score 5:1000
filterByRange.score on
filterLimits.score 0:1000

"""

track_comp = """
    track {sm}
    parent gene-conversion-by-sample
    compositeTrack on
    shortLabel {sm}-gc
    longLabel {sm} gene conversion
    type bigWig
    visibility full
        
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
        maxHeightPixels 100:30:5
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
autoScale on
visibility full

track Acceptor 
bigDataUrl all_candidate_windows_acceptor.bw
shortLabel Acceptor 
longLabel Acceptor
type bigWig
color 0,127,211
autoScale on
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


view_format_super = """
# SuperTrack declaration no type or visibility is required
# However "show" is needed to have a superTrack on by  default
track gene-conversion
longLabel gene conversion for HPRC samples
superTrack on show
shortLabel gene-conversion

"""
view_format_comp = """
    # Composite declaration, usually composite tracks are all of one type,
    #   and the type can be declared.
    # When a mixed type (some bigBeds, some bigWigs) you need to use the unusual
    #   'type bed 3' declaration.
    # The subGroup1 line will define groups,
    #   in this case the special 'view' group
    #   (a new subGroup2 could be metadata)
    # Later individual tracks will identify what 'subGroups id' they belong to.
    track gene-conversion-by-sample
    type bed 3
    longLabel gene conversion by sample
    parent gene-conversion
    compositeTrack on
    shortLabel gc-by-sample
    visibility full
    subGroup1 view Views bb=Colored_bigBed_items int=Interact_Data

"""
view_fromat_bb = """
        # This is the unexpected part about views,
        #    you need a separate parent to group the view
        # So this new view-specific stanza with "view id"
        #    can collect all tracks with some visibility settings
        track gene-conversion-by-sample-bb
        parent gene-conversion-by-sample on
        view bb
        visibility dense
        itemRgb on 
        maxItems 100000
        filter.score 5:1000
        filterByRange.score on
        filterLimits.score 0:1000

"""
view_format_bb_sm = """
            # Child bigBed in this view
            # The 'subGroups view=bb' shares this track belongs in a view,
            #    even though a parent declaration is also needed
            # All these tracks should be the same type of data
            track gene-conversion-by-sample-bb-{sm}
            parent gene-conversion-by-sample-bb
            type bigBed 9 +
            longLabel {sm} gene conversion bb 
            bigDataUrl gene-conversion/{sm}.bb
            shortLabel {sm}-gc-bb
            subGroups view=bb

"""
view_format_int = """
        # New View Stanza that collects all interact in this composite
        # This declares related bigInteract tracks    
        track gene-conversion-by-sample-interact
        parent gene-conversion-by-sample on
        view int
        visibility full
        maxHeightPixels 100:30:5
        maxItems 100000
        filter.score 5:1000
        filterByRange.score on
        filterLimits.score 0:1000

"""
view_format_int_sm = """
            # Child one Interact
            track gene-conversion-by-sample-interact-{sm}
            parent gene-conversion-by-sample-interact
            type bigInteract
            longLabel {sm} gene conversion interactions
            bigDataUrl gene-conversion/{sm}.interact.bb
            shortLabel {sm}-gc-interact
            subGroups view=int

"""

with open(snakemake.output.track, "w") as out:
    out.write(all_tracks)
    if False:
        out.write(track_db_header)
        # out.write(track_db_interact_header)
        [
            out.write(
                (track + track_interact).format(
                    sm=sm, pri=2 * idx + 1, pri2=2 * idx + 2
                )
            )  # pri=idx + 1, pri2=idx + 2))
            for idx, sm in enumerate(snakemake.params.samples)
        ]
    elif True:
        out.write(view_format_super)
        out.write(view_format_comp)
        # big beds
        out.write(view_fromat_bb)
        [out.write(view_format_bb_sm.format(sm=sm)) for sm in snakemake.params.samples]
        # bigInteract
        out.write(view_format_int)
        [out.write(view_format_int_sm.format(sm=sm)) for sm in snakemake.params.samples]
    else:
        out.write(track_super)
        [out.write(track_comp.format(sm=sm)) for sm in snakemake.params.samples]

open(snakemake.output.hub, "w").write(hub)

ref = snakemake.wildcards.ref
if ref == "CHM13_V1.1":
    print("changing ref")
    ref = "t2t-chm13-v1.1"
open(snakemake.output.genomes, "w").write(genomes.format(ref=ref))
