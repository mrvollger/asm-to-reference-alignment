table gene_conversion
"Windows with evidence of gene conversion"
(
string  chrom;		"Chromosome for original alignment"
uint    chromStart;	"Start position of original alignment"
uint    chromEnd;	"End position of original alignment"
string  name;		"Name with mismatch delta"
float    perID_by_all;		"Percent identity of alignment at original alignment location"
string  strand;		"strand NA"
uint    mismatches;		"mismatches at original alignment"
uint    donorMismatches;		"mismatches at donor alignment"
string  color;		"RGB color, blue = acceptor, orange = donor"
string  donorChrom;		"Chromosome for donor alignment"
uint    donorStart;	"Start position of donor alignment"
uint    donorEnd;	"End position of donor alignment"
float    donor_perID_by_all;		"Percent identity of alignment at donor location"
)