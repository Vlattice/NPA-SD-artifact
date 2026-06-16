"""Data sources and analysis parameters.

All defaults reproduce the sliding-window non-parental-allele (NPA) scan used in
the original report, so that the only intended change between the original and
corrected analyses is the segmental-duplication filter.
"""

# --- input data -------------------------------------------------------------

SD_TRACK_URL = (
    "https://raw.githubusercontent.com/abdeldar/SD_network/"
    "master/inputs/GRCh38GenomicSuperDup_sort.tab"
)

PEDIGREE_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt"
)

VCF_BASE_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/"
)


def vcf_filename(chrom: str) -> str:
    """Phased panel filename for one chromosome in the high-coverage release."""
    return (
        f"1kGP_high_coverage_Illumina.{chrom}.filtered."
        f"SNV_INDEL_SV_phased_panel.vcf.gz"
    )


# --- NPA detection (match the original scan) --------------------------------

WINDOW_SIZE = 60       # SNVs per sliding window
STEP = 20              # window advance, in SNVs
NPA_THRESHOLD = 5      # a window is positive for a trio at >= this many NPAs

# --- segmental-duplication handling -----------------------------------------

SD_BIN_SIZE = 500_000   # bin size for the genome-wide SD-vs-control comparison
SD_RICH_MIN = 5         # >= this many SD intervals in a bin -> "SD-rich"
SD_WINDOW_TOLERANCE = 2  # a window overlapping > this many SDs is dropped

# --- GRCh38 autosome lengths (bp) -------------------------------------------

CHROM_LENGTHS = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468,
}

# Loci used in the chromosome-3 reproduction.
CHR3_SD_HOTSPOT = ("chr3", 75_200_000, 75_700_000)   # strongest reported signal
CHR3_FRA3B = ("chr3", 59_700_000, 62_500_000)        # SD-free control locus

# The one apparent exception found genome-wide (subtelomeric 16q24.3).
CANDIDATE_LOCUS = ("chr16", 88_500_000, 89_000_000)
