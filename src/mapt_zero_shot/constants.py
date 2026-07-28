"""Canonical MAPT/Tau sequence and annotation constants."""

AMINO_ACIDS = tuple("ACDEFGHIKLMNPQRSTVWY")

# Adult CNS human Tau sequence, 441 aa, commonly used for 2N4R Tau/hTau40.
# UniProt P10636 canonical is PNS-tau/Big-tau (758 aa); this project uses the
# P10636-8 Tau-F isoform coordinate system for v1 of the atlas.
MAPT_441_SEQUENCE = (
    "MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKESPLQTPTEDGSEEPGSETSDAK"
    "STPTAEDVTAPLVDEGAPGKQAAAQPHTEIPEGTTAEEAGIGDTPSLEDEAAGHVTQARMVSKSKDGT"
    "GSDDKKAKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPKTPPSSGEPPKSGDRSGYSSPGSP"
    "GTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAKSRLQTAPVPMPDLKNVKSKIGSTENLKHQPG"
    "GGKVQIINKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEV"
    "KSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVYKSPVVSGDTSPRH"
    "LSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL"
)

REFERENCE = {
    "gene": "MAPT",
    "protein_name": "Tau",
    "isoform": "2N4R/Tau-F/hTau40",
    "uniprot": "P10636-8",
    "length": 441,
}

REGIONS = (
    ("N_terminal_projection", 1, 150),
    ("proline_rich_region", 151, 243),
    ("microtubule_repeat_R1", 244, 274),
    ("microtubule_repeat_R2_exon10", 275, 305),
    ("microtubule_repeat_R3", 306, 336),
    ("microtubule_repeat_R4", 337, 368),
    ("C_terminal_tail", 369, 441),
)

MOTIFS = (
    ("PHF6_star_VQIINK", 275, 280),
    ("PHF6_VQIVYK", 306, 311),
)

# Commonly discussed Tau phosphorylation sites; used for proximity annotation
# rather than as a complete PTM atlas.
PTM_SITES = (181, 202, 205, 212, 214, 217, 231, 235, 262, 356, 396, 404, 422)

KNOWN_PATHOGENIC_HOTSPOTS = (257, 279, 280, 301, 305, 315, 317, 320, 337, 389, 406)
