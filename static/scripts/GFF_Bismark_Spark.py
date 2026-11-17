"""
======================================================================================
GFF_Bismark_annotation — Apache Spark Version
======================================================================================
CHANGES FROM ORIGINAL:
1. SparkSession added with memory configs
2. Bismark reading: skiprows=1 replaced by filtering "Bismark" header line
3. add_all_methylation_counts: iterrows() loop replaced by Spark range join
4. CSV writing: save_as_flat_csv() produces identical flat files as pandas
5. Everything else (GFF3 parsing, region construction) UNCHANGED
======================================================================================
"""

from pyspark.sql import SparkSession
from pyspark.sql import functions as F
from pyspark.sql.types import StructType, StructField, StringType, LongType
import pandas as pd
import re
import os
import glob
import shutil

# ──────────────────────────────────────────────────────────────────────────────
# CHANGE 1: SparkSession with memory configuration
# local[*] auto-detects all cores given by PBS (ncpus=20)
# Memory must be set explicitly — default is only 1GB which crashes on 8GB files
# ──────────────────────────────────────────────────────────────────────────────
spark = SparkSession.builder \
    .appName("GFF_Bismark_Annotation") \
    .master("local[*]") \
    .config("spark.driver.memory", "250g") \
    .config("spark.executor.memory", "250g") \
    .config("spark.driver.maxResultSize", "20g") \
    .config("spark.sql.shuffle.partitions", "60") \
    .config("spark.memory.offHeap.enabled", "true") \
    .config("spark.memory.offHeap.size", "20g") \
    .getOrCreate()
spark.sparkContext.setLogLevel("WARN")


# ──────────────────────────────────────────────────────────────────────────────
# CHANGE 2: Flat CSV writer helper
# Produces single flat .csv file identical to pandas .to_csv()
# coalesce(1) → one part file → rename to final path → delete temp folder
# ──────────────────────────────────────────────────────────────────────────────
def save_as_flat_csv(sdf, final_path):
    tmp_path = final_path + "_tmp_spark"
    sdf.coalesce(1).write.csv(tmp_path, header=True, mode="overwrite")
    part_files = glob.glob(os.path.join(tmp_path, "part-*.csv"))
    if not part_files:
        raise FileNotFoundError(f"No part file found in {tmp_path}")
    shutil.move(part_files[0], final_path)
    shutil.rmtree(tmp_path)
    print(f"  Written: {final_path}")


def GFF_Bismark_annotation(GFF3, PROMOTER_LEN, output_dir,
                            CpG_OT, CHG_OT, CHH_OT,
                            CpG_OB, CHG_OB, CHH_OB):

    os.makedirs(output_dir, exist_ok=True)
   

    # ──────────────────────────────────────────────────────────────────────────
    # NO CHANGE: GFF3 parsing with pandas — identical to original
    # regex logic, protein sequence extraction all unchanged
    # ──────────────────────────────────────────────────────────────────────────
    df = pd.read_csv(
        GFF3,
        sep="\t",
        comment="#",
        header=None
    )
    df.columns = ['ID', 'Source', 'Portion', 'Start', 'End',
                  'score', 'Strand', 'Phase', 'Reference_id']

    def robust_gene_id(ref, portion):
        if pd.isna(ref):
            return ref
        # GTF style: gene_id "g1"
        m = re.search(r'gene_id "([^"]+)"', ref)
        if m:
            return m.group(1)
        # GFF3 gene line: ID=geneX
        if portion == "gene":
            m = re.search(r'ID=([^;]+)', ref)
            if m:
                return m.group(1)
        # GFF3 mRNA/transcript: Parent=geneX
        if portion in ["mRNA", "transcript"]:
            m = re.search(r'Parent=([^;]+)', ref)
            if m:
                return m.group(1)
        # Child features: Parent=mRNA → gene
        m = re.search(r'Parent=([^;]+)', ref)
        if m:
            return m.group(1).split(",")[0]
        # Transcript-like names: g1.t1 → g1
        if "." in ref:
            return ref.split(".")[0]
        return ref

    df["Reference_id"] = df.apply(
        lambda r: robust_gene_id(r["Reference_id"], r["Portion"]),
        axis=1
    )

    # NO CHANGE: Protein sequence extraction from GFF3 comment lines
    protein_dict = {}
    with open(GFF3) as f:
        current_gene = None
        collecting = False
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith("# start gene"):
                current_gene = line.split()[-1]
                seq = []
                collecting = False
            elif line.startswith("# protein sequence"):
                collecting = True
                seq.append(line.split("[", 1)[1])
            elif collecting:
                if "]" in line:
                    seq.append(line.replace("]", ""))
                    protein_dict[current_gene] = "".join(seq).replace(" ", "")
                    collecting = False
                else:
                    seq.append(line)

    if not protein_dict:
        print("WARNING: No protein sequences found in GFF comments")

    # ──────────────────────────────────────────────────────────────────────────
    # NO CHANGE: Strand splitting — identical to original
    # ──────────────────────────────────────────────────────────────────────────
    Positive_strand = df[df['Strand'] == '+'].copy()
    Negative_strand = df[df['Strand'] == '-'].copy()

    # ──────────────────────────────────────────────────────────────────────────
    # NO CHANGE: All genomic region construction — identical to original
    # Positive strand: transcript, start_codon, stop_codon, UTRs, CDS, intron,
    #                  exon, promoter
    # Negative strand: same regions
    # ──────────────────────────────────────────────────────────────────────────

    # ── Positive Strand ──
    Positive_strand_transcript  = Positive_strand[Positive_strand['Portion'] == 'transcript'].copy()
    Positive_strand_start_codon = Positive_strand[Positive_strand['Portion'] == 'start_codon'].copy()
    Positive_strand_stop_codon  = Positive_strand[Positive_strand['Portion'] == 'stop_codon'].copy()

    merged = Positive_strand_transcript.merge(
        Positive_strand_start_codon, on='Reference_id',
        suffixes=('_tx', '_sc'), how='inner'
    )
    five_prime_UTR_positive = pd.DataFrame({
        'ID':           merged['ID_tx'],
        'Source':       merged['Source_tx'],
        'Portion':      'five_prime_UTR',
        'Start':        merged['Start_tx'],
        'End':          merged['Start_sc'] - 1,
        'score':        merged['score_tx'],
        'Strand':       merged['Strand_tx'],
        'Phase':        '.',
        'Reference_id': merged['Reference_id']
    })

    merged_for_three = Positive_strand_transcript.merge(
        Positive_strand_stop_codon, on='Reference_id',
        suffixes=('_tx', '_sc'), how='inner'
    )
    Three_prime_UTR_positive = pd.DataFrame({
        'ID':           merged_for_three['ID_tx'],
        'Source':       merged_for_three['Source_tx'],
        'Portion':      'three_prime_UTR',
        'Start':        merged_for_three['End_sc'] + 1,
        'End':          merged_for_three['End_tx'],
        'score':        merged_for_three['score_tx'],
        'Strand':       merged_for_three['Strand_tx'],
        'Phase':        '.',
        'Reference_id': merged_for_three['Reference_id']
    })

    CDS_Positive_strand    = Positive_strand[Positive_strand['Portion'] == 'CDS'].copy()
    intron_positive_strand = Positive_strand[Positive_strand['Portion'] == 'intron'].copy()
    exon_positive_strand   = Positive_strand[Positive_strand['Portion'] == 'exon'].copy()

    CDS_Positive_strand['length']         = CDS_Positive_strand['End']     - CDS_Positive_strand['Start']     + 1
    intron_positive_strand['length']      = intron_positive_strand['End']   - intron_positive_strand['Start']   + 1
    exon_positive_strand['length']        = exon_positive_strand['End']     - exon_positive_strand['Start']     + 1
    Three_prime_UTR_positive['length']    = Three_prime_UTR_positive['End'] - Three_prime_UTR_positive['Start'] + 1
    Positive_strand_stop_codon['length']  = Positive_strand_stop_codon['End']  - Positive_strand_stop_codon['Start']  + 1
    Positive_strand_start_codon['length'] = Positive_strand_start_codon['End'] - Positive_strand_start_codon['Start'] + 1
    Positive_strand_transcript['length']  = Positive_strand_transcript['End']  - Positive_strand_transcript['Start']  + 1
    five_prime_UTR_positive['length']     = five_prime_UTR_positive['End']     - five_prime_UTR_positive['Start']     + 1

    Positive_strand_transcript["protein_sequence"] = \
        Positive_strand_transcript["Reference_id"].map(protein_dict)

    # ── Negative Strand ──
    Negative_strand_transcript  = Negative_strand[Negative_strand['Portion'] == 'transcript'].copy()
    Negative_strand_start_codon = Negative_strand[Negative_strand['Portion'] == 'start_codon'].copy()
    Negative_strand_stop_codon  = Negative_strand[Negative_strand['Portion'] == 'stop_codon'].copy()

    CDS_negative_strand    = Negative_strand[Negative_strand['Portion'] == 'CDS'].copy()
    intron_negative_strand = Negative_strand[Negative_strand['Portion'] == 'intron'].copy()
    exon_negative_strand   = Negative_strand[Negative_strand['Portion'] == 'exon'].copy()

    CDS_negative_strand['length']    = CDS_negative_strand['End']    - CDS_negative_strand['Start']    + 1
    intron_negative_strand['length'] = intron_negative_strand['End'] - intron_negative_strand['Start'] + 1
    exon_negative_strand['length']   = exon_negative_strand['End']   - exon_negative_strand['Start']   + 1

    # 5' UTR Negative strand
    Five_prime_UTR_negative_start = (Negative_strand_start_codon['End'] + 1).reset_index(drop=True).to_frame(name="Start")
    Five_prime_UTR_negative_End   = Negative_strand_transcript['End'].reset_index(drop=True).to_frame(name="End")
    Five_UTR_negative             = Negative_strand_transcript.drop(columns=['Start', 'End']).reset_index(drop=True)
    Five_UTR_negative_strand      = pd.concat([Five_UTR_negative, Five_prime_UTR_negative_start, Five_prime_UTR_negative_End], axis=1)
    Five_UTR_negative_strand['Portion'] = 'five_Prime_UTR'
    Five_UTR_negative_strand['length']  = Five_UTR_negative_strand['End'] - Five_UTR_negative_strand['Start'] + 1

    # 3' UTR Negative strand
    Three_prime_UTR_negative_end   = (Negative_strand_stop_codon['Start'] - 1).reset_index(drop=True).to_frame(name="End")
    Three_prime_UTR_negative_start = Negative_strand_transcript['Start'].reset_index(drop=True).to_frame(name="Start")
    Three_UTR_negative             = Negative_strand_transcript.drop(columns=['Start', 'End']).reset_index(drop=True)
    Three_UTR_negative_strand      = pd.concat([Three_UTR_negative, Three_prime_UTR_negative_start, Three_prime_UTR_negative_end], axis=1)
    Three_UTR_negative_strand['Portion'] = 'three_Prime_UTR'
    Three_UTR_negative_strand['length']  = Three_UTR_negative_strand['End'] - Three_UTR_negative_strand['Start'] + 1

    Negative_strand_transcript['length'] = Negative_strand_transcript['End'] - Negative_strand_transcript['Start'] + 1
    Negative_strand_transcript = Negative_strand_transcript.rename(columns={
        "Start": "Genomic_Start(tts)",
        "End":   "Genomic_End(tss)"
    })
    Negative_strand_transcript["protein_sequence"] = \
        Negative_strand_transcript["Reference_id"].map(protein_dict)

    # ── Promoters ──
    promoter_positive = pd.DataFrame({
        'ID':           five_prime_UTR_positive['ID'],
        'Source':       five_prime_UTR_positive['Source'],
        'Portion':      'promoter',
        'Start':        five_prime_UTR_positive['Start'] - PROMOTER_LEN,
        'End':          five_prime_UTR_positive['Start'] - 1,
        'score':        five_prime_UTR_positive['score'],
        'Strand':       five_prime_UTR_positive['Strand'],
        'Phase':        '.',
        'Reference_id': five_prime_UTR_positive['Reference_id']
    })
    promoter_positive['length'] = promoter_positive['End'] - promoter_positive['Start'] + 1

    Promoter_negative_start  = (Negative_strand_transcript['Genomic_End(tss)'] + 1).reset_index(drop=True).to_frame(name="Start")
    Promoter_negative_end    = (Negative_strand_transcript['Genomic_End(tss)'] + PROMOTER_LEN).reset_index(drop=True).to_frame(name="End")
    Promoter_negative        = Negative_strand_transcript.drop(columns=['Genomic_Start(tts)', 'Genomic_End(tss)'], errors="ignore").reset_index(drop=True)
    Promoter_negative_strand = pd.concat([Promoter_negative, Promoter_negative_start, Promoter_negative_end], axis=1)
    Promoter_negative_strand['Portion'] = 'Promoter'
    Promoter_negative_strand = Promoter_negative_strand.drop(columns=['protein_sequence'], errors='ignore')
    Promoter_negative_strand['length'] = Promoter_negative_strand['End'] - Promoter_negative_strand['Start'] + 1
    Promoter_negative_strand = Promoter_negative_strand[
        ['ID', 'Source', 'Portion', 'Start', 'End',
         'score', 'Strand', 'Phase', 'Reference_id', 'length']
    ]

    # NO CHANGE: Write transcript CSVs with pandas — same as original
    Positive_strand_transcript.to_csv(
        os.path.join(output_dir, "Positive_strand_transcript.csv"), index=False)
    Negative_strand_transcript.to_csv(
        os.path.join(output_dir, "Negative_strand_transcript.csv"), index=False)
    print("Transcript CSVs written.")

    # ──────────────────────────────────────────────────────────────────────────
    # CHANGE 3: Bismark file reading with Spark
    #
    # Original pandas logic (PRESERVED):
    #   - skiprows=1  → skips "Bismark methylation extractor version v0.25.1"
    #   - filter Methylation_condition == '+'
    #   - drop_duplicates on (Reference_id, Site)
    #
    # Bismark file format (from actual file):
    #   Line 1: "Bismark methylation extractor version v0.25.1"  ← skip this
    #   Col 1:  read_id       (Bismark_ID)
    #   Col 2:  strand        (Methylation_condition: + or -)
    #   Col 3:  chromosome    (Reference_id)
    #   Col 4:  position      (Site)
    #   Col 5:  context type  (Methyl_alphabet: z/Z/x/X/h/H)
    #
    # Spark change:
    #   skiprows=1 → filter out row where Bismark_ID starts with "Bismark"
    #   Everything else identical
    # ──────────────────────────────────────────────────────────────────────────
    bismark_schema = StructType([
        StructField("Bismark_ID",            StringType(), True),
        StructField("Methylation_condition", StringType(), True),
        StructField("Reference_id",          StringType(), True),
        StructField("Site",                  LongType(),   True),
        StructField("Methyl_alphabet",       StringType(), True),
    ])

    def read_bismark_spark(path):
        """
        Read Bismark methylation file into Spark DataFrame.
        Equivalent to original text_to_df() but distributed.

        Original:
            df = pd.read_csv(text, sep="\t", skiprows=1, header=None)
            df.columns = columns
            final_df = df[df['Methylation_condition']=='+'].copy()
            final_deduplicate_df = final_df.drop_duplicates(subset=["Reference_id","Site"])

        Spark equivalent below — same logic, distributed:
        """
        sdf = spark.read.csv(
            path,
            sep="\t",
            header=False,      # no column header in Bismark files
            schema=bismark_schema
        )

        # Replaces skiprows=1:
        # Skip the version header line "Bismark methylation extractor version..."
        sdf = sdf.filter(~F.col("Bismark_ID").startswith("Bismark"))

        # Identical to original: filter Methylation_condition == '+'
        sdf = sdf.filter(F.col("Methylation_condition") == "+")

        # Identical to original: drop_duplicates on Reference_id and Site
        sdf = sdf.dropDuplicates(["Reference_id", "Site"])

        return sdf

    print("Reading Bismark files...")
    df_CpG_OT = read_bismark_spark(CpG_OT)
    df_CHG_OT = read_bismark_spark(CHG_OT)
    df_CHH_OT = read_bismark_spark(CHH_OT)
    df_CpG_OB = read_bismark_spark(CpG_OB)
    df_CHG_OB = read_bismark_spark(CHG_OB)
    df_CHH_OB = read_bismark_spark(CHH_OB)
    print("Bismark files loaded.")

    # ── Region dictionaries — identical structure to original ──
    region_Positive_dfs = {
        "CDS_Positive":             CDS_Positive_strand,
        "intron_Positive":          intron_positive_strand,
        "exon_Positive":            exon_positive_strand,
        "Three_prime_UTR_Positive": Three_prime_UTR_positive,
        "Stop_codon_Positive":      Positive_strand_stop_codon,
        "Start_codon_Positive":     Positive_strand_start_codon,
        "Five_prime_UTR_Positive":  five_prime_UTR_positive,
        "Promoter_Positive":        promoter_positive,
    }
    site_Positive_dfs = {
        "CpG_OT": df_CpG_OT,
        "CHG_OT": df_CHG_OT,
        "CHH_OT": df_CHH_OT,
    }
    region_Negative_dfs = {
        "Three_prime_UTR_Negative": Three_UTR_negative_strand,
        "Five_prime_UTR_Negative":  Five_UTR_negative_strand,
        "exon_Negative":            exon_negative_strand,
        "intron_Negative":          intron_negative_strand,
        "CDS_Negative":             CDS_negative_strand,
        "Start_codon_Negative":     Negative_strand_start_codon,
        "Stop_codon_Negative":      Negative_strand_stop_codon,
        "Promoter_Negative":        Promoter_negative_strand,
    }
    site_Negative_dfs = {
        "CpG_OB": df_CpG_OB,
        "CHG_OB": df_CHG_OB,
        "CHH_OB": df_CHH_OB,
    }

    # ──────────────────────────────────────────────────────────────────────────
    # CHANGE 4: add_all_methylation_counts replaced with Spark range join
    #
    # Original logic (SLOW — triple nested loop):
    #   for each region row:
    #       count = (
    #           (site_df["Reference_id"] == row["ID"]) &
    #           (site_df["Site"] >= row["Start"]) &
    #           (site_df["Site"] <= row["End"])
    #       ).sum()
    #
    # Spark logic (FAST — distributed range join):
    #   Same conditions expressed as a join:
    #       region.ID == site.Reference_id     ← chromosome match
    #       site.Site >= region.Start          ← lower bound
    #       site.Site <= region.End            ← upper bound
    #   Then groupBy all region columns + count matching sites
    #
    # Result: identical counts, computed in parallel across all 20 cores
    # ──────────────────────────────────────────────────────────────────────────
    def add_all_methylation_counts_spark(region_dfs, site_dfs, output_dir,
                                          start_col="Start", end_col="End",
                                          site_col="Site"):
        output = {}

        for region_name, region_pd in region_dfs.items():
            print(f"\nProcessing: {region_name}")

            # Convert pandas region DataFrame to Spark
            region_sdf = spark.createDataFrame(region_pd)
            region_cols = region_pd.columns.tolist()

            for site_name, site_sdf in site_dfs.items():
                count_col = f"{site_name}_Count"
                print(f"  Counting {site_name}...")

                # Alias to avoid ambiguous column names
                r = region_sdf.alias("r")
                s = site_sdf.alias("s")

                # Range join — replaces iterrows() loop
                # Identical condition to original:
                #   site.Reference_id == region.ID         (chromosome match)
                #   site.Site >= region.Start              (within region)
                #   site.Site <= region.End                (within region)
                joined = r.join(
                    s,
                    on=(
                        (F.col("r.ID")             == F.col("s.Reference_id")) &
                        (F.col(f"s.{site_col}")    >= F.col(f"r.{start_col}")) &
                        (F.col(f"s.{site_col}")    <= F.col(f"r.{end_col}"))
                    ),
                    how="left"  # keep all regions even with zero sites (count=0)
                )

                # groupBy all region columns → count matched sites
                # F.count on site column: null (no match) → 0 automatically
                region_sdf = joined.groupBy(
                    *[F.col(f"r.{c}").alias(c) for c in region_cols]
                ).agg(
                    F.count(F.col(f"s.{site_col}")).alias(count_col)
                )

                # Update region_cols to include new count column
                region_cols = region_sdf.columns

            # CHANGE: Write as flat CSV using helper
            # Identical output to pandas .to_csv(path, index=False)
            csv_path = os.path.join(output_dir, f"{region_name}_methylation_counts.csv")
            save_as_flat_csv(region_sdf, csv_path)

            output[region_name] = region_sdf

        return output

    # ── Run for both strands ──
    print("\n=== Processing Positive Strand ===")
    positive_results = add_all_methylation_counts_spark(
        region_Positive_dfs, site_Positive_dfs, output_dir
    )

    print("\n=== Processing Negative Strand ===")
    negative_results = add_all_methylation_counts_spark(
        region_Negative_dfs, site_Negative_dfs, output_dir
    )

    print("\nAll done. Files written to:", output_dir)

    # Return structure identical to original
    # Values are Spark DataFrames instead of pandas
    # To convert: positive_results["CDS_Positive"].toPandas()
    return {
        "Positive": positive_results,
        "Negative": negative_results
    }


# ──────────────────────────────────────────────────────────────────────────────
# Entry point — identical arguments to original
# ──────────────────────────────────────────────────────────────────────────────
GFF_Bismark_annotation(
    "/cabinfs/home/subham9/Nipponbare/BS_seq/MOCK/SRR11777285/Output/nipponbare_annotation_with_UTR.gff3",
    1000,
    "/cabinfs/home/subham9/Nipponbare/BS_seq/MOCK/SRR11777285/Output/Spark_sample/",
    "/cabinfs/home/subham9/Nipponbare/BS_seq/MOCK/SRR11777285/Output/CpG_OT_SRR11777285__bismark_bt2_pe_deduplicated_namesorted.txt",
    "/cabinfs/home/subham9/Nipponbare/BS_seq/MOCK/SRR11777285/Output/CHG_OT_SRR11777285__bismark_bt2_pe_deduplicated_namesorted.txt",
    "/cabinfs/home/subham9/Nipponbare/BS_seq/MOCK/SRR11777285/Output/CHH_OT_SRR11777285__bismark_bt2_pe_deduplicated_namesorted.txt",
    "/cabinfs/home/subham9/Nipponbare/BS_seq/MOCK/SRR11777285/Output/CpG_OB_SRR11777285__bismark_bt2_pe_deduplicated_namesorted.txt",
    "/cabinfs/home/subham9/Nipponbare/BS_seq/MOCK/SRR11777285/Output/CHG_OB_SRR11777285__bismark_bt2_pe_deduplicated_namesorted.txt",
    "/cabinfs/home/subham9/Nipponbare/BS_seq/MOCK/SRR11777285/Output/CHH_OB_SRR11777285__bismark_bt2_pe_deduplicated_namesorted.txt"
)
