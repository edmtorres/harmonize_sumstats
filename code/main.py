from io_utils import load_ref_data, load_gwas_data, write_output
from liftover import perform_liftover
from join_utils import join_and_filter
from transformation import apply_transformations
from constants import OUTPUT_FILENAME
import argparse, os
import polars as pl

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome-build", type=str, required=True, choices=["hg19", "hg38"], help="Genome build of the GWAS file (e.g: hg19, hg38)")
    parser.add_argument("--ip-path", type=str, required=True, help="File path to input GWAS dataset")
    parser.add_argument("--op-path", type=str, required=True, help="File path to save output file")

    # Column name arguments
    parser.add_argument("--col-chr", required=True, help="Column name for chromosome")
    parser.add_argument("--col-pos", required=True, help="Column name for base pair location")
    parser.add_argument("--col-ea", required=True, help="Column name for effect allele")
    parser.add_argument("--col-oa", required=True, help="Column name for other allele")
    parser.add_argument("--col-effect", required=True, help="Column name for beta/OR")
    parser.add_argument("--col-eaf", required=True, help="Column name for effect allele frequency")
    parser.add_argument("--col-vid", required=True, help="Column name for variant ID")

    # Effect size type: beta or odds ratio
    parser.add_argument("--effect-type", type=str, required=True, choices=["beta", "OR"], help="'beta' or 'OR' (odds ratio) present")

    args = parser.parse_args()

    # Load the reference dataset
    ref_df = load_ref_data()
    print(f"Reference SNP list preview: {ref_df.head()}")

    # Load the GWAS dataset
    gwas_df = load_gwas_data(args.ip_path)
    print(f"GWAS summary statistics preview: {gwas_df.head()}")
    #column_names = gwas_df.columns
    #print(column_names)

    # If odds ratio is present, convert it to beta
    if args.effect_type == "OR":
        print("Converting to odds ratio to beta...")
        if "beta" in gwas_df.columns:
            print("Existing 'beta' column found. Dropping it before conversion...")
            # Drop the existing 'beta' column
            gwas_df = gwas_df.drop("beta")
        gwas_df = gwas_df.with_columns(
            (pl.col(args.col_effect).log().round(decimals=3)).alias("beta")  # log(OR) = beta
        )
        args.col_effect="beta"

    # Prepare for Join and Liftover: Cast 'chromosome' column in 'gwas_df' to Utf8 
    # This ensures matching data types for join keys between gwas_df and ref_df.
    gwas_df = gwas_df.with_columns(pl.col(args.col_chr).cast(pl.Utf8))

    # If the input genome build is not hg38, perform the liftover
    if args.genome_build != "hg38":
        gwas_df = perform_liftover(gwas_df, args.genome_build, args.col_chr, args.col_pos)
        print(f"GWAS DataFrame after liftover {gwas_df.head()}")
    else:
        # If no liftover, ensure '_hg38' column is present and same as 'base_pair_location'
        gwas_df = gwas_df.with_columns(pl.col(args.col_pos).alias("_hg38"))

    # Perform inner join and filtering of the DataFrames
    joined_df = join_and_filter(gwas_df, ref_df, args.col_chr, args.col_ea, args.col_oa)
    #print(gwas_df[args.col_vid, args.col_effect, args.col_effect])

    # Apply the transformations to beta and frequency when necessary
    joined_df = apply_transformations(joined_df, args.col_ea, args.col_oa, args.col_eaf, args.col_effect, args.col_vid)

    # Drop unnecessary columns
    joined_df = joined_df.drop([args.col_pos, "ref", "alt", "id"])
    print(f"Updated GWAS DataFrame preview: {joined_df.head()}")
    #print(gwas_df[args.col_vid, args.col_effect, args.col_effect])

    # Write the Output DataFrame to a .tsv file 
    # Ensure the output directory exists
    os.makedirs(args.op_path, exist_ok=True)
    write_output(joined_df, os.path.join(args.op_path, OUTPUT_FILENAME))

if __name__ == "__main__":
    main()