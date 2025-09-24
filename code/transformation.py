import polars as pl

def apply_transformations(df, col_ea, col_oa, col_eaf, col_effect, col_vid):
    '''
    Apply transformations (frequency, beta, allele swapping) only if the alleles need flipping.
    df: DataFrame to perform the transformations on
    '''
    return df.with_columns([
        # Changing 'effect_allele_frequency' based on the specified condition
        pl.when(
        (pl.col(col_ea) == pl.col("ref")) & # If GWAS effect_allele matches ref_df's ref allele
        (pl.col(col_oa) == pl.col("alt"))    # AND GWAS other allele matches ref_df's alt allele
    )
    .then(1 - pl.col(col_eaf))  # THEN, change the frequency
    .otherwise(pl.col(col_eaf)) # ELSE, keep original frequency
    .alias(col_eaf),

    # Changing 'beta' based on the specified condition
    pl.when(
        (pl.col(col_ea) == pl.col("ref")) &
        (pl.col(col_oa) == pl.col("alt"))
    )
    .then(-pl.col(col_effect))      # THEN, flip the sign of beta
    .otherwise(pl.col(col_effect))  # ELSE, keep original beta
    .alias(col_effect),

    # Changing 'effect_allele' based on the specified condition
    pl.when(
        (pl.col(col_ea) == pl.col("ref")) &
        (pl.col(col_oa) == pl.col("alt"))
    )
    .then(pl.col("alt"))                # THEN, set effect_allele to ref_df's alt
    .otherwise(pl.col(col_ea)) # ELSE, keep original effect_allele
    .alias(col_ea),

    # Changing 'other_allele' based on the specified condition
    pl.when(
        (pl.col(col_ea) == pl.col("ref")) &
        (pl.col(col_oa) == pl.col("alt"))
    )
    .then(pl.col("ref"))                # THEN, set other_allele to ref_df's ref
    .otherwise(pl.col(col_oa))  # ELSE, keep original other_allele
    .alias(col_oa),

    # Assign 'variant_id' directly from 'id' column from ref_df
    pl.col("id").alias(col_vid)
    ])