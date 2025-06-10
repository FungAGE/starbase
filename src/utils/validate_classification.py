import pandas as pd
from src.database.sql_manager import (
    fetch_captains,
    fetch_ships,
    fetch_meta_data,
)
from src.utils.seq_utils import write_temp_fasta
from src.utils.classification_utils import (
    check_exact_match,
    check_contained_match,
    check_similar_match,
    classify_family,
    classify_navis,
    classify_haplotype,
)


def run_checks(fasta_file, ships_df, accession_tag, seq_type, meta_dict, captains_df):
    (
        exact_match_result,
        contained_match_result,
        similar_match_result,
        similarities,
        family_dict,
        protein_file,
        classify_navis_result,
        classify_haplotype_result,
    ) = None, None, None, None, None, None, None, None

    exact_match_result = check_exact_match(fasta=fasta_file, existing_ships=ships_df)

    if exact_match_result is None:
        print(f"No exact match found for {accession_tag}")

        contained_match_result = check_contained_match(
            fasta=fasta_file,
            existing_ships=ships_df,
            min_coverage=0.95,
            min_identity=0.95,
        )

    if exact_match_result is None and contained_match_result is None:
        print(f"No contained match found for {accession_tag}")

        similar_match_result, similarities = check_similar_match(
            fasta=fasta_file,
            existing_ships=ships_df,
            threshold=0.95,
        )

    if (
        exact_match_result is None
        and contained_match_result is None
        and similar_match_result is None
    ):
        family_dict, protein_file = classify_family(
            fasta=fasta_file,
            seq_type=seq_type,
            meta_dict=meta_dict,
            pident_thresh=90,
            input_eval=0.001,
            threads=1,
        )

        if family_dict is not None:
            classify_navis_result = classify_navis(
                fasta=fasta_file,
                existing_captains=captains_df,
                threads=1,
            )

            if classify_navis_result is not None:
                classify_haplotype_result = classify_haplotype(
                    fasta=fasta_file,
                    existing_ships=ships_df,
                    navis=classify_navis_result,
                    similarities=similarities,
                )
    return (
        exact_match_result,
        contained_match_result,
        similar_match_result,
        similarities,
        family_dict,
        protein_file,
        classify_navis_result,
        classify_haplotype_result,
    )


def validate_classification():
    """
    Check if our classification pipeline can recreate the existing classification
    """
    meta_data = fetch_meta_data()
    ships_df = fetch_ships(curated=False, dereplicate=False, with_sequence=True)
    captains_df = fetch_captains(curated=False, dereplicate=False, with_sequence=True)
    # accession_tags = ships_df["accession_tag"].values

    if len(ships_df) == 0:
        return None

    results_df = pd.DataFrame(
        columns=[
            "accession_tag",
            "existing_family",
            "existing_navis",
            "existing_haplotype",
            "exact_match_result",
            "contained_match_result",
            "similar_match_result",
            "similarities",
            "family_dict",
            "protein_file",
            "classify_navis_result",
            "classify_haplotype_result",
        ]
    )

    for index, row in ships_df.iterrows():
        accession_tag = row["accession_tag"]

        meta_data_subset = meta_data[meta_data["accession_tag"] == accession_tag]
        meta_dict = meta_data_subset.to_dict(orient="records")

        # get existing classification information
        existing_family = meta_data_subset["familyName"].values[0]
        existing_navis = meta_data_subset["starship_navis"].values[0]
        existing_haplotype = meta_data_subset["starship_haplotype"].values[0]

        sequence = row["sequence"]
        fasta_file = write_temp_fasta(accession_tag, sequence)
        seq_type = "nucl"
        # seq_type = guess_seq_type(sequence)

        (
            exact_match_result,
            contained_match_result,
            similar_match_result,
            similarities,
            family_dict,
            protein_file,
            classify_navis_result,
            classify_haplotype_result,
        ) = run_checks(
            fasta_file,
            ships_df,
            accession_tag,
            seq_type,
            meta_dict,
            captains_df,
        )
        # append results to dataframe
        results_df = results_df.append(
            {
                "accession_tag": accession_tag,
                "existing_family": existing_family,
                "existing_navis": existing_navis,
                "existing_haplotype": existing_haplotype,
                "exact_match_result": exact_match_result,
                "contained_match_result": contained_match_result,
                "similar_match_result": similar_match_result,
                "similarities": similarities,
                "family_dict": family_dict,
                "protein_file": protein_file,
                "classify_navis_result": classify_navis_result,
                "classify_haplotype_result": classify_haplotype_result,
            },
            ignore_index=True,
        )

    return results_df


if __name__ == "__main__":
    validate_classification()
