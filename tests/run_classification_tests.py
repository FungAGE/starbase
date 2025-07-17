#!/usr/bin/env python3
"""
Example script for running classification workflow tests.
This demonstrates different ways to run the tests with various options.
"""

import os
import sys
from pathlib import Path

from test_classification_workflow import ClassificationTester

# Add the main project directory to the path so we can import src modules
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)


def run_quick_test():
    """Run a quick test with just a few sequences"""
    print("Running quick test (10 sequences only)...")

    tester = ClassificationTester(
        output_dir="test_results_quick", curated_only=False, max_sequences=10
    )

    tester.run_tests()
    print("Quick test completed!")


def run_curated_test():
    """Run test on curated sequences only"""
    print("Running test on curated sequences only...")

    tester = ClassificationTester(
        output_dir="test_results_curated",
        curated_only=True,
        max_sequences=50,  # Limit to 50 for testing
    )

    tester.run_tests()
    print("Curated test completed!")


def run_full_test():
    """Run the complete test on all sequences (WARNING: This could take a very long time!)"""
    print("WARNING: This will run tests on ALL sequences in the database.")
    print("This could take hours or days depending on your database size!")

    response = input("Are you sure you want to continue? (yes/no): ")
    if response.lower() != "yes":
        print("Test cancelled.")
        return

    tester = ClassificationTester(
        output_dir="test_results_full",
        curated_only=False,
        max_sequences=None,  # No limit
    )

    tester.run_tests()
    print("Full test completed!")


def analyze_previous_results(results_dir: str):
    """Analyze results from a previous test run"""
    results_path = Path(results_dir)
    if not results_path.exists():
        print(f"Results directory {results_dir} does not exist!")
        return

    # Find the most recent results
    csv_files = list(results_path.glob("detailed_results_*.csv"))
    if not csv_files:
        print(f"No result files found in {results_dir}")
        return

    # Sort by modification time and get the most recent
    latest_file = max(csv_files, key=lambda x: x.stat().st_mtime)

    import pandas as pd

    # Load and analyze results
    results_df = pd.read_csv(latest_file)

    print(f"\nAnalyzing results from: {latest_file}")
    print(f"Total sequences tested: {len(results_df)}")
    print(f"Successful classifications: {results_df['found_match'].sum()}")

    # Check if we have the new status columns
    has_status_cols = all(
        col in results_df.columns
        for col in ["family_status", "navis_status", "haplotype_status"]
    )

    if has_status_cols:
        print("\nStage Status Summary:")
        for stage in ["family", "navis", "haplotype"]:
            status_col = f"{stage}_status"
            if status_col in results_df.columns:
                status_counts = results_df[status_col].value_counts()
                predicted = status_counts.get("predicted", 0)
                failed = status_counts.get("failed", 0)
                skipped = status_counts.get("skipped", 0)
                error = status_counts.get("error", 0)
                print(
                    f"  {stage.title()}: {predicted} predicted, {failed} failed, {skipped} skipped, {error} error"
                )
    else:
        # Legacy format - try to infer from match_stage
        if "match_stage" in results_df.columns:
            print("\nClassification by stage:")
            stage_counts = results_df["match_stage"].value_counts()
            for stage, count in stage_counts.items():
                print(f"  {stage}: {count}")

    # Accuracy by type
    print("\nAccuracy by classification type:")
    print("(Only calculated for stages that made predictions)")
    for stage in ["family", "navis", "haplotype"]:
        correct_col = f"{stage}_correct"
        status_col = f"{stage}_status"

        if correct_col in results_df.columns:
            if has_status_cols and status_col in results_df.columns:
                # New format: only count sequences where stage made a prediction
                predicted_results = results_df[results_df[status_col] == "predicted"]
                testable_results = predicted_results[
                    predicted_results[correct_col].notna()
                ]
            else:
                # Legacy format: count all non-null results
                testable_results = results_df[results_df[correct_col].notna()]

            if len(testable_results) > 0:
                accuracy = testable_results[correct_col].mean()
                total_predicted = (
                    len(results_df[results_df[status_col] == "predicted"])
                    if has_status_cols and status_col in results_df.columns
                    else len(testable_results)
                )
                print(
                    f"  {stage.title()}: {accuracy:.2%} ({len(testable_results)} testable / {total_predicted} predicted)"
                )
            else:
                print(f"  {stage.title()}: No testable results")

    # Overall accuracy
    if "any_correct" in results_df.columns:
        overall_accuracy = results_df["any_correct"].mean()
        print(f"\nOverall accuracy: {overall_accuracy:.2%}")
    else:
        print("\nOverall accuracy: Not available in this results file")


def main():
    """Main menu for running different types of tests"""
    print("Classification Workflow Testing")
    print("=" * 40)
    print("1. Quick test (10 sequences)")
    print("2. Curated sequences test (50 sequences)")
    print("3. Full database test (WARNING: Very slow!)")
    print("4. Analyze previous results")
    print("5. Exit")

    while True:
        choice = input("\nSelect an option (1-5): ").strip()

        if choice == "1":
            run_quick_test()
            break
        elif choice == "2":
            run_curated_test()
            break
        elif choice == "3":
            run_full_test()
            break
        elif choice == "4":
            results_dir = input("Enter results directory path: ").strip()
            analyze_previous_results(results_dir)
            break
        elif choice == "5":
            print("Goodbye!")
            break
        else:
            print("Invalid choice. Please select 1-5.")


if __name__ == "__main__":
    main()
