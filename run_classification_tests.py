#!/usr/bin/env python3
"""
Example script for running classification workflow tests.
This demonstrates different ways to run the tests with various options.
"""

import os
import sys
from pathlib import Path

# Add the main directory to the path
sys.path.append(os.path.dirname(__file__))

from test_classification_workflow import ClassificationTester


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
    print(f"Overall accuracy: {results_df['overall_correct'].mean():.2%}")

    # Stage breakdown
    print("\nClassification by stage:")
    stage_counts = results_df["predicted_stage"].value_counts()
    for stage, count in stage_counts.items():
        print(f"  {stage}: {count}")

    # Accuracy by type
    print("\nAccuracy by classification type:")
    for col in [
        "family_correct",
        "navis_correct",
        "haplotype_correct",
        "accession_correct",
    ]:
        if col in results_df.columns:
            correct_results = results_df[results_df[col].notna()]
            if len(correct_results) > 0:
                accuracy = correct_results[col].mean()
                print(
                    f"  {col.replace('_correct', '').title()}: {accuracy:.2%} ({len(correct_results)} tested)"
                )


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
