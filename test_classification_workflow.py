#!/usr/bin/env python3
"""
Test script for the classification workflow.
Runs all sequences through the classification pipeline and compares results to existing classifications.
"""

import os
import sys
import tempfile
import pandas as pd
from pathlib import Path
from typing import Dict, Optional
from dataclasses import dataclass
from collections import defaultdict
import json
import argparse
from datetime import datetime

# Add src to path
sys.path.append(os.path.join(os.path.dirname(__file__), "src"))

from src.utils.test_helpers import (
    test_fetch_ships as fetch_ships,
    test_fetch_meta_data as fetch_meta_data,
    test_run_classification_workflow as run_classification_workflow,
)
from src.config.logging import get_logger

logger = get_logger(__name__)


@dataclass
class FetchParams:
    """Parameters for fetching data from database"""

    curated: bool = False
    with_sequence: bool = True
    dereplicate: bool = True


@dataclass
class UploadData:
    """Mock upload data structure for testing"""

    fasta_file: str
    seq_type: str = "nucl"
    fetch_ship_params: FetchParams = None
    fetch_captain_params: FetchParams = None
    blast_df: Optional[pd.DataFrame] = None

    def __post_init__(self):
        if self.fetch_ship_params is None:
            self.fetch_ship_params = FetchParams()
        if self.fetch_captain_params is None:
            self.fetch_captain_params = FetchParams()


class ClassificationTester:
    """Test the classification workflow accuracy for family, navis, and haplotype predictions"""

    def __init__(
        self,
        output_dir: str = None,
        curated_only: bool = False,
        max_sequences: int = None,
    ):
        self.curated_only = curated_only
        self.max_sequences = max_sequences
        self.output_dir = (
            Path(output_dir) if output_dir else Path("classification_test_results")
        )
        self.output_dir.mkdir(exist_ok=True)

        # Results tracking
        self.results = []
        self.classification_stats = defaultdict(int)
        self.confusion_matrices = {}

        # Create timestamp for this test run
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        # classification stages we want to test
        self.stages = {"family", "navis", "haplotype"}

    def fetch_test_data(self) -> pd.DataFrame:
        """Fetch sequences and their existing classifications for testing"""
        logger.info("Fetching test data from database...")

        # Fetch ships with sequences and metadata
        ships_df = fetch_ships(
            curated=self.curated_only, with_sequence=True, dereplicate=True
        )

        if ships_df.empty:
            raise ValueError("No ships data found in database")

        logger.info(f"Fetched {len(ships_df)} ships for testing")

        # Fetch additional metadata
        accession_tags = ships_df["accession_tag"].tolist()
        meta_df = fetch_meta_data(
            curated=self.curated_only, accession_tag=accession_tags
        )

        # Merge ships and metadata
        test_df = ships_df.merge(
            meta_df[["accession_tag", "familyName", "navis_name", "haplotype_name"]],
            on="accession_tag",
            how="left",
            suffixes=("", "_meta"),
        )

        # Remove rows without sequences
        test_df = test_df.dropna(subset=["sequence"])

        if self.max_sequences:
            test_df = test_df.head(self.max_sequences)
            logger.info(f"Limited to {len(test_df)} sequences for testing")

        return test_df

    def create_temp_fasta(self, sequence: str, accession: str) -> str:
        """Create temporary FASTA file for a sequence"""
        temp_file = tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False)
        temp_file.write(f">{accession}\n{sequence}\n")
        temp_file.close()
        return temp_file.name

    def extract_classification_result(self, workflow_result: Dict) -> Dict:
        """Extract classification information from workflow result"""
        classification = {
            "found_match": workflow_result.get("found_match", False),
            "match_stage": workflow_result.get("match_stage"),
            "match_result": workflow_result.get("match_result"),
            "complete": workflow_result.get("complete", False),
            "error": workflow_result.get("error"),
            "is_classification": False,
            "predicted_family": None,
            "predicted_navis": None,
            "predicted_haplotype": None,
        }

        # Check if we have classification results in class_dict
        class_dict = workflow_result.get("class_dict", {})
        if class_dict:
            classification["is_classification"] = True

            # Extract family classification
            if "family" in class_dict:
                classification["predicted_family"] = class_dict["family"]

            # Extract navis classification
            if "navis" in class_dict:
                classification["predicted_navis"] = class_dict["navis"]

            # Extract haplotype classification
            if "haplotype" in class_dict:
                classification["predicted_haplotype"] = class_dict["haplotype"]

            # Set the stage to the most specific classification found
            if classification["predicted_haplotype"]:
                classification["match_stage"] = "haplotype"
            elif classification["predicted_navis"]:
                classification["match_stage"] = "navis"
            elif classification["predicted_family"]:
                classification["match_stage"] = "family"

        return classification

    def compare_classifications(self, actual: Dict, predicted: Dict) -> Dict:
        """Compare actual vs predicted classifications"""
        comparison = {
            "is_testable": predicted.get("is_classification", False),
            "family_correct": None,
            "navis_correct": None,
            "haplotype_correct": None,
            "any_correct": False,
        }

        # Only evaluate if this is an classification
        if not comparison["is_testable"]:
            return comparison

        # Compare family
        actual_family = actual.get("familyName") or actual.get("familyName_meta")
        predicted_family = predicted.get("predicted_family")
        if actual_family and predicted_family:
            comparison["family_correct"] = actual_family == predicted_family

        # Compare navis
        actual_navis = actual.get("navis_name") or actual.get("navis_name_meta")
        predicted_navis = predicted.get("predicted_navis")
        if actual_navis and predicted_navis:
            comparison["navis_correct"] = actual_navis == predicted_navis

        # Compare haplotype
        actual_haplotype = actual.get("haplotype_name") or actual.get(
            "haplotype_name_meta"
        )
        predicted_haplotype = predicted.get("predicted_haplotype")
        if actual_haplotype and predicted_haplotype:
            comparison["haplotype_correct"] = actual_haplotype == predicted_haplotype

        # Any classification correct
        comparison["any_correct"] = any(
            [
                comparison["family_correct"] is True,
                comparison["navis_correct"] is True,
                comparison["haplotype_correct"] is True,
            ]
        )

        return comparison

    def run_single_test(self, row: pd.Series) -> Dict:
        """Run classification workflow on a single sequence"""
        accession = row["accession_tag"]
        sequence = row["sequence"]

        logger.debug(f"Testing sequence {accession}")

        # Create temporary FASTA file
        fasta_path = self.create_temp_fasta(sequence, accession)

        try:
            # Create upload data object
            upload_data = UploadData(
                fasta_file=fasta_path,
                seq_type="nucl",
                fetch_ship_params=FetchParams(
                    curated=self.curated_only, with_sequence=True
                ),
                fetch_captain_params=FetchParams(
                    curated=self.curated_only, with_sequence=True
                ),
            )

            # Run classification workflow
            workflow_result = run_classification_workflow(
                upload_data=upload_data, meta_dict=None, stages=self.stages
            )

            # Extract classification results
            predicted = self.extract_classification_result(workflow_result)
            actual = row.to_dict()

            # Compare results
            comparison = self.compare_classifications(actual, predicted)

            # Compile test result
            test_result = {
                "accession_tag": accession,
                "sequence_length": len(sequence),
                "actual_family": actual.get("familyName")
                or actual.get("familyName_meta"),
                "actual_navis": actual.get("navis_name")
                or actual.get("navis_name_meta"),
                "actual_haplotype": actual.get("haplotype_name")
                or actual.get("haplotype_name_meta"),
                "match_stage": predicted.get("match_stage"),
                "predicted_family": predicted.get("predicted_family"),
                "predicted_navis": predicted.get("predicted_navis"),
                "predicted_haplotype": predicted.get("predicted_haplotype"),
                "found_match": predicted.get("found_match"),
                "is_classification": predicted.get("is_classification"),
                "workflow_complete": predicted.get("complete"),
                "workflow_error": predicted.get("error"),
                "family_correct": comparison["family_correct"],
                "navis_correct": comparison["navis_correct"],
                "haplotype_correct": comparison["haplotype_correct"],
                "any_correct": comparison["any_correct"],
            }

            return test_result

        except Exception as e:
            logger.error(f"Error testing sequence {accession}: {str(e)}")
            return {
                "accession_tag": accession,
                "sequence_length": len(sequence),
                "actual_family": None,
                "actual_navis": None,
                "actual_haplotype": None,
                "match_stage": None,
                "predicted_family": None,
                "predicted_navis": None,
                "predicted_haplotype": None,
                "found_match": False,
                "is_classification": False,
                "workflow_complete": False,
                "workflow_error": str(e),
                "family_correct": None,
                "navis_correct": None,
                "haplotype_correct": None,
                "any_correct": False,
            }

        finally:
            # Clean up temporary file
            try:
                os.unlink(fasta_path)
            except OSError:
                pass

    def calculate_metrics(self) -> Dict:
        """Calculate performance metrics for classifications only"""
        if not self.results:
            return {}

        results_df = pd.DataFrame(self.results)

        # Filter to only classifications
        results = results_df[results_df["is_classification"]]

        metrics = {
            "total_sequences": len(results_df),
            "classified_sequences": len(results),
            "non_sequences": len(results_df) - len(results),
            "workflow_errors": len(results_df[results_df["workflow_error"].notna()]),
            "failed_classifications": len(results_df[~results_df["found_match"]]),
        }

        if len(results) > 0:
            # ML-specific metrics
            metrics["any_correct_accuracy"] = results["any_correct"].mean()

            # Stage-specific metrics for classifications
            stage_counts = results["match_stage"].value_counts(dropna=False)
            metrics["classification_by_stage"] = stage_counts.to_dict()

            # Accuracy by classification type (only for classifications)
            for classification_type in ["family", "navis", "haplotype"]:
                correct_col = f"{classification_type}_correct"
                if correct_col in results.columns:
                    correct_results = results[results[correct_col].notna()]
                    if len(correct_results) > 0:
                        metrics[f"{classification_type}_accuracy"] = correct_results[
                            correct_col
                        ].mean()
                        metrics[f"{classification_type}_tested"] = len(correct_results)
        else:
            metrics["any_correct_accuracy"] = 0
            metrics["classification_by_stage"] = {}

        return metrics

    def generate_confusion_matrix(
        self, actual_col: str, predicted_col: str
    ) -> pd.DataFrame:
        """Generate confusion matrix for classifications only"""
        if not self.results:
            return pd.DataFrame()

        results_df = pd.DataFrame(self.results)

        # Filter to only classifications with both actual and predicted values
        results = results_df[results_df["is_classification"]]
        valid_results = results[
            (results[actual_col].notna()) & (results[predicted_col].notna())
        ]

        if valid_results.empty:
            return pd.DataFrame()

        return pd.crosstab(
            valid_results[actual_col],
            valid_results[predicted_col],
            margins=True,
            margins_name="Total",
        )

    def save_results(self):
        """Save classification test results to files"""
        # Save detailed results
        results_df = pd.DataFrame(self.results)
        results_file = self.output_dir / f"detailed_results_{self.timestamp}.csv"
        results_df.to_csv(results_file, index=False)
        logger.info(f"Saved detailed results to {results_file}")

        # Save ML-only results
        results = results_df[results_df["is_classification"]]
        results_file = self.output_dir / f"only_results_{self.timestamp}.csv"
        results.to_csv(results_file, index=False)
        logger.info(f"Saved ML-only results to {results_file}")

        # Calculate and save metrics
        metrics = self.calculate_metrics()
        metrics_file = self.output_dir / f"metrics_{self.timestamp}.json"
        with open(metrics_file, "w") as f:
            json.dump(metrics, f, indent=2, default=str)
        logger.info(f"Saved metrics to {metrics_file}")

        # Generate and save confusion matrices for classifications
        confusion_matrices = {}
        for actual_col, predicted_col in [
            ("actual_family", "predicted_family"),
            ("actual_navis", "predicted_navis"),
            ("actual_haplotype", "predicted_haplotype"),
        ]:
            cm = self.generate_confusion_matrix(actual_col, predicted_col)
            if not cm.empty:
                classification_type = actual_col.replace("actual_", "")
                confusion_matrices[classification_type] = cm
                cm_file = (
                    self.output_dir
                    / f"confusion_matrix_{classification_type}_{self.timestamp}.csv"
                )
                cm.to_csv(cm_file)
                logger.info(
                    f"Saved {classification_type} confusion matrix to {cm_file}"
                )

        # Save summary report
        self.generate_summary_report(metrics, confusion_matrices)

    def generate_summary_report(self, metrics: Dict, confusion_matrices: Dict):
        """Generate a human-readable summary report for classifications"""
        report_file = self.output_dir / f"summary_report_{self.timestamp}.txt"

        with open(report_file, "w") as f:
            f.write("CLASSIFICATION WORKFLOW TEST REPORT\n")
            f.write("=" * 50 + "\n")
            f.write(f"Test run: {self.timestamp}\n")
            f.write(f"Curated only: {self.curated_only}\n")
            f.write(f"Max sequences: {self.max_sequences or 'All'}\n\n")

            f.write("OVERALL METRICS\n")
            f.write("-" * 20 + "\n")
            f.write(f"Total sequences tested: {metrics.get('total_sequences', 0)}\n")
            f.write(f"classified sequences: {metrics.get('classified_sequences', 0)}\n")
            f.write(
                f"Non-sequences (exact/contained/similar): {metrics.get('non_sequences', 0)}\n"
            )
            f.write(
                f"Failed classifications: {metrics.get('failed_classifications', 0)}\n"
            )
            f.write(f"Workflow errors: {metrics.get('workflow_errors', 0)}\n")
            f.write(
                f"classification accuracy: {metrics.get('any_correct_accuracy', 0):.2%}\n\n"
            )

            f.write("CLASSIFICATION BY STAGE\n")
            f.write("-" * 30 + "\n")
            stage_counts = metrics.get("classification_by_stage", {})
            for stage, count in stage_counts.items():
                f.write(f"{stage}: {count}\n")
            f.write("\n")

            f.write("ACCURACY BY CLASSIFICATION TYPE\n")
            f.write("-" * 40 + "\n")
            for classification_type in ["family", "navis", "haplotype"]:
                accuracy_key = f"{classification_type}_accuracy"
                tested_key = f"{classification_type}_tested"
                if accuracy_key in metrics:
                    accuracy = metrics[accuracy_key]
                    tested = metrics.get(tested_key, 0)
                    f.write(
                        f"{classification_type.title()}: {accuracy:.2%} ({tested} tested)\n"
                    )
            f.write("\n")

            f.write("CONFUSION MATRICES\n")
            f.write("-" * 25 + "\n")
            for classification_type, cm in confusion_matrices.items():
                f.write(f"\n{classification_type.title()} Classification:\n")
                f.write(str(cm))
                f.write("\n")

        logger.info(f"Saved summary report to {report_file}")

    def run_tests(self):
        """Run the complete classification test suite"""
        logger.info("Starting classification workflow tests...")

        # Fetch test data
        test_df = self.fetch_test_data()

        logger.info(f"Running classification tests on {len(test_df)} sequences...")

        # Run tests
        for idx, row in test_df.iterrows():
            try:
                result = self.run_single_test(row)
                self.results.append(result)

                if (idx + 1) % 10 == 0:
                    logger.info(f"Completed {idx + 1}/{len(test_df)} tests")

            except Exception as e:
                logger.error(f"Error in test {idx + 1}: {str(e)}")
                continue

        logger.info(f"Completed all tests. Processed {len(self.results)} sequences.")

        # Save results
        self.save_results()

        # Print summary
        metrics = self.calculate_metrics()
        print("\n" + "=" * 50)
        print("CLASSIFICATION TEST SUMMARY")
        print("=" * 50)
        print(f"Total sequences: {metrics.get('total_sequences', 0)}")
        print(f"classified sequences: {metrics.get('classified_sequences', 0)}")
        print(
            f"Non-sequences (exact/contained/similar): {metrics.get('non_sequences', 0)}"
        )
        print(f"classification accuracy: {metrics.get('any_correct_accuracy', 0):.2%}")
        print(f"Results saved to: {self.output_dir}")


def main():
    parser = argparse.ArgumentParser(description="Test classification workflow")
    parser.add_argument(
        "--output-dir",
        default="classification_test_results",
        help="Output directory for results",
    )
    parser.add_argument(
        "--curated-only", action="store_true", help="Only test curated sequences"
    )
    parser.add_argument(
        "--max-sequences", type=int, help="Maximum number of sequences to test"
    )
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")

    args = parser.parse_args()

    if args.verbose:
        import logging

        logging.getLogger().setLevel(logging.DEBUG)

    # Run classification tests
    tester = ClassificationTester(
        output_dir=args.output_dir,
        curated_only=args.curated_only,
        max_sequences=args.max_sequences,
    )

    tester.run_tests()


if __name__ == "__main__":
    main()
