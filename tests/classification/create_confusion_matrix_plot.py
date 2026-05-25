#!/usr/bin/env python3
"""
Create a confusion matrix visualization from CSV data.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os


def load_confusion_matrix(csv_path):
    """Load confusion matrix from CSV file."""
    df = pd.read_csv(csv_path, index_col=0)

    # Remove the 'Total' row and column if they exist
    if "Total" in df.index:
        df = df.drop("Total", axis=0)
    if "Total" in df.columns:
        df = df.drop("Total", axis=1)

    return df


def calculate_metrics(confusion_matrix):
    """Calculate precision, recall, and F1-score from confusion matrix."""
    # Convert to numpy array for easier calculation
    cm = confusion_matrix.values

    # Calculate per-class metrics
    precision = np.diag(cm) / np.sum(cm, axis=0)
    recall = np.diag(cm) / np.sum(cm, axis=1)
    f1_score = 2 * (precision * recall) / (precision + recall)

    # Handle division by zero
    precision = np.nan_to_num(precision)
    recall = np.nan_to_num(recall)
    f1_score = np.nan_to_num(f1_score)

    # Calculate overall accuracy
    accuracy = np.trace(cm) / np.sum(cm)

    return {
        "precision": precision,
        "recall": recall,
        "f1_score": f1_score,
        "accuracy": accuracy,
        "class_names": confusion_matrix.index.tolist(),
    }


def create_confusion_matrix_plot(
    confusion_matrix, save_path=None, title="Confusion Matrix"
):
    """Create a beautiful confusion matrix heatmap."""

    # Set up the matplotlib figure
    plt.figure(figsize=(12, 10))

    # Create the heatmap using matplotlib imshow
    im = plt.imshow(confusion_matrix.values, interpolation="nearest", cmap="Blues")

    # Add colorbar
    cbar = plt.colorbar(im, shrink=0.8)
    cbar.set_label("Count", rotation=270, labelpad=20)

    # Add text annotations
    thresh = confusion_matrix.values.max() / 2.0
    for i in range(confusion_matrix.shape[0]):
        for j in range(confusion_matrix.shape[1]):
            plt.text(
                j,
                i,
                format(confusion_matrix.values[i, j], "d"),
                ha="center",
                va="center",
                color="white" if confusion_matrix.values[i, j] > thresh else "black",
                fontsize=10,
            )

    # Set ticks and labels
    plt.xticks(range(len(confusion_matrix.columns)), confusion_matrix.columns)
    plt.yticks(range(len(confusion_matrix.index)), confusion_matrix.index)

    # Customize the plot
    plt.title(title, fontsize=16, fontweight="bold", pad=20)
    plt.xlabel("Predicted Family", fontsize=14, fontweight="bold")
    plt.ylabel("Actual Family", fontsize=14, fontweight="bold")

    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)

    # Adjust layout to prevent label cutoff
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Confusion matrix plot saved to: {save_path}")


def create_metrics_plot(metrics, save_path=None):
    """Create a bar plot showing precision, recall, and F1-score for each class."""

    class_names = metrics["class_names"]
    precision = metrics["precision"]
    recall = metrics["recall"]
    f1_score = metrics["f1_score"]

    # Set up the figure
    fig, ax = plt.subplots(figsize=(14, 8))

    # Set the width of bars and positions
    bar_width = 0.25
    r1 = np.arange(len(class_names))
    r2 = [x + bar_width for x in r1]
    r3 = [x + bar_width for x in r2]

    # Create bars
    bars1 = ax.bar(
        r1, precision, bar_width, label="Precision", alpha=0.8, color="skyblue"
    )
    bars2 = ax.bar(r2, recall, bar_width, label="Recall", alpha=0.8, color="lightgreen")
    bars3 = ax.bar(r3, f1_score, bar_width, label="F1-Score", alpha=0.8, color="salmon")

    # Customize the plot
    ax.set_xlabel("Family", fontweight="bold", fontsize=12)
    ax.set_ylabel("Score", fontweight="bold", fontsize=12)
    ax.set_title(
        f"Classification Metrics by Family (Overall Accuracy: {metrics['accuracy']:.3f})",
        fontweight="bold",
        fontsize=14,
    )
    ax.set_xticks([r + bar_width for r in range(len(class_names))])
    ax.set_xticklabels(class_names, rotation=45, ha="right")
    ax.legend()
    ax.set_ylim(0, 1.1)

    # Add value labels on bars
    def add_value_labels(bars):
        for bar in bars:
            height = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                height + 0.01,
                f"{height:.2f}",
                ha="center",
                va="bottom",
                fontsize=8,
            )

    add_value_labels(bars1)
    add_value_labels(bars2)
    add_value_labels(bars3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Metrics plot saved to: {save_path}")


def main():
    """Main function to create confusion matrix visualization."""
    parser = argparse.ArgumentParser(
        description="Create confusion matrix visualization"
    )
    parser.add_argument("input_csv", help="Path to the confusion matrix CSV file")
    parser.add_argument("--output-dir", default=".", help="Output directory for plots")
    parser.add_argument(
        "--title",
        default="Starship Family Classification Confusion Matrix",
        help="Title for the confusion matrix plot",
    )

    args = parser.parse_args()

    # Load the confusion matrix
    print(f"Loading confusion matrix from: {args.input_csv}")
    cm_df = load_confusion_matrix(args.input_csv)
    print(f"Confusion matrix shape: {cm_df.shape}")
    print(f"Classes: {cm_df.index.tolist()}")

    # Calculate metrics
    metrics = calculate_metrics(cm_df)
    print(f"Overall accuracy: {metrics['accuracy']:.3f}")

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Generate base filename from input
    base_name = os.path.splitext(os.path.basename(args.input_csv))[0]

    # Create confusion matrix plot
    cm_output_path = os.path.join(args.output_dir, f"{base_name}_heatmap.png")
    create_confusion_matrix_plot(cm_df, cm_output_path, args.title)

    # Create metrics plot
    metrics_output_path = os.path.join(args.output_dir, f"{base_name}_metrics.png")
    create_metrics_plot(metrics, metrics_output_path)

    # Print summary
    print("\nClassification Summary:")
    print("-" * 50)
    for i, class_name in enumerate(metrics["class_names"]):
        print(
            f"{class_name:12} | Precision: {metrics['precision'][i]:.3f} | "
            f"Recall: {metrics['recall'][i]:.3f} | F1-Score: {metrics['f1_score'][i]:.3f}"
        )
    print("-" * 50)
    print(f"Overall Accuracy: {metrics['accuracy']:.3f}")

    # Show plots
    plt.show()


if __name__ == "__main__":
    main()
