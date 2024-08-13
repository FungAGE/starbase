import subprocess
import pandas as pd
import plotly.io as pio
import plotly.graph_objects as go


def run_lastz(sequence1, sequence2, output_file):
    """
    Runs LASTZ to align two sequences and writes the output to a specified file.
    """
    command = f"lastz {sequence1} {sequence2} --format=general --output={output_file}"
    subprocess.run(command, shell=True, check=True)


def parse_lastz_output(output_file):
    """
    Parses the LASTZ output file and returns a DataFrame with relevant data.
    """
    columns = [
        "#score",
        "name1",
        "strand1",
        "size1",
        "zstart1",
        "end1",
        "name2",
        "strand2",
        "size2",
        "zstart2",
        "end2",
        "identity",
        "idPct",
        "coverage",
        "covPct",
    ]
    df = pd.read_csv(output_file, sep="\t", header=0, comment=None)
    return df


def create_alignment_plot(df, output_image_file):
    """
    Creates a Plotly scatter plot from the alignment DataFrame and saves it as a PNG.
    """
    # Extract individual columns and prepare the data for plotting
    x_values = []
    y_values = []

    for _, row in df.iterrows():
        x_values.append(row["zstart1"])
        y_values.append(row["zstart2"])
        x_values.append(row["end1"])
        y_values.append(row["end2"])

    print(x_values)
    print(y_values)
    # Create the scatter plot
    fig = go.Figure(
        data=go.Scatter(
            x=x_values,
            y=y_values,
        )
    )

    # Add layout details
    fig.update_layout(
        title="LASTZ Alignment",
        xaxis_title="Sequence Start",
        yaxis_title="Sequence End",
    )


# Example usage:
sequence1 = "tmp/seq1.fa"
sequence2 = "tmp/seq2.fa"
output_file = "tmp/lastz_output.txt"

# Run LASTZ alignment
run_lastz(sequence1, sequence2, output_file)

# Parse LASTZ output
alignment_df = parse_lastz_output(output_file)

# Create Plotly plot and save as PNG
create_alignment_plot(alignment_df, output_image_file)
