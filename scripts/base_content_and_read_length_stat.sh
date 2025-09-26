import os
import gzip
import pandas as pd
from Bio import SeqIO
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import traceback

def parse_fastq(file_path, n_bases, reverse=False):
    """
    Parses a FASTQ file and calculates base frequencies for the first n_bases or last n_bases (if reverse=True).
    """
    base_counts = {base: [0] * n_bases for base in 'ATCGN'}
    total_counts = [0] * n_bases
    sequence_lengths = []

    # Detect if the file is gzipped
    if file_path.endswith('.gz'):
        handle = gzip.open(file_path, 'rt')  # Open as text
    else:
        handle = open(file_path, 'r')

    with handle:
        sequence_found = False
        for record in SeqIO.parse(handle, 'fastq'):
            sequence_found = True
            sequence = str(record.seq)
            seq_len = len(sequence)
            sequence_lengths.append(seq_len)

            if not sequence.strip('N'):  
                continue

            if reverse:
                # Reverse indexing: count from the last base
                for i in range(1, n_bases + 1):
                    if i > seq_len:  # Skip if the sequence is shorter than the index
                        break
                    base = sequence[-i]
                    if base in base_counts:
                        base_counts[base][i - 1] += 1
                        total_counts[i - 1] += 1
            else:
                # Forward indexing: count from the first base
                for i in range(min(n_bases, seq_len)):  
                    base = sequence[i]
                    if base in base_counts:  
                        base_counts[base][i] += 1
                        total_counts[i] += 1

        if not sequence_found:
            raise ValueError("The FASTQ file is empty or contains no valid sequences.")

    if all(count == 0 for count in total_counts):
        raise ValueError("Total counts for all positions are zero. Check the input FASTQ file for valid sequences.")

    # Calculate frequencies
    base_freqs = {
        base: [count / float(total) if total > 0 else 0 for count, total in zip(counts, total_counts)]
        for base, counts in base_counts.items()
    }

    return pd.DataFrame(base_freqs), sequence_lengths


def plot_with_length_histogram(forward_df, reverse_df, lengths, output_plot, file_name, n_bases):
    """
    Plots forward and reverse base frequencies in separate rows and adds a length frequency histogram with updated style.
    """
    # Ensure lengths are integers
    lengths = [int(length) for length in lengths]
    length_counts = pd.Series(lengths).value_counts().sort_index()

    if length_counts.empty:
        raise ValueError("No valid sequence lengths found.")

    fig = make_subplots(
        rows=3, cols=1,
        row_heights=[0.4, 0.4, 0.2],
        shared_xaxes=False,
        vertical_spacing=0.15,
        subplot_titles=[
            "Base Frequencies (Forward)",
            "Base Frequencies (Reverse)",
            "Sequence Length Distribution"
        ],
        specs=[[{"secondary_y": False}], [{"secondary_y": False}], [{"secondary_y": False}]]
    )

    colors = ['blue', 'orange', 'green', 'red', 'purple']

    # 1. Forward base frequencies
    for base, color in zip(forward_df.columns, colors):
        fig.add_trace(
            go.Scatter(
                x=forward_df.index + 1,
                y=forward_df[base],
                mode='lines+markers',
                name=f'{base} (Forward)',
                line=dict(color=color),
                legendgroup="forward"
            ),
            row=1, col=1
        )

    # 2. Reverse base frequencies
    for base, color in zip(reverse_df.columns, colors):
        fig.add_trace(
            go.Scatter(
                x=-(reverse_df.index + 1),
                y=reverse_df[base],
                mode='lines+markers',
                name=f'{base} (Reverse)',
                line=dict(color=color, dash='dot'),
                legendgroup="reverse"
            ),
            row=2, col=1
        )

    # 3. Length frequency histogram
    fig.add_trace(
        go.Bar(
            x=length_counts.index,
            y=length_counts.values,
            name="Sequence Length",
            marker=dict( color="red"),
        ),
        row=3, col=1
    )

    # Add gray background rectangle
    fig.add_shape(
        type="rect",
        x0=min(length_counts.index) - 0.5,
        x1=max(length_counts.index) + 0.5,
        y0=0,
        y1=max(length_counts.values),
        fillcolor="lightgray",
        opacity=0.2,
        layer="below"
    )

    # Update axes
    fig.update_xaxes(title_text="Base Position (Forward: 1+)", row=1, col=1)
    fig.update_xaxes(title_text="Base Position (Reverse: -1-)", row=2, col=1)
    fig.update_xaxes(
        title_text="Sequence Length",
        tickmode='linear',
        dtick=1,
        row=3, col=1
    )
    fig.update_yaxes(title_text="Frequency", range=[0, 1], row=1, col=1)
    fig.update_yaxes(title_text="Frequency", range=[0, 1], row=2, col=1)
    fig.update_yaxes(
        title_text="Frequency",
        range=[0, max(length_counts.values) * 1.1],  # 修正 Y 轴范围
        row=3, col=1
    )

    # Update layout
    fig.update_layout(
        title=f'Base Frequencies and Sequence Length Distribution in {file_name}',
        height=900,
        template="plotly_white",
        legend_title="Base",
        font=dict(size=12)
    )

    # Save plot
    fig.write_html(output_plot)


def main(input_file,  n_bases, out_dir):
    try:
        # Determine file name and output paths
        file_name = os.path.splitext(os.path.basename(input_file))[0]
        output_html = os.path.join(out_dir, "{}_frequencies.html".format(file_name))

        print(f"Input file: {input_file}")
        print(f"Output HTML: {output_html}")

        # Parse forward and reverse base frequencies and collect sequence lengths
        forward_frequencies, lengths = parse_fastq(input_file, n_bases, reverse=False)
        reverse_frequencies, _ = parse_fastq(input_file, n_bases, reverse=True)

        print("Parsed FASTQ data successfully.")
        print(f"Forward frequencies: {forward_frequencies.shape}")
        print(f"Reverse frequencies: {reverse_frequencies.shape}")
        print(f"Lengths: {len(lengths)} sequences")

        # Generate interactive plot

        plot_with_length_histogram(forward_frequencies, reverse_frequencies, lengths, output_html, file_name, n_bases)
        print(f"Plot generated successfully and saved to {output_html}")

    except Exception as e:
        print(f"An error occurred: {e}")
        traceback.print_exc()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Calculate and plot base frequencies in a FASTQ file.')
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to the FASTQ file (supports .gz compression).')
    parser.add_argument('-n', '--num_bases', type=int, default=85, help='Number of bases to analyze (default: 85).')
    parser.add_argument('-o', '--out_dir', type=str, required=True, help='Generate an interactive HTML plot in out_dirt.')

    args = parser.parse_args()

    main(args.input_file,  args.num_bases, args.out_dir)