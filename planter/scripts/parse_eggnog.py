#!/opt/venv/bin/python
import argparse
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Define a COG category mapping dictionary
cog_mapping = {
    "J": "TRANSLATION, RIBOSOMAL STRUCTURE AND BIOGENESIS",
    "A": "RNA PROCESSING AND MODIFICATION",
    "K": "TRANSCRIPTION",
    "L": "REPLICATION, RECOMBINATION AND REPAIR",
    "B": "CHROMATIN STRUCTURE AND DYNAMICS",
    "D": "CELL CYCLE CONTROL, CELL DIVISION, CHROMOSOME PARTITIONING",
    "Y": "NUCLEAR STRUCTURE",
    "V": "DEFENSE MECHANISMS",
    "T": "SIGNAL TRANSDUCTION MECHANISMS",
    "M": "CELL WALL/MEMBRANE/ENVELOPE BIOGENESIS",
    "N": "CELL MOTILITY",
    "Z": "CYTOSKELETON",
    "W": "EXTRACELLULAR STRUCTURES",
    "U": "INTRACELLULAR TRAFFICKING, SECRETION, AND VESICULAR TRANSPORT",
    "O": "POSTTRANSLATIONAL MODIFICATION, PROTEIN TURNOVER, CHAPERONES",
    "X": "MOBILOME: PROPHAGES, TRANSPOSONS",
    "C": "ENERGY PRODUCTION AND CONVERSION",
    "G": "CARBOHYDRATE TRANSPORT AND METABOLISM",
    "E": "AMINO ACID TRANSPORT AND METABOLISM",
    "F": "NUCLEOTIDE TRANSPORT AND METABOLISM",
    "H": "COENZYME TRANSPORT AND METABOLISM",
    "I": "LIPID TRANSPORT AND METABOLISM",
    "P": "INORGANIC ION TRANSPORT AND METABOLISM",
    "Q": "SECONDARY METABOLITES BIOSYNTHESIS, TRANSPORT AND CATABOLISM",
    "R": "GENERAL FUNCTION PREDICTION ONLY",
    "S": "FUNCTION UNKNOWN",
}


# Function to expand COG categories
def expand_cog_categories(cog_string):
    if pd.isna(cog_string):
        return ["No COG assigned"]
    return [cog_mapping.get(cat, cat) for cat in cog_string]


def main():
    parser = argparse.ArgumentParser(
        description="Generate table of counts and a plot for eggNOG annotation file."
    )
    parser.add_argument("eggnog_file", help="Path to eggNOG annotation file")
    parser.add_argument(
        "output_dir", help="Directory to save the outputs (counts table and plot)"
    )

    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Read the eggNOG annotation file
    df = pd.read_csv(args.eggnog_file, sep="\t", comment="#", header=None)

    # Assign column names (update if necessary)
    df.columns = [
        "query",
        "seed_ortholog",
        "evalue",
        "score",
        "eggNOG_OGs",
        "max_annot_lvl",
        "COG_category",
        "Description",
        "Preferred_name",
        "GOs",
        "EC",
        "KEGG_ko",
        "KEGG_Pathway",
        "KEGG_Module",
        "KEGG_Reaction",
        "KEGG_rclass",
        "BRITE",
        "KEGG_TC",
        "CAZy",
        "BiGG_Reaction",
        "PFAMs",
    ]

    # Expand COG categories
    expanded_cogs = df["COG_category"].apply(expand_cog_categories)

    # Count the occurrences of each COG category
    cog_counts = expanded_cogs.explode().value_counts()

    # Save the counts to a CSV file
    counts_file = os.path.join(args.output_dir, "cog_category_counts.csv")
    cog_counts.to_csv(counts_file, header=["Count"], index_label="COG Category")

    # Create a bar plot with improved aesthetics
    plt.figure(figsize=(20, 12))
    ax = sns.barplot(x=cog_counts.index, y=cog_counts.values, palette="viridis")

    # Rotate x-axis labels and adjust their alignment
    plt.xticks(rotation=90, ha="right")

    # Add value labels on top of each bar
    for i, v in enumerate(cog_counts.values):
        ax.text(i, v + 0.5, str(v), ha="center", va="bottom")

    # Customize the plot
    plt.xlabel("COG Category", fontsize=14, fontweight="bold")
    plt.ylabel("Count", fontsize=14, fontweight="bold")
    plt.title("Distribution of COG Categories", fontsize=18, fontweight="bold")

    # Remove top and right spines
    sns.despine()

    # Adjust layout and save the figure
    plot_file = os.path.join(args.output_dir, "cog_distribution_plot.png")
    plt.tight_layout()
    plt.savefig(plot_file, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"COG category distribution table saved as '{counts_file}'")
    print(f"COG category distribution plot saved as '{plot_file}'")


if __name__ == "__main__":
    main()
