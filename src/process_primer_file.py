import pandas as pd
from Bio import SeqIO

def load_reference_sequence(fasta_file):
    """
    Load the reference sequence from a FASTA file and return both the forward and reverse strands.
    
    Args:
        fasta_file (str): Path to the FASTA file.

    Returns:
        tuple: (forward_strand, reverse_strand)
    """
    for record in SeqIO.parse(fasta_file, "fasta"):
        forward_strand = str(record.seq)
        reverse_strand = str(record.seq.reverse_complement())
        return forward_strand, reverse_strand

# File paths
bed_file_path = "/home/mahmadi/mpox_seqs/primer_output.bed"
fasta_file_path = "/home/mahmadi/mpox_seqs/reference.fasta"
output_file_path = "/home/mahmadi/mpox_seqs/primer.bed"

# Column names
col_names = ["ref", "start", "end", "name"]

# Read the primer bed file
primer_bed = pd.read_csv(bed_file_path, sep="\t", names=col_names)

# Process the DataFrame
primer_bed["pool"] = 1
primer_bed['handedness'] = primer_bed['name'].apply(lambda x: '+' if 'F' in x else '-')
primer_bed['name'] = "value_" + primer_bed['name']
primer_bed['name'] = primer_bed['name'].apply(lambda x: x.replace('_F', '_LEFT'))
primer_bed['name'] = primer_bed['name'].apply(lambda x: x.replace('_R', '_RIGHT'))

# Load reference sequences
reference_sequence, rev_ref = load_reference_sequence(fasta_file_path)

# Function to get sequence from the correct strand
def get_sequence(row, forward_seq, reverse_seq):
    try:
        if row['handedness'] == '+':
            # Extract from forward strand
            return forward_seq[row['start']-1:row['end']]
        else:
            # Extract from reverse strand
            return reverse_seq[row['end']-1:row['start']]
    except IndexError:
        # Handle out-of-bounds cases
        return "Out of bounds"

# Apply the function to the DataFrame
primer_bed['sequence'] = primer_bed.apply(lambda row: get_sequence(row, reference_sequence, rev_ref), axis=1)

# Debugging: Print rows with issues
print(primer_bed[primer_bed['sequence'] == "Out of bounds"])

# Save the DataFrame to a file
primer_bed.to_csv(output_file_path, sep="\t", index=False, header=False)
