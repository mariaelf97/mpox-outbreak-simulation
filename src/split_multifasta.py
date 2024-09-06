from Bio import SeqIO
import os

def split_multifasta_to_individual_files(multifasta_path, output_dir):
    """
    Splits a multi-record FASTA file into individual FASTA files.
    
    Args:
        multifasta_path (str): Path to the multi-record FASTA file.
        output_dir (str): Directory where individual FASTA files will be saved.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Parse the multi-record FASTA file
    with open(multifasta_path, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # Generate the output file path
            output_file_path = os.path.join(output_dir, f"{record.id}.fasta")
            
            # Write the individual record to a new FASTA file
            with open(output_file_path, "w") as output_handle:
                SeqIO.write(record, output_handle, "fasta")
            print(f"Record {record.id} written to {output_file_path}")

if __name__ == "__main__":
    # Replace these paths with your actual file and directory paths
    multifasta_file_path = "/home/mahmadi/mpox_seqs/sequence.fasta"
    output_directory = "/home/mahmadi/mpox_seqs/sequences/"
    
    split_multifasta_to_individual_files(multifasta_file_path, output_directory)
