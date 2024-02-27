import subprocess
import importlib
import os
import re
import sys

def run():
    print("HGVStoFASTA: HGVS missense variants to FASTA format converter v0.4.2")
    print(" ")
    print(" ")

if __name__ == "__main__":
    main()

# Check if Biopython is already installed
try:
    importlib.import_module("Bio")
    biopython_installed = True
    print("Biopython is already installed.")
    input("To continue, press any button...")
except ImportError:
    biopython_installed = False

# If Biopython is not installed, download and install it using pip
if not biopython_installed:
    user_input = input("Biopython is not installed. Do you want to install it? (Y/N): ").strip().lower()
    if user_input == 'y':
        print("Installing Biopython...")
        try:
            subprocess.run(["pip", "install", "biopython"], check=True)
            print("Biopython installation complete.")
            input("To continue, press any button...")
            # Check if Biopython is successfully installed
            try:
                importlib.import_module("Bio")
                biopython_installed = True
            except ImportError:
                biopython_installed = False
        except subprocess.CalledProcessError as e:
            print(f"Error installing Biopython: {e}")
    else:
        print("Biopython installation skipped.")
        input("To continue, press any button...")

# If Biopython is installed or successfully installed, continue with processing fasta sequences
if biopython_installed or importlib.import_module("Bio"):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO


# Protein name
protein_name = input("Enter the protein/gene name for FASTA header line: >").strip().upper()

# Prompt the user to choose input method
input_method = input("Choose input method for sequences (1 - manual input, 2 - input from file): ").strip()

# Initialize wild-type nucleotide sequence
wild_type_sequence = ""

# Input from file
if input_method == '2':
    file_path = input("Enter the path to the file containing the nucleotide / amino acid sequence (eg. C:\\wild_type.fasta): ").strip()
    try:
        with open(file_path, "r") as file:
            # Parse the FASTA file to extract the nucleotide sequence
            for record in SeqIO.parse(file, "fasta"):
                if not wild_type_sequence:
                    wild_type_sequence = str(record.seq).upper()
                else:
                    print("Error: Multiple sequences found in the input file. Please ensure the file contains only one sequence.")
                    exit()
    except FileNotFoundError:
        print("File not found. Please make sure the file path is correct.")
        exit()

# Manual input
elif input_method == '1':
    wild_type_sequence = input("Enter the wild-type nucleotide sequence of the protein: ").strip().upper()

else:
    print("Invalid input method. Please choose either 1 or 2.")
    exit()
    
# Determine if the sequence is nucleotide or amino acid
sequence_type = ""
valid_nucleotides = set("ACGTU")
valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")

if all(char in valid_nucleotides for char in wild_type_sequence):
    sequence_type = "nucleotide"
elif all(char in valid_amino_acids for char in wild_type_sequence):
    sequence_type = "amino acid"
else:
    print("Invalid sequence format. The sequence must consist of valid nucleotides (ACGTU) or amino acids (ACDEFGHIKLMNPQRSTVWY).")
    exit()

if sequence_type == "nucleotide":
# Prompt the user to choose input method
    input_method_var = input("Choose input method for missense variants in HGVS format (1 - manual input, 2 - input from file): ").strip()

# Initialize variant_notations list
    variant_notations = []

# Input from file
    if input_method_var == '2':
        print(" ")
        print("The text file must contain variant notations in the following format: ")
        print("c.92A>T")
        print("c.150G>A")
        print("...")    
        print(" ") 
        file_path = input("Enter the path to the file containing variant notations (eg.C:\\variantHGVS.txt): ").strip()
        try:
            with open(file_path, "r") as file:
                variant_notations = [line.strip() for line in file.readlines()]
        except FileNotFoundError:
            print("File not found. Please make sure the file path is correct.")
            exit()

# Manual input
    elif input_method_var == '1':
        variant_input = input("Enter all missense variant notations in HGVS, separated by commas (eg. c.107A>C, c.112A>C, ...): ").strip()
        variant_notations = [notation.strip() for notation in variant_input.split(",")]

    else:
        print("Invalid input method. Please choose either 1 or 2.")
        exit()
        
if sequence_type == "amino acid":
    # Prompt the user to choose input method for amino acid variants
    input_method_var = input("Choose input method for missense variants in HGVS format (1 - manual input, 2 - input from file): ").strip()

    # Initialize variant_notations list
    variant_notations = []

    # Input from file
    if input_method_var == '2':
        print(" ")
        print("The text file must contain variant notations in the following format: ")
        print("p.Arg107Cys")
        print("p.Arg112Cys")
        print("...")    
        print(" ") 
        file_path = input("Enter the path to the file containing variant notations (e.g., C:\\variantHGVS.txt): ").strip()
        try:
            with open(file_path, "r") as file:
                variant_notations = [line.strip() for line in file.readlines()]
        except FileNotFoundError:
            print("File not found. Please make sure the file path is correct.")
            exit()

    # Manual input
    elif input_method_var == '1':
        variant_input = input("Enter all missense variant notations in HGVS, separated by commas (eg. p.Arg107Cys, p.Glu112Cys, ...): ").strip()
        variant_notations = [notation.strip() for notation in variant_input.split(",")]

    else:
        print("Invalid input method. Please choose either 1 or 2.")
        exit()


# Print the wild-type sequence and the list of variants entered by the user for verification
print("\nEntered Sequence and Variants:")
print("")
print("Protein/Gene Name: >", protein_name)
print("")
print("Wild-type nucleotide/aminoacid sequence:", wild_type_sequence)
print("")
print("Variants:")
for variant in variant_notations:
    print(variant)
print("")
print("")
print("\nWarning: Please check the information thoroughly before proceeding.")
input("To continue, press any button...")
print("")

# Initialize variant_sequences dictionary
variant_sequences = {}

# Iterate over each variant notation
for notation in variant_notations:
    # Parse position and variant from the notation
    if sequence_type == "nucleotide":
        position = int(notation.split('.')[1][:-3])  # Extract position (e.g., 92)
        ref_nucleotide = notation[-3]  # Reference nucleotide (e.g., A)
        var_nucleotide = notation[-1]  # Variant nucleotide (e.g., T)
        
        # Apply nucleotide change to wild-type sequence
        variant_sequence = wild_type_sequence[:position-1] + var_nucleotide + wild_type_sequence[position:]
    elif sequence_type == "amino acid":
        # For amino acid sequences in HGVS format like p.Arg3Tyr
        # Parse position and variant from the notation
        parts = notation.split('.')
        if len(parts) != 2:
            print("Invalid HGVS format for amino acid variants. Please use the format like p.Arg3Tyr.")
            exit()
        amino_acid_change = parts[1]
        ref_aa, position, var_aa = re.match(r'([A-Za-z]+)(\d+)([A-Za-z]+)', amino_acid_change).groups()
        position = int(position)
        
        # Map amino acid abbreviations to single-letter codes
        aa_mapping = {
            'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Glu': 'E', 'Gln': 'Q',
            'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F',
            'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', 'Ter': '*'
        }
        ref_aa_code = aa_mapping.get(ref_aa, '?')  # Use '?' if abbreviation is not found
        var_aa_code = aa_mapping.get(var_aa, '?')  # Use '?' if abbreviation is not found

        
        # Apply amino acid change to wild-type sequence
        variant_sequence = wild_type_sequence[:position-1] + var_aa_code + wild_type_sequence[position:]

    # Store variant sequence in dictionary
    variant_sequences[notation] = variant_sequence

# Print variant sequences
for notation, sequence in variant_sequences.items():
    print(f"{notation}: {sequence}")

# Get the output file name from user input
output_file_path = input("Enter the name for the FASTA file (e.g., variants.fasta): ").strip()

# Save variant sequences to the output file in FASTA format
with open(output_file_path, "w") as f:
    f.write(f">{protein_name}\n")
    f.write(wild_type_sequence + "\n")

    for i, (notation, sequence) in enumerate(variant_sequences.items(), start=1):
        f.write(f">{protein_name}_variant_{i}\n")
        f.write(sequence + "\n")

print(f"\nVariant sequences saved to {output_file_path}")

# Open the FASTA file using the default system application
try:
    if os.name == "posix":  # Linux or macOS
        subprocess.run(["xdg-open", output_file_path], check=True)
    elif os.name == "nt":  # Windows
        os.startfile(output_file_path)
    else:
        print(f"Unsupported operating system. Please navigate to {output_file_path} manually.")
except Exception as e:
    print(f"Error opening the output file: {e}")

input("\nPress Enter to exit...")