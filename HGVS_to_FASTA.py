import subprocess
import importlib
import os

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
protein_name = input("Enter the protein name for FASTA header line: >").strip().upper()

# Wild-type protein sequence of the protein
wild_type_sequence = input("Enter the wild-type protein sequence of the protein. (eg. MQRVNMIMAESPGLITICLLGYLLSAECTVFLDHENANKILNRPKRY...) = ").strip()

# Prompt the user to input all variant notations, separated by commas
variant_input = input("Enter all missense variant notations in HGVS, separated by commas (eg. c.107A>C, c.112A>C, ...) = ").strip()
variant_notations = [notation.strip() for notation in variant_input.split(",")]

# Dictionary to store variant protein sequences
variant_sequences = {}

# Print the wild-type sequence and the list of variants entered by the user for verification
print("\nEntered Sequence and Variants:")
print("")
print("Protein Name: >", protein_name)
print("")
print("Wild-type sequence:", wild_type_sequence)
print("")
print("Variants:")
for variant in variant_notations:
    print(variant)
print("")
print("")
print("\nWarning: Please check the information thoroughly before proceeding.")
input("To continue, press any button...")
print("")

# Iterate over each variant notation
for notation in variant_notations:
    # Parse nucleotide change from the variant notation
    position = int(notation.split('.')[1][:-3])  # Extract position (e.g., 92)
    ref_nucleotide = notation[-3]  # Reference nucleotide (e.g., A)
    var_nucleotide = notation[-1]  # Variant nucleotide (e.g., T)
    
    # Translate reference and variant nucleotides to amino acids
    ref_aa = Seq(ref_nucleotide).translate()
    var_aa = Seq(var_nucleotide).translate()

    # Apply amino acid change to wild-type sequence
    variant_sequence = wild_type_sequence[:position-1] + str(var_aa) + wild_type_sequence[position:]

    # Store variant protein sequence in dictionary
    variant_sequences[notation] = variant_sequence

# Print variant protein sequences
for notation, sequence in variant_sequences.items():
    print(f"{notation}: {sequence}")
    
# Get the log file name from user input
log_file_path = input("Enter the name for the FASTA file (eg. variants.fasta): ").strip()

# Save variant protein sequences to the log file in FASTA format
with open(log_file_path, "w") as f:
    f.write(f"> Wild-type_{protein_name}\n")
    f.write(wild_type_sequence + "\n")

    for i, (notation, sequence) in enumerate(variant_sequences.items(), start=1):
        f.write(f"> Variant_{i}_{protein_name}_{notation}\n")
        f.write(sequence + "\n")

print(f"\nVariant protein sequences saved to {log_file_path}")

# Open the FASTA file using the default system application
try:
    if os.name == "posix":  # Linux or macOS
        subprocess.run(["xdg-open", log_file_path], check=True)
    elif os.name == "nt":  # Windows
        os.startfile(log_file_path)
    else:
        print(f"Unsupported operating system. Please navigate to {log_file_path} manually.")
except Exception as e:
    print(f"Error opening the log file: {e}")

input("\nPress Enter to exit...")

