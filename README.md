# HGVS missense variants to FASTA

This project is a Python script that allows researchers to generate variant protein sequences in FASTA format based on wild-type protein sequences and variant notations in HGVS format. This tool is particularly useful for researchers and bioinformaticians working with protein sequences to study the effects of missense variants.


### Requirements
Ensure you have Python 3.x installed on your system. 
If Biopython is not already installed, the script will prompt you to install it automatically. 
If the script fails to automatically install Biopython, you can install it manually using the following steps:

- Open a terminal or command prompt after installing Python.
- Run the following command to install Biopython using pip:
  
   ```
   pip install biopython
   ```

### Usage
1. Running the Script: Execute the HGVS_to_FASTA.py script using Python. Follow the prompts to input the required information, including the protein/gene name, wild-type sequence, variant notations, and the desired filename for the output FASTA file.
   - ! NEW UPDATE Input types for variants
   - 1-Manual Input
   - 2-File Input: Users can input variant notations from a text file, making it convenient to handle large lists of variants.
     If you choose to input variant notations from a text file, the file should contain the variant notations in the following format:
     
      ```
      c.92A>T
      c.150G>A
      ...
      ```


3. Review and Save: Verify the entered information before proceeding. Once confirmed, the script will generate the variant protein sequences and save them to the specified file in FASTA format.

4. Viewing Results: After execution, the script will automatically open the generated FASTA file using the default system application for viewing.


### Questions
If you have any questions or encounter issues, don't hesitate to reach out.

### License
This project is licensed under the  GPL-3.0 License - see the LICENSE file for details.
