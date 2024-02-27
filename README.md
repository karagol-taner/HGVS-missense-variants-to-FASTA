# HGVS missense variants to FASTA

This project is a Python script that allows researchers to generate variant protein and nucleotide sequences in FASTA format based on wild-type protein/nucleotide sequences and variant notations in HGVS format. This tool is particularly useful for researchers and bioinformaticians working with protein sequences to study the effects of missense variants.


### Requirements
If you want to run the script natively on your local computer, ensure you have Python 3.x installed on your system. 
If Biopython is not already installed, the script will prompt you to install it automatically. 
If the script fails to automatically install Biopython, you can install it manually using the following steps:

- Open a terminal or command prompt after installing Python.
- Run the following command to install Biopython using pip:
  
   ```
   pip install biopython
   ```

### Usage
1. Follow the prompts to input the required information, including the protein/gene name, wild-type sequence, variant notations, and the desired filename for the output FASTA file. For running the Python script:
    - Option 1: Install the package via pip, and use the command HGSVtoFASTA directly:
      ```
      pip install HGSVtoFASTA
      HGSVtoFASTA
      ```
      
    - Option 2: You can run the script directly from Google Colab: [HGVS-to-FASTA.ipynb on Google Colab](https://colab.research.google.com/drive/1yiqgo0joTMsBdOz57pTI6i5LgBDWO3zw?usp=sharing).
    
    - Option 3: Execute the HGVS_to_FASTA.py script using Python. 
    For runing the script you can download the [.py file](HGVS_to_FASTA.py) and open a terminal/command prompt. Navigate to the directory where the .py file is located, and then execute the script using the following command:
      ```
      python HGVS_to_FASTA.py
      ```
  
3. Input types for variants
   - 1-Manual Input
   - 2-File Input: Users can input variant notations from a text file, making it convenient to handle large lists of variants.
     If you choose to input variant notations from a text file, the file should contain the variant notations in the following format for nucleotide:
     
      ```
      c.92A>T
      c.150G>A
      ...
      ```
      
      the following format for aminoacid:
       ```
      p.Arg107Cys
      p.Arg112Cys
      ...
      ```

4. Review and Save: Verify the entered information before proceeding. Once confirmed, the script will generate the variant protein/nucleotide sequences and save them to the specified file in FASTA format.

5. Viewing Results: After execution, the script will automatically open the generated FASTA file using the default system application for viewing.


### Questions
If you have any questions or encounter issues, don't hesitate to reach out.

### License
This project is licensed under the  GPL-3.0 License - see the LICENSE file for details.
