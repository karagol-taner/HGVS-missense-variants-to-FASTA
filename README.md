# HGVS missense variants to FASTA

This project is a Python script that allows researchers to generate variant protein and nucleotide sequences in FASTA format based on wild-type protein/nucleotide sequences and variant notations in HGVS format. This tool is particularly useful for researchers and bioinformaticians working with protein sequences to study the effects of missense variants.

### Installation

- Option 1 (RECOMMENDED, for all platforms): You can install the package via pip (if you have Python 3.x on your system), and use the command HGSVtoFASTA directly:
```
pip install HGVStoFASTA      
```
```
HGVStoFASTA
```


- Option 2 (for all platforms): You can run the script directly from Google Colab: [HGVS-to-FASTA.ipynb on Google Colab](https://colab.research.google.com/drive/1yiqgo0joTMsBdOz57pTI6i5LgBDWO3zw?usp=sharing).
  
- Option 3 (for Windows): You can directly run this .EXE file: [Download .EXE file](https://drive.google.com/file/d/1rrDwS52b_H1F8sdb93SGLpxQ4JhbRyBB/view?usp=sharing).

Latest version on PyPI:

[![PyPI version](https://badge.fury.io/py/HGVStoFASTA.svg)](https://badge.fury.io/py/HGVStoFASTA)

### Requirements
If you want to run the script natively on your local computer (Option 1), ensure you have Python 3.x installed on your system. 

If Biopython is not already installed, the script will prompt you to install it automatically. 
If the script fails to automatically install Biopython, you can install it manually using the following steps:

- Open a terminal or command prompt after installing Python.
- Run the following command to install Biopython using pip:
  
   ```
   pip install biopython
   ```

### Usage/Features
1. Cross-Platform: Supports Windows, Mac, and Linux. Compatible with High Performance Computing (HPC).

2. Quick Start: Follow the prompts to input the required information, including the protein/gene name, wild-type sequence, variant notations, and the desired filename for the output FASTA file. 


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
