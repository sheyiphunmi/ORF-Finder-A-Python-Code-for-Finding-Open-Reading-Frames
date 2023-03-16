# ORF Finder: A Python Code for Finding Open Reading Frames

This is a Python script that finds all the open reading frames (ORFs) in a given genome sequence, and writes the results to a file. It also writes the top N (default: 20) ORFs by length for each start codon into the same file. The script takes in a genetic code file, a genome file, and the minimum and maximum length threshold as input parameters.

The script contains several functions that perform different tasks:
* **read_genetic_code**: Reads a file containing the genetic code and returns the start and stop codons.
*	**read_genome**: Reads a file containing the genome sequence and returns the genome.
*	**reverse_complement**: Returns the reverse complement of a DNA sequence.
*	**find_sort_orfs**: Finds all the open reading frames (ORFs) in a genome sequence and writes the results to a file. Also writes the top N orfs by length for each start codon into the same file.
*	**top_n_orfs**: Sorts the ORFs by length and writes the top N (default: 13,15,18,21) ORFs by length to a text file. Takes in a list of dictionaries containing the ORFs for each start codon, the name of the output file, and the number of top ORFs to write to the output file.


# Usage

1. Clone the repository:

```
git clone https://github.com/sheyiphunmi/ORF-Finder-A-Python-Code-for-Finding-Open-Reading-Frames.git
```

2. Run the Scripts
 
```
python orf_finder.py <genome_file_name> <genetic_code_file_name> <min_len_threshold> <max_len_threshold>
```

where:

* <genome_file_name>: The name of the FASTA file containing the genome sequence.
* <genetic_code_file_name>: The name of the text file containing the genetic code table. The first line of the file is expected to have the start codons seperated by space and the third line is expected to have the stop codons seperated by space
* <min_len_threshold>: The minimum ORF length threshold (optional).
* <max_len_threshold>: The maximum ORF length threshold (optional).

# Outputs

* **orf_all.txt**: This file contains information about all the open reading frames (ORFs) found in the genome for each start codon and their corresponding top 20 ORFs sorted based on length. The file includes the reading frame, start and stop positions, length, and sequence for each ORF.
* **top_N_orf_all.txt**: This file contains the top N ORFs by length after combining the ORFs found for each start codon. By default, this will generate four text files: top_13_orfs.txt, top_15_orfs.txt, top_18_orfs.txt, and top_21_orfs.txt. These files contain the top N ORFs across all start codons, sorted by length.

# Contributions

Contributions to the ORF Finder are welcome. Please create a pull request or open an issue to suggest new features, report bugs, or ask questions.

# License

This project is licensed under the [MIT License](https://opensource.org/license/mit/)
