import sys
import os
import logging

# # define the start codons 
# start_codons_stdcode  = "ATG"
# # define the stop codons
# stop_codons_stdcode  = "TGA TAA TAG"

len_codon = 3
# set the threshold for ORF length
#min_len_threshold, max_len_threshold = 200, 1850

# Set the number of orfs to be writing to file after sorting by length.
# a list of different. 13 means write the top 13 ORFS from the combined ORFs for all start codons.
top_n_orf = [13,15,18,21]  # Must be a list; length do not matter.

logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

def read_genetic_code(genetic_code_file_name):
    """
    Reads the genetic code file and returns the start and stop codons

    Args:
    - genetic_code_file_name (str): the name of the file containing the genetic code

    Returns:
    - starts (list of str): the start codons
    - stops (list of str): the stop codons
    """
    with open(genetic_code_file_name) as file:
        file.readline()
        starts = file.readline().strip().split()
        stops = file.readline().strip().split()
    return starts, stops


def read_genome(genome_file_name):
    """
    Reads the genome file and returns the genome

    Args:
    - genome_file_name (str): the name of the file containing the genome

    Returns:
    - genome (str): the genome
    """
    with open(genome_file_name) as file:
        file.readline() #next(file)
        genome = ''.join(file.read().strip().split())
    return genome

def reverse_complement(genome):
    """
    Returns the reverse complement of a DNA sequence

    Args:
    - genome (str): the DNA sequence

    Returns:
    - rev_seq (str): the reverse complement of the DNA sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
    # Reverse the sequence and complement each nucleotide
    rev_seq = ''.join([complement[nuc] for nuc in genome[::-1]])  # genome[::-1]] reverses the order of the nucleotides in the input sequence genome
    # rev_seq = []
    # # for nuc in genome[::-1]:
    # #     "".join(rev_seq.append(complement[nuc]))
    return rev_seq


def find_sort_orfs(genome, starts, stops, output_file_name, min_len_threshold = None, max_len_threshold = None):
    """
    Finds all the open reading frames (ORFs) in a genome sequence and writes the results to a file. 
    Also write the top N orfs by length for each start codon into the same file


    Args:
    - genome (str): the genome sequence
    - starts (list of str): the start codons
    - stops (list of str): the stop codons
    - output_file_name (str): the name of the output file
    - min_len_threshold (int): Minimum length threshold of ORF (Default = None)
    - max_len_threshold (int): Maximum length threshold of ORF (Default = None)

    Returns:
    - orfs (list of dictionaries): orfs for all start codons and their information
    """
    orfs = []
    with open(output_file_name, 'w') as output_file:
        for start_codon in starts:
            start_codon_orfs = []
            logging.info(f"\nSTART CODON USED: {start_codon}")
            output_file.write(f"\n\nSTART CODON USED: {start_codon}\n")
            for reading_frame in range(3):
                logging.info(f"\nSearching for ORFs in reading frame {reading_frame + 1}")
                output_file.write(f"\n\nORFS IN READING FRAME {reading_frame + 1}\n")
                for i in range(reading_frame, len(genome) - 2, 3):
                    if genome[i:i+3] == start_codon:  # if a start codon is found at the current position
                        for j in range(i+3, len(genome) - 2, 3): # starting at the next codon, iterate over the genome in increments of three
                            if genome[j:j+3] in stops: # if a stop codon is found at the current position
                                len_orf = j - i + 3
                                start_orf, stop_orf = i, j
                                if len_orf % 3 == 0: # Checks if ORF is a potential coding sequence
                                    if len_orf >= min_len_threshold and len_orf <= max_len_threshold:
                                        # Check if the ORF is already included in a longer ORF
                                        contained = False
                                        for prev_orf in orfs:
                                            if prev_orf["reading_frame"] == reading_frame+1 and prev_orf["start"] <= start_orf and prev_orf["stop"] >= stop_orf:
                                                contained = True
                                                break
                                        if not contained:
                                            orf = genome[start_orf:stop_orf + 3]
                                            # Display information about the ORF that was found
                                            logging.info(f"ORF found at {start_orf} {stop_orf}")
                                            output_file.write(f"\nORF found at start: {start_orf} stop: {stop_orf}. The length is {len_orf}\n")
                                            orf_dict = {"reading_frame": reading_frame+1,"start": i+1, "stop": stop_orf+3, "length": len_orf, "sequence": orf}
                                            output_file.write(f"{orf}\n")
                                            orfs.append(orf_dict)
                                            start_codon_orfs.append(orf_dict)  ## save information about each orf in a list
                                    break

            # Prints top 20 ORFS for each start codon based on length
            output_file.write(f"\n\nPRINTING THE RESULT FOR THE TOP 20 ORFS FOR {start_codon}\n")
            start_codon_orfs.sort(key=lambda x: x["length"], reverse=True)
            for i, gene in enumerate(start_codon_orfs[:20]):
                output_file.write(f"\nORF found in reading frame {gene['reading_frame']} at start: {gene['start']} stop: {gene['stop']}. The length is {gene['length']}\n")
                output_file.write(f"{gene['sequence']}\n")
    return orfs


def top_n_orfs(orfs_list, output_file_name, top_n_orf):
    """
    Sorts the ORFs by length and writes the top N ORFs by length to a text file.

    Args:
    - orfs_list (list): A list of dictionaries containing the ORFs for each start codon.
    - output_file_name (str): The name of the output file.
    - top_n_orf (int): The number of top ORFs to write to the output file.

    Returns:
    - None
    """
    if not orfs_list:
        raise ValueError("Empty list of ORFs.")
    if top_n_orf > len(orfs_list):
        raise ValueError("Number of top ORFs is greater than the total number of ORFs.")

    with open(output_file_name, 'w') as output_file:
        # Sort the list of all ORFs found in the genome by length, in descending order
        orfs_list.sort(key=lambda x: x["length"], reverse=True)

        # Write the top N ORFs to the output file
        output_file.write(f"\n\nPRINTING THE RESULT FOR THE TOP {top_n_orf} ORFs\n")
        for i, gene in enumerate(orfs_list[:top_n_orf]):
            output_file.write(f"\nORF found in reading frame {gene['reading_frame']} at start: {gene['start']} stop: {gene['stop']}. The length is {gene['length']}\n")
            output_file.write(f"{gene['sequence']}\n")


def main():
    genome_file_name = sys.argv[1]
    genetic_code_file_name = sys.argv[2]
    if len(sys.argv) == 3:
        print("Using default minimum and maximum length thresholds: 200 and 1850, respectively")
        min_len_threshold = 200
        max_len_threshold = 1850
    elif len(sys.argv) == 4:
        min_len_threshold = int(sys.argv[3])
        max_len_threshold = 1850
        print("Using the minimum length threshold as {} and the maximum as the default (1850).".format(sys.argv[3]))
    elif len(sys.argv) == 5:
        min_len_threshold = int(sys.argv[3])
        max_len_threshold = int(sys.argv[4])
        print("Using the minimum length threshold as {} and the maximum as {}.".format(sys.argv[3], sys.argv[4]))
    else:
        print("Invalid number of arguments")
        return

    # Read the genetic code
    starts, stops = read_genetic_code(genetic_code_file_name)

    # Read the genome
    genome = read_genome(genome_file_name)

    # Get the reverse complement of the genome sequence, then
    # modify the genome sequence in place by concatenating the reverse complement to it
    # to form a double strands
    dsgenome = genome + reverse_complement(genome)

    # Find ORFs for each start codon in all reading frames
    orfs = find_sort_orfs(dsgenome, starts, stops, "orf_all.txt", min_len_threshold, max_len_threshold)

    # Print top N for all start codons
    for num in top_n_orf:
        top_n_orfs(orfs,f"top_{num}_orf_all.txt", num)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python find_orfs.py genome_file genetic_code_file [min_length_threshold] [max_length_threshold]")
        sys.exit()
    main()
    logging.info("Done searching for ORFs")


   
