### IMPORT LIBRARIES ###


import requests, sys, csv, os
#from Bio.Blast import NCBIWWW as NCBI


### FUNCTION DEFINITIONS ###


def find_seq(id_csv):
    # Convert transcript IDs in a CSV file into gene/transcript sequences using
    # Ensembl web addresses

    # Read ID from csv
    seq_list=list(csv.reader(open(id_csv,encoding = "utf-8-sig")))

    # Formulate Ensembl web addresses for transcripts and genes
    addr = "https://rest.ensembl.org/sequence/id/"
    header_gene="?content-type=application/json"
    header_transcript="?content-type=application/json;type=cds"

    # Search for transcript and gene sequences on Ensembl
    for  id_transcript in seq_list:
        print('Searching for sequences : '+id_transcript[0])
        row = seq_list.index(id_transcript)
        req_gene = requests.get(addr+id_transcript[0]+header_gene)
        req_transcript = requests.get(addr+id_transcript[0]+header_transcript)

        # Annotate list for unidentifiable IDs
        if (req_gene.ok == 0) or (req_transcript.ok == 0):
            seq_list[row][1:3]=["Invalid","NA","NA"]
            print('Transcript or gene sequence NOT found.')
            continue

        # Annotate list for discovered IDs
        seq_gene = req_gene.json()['seq']
        seq_transcript = req_transcript.json()['seq']
        seq_list[row][1:3] = ["Valid",seq_gene,seq_transcript]
        print('Transcript and gene sequence found.')

    print('Sequence search complete.')
    return seq_list

def find_all(a_str, sub):
    # Find a subsequence (like a PAM site) in a sequence string

    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += 1

def reverse_complement(dna):
    # Generate a reverse complement of a sequence string to search the opposite
    # strand

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])

def build_guide(sequence,PAM_locations,orientation):
    # Take identified PAM locations in a sequence and build a guide from it

    # NEED TO REPLACE PAMS FAR FROM START SITE WITH 'INVALID'

    guide_list=list()

    # Was a PAM sequence detected?
    if sum(site > 0 for site in PAM_locations) > 0:

        # Sense strand: Build guide sequence if PAM is close to start site
        if orientation > 0:
            sites = (location for location in PAM_locations if location < 36 and location > 9)
            for location in sites:
                guide_list.append([location-31, 'sense', sequence[location-1:location+22]])

        # Antisense strand: Build guide sequence if PAM is close to start site
        else:
            sites = (location for location in PAM_locations if location < 51 and location > 24)
            for location in sites:
                guide_list.append([location-28, 'antisense', reverse_complement(sequence[location-20:location+3])])

    else:
        guide_list.append('NA')

    print(guide_list)
    return guide_list

def find_PAM(seq_list):
    # Find 'NGG' PAMs near start site in list of sequences

    print('Searching for PAMs.')

    # Identify indices of valid sequences
    for sequence in seq_list:
        if sequence[1] == 'Valid':
            row = seq_list.index(sequence)

            # Find where the transcript starts and ends in the gene sequence
            TSS = seq_list[row][2].find(seq_list[row][3][0:9])
            TES = seq_list[row][2].find(seq_list[row][3][-10:-1])

            # Copy gene sequence around start and end sites
            seq_start = seq_list[row][2][TSS-30:TSS+30]
            seq_end = seq_list[row][2][TES-30:TES+30]

            # Find all NGG sequences on both strands around start site
            PAM_start_sense = list(find_all(seq_start,'GG'))
            PAM_start_antisense = list(find_all(seq_start,'CC'))
            print(seq_list[row][0]+': a total of '+str(sum(i > 0 for i in PAM_start_sense) + sum(i > 0 for i in PAM_start_antisense))+' NGG PAMs found around the start site.')
            print('Viable PAMs at start site :')

            # Construct guide sequences for PAMs within proximity to start site
            build_guide(seq_start,PAM_start_sense,1)
            build_guide(seq_start,PAM_start_antisense,-1)

            # Find all NGG sequences on both strands around start site
            PAM_end_sense = list(find_all(seq_end,'GG'))
            PAM_end_antisense = list(find_all(seq_end,'CC'))
            print(seq_list[row][0]+': a total of '+str(sum(i > 0 for i in PAM_end_sense) + sum(i > 0 for i in PAM_end_antisense))+' NGG PAMs found around the end site.')
            print('Viable PAMs at the end site :')

            # Construct guide sequences for PAMs within proximity to start site
            build_guide(seq_start,PAM_end_sense,1)
            build_guide(seq_start,PAM_end_antisense,-1)

            # Append TSS, TES, PAM sites and guide sequences to sequence list
            seq_list[row][4:9] = [seq_start,PAM_start_sense, PAM_start_antisense,seq_end,PAM_end_sense,PAM_end_antisense]

    print('PAM seach complete.')
    return seq_list


### MAIN PROGRAM ###


id_csv = input('Transcript ID file : ')
print('Loading : '+id_csv)
seq_list = find_seq(id_csv)
seq_list_with_guides = find_PAM(seq_list)

# NEED TO IDENTIFY SEQUENCES WITH OFFTARGET BINDING

outputtext=NCBI.qblast('blastn','nt','CGGCCGCGCAATGGGCACCCGCG')
print(outputtext.read())
