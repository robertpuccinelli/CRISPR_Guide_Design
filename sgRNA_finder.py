import requests, sys, csv, os

def find_seq(id_csv):
    seq_list=list(csv.reader(open(id_csv,encoding = "utf-8-sig")))
    addr = "https://rest.ensembl.org/sequence/id/"
    header_gene="?content-type=application/json"
    header_transcript="?content-type=application/json;type=cds"

    for  id_transcript in seq_list:
        print('Searching for sequences : '+id_transcript[0])
        row = seq_list.index(id_transcript)
        req_gene = requests.get(addr+id_transcript[0]+header_gene)
        req_transcript = requests.get(addr+id_transcript[0]+header_transcript)

        if (req_gene.ok == 0) or (req_transcript.ok == 0):
            seq_list[row][1:3]=["Invalid","NA","NA"]
            print('Transcript or gene sequence NOT found.')
            continue
        seq_gene = req_gene.json()['seq']
        seq_transcript = req_transcript.json()['seq']
        seq_list[row][1:3] = ["Valid",seq_gene,seq_transcript]
        print('Transcript and gene sequence found.')

    print('Sequence search complete.')
    return seq_list

def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += 1

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])

def build_guide(sequence,PAMs,orientation):
    guide=list()
    if sum(site > 0 for site in PAMs) > 0:
        if orientation > 0:
            sites = (location for location in PAMs if location < 51 and location > 31)
            for location in sites:
                guide.append([location-31, 'sense', sequence[location-22:location+1]])
        else:
            sites = (location for location in PAMs if location < 32 and location > 12)
            for location in sites:
                guide.append([location-28, 'antisense', reverse_complement(sequence[location:location+23])])
    else:
        guide.append('NA')

    print(guide)

    return

def find_PAM(PAM_list):
    print('Searching for PAMs.')
    for s in PAM_list:
        if s[1] == 'Valid':
            row = PAM_list.index(s)
            TSS = PAM_list[row][2].find(PAM_list[row][3][0:9])
            TES = PAM_list[row][2].find(PAM_list[row][3][-10:-1])
            seq_start = PAM_list[row][2][TSS-30:TSS+30]
            seq_end = PAM_list[row][2][TES-30:TES+30]

            PAM_start_sense = list(find_all(seq_start,'GG'))
            PAM_start_antisense = list(find_all(seq_start,'CC'))
            print(PAM_list[row][0]+': a total of '+str(sum(i > 0 for i in PAM_start_sense) + sum(i > 0 for i in PAM_start_antisense))+' NGG PAMs found around the start site.')
            print('Viable PAMs at start site :')
            build_guide(seq_start,PAM_start_sense,1)
            build_guide(seq_start,PAM_start_antisense,-1)

            PAM_end_sense = list(find_all(seq_end,'GG'))
            PAM_end_antisense = list(find_all(seq_end,'CC'))
            print(PAM_list[row][0]+': a total of '+str(sum(i > 0 for i in PAM_end_sense) + sum(i > 0 for i in PAM_end_antisense))+' NGG PAMs found around the end site.')
            print('Viable PAMs at the end site :')
            build_guide(seq_start,PAM_end_sense,1)
            build_guide(seq_start,PAM_end_antisense,-1)

            PAM_list[row][4:9] = [seq_start,PAM_start_sense, PAM_start_antisense,seq_end,PAM_end_sense,PAM_end_antisense]
    print('PAM seach complete.')
    return PAM_list
ENST00000261890.6
print()
print(PAM_list)
id_csv = input('Transcript ID file : ')
print('Loading : '+id_csv)
seq_list = find_seq(id_csv)
PAM_list = find_PAM(seq_list)
outputtext=NCBI.qblast('blastn','nt','CGGCCGCGCAATGGGCACCCGCG')
print(outputtext.read())
