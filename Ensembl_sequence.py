import requests, sys, csv, os

def get_seq(id_csv):
    id_list=list(csv.reader(open(id_csv,encoding = "utf-8-sig")))
    addr = "https://rest.ensembl.org/sequence/id/"
    header_gene="?content-type=application/json"
    header_transcript="?content-type=application/json;type=cds"

    for  id_transcript in id_list:
        row = id_list.index(id_transcript)
        req_gene = requests.get(addr+id_transcript[0]+header_gene)
        req_transcript = requests.get(addr+id_transcript[0]+header_transcript)

        if (req_gene.ok == 0) or (req_transcript.ok == 0):
            id_list[row][1:3]=["Invalid","NA","NA"]
            continue
        seq_gene = req_gene.json()['seq']
        seq_transcript = req_transcript.json()['seq']
        id_list[row][1:3] = ["Valid",seq_gene,seq_transcript]
    return id_list

def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += 1

def get_PAM(id_list,PAM):
    for s in id_list:
        if s[1] == 'Valid':
            row = id_list.index(s)

            print(s[3])

id_list=get_seq("transcripts.csv")

TSS=seq_gene.find(seq_transcript[0:9])
TES=seq_gene.find(seq_transcript[-10:-1])
seq_start=seq_gene[TSS-20:TSS+20]
seq_end=seq_gene[TES-20:TES+20]
print(list(find_all(seq_start,'GG')))
print(list(find_all(seq_end,'GG')))
