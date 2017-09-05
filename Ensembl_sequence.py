import requests, sys, csv

def get_seq(id_csv):

    id_list=list(csv.reader(open("transcripts.csv",encoding = "utf-8-sig")))
    print(id_list)
    addr = "https://rest.ensembl.org/sequence/id/"
    header_gene="?content-type=application/json"
    header_transcript="?content-type=application/json;type=cds"

    item=0
    for  id_transcript in id_list:
        req_gene = requests.get(addr+id_transcript[0]+header_gene)
        req_transcript = requests.get(addr+id_transcript[0]+header_transcript)

        if (req_gene.ok == 0) or (req_transcript.ok == 0):
            id_list[item][1:3]=["Invalid","NA","NA"]
            item=+ 1
            continue
        seq_gene = req_gene.json()['seq']
        seq_transcript = req_transcript.json()['seq']
        id_list[item][1:3] = ["Valid",seq_gene,seq_transcript]
        item =+ 1
    return id_list

def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += 1

id_list=get_seq("transcripts.csv")

TSS=seq_gene.find(seq_transcript[0:9])
TES=seq_gene.find(seq_transcript[-10:-1])
seq_start=seq_gene[TSS-20:TSS+20]
seq_end=seq_gene[TES-20:TES+20]
print(list(find_all(seq_start,'GG')))
print(list(find_all(seq_end,'GG')))
