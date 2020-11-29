import csv
import requests
import json
from time import sleep

from urllib.parse import urlparse


def filter_disease(input_filepath, disease_id, output_filepath):
    tsv_in = open(input_filepath, 'r')
    tsv_out = open(output_filepath, 'w', newline='')

    reader = csv.reader(tsv_in, delimiter='\t')
    writer = csv.writer(tsv_out, delimiter='\t')

    header = reader.__next__()
    disease_id_idx = header.index('diseaseId')
    writer.writerow(header)
    
    print("Filtering genes for disease {}".format(disease_id))
    cnt = 0
    for line in reader:
        temp_disease_id = line[disease_id_idx].strip()
        if temp_disease_id == disease_id:
            writer.writerow(line)
            cnt += 1

    tsv_in.close()
    tsv_out.close()
    print("Filtering finished")
    print("{} genes found".format(cnt))


def get_gene_HGNC(gene_symb):
    headers = {
    'Accept': 'application/json',
    }

    uri = 'http://rest.genenames.org'
    path = '/search/symbol/' + gene_symb + '+AND+status:Approved'
    url = urlparse(uri+path).geturl()

    response = requests.get(url=url, headers=headers)

    if response.status_code == 200: 
        data = json.loads(response.text)
        return data['response']
    else:
        print('Error detected: ' + response.status_code)
        return None


def fetch_HGNC(input_filepath, output_filepath):
    tsv_in = open(input_filepath, 'r')
    csv_out = open(output_filepath, 'w', newline='')

    reader = csv.reader(tsv_in, delimiter='\t')
    writer = csv.writer(csv_out, delimiter=',')

    header = reader.__next__()
    gene_id_idx = header.index('geneId')
    gene_symb_idx = header.index('geneSymbol')
    writer.writerow(['geneId', 'geneSymbol', 'uniprotAC',
                     'proteinName', 'description', 'notes'])

    cnt = 0
    approved = 0
    not_approved = []
    for line in reader:
        gene_id = line[gene_id_idx]
        gene_symb = line[gene_symb_idx]
        res = get_gene_HGNC(gene_symb)
        cnt += 1
        
        if res['numFound'] == 0:
            not_approved.append(gene_symb)
        else:
            writer.writerow([gene_id, gene_symb, '', '', '', ''])
            approved += 1

        if cnt % 10 == 0:
            print(cnt)
            sleep(1.5)


    tsv_in.close()
    csv_out.close()

    print("Found {} NOT APPROVED genes:".format(len(not_approved)))
    print(not_approved)


if __name__ == "__main__":
    # filter_disease('data/curated_gene_disease_associations.tsv',
    #                'C0019208',
    #                'data/filtered_curated_gene_disease_associations.tsv')

    fetch_HGNC('data/filtered_curated_gene_disease_associations.tsv',
               'data/approved_genes.csv')
