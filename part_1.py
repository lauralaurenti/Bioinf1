import csv
import requests
import json
import pandas as pd

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
    tsv_out = open(output_filepath, 'w', newline='')

    reader = csv.reader(tsv_in, delimiter='\t')
    writer = csv.writer(tsv_out, delimiter='\t')

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
            writer.writerow([gene_id, gene_symb, '-', '-', '-', '-'])
            approved += 1

        if cnt % 10 == 0:
            print(cnt)
            sleep(1.5)

    tsv_in.close()
    tsv_out.close()

    print("Found {} NOT APPROVED genes:".format(len(not_approved)))
    print(not_approved)


def get_genes_ids(filepath):
    with open(filepath, 'r') as tsv_in:
        reader = csv.reader(tsv_in, delimiter='\t')
        _ = reader.__next__()

        genes_ids = []
        for line in reader:
            genes_ids.append(int(line[0].strip()))

    return genes_ids


def collect_interactions(input_filepath, seed_genes, output_filepath):
    print("Loading Biogrid table...")
    biogrid = pd.read_csv(input_filepath, sep='\t')

    # take only interactions between human genes
    print("Filtering human genes interactions...")
    biogrid = biogrid.loc[(biogrid['Organism ID Interactor A'] == 9606)
                          & (biogrid['Organism ID Interactor B'] == 9606)]

    # filter columns containint interactors IDs
    biogrid = biogrid[['Entrez Gene Interactor A', 'Entrez Gene Interactor B']]
    biogrid.columns = ['InteractorA', 'InteractorB']
    biogrid = biogrid.apply(pd.to_numeric)

    # filter interactions containing at leas 1 seed gene
    print("Filtering seed genes interactions...")
    seed_genes_inters = biogrid.loc[biogrid['InteractorA'].isin(seed_genes)
                                    | biogrid['InteractorB'].isin(seed_genes)]

    # create a list with all non seed genes interacting with seed genes
    all_genes = set(
        list(seed_genes_inters['InteractorA']) + list(seed_genes_inters['InteractorB']))
    non_seed_genes = [gene for gene in all_genes if gene not in seed_genes]

    # filter interaction between non seed genes
    print("Filtering non seed genes interactions...")
    non_seed_genes_inters = biogrid.loc[biogrid['InteractorA'].isin(non_seed_genes)
                                        & biogrid['InteractorB'].isin(non_seed_genes)]

    # concat all the interactions found
    interactions = pd.concat([seed_genes_inters, non_seed_genes_inters],
                             ignore_index=True)

    # remove duplicates
    print("Removing duplicates...")
    interactions['check_string'] = interactions.apply(lambda row: ''.join(
        sorted([str(row['InteractorA']), str(row['InteractorB'])])), axis=1)
    interactions.drop_duplicates('check_string', inplace=True)
    interactions.drop(columns=['check_string'], inplace=True)

    # save the interaction tsv file
    print("Sorting table...")
    interactions.sort_values(by=['InteractorA', 'InteractorB'], inplace=True)
    interactions.to_csv(output_filepath, sep='\t', index=False)
    print("### Finished ###")


def stats_summary(seed_genes, interactions_path):
    interactions = pd.read_csv(interactions_path, sep='\t')

    biogrid_seed_genes = set()
    all_biogrid_genes = set()
    for _, row in interactions.iterrows():
        geneA, geneB = row
        if geneA in seed_genes:
            biogrid_seed_genes.add(geneA)
        if geneB in seed_genes:
            biogrid_seed_genes.add(geneB)
        all_biogrid_genes.add(geneA)
        all_biogrid_genes.add(geneB)

    print("No. of Disgenet seed genes: ", len(seed_genes))
    print("No. of Biogrid seed genes: ", len(biogrid_seed_genes))
    print("Total no. of interacting genes: ", len(all_biogrid_genes))
    print("Total no. of interactions: ", interactions.shape[0])


if __name__ == "__main__":
    # filter_disease('data/curated_gene_disease_associations.tsv',
    #                'C0019208',
    #                'data/filtered_curated_gene_disease_associations.tsv')

    # fetch_HGNC('data/filtered_curated_gene_disease_associations.tsv',
    #            'data/approved_genes.tsv')

    genes = get_genes_ids('data/approved_genes.tsv')
    # collect_interactions('data/BIOGRID-ALL-4.2.191.tsv',
    #                      genes, 'data/interactions.tsv')

    stats_summary(genes, 'data/interactions.tsv')
