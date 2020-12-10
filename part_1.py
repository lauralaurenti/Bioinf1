import csv
import requests
import json
import pandas as pd
import re

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


def add_uniprot_data(genes_path, uniprot_data_path):
    uniprot_data = {}

    with open(uniprot_data_path, "r") as fp_in:
        reader = csv.reader(fp_in, delimiter='\t')
        _ = reader.__next__()

        for line in reader:
            gene_sym = line[2]
            uniprot_data[gene_sym] = {
                "uniprotAC": line[0],
                "proteinName": line[3]
            }

    genes = pd.read_csv(genes_path, sep="\t")

    for idx, row in genes.iterrows():
        gene_sym = row['geneSymbol']
        try:
            gene_info = uniprot_data[gene_sym]
            genes.at[idx, 'uniprotAC'] = gene_info['uniprotAC']

            try:
                proteinName = re.match(r".+?(?=[\(,])",
                                       gene_info['proteinName']).group()
            except AttributeError:
                proteinName = gene_info['proteinName']

            genes.at[idx, 'proteinName'] = proteinName
        except KeyError:
            continue

    genes.to_csv(genes_path, sep="\t", index=False)


def get_genes_ids(filepath):
    with open(filepath, 'r') as tsv_in:
        reader = csv.reader(tsv_in, delimiter='\t')
        _ = reader.__next__()

        genes_ids = set()
        for line in reader:
            genes_ids.add(int(line[0].strip()))

    return genes_ids


def collect_interactions(input_filepath, seed_genes, output_filepath):
    print("Loading Biogrid table...")
    biogrid = pd.read_csv(input_filepath, sep='\t')

    # take only interactions between human genes
    print("Filtering human genes interactions...")
    biogrid = biogrid.loc[(biogrid['Organism ID Interactor A'] == 9606)
                          & (biogrid['Organism ID Interactor B'] == 9606)]

    # filter columns containint interactors IDs
    biogrid = biogrid[['Entrez Gene Interactor A', 'Entrez Gene Interactor B',
                       'Official Symbol Interactor A', 'Official Symbol Interactor B']]
    biogrid.columns = ['InteractorA', 'InteractorB', 'Sym_A', 'Sym_B']
    biogrid[['InteractorA', 'InteractorB']] = biogrid[['InteractorA',
                                                       'InteractorB']].apply(pd.to_numeric)

    # filter interactions containing at least 1 seed gene
    print("Filtering seed genes interactions...")
    seed_genes_inters = biogrid.loc[biogrid['InteractorA'].isin(seed_genes)
                                    | biogrid['InteractorB'].isin(seed_genes)]

    # create a list with all non seed genes interacting with seed genes
    all_genes = set(
        list(seed_genes_inters['InteractorA']) + list(seed_genes_inters['InteractorB']))
    non_seed_genes = set(
        [gene for gene in all_genes if gene not in seed_genes])

    with open('data/uniprot_non_seed_ids.txt', "w") as fp_in:
        for el in non_seed_genes:
            fp_in.write(str(el) + "\n")

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
        geneA, geneB, _, _ = row
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


def load_uniprot_non_seed():
    res = {}
    with open('data/uniprot_non_seed_genes_reviewed.tab', 'r') as fp_in:
        reader = csv.reader(fp_in, delimiter='\t')
        for line in reader:
            entry, sym = line
            res[sym] = entry
    return res


def arrange_interaction_data(approved_genes_path, interactions_path,
                             seed_genes_interactome_path, disease_interactome_path):
    approved_genes = pd.read_csv(approved_genes_path, sep='\t')
    seed_genes = set(approved_genes['geneId'])
    interactions = pd.read_csv(interactions_path, sep='\t')
    non_seed_genes = load_uniprot_non_seed()

    seed_out = open(seed_genes_interactome_path, 'w', newline='')
    writer_seed = csv.writer(seed_out, delimiter='\t')
    writer_seed.writerow(['interactor A gene symbol', 'interactor B gene symbol',
                          'interactor A Uniprot AC', 'interactor B Uniprot AC'])

    disease_out = open(disease_interactome_path, 'w', newline='')
    writer_disease = csv.writer(disease_out, delimiter='\t')
    writer_disease.writerow(['interactor A gene symbol', 'interactor B gene symbol',
                             'interactor A Uniprot AC', 'interactor B Uniprot AC'])

    for _, row in interactions.iterrows():
        geneA_id, geneB_id, geneA_sym, geneB_sym = row

        geneA_is_seed = geneA_id in seed_genes
        geneB_is_seed = geneB_id in seed_genes

        if geneA_is_seed or geneB_is_seed:

            # populate both seed genes interactome and disease interactome
            if geneA_is_seed and geneB_is_seed:
                geneA = approved_genes.loc[approved_genes['geneId']
                                           == geneA_id].iloc[0]
                geneB = approved_genes.loc[approved_genes['geneId']
                                           == geneB_id].iloc[0]

                geneA_AC = geneA['uniprotAC']
                geneB_AC = geneB['uniprotAC']

                writer_seed.writerow(
                    [geneA_sym, geneB_sym, geneA_AC, geneB_AC])
                writer_disease.writerow(
                    [geneA_sym, geneB_sym, geneA_AC, geneB_AC])

            # populate only disease interactome
            if geneA_is_seed ^ geneB_is_seed:
                if not geneA_is_seed:
                    geneB = approved_genes.loc[approved_genes['geneId']
                                               == geneB_id].iloc[0]
                    geneB_AC = geneB['uniprotAC']
                    geneA_AC = non_seed_genes.get(geneA_sym, "-")
                else:
                    geneA = approved_genes.loc[approved_genes['geneId']
                                               == geneA_id].iloc[0]
                    geneA_AC = geneA['uniprotAC']
                    geneB_AC = non_seed_genes.get(geneB_sym, "-")
                writer_disease.writerow(
                    [geneA_sym, geneB_sym, geneA_AC, geneB_AC])

    seed_out.close()
    disease_out.close()


def get_disease_interactome_gene_symbols(filepath):
    with open(filepath, 'r') as tsv_in:
        reader = csv.reader(tsv_in, delimiter='\t')
        _ = reader.__next__()

        genes_symbols = set()
        for line in reader:
            genes_symbols.add(line[0].strip())
            genes_symbols.add(line[1].strip())

    return genes_symbols


def get_Enrichr_userId(gene_list):
    ENRICHR_URL = 'http://maayanlab.cloud/Enrichr/addList'
    genes_str = '\n'.join(gene_list)
    payload = {
        'list': (None, genes_str),
    }

    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')

    data = json.loads(response.text)
    return data['userListId']


def fetch_Enrichr(user_list_id, gene_set_library):

    ENRICHR_URL = 'http://maayanlab.cloud/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'
    response = requests.get(
        ENRICHR_URL + query_string % (user_list_id, gene_set_library)
    )
    if not response.ok:
        raise Exception('Error fetching enrichment results')

    data = json.loads(response.text)
    temp = data[gene_set_library]
    temp = temp[:10]
    data[gene_set_library] = temp

    filepath = 'data/' + gene_set_library + '.json'
    with open(filepath, "w") as fp:
        json.dump(data, fp)


def parse_Enrichr_data(gst):
    print("Parsing " + gst)
    filepath = 'data/' + gst + '.json'
    with open(filepath, "r") as fp:
        data = json.load(fp)
    data = data[gst]

    cols = ["index", "name", "p-value",
            "adj p-value", "odds ratio",
            "combined score"]
    table_data = []

    for el in data:
        row = el[0:3] + [el[6]] + el[3:5]
        table_data.append(row)

    df = pd.DataFrame(table_data, columns=cols)
    df.set_index("index", inplace=True)
    tsv_path = 'data/' + gst + '.tsv'
    df.to_csv(tsv_path, sep="\t")
    print(df)
    print()


def cut_descriptions(filepath):
    from tempfile import NamedTemporaryFile
    import shutil

    tempfile = NamedTemporaryFile(mode='w', delete=False)

    with open(filepath, 'r') as tsvfile, tempfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        writer = csv.writer(tempfile, delimiter='\t')

        header = reader.__next__()
        writer.writerow(header)

        for row in reader:
            desc = row[4].split()
            new_desc = " ".join(desc[:20])
            new_row = row[:4] + [row[4]] + row[5:]
            writer.writerow(row)

    shutil.move(tempfile.name, filepath)


if __name__ == "__main__":
    # filter_disease('data/curated_gene_disease_associations.tsv',
    #                'C0019208',
    #                'data/filtered_curated_gene_disease_associations.tsv')

    # create a set of approved seed genes
    # fetch_HGNC('data/filtered_curated_gene_disease_associations.tsv',
    #            'data/approved_genes.tsv')

    # save genes IDs in a file for Uniprot
    # genes = get_genes_ids('data/approved_genes.tsv')
    # with open('data/uniprot_gene_ids.txt', 'w', newline="\n") as fp_in:
    #     for gene in genes:
    #         fp_in.write(str(gene) + "\n")

    # add uniprot information to approved genes
    # add_uniprot_data('data/approved_genes.tsv',
    #                  'data/uniprot_seed_genes_reviewed.tab')

    # genes = get_genes_ids('data/approved_genes.tsv')
    # collect_interactions('data/BIOGRID-ALL-4.2.191.tsv',
    #                      genes, 'data/interactions.tsv')

    # stats_summary(genes, 'data/interactions.tsv')

    # arrange_interaction_data('data/approved_genes.tsv',
    #                          'data/interactions.tsv',
    #                          'data/seed_genes_interactome.tsv',
    #                          'data/disease_interactome.tsv')

    # genes_symbols = get_disease_interactome_gene_symbols(
    #     'data/disease_interactome.tsv')
    # user_list_id = get_Enrichr_userId(genes_symbols)

    gene_set_libraries = [
        'KEGG_2019_Human',
        'GO_Biological_Process_2018',
        'GO_Molecular_Function_2018',
        'GO_Cellular_Component_2018'
    ]

    # for gst in gene_set_libraries:
    #     fetch_Enrichr(user_list_id, gst)

    for gst in gene_set_libraries:
        parse_Enrichr_data(gst)

    # with open('data/KEGG_2019_Human.json', "r") as fp:
    #     data = json.load(fp)

    # print
