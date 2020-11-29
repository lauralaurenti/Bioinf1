import csv


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


if __name__ == "__main__":
    filter_disease('data/curated_gene_disease_associations.tsv',
                   'C0019208',
                   'data/filtered_curated_gene_disease_associations.tsv')
