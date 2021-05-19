#!/usr/bin/env python3

import os
import csv
import sys


DELIMITER = '\t'
BGC_ID_COLUMN_NAME = 'BGC'
CTG_ID_COLUMN_NAME = 'CONTIG'
ORF_ID_COLUMN_NAME = 'ORF'
DOMAIN_TYPE_COLUMN_NAME = 'A-ID'
PRED_COLUMN_NAME = 'PRED_TOP5'


def main():
    if len(sys.argv) != 4:
        print("Usage: %s original_full.tsv  updated_predictions.tsv  new_full.tsv" % sys.argv[0])
        sys.exit(1)
    original_full_tsv = sys.argv[1]
    updated_predictions_tsv = sys.argv[2]
    new_full_tsv = sys.argv[3]
    if not os.path.isfile(original_full_tsv) or not os.path.isfile(updated_predictions_tsv):
        print(f"Error! {original_full_tsv} and {updated_predictions_tsv} should exists!")
        sys.exit(1)
    if os.path.exists(new_full_tsv):
        print(f"Warning! {new_full_tsv} already exists! Will rewrite it")

    header = None
    original_content = []
    with open(original_full_tsv, 'r') as f:
        reader = csv.reader(f, delimiter=DELIMITER)
        for row in reader:
            if header is None:
                header = row
                bgc_id_column = header.index(BGC_ID_COLUMN_NAME)
                ctg_id_column = header.index(CTG_ID_COLUMN_NAME)
                orf_id_column = header.index(ORF_ID_COLUMN_NAME)
                domain_type_column = header.index(DOMAIN_TYPE_COLUMN_NAME)
                pred_column = header.index(PRED_COLUMN_NAME)
            else:
                original_content.append(row)

    updated_predictions_content = dict()  # BGC_ID --> list of rows
    updates_header = None
    with open(updated_predictions_tsv, 'r') as f:
        reader = csv.reader(f, delimiter=DELIMITER)
        for row in reader:
            if updates_header is None:
                updates_header = row
                updates_bgc_id_column = updates_header.index(BGC_ID_COLUMN_NAME)
                updates_ctg_id_column = updates_header.index(CTG_ID_COLUMN_NAME)
                updates_orf_id_column = updates_header.index(ORF_ID_COLUMN_NAME)
                updates_domain_type_column = updates_header.index(DOMAIN_TYPE_COLUMN_NAME)
                updates_pred_column = updates_header.index(PRED_COLUMN_NAME)
            else:
                bgc_id = row[updates_bgc_id_column]
                if bgc_id not in updated_predictions_content:
                    updated_predictions_content[bgc_id] = []
                updated_predictions_content[bgc_id].append(row)

    prev_bgc_id = None
    prev_uniq_domain_id = None
    bgc_id_row_counter = 0
    bgc_id_aux_row_counter = 0
    original_content.sort(key=lambda x: (x[bgc_id_column],
                                         x[ctg_id_column] if x[ctg_id_column].startswith('ctg') else 'zzz',
                                         x[orf_id_column],
                                         x[domain_type_column][0] if x[domain_type_column] else '',  # 'A' for A-domains and e.g. 'T' for 'TE'
                                         int(x[domain_type_column][1:]) if len(x[domain_type_column]) > 1 and
                                                                           x[domain_type_column][1:].isdigit() else 0))
    tmp = [row for row in original_content if row[bgc_id_column] == 'BGC0001214']
    tmp2 = updated_predictions_content['BGC0001214']
    with open(new_full_tsv, 'w') as f:
        csv_writer = csv.writer(f, delimiter=DELIMITER)
        csv_writer.writerow(header)
        for row in original_content:
            bgc_id = row[bgc_id_column]
            if bgc_id == prev_bgc_id:
                bgc_id_row_counter += 1
            else:
                prev_bgc_id = bgc_id
                bgc_id_row_counter = 0
                bgc_id_aux_row_counter = 0

            if bgc_id == ''.join(row):
                print(f"Warning! All values are empty for {bgc_id}. Skipping entire line")
                bgc_id_aux_row_counter += 1
                continue

            row[orf_id_column] += '_old'  # to distinguish from rows with updated ORF ID
            if bgc_id not in updated_predictions_content:
                if bgc_id_row_counter == 0:
                    print(f"Warning! {bgc_id} is absent in the update, will write previous lines for it as is")
                csv_writer.writerow(row)
            else:
                domain = row[domain_type_column]
                if not domain.startswith('A'):  # aux domain, e.g. 'TE'
                    bgc_id_aux_row_counter += 1
                    csv_writer.writerow(row)
                else:
                    uniq_domain_id = '_'.join([bgc_id, row[ctg_id_column], row[orf_id_column], domain])
                    if uniq_domain_id == prev_uniq_domain_id:  # special case: duplicating previous row (iterative NRP)
                        bgc_id_aux_row_counter += 1
                    else:
                        prev_uniq_domain_id = uniq_domain_id

                    if bgc_id_row_counter - bgc_id_aux_row_counter < 0 or \
                            bgc_id_row_counter - bgc_id_aux_row_counter >= len(updated_predictions_content[bgc_id]):
                        updates_domain = None
                    else:
                        updates_row = updated_predictions_content[bgc_id][bgc_id_row_counter - bgc_id_aux_row_counter]
                        updates_domain = updates_row[updates_domain_type_column]
                    if domain != updates_domain:
                        if not row[pred_column]:  # special case
                            print(f"Warning on {bgc_id}! No predictions in the row {bgc_id_row_counter} with A domain. "
                                  f"Keeping the row as is")
                            csv_writer.writerow(row)
                            bgc_id_aux_row_counter += 1
                        else:
                            print(f"Critical error on {bgc_id}! Domains does not match in the row {bgc_id_row_counter} "
                                  f"of the origin: {domain} vs {updates_domain}")
                            sys.exit(1)
                    else:
                        # bgc_id, ctg_id, domain should be the same for original and updated, so updating only orf_id, pred
                        row[orf_id_column] = updates_row[updates_orf_id_column]
                        row[pred_column] = updates_row[updates_pred_column]
                        csv_writer.writerow(row)

    print(f"Done! See results in\n{new_full_tsv}")


if __name__ == "__main__":
    main()
