import os
import csv

#BGC, CONTIG, ORF, A-ID, PRED_TOP5, START_POS, END_POS, STRAND
path_to_pred = "/home/olga/CAB/NRP/data/DataBase/MiBig_2019/MiBig2019BGCFilter/prediction/MINOWA/"

#создать и открыть для записи csv файл
csv_file = "out/mibig_detail.csv"
with open(csv_file, mode='w') as wf:
    csv_writer = csv.writer(wf, delimiter=',')

    csv_writer.writerow(['BGC', 'CONTIG', 'ORF', 'A-ID', 'PRED_TOP', 'START_POS', 'END_POS', 'STRAND'])
    #обойти все папочки которые есть в данной папке, вытащить их имена
    for BGC in os.listdir(path_to_pred):
        #загрузить в мапу информацию о начале, конце и стренде
        path_to_coord = path_to_pred + BGC + "/txt/" + BGC + "_gene.txt"
        bg_coord = dict()
        ed_coord = dict()
        strand = dict()
        if os.path.isfile(path_to_coord):
            with open(path_to_coord, mode='r') as rf:
                csv_reader = csv.reader(rf, delimiter='\t')
                for row in csv_reader:
                    bg_coord[row[0]] = row[1]
                    ed_coord[row[0]] = row[2]
                    strand[row[0]] = row[3]    

            #для каждого элемента посмотреть какие есть контиги, выписать все номера контиков по соотвествующим файликам
            nrpspred2 = path_to_pred + BGC + "/nrpspks_predictions_txt/"
            for pred_file in os.listdir(nrpspred2):
                if "nrpspredictor2_codes.txt" in pred_file:
                    with open(nrpspred2 + "/" + pred_file, "r") as rf:
                        for row in rf:
                            parts = row.split();
                            cur_contig = parts[0].split('_')[0]
                            cur_orf = parts[0].split('_')[1]
                            cur_aa = parts[0].split('_')[2]
                            cur_top5 = ';'.join(parts[-1].split(';')[:5])
                            cur_co = cur_contig + "_" + cur_orf
                            csv_writer.writerow([BGC, cur_contig, cur_orf, cur_aa, parts[-1], bg_coord[cur_co], ed_coord[cur_co], strand[cur_co]])
        else:
            print(BGC)
