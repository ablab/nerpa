import csv
path="/media/hosein/My Passport/hosein/Desktop/project/sequence_data/bacteria_all/"

with open('ball.csv') as csvfile:
    with open('ball_org.csv', 'w') as csvw:
        reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        writer = csv.writer(csvw, delimiter=',', quotechar='"')
        for row in reader:
            print(row)
            nrow = row
            if (nrow[0] == 'Score'):
                nrow.append('Peptide predicted organism (from the corresponding NRPSprediction file)')
            else:
                nrow[-2] = nrow[-2].split('/')[0]
                filename = nrow[-2] + '.fna'
                line = open(path + filename).readline()
                nrow.append(' '.join((line.split('genome')[0].split())[1:]) )
            print(nrow)
            writer.writerow(nrow)
        
                
                
                
        
    
