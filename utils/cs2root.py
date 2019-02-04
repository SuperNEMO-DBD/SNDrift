import csv
import ROOT as root

def read_field(filename, wt):
    with open(filename) as f:
        reader = csv.reader(f, delimiter=' ') 
        data = []
        numbersflag = False
        block = 0
        for row in reader:
            if row[0] == '#':
                if numbersflag:
                    block += 1
                    numbersflag = False
                continue
            else:
                numbersflag = True
                data.append((wt[block], float(row[0]), float(row[2]), float(row[4])))
    return data

weights = (0.95, 0.04, 0.01) # 95% He, 4% Ethanol, 1% Ar
alldata = read_field('data/trackergasCS.txt', weights)
nt = root.TNtupleD("cs","MagBoltz cross sections","weight:energy:csel:csinel")
n = len(alldata)
print ('read: ',n)

ff = root.TFile('data/trackergasCS.root','recreate')
for (weight, en,csel,csin) in alldata:
    nt.Fill(weight, en, csel, csin)

nt.Write()
ff.Close()

