import csv
import itertools
import ROOT as root

def read_field(filename):
    with open(filename,"r") as f:
        lines = itertools.islice(f,9,None)
        reader = csv.reader(lines)
        data = []
        for line in reader:
            point = []
            for col in line:
                point.append(float(col))
            data.append(tuple(point))
    return data

data = read_field('sntracker.csv')
nt = root.TNtupleD("drift","Drift Field","x:y:ex:ey")
n = len(data)
print 'read: ',n

ff = root.TFile('sntracker_driftField.root','recreate')
for (x,y,ex,ey) in data:
    nt.Fill(x,y,ex,ey)

nt.Write()
ff.Close()

