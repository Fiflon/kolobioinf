from scipy.stats import genextreme
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from random import randint
from Bio.SeqRecord import SeqRecord
from io import StringIO


ax = plt.subplot()

c=-0.2

r = genextreme.rvs(c,size=10,loc=8000,scale=500)

r=r.astype(int)

print(r)
print(min(r))

ax.hist(r,bins=50, density=True, histtype="stepfilled")

caly_genom = SeqIO.read("sequence.fasta", "fasta")
i = 0
fragmenty = []
limit = len(caly_genom.seq)
for dlugosc in r:
    if dlugosc < 0:
        dlugosc = dlugosc * (-1)

    if dlugosc >= limit:
        dlugosc = randint(0, limit)
    start = randint(0, limit - dlugosc)
    end = start + dlugosc
    fragment = caly_genom.seq[start:end]
    record = SeqRecord(fragment, "odczyt_%i_s%i_d%i.fasta" % (i + 1,start,dlugosc), "", "")
    fragmenty.append(record)
    i = i+1

SeqIO.write(fragmenty, "fragmenty.fasta", "fasta")