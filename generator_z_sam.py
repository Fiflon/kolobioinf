from scipy.stats import genextreme
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from random import randint
from Bio.SeqRecord import SeqRecord
from io import StringIO


ax = plt.subplot()

c=-0.2

r = genextreme.rvs(c,size=10000,loc=8000,scale=500)

r=r.astype(int)
ax.hist(r,bins=50, density=True, histtype="stepfilled")

caly_genom = SeqIO.read("kng_reg01.fasta", "fasta")
i = 0
suma = 0
fragmenty = []
limit = len(caly_genom.seq)
plik_sam = open("fragmenty.sam","w")
tmp=""
plik_sam.write("@SQ\tSN:kng_reg01\tLN:1000000\n")
for dlugosc in r:
    if dlugosc < 0:
        dlugosc = dlugosc * (-1)

    if dlugosc >= limit:
        dlugosc = randint(0, limit)
    start = randint(0, limit - dlugosc)
    end = start + dlugosc
    fragment = caly_genom.seq[start:end]
    losowanieCzyOdwrotniekom = randint(0,1)
    if losowanieCzyOdwrotniekom == 1:
        fragment = caly_genom.seq[start:end].reverse_complement()

    record = SeqRecord(fragment, "odczyt_%i_s%i_d%i_rc%i.fasta" % (i + 1,start,dlugosc,losowanieCzyOdwrotniekom), "", "")
    fragmenty.append(record)
    i = i+1
    tmp="odczyt_"+str(i)+"_s"+str(start)+"_d"+str(dlugosc)+"_rc"+str(losowanieCzyOdwrotniekom)+".fasta"
    plik_sam.write(">"+tmp+"\t0\tkng_reg01\t"+str(start)+"\t60\t"+str(dlugosc)+"M\t*\t0\t0\t"+str(fragment)+"\t*\n")
    suma+=dlugosc
    if suma>20000000:
        break

SeqIO.write(fragmenty, "fragmenty.fasta", "fasta")