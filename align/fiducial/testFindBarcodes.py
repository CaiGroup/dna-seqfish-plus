from findBarcodes import psfDecoder
import os
import networkx as nx

#os.chdir('/Users/jonathanwhite/Documents/NIH3T3')
os.chdir('C:\\Users\\jonat\\Documents\\NIH3T3')
#os.chdir('D:\decoding')
decoder = psfDecoder()

nRounds = 4

alignment = 'NIH3T3_bnd1000_thr1000_Pos1_z1_ch1_Dots_aligned.csv'

decoder.findRoundGraph(alignment, nRounds, 2)

decoder.saveRounds('NIH3T3_bnd1000_thr1000_Pos1_z1_ch1_Dots_barcodes_mError_2p.pkl')

nPSFs = []

#for bc in decoder.rounds.barcodes:
#    nPSFs.append(sum(bc.hybs))
'''
import matplotlib.pyplot as plt

#plt.hist(nPSFs)
#plt.show()

#for g in nx.weakly_connected_components(decoder.roundMatchDotGraph): print(len(g))
print(decoder.roundMatchDotGraph.number_of_nodes())
wcc = list(nx.weakly_connected_components(decoder.roundMatchDotGraph))

lwcc = [len(g) for g in wcc]

plt.hist(lwcc)
plt.show()
'''
decoder.roundMatchDotGraph.number_of_nodes()