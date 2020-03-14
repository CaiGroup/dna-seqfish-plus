import pandas as pd
from bisect import bisect
import numpy as np
import pickle
import networkx as nx
from pybktree import BKTree
from scipy.spatial import cKDTree as KDTree
from scipy.stats import norm

class psfDecoder:
    '''
    Decodes image from psf
    '''
    def __init__(self):
        self.alignedPSFs = None
        self.rounds = None
        # columns: Gene, Likelihood, Avecoords, AveAmp, rnpch, rnrow, rncol, rnamp
        self.acceptedBarcodes = pd.DataFrame(columns=('gene', 'row_ave', 'col_ave', 'amp_ave',
                                                      'r1pchn', 'r2pchn', 'r3pchan', 'r4pchan',
                                                      'r1row', 'r2row', 'r3row', 'r4row',
                                                      'r1col', 'r2col', 'r3col', 'r4col',
                                                      'r1amp', 'r2amp', 'r3amp', 'r4amp'))
        self.geneCounts = {}

    def load_barcode_key(self, filename):
        self.key = barcodeKey(filename, 0)

    def saveGenes(self, filename):
        self.acceptedBarcodes.to_csv('barcodes_' + filename)

        #make dataframe of gene counts in same order as key
        Gene = self.key.key['Gene']
        geneCountsDF = pd.DataFrame(index=Gene)
        geneCountsArray = []
        for i, keyRow in self.key.key.iterrows():
            geneName = keyRow['Gene']
            if geneName in self.geneCounts:
                geneCountsArray.append(self.geneCounts[geneName])
            else:
                geneCountsArray.append(0)
        geneCountsDF['Counts'] = geneCountsArray
        geneCountsDF.to_csv('counts_' + filename)

    def findRoundGraph(self, alignedPSFsFile, nRounds, maxError):
        '''
        Initializes objects, makes barcodes from aligned psfs dataframe
        :param alignedPSFs:
        '''
        self.alignedPSFs = pd.read_csv(alignedPSFsFile, sep=',', header=0, index_col=('Hyb', 'Dot'))
        self.roundMatchDotGraph = nx.DiGraph()
        self.rounds = []
        previousRound = None

        for round in range(nRounds):
            self.rounds.append(roundDots(round, self.roundMatchDotGraph, self.alignedPSFs, maxError, previousRound))
            previousRound = self.rounds[round]

    def estimate_Barcodes(self):
        return

    def saveRounds(self, filename):
        with open(filename, 'wb') as f:
            pickle.dump(self.rounds, f)

    def loadRounds(self, filename):
        with open(filename, 'rb') as f:
            self.rounds = pickle.load(f)
        self.roundMatchDotGraph = self.rounds[0].roundDotGraph

    def trimMG(self):
        print('trimming')
        toRemove = [g for g in nx.weakly_connected_components(self.roundMatchDotGraph) if len(g) < 3 or len(g) > 20]
        for g in toRemove:
            for node in list(g,):
                self.roundMatchDotGraph.remove_node(node)

    def decode(self):
        #print('decoding')

        for i, subGraph in enumerate(nx.weakly_connected_component_subgraphs(self.roundMatchDotGraph)):
        #for i, subGraph in enumerate(self.roundMatchDotGraph.weakly_connected_component_subgraphs()):
            #print(i)
            # initialize list for matching barcodes. After all matching barcodes in the graph are found, conflicts will
            # be resolved
            putative_bcs = []
            p_bc_info = []
            #paths = nx.shortest_path(subGraph, weight='loglik')
            #for path in paths:
            for source, (tldict, tpdict) in nx.all_pairs_dijkstra(subGraph):#,weight='loglik'):
                #print(source)

                for target in tldict:
                    if tldict[target] > 2:
                        #check if path corresponds to barcode
                        #print('checking')
                        barcode = self.getBarcodeFromNodes(tpdict[target])
                        match = self.key.findCloseBarcode(barcode)
                        if match:
                            print('matched')
                            putative_bcs.append(set(tpdict[target]))
                            p_bc_info.append(match)
                            
            #print('len putative barcodes:', len(p_bc_info))

            # find conflicting dots and barcodes
            if putative_bcs:
                conflicting_barcodes = []
                conflicts = set()
                union = set()
                for bc in putative_bcs:
                    conflicts |= (union & bc)
                    union |= bc

                for i, bc in enumerate(putative_bcs):
                    conflicting = bc & conflicts
                    if conflicting:
                        conflicting_barcodes.append(bc)
                    else: # accept
                        print("Accept:", p_bc_info[0])
                        self.addAcceptedBarcode(p_bc_info[0], list(bc))

        #subGraphNodes = [list(g) for g in nx.weakly_connected_component_subgraphs(self.roundMatchDotGraph)]

    def getBarcodeFromNodes(self, nodes):
        '''
        From a list of the points in a path, get teh barcode ready to look up in the barcode key
        :param nodes: list of nodes in a putative barcode.
        :return:
        4 element tuple containing the psuedochannels of each point in the barcode.
        '''
        if len(nodes) < 3:
            return None
        bc = [None]*4
        for node in nodes:
            bc[self.roundMatchDotGraph.node[node]['rnd']] = self.roundMatchDotGraph.node[node]['pseudochannel']

        return tuple(bc)

    def addAcceptedBarcode(self, info, bcl):
        '''
        Add accepted barcode to accepted PSF dataframe
        :param gene: String name of gene
        :param bcl: list of tuples containing the coordinates of the point.
        :return:
        '''
        gene = info[0]
        barcode = info[1]
        hd = info[2]
        row_ave = 0
        col_ave = 0
        amp_ave = 0
        g = self.roundMatchDotGraph.node
        for p in bcl:
            row_ave += p[0]
            col_ave += p[1]
            amp_ave += g[p]['amp']
        row_ave = row_ave / len(bcl)
        col_ave = col_ave / len(bcl)
        amp_ave = amp_ave / len(bcl)
        bcn = [None]*4

        #order barcodes
        bcn[g[bcl[0]]['rnd']] = bcl[0]
        bcn[g[bcl[1]]['rnd']] = bcl[1]
        bcn[g[bcl[2]]['rnd']] = bcl[2]
        bcn[g[bcl[3]]['rnd']] = bcl[3]

        self.acceptedBarcodes = self.acceptedBarcodes.append(pd.DataFrame([{
            'gene': gene, 'row_ave': row_ave, 'col_ave': col_ave, 'amp_ave': amp_ave,
            'r1pchn':  g[bcn[0]]['pseudochannel'], 'r2pchn': g[bcn[1]]['pseudochannel'],
            'r3pchan': g[bcn[2]]['pseudochannel'], 'r4pchan': g[bcn[3]]['pseudochannel'],
            'r1row': bcn[0][0], 'r2row': bcn[1][0], 'r3row': bcn[2][0], 'r4row': bcn[3][0],
            'r1col': bcn[0][1], 'r2col': bcn[1][1], 'r3col': bcn[2][1], 'r4col': bcn[3][1],
            'r1amp': g[bcn[0]]['amp'], 'r2amp': g[bcn[1]]['amp'], 'r3amp': g[bcn[2]]['amp'], 'r4amp': g[bcn[3]]['amp']
        }]), ignore_index=True)

        if gene in self.geneCounts:
            self.geneCounts[gene] += 1
        else:
            self.geneCounts[gene] = 1

class roundDots:
    '''
    Inserts psfs in a round into roundDotGraph and provides searching for PSfs with a KDTree within the round.
    '''
    def __init__(self, rnd, roundDotGraph, psfData, maxError = 1,  previousRound = None):
        '''
        Contsructor. Saves subset of psf dataframe for round, which is made searchable by a KDTree, populates
        roundDOtGraph with dots from this round, and draws edges from close dots in the last round to dots in this round.
        :param rnd: int round number. Used to slice psfData DataFrame (see below)
        :param roundDotGraph: The graph connecting PSFs that are close to one another between round to allow
            for matching barcodes.
        :param psfData: Pandas Dataframe containing information on all PSFs
        :param previousRound: If this is not the first round, a reference to the previous round to draw edges in roundDotGraph.
        '''
        self.maxError = maxError
        self.maxErrorSq = maxError**2
        self.round = rnd
        self.roundDotGraph = roundDotGraph
        startHyb = rnd * 20
        stopHyb = startHyb+19
        self.psfs = psfData.loc[startHyb:stopHyb]
        self.dotTree = KDTree(self.psfs.loc[:,'row':'col'])
        self.previousRound = previousRound
        for i, psf in self.psfs.iterrows():
            coords = (psf['row'], psf['col'])
            hyb = i[0] # since PSFs dataframe is multiIndexed; hyb then index of dot in hyb
            self.roundDotGraph.add_node(coords, amp=psf.amp, RMSE=psf.RMSE, rnd=rnd, pseudochannel=hyb+1-rnd*20)
            if previousRound:
                # add edges connecting psfs from previous round within search radius
                matches = previousRound.matchPSF(psf)
                for j, match in matches.iterrows():
                    # Calculate Log Likelihood that both PSFs were generated by the same object based and position and amplitude
                    # w = (match['row']-psf['row'])**2 + (match['col']-psf['col'])**2
                    sigma = psf['amp']/7.5
                    likelihood_amp = (1 - norm.cdf(match['amp'], psf['amp'], sigma))*2
                    #if likelihood_amp > 0.001:
                    llikelihood_amp = np.log2(likelihood_amp)
                    matchCoords = (float(match['row']), float(match['col']))
                    # give set loglik as negative, because path finding algorithms minimize path
                    # so they can calculate maxlikelihood if set negative
                    self.roundDotGraph.add_edge(matchCoords, coords, loglik=-llikelihood_amp)

    def matchPSF(self, psf):
        '''
        Checks for dots in round that match input dot to within error
        :param psf:
        :return:
        list of matching dots from round
        '''
        matches = self.dotTree.query_ball_point((psf.row, psf.col), self.maxError)
        return self.psfs.iloc[matches]


class barcodeKey:
    '''
    Searchable datastructure for matching barcodes
    '''

    def __init__(self, barcode_key_file, maxHamDist):
        self.maxHamDist = maxHamDist
        self.key = pd.read_csv(barcode_key_file, header=0)#, names=('gene', 'id', 'r1', 'r2', 'r3', 'r4'))#0)
        self.bktree = BKTree(hamming_distance)
        self.barcodeDict = {}
        for i, row in self.key.iterrows():
            barcode = (row['R1Pchan'], row['R2Pchan'], row['R3Pchan'], row['R4Pchan'])
            self.bktree.add(barcode)
            self.barcodeDict[barcode] = row['Gene']


    def findCloseBarcode(self, barcode):
        ret = self.bktree.find(barcode, self.maxHamDist)
        if not ret:
            return None
        hd=ret[0][0]
        bcMatch = ret[0][1]
        #hd = hamming_distance(barcode, closestMatch)
        gene = self.barcodeDict[bcMatch]
        return gene, bcMatch, hd


def hamming_distance(bc1, bc2):
    '''
    finds the hamming distance between two possible barcodes
    :param bc1: 4 element np array of int pseudocolors
    :param bc2: 4 element np array of int pseudocolors
    :return:
    hamming distance
    '''
    hd = 0
    for i in range(4):
        # increment hamming distance if either element is empty, or if they are unequal
        if bc1[i] == -1 or bc2[i] == -1 or bc1[i] != bc2[i]:
            hd += 1
    return hd
