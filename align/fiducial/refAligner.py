import pandas as pd
import numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree as KDTree
from copy import copy

class RefAligner:
    """
    This Class aligns dots in a seqFISH experiment image containing fiducial markers, to a reference image containing
    only the fiducial markers
    """
    def __init__(self, ro, ref, no_hybs = False):
        """
        Initialize RefAligner object. Save pandas DataFrames of reference and readout points
        :param ro: string name of reodout points csv file.
        :param ref: string name of reference points csv file.
        """

        if no_hybs:
            self.ro = pd.read_csv(ro)
            self.ro['Hyb'] = 1
            self.ro.set_index('Hyb', inplace=True)
        else:
            self.ro = pd.read_csv(ro, index_col='Hyb')
        self.ref = pd.read_csv(ref)
        self.ref['n_trav_matched'] = 0
        self.ref['n_unmatched'] = 0

        self.ref.sort_values(by='row', inplace=True)

        # define data members that will be used in alignment
        self.matchDict = {}
        self.edgeMatches = set()
        self.trav_matches = set()
        self.ref_dots = None
        self.reserved_ref_dots = None
        self.ro_hyb = None
        self.ro_lat_tree = None
        self.trav_ind = 0
        self.traversal_queue = []
        self.n_ambiguous = 0
        self.n_trav_matched = 0
        self.offsets = None
        self.matchesDF = pd.DataFrame(columns=('ref_row', 'ref_col', 'ref_z', 'comp_row', 'comp_col', 'comp_z',
                                               'aligned_row', 'aligned_col', 'aligned_z'))

        # define threshold parameters set by their own methods
        self.xyse = None
        self.xyse_sq = None
        self.zse = None
        self.xyte = None
        self.xyte_sq = None
        self.zte = None
        self.xyme = None
        self.xyme_sq = None
        self.zme = None
        self.min_edge_matches = None
        self.min_bright_prop = None
        self.max_bright_prop = None
        self.n_umnatch_give_up = None
        self.min_dot_matches = None
        self.n_longest_edges = None

    def _find_all_ref_edges(self):
        """
        returns sorted dataframe describing edges in reference graph.
        :return:
        """
        print('Finding reference edges')
        n_dots = len(self.ref)
        #edges = []
        n_edges = int(((n_dots-1)**2-(n_dots-1))/2)
        edge_array = np.zeros([n_edges, 12])
        i = 0
        for u_i in range(n_dots):
            for v_i in range(u_i + 1, n_dots):
                #edge = []
                # add row, col, z, amp of each node to edge list
                edge_array[i,0:4] = np.array(self.ref.iloc[u_i]['row':'amp'])
                #edge = list(self.ref.iloc[u_i])
                #edge += list(self.ref.iloc[v_i])
                edge_array[i, 4:8] = np.array(self.ref.iloc[v_i]['row':'amp'])

                #now add col_dist, row_dist and absolute distance
                rdist = self.ref.iloc[v_i]['row'] -self.ref.iloc[u_i]['row']
                cdist = self.ref.iloc[v_i]['col'] - self.ref.iloc[u_i]['col']
                zdist = self.ref.iloc[v_i]['z'] - self.ref.iloc[u_i]['z']
                length = np.sqrt(rdist**2 + cdist**2)
                #edge += [rdist, cdist, zdist, length]
                edge_array[i, 8:12] = np.array([rdist, cdist, zdist, length])
                #edges.append(edge)
        edgedf = pd.DataFrame(data=edge_array, columns=['u_row', 'u_col', 'u_z', 'u_amp', 'v_row', 'v_col', 'v_z', 'v_amp',
                                                   'row_dist', 'col_dist', 'z_dist', 'length'])
        edgedf.sort_values(by='length', ascending=False, inplace=True)

        return edgedf

    def _find_ref_edges(self, n_longest):
        n_corner = int(np.ceil(np.sqrt(n_longest)))
        # find dots closest to each corner
        print('Finding Reference Edges')
        ref_dists = copy(self.ref)
        ref_dists['ul_dist'] = 0
        ref_dists['ur_dist'] = 0
        ref_dists['bl_dist'] = 0
        ref_dists['br_dist'] = 0
        for i, dot in self.ref.iterrows():
            ref_dists['ul_dist'].loc[i] = dot.row**2 + dot.col**2
            ref_dists['ur_dist'].loc[i] = dot.row**2 + (2048 - dot.col)**2
            ref_dists['bl_dist'].loc[i] = (2048 - dot.row)**2 + dot.col**2
            ref_dists['br_dist'].loc[i] = (2048 - dot.row)**2 + (2048 - dot.col)**2

        ul_dots = ref_dists.sort_values(by='ul_dist').iloc[:n_corner]
        ur_dots = ref_dists.sort_values(by='ur_dist').iloc[:n_corner]
        bl_dots = ref_dists.sort_values(by='bl_dist').iloc[:n_corner]
        br_dots = ref_dists.sort_values(by='br_dist').iloc[:n_corner]

        edges_array = np.zeros([2*n_corner**2, 12])#([2*n_corner**2, 12])
        i = 0
        for j, ul_dot in ul_dots.iterrows():
            for k, br_dot in br_dots.iterrows():
                edges_array[i, 0:4] = np.array(ul_dot['row':'amp'])
                edges_array[i, 4:8] = np.array(br_dot['row':'amp'])

                # now add col_dist, row_dist and absolute distance
                r_dist = br_dot['row'] - ul_dot['row']
                c_dist = br_dot['col'] - ul_dot['col']
                z_dist = br_dot['z'] - ul_dot['z']
                length = np.sqrt(r_dist ** 2 + c_dist ** 2)
                edges_array[i, 8:12] = np.array([r_dist, c_dist, z_dist, length])
                i += 1

        for j, ur_dot in ur_dots.iterrows():
            for k, bl_dot in bl_dots.iterrows():
                edges_array[i, 0:4] = np.array(ur_dot['row':'amp'])
                edges_array[i, 4:8] = np.array(bl_dot['row':'amp'])

                # now add col_dist, row_dist and absolute distance
                r_dist = bl_dot['row'] - ur_dot['row']
                c_dist = bl_dot['col'] - ur_dot['col']
                z_dist = bl_dot['z'] - ur_dot['z']
                length = np.sqrt(r_dist ** 2 + c_dist ** 2)
                edges_array[i, 8:12] = np.array([r_dist, c_dist, z_dist, length])
                i += 1

        edges_df = pd.DataFrame(data=edges_array,
                              columns=['u_row', 'u_col', 'u_z', 'u_amp', 'v_row', 'v_col', 'v_z', 'v_amp',
                                       'row_dist', 'col_dist', 'z_dist', 'length'])
        edges_df.sort_values(by='length', ascending=False, inplace=True)

        edges_df = edges_df.iloc[:n_longest]

        return edges_df

    # Parameter setting methods
    def set_n_longest_edges(self, nle, all_pairs=False):
        self.n_longest_edges = nle
        if all_pairs:
            self.ref_edges = self._find_all_ref_edges()
        else:
            self.ref_edges = self._find_ref_edges(nle)

    def set_xy_search_error(self, xyse):
        self.xyse = xyse
        self.xyse_sq = xyse**2

    def set_z_search_error(self, zse):
        self.zse = zse

    def set_xy_traversal_error(self, xyte):
        self.xyte = xyte

    def set_z_traversal_error(self, zte):
        self.zte = zte

    def set_xy_match_error(self, xyme):
        self.xyme = xyme
        self.xyme_sq = xyme**2

    def set_z_match_error(self, zme):
        self.zme = zme

    def set_min_edge_match(self, min_edge_matches):
        self.min_edge_matches = min_edge_matches

    def set_min_bright_prop(self, mbp):
        self.min_bright_prop = mbp

    def set_max_bright_prop(self, mbp):
        self.max_bright_prop = mbp

    def set_n_unmatch_give_up(self, nugu):
        self.n_umnatch_give_up = nugu

    def set_min_dot_matches(self,mdm):
        self.min_dot_matches = mdm

    def align(self, hybs=None, saveMatches=False, drop=True):
        """
        aligns hybridizations in the readout points to the reference
        :param hybs: list of hybridizations to align. If None, aligns all hybridizations.
        :return:
        """
        if not hybs:
            hybs = list(range(1, 1 + max(self.ro.index)))

        offsets = []
        for hyb in hybs:
            hyb_offsets = self.align_hyb(hyb)
            if saveMatches:
                self.add_hyb_to_match_df(hyb, hyb_offsets)
            hyb_offsets = [hyb] + hyb_offsets
            offsets.append(hyb_offsets)
        
        self.offsets = pd.DataFrame(offsets, columns=('hyb', 'row', 'col', 'z', 'row_SE', 'col_SE', 'z_SE', 'n_matches'))
        return self.offsets
        #self.offsets.set_index('hyb')
        #if __name__ == '__main__':
            #pool = Pool(n_processes)
            #pool.map(self.align_hyb, hybs)

    def align_par(self, hybs=None, n_processes=1):

        if not hybs:
            hybs = list(range(1, 1 + max(self.ro.index)))

        if __name__ == '__main__':
            pool = Pool(n_processes)
            hyb_offsets = pool.map(self.align_hyb, hybs)
            pool.close()
            pool.join()
            self.offsets.DataFrame(hyb_offsets, columns=('row', 'col', 'z', 'row_SE', 'col_SE', 'z_SE', 'n_matches', 'hyb'))
            self.offsets.set_index('hyb', inplace=True)

    def align_hyb(self, hyb):
        """
        Finds the alignment of a readout hybridization to the reference
        :param hyb: integer number of the hybridization to align
        :return:
        """
        print('Aligning hyb', hyb)
        self.ro_hyb = copy(self.ro.loc[hyb])
        self.ro_hyb['matches'] = 0

        #clear the dictionaries of matched dots and edges before aligning each hybridization
        self.matchDict = {}
        self.edgeMatches = set()
        self.ref_dots = copy(self.ref)
        self.ref_dots.set_index(['row', 'col', 'z'], drop=False, inplace=True)
        self.reserved_ref_dots = pd.DataFrame(columns=self.ref.columns)
        self.ro_lat_tree = KDTree(self.ro_hyb.loc[:, 'row':'col'])
        self.trav_ind = 0

        # find a pair of dots in the hybridization readout that match a pair of the reference fiducial markers
        i = 0
        for edge_ind, edge in self.ref_edges.iterrows():
            #print(edge)
            #print('edge number:', edge_ind)
            if self._find_ro_matching_pair(edge):
                break
            i += 1
            if i >= self.n_longest_edges:
                print("Giving up search for matching edges")
                break
        offsets_ses = self.est_offsets()
        print('offsets and SEs:', offsets_ses)
        return offsets_ses

    def _find_ro_matching_pair(self, edge):
        """
        Given an edge in the reference graph, look for the corresponding matching edge in the dataframe
        of dots in the readout image.
        :param edge: Pandas series with data representing an edge: 'u_row', 'u_col', 'u_z', 'u_amp', 'u_row', 'u_col',
        'u_z', 'u_amp', 'row_dist', 'col_dist', 'length'
        :return:
        """
        #print('subsetting readout...')
        # ToDo: in C or Cython, subset with a single for loop checking each element for each threshold, and passing if one is failed

        # subset points in read out image to search through by excluding dots with row values that preclude
        # them from being matched to the reference edge
        # u is the dot with lower row, v is the dot with higher row
        u_r_max = 2048 - edge['row_dist'] + self.xyse
        v_r_min = edge['row_dist'] - self.xyse

        ro_u_candidates = self.ro_hyb.iloc[list(self.ro_hyb['row'] < u_r_max)]
        ro_v_candidates = self.ro_hyb.iloc[list(self.ro_hyb['row'] > v_r_min)]

        # exclude dots whose column values preclude them from being matched to the reference edge
        if float(edge['u_col']) < float(edge['v_col']):
            # u has a lower column than v, edge column distance is positive
            u_c_max = 2048 - edge['col_dist'] + self.xyse
            v_c_min = edge['col_dist'] - self.xyse

            ro_u_candidates = ro_u_candidates.iloc[list(ro_u_candidates['col'] < u_c_max)]
            ro_v_candidates = ro_v_candidates.iloc[list(ro_v_candidates['col'] > v_c_min)]

        else:
            # u has a larger column than v, column distance is negative
            u_c_min = -edge['col_dist'] - self.xyse
            v_c_max = 2048 + edge['col_dist'] + self.xyse

            ro_u_candidates = ro_u_candidates.iloc[list(ro_u_candidates['col'] > u_c_min)]
            ro_v_candidates = ro_v_candidates.iloc[list(ro_v_candidates['col'] < v_c_max)]

        # subset by brightness
        ucands_amp_gt_lbnd = np.greater(np.array(ro_u_candidates['amp']), float(edge['u_amp']) * self.min_bright_prop)
        ucands_amp_lt_ubnd = np.less(np.array(ro_u_candidates['amp']), float(edge['u_amp']) * self.max_bright_prop)
        ro_u_candidates = ro_u_candidates.loc[list(np.logical_and(ucands_amp_gt_lbnd, ucands_amp_lt_ubnd))]

        v_cands_amp_gt_lbnd = np.greater(np.array(ro_v_candidates['amp']), float(edge['v_amp']) * self.min_bright_prop)
        v_cands_amp_lt_ubnd = np.less(np.array(ro_v_candidates['amp']), float(edge['v_amp']) * self.max_bright_prop)
        ro_v_candidates = ro_v_candidates.loc[list(np.logical_and(v_cands_amp_gt_lbnd, v_cands_amp_lt_ubnd))]

        ro_u_candidates.sort_values(by='row', inplace=True)
        ro_v_candidates.sort_values(by='row', inplace=True)

        min_v = 0
        # print('n ucands:', len(u_cands))
        # print('n vcands:', len(v_cands))
        # keep track of closest matching edge in case none match within error

        # when reference graph is sufficiently traversed, set true to break
        done = False
        #print('n u candidates:', len(ro_u_candidates), 'n v candidates:', len(ro_v_candidates))

        #print('searching readout...')
        # search through candidates for reference edge
        for i, udot in ro_u_candidates.iterrows():
            # print(minV, 'of', len(v_cands))
            for j in range(min_v, len(ro_v_candidates)):
                rdist = ro_v_candidates.iloc[j]['row'] - udot['row']
                cdist = ro_v_candidates.iloc[j]['col'] - udot['col']
                zdist = ro_v_candidates.iloc[j]['z'] - udot['z']

                if rdist < edge['row_dist'] - self.xyse:
                    min_v += 1
                    # print('Increment minV to', minV)

                elif rdist > edge['row_dist'] + self.xyse:
                    # print('rdist to vcand out of bounds. Break')
                    break

                elif (rdist - edge['row_dist']) ** 2 + (cdist - edge['col_dist']) ** 2 < self.xyse_sq and \
                        np.abs(zdist - edge['z_dist']) < self.zse:
                    # found match! run graph traversal
                    done = self.traverse_reference(edge, udot, ro_v_candidates.iloc[j])
                    if done:
                        break
            if done:
                break
        return done

    def traverse_reference(self, matched_ref_edge, ro_u, ro_v):
        """
        From an initial matched edge dots in the reference image and the reodout image, traverse the
        reference graph in search of the rest of the fiducial markers in the readout image.
        :param matched_ref_edge: Pandas series with data representing an edge: 'u_row', 'u_col', 'u_z', 'u_amp', 'u_row', 'u_col',
        'u_z', 'u_amp', 'row_dist', 'col_dist', 'length'
        :param ro_u: Pandas series representing u dot matched to reference edge. Has parameters: row, col, z, amp
        :param ro_v: Pandas series representing v dot matched to reference edge. Has parameters: row, col, z, amp
        :return:
        """
        #print("Starting Traversal")
        self.ref_dots = copy(self.ref)
        self.ref_dots.set_index(['row', 'col', 'z'], drop=False, inplace=True)
        ref_u_in = matched_ref_edge['u_row':'u_amp']
        ref_u_in.rename({'u_row': 'row', 'u_col': 'col', 'u_z': 'z', 'u_amp': 'amp'}, inplace=True)
        ref_v_in = matched_ref_edge['v_row':'v_amp']
        ref_v_in.rename({'v_row': 'row', 'v_col': 'col', 'v_z': 'z', 'v_amp': 'amp'}, inplace=True)

        # each entry in the traversal queue is a tuple. The first entry is a pd.Series with info on the ref dot,
        # and the second entry is a pd.Series with info in the readout dot
        self.traversal_queue = [(ref_u_in, ro_u), (ref_v_in, ro_v)]
        self.n_ambiguous = 0
        self.n_trav_matched = 0
        self.total_unmatched = 0
        self.trav_ind = 0

        while self.traversal_queue:
            #print('Traversal Queue length:', len(self.traversal_queue), '; ref dot length: ', len(self.ref_dots))
            matched_dot = self.traversal_queue.pop(0)
            if tuple(matched_dot[0][:3]) in self.ref_dots.index:

                self.ref_dots.drop(tuple(matched_dot[0][:3]), inplace=True)

                self.find_ro_neighbors(matched_dot)
                #print(n_trav_matched, n_unmatched)

            
            matched_dot_key = to_tuples(matched_dot)
            if matched_dot_key in self.matchDict and self.matchDict[matched_dot_key] < self.min_edge_matches:
                self.ref_dots.append(matched_dot[0])

        n_well_matched = self._n_well_matched()
        found_enough = n_well_matched > self.min_dot_matches
            

        if not found_enough:
            # if more searching is required, add the initial dots back to the reference dataframe
            u_match_key = to_tuples((ref_u_in, ro_u))
            v_match_key = to_tuples((ref_v_in, ro_v))

            if u_match_key not in self.matchDict or self.matchDict[u_match_key] < self.min_dot_matches:
                self.ref_dots.append(ref_u_in)
            if v_match_key not in self.matchDict or self.matchDict[v_match_key] < self.min_dot_matches:
                self.ref_dots.append(ref_v_in)
            #print('Consider implementing looking through reserved ref_dots')

        print('n well matched:', n_well_matched, '; n trav matched:', self.n_trav_matched,
              '; n unmatched:', self.total_unmatched, '; n_ambiguous:', self.n_ambiguous)

        return found_enough

    def find_ro_neighbors(self, matched_dot):
        """
        For a dot matched in the reference and the readout image, search for neighbors in the readout image.
        :param matched_dot: tuple of pandas series representing dots in the reference frame and the readout frame:
            (ref_dot = pd.Series(row, col, z, amp), ro_dot = pd.Series(row, col, z, amp))
        :return: tuple of pandas series representing dots in the reference frame and the readout frame:
            (ref_dot = pd.Series(row, col, z, amp), ro_dot = pd.Series(row, col, z, amp))
        """
        self.dot_neighbor_matched = 0
        self.dot_neighbor_unmatched = 0
        matched_dot_key = to_tuples(matched_dot)

        #print('Searching through', len(self.ref_dots), 'neighbors in reference.')
        i = 0
        keep_going = True
        while i < len(self.ref_dots):
            ref_neighbor = self.trav_next()
            keep_going = self._compare_dots(matched_dot, matched_dot_key, ref_neighbor)
            if not keep_going:
                break
            i += 1
        # if made it to here without making enough matches, see if we can find enough matches in the reserved dots
        if keep_going:
            #print('looking through reserved dots.')
            for i, reserved_ref_dot in self.reserved_ref_dots.iterrows():
                keep_going = self._compare_dots(matched_dot, matched_dot_key, reserved_ref_dot, reserved=True)
                if not keep_going:
                    break

    def _compare_dots(self, matched_dot, matched_dot_key, ref_neighbor, reserved=False):

        # find position vector separating matched fiducial dot from neighbor in reference image
        rdist = ref_neighbor.row - matched_dot[0].row
        cdist = ref_neighbor.col - matched_dot[0].col
        zdist = ref_neighbor.z - matched_dot[0].z

        # search for neighbor in read out hyb
        search_row = matched_dot[1].row + rdist
        search_col = matched_dot[1].col + cdist
        search_point = [search_row, search_col]
        nn_dists, nn_inds = self.ro_lat_tree.query(search_point, 2)
        within_trav_error_inds = [nn_inds[i] for i in (0, 1) if nn_dists[i] < self.xyte]

        if within_trav_error_inds:
            within_error_inds = [i for i in within_trav_error_inds if
                                 np.abs(self.ro_hyb.iloc[i]['z'] - zdist - matched_dot[1].z) < self.zte]
            matches = self.ro_hyb.iloc[within_error_inds]
        else:
            matches = ()

        if len(matches) == 1:  # we found a match!
            matches = matches.iloc[0]
            new_match = (ref_neighbor['row':'z'], matches['row':'z'])
            if to_tuples(new_match) not in self.edgeMatches:
                self.n_trav_matched += 1
                self.dot_neighbor_matched += 1
                dropped_neighbor = self.process_match(new_match, matched_dot)
                if dropped_neighbor:
                    return True
            ref_neighbor_key = tuple(ref_neighbor['row':'z'])
            if not reserved:
                self.ref_dots.loc[ref_neighbor_key, 'n_trav_matched'] += 1

        elif len(matches) > 1:  # ambiguous, add to counter and move on
            self.n_ambiguous += 1
        else:
            self.total_unmatched += 1
            self.dot_neighbor_unmatched += 1
            ref_neighbor_key = tuple(ref_neighbor['row':'z'])
            if not reserved:
                self.ref_dots.loc[ref_neighbor_key, 'n_unmatched'] += 1
                try:
                    no_matches = np.equal(self.ref_dots.loc[ref_neighbor_key, 'n_trav_matched'], 0)[0]
                except:
                    no_matches = np.equal(self.ref_dots.loc[ref_neighbor_key, 'n_trav_matched'], 0)
                try:
                    too_many_failed_searches = np.equal(self.ref_dots.loc[ref_neighbor_key, 'n_unmatched'], self.n_umnatch_give_up)[0]
                except:
                    too_many_failed_searches = np.equal(self.ref_dots.loc[ref_neighbor_key, 'n_unmatched'],
                                                        self.n_umnatch_give_up)

                if no_matches and too_many_failed_searches:
                    #drop_ind = self.ref_dots.index.get_loc(self.ref_dots.loc[ref_neighbor_key].name)
                    drop_ind = self.ref_dots.index.get_loc(ref_neighbor_key)
                    if type(drop_ind) == slice:
                        drop_ind = drop_ind.start
                    self.ref_dots = self.ref_dots.drop(ref_neighbor_key)
                    if drop_ind < self.trav_ind:
                        self.trav_ind -= 1
                    # print('removing', ref_neighbor_key, 'too many failed searches.')
                    return True

            if self.dot_neighbor_matched == 0 and self.dot_neighbor_unmatched >= self.n_umnatch_give_up:
                # print(self.n_umnatch_give_up, 'unmatched: breaking.')
                # return n_matched, n_unmatched
                return False
        if matched_dot_key in self.matchDict and self.matchDict[matched_dot_key] >= self.min_edge_matches:
            # print('Matched enough after looking through', i, 'neighbors. Moving on.')
            # return n_matched, n_unmatched
            return False
        return True

    def trav_next(self):
        """
        Get the next fiducial dot from the reference image to look for in the readout image
        :return:
        """
        if self.trav_ind >= len(self.ref_dots):
            self.trav_ind = 0
        to_return = self.ref_dots.iloc[self.trav_ind]
        self.trav_ind += 1
        return to_return

    def process_match(self, new_match, old_match, reserved=False):
        """
        Takes dots matched by a KDTree look up to within self.xyte and within self.zte, does the necessary accounting:
        Adds the new dot the traveral queue if it has not already been added.
        If the match is within the more stringent self.xyme and self.zme, adds an edge match point to each match
        towards being included in the offset calculation
        :param new_match: tuple of pandas series representing dots in the reference frame and the readout frame for the new match:
            (ref_dot = pd.Series(row, col, z, amp), ro_dot = pd.Series(row, col, z, amp))
        :param old_match: tuple of pandas series representing dots in the reference frame and the readout frame for the \
            previously known match:
            (ref_dot = pd.Series(row, col, z, amp), ro_dot = pd.Series(row, col, z, amp))
        :param reserved:
        :return: True if drops neighbor upon neighbor reaching min_edge_matches
        """
        
        new_match_key = to_tuples(new_match)
        old_match_key = to_tuples(old_match)
        
        # if match has reached this point, is within traveral error; add to traveral queue is not
        if new_match_key not in self.trav_matches:
            self.traversal_queue.append(new_match)
            self.trav_matches |= {new_match_key}

        if self.within_error(new_match_key, old_match_key) and (new_match_key[0], old_match_key[0]) not in self.edgeMatches:
            if new_match_key in self.matchDict:
                self.matchDict[new_match_key] += 1
            else:
                self.matchDict[new_match_key] = 1
            if old_match_key in self.matchDict:
                self.matchDict[old_match_key] += 1
            else:
                self.matchDict[old_match_key] = 1
            self.edgeMatches |= {(old_match_key[0], new_match_key[0]), (new_match_key[0], old_match_key[0])}

            # remove from ref_dots dataframe if reached minimum number of matches to save time searching later
            if self.matchDict[old_match_key] >= self.min_edge_matches and old_match_key[0] in self.ref_dots.index:
                self.reserved_ref_dots = self.reserved_ref_dots.append(self.ref_dots.loc[old_match[0]])
                #drop_ind = self.ref_dots.index.get_loc(self.ref_dots.loc[old_match_key[0]].name)
                drop_ind = self.ref_dots.index.get_loc(old_match_key[0])
                if type(drop_ind) == slice:
                    drop_ind = drop_ind.start
                #print('old_match drop ind:', drop_ind)
                self.ref_dots = self.ref_dots.drop(old_match_key[0])
                if drop_ind < self.trav_ind:
                    self.trav_ind -= 1
                #print('dropped old match:', old_match_key)
            if self.matchDict[new_match_key] >= self.min_edge_matches and new_match_key[0] in self.ref_dots.index:
                self.reserved_ref_dots = self.reserved_ref_dots.append(self.ref_dots.loc[new_match[0]])
                try:
                    #drop_ind = self.ref_dots.index.get_loc(self.ref_dots.loc[new_match_key[0]].name)
                    drop_ind = self.ref_dots.index.get_loc(new_match_key[0])
                    if type(drop_ind) == slice:
                        #print('start:', drop_ind.start, 'stop:', drop_ind.stop)
                        drop_ind = drop_ind.start
                    #print('new_match drop ind:', drop_ind)
                    if drop_ind < self.trav_ind:
                        self.trav_ind -= 1
                except AttributeError as e:
                    print(e)

                self.ref_dots = self.ref_dots.drop(new_match_key[0])
                #print('dropped new match:', new_match_key)
                return True
            # reset self.trav_ind if now out of bounds
            if self.trav_ind >= len(self.ref_dots):
                self.trav_ind = 0
        return False



    def within_error(self, new_match, old_match):
        """
        Calculte the error in the edge difference the matched dots in the reference image and the readout image.
        If the error is less than that allowed by self.xyme and self.zme, then return True. Otherwise return False.
        :param new_match: tuple of pandas series representing dots in the reference frame and the readout frame for the new match:
            (ref_dot = pd.Series(row, col, z, amp), ro_dot = pd.Series(row, col, z, amp))
        :param old_match: tuple of pandas series representing dots in the reference frame and the readout frame for the \
            previously known match:
        :return:
        """
        ro_row_diff = new_match[1][0] - old_match[1][0]
        ref_row_diff = new_match[0][0] - old_match[0][0]

        ro_col_diff = new_match[1][1] - old_match[1][1]
        ref_col_diff = new_match[0][1] - old_match[0][1]

        ro_z_diff = new_match[1][2] - old_match[1][2]
        ref_z_diff = new_match[0][2] - old_match[0][2]
        xy_error_sq = (ro_row_diff - ref_row_diff) ** 2 + (ro_col_diff - ref_col_diff) ** 2
        z_error_abs = np.abs(ro_z_diff - ref_z_diff)
        #try:
        within_error = xy_error_sq < self.xyme_sq and z_error_abs < self.zme
        #except ValueError as ve:
            #print('what?')
            #raise ve
        return within_error

    def _n_well_matched(self):
        """
        Find the number of dots that have been matched by at least the minimum number of edges to be included in
        the offset calculation (self.min_edge_matches).
        :return: integer number edges meeting the threshold to be included in the offset calculation.
        """
        # Find number of well matched dots
        n_well_matched = len([match for match in self.matchDict if self.matchDict[match] >= self.min_edge_matches])
        return n_well_matched

    def est_offsets(self):
        """

        :return:
        """
        rdisp = []
        cdisp = []
        zdisp = []
        nmatches = []
        for match in self.matchDict:
            # only consider dots that were matched more than the minimum allowed number of times
            if self.matchDict[match] >= self.min_edge_matches:
                rdisp.append(match[1][0] - match[0][0])
                cdisp.append(match[1][1] - match[0][1])
                zdisp.append(match[1][2] - match[0][2])
                nmatches.append(self.matchDict[match])

        dot_r_offsets = np.array(rdisp)
        dot_c_offsets = np.array(cdisp)
        dot_z_offsets = np.array(zdisp)
        n_dot_edges = np.array(nmatches)

        # remove any outliers in XY that my have been introduced by mismatches. Don't consider z since z measurements are poor.
        outlier = True
        while outlier:
            # compute summary stats
            r_mean = np.mean(dot_r_offsets)
            c_mean = np.mean(dot_c_offsets)
            rstdv = np.std(dot_r_offsets)
            cstdv = np.std(dot_c_offsets)

            # check if there are lateral outliers
            outliers = np.logical_or(abs(dot_r_offsets - r_mean) > 2 * rstdv, abs(dot_c_offsets - c_mean) > 2 * cstdv)
            outlier = np.any(outliers)

            # remove outliers
            dot_r_offsets = dot_r_offsets[np.logical_not(outliers)]
            dot_c_offsets = dot_c_offsets[np.logical_not(outliers)]
            n_dot_edges = n_dot_edges[np.logical_not(outliers)]

        z_mean = np.mean(dot_z_offsets)
        zstdv = np.std(dot_z_offsets)

        r_mean_se = rstdv / np.sqrt(len(dot_r_offsets))
        c_mean_se = cstdv / np.sqrt(len(dot_c_offsets))
        z_mean_se = zstdv / np.sqrt(len(dot_z_offsets))

        return [r_mean, c_mean, z_mean, r_mean_se, c_mean_se, z_mean_se, len(dot_r_offsets)]

    def save_offsets(self, filename):
        self.offsets.to_csv(filename, index=False)

    def add_hyb_to_match_df(self, hyb_num, offsets):
        '''
        returns dataframe of fidicuial coordinates in reference image, matched coordinates in hyb images,
        and aligned coordinates off in hyb images. Columns indexed by hyb and then number in hyb
        :param hyb: integer current hyb number
        :return:
        Dataframe as described above
        '''
        r_offset, c_offset, z_offset, r_mean_se, c_mean_se, z_mean_se, n_dots = offsets
        ref_row = []
        ref_col = []
        ref_z = []
        comp_row = []
        comp_col = []
        comp_z = []
        for match in self.matchDict:
            if self.matchDict[match] >= self.min_edge_matches:
                ref_row.append(match[0][0])
                ref_col.append(match[0][1])
                ref_z.append(match[0][2])
                comp_row.append(match[1][0])
                comp_col.append(match[1][1])
                comp_z.append(match[1][2])

        matchesDF = pd.DataFrame()#index=index)
        matchesDF['Hyb'] = [hyb_num]*len(ref_row)
        matchesDF['ref_row'] = ref_row
        matchesDF['ref_col'] = ref_col
        matchesDF['ref_z'] = ref_z
        matchesDF['comp_row'] = comp_row
        matchesDF['comp_col'] = comp_col
        matchesDF['comp_z'] = comp_z
        matchesDF['aligned_row'] = matchesDF['comp_row'] - r_offset
        matchesDF['aligned_col'] = matchesDF['comp_col'] - c_offset
        matchesDF['aligned_z'] = matchesDF['comp_z'] - z_offset

        self.matchesDF = self.matchesDF.append(matchesDF)

        return matchesDF

    def save_matches(self, save_name):
        self.matchesDF.to_csv(save_name, index=False)

def to_tuples(match):
    """
    :param match:tuple of pandas series representing dots in the reference frame and the readout frame for the new match:
        (ref_dot = pd.Series(row, col, z, amp), ro_dot = pd.Series(row, col, z, amp))
    :return:
    """
    if type(match[0]) == pd.Series:
        ref_tup = tuple(match[0][:3])
    else:
        ref_tup = tuple(match[0].index[0])
    if type(match[1]) == pd.Series:
        ro_tup = tuple(match[1][:3])
    else:
        ro_tup = tuple(match[1].index[0])
    return (ref_tup, ro_tup)