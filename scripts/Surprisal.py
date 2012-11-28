#!/usr/bin/env python
"""
SurprisalCalculator objects are used to calculate the surprisal
values for two sparse transition count matrices. The variance of 
the surprisal values can also be calculated.
"""

import sys
import os 
import numpy as np
from scipy.io import mmread

class SurprisalCalculator:

    def __init__(self, sparse1, sparse2):
        self.sparse1 = sparse1
        self.sparse2 = sparse2

    def prepare_sparse_matrices(self, i):
        """
        Prepares count arrays from sparse matrices within a  
        SurprisalCalculator object for surprisal calculations
        for every state i. 
        """
        sparse1 = self.sparse1
        sparse2 = self.sparse2
    
        j1 = sparse1.col[sparse1.row == i] 
        j2 = sparse2.col[sparse2.row == i] 
        j_states = np.unique(np.append(j1, j2))
        j_states.sort()

        counts1 = sparse1.data[sparse1.row == i]
        counts2 = sparse2.data[sparse2.row == i] 

        for index in range(len(j_states)):
            if not j1.__contains__(j_states[index]):  
                counts1 = np.insert(counts1, index, 0)     
            if not j2.__contains__(j_states[index]):  
                counts2 = np.insert(counts2, index, 0)
        

        return counts1.astype(np.float64), counts2.astype(np.float64)

    def calculate_entropy(self, counts):
        """Calculates the entropy of a given set of transition counts"""
        entropy_sum = 0.0
        total_counts = np.sum(counts).astype(np.float64)

        for index in range(len(counts)):
            if counts[index] != 0:
                entropy_sum += counts[index]*np.log(counts[index])

        if total_counts == 0:
            entropy = 0.0
        else:
            entropy = (1.0/total_counts)*(total_counts*np.log(total_counts) - entropy_sum)
             
        return entropy

    def calculate_surprisal(self, counts1, counts2, weighted=False):
        """
        Calculates surprisal for any given set of counts. If weighted, 
        the surprisal value is divided by the total number of counts. 
        """
        combined_counts = counts1 + counts2
        total_combined_counts = np.sum(combined_counts)
        total_counts1 = np.sum(counts1)
        total_counts2 = np.sum(counts2)

        combined_entropy = self.calculate_entropy(combined_counts)
        entropy1 = self.calculate_entropy(counts1)
        entropy2 = self.calculate_entropy(counts2)

        surprisal = ((total_combined_counts * combined_entropy) -
                        (total_counts2 * entropy2) - (total_counts1 * entropy1))
        if (weighted):
            surprisal /= float(total_combined_counts)

        return surprisal
 
    def calculate_all_surprisal(self, verbose=False, weighted=False):
        """
        Calculates surprisal for all states of the sparse matrices
        in the SurprisalCalculator object. Returns in array where
        array[i] is the surprisal of state i. 
        """
        surprisals = []
        for i in range(self.sparse1.shape[0]):
            if (verbose):
                print("Working on state %d"%i) 

            counts1, counts2 = self.prepare_sparse_matrices(i) 

            if(not weighted):
                surprisals.append(self.calculate_surprisal(counts1, counts2))
            if(weighted):
                surprisals.append(self.calculate_surprisal(counts1, counts2, True))
        return surprisals

    def get_covariance_matrix(self, sparse):
        """
        Returns the covariance matrix of a multinomial distribution for 
        p_i ~ counts
        """
        covariance_matrix = np.zeros((sparse.shape[0], sparse.shape[1]))
        
        for i in range(sparse.shape[0]):
            j_states = sparse.col[sparse.row == i]
            counts = sparse.data[sparse.row == i].astype(np.float64)
            total_counts = np.sum(counts).astype(np.float64)
            self_transition_counts = sparse.data[(sparse.row == i)*(sparse.col == i)][0]

            for index in range(len(j_states)):
                j_state = j_states[index]
                covariance_matrix[i][j_state] = -self_transition_counts*counts[index]/total_counts
                if i == j_state:
                    covariance_matrix[i][i] += self_transition_counts
        return covariance_matrix

    def estimate_surprisal_variance_analytical(self, istate, weighted=True):
        """
        Analytical estimation of variance. Still not quite working, use
        bootstrap method instead. 
        """
        covariance1 = self.get_covariance_matrix(self.sparse1)
        covariance2 = self.get_covariance_matrix(self.sparse2)
        size = self.sparse1.shape[1] #Corresponds to number of j-states
                                     
        diagonal_covariances = np.zeros((2*size, 2*size), np.float)
        diagonal_covariances[0:size, 0:size] = covariance1
        diagonal_covariances[size:2*size, size:2*size] = covariance2
        sensitivities1 = np.zeros(size)
        sensitivities2 = np.zeros(size)

        j_states1 = self.sparse1.col[self.sparse1.row == istate]
        j_states2 = self.sparse2.col[self.sparse2.row == istate]
        j_states = np.unique(np.append(j_states1, j_states2))
        counts1, counts2 = self.prepare_sparse_matrices(istate)
        total_counts1 = np.sum(counts1).astype(np.float64)
        total_counts2 = np.sum(counts2).astype(np.float64)
        total_counts = total_counts1 + total_counts2

        for i in range(len(j_states)):
            j_state = j_states[i]
            if counts1[i] != 0:
                sensitivities1[j_state] = (np.log(total_counts1 / total_counts) - 
                                          np.log(counts1[i] / (counts1[i] + counts2[i])))

            if counts2[i] != 0:
                sensitivities2[j_state] = (np.log(total_counts2 / total_counts) -
                                          np.log(counts2[i] / (counts1[i] + counts2[i])))

            if weighted:
                sensitivities1[j_state] /= total_counts
                sensitivities2[j_state] /= total_counts

            sensitivities = np.append(sensitivities1, sensitivities2)
            #print 'sensitivities.shape', sensitivities.shape
            #print 'diagonal_covariances.shape', diagonal_covariances.shape
            variance = (sensitivities.T).dot(diagonal_covariances).dot(sensitivities)

        return variance

    def estimate_surprisal_variance_bootstrap(self, counts1, counts2, n_bootstraps=100, weighted=False):
        """
        Calculates the variance of two sets of transition counts by 
        creating pseudo-count matrices from transition probabilities
        and then estimating the variance of the resulting array.    
        """
        total_counts1 = np.sum(counts1)
        prob1 = np.divide(counts1, total_counts1)

        total_counts2 = np.sum(counts2)
        prob2 = np.divide(counts2, total_counts2)

        surprisals = []
        sampled_counts1 = np.random.multinomial(total_counts1, prob1,  
                                                size = n_bootstraps)
        sampled_counts2 = np.random.multinomial(total_counts2, prob2,
                                                size = n_bootstraps) 
        for trial in range(n_bootstraps):
            surprisal = self.calculate_surprisal(sampled_counts1[trial,:],
                                                 sampled_counts2[trial,:], weighted=weighted)
            surprisals.append(surprisal)
        return np.array(surprisals).var()

if __name__ == '__main__':
    sparse1 = mmread(sys.argv[1])
    sparse2 = mmread(sys.argv[2])
    obj = SurprisalCalculator(sparse1, sparse2) 
    print 'state\tsurprisal\tVariances (analytical)\tVariances (bootstrap)\tN_row_counts'
    for istate in range(72):
        counts1, counts2 = obj.prepare_sparse_matrices(istate)
        surprisal = obj.calculate_surprisal(counts1, counts2, weighted=True)
        var_analytical = obj.estimate_surprisal_variance_analytical(istate)
        var_bootstrap  = obj.estimate_surprisal_variance_bootstrap(counts1, counts2, n_bootstraps=1000, weighted=True)
        num_row_counts = np.sum(counts1) + np.sum(counts2)
        print '%d\t%e\t%e\t%e\t%d'%(istate, surprisal, var_analytical, var_bootstrap, num_row_counts)
    #print "Weighted Surprisals", obj.calculate_all_surprisal(weighted=True)                
