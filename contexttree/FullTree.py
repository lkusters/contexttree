# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:05:40 2016

@author: Lieneke Kusters
"""

import copy
from contexttree.TreeCounts import TreeCounts

ALPHABET = "ACGT"
base2bin = {'A': 0b0, 'C': 0b1, 'G': 0b10, 'T': 0b11}
bin2base = ['A', 'C', 'G', 'T']


class FullTree(TreeCounts):
    """ Inherits from TreeCounts, adds functionality of rates
    probabilities and rates calculation
    """

    _symbollogprobs = []
    _rself = None

    def getprobs(self):
        """ Use the counts to compute the corresponding probabilities of the
        symbols in log2-space and return (also store achievable compression
        rate, rself)
        """

        if not (self._rself is None):
            # the value is already there
            return self._symbollogprobs

        if self._sequencelength == 0:
            self._rself = None
            return []

        symbollogprobs = []
        rself = 0
        for idx in range(len(self._symbolcounts)):
            counts = self._symbolcounts[idx]
            logprobs = self._counts2logprobs(counts)
            symbollogprobs.append(logprobs)
            rself -= sum([a*b for a, b in zip(counts, logprobs)])
        self._rself = rself/self._sequencelength
        self._symbollogprobs = symbollogprobs
        return symbollogprobs

    def getrself(self):
        """ return rself or calculate rself if not calculated yet"""

        if self._rself is None:
            self.getprobs()

        return self._rself

    def getratetree(self, tree):
        """ estimate the achievable compression rate when applying this
        tree model for compression of the sequence corresponding to the input
        contexttree
        """

        self._verifysametype(tree)

        if tree._sequencelength == 0:  # the tree is empty
            return 0
        if self._sequencelength == 0:  # the model is empty
            return 2

        symbolprobs = self.getprobs()

        rate = 0
        for counts, probs in zip(tree._symbolcounts, symbolprobs):
            rate -= sum([c*p for c, p in zip(counts, probs)])

        return rate/tree._sequencelength

    def getratesequence(self, sequence):
        """ apply the model to a sequence and return the list of corresponding
        rates corresponding to each symbol in the sequence
        """
        raise ValueError(
                         "Sorry this functionality has not been " +
                         "implemented yet")

    def getdivergence(self, tree):
        """ This function returns the estimated KL-Divergence of the
        probability distribution of the own tree in comparison to other tree

        We use D(q_z||p_x) ~ lim_{n-> inf} 1/n log2(q_z(Z)/p_x(Z))
         = Rother - Rself  (for sequence Z)

        Here 'tree' corresponds to the (model of the) input sequence Z
        Here q_z is the probability distribution of this tree and p_x
        corresponds to the other tree
        """

        rother = self.getratetree(tree)
        rself = tree.getrself()

        divergence = rother-rself
        return divergence

    def getdistance(self, tree):
        """ Use divergence as a distance metric, by estimating it in
        both directions """

        div1 = self.getdivergence(tree)
        div2 = tree.getdivergence(self)

        return (div1+div2)/2

    def getcopy(self):
        """ Make a copy of this tree and return it
        """

        tree = FullTree(self._maximumdepth)
        for attr in vars(self):
            setattr(tree, attr, copy.deepcopy(getattr(self, attr)))

        return tree
