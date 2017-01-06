# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:05:40 2016

@author: Lieneke Kusters
"""

import copy
from contexttree.TreeCounts import TreeCounts
import warnings

ALPHABET = "ACGT"
base2bin = {'A': 0b0, 'C': 0b1, 'G': 0b10, 'T': 0b11}

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
        rself = 0
        if not (self._rself is None):
            # the value is already there
            return self._symbollogprobs

        if self._sequencelength == 0:
            self._rself = None
            return dict()

        symbollogprobs = dict()
        for key, counts in self._symbolcounts.items():
            logprobs = self._counts2logprobs(counts)
            symbollogprobs[key] = logprobs
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

        if len(tree._symbolcounts) == 0:  # the tree is empty
            return 0

        symbolprobs = self.getprobs()

        rate = 0
        for key, val in tree._symbolcounts.items():
            if key in symbolprobs:
                rate -= sum([a*b for a, b in zip(val, symbolprobs[key])])
            else:
                rate -= sum([a*-2 for a in val])
        return rate/tree._sequencelength

    def getratesequence(self, sequence):
        """ apply the model to a sequence and return the list of corresponding
        rates corresponding to each symbol in the sequence
        """
        raise ValueError(
                         "Sorry this functionality has not been " +
                         "implemented yet")

        # first verify if input is valid
        if len(sequence) <= self._maximumdepth:
            warnings.warn("sequence length {0}, is too short\nsequence: {1}"
                          .format(str(len(sequence)), sequence)
                          )
            sequences = []
        elif not(self._verifyinputsequence(sequence)):
            sequences = filter(lambda s: len(s) > self._maximumdepth,
                               sequence.split('N'))
        else:
            sequences = [sequence]

        # Now initialize the rate
        rate = 0
        seqlen = 0
        # Prepare conversion table
        keys = dict(zip(ALPHABET, range(len(ALPHABET))))
        symbolprobs = self.getprobs()

        for sequence in sequences:
            sequence = sequence.upper()  # upper case

            # Special case, tree of depth 0
            if self._maximumdepth == 0:
                rate -= sum([sequence.count(ALPHABET[index]) *
                             symbolprobs[ALPHABET[index]] for
                            index in range(4)])
                seqlen += len(sequence)

            else:
                initcontext = ''
                for i in range(self._maximumdepth):
                    initcontext += sequence[i]
                sequence = sequence[self._maximumdepth:]
                initcontext = initcontext[::-1]

                # we start with initial context initcontext
                context = initcontext
                # now each next state is just a shift
                for symbol in sequence:
                    if context in symbolprobs:
                        rate += symbolprobs[context][keys[symbol]]
                        seqlen += 1
                    else:
                        rate += 1/4  # default value
                        seqlen += 1

                    context = symbol+context[:-1]
        if seqlen > 0:
            return rate/seqlen
        else:
            return 0

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
