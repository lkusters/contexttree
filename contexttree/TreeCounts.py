# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 10:19:12 2016

@author: Lieneke Kusters
"""

from __future__ import division
import numpy as np
import warnings
import copy

ALPHABET = "ACGT"
base2bin = {'A': 0b0, 'C': 0b1, 'G': 0b10, 'T': 0b11}
bin2base = ['A', 'C', 'G', 'T']
# TODO Rself support!


class TreeCounts:
    """ Parent class for all the tree objects, the class object can be
    initialized with a depth and an input string for which the counts are
    stored. This class stores a dictionary which represents the leafs and
    corresponding counts of the context tree structure.
    """

    def __init__(self, depth, sequence=None):
        # Usually the tree should be initialized with a sequence,
        # however we support initialization without sequence
        # for calling from self.copy()
        if not isinstance(depth, int):
            raise ValueError("Invalid maximum depth", depth)

        self._initialcontext = []
        self._symbolcounts = [[0, 0, 0, 0] for idx in range(4**depth)]
        self._sequencelength = 0
        self._maximumdepth = depth

        if sequence is not None:
            self._countsymbols(sequence)

    def __str__(self):
        depth = self._maximumdepth
        return "tree: {2} of depth {0}, initcontext {1}".format(
           depth, self._initialcontext, type(self)) +\
           "\n symbolcounts: \n" + str(
            [''.join([bin2base[idx & 4**d] for d in range(depth)]) + ' ' +
             str(self._symbolcounts[idx])+'\n' for idx in range(4**depth)])

    # Verification of input sequence / tree
    def _verifyinputsequence(self, sequence):
        """ verify that the sequence has only symbols in the alphabet
        """
        # alphabet in sequence
        no_valid_symbols = 0  # number of valid symbols
        for symbol in ALPHABET:
            no_valid_symbols += sequence.count(symbol)
        if not no_valid_symbols == len(sequence):
            # invalid symbols are in the sequence
            countN = sequence.count('N')
            if -no_valid_symbols + len(sequence) == countN:
                warnings.warn(
                    "Sequence has values that are not in alphabet ({0}): "
                    "{1} N's were found".format(ALPHABET, str(countN))
                )
                return False
            else:  # it was not N, so there is no 'support' for this symbol
                raise ValueError(
                            "Sequence has values that are not in alphabet"
                            " ({0}): {1} valid of total {2} symbols \n"
                            "Also, only {3} are 'N' \nsequence: {4}"
                            .format(ALPHABET, str(no_valid_symbols),
                                    str(len(sequence)),
                                    str(sequence.count('N')),
                                    str(sequence)
                                    )
                             )
        else:
            return True  # valid sequence: no invalid symbols found

    def _verifytreedephts(self, tree):
        """ verify that the input tree has the same depth as the source
        """
        if not tree._maximumdepth == self._maximumdepth:
            raise ValueError("trees cannot interact, as their depth is " +
                             "not matching ",
                             self._maximumdepth, tree._maximumdepth)

    def _verifysametype(self, tree):
        """ verify that both trees are same type and depth"""
        self._verifytreedephts(tree)
        if not(type(tree) == type(self)):
            raise ValueError("both trees should be of same type", type(self),
                             type(tree))

    # Helper functions
    def _counts2logprobs(self, counts):
        """ convert symbolcounts to probabilities
        """
        denum = np.log2(sum(counts)+len(ALPHABET)/2)
        logprobs = [np.log2(c+1/2) - denum for c in counts]
        return logprobs

    # Initialize the tree, by counting symbol occurences in the sequence
    def _countsymbols(self, sequence):
        """ Count the symbol occurences for a context tree, for the contexts at
            maximum depth

        Keyword arguments:
        sequence:   (str) The sequence for which we count the symbolcounts

        Finds and stores:
        counts: (list) idxs correspond to occuring contexts, counts are
                symbol counts for symbols of alphabet given context
        """

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

        mask = 4**depth - 1
        for sequence in sequences:
            sequence = sequence.upper()  # upper case
            depth = 2 #self._maximumdepth
            _symbolcounts = [[0, 0, 0, 0] for i in range(4**depth)]
            sequence = 'CGTAA'
            initcontext = 0b0
            for i in range(depth):
                initcontext = initcontext << 2
                initcontext += base2bin[sequence[i]]

            sequence = sequence[depth:]
            for symbol in sequence:
                _symbolcounts[initcontext][base2bin[symbol]] += 1
                initcontext = (initcontext << 2) & mask
                initcontext += base2bin[symbol]
            print(_symbolcounts)

            self._sequencelength += len(sequence)

        del sequences

    # Functions for updating the counts in the tree or combining trees
    def updatesymbolcounts(self, sequence):
        """ update symbol counts: call __countsymbols()

        Keyword arguments:
        sequence:   (str) The sequence for which we count the symbolcounts
        """

        if len(self._symbolcounts)==0:
            warnings.warn("cannot update, since no symbols were counted " +
                          "yet, initializing with current input sequence " +
                          "instead")

        self._countsymbols(sequence)

    def combine(self, tree):
        """ add the input tree to this tree
        """

        self._verifytreedephts(tree)

        if len(self._symbolcounts)==0 and len(tree._symbolcounts)==0:
            warnings.warn("combining with an empty tree")
            # do nothing
        elif len(tree._symbolcounts)==0:
            warnings.warn("combining with an empty tree")
            # do nothing
        elif len(self._symbolcounts)==0:
            warnings.warn("combining with an empty tree")
            # copy tree to self
            for attr in vars(tree):
                setattr(self, attr, getattr(tree, attr))
        else:
            for key, val in tree._symbolcounts.items():
                if key in self._symbolcounts:
                    self._symbolcounts[key] = [a+b for (a, b) in
                                               zip(self._symbolcounts[key],
                                                   val)]
                else:
                    self._symbolcounts[key] = val
            self._sequencelength += tree._sequencelength
            self._initialcontext += [tree._initialcontext]

    def getcopy(self):
        """ Make a copy of this tree and return it
        """

        tree = TreeCounts(self._maximumdepth)
        for attr in vars(self):
            setattr(tree, attr, copy.deepcopy(getattr(self, attr)))

        return tree

    def getdepth(self):
        return self._maximumdepth
