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
        self._symbolcounts = None
        self._sequencelength = 0
        self._maximumdepth = depth
        self._rself = None  # achievable compression rate of full source tree

        if sequence is not None:
            self._countsymbols(sequence)

    def __str__(self):
        return "tree: {2} of depth {0}, initcontext {1}".format(
           self._maximumdepth, self._initialcontext, type(self)) +\
           "\n symbolcounts: \n {0}".format(str(self._symbolcounts))

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
        counts: (dict) keys are occuring contexts (tuple), counts are
                symbol counts for symbols of alphabet given context
        """

        # first verify if input is valid
        if not(self._verifyinputsequence(sequence)):
            # warnings.warn("For now, we remove invalid symbol, N, however, "
            #              "in the future, we want to split the sequences "
            #              "around the invalid characters, before model "
            #              "construction in order to prevent invalid contexts."
            #              )
            # sequence = sequence.replace('N', '')
            sequences = filter(None, sequence.split('N'))
        else:
            sequences = [sequence]

        for sequence in sequences:
            if len(sequence) <= self._maximumdepth:
                # we need a sequence of at least length > self._maximumdepth
                warnings.warn("sequence length {0}, is too short\nsequence: {1}"
                              .format(str(len(sequence)), sequence)
                              )
                continue

            # Now prepare the data
            if self._symbolcounts is None:
                counts = dict()
            else:
                counts = self._symbolcounts

            keys = dict(zip(ALPHABET, range(len(ALPHABET))))  # Conversion table
            sequence = sequence.upper()  # upper case

            # Special case, tree of depth 0
            if self._maximumdepth == 0:
                counts[''] = [sequence.count(ALPHABET[index]) for
                              index in range(4)]
                return counts

            initcontext = ''
            for i in range(self._maximumdepth):
                initcontext += sequence[i]
            sequence = sequence[self._maximumdepth:]
            initcontext = initcontext[::-1]

            # we start with initial context initcontext
            context = initcontext
            # now each next state is just a shift
            for symbol in sequence:
                if context in counts:
                    counts[context][keys[symbol]] += 1
                else:
                    # context has not occured before, so initialize new context
                    # with 0 counts
                    counts[context] = [0 for sym in range(len(ALPHABET))]
                    counts[context][keys[symbol]] += 1
                context = symbol+context[:-1]
            self._initialcontext += [initcontext]
            self._sequencelength += len(sequence)

            self._symbolcounts = counts
            self._rself = None  # just became invalid in case it was set before

    # Functions for updating the counts in the tree or combining trees
    def updatesymbolcounts(self, sequence):
        """ update symbol counts: call __countsymbols()

        Keyword arguments:
        sequence:   (str) The sequence for which we count the symbolcounts
        """

        if self._symbolcounts is None:
            warnings.warn("cannot update, since no symbols were counted " +
                          "yet, initializing with current input sequence " +
                          "instead")

        self._countsymbols(sequence)

    def combine(self, tree):
        """ add the input tree to this tree
        """

        self._verifytreedephts(tree)

        if self._symbolcounts is None and tree._symbolcounts is None:
            warnings.warn("combining with an empty tree")
            # do nothing
        elif tree._symbolcounts is None:
            warnings.warn("combining with an empty tree")
            # do nothing
        elif self._symbolcounts is None:
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
            self._rself = None  # just became invalid and should be updated

    def getcopy(self):
        """ Make a copy of this tree and return it
        """

        tree = TreeCounts(self._maximumdepth)
        for attr in vars(self):
            setattr(tree, attr, copy.deepcopy(getattr(self, attr)))

        return tree
    def getdepth(self):
        
        return self._maximumdepth
