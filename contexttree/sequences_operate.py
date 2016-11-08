# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:52:50 2016

@author: Lieneke Kusters

functions that use FullTree class to operate on DNA sequences
"""

from contexttree.FullTree import FullTree
from Bio.Seq import Seq


def seqdistance(depth, seq1, seq2, rev):
    """
    Compare 2 DNA sequences, by constructing the trees and computing distance,
    based on KL-Divergence
    rev = True means that reverse complement of the sequences should be
    included in the models
    """

    tree1 = FullTree(depth, seq1)
    tree2 = FullTree(depth, seq2)
    if rev:
        tree1.updatesymbolcounts(str(Seq(seq1).reverse_complement()))
        tree2.updatesymbolcounts(str(Seq(seq2).reverse_complement()))

    return tree1.getdistance(tree2)


def seqsmodel(depth, seqgenlist, rev):
    """
    Construct the full tree model, corresponding to a generator list of
    sequences. (note that the list should be generator)
    rev = True means that reverse complement of the sequences should be
    included in the models
    """

    seq = next(seqgenlist)
    tree = FullTree(depth, seq)
    if rev:
        tree.updatesymbolcounts(str(Seq(seq).reverse_complement()))
        for seq in seqgenlist:
            tree.updatesymbolcounts(seq)
            tree.updatesymbolcounts(str(Seq(seq).reverse_complement()))
    else:
        for seq in seqgenlist:
            tree.updatesymbolcounts(seq)
    return tree


def modelapply(model, seq):
    """
    Apply the full tree model to a sequence and get resulting compression rate.
    """

    return model.getratetree(FullTree(model.getdepth(), seq))
