"""@package colortrees
Functions for rendering and coloring ETE and NetworkX trees

The coloring is based on attributes attached to the nodes/branches
of the tree.
"""

from __future__ import (print_function, division)
from builtins import (str, zip)
import sys ## for finding python version
import ete3
import ete3.treeview
import networkx as nx
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import rgb2hex
from mhcshrubs import auxiliary as aux
from mhcshrubs import (mhctools, treex)



def mkEteHlaTree(newick_str, midpoint=False, colorfun=(lambda x: 'black')):
    """
    Transform a Newick tree into the ETE format.

    Args:
        newick_str (str): a string representing the tree in the Newick format

    Kwargs:
        midpoint (bool): if True, find the midpoint of the tree and re-root.
        colorfun (function): a rule for coloring faces (TODO)

    Returns:
        A tuple consisting of an ete3.Tree object and a list of integer node names
    """
    tree = ete3.Tree(newick_str)
    if midpoint:
        tree.set_outgroup(tree.get_leaves()[0]) ## make arbitrary leaf the outgroup
        outgrp_node = tree.get_midpoint_outgroup() ## now use ete algorithm
        tree.set_outgroup(outgrp_node) ## and re-define the outgroup
    ## give nodes a name to be used in Bayesian model
    nodeNames = []
    for i, node in enumerate(tree.traverse()):
        if node.name != '':
            hla = mhctools.MhcObject(node.name, fmt="Newick")
            node.add_feature('hla', str(hla)) ## FIXME: use MhcObject
        node.name = str(i)
        nodeNames.append(i)
        ## set width and color of edges
        size = 10
        nstyle = ete3.treeview.NodeStyle()
        nstyle["vt_line_width"] = size
        nstyle["hz_line_width"] = size ## FIXME make sure that ultrametric trees are rendered properly
        ## make sure that ultrametric trees are rendered correctly
        nstyle["size"] = 0 ## do not render nodes
        node.set_style(nstyle)
        ## NB: colorEteTree overrides set_style
    ## set colors of the HLA tags
    for leaf in tree:
        N = ete3.AttrFace("hla", fsize=50, fgcolor=colorfun(leaf.hla))
        leaf.add_face(N, 0)
    return (tree, nodeNames)


def getPathsFromEteTree(tree):
    paths = []
    for leaf in tree:
        ## create a path
        node = leaf
        path = []
        while node.up:
            path.append(node)
            node = node.up
        ## add paths to path
        paths += [(leaf.hla, path)]
    return paths


def addNodeFeatures(tree, nodeNames, fname, fvalues):
    for i in nodeNames:
        node = tree&str(i)
        node.add_feature(fname, fvalues[i])


def addLeafBars(tree, hlaAlleles, counts, represDict=None):
    """
    input:
    - tree -- an ete object with hlaAlleles in the leafs
    - hlaAlleles -- ordered HLA alleles, indexed by locus
    - counts -- counts in the same order as hlaAlleles
    - represDict -- if some of the hlaAlleles are not in the tree,
      then they must be represented by another alleles. (default: None)
    output: None, the tree is modified
    """
    ## loop over leaves of the tree
    if represDict is not None:
        equivClassDict = aux.invertRepresDict(represDict)
    for leaf in tree:
        ## find the index of the HLA allele
        mhc = mhctools.MhcObject(leaf.hla)
        if represDict is not None:
            ## note that the represDict can contain more alleles than listed in hlaAlleles
            mhcs = [x for x in equivClassDict[mhc] if x in hlaAlleles[mhc.locus]]
        else:
            mhcs = [mhc]
        total_count = 0
        for mhc in mhcs:
            hlaIdx = hlaAlleles[mhc.locus].index(mhc)
            total_count += counts[mhc.locus][hlaIdx]
        ## add a bar to the leaf
        F = ete3.RectFace(width=total_count*2.5, height=10, fgcolor='k', bgcolor='k')
        leaf.add_face(F, 0)


## function that lists the relevant values and nodes
def getNodesAndVals(tree, fname, fun=(lambda x: x), f_scale_dist=False):
    ## FIXME: pre or post function application?
    fvalues = [fun(node.__dict__[fname]) * (node.dist if f_scale_dist else 1.0)
               for node in tree.traverse()
               if fname in node.features]
    fnodeNames = [node.name for node in tree.traverse()
                  if fname in node.features]
    return (fvalues, fnodeNames)


def colorEteTree(tree, cfname=None, wfname=None, cffun=(lambda x: x), wffun=(lambda x: x),
                 cf_scale_dist=False, wf_scale_dist=False, center=0.0, cax=None, wax=None):
    """
    Color the edges of an ETE tree, using attributes of the nodes.

    The function returns nothing, but the tree in the argument is modified.

    Args:
        tree (ete3.Tree): an ETE tree

    Kwargs:
        cfname (str): the name of the feature that is used for color of the branch
        wfname (str): the name of the feature that is used for the width of the branch
        cffun (function):  a transformation of the color feature
        wffun (function): a transformation of the width feature
        cf_scale_dist (bool): boolean indicating that the color feature should
            be scaled with the branch length
        wf_scale_dist (bool): boolean indicating that the width feature should
            be scaled with the branch length
        center (float): defines the midpoint of the colormap
        cax: matplotlib.pyplot axis object for a color bar (ignored if None)
        wax: matplotlib.pyplot axis object for a width scale (ignore if None)

    Returns:
        ranges for the color scale and width scale
    """
    ## a colormap
    cmap = get_cmap("coolwarm") ## TODO: parameter
    cmin = 0; cmax = 0; wmin = 0; wmax = 0
    defaultcol = 'k'
    defaultsize = 10; sizescale = 40
    ## collect fvalues and determine the scale
    cscaledvals = {}
    wscaledvals = {}
    if cfname is not None:
        cfvalues, cfnodeNames = getNodesAndVals(tree, cfname, cffun, f_scale_dist=cf_scale_dist)
        if len(cfnodeNames) > 0:
            cfvalues = [x-center for x in cfvalues]
            maxval = np.max(cfvalues)
            minval = np.min(cfvalues)
            rangeval = max(abs(maxval), abs(minval))
            cmax = rangeval
            cmin = -rangeval
            scaledvals = [(x+rangeval) / (2*rangeval) for x in cfvalues]
            cscaledvals = dict(zip(cfnodeNames, scaledvals))
    if wfname is not None:
        wfvalues, wfnodeNames = getNodesAndVals(tree, wfname, wffun, f_scale_dist=wf_scale_dist)
        if len(wfnodeNames) > 0:
            maxval = np.max(wfvalues)
            minval = np.min(wfvalues)
            wmax = maxval
            wmin = minval
            if maxval > minval:
                scaledvals = [(x - minval)/(maxval - minval) for x in wfvalues]
                wscaledvals = dict(zip(wfnodeNames, scaledvals))
    ## add color to nodes and change width
    for node in tree.traverse():
        i = node.name
        nstyle = ete3.treeview.NodeStyle()
        nstyle["size"] = 0 ## do not render nodes
        nstyle["vt_line_color"] = defaultcol
        nstyle["hz_line_color"] = defaultcol
        nstyle["vt_line_width"] = defaultsize
        nstyle["hz_line_width"] = defaultsize
        if i in list(cscaledvals.keys()):
            col = cmap(cscaledvals[i])
            nstyle["vt_line_color"] = rgb2hex(col)
            nstyle["hz_line_color"] = rgb2hex(col)
        if i in list(wscaledvals.keys()):
            hsize = sizescale * wscaledvals[i] + defaultsize
            vsize = defaultsize
            nstyle["vt_line_width"] = vsize
            nstyle["hz_line_width"] = hsize
        node.set_style(nstyle)
    if cax is not None:
        norm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)
        cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')
        tick_locator = mpl.ticker.MaxNLocator(nbins=6)
        cb.locator = tick_locator
        cb.update_ticks()
        cax.spines['top'].set_visible(False)
        cax.spines['right'].set_visible(False)
        cax.spines['bottom'].set_visible(True)
        cax.spines['left'].set_visible(False)
        cax.set_xlabel("$\\Sigma\\alpha$")
    if wax is not None:
        xs = [wmin, wmax]
        ys = [sizescale * x + defaultsize for x in [0,1]]
        tops = [0.5*y for y in ys]
        bottoms = [-0.5*y for y in ys]
        wax.fill_between(xs, bottoms, tops, color='k', lw=0)
        yrange = np.max(tops) - np.min(bottoms)
        wax.set_ylim(np.min(bottoms) - 0.2 * yrange, np.max(tops) + 0.2 * yrange)
        wax.spines['top'].set_visible(False)
        wax.spines['right'].set_visible(False)
        wax.spines['bottom'].set_visible(True)
        wax.spines['left'].set_visible(False)
        wax.set_yticks([])
        wax.set_xlabel("$|\\alpha|$")



def addCumulativeFeature(node, fname, cfname=None, cfval=None, ffun=(lambda x: x),
                         scale_dist=False, dfun=(lambda x: x)):
    """
    Compute the cumulative value of a trait along the tree.

    The function is called recursively until a leaf is found

    Todo:
        Skip edges that don't have the feature?

    Args:
        node (ete3.Node): the current node
        fname (str): the name of the feature, used to get the value
        cfname (str): the name of the cumulative feature.
            If None, '_cumulative' is added to fname (default: None)
        cfval: the cumulative value so far
        ffun (function): an optinal transformation applied before adding (default: identity)
        scale_dist (bool): a boolean indicating that the cumulative value
            should be scaled with the branch length
        dfun (function): if scale_dist is True, then the branch length (node.dist)
            is transformed with dfun (default: identity)
    """
    if cfname is None:
        cfname = fname + "_cumulative"
    scaled_fval = ffun(node.__dict__[fname]) * dfun(node.dist) if scale_dist else 1.0
    ## add the cumulative feature to the node
    if cfval is None:
        cfval_next = scaled_fval
    else:
        cfval_next = cfval + scaled_fval
    node.add_feature(cfname, cfval_next)
    ## call this function for node's children
    for child in node.children:
        addCumulativeFeature(child, fname, cfname=cfname, cfval=cfval_next, ffun=ffun,
                             scale_dist=scale_dist, dfun=dfun)


def getTreeStyle():
    """
    Define an acceptible TreeStyle for rendering

    Returns:
        A ete3 TreeStyle object.
    """
    ts = ete3.treeview.TreeStyle()
    ts.show_branch_length = False
    ts.mode = "c"
    ts.show_leaf_name = False
    ts.optimal_scale_level = "full"
    ts.root_opening_factor = 0.25
    return ts


def mkColorBar(filename): ## TODO: add a colorbar to a figure
    # Make a figure and axes with dimensions as desired.
    fig = plt.figure(figsize=(4, 0.3))
    ax1 = fig.add_subplot(111)
    # Set the colormap and norm to correspond to the data for which
    # the colorbar will be used.
    cmap = mpl.cm.get_cmap("coolwarm")
    norm = mpl.colors.Normalize(vmin=-1, vmax=1)
    # ColorbarBase derives from ScalarMappable and puts a colorbar
    # in a specified axes, so it has everything needed for a
    # standalone colorbar.  There are many more kwargs, but the
    # following gives a basic continuous colorbar with ticks
    # and labels.
    ticks = [-0.75, 0, 0.75]
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm,
                                    orientation='horizontal', ticks=ticks)
    ax1.set_xticklabels(["protective", "neutral", "deleterious"])
    fig.savefig(filename, dpi=200, bbox_inches="tight")



## NetworkX-related functions


def eteToNetworkX(tree):
    """
    Transform an ETE MHC tree into a NetworkX graph

    The colors and widths of the branches are copied to the NetworkX tree.
    The resulting tree can then be rendered with renderTreeNetworkX

    Args:
        tree (ete3.Tree): an ETE tree object

    Returns:
        The same tree as a networkx.Graph object
    """
    ## construct the graph
    G = nx.Graph()
    edge_colors = {}
    edge_widths = {}
    node_labels = {}
    for node in tree.traverse():
        G.add_node(node.name)
        if node.up:
            G.add_edge(node.name, node.up.name, weight=node.dist)
            edge = (node.name, node.up.name)
            edge_color = node.img_style["hz_line_color"]
            edge_width = node.img_style["hz_line_width"]
            edge_colors[edge] = edge_color
            edge_widths[edge] = edge_width
        if node.is_leaf():
            node_label = node.hla
            node_labels[node.name] = node_label
    ## add edge color and width attributes
    nx.set_edge_attributes(G, edge_colors, name="color")
    ## requires networkx high enough verion
    nx.set_edge_attributes(G, edge_widths, name="width")
    nx.set_node_attributes(G, node_labels, name="label")
    return G


def renderTreeNetworkX(G, filename, colorfun=(lambda x: "black")):
    fig, (ax) = plt.subplots(1, 1, figsize=(7.5,7.5))
    edge_colors = nx.get_edge_attributes(G, "color")
    edge_width_scale = 0.1
    edge_widths = nx.get_edge_attributes(G, "width")
    node_labels = nx.get_node_attributes(G, "label")
    node_label_colors = {n : colorfun(x) for n, x in node_labels.items()}
    positions = treex.tree_layout(G)
    nx.draw_networkx(G, pos=positions, ax=ax, node_size=0, with_labels=False,
                     edge_color=[edge_colors[e] for e in G.edges()],
                     width=[edge_widths[e]*edge_width_scale for e in G.edges()])
    treex.rotated_leaf_labels(ax, G, positions, node_labels, color=node_label_colors)
    ax.axis('off')
    fig.savefig(filename, dpi=200, bbox_inches='tight')
    plt.close(fig) ## prevent warnings when running in parallel
