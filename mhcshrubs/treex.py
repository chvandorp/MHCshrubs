"""@package treex
Module for drawing NetworkX graphs in a tree layout
"""
from __future__ import (print_function, division)
from builtins import zip
import networkx as nx
import numpy as np
import matplotlib as plt
import operator

def angleToAlign(angle):
    """
    Function for aligning text. Transform an angle to text alignment properties

    DEPRECATED: this issue is solved by passing the kwarg rotation_mode="anchor"

    Args:
        angle (float): the desired angle of the rotated text

    Returns:
        A dictionary of alignment codes
    """
    if angle < 0.5*np.pi:
        props = {'ha': 'left', 'va': 'bottom'}
    elif angle < np.pi:
        props = {'ha': 'right', 'va': 'bottom'}
        angle += np.pi
    elif angle < 1.5*np.pi:
        props = {'ha': 'right', 'va': 'top'}
        angle += np.pi
    else:
        props = {'ha': 'left', 'va': 'top'}
    return props

def distribute_arclengths(G, n_prev, n, a0, da, arclengths):
    """
    Auxiliary function for tree_layout. Determines the directions of the edges.
    """
    d = nx.degree(G, n)
    a1 = a0
    if d > 1:
        for n_next in nx.all_neighbors(G, n):
            if n_next != n_prev:
                a1 = distribute_arclengths(G, n, n_next, a1, da, arclengths)
    else:
        a1 += da
    arclengths[(n_prev, n)] = (a0, a1)
    return a1


def compute_positions(G, n_prev, n, arclengths, positions):
    """
    Auxiliary function for tree_layout. Determines the positions,
    using pre-computed directions.
    """
    x_prev, y_prev = positions[n_prev]
    a0, a1 = arclengths[(n_prev, n)]
    direction = (a0 + a1)/2.0
    weight = G[n_prev][n]['weight']
    x = x_prev + weight * np.cos(direction)
    y = y_prev + weight * np.sin(direction)
    positions[n] = np.array([x, y])
    if G.degree(n) > 1:
        for n_next in nx.all_neighbors(G, n):
            if n_next != n_prev:
                compute_positions(G, n, n_next, arclengths, positions)


def find_middle_edge(G):
    """
    Find an edge in the middle of a tree.
    """
    leafs = [n for n in G if G.degree(n) == 1]
    paths = [nx.shortest_path(G, l1, l2, weight='weight') for l1 in leafs for l2 in leafs]
    path_lengths = [nx.shortest_path_length(G, l1, l2, weight='weight') for l1 in leafs for l2 in leafs]
    max_index, max_value = max(enumerate(path_lengths), key=operator.itemgetter(1))
    longest_path = paths[max_index]
    cumul_edge_weight = 0.0
    middle_edge = None
    for edge in zip(longest_path[:-1], longest_path[1:]):
        cumul_edge_weight += G[edge[0]][edge[1]]['weight']
        if cumul_edge_weight > max_value/2.0:
            middle_edge = edge
            break
    return middle_edge


def tree_layout(G):
    """
    Draw an unrooted tree layout. Use the weights as edge length.

    Args:
        G: a networkx graph. If G is not a tree, then an exception is raised

    Returns:
        a dictionary of node positions
    """
    pos = {}
    ## the graph must be a tree
    if not nx.is_tree(G): raise nx.NetworkXUnfeasible()
    ## identify leafs and arc length unit
    leafs = [n for n in G if G.degree(n) == 1]
    da = 1.0/len(leafs) * 2.0 * np.pi
    ## determine arclengths, positions
    n0 = leafs[0]
    n1 = next(nx.all_neighbors(G, n0)) ## a leaf only has one neighbor
    arclengths = {}
    positions = {n0 : np.array([0.0, 0.0])}
    distribute_arclengths(G, n0, n1, da, da, arclengths) ## first da assigned to leaf n0
    compute_positions(G, n0, n1, arclengths, positions)
    return positions

def rotated_leaf_labels(ax, G, positions, labels, color="black", padding=3, fontsize=10):
    """
    Put labels at the right position and angle, so that they are aligned with the leafs.

    Args:
        ax: a matplotlib.pyplot axis object
        G: a networkx Graph
        positions: a dictionary of positions of the nodes
        labels: a dictionary of labels (str) for the nodes

    Kwargs:
        color: either a string (str) or a dictionary (dict). If a string is given,
            all labels get the same color (default: "black"), otherwise,
            for each node a different color can be given.
        padding: space between node and label
        fontsize: the size of the text
    """
    for n1, label in labels.items():
        x1, y1 = positions[n1]
        ## find neighbor
        n0 = next(nx.all_neighbors(G, n1)) ## assume that there is one neighbor
        x0, y0 = positions[n0]
        h = (x1-x0) / G[n0][n1]["weight"]
        h = -1.0 if h < -1.0 else 1.0 if h > 1.0 else h ## h might be a bit "off"
        angle = np.arccos(h)
        if y1 < y0: angle = 2.0*np.pi - angle
        space = " " * padding ## a number of spaces
        col = color[n1] if isinstance(color, dict) else color
        ax.text(x1, y1, space + label, rotation=angle/(2*np.pi)*360.0,
                rotation_mode="anchor", va="center", fontsize=fontsize, color=col)
