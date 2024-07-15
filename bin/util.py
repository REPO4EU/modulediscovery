import graph_tool.all as gt
from pathlib import Path


def load_graph(path):
    """
    Load a graph-tool graph from a file. The file format is determined by the file extension.
    """
    extension = Path(path).suffix
    if extension in [".gt", ".graphml", ".xml", ".dot", ".gml"]:
        return gt.load_graph(path)
    else:
        return gt.load_graph_from_csv(path)


def read_seeds(path):
    """
    Loads a list of seeds from a file containing one line per seed gene.
    """
    with open(path, "r") as file:
        seeds = [line.strip() for line in file.readlines() if line.strip()]
    return seeds


def name2index(g):
    """
    Create a mapping from gene name to vertex index.
    """
    index2name = g.vertex_properties["name"]
    return {index2name[v]: v for v in g.iter_vertices()}
