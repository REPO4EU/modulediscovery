import graph_tool.all as gt
import pandas as pd
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


def vp2df(g):
    """Convert the vertex properties of a graph to a pandas DataFrame. 'name' will be used as index."""
    # Get all vertex properties
    vertex_props = g.vertex_properties

    # Prepare a dictionary to hold the data for the DataFrame
    data = {
        prop_name: [prop[v] for v in g.vertices()]
        for prop_name, prop in vertex_props.items()
    }

    # Create and return the DataFrame
    df = pd.DataFrame(data)
    df.set_index("name", inplace=True)
    return df


def ep2df(g):
    """Convert the edge properties of a graph to a pandas DataFrame. 'source' and 'target' will be used as index."""
    # Get all edge properties
    edge_props = g.edge_properties

    # Prepare a dictionary to hold the data for the DataFrame
    data = {
        prop_name: [prop[e] for e in g.edges()]
        for prop_name, prop in edge_props.items()
    }
    source_list = []
    target_list = []
    for e in g.iter_edges():
        source_list.append(g.vp["name"][e[0]])
        target_list.append(g.vp["name"][e[1]])

    data["source"] = source_list
    data["target"] = target_list

    # Create and return the DataFrame
    df = pd.DataFrame(data)
    df.set_index(["source", "target"], inplace=True)
    return df
