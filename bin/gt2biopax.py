#!/usr/bin/env python

import argparse
import logging
import sys
from pathlib import Path
from subprocess import call

import graph_tool.all as gt
from pybiopax import biopax, model_to_owl_file

logger = logging.getLogger()


def get_uniprot_from_entrez(entrez_ids: list[str]) -> dict[str, list[str]]:
    import requests
    response = requests.post('https://api.nedrex.net/open/relations/get_encoded_proteins', json={'nodes': entrez_ids})
    response.raise_for_status()
    return response.json()


class BioPAXFactory():
    def __init__(self, input_path: Path, id_space: str = 'entrez'):
        self.input_path = input_path
        self.id_space = id_space
        self.g = None
        self.biopaxmodel = None
        self.xRefs: dict[str, biopax.Xref] = {}
        self.entityRefs: dict[str, biopax.EntityReference] = {}
        self.entities: dict[str, biopax.BioPaxObject] = {}

    def get_owl_path(self) -> Path:
        return self.input_path.with_suffix('.owl')

    def load_graph(self):
        self.g = gt.load_graph(str(self.input_path))
        logger.debug(f"{self.g=}")

    def create_biopax_model(self):
        if not self.g:
            self.load_graph()
        if self.id_space == 'entrez':
            # entrez_ids = [f"entrez.{self.g.vp['name'][v_i]}" for v_i in self.g.get_vertices()]
            # gene2prot = get_uniprot_from_entrez(entrez_ids)
            # for k, v in gene2prot.items():
            #     [self.add_protein(v_i) for v_i in v]
            #     self.add_gene(k)
            self.add_gene_info()
        elif self.id_space == 'uniprot':
            self.add_protein_info()
        self.biopaxmodel = biopax.BioPaxModel(objects=self.get_bioax_objects())

    def add_gene_info(self):
        for v_i in self.g.get_vertices():
            entrez_id = self.g.vp['name'][v_i]
            self.add_gene(entrez_id)
        for e_i, e_j in self.g.get_edges():
            entrez_id1 = self.g.vp['name'][e_i]
            entrez_id2 = self.g.vp['name'][e_j]
            self.add_GGI(entrez_id1, entrez_id2)

    def add_protein_info(self):
        for v_i in self.g.get_vertices():
            uniprot_id = self.g.vp['name'][v_i]
            self.add_protein(uniprot_id)
        for e_i, e_j in self.g.get_edges():
            uniprot_id1 = self.g.vp['name'][e_i]
            uniprot_id2 = self.g.vp['name'][e_j]
            self.add_PPI(uniprot_id1, uniprot_id2)

    def write(self):
        if not self.biopaxmodel:
            self.create_biopax_model()
        model_to_owl_file(self.biopaxmodel, self.get_owl_path())
        # if validate:
        #     self.validate(self.get_owl_path())

    # def validate(self, owl_path: Path):
    #     call(
    #         ["java", "-jar", "paxtools-5.3.0.jar", "validate", owl_path.name,
    #          f"{owl_path.stem}_validation.html", "notstrict"])

    def add_protein(self, uniprot_id):
        # TODO: add all the other properties

        # uniprot_id = self.g.vp['name'][v]
        uniXRef = self.xRefs.setdefault(uniprot_id,
                                        biopax.UnificationXref(uid=f"{uniprot_id}.XREF", db='uniprot', id=uniprot_id))
        entityRef = self.entityRefs.setdefault(uniprot_id,
                                               biopax.ProteinReference(uid=f"{uniprot_id}.REF", xref=uniXRef))
        self.entities[uniprot_id] = biopax.Protein(uid=uniprot_id, entity_reference=entityRef)

    def get_bioax_objects(self):
        return list(self.xRefs.values()) + list(self.entityRefs.values()) + list(self.entities.values())

    def add_PPI(self, uniprot_id1, uniprot_id2):
        # uniprot_id1 = self.g.vp['name'][e_i]
        # uniprot_id2 = self.g.vp['name'][e_j]
        interaction_id = f"{uniprot_id1}_{uniprot_id2}"
        self.entities[interaction_id] = biopax.MolecularInteraction(uid=interaction_id,
                                                                    participant=[self.entities[uniprot_id1],
                                                                                 self.entities[
                                                                                     uniprot_id2]])

    def add_gene(self, entrez_id):
        uniXRef = self.xRefs.setdefault(entrez_id,
                                        biopax.UnificationXref(uid=f"{entrez_id}.XREF", db='ncbi gene', id=entrez_id))
        self.entities[entrez_id] = biopax.Gene(uid=entrez_id, xref=uniXRef, organism=None)

    def add_GGI(self, entrez_id1, entrez_id2):
        interaction_id = f"{entrez_id1}_{entrez_id2}"
        self.entities[interaction_id] = biopax.GeneticInteraction(uid=interaction_id,
                                                                  participant=[self.entities[entrez_id1],
                                                                               self.entities[
                                                                                   entrez_id2]])


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse network files to different formats.",
        epilog="Example: python gt2biopax.py network.gt --namespace entrez",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Input network.",
    )
    # parser.add_argument(
    #     "-f",
    #     "--format",
    #     help="Output format (default gt).",
    #     choices=("gt","diamond", "domino", "robust"),
    #     default="gt",
    # )
    parser.add_argument(
        "-i",
        "--idspace",
        help="ID space of the given network.",
        type=str,
        choices=['entrez', 'uniprot'],
        default="entrez",
    )
    # parser.add_argument(
    #     "--no_validate",
    #     help="The resulting biopax file is not validated.",
    #     action="store_false",
    # )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    logger.debug(f"{args=}")
    biopax = BioPAXFactory(args.file_in, args.idspace)
    biopax.write()


if __name__ == "__main__":
    sys.exit(main())
