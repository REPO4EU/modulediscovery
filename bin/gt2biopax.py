#!/usr/bin/env python
import argparse
import logging
import sys
from pathlib import Path
import requests

from nedrex.core import iter_edges
from nedrex.core import get_nodes

import graph_tool.all as gt
from pybiopax import biopax, model_to_owl_file

import nedrex

open_url = "http://10.162.163.34:7123/open"
nedrex.config.set_url_base(open_url)
from nedrex.core import api_keys_active, get_api_key

if api_keys_active():
    api_key = get_api_key(accept_eula=True)
    nedrex.config.set_api_key(api_key)

logger = logging.getLogger()


def get_uniprot_from_entrez(entrez_ids: list[str]) -> dict[str, list[str]]:

    # returns the ids without the prefixes uniprot. / entrez.
    response = requests.post(
        "https://api.nedrex.net/open/relations/get_encoded_proteins",
        json={"nodes": entrez_ids},
    )
    response.raise_for_status()
    return response.json()


def get_proteins(uniprot_ids: list[str]) -> any:
    proteins = list({id: "" for id in uniprot_ids})
    result = getNodeDict(proteins, "protein")
    return result


def get_genes_to_proteins(uniprot_ids):
    edges = getEdges("protein_encoded_by_gene", target_domain_ids=uniprot_ids)
    prot2gene = {}
    for uniprot_id in uniprot_ids:
        prot2gene[uniprot_id] = ""
    genes = []
    for e in edges:
        if e["sourceDomainId"] in prot2gene:
            genes.append(e["targetDomainId"])
            prot2gene[e["sourceDomainId"]] = e["targetDomainId"]
    return genes, prot2gene


def getEdges(
    type: str, source_domain_ids=[], target_domain_ids=[], extra_attributes=[]
):
    attributes = ["sourceDomainId", "targetDomainId", "dataSources", "type"]
    attributes.extend(extra_attributes)
    body = {
        "source_domain_id": source_domain_ids,
        "target_domain_id": target_domain_ids,
        "attributes": attributes,
    }
    try:
        response = requests.post(url=f"{open_url}/{type}/attributes/json", json=body)
        response.raise_for_status()
        data = response.json()
    except requests.exceptions.RequestException as e:
        print(f"HTTP Anfrage fehlgeschlagen: {e}")
        return None
    return data


def getNodeDict(ids, node_type, extra_attributes=[]):
    attributes = ["primaryDomainId", "type", "displayName"]
    attributes.extend(extra_attributes)
    body = {
        "node_ids": ids,
        "attributes": attributes,
    }
    try:
        response = requests.post(
            url=f"{open_url}/{node_type}/attributes/json", json=body
        )
        response.raise_for_status()
        data = response.json()
    except requests.exceptions.RequestException as e:
        print(f"HTTP Anfrage fehlgeschlagen: {e}")
        return None
    result = {res["primaryDomainId"]: res for res in data}
    return result


def get_nedrex_data(entrez_ids: list[str], uniprot_ids=None, protein2gene=None) -> any:

    proteins = []

    gene2prot = {}

    if uniprot_ids is None:
        # in gene2prot: without prefix uniprot. / entrez.
        gene2prot = get_uniprot_from_entrez(entrez_ids)
        uniprot_ids = [
            "uniprot." + protein
            for proteins_list in gene2prot.values()
            for protein in proteins_list
        ]
        uniprot_ids = list(set(uniprot_ids))
        # dict with encoding gene for each protein
        protein2gene = {
            "uniprot." + protein: "entrez." + gen
            for gen, proteins in gene2prot.items()
            for protein in proteins
        }

    # get all needed protein nodes as dict
    proteins = get_proteins(uniprot_ids)

    # get all needed gene nodes + create dict for faster access
    genes = getNodeDict(entrez_ids, "gene")

    # list for all needed edges
    edges = []

    edges_new = getEdges(
        "gene_associated_with_disorder", source_domain_ids=list(genes.keys())
    )
    edges.extend(edges_new)

    nodes_to_get = list({e["targetDomainId"] for e in edges_new})
    disorders = getNodeDict(nodes_to_get, "disorder")

    # get drugs that target a protein (encoded by our genes)
    edges_new = getEdges("drug_has_target", target_domain_ids=list(proteins.keys()))
    edges.extend(edges_new)

    nodes_to_get = list({e["sourceDomainId"] for e in edges_new})
    drugs = getNodeDict(nodes_to_get, "drug")

    # get edges for go annotations; no nodes needed
    edges_new = getEdges(
        "protein_has_go_annotation",
        source_domain_ids=list(proteins.keys()),
        extra_attributes=["qualifiers"],
    )
    edges_relevant = []
    for e in edges_new:
        if "is_active_in" in e["qualifiers"]:
            edges_relevant.append(e)

    edges.extend(edges_relevant)

    # get side effects for drugs
    edges_new = getEdges("drug_has_side_effect", source_domain_ids=list(drugs.keys()))
    edges.extend(edges_new)

    nodes_to_get = list({e["targetDomainId"] for e in edges_new})
    # get side effect nodes for drugs
    sideeffects = getNodeDict(nodes_to_get, "side_effect")

    edges_new = getEdges(
        "drug_has_indication",
        source_domain_ids=list(drugs.keys()),
        target_domain_ids=list(disorders.keys()),
    )
    edges.extend(edges_new)

    edges_new = getEdges(
        "drug_has_contraindication",
        source_domain_ids=list(drugs.keys()),
        target_domain_ids=list(disorders.keys()),
    )
    edges.extend(edges_new)

    return (
        genes,
        disorders,
        edges,
        drugs,
        proteins,
        protein2gene,
        gene2prot,
        sideeffects,
    )


# switched mapping decides on the direction of the mapping: True -> target to source, False -> source to target
# type refers to the type of the edge
def create_dict_mapping(edges, type, switched_mapping=False):
    mapping = {}
    for edge in edges:
        if edge["type"] == type:
            source_domain_id = edge["sourceDomainId"]
            target_id = edge["targetDomainId"]
            if switched_mapping:
                mapping.setdefault(target_id, []).append(
                    {"id": source_domain_id, "dataSources": edge["dataSources"]}
                )
            else:
                mapping.setdefault(source_domain_id, []).append(
                    {"id": target_id, "dataSources": edge["dataSources"]}
                )
    return mapping


class BioPAXFactory:
    def __init__(self, input_path: Path, id_space: str = "entrez"):
        self.input_path = input_path
        self.id_space = id_space
        self.g = None
        self.biopaxmodel = None
        self.xRefs: dict[str, biopax.Xref] = {}
        self.entityRefs: dict[str, biopax.EntityReference] = {
            "gene_associated_with_disorder.vocab": biopax.UnificationXref(
                uid="gene_associated_with_disorder.XREF",
                db="PSI-MI",
                id="MI:0361",
                comment="gene associated with disorder",
            ),
            "drug_has_side_effect.vocab": biopax.UnificationXref(
                uid="drug_has_side_effect.XREF",
                db="PSI-MI",
                id="MI:0361",
                comment="drug has side effect",
            ),
            "drug_has_indication.vocab": biopax.UnificationXref(
                uid="drug_has_indication.XREF",
                db="PSI-MI",
                id="MI:0361",
                comment="drug has indication",
            ),
            "drug_has_contraindication.vocab": biopax.UnificationXref(
                uid="drug_has_contraindication.XREF",
                db="PSI-MI",
                id="MI:0361",
                comment="drug has contraindication",
            ),
            "gene_product.vocab": biopax.UnificationXref(
                uid="gene_product.XREF",
                db="PSI-MI",
                id="MI:0251",
                comment="gene product",
            ),
            "cellular_component.vocab": biopax.UnificationXref(
                uid="cellular_component.XREF",
                db="PSI-MI",
                id="MI:0354",
                comment="cellular component",
            ),
        }
        self.entities: dict[str, biopax.BioPaxObject] = {}
        self.edgeTypes: dict[str, biopax.RelationshipTypeVocabulary] = {
            "gene_associated_with_disorder": biopax.RelationshipTypeVocabulary(
                term=["additional information"],
                comment="gene associated with disorder",
                uid="gene_associated_with_disorder.vocab",
                xref=self.entityRefs["gene_associated_with_disorder.vocab"],
            ),
            # "variant_associated_with_disorder": biopax.RelationshipTypeVocabulary(term = ["variant associated with disorder"], uid = "variant_associated_with_disorder.XREF"),
            # "variant_affects_gene": biopax.RelationshipTypeVocabulary(term = ["variant affects gene"], uid = "variant_affects_gene.XREF"),
            "gene_product": biopax.RelationshipTypeVocabulary(
                term=["gene product"],
                uid="gene_product.vocab",
                xref=self.entityRefs["gene_product.vocab"],
            ),
            "cellular_component": biopax.RelationshipTypeVocabulary(
                term=["cellular component"],
                uid="cellular_component.vocab",
                xref=self.entityRefs["cellular_component.vocab"],
            ),
            "drug_has_indication": biopax.RelationshipTypeVocabulary(
                term=["drug has indication"], uid="drug_has_indication.XREF"
            ),
            "drug_has_contraindication": biopax.RelationshipTypeVocabulary(
                term=["drug has contraindication"], uid="drug_has_contraindication.XREF"
            ),
            "drug_has_side_effect": biopax.RelationshipTypeVocabulary(
                term=["additional information"],
                comment="drug has side effect",
                uid="drug_has_side_effect.vocab",
                xref=self.entityRefs["drug_has_side_effect.vocab"],
            ),
        }
        self.organism: dict[str, biopax.BioSource] = {
            "human": biopax.BioSource(
                uid="human", display_name="Homo sapiens", standard_name="human"
            )
        }

    def get_owl_path(self) -> Path:
        return self.input_path.with_suffix(".owl")

    def load_graph(self):
        self.g = gt.load_graph(str(self.input_path))
        logger.debug(f"{self.g=}")

    def create_biopax_model(self):
        if not self.g:
            self.load_graph()
        if self.id_space == "entrez":
            # TODO: entrez prefix will already be added in the future
            entrez_ids = []
            for v_i in self.g.get_vertices():
                präfix = (
                    ""
                    if str(self.g.vp["name"][v_i]).startswith("entrez.")
                    else "entrez."
                )
                entrez_ids.append(präfix + self.g.vp["name"][v_i])
            self.add_info(entrez_ids)
        elif self.id_space == "uniprot":
            # TODO: entrez prefix will already be added in the future
            uniprot_ids = []
            for v_i in self.g.get_vertices():
                präfix = (
                    ""
                    if str(self.g.vp["name"][v_i]).startswith("uniprot.")
                    else "uniprot."
                )
                uniprot_ids.append(präfix + self.g.vp["name"][v_i])
            self.add_info(uniprot_ids, True)

        self.biopaxmodel = biopax.BioPaxModel(objects=self.get_bioax_objects())

    def add_info(self, ids, protein=False):
        if protein:
            entrez_ids, prot2gene = get_genes_to_proteins(ids)
            (
                genes,
                disorders,
                edges,
                drugs,
                proteins,
                protein2gene,
                gene2prot,
                sideeffects,
            ) = get_nedrex_data(entrez_ids, ids, prot2gene)
        else:
            (
                genes,
                disorders,
                edges,
                drugs,
                proteins,
                protein2gene,
                gene2prot,
                sideeffects,
            ) = get_nedrex_data(ids)

        # TODO: add disorders + sideeffects as soon as Biopax supports it

        protein2go = create_dict_mapping(edges, "ProteinHasGOAnnotation")

        if protein:
            self.add_protein_info(proteins, protein2gene, protein2go, True)
        else:
            self.add_protein_info(proteins, protein2gene, protein2go, False, gene2prot)

        gene2disorder = create_dict_mapping(edges, "GeneAssociatedWithDisorder")
        drug2disorderIndication = create_dict_mapping(edges, "DrugHasIndication")
        drug2disorderContraindication = create_dict_mapping(
            edges, "DrugHasContraindication"
        )
        protein2drug = create_dict_mapping(edges, "DrugHasTarget", True)
        drug2sideeffect = create_dict_mapping(edges, "DrugHasSideEffect")

        self.add_drug_info(
            protein2drug,
            drugs,
            drug2sideeffect,
            drug2disorderIndication,
            drug2disorderContraindication,
        )
        self.add_gene_info(gene2disorder, genes)

    def add_gene_info(self, gene2disorder, genes):
        for entrez_id, gene in genes.items():
            if entrez_id in gene2disorder:
                self.add_gene(entrez_id, gene, gene2disorder[entrez_id])
            else:
                self.add_gene(entrez_id, gene)

    def add_drug_info(
        self,
        protein2drug,
        drug_nodes,
        drug2sideeffect,
        drug2disorderIndication,
        drug2disorderContraindication,
    ):
        for p, drugs in protein2drug.items():
            for drug in drugs:
                id = drug["id"]
                drug_node = drug_nodes[id]
                sides = []
                indications = []
                contraindications = []
                if id in drug2sideeffect:
                    sides = drug2sideeffect[id]
                if id in drug2disorderIndication:
                    indications = drug2disorderIndication[id]
                if id in drug2disorderContraindication:
                    contraindications = drug2disorderContraindication[id]
                self.add_drug(
                    drug,
                    p,
                    drug_node["displayName"],
                    sides,
                    indications,
                    contraindications,
                )

    def add_drug(
        self,
        drug,
        uniprot_id,
        display_name,
        sideeffects,
        indications,
        contraindications,
    ):
        drug_id = drug["id"]
        # split drug_id by "." into db and id
        this_db, this_id = drug_id.split(".")
        uniXRef = [
            self.xRefs.setdefault(
                drug_id,
                biopax.UnificationXref(
                    uid=f"{drug_id}.XREF",
                    db=this_db,
                    id=this_id,
                ),
            )
        ]

        for sideeffect in sideeffects:
            id = sideeffect["id"]
            uniXRef.append(
                self.xRefs.setdefault(
                    id,
                    biopax.RelationshipXref(
                        uid=f"{id}.XREF",
                        db="sider",
                        id=id,
                        comment=sideeffect["dataSources"],
                        relationship_type=self.edgeTypes["drug_has_side_effect"],
                    ),
                )
            )
        for indication in indications:
            id = indication["id"]
            uniXRef.append(
                self.xRefs.setdefault(
                    id,
                    biopax.RelationshipXref(
                        uid=f"{id}.XREF",
                        db="drugbank",
                        id=id,
                        comment=indication["dataSources"],
                        relationship_type=self.edgeTypes["drug_has_indication"],
                    ),
                )
            )
        for contraindication in contraindications:
            id = contraindication["id"]
            uniXRef.append(
                self.xRefs.setdefault(
                    id,
                    biopax.RelationshipXref(
                        uid=f"{id}.XREF",
                        db="drugbank",
                        id=id,
                        comment=contraindication["dataSources"],
                        relationship_type=self.edgeTypes["drug_has_contraindication"],
                    ),
                )
            )
        entityRef = self.entityRefs.setdefault(
            drug_id,
            biopax.SmallMoleculeReference(
                uid=f"{drug_id}.REF", xref=uniXRef, display_name=display_name
            ),
        )
        self.entities[drug_id] = biopax.SmallMolecule(
            uid=drug_id, entity_reference=entityRef, display_name=display_name
        )

        if uniprot_id:
            uniprot_id = uniprot_id.lstrip("uniprot.")
            self.add_drug_protein_interaction(uniprot_id, drug_id)

    def add_protein_info(
        self, proteins, protein2gene, protein2go, uniprot_ids=True, gene2prot=None
    ):
        for uniprot_id, protein in proteins.items():
            encoding_gene = protein2gene[uniprot_id]
            go_to_protein = []
            if uniprot_id in protein2go:
                go_to_protein = protein2go[uniprot_id]
            if protein:
                self.add_protein(uniprot_id, encoding_gene, go_to_protein, protein)
            else:
                self.add_protein(uniprot_id, encoding_gene, go_to_protein)

        if uniprot_ids:
            for e_i, e_j in self.g.get_edges():
                uniprot_id1 = self.g.vp["name"][e_i].lstrip("uniprot.")
                uniprot_id2 = self.g.vp["name"][e_j].lstrip("uniprot.")
                self.add_PPI(uniprot_id1, uniprot_id2)
        else:
            for e_i, e_j in self.g.get_edges():
                entrez_id1 = self.g.vp["name"][e_i].lstrip("entrez.")
                entrez_id2 = self.g.vp["name"][e_j].lstrip("entrez.")
                # gene2prot is not annotated with the id prefixes uniprot. / entrez.
                uniprot_ids1 = gene2prot[entrez_id1]
                uniprot_ids2 = gene2prot[entrez_id2]
                for uniprot_id1 in uniprot_ids1:
                    for uniprot_id2 in uniprot_ids2:
                        self.add_PPI(uniprot_id1, uniprot_id2)

    def write(self):
        if not self.biopaxmodel:
            self.create_biopax_model()
        model_to_owl_file(self.biopaxmodel, self.get_owl_path())

    def add_protein(self, protein_id, gene_id, go_s, protein=None):
        # in the biopax file: id-prefix uniprot./entrez. should be removed
        gene_id = gene_id.lstrip("entrez.")
        uniprot_id = protein_id.lstrip("uniprot.")

        uniXRef = [
            self.xRefs.setdefault(
                uniprot_id,
                biopax.UnificationXref(
                    uid=f"{uniprot_id}.XREF", db="uniprot", id=uniprot_id
                ),
            )
        ]
        if gene_id:
            uniXRef.append(
                self.xRefs.setdefault(
                    gene_id,
                    biopax.RelationshipXref(
                        uid=f"{gene_id}.XREF",
                        db="NCBI GENE",
                        id=gene_id,
                        relationship_type=self.edgeTypes["gene_product"],
                    ),
                )
            )
        for go in go_s:
            id = go["id"].replace("go.", "GO:")
            uniXRef.append(
                self.xRefs.setdefault(
                    id,
                    biopax.RelationshipXref(
                        uid=id,
                        db=go["dataSources"],
                        id=id,
                        relationship_type=self.edgeTypes["cellular_component"],
                    ),
                )
            )
        if protein:
            displayName = [protein["displayName"]]
        else:
            displayName = []

        entityRef = self.entityRefs.setdefault(
            uniprot_id,
            biopax.ProteinReference(
                uid=f"{uniprot_id}.REF",
                xref=uniXRef,
                display_name=displayName,
                organism=self.organism["human"],
            ),
        )
        self.entities[uniprot_id] = biopax.Protein(
            uid=uniprot_id, entity_reference=entityRef, display_name=displayName
        )

    def get_bioax_objects(self):
        return (
            list(self.xRefs.values())
            + list(self.entityRefs.values())
            + list(self.entities.values())
            + list(self.edgeTypes.values())
            + list(self.organism.values())
        )

    def add_PPI(self, uniprot_id1, uniprot_id2):
        interaction_id = f"{uniprot_id1}_{uniprot_id2}"
        self.entities[interaction_id] = biopax.MolecularInteraction(
            uid=interaction_id,
            participant=[self.entities[uniprot_id1], self.entities[uniprot_id2]],
            display_name=[f"{uniprot_id1} {uniprot_id2}"],
        )

    def add_drug_protein_interaction(self, uniprot_id, drug_id):
        interaction_id = f"{drug_id}_{uniprot_id}"
        self.entities[interaction_id] = biopax.MolecularInteraction(
            uid=interaction_id,
            participant=[self.entities[drug_id], self.entities[uniprot_id]],
            display_name=[f"{drug_id} {uniprot_id}"],
            comment="drug has target",
        )

    def add_gene(self, entrez_id, gene, associated_disorders=None):
        # in the biopax file: id-prefix entrez. should be removed
        entrez_id = entrez_id.lstrip("entrez.")
        uniXRef = [
            self.xRefs.setdefault(
                entrez_id,
                biopax.UnificationXref(
                    uid=f"{entrez_id}.XREF", db="NCBI GENE", id=entrez_id
                ),
            )
        ]
        if associated_disorders is not None:
            for disorder in associated_disorders:
                id = disorder["id"]
                uniXRef.append(
                    self.xRefs.setdefault(
                        id,
                        biopax.RelationshipXref(
                            uid=f"{id}.XREF",
                            db="MONDO",
                            # remove mondo. prefix from id
                            id=id.replace("mondo.", ""),
                            comment=disorder["dataSources"],
                            relationship_type=self.edgeTypes[
                                "gene_associated_with_disorder"
                            ],
                        ),
                    )
                )
        self.entities[entrez_id] = biopax.Gene(
            uid=entrez_id,
            xref=uniXRef,
            organism=None,
            display_name=[gene["displayName"]],
        )


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

    parser.add_argument(
        "-i",
        "--idspace",
        help="ID space of the given network.",
        type=str,
        choices=["entrez", "uniprot"],
        default="entrez",
    )

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
