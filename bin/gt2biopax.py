#!/usr/bin/env python
import json
import argparse
import logging
import sys
from pathlib import Path

from nedrex.core import iter_edges
from nedrex.core import get_nodes

import graph_tool.all as gt
from pybiopax import biopax, model_to_owl_file

import nedrex
open_url = "https://api.nedrex.net/open/"
nedrex.config.set_url_base(open_url)
from nedrex.core import api_keys_active, get_api_key

if api_keys_active():
    api_key = get_api_key(accept_eula=True)
    nedrex.config.set_api_key(api_key)

logger = logging.getLogger()


def get_uniprot_from_entrez(entrez_ids: list[str]) -> dict[str, list[str]]:
    import requests

    response = requests.post("https://api.nedrex.net/open/relations/get_encoded_proteins", json={"nodes": entrez_ids})
    response.raise_for_status()
    return response.json()


def get_nedrex_data(entrez_ids: list[str], gene2prot: dict[str, list[str]]) -> any:
    genes = get_nodes(node_type = "gene", node_ids=entrez_ids)

    genes = {gene['primaryDomainId']: gene for gene in genes}

    #file_path = 'relevant_edges.json'

    edge_types_to_get = ['gene_associated_with_disorder'] #'variant_affects_gene'
    edges = [e for edge_type in edge_types_to_get for e in iter_edges(edge_type)]

    edges_relevant = []
    for e in edges:
        if e["sourceDomainId"] in genes or e["targetDomainId"] in genes:
            edges_relevant.append(e)
    
    edges = edges_relevant

    # variants_to_get = set()
    # for e in edges:
    #     if e["type"] == "VariantAffectsGene":
    #         variants_to_get.add(e["sourceDomainId"])
    # variants_to_get = list(variants_to_get)

    # batch_size = 300  # Die Anzahl der IDs in jeder Gruppe

    # variants = []
    # for i in range(0, len(variants_to_get), batch_size):
    #     batch_ids = variants_to_get[i:i+batch_size]

    #     # get_nodes für die aktuelle Gruppe von IDs aufrufen
    #     variant_nodes = get_nodes(node_type="genomic_variant", node_ids=batch_ids)
    #     variants.extend(variant_nodes)
        
    # edges_new = [e for e in iter_edges("variant_associated_with_disorder")]
    
    # # create dictionary for faster access
    # variants = {variant['primaryDomainId']: variant for variant in variants}

    # edges_relevant = []
    # for e in edges_new:
    #     if e["sourceDomainId"] in variants:
    #         edges_relevant.append(e)
    
    # edges.extend(edges_relevant)

    # # Liste als JSON in die Datei schreiben
    # with open(file_path, 'w') as json_file:
    #     json.dump(edges, json_file)

    # just for testing purposes: get saved edges from file
    #with open(file_path, 'r') as json_file:
    #    edges = json.load(json_file)

    # variants_to_get = set()
    # for e in edges:
    #     if e["type"] == "VariantAffectsGene":
    #         variants_to_get.add(e["sourceDomainId"])
    # variants_to_get = list(variants_to_get)

    
    batch_size = 300  # Die Anzahl der IDs in jeder Gruppe

    variants = []
    # for i in range(0, len(variants_to_get), batch_size):
    #     batch_ids = variants_to_get[i:i+batch_size]

    #     # get_nodes für die aktuelle Gruppe von IDs aufrufen
    #     variant_nodes = get_nodes(node_type="genomic_variant", node_ids=batch_ids)
    #     variants.extend(variant_nodes)

    variants = {variant['primaryDomainId']: variant for variant in variants}

    disorders_to_get = set()

    for e in edges:
        if e["type"] == "GeneAssociatedWithDisorder": # or e["type"] == "VariantAssociatedWithDisorder":
            disorders_to_get.add(e["targetDomainId"])
    disorders_to_get = list(disorders_to_get)

    disorders = []
    for i in range(0, len(disorders_to_get), batch_size):
        batch_ids = disorders_to_get[i:i+batch_size]

        # get_nodes für die aktuelle Gruppe von IDs aufrufen
        disorder_nodes = get_nodes(node_type="disorder", node_ids=batch_ids)
        disorders.extend(disorder_nodes)

    # Dictionary erstellen, wobei die primaryDomainId als Schlüssel verwendet wird für schnellen Zugriff
    disorders = {disorder['primaryDomainId']: disorder for disorder in disorders}

    # get drugs that target our associated disorders
    # note: not used until now bc biopax does not support disorders
    edges_drugs = [e for e in iter_edges("drug_has_indication")]
    edges_relevant = []
    for e in edges_drugs:
        if e["targetDomainId"] in disorders:
            edges_relevant.append(e)
    edges_drugs = edges_relevant

    protein2gene = {"uniprot."+protein: gen for gen, proteins in gene2prot.items() for protein in proteins}

    # get drugs that target a protein encoded by our genes
    edges_drugs = [e for e in iter_edges("drug_has_target")]
    edges_relevant = []
    for e in edges_drugs:
        if e["targetDomainId"] in protein2gene:
            edges_relevant.append(e)
    edges_drugs = edges_relevant
    
    drugs_to_get = set()
    for e in edges_drugs:
        if e["type"] == "DrugHasIndication" or e["type"] == "DrugHasTarget":
            drugs_to_get.add(e["sourceDomainId"])
    drugs_to_get = list(drugs_to_get)

    drugs = []
    for i in range(0, len(drugs_to_get), batch_size):
        batch_ids = drugs_to_get[i:i+batch_size]

        # get_nodes für die aktuelle Gruppe von IDs aufrufen
        drug_nodes = get_nodes(node_type="drug", node_ids=batch_ids)
        drugs.extend(drug_nodes)

    edges.extend(edges_drugs)

    return genes, disorders, edges, variants, drugs
    

def get_associated_disorders_for_genes(genes, disorders, edges) -> dict[str, any]:
    if genes and disorders and edges:
        gene2disorder = {}
        for edge in edges:
            if edge["type"] == "GeneAssociatedWithDisorder":
                gene_id = edge["sourceDomainId"]
                disorder_id = edge["targetDomainId"]
                gene2disorder.setdefault(gene_id, []).append({"disorder":disorder_id, "dataSources": edge["dataSources"]})
        return gene2disorder
    return None

# note: not used yet
def get_associated_disorders_for_variants(variants, disorders, edges) -> dict[str, any]:
    if variants and disorders and edges:
        variant2disorder = {}
        for edge in edges:
            if edge["type"] == "VariantAssociatedWithDisorder":
                variant_id = edge["sourceDomainId"]
                disorder_id = edge["targetDomainId"]
                variant2disorder.setdefault(variant_id, []).append({"disorder":disorder_id, "dataSources": edge["dataSources"]})
        return variant2disorder
    return None

# note: not used yet
def get_variants_to_affect_gene(variants, genes, edges) -> dict[str, any]:
    if variants and genes and edges:
        gene2variant = {}
        for edge in edges:
            if edge["type"] == "VariantAffectsGene":
                variant_id = edge["sourceDomainId"]
                gene_id = edge["targetDomainId"]
                gene2variant.setdefault(gene_id, []).append({"variant":variant_id, "dataSources": edge["dataSources"]})
        return gene2variant
    return None

def get_drugs_targeting_protein(drugs, edges) -> dict[str, any]:
    if drugs and edges:
        protein2drug = {}
        for edge in edges:
            if edge["type"] == "DrugHasTarget":
                drug_id = edge["sourceDomainId"]
                protein_id = edge["targetDomainId"]
                protein2drug.setdefault(protein_id, []).append({"drug":drug_id, "dataSources": edge["dataSources"]})
        return protein2drug
    return None
    


class BioPAXFactory:
    def __init__(self, input_path: Path, id_space: str = "entrez"):
        self.input_path = input_path
        self.id_space = id_space
        self.g = None
        self.biopaxmodel = None
        self.xRefs: dict[str, biopax.Xref] = {}
        self.entityRefs: dict[str, biopax.EntityReference] = {}
        self.entities: dict[str, biopax.BioPaxObject] = {}
        self.edgeTypes: dict[str, biopax.RelationshipTypeVocabulary] = {
            "gene_associated_with_disorder": biopax.RelationshipTypeVocabulary(term = ["gene associated with disorder"], uid = "gene_associated_with_disorder.vocab"),
            #"variant_associated_with_disorder": biopax.RelationshipTypeVocabulary(term = ["variant associated with disorder"], uid = "variant_associated_with_disorder.XREF"),
            #"variant_affects_gene": biopax.RelationshipTypeVocabulary(term = ["variant affects gene"], uid = "variant_affects_gene.XREF"),
            "drug_has_target": biopax.RelationshipTypeVocabulary(term = ["drug has target"], uid = "drug_has_target.vocab"),
            "gene_product": biopax.RelationshipTypeVocabulary(term = ["gene product"], uid = "gene_product.vocab"),
            #"drug_has_indication": biopax.RelationshipTypeVocabulary(term = ["drug has indication"], uid = "drug_has_indication.XREF")
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
            entrez_ids = [f"entrez.{self.g.vp['name'][v_i]}" for v_i in self.g.get_vertices()]
            gene2prot = get_uniprot_from_entrez(entrez_ids)

            # get all data from nedrex
            genes, disorders, edges, variants, drugs = get_nedrex_data(entrez_ids, gene2prot)

            gene2disorder= get_associated_disorders_for_genes(genes, disorders, edges)
            #variant2disorder = get_associated_disorders_for_variants(variants, disorders, edges)
            #gene2variant = get_variants_to_affect_gene(variants, genes, edges)
            protein2drug = get_drugs_targeting_protein(drugs, edges)
            # TODO: add disorders as soon as Biopax supports it

            for k, v in gene2prot.items():
                [self.add_protein(v_i, k) for v_i in v]
            self.add_drug_info(protein2drug)
            self.add_gene_info(gene2disorder, genes)
        elif self.id_space == "uniprot":
            uniprot_ids = [f"entrez.{self.g.vp['name'][v_i]}" for v_i in self.g.get_vertices()]
            self.add_protein_info()
        self.biopaxmodel = biopax.BioPaxModel(objects=self.get_bioax_objects())

    def add_gene_info(self, gene2disorder, genes):
        for v_i in self.g.get_vertices():
            entrez_id = self.g.vp["name"][v_i]
            if "entrez." + entrez_id in gene2disorder:
                self.add_gene(entrez_id, genes["entrez." + entrez_id] ,gene2disorder["entrez." + entrez_id])
            else:
                self.add_gene(entrez_id, genes["entrez." + entrez_id])
    
    def add_drug_info(self, protein2drug):
        for p, drugs in protein2drug.items():
            for drug in drugs:
                self.add_drug(drug, p)

    def add_drug(self, drug, uniprot_id):
        drug_id = drug["drug"]
        uniXRef = [self.xRefs.setdefault(
            drug_id, biopax.UnificationXref(uid=f"{drug_id}.XREF", db=drug["dataSources"], id=drug_id)
        )]
        if uniprot_id:
            uniXRef.append(self.xRefs.setdefault(
                uniprot_id, biopax.RelationshipXref(uid=f"{uniprot_id}.XREF", db="uniprot", id=uniprot_id, relationship_type=self.edgeTypes["drug_has_target"])
            ))
        entityRef = self.entityRefs.setdefault(
            drug_id, biopax.SmallMoleculeReference(uid=f"{drug_id}.REF", xref=uniXRef)
        )
        self.entities[drug_id] = biopax.SmallMolecule(uid=drug_id, entity_reference=entityRef)

    def add_protein_info(self):
        for v_i in self.g.get_vertices():
            uniprot_id = self.g.vp["name"][v_i]
            self.add_protein(uniprot_id)
        for e_i, e_j in self.g.get_edges():
            uniprot_id1 = self.g.vp["name"][e_i]
            uniprot_id2 = self.g.vp["name"][e_j]
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

    def add_protein(self, uniprot_id, gene_id = None):
        # TODO: add all the other properties

        # uniprot_id = self.g.vp['name'][v]
        uniXRef = [self.xRefs.setdefault(
            uniprot_id, biopax.UnificationXref(uid=f"{uniprot_id}.XREF", db="uniprot", id=uniprot_id)
        )]
        if gene_id:
            uniXRef.append(self.xRefs.setdefault(
                gene_id, biopax.RelationshipXref(uid=f"{gene_id}.XREF", db="NCBI GENE", id=gene_id, relationship_type=self.edgeTypes["gene_product"])
            ))
        entityRef = self.entityRefs.setdefault(
            uniprot_id, biopax.ProteinReference(uid=f"{uniprot_id}.REF", xref=uniXRef)
        )
        self.entities[uniprot_id] = biopax.Protein(uid=uniprot_id, entity_reference=entityRef)

    def get_bioax_objects(self):
        return list(self.xRefs.values()) + list(self.entityRefs.values()) + list(self.entities.values()) + list(self.edgeTypes.values())

    def add_PPI(self, uniprot_id1, uniprot_id2):
        interaction_id = f"{uniprot_id1}_{uniprot_id2}"
        self.entities[interaction_id] = biopax.MolecularInteraction(
            uid=interaction_id, participant=[self.entities[uniprot_id1], self.entities[uniprot_id2]]
        )

    def add_gene(self, entrez_id, gene, associated_disorders = None):
        uniXRef = [self.xRefs.setdefault(
            entrez_id, biopax.UnificationXref(uid=f"{entrez_id}.XREF", db="NCBI GENE", id=entrez_id)
        )]
        if associated_disorders is not None:
            for disorder in associated_disorders:
                id = disorder["disorder"]
                uniXRef.append(self.xRefs.setdefault(
                    id, biopax.RelationshipXref(uid=f"{id}.XREF", db="MONDO", id=id, comment=disorder["dataSources"], relationship_type= self.edgeTypes["gene_associated_with_disorder"])
                ))
        self.entities[entrez_id] = biopax.Gene(uid=entrez_id, xref=uniXRef, organism=None, display_name=[gene["displayName"]])


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
        choices=["entrez", "uniprot"],
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
