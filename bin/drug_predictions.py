#!/usr/bin/env python
import argparse
from pathlib import Path
import logging
import sys
import pandas as pd
import drugstone
from drugstone import new_task

drugstone.print_license()
drugstone.accept_license()


logger = logging.getLogger()


class DrugPredictions:
    def __init__(self, input_path: Path, id_space: str = "entrez"):
        self.input_path = input_path
        self.id_space = id_space
        self.df = None
        self.nodes = None

    def load_file(self, input_path):
        logger.debug("Loading file from %s", input_path)
        self.df = pd.read_csv(input_path, sep="\t")
        self.nodes = set(self.df["name"])

    # TODO: adjust output path
    def get_drug_csv_path(self) -> Path:
        return str(self.input_path).replace(".tsv", "") + "_trustrank_drugs"

    def get_drug_set_trustrank(self, id_set, identifier, filename):
        parameters = {
            "identifier": identifier,
            "algorithm": "trustrank",
            "target": "drug",
            "includeIndirectDrugs": False,
            "includeNonApprovedDrugs": False,
            "result_size": 50,
        }
        task = new_task(id_set, parameters)
        r = task.get_result()
        detailsForNodes = r.get_raw_result()["nodeAttributes"]["details"]
        r.download_drugs_csv(name=filename)
        drugs = r.get_drugs()
        return [
            {
                "id": entry["drugId"],
                "score": entry["score"],
                "proteins_targeted": entry["hasEdgesTo"],
                "status": entry["status"],
                "label": entry["label"],
            }
            for entry in drugs.values()
        ], detailsForNodes

    def parse_drug_predictions(self, drugs, nodeDetails):
        expanded_df_rows = []

        for index, row in self.df.iterrows():
            node_id = row["name"]
            details = nodeDetails.get(node_id, None)
            if details:
                symbol = details.get("symbol", None)
                proteins_symbols_str = " & ".join(symbol)
                row["symbol"] = proteins_symbols_str

            drug_rows = []

            for drug in drugs:
                if node_id in drug["proteins_targeted"]:
                    new_row = row.copy()

                    new_row["drug_id"] = drug["id"]
                    new_row["score"] = drug["score"]
                    new_row["drug_name"] = drug["label"]
                    new_row["status"] = drug["status"]

                    drug_rows.append(new_row)

            if drug_rows:
                expanded_df_rows.extend(drug_rows)
            else:
                expanded_df_rows.append(row)

        expanded_df = pd.DataFrame(expanded_df_rows)

        df_sorted_within_groups = (
            expanded_df.groupby("name", group_keys=False)
            .apply(lambda x: x.sort_values(by=["score"], ascending=[False]))
            .reset_index(drop=True)
        )

        self.df = df_sorted_within_groups

    def create_drug_predictions(self):
        if self.df is None:
            self.load_file(self.input_path)
        logger.debug("Creating drug predictions")
        drugs, nodeDetails = self.get_drug_set_trustrank(
            self.nodes, self.id_space, self.get_drug_csv_path()
        )
        self.parse_drug_predictions(drugs, nodeDetails)
        self.df.to_csv(self.get_drug_csv_path() + ".tsv", sep="\t", index=False)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse tsv file and create drug predictions.",
        epilog="Example: python3 drug_predictions.py --namespace entrez network.nodes.tsv",
    )

    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="TSV file of the nodes of the network.",
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
    predictor = DrugPredictions(args.file_in, args.idspace)
    predictor.create_drug_predictions()


if __name__ == "__main__":
    sys.exit(main())
