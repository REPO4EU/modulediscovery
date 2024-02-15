#!/usr/bin/env python

import sys
import argparse


def add_header(input_file, header, output_file):

    with open(input_file, "r") as file:
        content = file.read()

    if not content.startswith(header):
        with open(output_file, "w") as output:
            output.write(f"{header}\n{content}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Checks if file has correct input for gprofiler."
    )
    parser.add_argument("--input_file", required=True, help="Input GT file path")
    parser.add_argument("--header", required=True, help="for examplea 'gene id'")
    parser.add_argument("--output_file", required=True, help="Output file path")

    args = parser.parse_args()
    input_file = args.input_file
    header = args.header
    output_file = args.output_file

    add_header(input_file, header, output_file)
