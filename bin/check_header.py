#!/usr/bin/env python

import sys
import argparse

def process_file(input_file):

    expected_content = "gene_id"

    with open(input_file, 'r+') as file:
        content = file.read()

        if not content.startswith(expected_content):
           
            file.seek(0, 0)

            file.write(f"{expected_content}\n{content}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Checks if file has correct input for gprofiler.')
    parser.add_argument('--input', required=True, help='Input GT file path') 
    
    args = parser.parse_args()
    input_file = args.input 

    process_file(input_file)