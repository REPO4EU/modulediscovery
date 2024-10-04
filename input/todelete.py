import pandas as pd

def modify_target_domain_id(input_file, output_file, mapping):
    """Modify targetDomainId in the input file using the mapping."""
    # Read the input TSV file
    df = pd.read_csv(input_file, sep='\t')
    # Modify targetDomainId based on mapping
    df['targetDomainId'] = df['targetDomainId'].apply(lambda x: mapping.get(x, x))
    # Write the modified DataFrame to a new TSV file
    df.to_csv(output_file, sep='\t', index=False)


if __name__ == "__main__":
    input_file = 'nedrex_drug_target_genename_20240208_nedrex.tsv'
    output_file = 'nedrex_drug_target_geneid_20240208_nedrex.tsv'
    mapping_file = 'mapping.tsv'
    with open(mapping_file) as mapping:
        mapping = {gene: k for k, gene in (l.strip().split("\t") for l in mapping)}
    modify_target_domain_id(input_file, output_file, mapping)

