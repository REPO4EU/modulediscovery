# Including a new active module detection tool

1. Create a new branch for your tool.
2. Add a function to the bin/graph_tool_parser.py script for preparing the tool-specific network input format. The script is built around the [graph-tool](https://graph-tool.skewed.de/) Python package. An example is the safe_diamond() function, which saves a simple edge list in CSV format. Add the function as an option in the save() function and a command line option for `--format` in parse_args(). The output file name has to include the option specified with `--format` since nextflow uses this pattern to check whether the output file was successfully generated. The script expects a .gt file as input. Run the pipeline with the "test" profile to generate a .gt example file in `<OUTDIR>/graphtoolparser`, which you can use to test the parsing function by executing the parsing script directly via the command line.
3. Create a module for the tool. (Example with comments: `modules/local/diamond/main.nf` and `modules/local/domino/`)
4. Create a subworkflow wrapping the tool together with the input parser. (Example with comments: `subworkflows/local/gt_diamond/main.nf` and `subworkflows/local/gt_domino/main.nf`)
5. Include the subworkflow in the `workflows/modulediscovery.nf` file. Again, DIAMOnD and DOMINO are included as examples.
6. Test checks locally:
   1. Run tests via, e.g., `nextflow run main.nf -profile singularity,test --outdir results`.
   2. Run `nf-core lint`.
   3. Check your code style. This will automatically happen before you commit, if you use pre-commit, which can be set up with: `pre-commit install`. After each commit, it will automatically check your code style and fix it where possible. If changes were made, you have to commit again.
7. Create a pull request against the dev branch.

Further information

- [FAQ sheet](https://docs.google.com/document/d/1WgBIFrrcxFKN0I-zJbuS7PUCmyCLPTWx6xAHg1zi4FA/edit?usp=sharing)
- [Workflow schema](https://docs.google.com/drawings/d/1X7U79dAZaeRdGdIsXoEKw74MNqjxCHq3RuNASBYCiB4/edit?usp=sharing)
