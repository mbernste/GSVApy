## GSVApy: A Python wrapper for the GSVA algorithm

This is a Python script for running the [Gene Set Variation Analysis (GSVA) algorithm](https://doi.org/10.1186/1471-2105-14-7). Specifically, it is a wrapper around the [GSVA R package](https://www.bioconductor.org/packages/release/bioc/html/GSVA.html).

Specifically, given a genes-by-samples gene expression matrix, this algorithm outputs a gene\_set-by-sample matrix with an enrichment score for each gene set within each sample. 

### Dependencies

This script requires the following Python packages as described in `requirements.txt`:
* rpy2
* optparse
* pandas

This script requires the following R packages:
* [GSVA](https://www.bioconductor.org/packages/release/bioc/html/GSVA.html)

I recommend installing these packages in a conda virtual environment. I have found the interaction between R and Python through rpy2 to be a bit fragile when using global installations. 

### Running on the command line

The `run_gsva.py` script can be run on the command as follow:

```
Usage: python run_gsva.py <input_expression_data> <input_GMT_gene_set_file>

Options:
  -h, --help            show this help message and exit
  -t, --transpose       Take transpose of input
  -d DISTRIBUTION, --distribution=DISTRIBUTION
                        Distribution to use in GSVA {'Poisson' or 'Gaussian'}
  -o OUT_FILE, --out_file=OUT_FILE
                        Output file
```

We provide an example gene expression dataset in `example_data/example.tsv`. We also provide gene sets corresponding to Gene Ontology (GO) terms in `gene_sets/c5.bp.v7.1.symbols.gmt`. To run GSVA on this example data using these gene sets, one would run the command:

`python run_gsva.py -t -d Poisson -o example_GSVA_output.tsv example_data/example.tsv gene_sets/c5.bp.v7.1.symbols.gmt`

Note, the `-t` parameter is required here because the input expression matrix is a samples-by-genes matrix rather than a genes-by-samples matrix and thus, we must take the transpose of the matrix before passing it to GSVA.
