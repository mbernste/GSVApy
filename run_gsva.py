##############################################################################
#   Given a TSV file of expression data, run the GSVA algorithm
##############################################################################

from optparse import OptionParser
import json
import pandas as pd
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.conversion import localconverter
import rpy2.robjects as ro
import sys

def run_GSVA(df_expr, gene_set_to_genes, distr='Gaussian'):
    """
    Parameters
    ----------
    df_expr : DataFrame
        A gene-by-sample Pandas DataFrame storing the gene expression 
        matrix.
    gene_set_to_genes: dictionary
        A dictionary mapping each gene set name (string) to list of gene
        names (list of strings).
    distr: string
        The distribution to use in the GSVA algorithm. Must be either 
        'Poisson' or 'Gaussian'.

    Returns
    -------
    df_gsva: DataFrame
        A gene set-by-sample dataframe storing the GSVA scores for each
        gene set within each sample.
    """
    # Re-format gene-set information to make compatible with GSVA
    gene_set_names = []
    gene_lists = []
    for gene_set, genes in gene_set_to_genes.items():
        gene_set_names.append(gene_set)
        gene_lists.append(genes)

    # Convert Python objects to rpy2 objects to pass to R
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_expr = ro.conversion.py2rpy(df_expr)
        r_genes = ro.conversion.py2rpy(list(df_expr.index))
        r_gene_set_names = ro.conversion.py2rpy(gene_set_names)
        r_gene_lists = ro.conversion.py2rpy(gene_lists)
        r_gs = ro.vectors.ListVector(gene_set_to_genes)

    # The R function for running GSVA
    rstring="""
        function(expr, genes, gs) {{
            library(GSVA)
            gs <- lapply(gs, as.character)
            expr <- as.matrix(expr)
            rownames(expr) <- unlist(genes)
            colnames(expr) <- colnames(expr)
            res <- gsva(expr, gs, kcdf="{}", mx.diff=TRUE)
            df <- data.frame(res)
            df
        }}
    """.format(distr)

    # Run the R code
    r_func = ro.r(rstring)
    r_res = r_func(r_expr, r_genes, r_gs)

    # Convert the R output  to Python objects
    with localconverter(ro.default_converter + pandas2ri.converter):
        df_gsva = ro.conversion.rpy2py(r_res)

    return df_gsva


def _parse_gene_sets(gene_sets_f):
    gene_set_to_genes = {}
    with open(gene_sets_f, 'r') as f:
        for l in f:
            toks = l.split('\t')
            gene_set = toks[0]
            genes = toks[2:]
            gene_set_to_genes[gene_set] = genes
    return gene_set_to_genes


def main():
    usage = "python run_gsva.py <input_expression_data> <input_GMT_gene_set_file>"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--transpose", action="store_true", help="Take transpose of input")
    parser.add_option("-d", "--distribution", help="Distribution to use in GSVA {'Poisson' or 'Gaussian'}")
    parser.add_option("-o", "--out_file", help="Output file")
    (options, args) = parser.parse_args()

    data_f = args[0]
    gene_sets_f = args[1]

    if options.distribution is None:
        print("Warning! No distribution was specified (see the '--distribution' flag). Using 'Guassian' by default.")
        distr = 'Gaussian'
    else:
        distr = options.distribution
        assert distr in {'Poisson', 'Gaussian'}, f"The `--distribution` argument must be either 'Poisson' or 'Gaussian'. '{distr}' is invalid."
    out_f = options.out_file

    gene_set_to_genes = _parse_gene_sets(gene_sets_f)

    df = pd.read_csv(data_f, sep='\t', index_col=0)
    if options.transpose:
        df = df.transpose()

    res_df = run_GSVA(df, gene_set_to_genes, distr=distr)
    res_df.to_csv(out_f, sep='\t')

if __name__ == "__main__":
    main()
