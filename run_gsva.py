##############################################################################
#   Given a TSV file of expression data, run the GSVA algorithm
##############################################################################

from optparse import OptionParser
import json
import numpy as np
import pandas as pd
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.conversion import localconverter
import rpy2.robjects as ro
import sys

def run_GSVA(df_counts, gene_set_to_genes):
    gene_set_names = []
    gene_lists = []
    for gene_set, genes in gene_set_to_genes.items():
        gene_set_names.append(gene_set)
        gene_lists.append(genes)
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_cts = ro.conversion.py2rpy(df_counts)
        r_genes = ro.conversion.py2rpy(list(df_counts.index))
        r_gene_set_names = ro.conversion.py2rpy(gene_set_names)
        r_gene_lists = ro.conversion.py2rpy(gene_lists)
        r_gs = ro.vectors.ListVector(gene_set_to_genes)
    rstring="""
        function(counts, genes, gs) {
            library(GSVA)
            gs <- lapply(gs, as.character)
            cts <- as.matrix(counts)
            rownames(cts) <- unlist(genes)
            colnames(cts) <- colnames(counts)
            res <- gsva(cts, gs, kcdf="Poisson")
            df <- data.frame(res)
            df
        }
    """
    r_func = ro.r(rstring)
    r_res = r_func(r_cts, r_genes, r_gs)
    with localconverter(ro.default_converter + pandas2ri.converter):
        df_res = ro.conversion.rpy2py(r_res)
    print(df_res)
    return df_res


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
    usage = "" # TODO
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--transpose", action="store_true", help="Take transpose of input")
    parser.add_option("-o", "--out_file", help="Output file")
    (options, args) = parser.parse_args()

    data_f = args[0]
    gene_sets_f = args[1]
    out_f = options.out_file

    gene_set_to_genes = _parse_gene_sets(gene_sets_f)

    df = pd.read_csv(data_f, sep='\t', index_col=0)
    if options.transpose:
        df = df.transpose()

    res_df = run_GSVA(df, gene_set_to_genes)
    res_df.to_csv(out_f, sep='\t')

if __name__ == "__main__":
    main()
