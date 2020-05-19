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
            #print(gs)
            #print(sessionInfo())
            res <- gsva(cts, gs, kcdf="Poisson")
            print(res)
            print(rownames(res))
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


def main():

    #cells = ['PJ016_4', 'PJ016_5', 'PJ016_6']
    #counts = load_GSE103224.counts_matrix_for_cells(cells)

    data_f = sys.argv[1]
    df = pd.read_csv(data_f, sep='\t', index_col=0)
    df = df.transpose()
    with open('gene_set_to_genes.json', 'r') as f:
        gene_set_to_genes = json.load(f)
    run_GSVA(df, gene_set_to_genes)


if __name__ == "__main__":
    main()
