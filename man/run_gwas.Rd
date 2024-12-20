% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gwas.R
\name{run_gwas}
\alias{run_gwas}
\alias{run_gwas.bed}
\alias{run_gwas.gds}
\title{Run Genome-Wide Association Study (GWAS)}
\usage{
run_gwas(file, model_params, output_file, chunk_size = 10000, min_maf = NULL)

\method{run_gwas}{bed}(file, model_params, output_file, chunk_size = 10000, min_maf = NULL)

\method{run_gwas}{gds}(file, model_params, output_file, chunk_size = 10000, min_maf = NULL)
}
\arguments{
\item{file}{A character string specifying the path to the genotype file.}

\item{model_params}{A \code{model_params} object.}

\item{output_file}{A character string specifying the base name for the output
files.}

\item{chunk_size}{An integer specifying the number of SNPs to process in each
chunk. Default is 10000.}

\item{min_maf}{A numeric value specifying the minimum minor allele frequency
(MAF) threshold. Default is NULL.}
}
\value{
Invisibly returns \code{TRUE} if the analysis completes successfully.
}
\description{
This function performs a GWAS analysis on genotype data stored in various
file formats.
}
\details{
Now supports the following file formats:
\itemize{
\item plink \code{BED} file
\item \code{GDS} file
}

The results are written to two CSV files: one for valid results and one for
invalid results (e.g., SNPs with MAF below the threshold).
}
\examples{
\dontrun{
todo!
}

}
