# Huntsman Lab Single Cell RNA Sequencing Pipeline

The Huntsman Lab single cell RNA sequencing pipeline is implemented in R (4.0.4), and uses the `SingleCellExperiment` class to interact with the dataset. The pipeline has been designed with 10x Genomics output in mind. If other scRNA-seq technologies are to be fed into this pipeline, implementation of additional features, and refacting may be necessary.

The pipeline was inspired by [Orchestrating Single-Cell Analysis](https://bioconductor.org/books/release/OSCA/index.html).

## Installation

This packages can be installed directly from GitHub. Do not submit this package to CRAN, or Bioconductor.

```{R}
devtools::install_github('Huntsmanlab/blackbox)
```

## Background + Expected practice

The pipeline has been designed as a R package, and as such any changes to the pipeline should comply with R package development conventions. For information on R package development you can read the [R packages](https://r-pkgs.org/) book by Hadley Wickham.

Below is the current directory struture for our the `blackbox` package. All scripts related to processing and analyzing scRNA-seq data is saved in the `R/` directory. All functions are to be documented using [`roxygen2`](https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html).

If any new libraries are to be added to the pipeline, be sure to include the in the `DESCRIPTION` file, along with the version used when implemented in the pipeline.

Do **not** push any changes to the `master` branch without having another review your code, and being sure it does not break any currently implemented protocols. Make your own branch by cloning this repo, and executing the following command in your terminal: `$ git branch <BRANCH-NAME>`. 

After any changes to the pipeline, always run `devtools::document(roclets = c('rd', 'collate', 'namespace'))`, or `CMD + SHIFT + D` to update the `NAMESPACE` of the pipeline. All new protocols are to be added to the dispatcher functions. Protocols should not be public. 

```
.
├── DESCRIPTION
├── LICENSE
├── LICENSE.md
├── NAMESPACE
├── R/
│   ├── clustering.R
│   ├── countQC.R
│   ├── diffExp.R
│   ├── dimReduced.R
│   ├── featureSelect.R
│   ├── normCounts.R
│   ├── prep.R
│   └── utils.R
├── README.md
├── blackbox.Rproj
├── man/
```

### Conventions

Lets be sure to keep true to these conventions. If these conventions are not maintained with an evergrowing codebase... entropy will increase making it really challenging for debugging and increase the amount of training for future colleagues.

1. All function which are not dispatchers should be prefixed by a `.` to indicate that they are private.

2. Dispatcher functions are to be named in camelCase.

3. Helper, or additional functions, are to be named in lower case and separated by `.`.

4. Constants must be in uppercase.

5. Functions and variable must have clear and relevant names. Autocomplete exists.

6. Avoid using `for` and `while` loops for large data structures. 

7. Adding `TODO` to uncompleted fuctions.

## Project structure (TODO: implement command using python and click)

To create a new scRNA-seq analysis, run the following command.

```{zsh}
$ blackbox create DH
```

The above command will generate the directory tree shown below.

```
.
├── README.md
├── data
│   ├── diagnostic
│   ├── processed
│   ├── raw
│   │   └── DH
│   │       ├── barcodes.tsv
│   │       ├── genes.tsv
│   │       └── matrix.mtx
│   └── temp
├── dh.Rproj
├── output
├── run_analysis.R
```

## Workflow

The main function for each stage of the workflow is has the same name as the file it is sourced from, except for `source.R`. 

All data is saved in the `processed` directory as a `RDS` after the main function has completed. RDS files are named after the source script that generated them. The main functions will also return a `SingleCellExperiment` object.

1. `source.R`: Prepares environment of scRNA-seq analysis by importing various packages from CRAN and Bioconductor. `prep.data` function is **required** to run on the 10x Genomics dataset, prior to starting analysis.

2. `countQC`: Will perform quality control counts on `SingleCellExperiment` object using fixed or adaptive thresholds. The fixed threshold can be altered by the user. The adaptive thresholds are calculated by quartiles ranges.

3. `normCounts.R`: Normalized quality controlled counts using library size factors (`vanilla`) or deconvolution method. 

4. `featureSelect.R`: Selected highly variable genes using either `modelGeneVar` or `modelGeneByPoisson`, as defined by the user.

5. `dimReduced.R`: Performs dimensionality reduction through PCA, and PCA-initialized TSNE and UMAP.

6. `clustering.R`

## TODOs

*Prefix with `TODO` for these TODOs to appear in the TODO file.*

To generate a TODO file, run the following command in your terminal.

```{zsh}
$ rg --stats --regexp "TODO" --glob "\!man"
```

TODO: Find a better method for passing along `IDENTIFIER` rather than making it into a temporary global variable. `rm` will not delete the global variable made within the same scope.
TODO: Everything upto `batchCorrect.R` is operating on a single `SingleCellExperiment` object. Should I change everything to operator on a list of `SingleCellExperiment`?
TODO: Subset context params within a scope fo a function rather than as a paramter?
