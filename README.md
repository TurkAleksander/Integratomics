This repo contains all the necessary materials to run your own positional integration analysis, with methods based on Breitling et al. (2004) and Knijnenburg et al. (2009).
The analysis can be run either directly from R or from within a Docker for stability.
If you experience problems with the Docker, you may want to check the versions of the dependencies that are being installed for Alpine Linux.

File descriptions:
- **R_integration_code.r**: Code for running the analysis with placeholder directories. It is NOT recommended to run on a local PC, as the permutation tests are highly RAM-intensive. Both this code and the Docker version are quite comment-rich for clarity.
- **R_integration_code_dockerVersion.r**: Dockerized version of the integration code. Draws directories from arguments passed when running the Docker. Will not run without an input directory. If you do not specify an output, it will output results to your input directory.
- **Dockerfile_integration**: Dockerfile with instructions on making the image. Downloads Alpine Linux, installs dependencies, and necessary packages and clones this repository anew with every build. Uses specific versions of packages for stability, this can be changed here. If you do so, run the build command with --no-cache to reinstall everything.
- **Docker_build_template.sh**: How to build the Docker image - this is where you can also specify the Docker image's name.
- **Docker_run_template.sh**: How to run the Docker with placeholder directories - here you specify where the input files are on your system, as well as where you want the output to be.
- **hg38_UCSC_chrom_lengths.txt**: Chromosome lengths from the UCSC browser - needed in case a genomic backbone needs to be constructed from scratch (not recommended - it takes a while as it's not optimized at all).
- **locationBackbone.txt**: A table of 10,000 base pair intervals along each chromosome. Each interval has a 5,000bp overlap with its neighbouring intervals.
- **locationBackbone_R-version.txt**: Same as locationBackbone.txt, but constructed using the in-built code ("from scratch"). The only difference is that it uses a Unix line ending ("\n") instead of ("\r\n").
- **mart_export.txt**: Holds the location (chr, start, end) of protein-coding genes with a UCSC stable-ID on chromosomes 1-22, X, Y and MT (build GRCh38.p14) - 20,037 genes. This file is needed to both estimate gene density of intervals and to see which genes map to them, which is handy for the results.

For _details_ on the integration procedure itself see the **Supplementary Methods** section of the article.

Helpful Docker guides:
- https://docs.docker.com/get-started/docker-overview/
- https://docker-curriculum.com/
