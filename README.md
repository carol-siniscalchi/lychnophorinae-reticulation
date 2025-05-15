# lychnophorinae-reticulation

Data, scripts and figures from "Phylogenomics reveal widespread ancient hybridization in the evolution of Lychnophorinae (Vernonieae, Asteraceae)", under revision in TAXON.

## Repository structure:

- **alignments:** contains all the aligments used for analyses, separated into phyluce and HybPiper assemblies.
- **assembly_comparisons:** this folder contains the R script used to plot figures in the Appendix, with the necessary data files to run it. PDF of the figures are provided.
- **PhyloNetworks:** this folder has all the inputs and outputs of all the seven PhyloNetworks analyses. Files are named similarly in each folder, except for the input trees. The scripts folder has a general script that can be used to replicate analysis. The loglik_plot folder has the input and script to produce the plot in Fig 4I. TICR_results has the p-values of all seven tests.
- **QuartetSampling:** contains the species tree and concatenated alignment needed to run the QS analysis, plus all the output files. The figtree file was used as basis for Fig 3.
- **tree_comparisons:** contains script and data files used to generate the overlayed trees from Fig 2, plus some of the Supplemental Figures. PDFs are provided.
- **trees:** contains all trees generated, divided by assembly method.

A static copy of this repository is provided at Zenodo (10.5281/zenodo.15426471). 
