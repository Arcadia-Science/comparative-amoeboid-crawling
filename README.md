# Comparative amoeboid crawling

This repository contains analysis code to accompany the publication 'Comparative analyses reveal substantial variation among unicellular organisms' ([link](https://research.arcadiascience.com/pub/result-comparative-crawling/draft)).

Notebooks containing expanded methods and the code + analyses for generating all figures in the publication can be run in Binder:


[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ryanayork/comparative-amoeboid-crawling/main)

The repo is structured as follows:

'00_data': Contains cell shape matrices for all cell types/trials examined in the pub.
'01_utils': Contains .R code necessary to run the analyses presented in the pub and accompanying Jupyter notebooks.
'02_analysis_files': Contains cached intermediate files generated from the analyses contained in the repo, including the movement space embeddings used in the pub to facilitate replication of the published figures.
'03_notebooks': Contains three Jupyter notebooks outlining the full suite of analyses presented in the pub.

R packages and versions used in this repo:
umap v0.2.7.0
scales v1.1.1
MASS v7.3.55
RColorBrewer v1.1.2
colorRamps v2.3
vegan v2.5.7
igraph v1.2.11
entropy v1.3.1
jpeg v0.1.9
repr v1.1.4
swaRm v0.5.0
dunn.test v1.3.5
vioplot v0.3.7
lsa v0.73.3

---

![github_front_page_figure](https://user-images.githubusercontent.com/64554648/182730384-e9e60ead-8651-4c23-a484-b69bdd6ffe88.png)


