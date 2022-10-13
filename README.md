![github_front_page_figure](https://user-images.githubusercontent.com/64554648/182730384-e9e60ead-8651-4c23-a484-b69bdd6ffe88.png)

---

This repository contains analysis code to accompany the publication ['Distinct spatiotemporal movement properties reveal sub-modalities in crawling cell types'](https://research.arcadiascience.com/pub/result-comparative-crawling).

Notebooks containing expanded methods and the code + analyses for generating all figures in the publication can be run in Binder:
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Arcadia-Science/comparative-amoeboid-crawling/HEAD)

**Note for Binder users:** For repositories that rely on a lot of packages, the launch process may take some time and you may see the error message `Failed to connect to event stream`. This means that the connection between the logs service and the UI is broken. But the build is still happening in the background and the UI should fix itself upon refreshing the page. You can learn more about this issue [here](https://discourse.jupyter.org/t/failed-to-connect-to-event-stream/527).

## Directory structure

`00_data/` Cell shape data matrices for all cell types/trials examined in the pub.<br>
`01_utils/` .R code necessary to run the analyses presented in the pub and accompanying Jupyter notebooks.<br>
`02_analysis_files/` Cached intermediate files generated from the analyses contained in the repo, including the movement space embeddings used in the pub to facilitate replication of the published figures.<br>
`03_notebooks/` Jupyter notebooks outlining the full suite of analyses presented in the pub.<br>

## R packages and versions used in this repo:

`umap v0.2.7.0`<br>
`scales v1.1.1`<br>
`MASS v7.3-54`<br>
`RColorBrewer v1.1-2`<br>
`colorRamps v2.3`<br>
`vegan v2.5-7`<br>
`igraph v1.2.11`<br>
`entropy v1.3.1`<br>
`jpeg v0.1-9`<br>
`repr v1.1.4`<br>
`swaRm v0.5.0`<br>
`dunn.test v1.3.5`<br>
`vioplot v0.3.7`<br>
`lsa v0.73.2`<br>
`shape v1.4.6`<br>
`gtools v3.9.2`<br>
`pracma v2.3.6`<br>
`fractaldim v0.8-5`<br>
`here v1.0.0`<br>
`gplots v3.1.3`<br>
