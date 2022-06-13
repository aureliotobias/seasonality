## Seasonality
**R** script files to assess seasonality and the role of its potential drivers in environmental epidemiology studies.

* **main_analysis.R** reproduces the examples from the published manuscript. 
* It is necessary to upload the following **ancillary functions**: 
    * **cyclic.R** - function to generate the clyclic spline (adapted from the R packge DLNM by <a href="https://pubmed.ncbi.nlm.nih.gov/22003319/" target="_blank">Gasparrini 2017</a>).  
    * **findmax.R** - function to estimate peak of seasonality (adapted from the R function to estimate the minimum the minimum mortality temperature by <a href="https://pubmed.ncbi.nlm.nih.gov/27748681/" target="_blank">Tobías et al. 2017</a>). 
    * **findmin.R** - function to estimate trhough of seasonality (adapted from the R function to estimate the minimum the minimum mortality temperatureby <a href="https://pubmed.ncbi.nlm.nih.gov/27748681/" target="_blank">Tobías et al. 2017</a>). 
    * **attrs.R** - function to estimate attributable fraction (adapted from the R function to estimate attributable measures from DLNM by <a href="https://pubmed.ncbi.nlm.nih.gov/24758509/" target="_blank">Gasparrini and Leone 2014</a>). 
<br>
<b>Reference:</b> Madaniyazi L, Tobías A, Kim Y, Chung Y, Armstrong B, Hashizume M. <a href="https://academic.oup.com/ije/advance-article/doi/10.1093/ije/dyac115/6593248" target="_blank">Assessing seasonality and the role of its potential drivers in environmental epidemiology: a tutorial</a>. <b>International Journal of Epidemiology</b> 2022, doi.org/10.1093/ije/dyac115.
