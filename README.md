# MS-CleanR: A package for cleaning and annotating LC-MS data

The MS-CleanR package provides functions for feature filtering and annotation of LC-MS data.

See the publication and tutorials (pdf files included in the master branch) for more information.

Needs MS-DIAL (v4.00 or higher) and MS-FINDER (3.30 or higher): http://prime.psc.riken.jp/compms/index.html

### Short description:
![Image of MSC-scheme](https://github.com/eMetaboHUB/MS-CleanR/blob/master/MSC-scheme.png)
MS-CleanR use as input MS-DIAL peak list processed in data dependent analysis (DDA) or data independent analysis (DIA) using either positive ionization mode (PI) or negative ionization mode (NI) or both. First, MS-CleanR apply generic filters encompassing blank injection signal subtraction, background ions drift removal, unusual mass defect filtering, relative standard deviation threshold (RSD) based on sample class and relative mass defect (RMD) window filtering. All these options are tunable by the user. The second step involves a feature clustering method based on MS-DIAL peak character estimation algorithm followed by parental signal extraction using multi-level optimization of modularity algorithm. Optionally, MS-CleanR can merge PI and NI mode during this step. Then, all selected features are exported to MS-FINDER program for in silico-based annotation using hydrogen rearrangement rules (HRR) scoring system. At this step, multiple databases can be queried and each annotation results will be handled by MS-CleanR. The final step will merge annotation results to the filtered peak list by prioritizing database annotation depending on user choice. Optionally, all results can be exported as .msp file for mass spectral similarity networking purpose.   

### Installation
```
devtools::install_github("eMetaboHUB/MS-CleanR")
library(mscleanr)
runGUI() 
```
### Citation
Publication link: https://www.biorxiv.org/content/10.1101/2020.04.09.033308v2

### Credits
- [Université Toulouse III -- Paul Sabatier](https://www.univ-tlse3.fr)
- [MetaToul-AgromiX Platform](https://www.lrsv.ups-tlse.fr/metatoul-en/)

### Licence
GPL-3
