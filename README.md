# CatIB-TS

Repository for parsing, plotting and analysis of an _Escherichia coli_ strain library producing glucose dehydrogenase from _Bacillus subtilis_ (_Bs_GDH) as Catalytically active Inclusion Bodies (CatIBs). The overall goal is to perform consecutive screening experiments using a Bayesian process model and Thompson Sampling as a policy to select candidates.

This projects provides the raw data and data analysis notebooks for the manuscript __"High-Throughput Screening of Catalytically Active Inclusion Bodies Using Laboratory Automation and Bayesian Optimization"__ (2023) by Laura M. Helleckes*, Kira Küsters*, Christian Wagner, Rebecca Hamel, Ronja Saborowski, Wolfgang Wiechert, Marco Oldiges.

*These authors contributed equally.

## Structure
Raw data can be found in the `data` folder. Data analysis is conducted in `notebooks`, where plots for the accompanying paper can be found as well.
The results for individual experiments can be found by a unique idetentifier, their so-called Run ID.

The following runs were conducted:

| Preculture ID  | Main Culture ID | Assay ID | Description    |
| -------------- | --------------- | -------- | -------------- |
| D15DTK         | D19YYZ          | D1X8DM   | TS Round 1     |
| D95YC6         | D9A1YM          | D9XCLX   | TS Round 2     |
| DAFT9R         | DAMZ19          | DB984D   | TS Round 3     |
| DDA79K         | DDEBDB          | DE3MB4   | Manual Round 4 |

Thompson Sampling on the reaction rates provided by the process model was performed in three rounds, leading to >95 % certainty in identifying the best performer from the library. The fourth screening round was designed manually to screen the top-10 CatIB variants as well as two previously unseen variants.

## Citation of code
This repository and the corresponding Python package for data analysis (`catibts`) is licensed under the [GNU Affero General Public License v3.0](https://github.com/JuBiotech/Supplement-to-Helleckes-Kuesters-et-al.-2023/blob/main/LICENSE.md).
Head over to Zenodo to generate a BibTeX citation for the latest release.
