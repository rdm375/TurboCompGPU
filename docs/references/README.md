# Archived research papers

## Kovář et al. (2021)

The PDF in this directory is the publisher's article:

> Patrik Kovář, Adam Tater, Pavel Mačák, and Tomáš Vampola, “Searching for the Most Suitable Loss Model Set for Subsonic Centrifugal Compressors Using an Improved Method for Off-Design Performance Prediction,” *Energies* 14(24), 8545 (2021). https://doi.org/10.3390/en14248545

Publisher page: <https://www.mdpi.com/1996-1073/14/24/8545>

The publisher identifies the article as open access under the Creative Commons Attribution 4.0 International license (CC BY 4.0): <https://creativecommons.org/licenses/by/4.0/>. The article is included here with attribution under that license.

## Chen et al. (2025)

Hua Chen, Gang Li, Shengping Wang, Ning Wang, Lifeng Zhou, Hao Zhou, Yukang
Sun, and Lijun Liu, “The Influence of Variable Operating Conditions and
Components on the Performance of Centrifugal Compressors in Natural Gas Storage
Reservoirs,” *Energies* 18(15), 3930 (2025).
[DOI: 10.3390/en18153930](https://doi.org/10.3390/en18153930).

Publisher PDF: `chen-etal-2025-natural-gas-storage-compressor.pdf`. The publisher
identifies this article as open access under CC BY 4.0. Its Tables 2-4 are
transcribed, with source labels and units, in `../../data/validation/`.

## Lin et al. (2026)

Yuming Lin, Shuai Wang, Chuanyu Zhang, Yuhui Liu, Yuxuan He, Zhiyi Xiong, Yang
Xi, and Weichao Yu, “Healthy-State Performance Modeling of a Multistage Natural
Gas Centrifugal Compressor Using a CFD-Generated Baseline and Factory-Data
Correction,” *Processes* 14(13), 2154 (2026).
[DOI: 10.3390/pr14132154](https://doi.org/10.3390/pr14132154).

Publisher PDF: `lin-etal-2026-natural-gas-compressor.pdf`. The publisher
identifies this article as open access under CC BY 4.0. Its data availability
statement says the data are not public due to confidentiality restrictions.

Both PDFs are reproduced from the publisher under the Creative Commons
Attribution 4.0 International license. Attribution: the authors named above;
source publisher: MDPI; original article DOI links are given above; license:
<https://creativecommons.org/licenses/by/4.0/>. SHA-256:

- `chen-etal-2025-natural-gas-storage-compressor.pdf`: `31ebc5111ab405dcc0cd690afad3285e899f4218f05640aee014193087a620a2`
- `lin-etal-2026-natural-gas-compressor.pdf`: `68bb52341472770133d2e0d3c5be8424e9986d14e1723fce90d83fe0a2a9cf9e`

## Additional physics-validation references

### Yang, Liu, and Zhao (2023)

Xiuxin Yang, Yan Liu, and Guang Zhao, “Modification and Validation of 1D Loss
Models for the Off-Design Performance Prediction of Centrifugal Compressors
with Splitter Blades,” *Machines* 11(1), 118 (2023).
[DOI: 10.3390/machines11010118](https://doi.org/10.3390/machines11010118).

Publisher PDF: `yang-liu-zhao-2023-splitter-blade-loss-model-validation.pdf`.
The publisher article is open access under CC BY 4.0. Attribution: Yang, Liu,
and Zhao; *Machines* (2023); DOI above; license
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). SHA-256:
`548678222af5e783ab4b5d622843db2e667ca0e42f5a3eb3d2ba62311c863fc8`.

This is a 1D loss-model and splitter-blade comparison, including Eckardt-O and
SRV2-O examples. Its optimized-coefficient scores are calibration results and
must not be reported as independent validation scores for TurboCompGPU.

### Avramenko, Frolov, and Hämäläinen (2015)

Anna Avramenko, Alexey Frolov, and Jari Hämäläinen, “Simulation and Modeling of
Flow in a Gas Compressor,” *Journal of Applied Mathematics* (2015), Article
682469. [DOI: 10.1155/2015/682469](https://doi.org/10.1155/2015/682469).

The publisher makes the open-access paper available at the DOI link and states
that its license is a Creative Commons Attribution License. The article's
eight experimental regimes are transcribed, with units and pressure definition,
in `../../data/validation/avramenko_etal_2015_diffuser.csv`. The published
section-average static pressure ratios must not be treated as total-pressure
ratios. A local PDF copy is not included because the publisher download was
unavailable from the repository environment; the DOI remains the version of
record.

### Srinivasan, Kurz, and David (2024)

Anand Srinivasan, Rainer Kurz, and Timothy David, “A Practical Evaluation of
Equations of State and the New PTC-10 Polytropic Calculations for Centrifugal
Compressors,” Proceedings of ASME Turbo Expo 2024, paper GT2024-125064.
[DOI: 10.1115/GT2024-125064](https://doi.org/10.1115/GT2024-125064).

ASME's [paper record](https://asme-turboexpo.secure-platform.com/a/solicitations/223/sessiongallery/15718/application/125064)
describes its natural-gas EOS and PTC-10 comparisons. The paper is cited here;
its full text is not copied into this repository. Its test points are
multistage/whole-machine thermodynamic checks, not matched single-stage loss
map validation for TurboCompGPU. The paper excludes GERG and AGA from its EOS
comparison, so it is complementary to, not a replacement for, our GERG screen.
