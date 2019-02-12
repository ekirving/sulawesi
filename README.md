# Sulawesi
This repository contains R code for plotting Figure 4 from 
the paper [Synchronous diversification of Sulawesi's iconic artiodactyls driven by recent geological events](https://doi.org/10.1098/rspb.2017.2566)

![Figure 4](./pdf/Figure_4.png?raw=true)


If you reuse any of this code then please cite the paper:
> Frantz, L.A.F., Rudzinski, A., Nugraha, A.M.S., Evin, A., Burton, J., Hulme-Beaman, A., Linderholm, A., Barnett, R., 
> Vega, R., Irving-Pease, E.K., Haile, J., Allen, R., Leus, K., Shephard, J., Hillyer, M., Gillemot, S., Hurk, J. van 
> den, Ogle, S., Atofanei, C., Thomas, M.G., Johansson, F., Mustari, A.H., Williams, J., Mohamad, K., Damayanti, C.S., 
> Wiryadi, I.D., Obbles, D., Mona, S., Day, H., Yasin, M., Meker, S., McGuire, J.A., Evans, B.J., Rintelen, T. von, Ho, 
> S.Y.W., Searle, J.B., Kitchener, A.C., Macdonald, A.A., Shaw, D.J., Hall, R., Galbusera, P., Larson, G., 2018. 
> Synchronous diversification of Sulawesi’s iconic artiodactyls driven by recent geological events. *Proc. R. Soc. B* 285, 
> 20172566 https://doi.org/10.1098/rspb.2017.2566

## Installation

To reproduce the figure from the paper you will need to install the following dependencies.


### R

R ≥ 3.4 with the following modules:

* [devtools](https://cran.r-project.org/web/packages/devtools/)
* [mapdata](https://cran.r-project.org/web/packages/mapdata/)
* [maptools](https://cran.r-project.org/web/packages/maptools/)
* [maps](https://cran.r-project.org/web/packages/maps/)
* [tess3r](https://cran.r-project.org/web/packages/tess3r/)

```R
install.packages(c("devtools", "mapdata", "maptools", "maps"))
```

```R
devtools::install_github("bcm-uga/TESS3_encho_sen")
```

## Running the code

To reproduce the figure, run:

```bash
Rscript Sula_Tess3.R 
```

## Author

Evan K. Irving-Pease, [PalaeoBARN](https://www.palaeobarn.com/), University of Oxford 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
