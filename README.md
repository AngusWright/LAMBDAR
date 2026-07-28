# LAMBDAR

LAMBDAR (Lambda Adaptive Multi-Band Deblending Algorithm in R) measures
matched-aperture galaxy fluxes from FITS images. Given an aperture catalogue
and image data, it creates apertures, optionally convolves and deblends them,
and writes photometric fluxes and uncertainties.

The main entry point is `measure.fluxes()`. The package also provides helpers
for aperture construction, point-spread-function handling, astrometry, FITS
I/O, and creating simulated images.

## Installation

Install from GitHub with [remotes](https://remotes.r-lib.org/):

```r
install.packages(c("remotes", "FITSio"))
remotes::install_github("AngusWright/LAMBDAR", dependencies = TRUE)
library(LAMBDAR)
```

LAMBDAR depends on `doParallel`, `foreach`, `magicaxis`, `nortest`, `seqinr`,
`iterators`, `plotrix`, `RANN`, `matrixStats`, `data.table`, `FITSio`, and
[`astro`](https://github.com/AngusWright/astro). `astro` is declared as a
GitHub remote and is installed by `install_github()`. Optional functionality
uses `EBImage` and `reticulate`.

To install a local clone instead:

```sh
R CMD INSTALL LAMBDAR
```

Run this command from the directory that contains the clone.

## Quick start

The package includes an SDSS r-band FITS image and a matching aperture
catalogue (`SDSS.sample` and `ApCat.sample`). The catalogue has 274 rows and
the columns `CATAID`, `RAdeg`, `DECdeg`, `rotW2N`, `radminasec`,
`radmajasec`, and `rotN2E`.

The following reproduces the package's basic measurement workflow:

```r
library(LAMBDAR)

data("SDSS.sample", envir = environment())
data("ApCat.sample", envir = environment())

write.fits(file = "SampleImage.fits", SDSS.sample)
write.csv(
  ApCat.sample, "SampleCat.csv",
  row.names = FALSE, quote = FALSE
)

create.par.file(
  filename = "Lambdar_Sample.par",
  DataMap = "SampleImage.fits",
  Catalogue = "SampleCat.csv"
)
measure.fluxes("Lambdar_Sample.par")

results <- read.csv("LAMBDAR_Results.csv", stringsAsFactors = FALSE)
summary(results$DFAFlux_Jy)
```

`Lambdar_Sample.par` contains all run settings. Generate a fully annotated
default parameter file with:

```r
measure.fluxes("--makepar")
```

The principal result is `LAMBDAR_Results.csv`; with `do.return = TRUE`,
`measure.fluxes()` also returns a list including `SFAflux`, `SFAerror`,
`DFAflux`, and `DFAerror`.

## Function examples

The outputs below were produced by running the shown calls against this
checkout.

### Spherical coordinates to Cartesian coordinates

```r
sph.to.car(c(0, 90), c(0, 0))
```

```text
                x y z
[1,] 1.000000e+00 0 0
[2,] 6.123234e-17 1 0
```

### Pixel coordinates to right ascension and declination

```r
xy.to.radec(
  c(100, 110), c(200, 205),
  ra0 = 150, dec0 = 2, x0 = 100, y0 = 200,
  xscale = 0.001, yscale = 0.001
)
```

```text
         RA   DEC
[1,] 150.00 2.000
[2,] 150.01 2.005
```

### Gaussian FWHM estimate

```r
sd <- 5
x <- matrix(-10:10, 21, 21, byrow = TRUE)
zdist <- matrix(
  exp(-1 * (t(x)^2 / (2 * sd^2) + x^2 / (2 * sd^2))),
  21, 21
)
paste(2 * sqrt(2 * log(2)) * sd, "~", get.fwhm(zdist))
```

```text
[1] "11.7741002251547 ~ 12"
```

The analytic FWHM is approximately 11.77 pixels; `get.fwhm()` returns the
integer pixel estimate, 12.

### Sub-pixel elliptical aperture

```r
grid <- expand.grid(-2:2, -2:2)
aperture <- generate.aperture(
  grid[, 1], grid[, 2],
  xstep = 1, ystep = 1,
  axrat = 0.5, axang = 0, majax = 2,
  resample.iterations = 0, pixscale = TRUE
)

matrix(aperture[, 3], nrow = 5)
sum(aperture[, 3])
```

```text
     [,1] [,2] [,3] [,4] [,5]
[1,] 0.00  0.0  0.0  0.0 0.00
[2,] 0.25  0.5  0.5  0.5 0.25
[3,] 0.50  1.0  1.0  1.0 0.50
[4,] 0.25  0.5  0.5  0.5 0.25
[5,] 0.00  0.0  0.0  0.0 0.00
[1] 8
```

## Simulation

`create.sim()` produces a simulated FITS image and catalogue from a LAMBDAR
parameter file. It accepts observing parameters such as exposure time,
collecting area, effective wavelength, and filter width:

```r
create.sim(
  "Lambdar_Simulation.par",
  ObsParm = data.frame(exp = 53.9, area = 4.91, lamEff = 1111.2, Weff = 6165)
)
```

See `?measure.fluxes`, `?create.par.file`, and `?create.sim` for the full
parameter reference and complete examples.

## License

GPL-2. See [`DESCRIPTION`](DESCRIPTION).
