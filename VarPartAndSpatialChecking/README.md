# Variation partitioning and checking for missing environmental factor in spatial component

#### Mario Jose Marques-Azevedo

This text show how to run partition variance and check if spatial component cam be attributed to neutral process, rather than environmental factors that are missing from the analysis. The partition variance was implemented by [Danilo Rafael Mesquita Neves](http://lattes.cnpq.br/7825577355355814) and [Pedro V. Eisenlohr](http://pedroeisenlohr.webnode.com.br/) and following Anderson and Legendre (1999), Blanchet et al. (2008), Bocard and Legendre (2011), Dray et al. (2006), Legendre et al. (2012), Legendre and Gallagher (2001) and Peres-Neto and Legendre (2010). The protocol to check spatial component was implemented by [Mario Jose Marques-Azevedo](https://github.com/mariojose) and following Diniz-Filho et al. (2012).

## Installing packages

Before start with analysis, we need install some packages:

- vegan;
- spdep (spacemakeR dependence);
- ade4 (spacemakeR dependence);
- tripack (spacemakeR dependence);
- packfor;
- spacemakeR.

If packages are not installed, run the following commands:

```r
install.packages("vegan") 
install.packages("spdep") # spacemakeR dependence
install.packages("ade4") # spacemakeR dependence
install.packages("tripack") # spacemakeR dependence

install.packages("packfor", repos = "http://R-Forge.R-project.org")
install.packages("spacemakeR", repos = "http://R-Forge.R-project.org")
```

After installing the packages we need load them.

```r
library(vegan)
library(packfor)
library(spacemakeR)
library(spdep)
```

Do not forget to set your working directory with `setwd` command.

## Variation partitioning

Before starting the analysis, we need load and check structure of database. We will use three variables names to data:

- `ll`: coordinate data with sites in rows and longitude and latitude in columns;
- `spp`: community data with sites in rows and species in columns;
- `env`: environmental data with sites in rows and variables in columns.

Despite it is not necessary to remove collinear variables prior to variation partitioning, you can eliminate collinear or otherwise superfluous variables if you are interested to obtain a parsimonious model.

We will use Barro Colorado Island (`BCI`) database from Condit et al. (2002) and Zanne et al. (2014) available in `vegan` package as example.

```r
# Load data from vegan package
data(BCI)
BCI.env <- read.delim ('http://www.davidzeleny.net/anadat-r/lib/exe/fetch.php?media=data:bci.env.txt')

# Using our variable names to data
ll <- BCI.env[ ,2:3]
spp <- BCI
env <- BCI.env[ ,-(1:3)]

# Naming rows with sites names
rownames(ll) <- rownames(spp) <- rownames(env) <- paste("site", 1:dim(ll)[1], sep="_")

# Removing original data from dashboard
rm(BCI, BCI.env)
```

To load your data, you can use `read.table` command like that:

```r
# Coordinate data (rows = sites; columns = long, lat)
ll <- read.table(file.choose(), row.names = 1, header = T, sep = ",")

# Community data (rows = sites; columns = species)
spp <- read.table(file.choose(), row.names =  1, header = T, sep = ",")

# Environmental data (rows = sites; columns = variables)
env <- read.table(file.choose(), row.names = 1, header = T, sep = ",")
```

Check if your database file is using comma (",") as data separator character. If you is using other, like semicolon (";") or tabulation ("\\t"), change "sep" parameter of `read.table` command for your case.

Database files (for species abundance, environmental variables and coordinates data) must have the same rows names. For the three files, rows names are sites names of samples. Check your data using:

```r
View(ll)
View(spp)
View(env)
```

We need to prepare the data to analysis. For that, we need remove unicates (species that occur in one site) and standardize species and environmental data.

```r
# Removing unicates and Hellinger transformation
spp.std <- decostand(spp[ ,apply(spp > 0, 2, sum) > 1], "hellinger")

# Standardization of environmental data
env.std <- decostand(env, "standardize")
```

We will generate and test the significance of spatial variables (Moran's Eigenvector Maps - MEMs) to use them in our analysis.

```r
# Convert coordinates to neighbours list. See Borcard et al. (2011) for details.
ll.nb <- tri2nb(ll)

# Spatial weights to neighbour list. See Borcard et al. (2011) for details.
ll.wgh <- nb2listw(ll.nb, style = "B")

# Compute Moran's eigenvectors
mem <- scores.listw(ll.wgh)
```

To make data visualization easier, we will naming columns and rows of generated MEMs vectors. 

```r
# Name columns
colnames(mem$vectors) <- paste("MEM_", 1:ncol(mem$vectors), sep = "")

# Name rows
rownames(mem$vectors) <- rownames(ll)
```

Now, we can compute and test Moran's I for MEMs and select only MEMs that are significant at 0.05 confidence level.

```r
# Compute and test Moran's I for eigenvectors of spatial weighting matrices
mem.test <- test.scores(mem, ll.wgh, 999)

# Select significant MEMs at 0.05 confidence level
mem.sel <- list(values = mem$values[mem.test[ ,"pval"] < 0.05], 
                vectors = mem$vectors[ ,mem.test[ ,"pval"] < 0.05])

# Checking how many MEMs was selected (rows, columns)
dim(mem.sel$vectors)
```

```
[1] 50 39
```

Finishing data to variation partitioning, we will select significant MEMs for our species data using forward selection.

```r
# Selecting MEMs via forward selection
rda.spa <- rda(spp.std, mem.sel$vectors)
rda.spa.sel <- forward.sel(spp.std, mem.sel$vectors, adjR2thresh = RsquareAdj(rda.spa)$adj.r.squared)

#Checking selected MEMs for species data
rda.spa.sel
```

```
  variables order         R2      R2Cum   AdjR2Cum        F  pval
1     MEM_1     1 0.08306955 0.08306955 0.06396683 4.348572 0.001
2     MEM_7     7 0.05784531 0.14091486 0.10435804 3.164680 0.001
3     MEM_3     3 0.04802067 0.18893552 0.13604001 2.723520 0.001
4     MEM_5     5 0.04517433 0.23410986 0.16603073 2.654225 0.001
5     MEM_4     4 0.04403388 0.27814374 0.19611462 2.684040 0.001
6     MEM_2     2 0.04203541 0.32017915 0.22532042 2.658822 0.002
7     MEM_9     9 0.03778533 0.35796448 0.25095856 2.471801 0.001
```

The same procedure for environmental data.

```r
# Selecting environmental variables via forward selection
rda.env <- rda(spp.std, env.std)
rda.env.sel <- forward.sel(spp.std, env.std, adjR2thresh = RsquareAdj(rda.env)$adj.r.squared)

#Checking selected MEMs for environmental data
rda.env.sel
```

```
  variables order         R2     R2Cum   AdjR2Cum        F  pval
1     slope     3 0.10651196 0.1065120 0.08789763 5.722040 0.001
2 elevation     1 0.05324654 0.1597585 0.12400354 2.978414 0.001
3    convex     2 0.05551453 0.2152730 0.16409519 3.254213 0.001
4  aspectEW     4 0.03093983 0.2462129 0.17920957 1.847063 0.008
```

We will create a `spatial` variable with only MEMs selected by forwarding select procedure. The same will be done to environmental variables selected, `environment` variable. This variables will be used in variation partitioning.

```r
# Preparing selected variables for variance partitioning

# Filter MEMs from significant selected MEMs to MEMs selected in RDA
spatial <- list(values = mem.sel$values[colnames(mem.sel$vectors) %in% rda.spa.sel$variables],
                vectors = mem.sel$vectors[ ,colnames(mem.sel$vectors) %in% rda.spa.sel$variables])

# Checking selected MEMs
dim(spatial$vectors)
```

```
[1] 50  7
```

```r
# Filter environment data to selected environmental variables in RDA
environment <- env.std[ ,rda.env.sel$variables]

# Checking selected environmental variables
dim(environment)
```

```
[1] 50  4
```

Finally we will run variation partitioning.

```r
# Variance partitioning
all.varpart <- varpart(spp.std, environment, spatial$vectors)
all.varpart
```

```
Partition of variation in RDA

Call: varpart(Y = spp.std, X = environment, spatial$vectors)

Explanatory tables:
X1:  environment
X2:  spatial$vectors 

No. of explanatory tables: 2 
Total variation (SS): 12.203 
            Variance: 0.24904 
No. of observations: 50 

Partition table:
                     Df R.squared Adj.R.squared Testable
[a+b] = X1            4   0.24621       0.17921     TRUE
[b+c] = X2            7   0.35796       0.25096     TRUE
[a+b+c] = X1+X2      11   0.48340       0.33386     TRUE
Individual fractions                                    
[a] = X1|X2           4                 0.08290     TRUE
[b]                   0                 0.09631    FALSE
[c] = X2|X1           7                 0.15465     TRUE
[d] = Residuals                         0.66614    FALSE
---
Use function 'rda' to test significance of fractions of interest
```

```r
plot(all.varpart, Xnames = c("environment", "spatial"))
```

<div style="text-align: center">
  <img src="varpart.png" alt="Variation partitioning" />
</div>
<br />

Now we can testing environmental and spatial components significance.

```r
# Testing the environmental significance, after considering the effect of selected MEMs
rda.env.spa <- rda(spp.std, environment, spatial$vectors)
env.aov <- anova(rda.env.spa)
env.aov
```

```
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(X = spp.std, Y = environment, Z = spatial$vectors)
         Df Variance      F Pr(>F)    
Model     4 0.031238 2.3067  0.001 ***
Residual 38 0.128652                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

```r
# Testing the spatial significance, after considering the effect of selected 
#  environmental variables
rda.spa.env <- rda(spp.std, spatial$vectors, environment) 
spa.aov <- anova(rda.spa.env)
spa.aov
```

```
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(X = spp.std, Y = spatial$vectors, Z = environment)
         Df Variance      F Pr(>F)    
Model     7 0.059068 2.4924  0.001 ***
Residual 38 0.128652                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

## Checking for missing environmental factor in spatial component

Diniz et al. (2012) proposed a protocol to test if amount of variation predicted by pure spatial component [c] can be safely attributed to neutral process, rather than lack of spatially structured environmental predictor. If matrix correlation of species abundance (R) hasn't correlation with matrix of species spatial maps (M), then spatial component [c] can safely attributed do neutral process.

This protocol is valid only if you species data contain species abundance. For occurrence data this protocol do not work. 

Before start the analysis we need the download function `correlogI`. You can download the [`correlogI`](https://github.com/MarioJose/r-functions/tree/master/correlogI) function at [this link](https://raw.githubusercontent.com/MarioJose/r-functions/master/correlogI/correlogI.r). Save the function script at your work directory. Read [`correlogI`](https://github.com/MarioJose/r-functions/tree/master/correlogI) function page for more information.

```r
# Load correlogram function.
source("correlogI.r")
```

We will define the number of class for our correlogram using Sturges's rule

```r
# Number of classes: Sturges' rule
nc <- round(1 + (1/log10(2)) * log10(length(dist(ll)[upper.tri(dist(ll), diag = FALSE)])), 0)
```

We will predict species abundance from pure spatial variation [c]. For this we predict abundance from a full model [a + b + c] and subtract from abundance predicted from environmental model [a + b]. After that, we will create a correlation matrix (R) of abundance predicted by pure spatial component.

```r
# Predicting species abundance from pure spatial variation [c]
pred.spa <- predict(rda(spp.std, cbind(environment, spatial$vectors))) - predict(rda(spp.std, environment))

# Matrix with species correlation
R <- cor(pred.spa)
```

Now, we will calculate correlogram for each species.

```r
# Matrix with species correlogram
corr <- matrix(NA, nrow = nc, ncol = dim(pred.spa)[2], dimnames = list(1:nc, colnames(pred.spa)))

for(i in 1:dim(corr)[2]){
  corr[ ,i] <- correlogI(pred.spa[ ,i], ll)$result$I
} 
rm(i)
```

Then, we can create a matrix of pairwise Manhattan similarity (M) among species correlogram.

```r
# Matrix of pairwise similarity among correlogram (Manhattan distance)
M <- matrix(NA, nrow = dim(pred.spa)[2], dim(pred.spa)[2], dimnames = list(names(pred.spa), names(pred.spa)))

# 3D array with x and y to species and z to distance classes
tmp <- array(NA, c(dim(pred.spa)[2],dim(pred.spa)[2],nc))

for(k in 1:nc){
  tmat <- matrix(rep(corr[k, ], times = dim(pred.spa)[2]), nrow = dim(pred.spa)[2], ncol = dim(pred.spa)[2])
  tmp[, , k] <- abs(tmat - t(tmat)) / nc
}

# Sum distance classes difference (z axis)
M[] <- apply(tmp, c(1,2), sum, na.rm = TRUE)

rm(k, tmp)
```

Finally, we can check if neutral process can be safely attributed to spatial component.

```r
# Checking with R and M are correlated
mantel(M, R)
```

```
Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel(xdis = M, ydis = R) 

Mantel statistic r: -0.02087 
      Significance: 0.971 

Upper quantiles of permutations (null model):
   90%    95%  97.5%    99% 
0.0129 0.0167 0.0188 0.0214 
Permutation: free
Number of permutations: 999
```

If we have no correlation (Matel correlation) among M and R matrix, we can safely attribute neutral process to spatial component.

# Bibliography

Anderson, M. J. & P. Legendre. 1999. An empirical comparison of permutation methods for tests of partial regression coefficients in a linear model. [Journal of Statistical Computation and Simulation 62: 271-303](http://www.tandfonline.com/doi/abs/10.1080/00949659908811936).

Blanchet F. G., P. Legendre, and D. Borcard. 2008. Forward selection of explanatory variables. [Ecology 89: 2623-2632](http://www.esajournals.org/doi/abs/10.1890/07-0986.1).

Borcard, D., F. Gillet & P. Legendre. 2011. Numerical Ecology with R. [Springer, New York, 302p](http://www.springer.com/us/book/9781441979759).

Diniz-Filho, J. A. F. et al. 2012. Spatial autocorrelation analysis allows disentangling the balance between neutral and niche processes in metacommunities. [Oikos 121: 201–210](http://onlinelibrary.wiley.com/doi/10.1111/j.1600-0706.2011.19563.x/abstract).

Dray, S., P. Legendre and P. Peres-Neto. 2006. Spatial modelling: a comprehensive framework for principal coordinate analysis of neighbor matrices (PCNM). [Ecological Modelling 196: 483-493](http://www.sciencedirect.com/science/article/pii/S0304380006000925).

Legendre, P., D. Borcard and D. W. Roberts. 2012. Variation partitioning involving orthogonal spatial eigenfunction submodels. [Ecology 93: 1234-1240](http://www.esajournals.org/doi/abs/10.1890/11-2028.1).

Legendre, P. and E. Gallagher. 2001. Ecologically meaningful transformations for ordination of species data. [Oecologia 129: 271-280](http://link.springer.com/article/10.1007/s004420100716).

Peres-Neto, P. R. and P. Legendre. 2010. Estimating and controlling for spatial structure in the study of ecological communities. [Global Ecology and Biogeography 19: 174-184](http://onlinelibrary.wiley.com/doi/10.1111/j.1466-8238.2009.00506.x/abstract).

# Data bibliography

Condit, R, Pitman, N, Leigh, E.G., Chave, J., Terborgh, J., Foster, R.B., Nuñez, P., Aguilar, S., Valencia, R., Villa, G., Muller-Landau, H.C., Losos, E. & Hubbell, S.P. (2002). Beta-diversity in tropical forest trees. [Science 295, 666–669](http://www.sciencemag.org/content/295/5555/666).

Zanne A.E., Tank D.C., Cornwell, W.K., Eastman J.M., Smith, S.A., FitzJohn, R.G., McGlinn, D.J., O’Meara, B.C., Moles, A.T., Reich, P.B., Royer, D.L., Soltis, D.E., Stevens, P.F., Westoby, M., Wright, I.J., Aarssen, L., Bertin, R.I., Calaminus, A., Govaerts, R., Hemmings, F., Leishman, M.R., Oleksyn, J., Soltis, P.S., Swenson, N.G., Warman, L. & Beaulieu, J.M. (2014) Three keys to the radiation of angiosperms into freezing environments. [Nature 506, 89–92](http://www.nature.com/nature/journal/v506/n7486/full/nature12872.html), doi:10.1038/nature12872.
