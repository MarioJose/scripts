# Executing instructions in multi cores in R

#### Mario Jose Marques-Azevedo

Here we will see how to execute more than one [instruction](https://en.wikipedia.org/wiki/Process_%28computing%29) simultaneously on each [core](https://en.wikipedia.org/wiki/Multi-core_processor) of our processor. When we run a command or instruction, in short, this is translated to machine language and all needed calculus is processed by the processor core. Modern processors are assembled by more than one core, or core with [multithread](https://en.wikipedia.org/wiki/Hyper-threading), to execute more than one task simultaneously. This is a very little snapshot of the things.

Then, if my processor has more than one core, or thread, can I run a slow or complicated calculus more fast? More or less! A calculus is one instruction and, in general, you can not split a calculus in part to send them for each core. But, if you have to execute the same calculus more than one time, you can send each of this time to one core and run the all process more fast.

Let's see it in practice. First we need to load `parallel` package. This package is in default R installation.

```r
# Load package
library("parallel")
```

To detect the number of cores of our processor we can execute:

```r
# Check number of cores
nc <- detectCores()
nc
```

```
[1] 4
```

We have four cores. Then we can virtually run our instruction four times more fast.

Suppose that you have a function `fn` that execute some calculus and you need to run this function 1000 times and get the mean of this 1000 results. This function will get a vector of sampled abundance, for instance, re-sample it and calculate Shannon diversity index.

```r
# Function to re-sample abundance with informed size and calculate Shannon index
fn <- function(x, size){
  expand_abd <- rep(x, times=x)
  s <- table(sample(expand_abd, size = size, replace = TRUE))
  p <- s / sum(s)
  H <- - sum(p * log(p))
  return(H)
}
```

We will creating a sampled abundance. For this we will to use `set.seed` function for the results be the same of your execution of this tutorial.

```r
# Creating a abundance vector
set.seed(42)
abd <- round(rweibull(100, 1.5, 40), 0)
abd
```

```
  [1]   8   6  46  13  23  30  18  64  22  20  34  19   7  49  34   6   3  66  33  28   9  63
 [23]   2   6  74  30  38   9  35  13  18  14  39  21 125  13 116  54   9  25  39  35  88   4
 [45]  36   5  10  23   4  25  43  42  38  16  88  17  21  58  49  30  21   3  17  27  12  56
 [67]  48  13  20  51  86  63  53  33  55  19 114  39  30 139  27  60  41  23  16  28  51  72
 [89]  73  45  22 164  54   7   7  18  43  30  18  25
```

Running our function with abundance to calculate Shannon index for re-samples of size 80, for instance.

```r
# Calculating the Shannon index of re-sampled abundance
set.seed(42)
fn(abd, 80)
```

```
[1] 3.458882
```

We can execute this function 1000 times using `for` function:

```r
# Calculating the Shannon index of the re-sample abundance 1000 times
set.seed(42)
results <- c()
for(i in 1:1000){ results[i] <- fn(abd, 80) }
mean(results)
```

```
[1] 3.459867
```

Now we can check how long this 1000 re-samples takes? First we will create a function  with all this commands. This function will be used after to run in cluster.

```r
fn2 <- function(rand, abd, size){
  results <- c()
  for(i in 1:rand){ results[i] <- fn(abd, size) }
  return(mean(results))
}
```

We can see that the result of `fn2` is the same:

```r
set.seed(42)
fn2(rand = 1000, abd = abd, size = 80)
```

```
[1] 3.459867
```

Now checking how long it takes:

```r
# Checking how long the function takes
set.seed(42)
system.time(fn2(rand = 1000, abd = abd, size = 80))
```

```
   user  system elapsed 
  0.368   0.000   0.369 
```

It takes 0.369 milliseconds.

Now, we will see how to execute this second function in cluster. First we need to create a cluster object informing the number of cores. It is like to run a R environment in each core.

```r
# Creating cluster
cl <- makeCluster(nc)
```

The function `fn` exist only in the current environment. Then, we need to export this function to other virtual environment created for each core.

```r
# Exporting variables to virtual environment created for each core
clusterExport(cl, "fn")
```

Now, we want to execute 1000 times the function `fn`. For this we create the function `fn2` with `rand` parameter to inform how many randomization we want. How we have four cores, then we can run 250 time the function `fn`, through the `fn2` with parameter `rand = 250`, in each core. After that we will executed 1000 times the function `fn` through the `fn2`.

```r
# Running function in each core
set.seed(42)
results2 <- clusterApplyLB(cl = cl, x = c(250, 250, 250, 250), fun = fn2, abd = abd, size = 80)
results2
```

```
[[1]]
[1] 3.471441

[[2]]
[1] 3.45454

[[3]]
[1] 3.458044

[[4]]
[1] 3.458233
```

The `clusterApplyLB` needs a parameter `x` informing the value of `rand` parameter ofn `fn2` to be executed in each core. Then, we inform 250 for each of four cores (sum 1000). The `clusterApplyLB` return a list with four elements, one for each core. This is the mean of 250 results of `fn` executed in each core. We can to calculate the mean of this results:

```r
mean(unlist(results2))
```

To check how long the execute in cluster takes:

```r
set.seed(42)
system.time(clusterApplyLB(cl = cl, x = c(250, 250, 250, 250), fun = fn2, abd=abd, size=80))
```

```
   user  system elapsed 
  0.004   0.000   0.187 
```

It takes 0.187 milliseconds against 0.369 of non cluster execution.
