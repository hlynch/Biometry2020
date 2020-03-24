Week 6 Lab
=============
  
Today we will do a short exercise to illustrate the permutation method of dealing with multiple comparisons.

First, we will simulate 10 groups of data (say, heights) from the *same* distribution using the normal distribution and put the data for each group in a list for easy access:


```r
data.all<-list()
for (i in 1:10)
  {
    data.all[[i]]<-rnorm(10)  #Note the double brackets for a list
  }
```

Now we can compare each group pairwise using a t.test.


```r
p.values<-matrix(ncol=10,nrow=10)
for (i in 1:9)
  {
    for (j in (i+1):10)
      {
        p.values[i,j]<-t.test(data.all[[i]],data.all[[j]])$p.value 
      }
  }
p.values
```

```
##       [,1]      [,2]        [,3]       [,4]        [,5]       [,6]      [,7]
##  [1,]   NA 0.1729018 0.104235308 0.38671535 0.358751613 0.42611831 0.9536496
##  [2,]   NA        NA 0.005885461 0.65742856 0.003565057 0.03120549 0.2157585
##  [3,]   NA        NA          NA 0.01896172 0.247991507 0.35582546 0.1027472
##  [4,]   NA        NA          NA         NA 0.046718272 0.10066897 0.4382449
##  [5,]   NA        NA          NA         NA          NA 0.94091709 0.3473146
##  [6,]   NA        NA          NA         NA          NA         NA 0.4084777
##  [7,]   NA        NA          NA         NA          NA         NA        NA
##  [8,]   NA        NA          NA         NA          NA         NA        NA
##  [9,]   NA        NA          NA         NA          NA         NA        NA
## [10,]   NA        NA          NA         NA          NA         NA        NA
##             [,8]       [,9]      [,10]
##  [1,] 0.71052374 0.74899973 0.58206467
##  [2,] 0.05148685 0.04444874 0.04749968
##  [3,] 0.15364137 0.12714271 0.23492372
##  [4,] 0.18554442 0.18673959 0.15124713
##  [5,] 0.56090308 0.47100497 0.78802709
##  [6,] 0.61494470 0.55470902 0.78490235
##  [7,] 0.67513350 0.71032480 0.55517387
##  [8,]         NA 0.93915627 0.82394382
##  [9,]         NA         NA 0.76206778
## [10,]         NA         NA         NA
```

Now we can see how many of these p.values are "significant". We know these are false positives, because all the data were generated from the same distribution.


```r
false.positives<-sum(p.values<0.05,na.rm=T)
false.positives
```

```
## [1] 7
```

We could correct this using the Bonferonni method:


```r
k<-45
new.threshold.B<-0.05/k
new.threshold.B
```

```
## [1] 0.001111111
```

```r
false.positives.B<-sum(p.values<new.threshold.B,na.rm=T)
false.positives.B
```

```
## [1] 0
```

We could correct this using the Dunn-Sidak method:


```r
k<-45
new.threshold.DS<-1-((1-0.05)^(1/k))
new.threshold.DS
```

```
## [1] 0.001139202
```

```r
false.positives.DS<-sum(p.values<new.threshold.DS,na.rm=T)
false.positives.DS
```

```
## [1] 0
```

We could correct this using the randomization method. This requires simulating data under the null hypothesis to generate a null distribution of p-values.



```r
p.values.all<-c()
min.p.values.all<-c()
for (k in 1:1000)
  {
    data.null<-list()
    for (i in 1:10)
      {
        data.null[[i]]<-rnorm(10)  #Note the double brackets for a list
      }
    p.values.null<-matrix(ncol=10,nrow=10)
    for (i in 1:9)
      {
        for (j in (i+1):10)
          {
            p.values.null[i,j]<-t.test(data.null[[i]],data.null[[j]])$p.value 
          }
      }
    p.values.all<-c(p.values.all,c(p.values.null)[!is.na(c(p.values.null))])
    min.p.values.all<-c(min.p.values.all,min(c(p.values.null)[!is.na(c(p.values.null))]))
  }
new.threshold.R<-quantile(min.p.values.all,probs=c(0.05))
new.threshold.R
```

```
##          5% 
## 0.001562223
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
