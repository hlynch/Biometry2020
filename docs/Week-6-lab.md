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
##       [,1]      [,2]      [,3]       [,4]       [,5]      [,6]      [,7]
##  [1,]   NA 0.4889276 0.1169553 0.16248366 0.11367804 0.8376771 0.9767027
##  [2,]   NA        NA 0.3911329 0.08043810 0.46480534 0.6651530 0.5909423
##  [3,]   NA        NA        NA 0.02124393 0.81661687 0.2093530 0.2178840
##  [4,]   NA        NA        NA         NA 0.02205202 0.1505764 0.2548106
##  [5,]   NA        NA        NA         NA         NA 0.2359242 0.2517949
##  [6,]   NA        NA        NA         NA         NA        NA 0.8592987
##  [7,]   NA        NA        NA         NA         NA        NA        NA
##  [8,]   NA        NA        NA         NA         NA        NA        NA
##  [9,]   NA        NA        NA         NA         NA        NA        NA
## [10,]   NA        NA        NA         NA         NA        NA        NA
##            [,8]       [,9]     [,10]
##  [1,] 0.9271397 0.42442366 0.5215898
##  [2,] 0.5282700 0.98765484 0.9053410
##  [3,] 0.1727445 0.36119215 0.3016993
##  [4,] 0.2484457 0.06589004 0.0806624
##  [5,] 0.1955065 0.42581384 0.3487619
##  [6,] 0.8068625 0.62495362 0.7236930
##  [7,] 0.9629749 0.56020391 0.6359234
##  [8,]        NA 0.49060599 0.5680220
##  [9,]        NA         NA 0.8808201
## [10,]        NA         NA        NA
```

Now we can see how many of these p.values are "significant". We know these are false positives, because all the data were generated from the same distribution.


```r
false.positives<-sum(p.values<0.05,na.rm=T)
false.positives
```

```
## [1] 2
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
## 0.001934238
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
