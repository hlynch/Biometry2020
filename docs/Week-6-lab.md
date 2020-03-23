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
##       [,1]      [,2]       [,3]      [,4]       [,5]      [,6]      [,7]
##  [1,]   NA 0.9590248 0.01546473 0.3125127 0.45327821 0.3417626 0.2044221
##  [2,]   NA        NA 0.02888612 0.3775858 0.53370622 0.3913096 0.2691391
##  [3,]   NA        NA         NA 0.1000010 0.02035348 0.2382457 0.1003851
##  [4,]   NA        NA         NA        NA 0.66668390 0.9042687 0.8366508
##  [5,]   NA        NA         NA        NA         NA 0.6464573 0.4612553
##  [6,]   NA        NA         NA        NA         NA        NA 0.9744585
##  [7,]   NA        NA         NA        NA         NA        NA        NA
##  [8,]   NA        NA         NA        NA         NA        NA        NA
##  [9,]   NA        NA         NA        NA         NA        NA        NA
## [10,]   NA        NA         NA        NA         NA        NA        NA
##             [,8]       [,9]      [,10]
##  [1,] 0.95680469 0.49305497 0.48081454
##  [2,] 0.99993272 0.56432242 0.54876861
##  [3,] 0.01974107 0.03767728 0.05208686
##  [4,] 0.34828850 0.69343226 0.74018374
##  [5,] 0.50187658 0.99486065 0.96206806
##  [6,] 0.37137493 0.66273267 0.69847382
##  [7,] 0.23553304 0.51478304 0.57004048
##  [8,]         NA 0.53844886 0.52400245
##  [9,]         NA         NA 0.96132701
## [10,]         NA         NA         NA
```

Now we can see how many of these p.values are "significant". We know these are false positives, because all the data were generated from the same distribution.


```r
false.positives<-sum(p.values<0.05,na.rm=T)
false.positives
```

```
## [1] 5
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
## 0.001404029
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
