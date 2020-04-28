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
##       [,1]      [,2]      [,3]       [,4]      [,5]      [,6]      [,7]
##  [1,]   NA 0.8398241 0.8444179 0.09719167 0.7045398 0.3450280 0.3558726
##  [2,]   NA        NA 0.9905475 0.13138685 0.8385337 0.4529562 0.4619116
##  [3,]   NA        NA        NA 0.11853829 0.8268582 0.4293861 0.4402397
##  [4,]   NA        NA        NA         NA 0.2498909 0.3457580 0.3758734
##  [5,]   NA        NA        NA         NA        NA 0.6647154 0.6637051
##  [6,]   NA        NA        NA         NA        NA        NA 0.9841811
##  [7,]   NA        NA        NA         NA        NA        NA        NA
##  [8,]   NA        NA        NA         NA        NA        NA        NA
##  [9,]   NA        NA        NA         NA        NA        NA        NA
## [10,]   NA        NA        NA         NA        NA        NA        NA
##             [,8]       [,9]      [,10]
##  [1,] 0.91390210 0.86097618 0.54847291
##  [2,] 0.71949172 0.95103779 0.38151272
##  [3,] 0.71833717 0.96016673 0.36641698
##  [4,] 0.03853593 0.07351224 0.01198885
##  [5,] 0.59379339 0.77882531 0.33784437
##  [6,] 0.18203967 0.32642763 0.05231949
##  [7,] 0.20559294 0.34751246 0.06912830
##  [8,]         NA 0.70234099 0.47403128
##  [9,]         NA         NA 0.27864076
## [10,]         NA         NA         NA
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
## 0.001684734
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
