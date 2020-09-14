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
##       [,1]      [,2]      [,3]      [,4]      [,5]       [,6]       [,7]
##  [1,]   NA 0.6061354 0.7870907 0.3551829 0.2260448 0.54853508 0.47794898
##  [2,]   NA        NA 0.7510016 0.6160305 0.3979571 0.15594067 0.81672184
##  [3,]   NA        NA        NA 0.4058755 0.2231961 0.24718346 0.57439330
##  [4,]   NA        NA        NA        NA 0.7625149 0.05627899 0.78032096
##  [5,]   NA        NA        NA        NA        NA 0.01565354 0.54076896
##  [6,]   NA        NA        NA        NA        NA         NA 0.09313118
##  [7,]   NA        NA        NA        NA        NA         NA         NA
##  [8,]   NA        NA        NA        NA        NA         NA         NA
##  [9,]   NA        NA        NA        NA        NA         NA         NA
## [10,]   NA        NA        NA        NA        NA         NA         NA
##             [,8]       [,9]     [,10]
##  [1,] 0.43261644 0.74458082 0.8186225
##  [2,] 0.76152198 0.27507595 0.7534045
##  [3,] 0.50832263 0.42287781 0.9776742
##  [4,] 0.80880172 0.10854004 0.4348200
##  [5,] 0.54797063 0.03619238 0.2653002
##  [6,] 0.05939645 0.67943323 0.3229237
##  [7,] 0.95455809 0.17526550 0.5930923
##  [8,]         NA 0.12451415 0.5370544
##  [9,]         NA         NA 0.4995773
## [10,]         NA         NA        NA
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
## 0.001795757
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
