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
##       [,1]       [,2]       [,3]       [,4]      [,5]      [,6]       [,7]
##  [1,]   NA 0.04771869 0.73144332 0.50226598 0.6039063 0.4938343 0.20057068
##  [2,]   NA         NA 0.02794328 0.02408972 0.1647819 0.3102730 0.42841549
##  [3,]   NA         NA         NA 0.70541550 0.4212424 0.3550928 0.12237659
##  [4,]   NA         NA         NA         NA 0.2962339 0.2554017 0.09315468
##  [5,]   NA         NA         NA         NA        NA 0.8206854 0.49775114
##  [6,]   NA         NA         NA         NA        NA        NA 0.71841811
##  [7,]   NA         NA         NA         NA        NA        NA         NA
##  [8,]   NA         NA         NA         NA        NA        NA         NA
##  [9,]   NA         NA         NA         NA        NA        NA         NA
## [10,]   NA         NA         NA         NA        NA        NA         NA
##            [,8]      [,9]     [,10]
##  [1,] 0.5001249 0.8399840 0.7249921
##  [2,] 0.2826257 0.1005302 0.1299044
##  [3,] 0.3558087 0.6199902 0.5220542
##  [4,] 0.2546895 0.4367209 0.3670860
##  [5,] 0.8404449 0.7761182 0.8828273
##  [6,] 0.9769434 0.6349969 0.7229124
##  [7,] 0.6855999 0.3359122 0.4110419
##  [8,]        NA 0.6476574 0.7389737
##  [9,]        NA        NA 0.8920462
## [10,]        NA        NA        NA
```

Now we can see how many of these p.values are "significant". We know these are false positives, because all the data were generated from the same distribution.


```r
false.positives<-sum(p.values<0.05,na.rm=T)
false.positives
```

```
## [1] 3
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
## 0.001841243
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
