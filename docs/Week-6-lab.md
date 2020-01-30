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
##       [,1]      [,2]      [,3]       [,4]      [,5]      [,6]       [,7]
##  [1,]   NA 0.3849308 0.9704138 0.33200210 0.6133343 0.4122114 0.38308730
##  [2,]   NA        NA 0.4415239 0.08831059 0.1788280 0.1163827 0.06626621
##  [3,]   NA        NA        NA 0.41904852 0.6969883 0.4991969 0.51943005
##  [4,]   NA        NA        NA         NA 0.6053750 0.8823457 0.70034584
##  [5,]   NA        NA        NA         NA        NA 0.7182635 0.79718970
##  [6,]   NA        NA        NA         NA        NA        NA 0.84160204
##  [7,]   NA        NA        NA         NA        NA        NA         NA
##  [8,]   NA        NA        NA         NA        NA        NA         NA
##  [9,]   NA        NA        NA         NA        NA        NA         NA
## [10,]   NA        NA        NA         NA        NA        NA         NA
##            [,8]      [,9]      [,10]
##  [1,] 0.7376083 0.5230689 0.12167977
##  [2,] 0.3416987 0.7990741 0.01872876
##  [3,] 0.7821006 0.5673870 0.21576457
##  [4,] 0.6759653 0.1269892 0.70809391
##  [5,] 0.9763608 0.2564470 0.30534228
##  [6,] 0.7644202 0.1659325 0.58413051
##  [7,] 0.8428348 0.1035772 0.29773446
##  [8,]        NA 0.4317111 0.45977034
##  [9,]        NA        NA 0.02887969
## [10,]        NA        NA         NA
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
## 0.001933103
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
