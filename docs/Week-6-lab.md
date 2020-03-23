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
##       [,1]      [,2]       [,3]      [,4]       [,5]        [,6]      [,7]
##  [1,]   NA 0.0781242 0.74273144 0.8951424 0.92631169 0.398602331 0.7643041
##  [2,]   NA        NA 0.07513988 0.1261670 0.07629715 0.008767037 0.1381302
##  [3,]   NA        NA         NA 0.6669800 0.67681351 0.700472529 0.5587957
##  [4,]   NA        NA         NA        NA 0.96026802 0.351106649 0.8790090
##  [5,]   NA        NA         NA        NA         NA 0.327385173 0.8250740
##  [6,]   NA        NA         NA        NA         NA          NA 0.2494224
##  [7,]   NA        NA         NA        NA         NA          NA        NA
##  [8,]   NA        NA         NA        NA         NA          NA        NA
##  [9,]   NA        NA         NA        NA         NA          NA        NA
## [10,]   NA        NA         NA        NA         NA          NA        NA
##               [,8]        [,9]     [,10]
##  [1,] 0.1097414721 0.145854416 0.9831333
##  [2,] 0.0006256912 0.000914067 0.1696257
##  [3,] 0.3138918948 0.380788598 0.7973829
##  [4,] 0.1017078998 0.133233997 0.8986386
##  [5,] 0.0747195456 0.102416498 0.9250608
##  [6,] 0.4373831408 0.543832075 0.5200867
##  [7,] 0.0560814992 0.076446077 0.7955395
##  [8,]           NA 0.845268708 0.2259669
##  [9,]           NA          NA 0.2737796
## [10,]           NA          NA        NA
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
## [1] 2
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
## [1] 2
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
##        5% 
## 0.0022674
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 2
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
