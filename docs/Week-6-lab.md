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
##       [,1]     [,2]      [,3]      [,4]      [,5]      [,6]      [,7]      [,8]
##  [1,]   NA 0.562243 0.1570045 0.8939159 0.2001269 0.8679510 0.4576649 0.6892632
##  [2,]   NA       NA 0.4201393 0.5677350 0.5386224 0.6617155 0.9116861 0.8241854
##  [3,]   NA       NA        NA 0.2326374 0.7831910 0.1933503 0.4463640 0.2717204
##  [4,]   NA       NA        NA        NA 0.2896541 0.7963900 0.4962321 0.6678212
##  [5,]   NA       NA        NA        NA        NA 0.2483737 0.5803233 0.3528444
##  [6,]   NA       NA        NA        NA        NA        NA 0.5512613 0.8107414
##  [7,]   NA       NA        NA        NA        NA        NA        NA 0.7156467
##  [8,]   NA       NA        NA        NA        NA        NA        NA        NA
##  [9,]   NA       NA        NA        NA        NA        NA        NA        NA
## [10,]   NA       NA        NA        NA        NA        NA        NA        NA
##            [,9]     [,10]
##  [1,] 0.2154065 0.3134618
##  [2,] 0.5512992 0.6910959
##  [3,] 0.7909162 0.6657463
##  [4,] 0.2975460 0.3751635
##  [5,] 0.9997327 0.8438656
##  [6,] 0.2658479 0.3806858
##  [7,] 0.5948436 0.7519638
##  [8,] 0.3713651 0.5080340
##  [9,]        NA 0.8490839
## [10,]        NA        NA
```

Now we can see how many of these p.values are "significant". We know these are false positives, because all the data were generated from the same distribution.


```r
false.positives<-sum(p.values<0.05,na.rm=T)
false.positives
```

```
## [1] 0
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
##        5% 
## 0.0014825
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
