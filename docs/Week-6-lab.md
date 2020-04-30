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
##       [,1]      [,2]      [,3]      [,4]      [,5]      [,6]        [,7]
##  [1,]   NA 0.2744116 0.5206603 0.9754929 0.6155222 0.4195247 0.318021026
##  [2,]   NA        NA 0.7240330 0.3670720 0.0479781 0.7664169 0.006583227
##  [3,]   NA        NA        NA 0.5658674 0.2140186 0.9193384 0.077808893
##  [4,]   NA        NA        NA        NA 0.7041328 0.4884097 0.437047815
##  [5,]   NA        NA        NA        NA        NA 0.1237624 0.557093812
##  [6,]   NA        NA        NA        NA        NA        NA 0.030054969
##  [7,]   NA        NA        NA        NA        NA        NA          NA
##  [8,]   NA        NA        NA        NA        NA        NA          NA
##  [9,]   NA        NA        NA        NA        NA        NA          NA
## [10,]   NA        NA        NA        NA        NA        NA          NA
##            [,8]       [,9]     [,10]
##  [1,] 0.9233291 0.29643904 0.8745655
##  [2,] 0.1918788 0.83621128 0.2285712
##  [3,] 0.5119691 0.65718710 0.4389187
##  [4,] 0.9104909 0.35767753 0.9150853
##  [5,] 0.4660795 0.09817188 0.7639462
##  [6,] 0.3762577 0.68997455 0.3487491
##  [7,] 0.1714678 0.03233349 0.4444422
##  [8,]        NA 0.26318276 0.7861924
##  [9,]        NA         NA 0.2488797
## [10,]        NA         NA        NA
```

Now we can see how many of these p.values are "significant". We know these are false positives, because all the data were generated from the same distribution.


```r
false.positives<-sum(p.values<0.05,na.rm=T)
false.positives
```

```
## [1] 4
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
##         5% 
## 0.00210532
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
