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
##       [,1]     [,2]      [,3]      [,4]      [,5]       [,6]       [,7]
##  [1,]   NA 0.659118 0.5180533 0.7561279 0.8018007 0.07997576 0.76446051
##  [2,]   NA       NA 0.8345432 0.8439990 0.8059616 0.12760047 0.40451463
##  [3,]   NA       NA        NA 0.6485344 0.6210731 0.14776880 0.27053303
##  [4,]   NA       NA        NA        NA 0.9479414 0.07052512 0.46167404
##  [5,]   NA       NA        NA        NA        NA 0.07191243 0.51659585
##  [6,]   NA       NA        NA        NA        NA         NA 0.02802564
##  [7,]   NA       NA        NA        NA        NA         NA         NA
##  [8,]   NA       NA        NA        NA        NA         NA         NA
##  [9,]   NA       NA        NA        NA        NA         NA         NA
## [10,]   NA       NA        NA        NA        NA         NA         NA
##             [,8]       [,9]      [,10]
##  [1,] 0.99962920 0.66496699 0.95633785
##  [2,] 0.61242552 0.96570632 0.69538656
##  [3,] 0.44828080 0.78048164 0.54808257
##  [4,] 0.71268895 0.86463941 0.79976446
##  [5,] 0.76885936 0.82222736 0.84618033
##  [6,] 0.05095133 0.09706022 0.08264626
##  [7,] 0.73131844 0.38764540 0.71335626
##  [8,]         NA 0.61053815 0.95069823
##  [9,]         NA         NA 0.70342531
## [10,]         NA         NA         NA
```

Now we can see how many of these p.values are "significant". We know these are false positives, because all the data were generated from the same distribution.


```r
false.positives<-sum(p.values<0.05,na.rm=T)
false.positives
```

```
## [1] 1
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
## 0.002032609
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
