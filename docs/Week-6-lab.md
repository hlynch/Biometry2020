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
##       [,1]       [,2]      [,3]      [,4]       [,5]       [,6]       [,7]
##  [1,]   NA 0.08941968 0.5691592 0.1036412 0.02043682 0.01528556 0.14324170
##  [2,]   NA         NA 0.1491839 0.7629038 0.53716120 0.69772625 0.38422637
##  [3,]   NA         NA        NA 0.1766957 0.03290488 0.02040242 0.24106202
##  [4,]   NA         NA        NA        NA 0.32362750 0.41518906 0.51989273
##  [5,]   NA         NA        NA        NA         NA 0.74246620 0.10060366
##  [6,]   NA         NA        NA        NA         NA         NA 0.08271415
##  [7,]   NA         NA        NA        NA         NA         NA         NA
##  [8,]   NA         NA        NA        NA         NA         NA         NA
##  [9,]   NA         NA        NA        NA         NA         NA         NA
## [10,]   NA         NA        NA        NA         NA         NA         NA
##             [,8]       [,9]      [,10]
##  [1,] 0.05308700 0.68704062 0.47373679
##  [2,] 0.92789484 0.17642687 0.20504300
##  [3,] 0.08510151 0.93188646 0.82931981
##  [4,] 0.79622719 0.22301262 0.25862679
##  [5,] 0.42127745 0.04714574 0.05020011
##  [6,] 0.55513954 0.04432155 0.03972556
##  [7,] 0.30261992 0.35076481 0.41174391
##  [8,]         NA 0.13120071 0.14099193
##  [9,]         NA         NA 0.80033184
## [10,]         NA         NA         NA
```

Now we can see how many of these p.values are "significant". We know these are false positives, because all the data were generated from the same distribution.


```r
false.positives<-sum(p.values<0.05,na.rm=T)
false.positives
```

```
## [1] 7
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
## 0.001907274
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
