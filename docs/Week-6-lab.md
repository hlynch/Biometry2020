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
##       [,1]       [,2]      [,3]       [,4]       [,5]        [,6]      [,7]
##  [1,]   NA 0.02101172 0.2883720 0.31375819 0.52722608 0.936924381 0.3227857
##  [2,]   NA         NA 0.1363263 0.09404364 0.03968569 0.009125264 0.2202754
##  [3,]   NA         NA        NA 0.90871789 0.58826775 0.245559885 0.9543071
##  [4,]   NA         NA        NA         NA 0.65003362 0.264529450 0.8755873
##  [5,]   NA         NA        NA         NA         NA 0.504789237 0.6048047
##  [6,]   NA         NA        NA         NA         NA          NA 0.2993466
##  [7,]   NA         NA        NA         NA         NA          NA        NA
##  [8,]   NA         NA        NA         NA         NA          NA        NA
##  [9,]   NA         NA        NA         NA         NA          NA        NA
## [10,]   NA         NA        NA         NA         NA          NA        NA
##             [,8]       [,9]      [,10]
##  [1,] 0.25061009 0.80117784 0.70817836
##  [2,] 0.09483964 0.02761469 0.06545072
##  [3,] 0.99267534 0.39222619 0.54033131
##  [4,] 0.90413631 0.42941159 0.58941589
##  [5,] 0.54522668 0.69958538 0.85802958
##  [6,] 0.18938869 0.83419380 0.72752861
##  [7,] 0.94447625 0.42454318 0.55045645
##  [8,]         NA 0.34875903 0.51297850
##  [9,]         NA         NA 0.87902643
## [10,]         NA         NA         NA
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
##          5% 
## 0.002092462
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
