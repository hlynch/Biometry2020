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
##       [,1]       [,2]       [,3]       [,4]       [,5]       [,6]       [,7]
##  [1,]   NA 0.06408012 0.93575444 0.51290532 0.94419780 0.80584253 0.91727251
##  [2,]   NA         NA 0.06413432 0.09851288 0.01098751 0.04303638 0.04921738
##  [3,]   NA         NA         NA 0.47475404 0.97209059 0.74349412 0.85162495
##  [4,]   NA         NA         NA         NA 0.27130473 0.59467123 0.53300435
##  [5,]   NA         NA         NA         NA         NA 0.64198111 0.82147027
##  [6,]   NA         NA         NA         NA         NA         NA 0.87817793
##  [7,]   NA         NA         NA         NA         NA         NA         NA
##  [8,]   NA         NA         NA         NA         NA         NA         NA
##  [9,]   NA         NA         NA         NA         NA         NA         NA
## [10,]   NA         NA         NA         NA         NA         NA         NA
##             [,8]       [,9]       [,10]
##  [1,] 0.80482955 0.94891162 0.047736634
##  [2,] 0.03108641 0.06538743 0.839085811
##  [3,] 0.87674706 0.98703127 0.048326257
##  [4,] 0.32862453 0.48470903 0.069802711
##  [5,] 0.80251662 0.98877014 0.007911721
##  [6,] 0.58226465 0.75711687 0.030401358
##  [7,] 0.70421164 0.86537209 0.035491618
##  [8,]         NA 0.86271257 0.022728811
##  [9,]         NA         NA 0.049234619
## [10,]         NA         NA          NA
```

Now we can see how many of these p.values are "significant". We know these are false positives, because all the data were generated from the same distribution.


```r
false.positives<-sum(p.values<0.05,na.rm=T)
false.positives
```

```
## [1] 11
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
## 0.001629932
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
