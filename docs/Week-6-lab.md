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
##       [,1]      [,2]      [,3]      [,4]       [,5]      [,6]      [,7]
##  [1,]   NA 0.9241858 0.3991918 0.3610815 0.44931216 0.8493008 0.6825321
##  [2,]   NA        NA 0.3917830 0.3215554 0.32865141 0.9120191 0.7129262
##  [3,]   NA        NA        NA 0.9089897 0.08659679 0.4562908 0.5092367
##  [4,]   NA        NA        NA        NA 0.04080235 0.4002640 0.4013512
##  [5,]   NA        NA        NA        NA         NA 0.2853695 0.1382617
##  [6,]   NA        NA        NA        NA         NA        NA 0.8177154
##  [7,]   NA        NA        NA        NA         NA        NA        NA
##  [8,]   NA        NA        NA        NA         NA        NA        NA
##  [9,]   NA        NA        NA        NA         NA        NA        NA
## [10,]   NA        NA        NA        NA         NA        NA        NA
##            [,8]      [,9]     [,10]
##  [1,] 0.9903510 0.8909226 0.2680673
##  [2,] 0.8975577 0.9638710 0.2233531
##  [3,] 0.3165321 0.3952444 0.8871591
##  [4,] 0.2271254 0.3113674 0.7012903
##  [5,] 0.3718764 0.2869577 0.0269793
##  [6,] 0.8091398 0.9433219 0.2838182
##  [7,] 0.5857554 0.7341510 0.2575689
##  [8,]        NA 0.8548518 0.1520997
##  [9,]        NA        NA 0.2110243
## [10,]        NA        NA        NA
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
## 0.001309026
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
