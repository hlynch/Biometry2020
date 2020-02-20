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
##       [,1]       [,2]       [,3]       [,4]      [,5]       [,6]      [,7]
##  [1,]   NA 0.03142679 0.07324319 0.58334044 0.3042353 0.89719896 0.2162605
##  [2,]   NA         NA 0.51904647 0.01555602 0.1632379 0.03349189 0.4882143
##  [3,]   NA         NA         NA 0.03505411 0.3822214 0.07824348 0.8509301
##  [4,]   NA         NA         NA         NA 0.1447291 0.49389663 0.1136429
##  [5,]   NA         NA         NA         NA        NA 0.34080105 0.6298379
##  [6,]   NA         NA         NA         NA        NA         NA 0.2394194
##  [7,]   NA         NA         NA         NA        NA         NA        NA
##  [8,]   NA         NA         NA         NA        NA         NA        NA
##  [9,]   NA         NA         NA         NA        NA         NA        NA
## [10,]   NA         NA         NA         NA        NA         NA        NA
##             [,8]       [,9]      [,10]
##  [1,] 0.39098039 0.70822867 0.75783546
##  [2,] 0.08871733 0.01164934 0.08118244
##  [3,] 0.21836103 0.02392675 0.18478357
##  [4,] 0.18199958 0.79801729 0.43455135
##  [5,] 0.77072705 0.13659902 0.54400066
##  [6,] 0.44150176 0.59061055 0.83676673
##  [7,] 0.47704191 0.12141730 0.35877976
##  [8,]         NA 0.16831456 0.68006345
##  [9,]         NA         NA 0.51224182
## [10,]         NA         NA         NA
```

Now we can see how many of these p.values are "significant". We know these are false positives, because all the data were generated from the same distribution.


```r
false.positives<-sum(p.values<0.05,na.rm=T)
false.positives
```

```
## [1] 6
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
## 0.001989079
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
