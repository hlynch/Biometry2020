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
##       [,1]      [,2]       [,3]       [,4]      [,5]      [,6]      [,7]
##  [1,]   NA 0.1797745 0.09348012 0.03907419 0.1540715 0.3158670 0.1329998
##  [2,]   NA        NA 0.73577804 0.29139298 0.8873464 0.6922953 0.8479872
##  [3,]   NA        NA         NA 0.39613388 0.5450237 0.4481149 0.8950661
##  [4,]   NA        NA         NA         NA 0.1912115 0.1717624 0.3658186
##  [5,]   NA        NA         NA         NA        NA 0.7273158 0.7025377
##  [6,]   NA        NA         NA         NA        NA        NA 0.5574168
##  [7,]   NA        NA         NA         NA        NA        NA        NA
##  [8,]   NA        NA         NA         NA        NA        NA        NA
##  [9,]   NA        NA         NA         NA        NA        NA        NA
## [10,]   NA        NA         NA         NA        NA        NA        NA
##             [,8]      [,9]     [,10]
##  [1,] 0.04547566 0.3082292 0.1269600
##  [2,] 0.39086058 0.6569169 0.7258572
##  [3,] 0.54297678 0.4024009 0.9361898
##  [4,] 0.75323327 0.1511230 0.5022768
##  [5,] 0.23985430 0.6790090 0.5964511
##  [6,] 0.22167818 0.9785759 0.4851340
##  [7,] 0.49583453 0.5183113 0.8549966
##  [8,]         NA 0.1907083 0.6688448
##  [9,]         NA        NA 0.4529850
## [10,]         NA        NA        NA
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
## 0.001621305
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
