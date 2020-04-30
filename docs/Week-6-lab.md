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
##       [,1]      [,2]       [,3]       [,4]      [,5]       [,6]      [,7]
##  [1,]   NA 0.1465921 0.77257153 0.55059121 0.8343373 0.55113878 0.7600500
##  [2,]   NA        NA 0.04296811 0.06864749 0.3094667 0.02552689 0.1696755
##  [3,]   NA        NA         NA 0.68439705 0.6350998 0.71196652 0.4891795
##  [4,]   NA        NA         NA         NA 0.4668270 0.90528919 0.3667198
##  [5,]   NA        NA         NA         NA        NA 0.46272396 0.9732882
##  [6,]   NA        NA         NA         NA        NA         NA 0.3105018
##  [7,]   NA        NA         NA         NA        NA         NA        NA
##  [8,]   NA        NA         NA         NA        NA         NA        NA
##  [9,]   NA        NA         NA         NA        NA         NA        NA
## [10,]   NA        NA         NA         NA        NA         NA        NA
##            [,8]       [,9]      [,10]
##  [1,] 0.6121629 0.45113128 0.89174636
##  [2,] 0.2751801 0.03323231 0.06039086
##  [3,] 0.3631825 0.57076272 0.85592070
##  [4,] 0.2885187 0.92738124 0.59013761
##  [5,] 0.8313293 0.38432153 0.72996147
##  [6,] 0.2266276 0.80872743 0.58766456
##  [7,] 0.8081810 0.26475106 0.60414532
##  [8,]        NA 0.19929686 0.45640230
##  [9,]        NA         NA 0.47478440
## [10,]        NA         NA         NA
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
## 0.001710507
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
