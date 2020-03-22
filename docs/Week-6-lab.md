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
##       [,1]      [,2]      [,3]       [,4]       [,5]       [,6]      [,7]
##  [1,]   NA 0.4418442 0.5226317 0.14758280 0.18351401 0.39785340 0.7653187
##  [2,]   NA        NA 0.2128664 0.04447491 0.54802569 0.96267854 0.6467665
##  [3,]   NA        NA        NA 0.61324528 0.08797988 0.18751276 0.3896532
##  [4,]   NA        NA        NA         NA 0.01705245 0.03250211 0.1052154
##  [5,]   NA        NA        NA         NA         NA 0.56706669 0.3024326
##  [6,]   NA        NA        NA         NA         NA         NA 0.6030946
##  [7,]   NA        NA        NA         NA         NA         NA        NA
##  [8,]   NA        NA        NA         NA         NA         NA        NA
##  [9,]   NA        NA        NA         NA         NA         NA        NA
## [10,]   NA        NA        NA         NA         NA         NA        NA
##               [,8]       [,9]        [,10]
##  [1,] 0.0366422437 0.96862989 1.211869e-02
##  [2,] 0.2713292459 0.44930918 1.358575e-01
##  [3,] 0.0210753440 0.56898770 9.061772e-03
##  [4,] 0.0003116687 0.20558686 4.909146e-05
##  [5,] 0.7422292580 0.19539654 4.897307e-01
##  [6,] 0.2750550104 0.40898265 1.329668e-01
##  [7,] 0.0996931003 0.75260798 4.152423e-02
##  [8,]           NA 0.05172595 6.039662e-01
##  [9,]           NA         NA 2.049962e-02
## [10,]           NA         NA           NA
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
## [1] 2
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
## [1] 2
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
## 0.00185616
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 2
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
