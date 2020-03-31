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
##       [,1]     [,2]      [,3]      [,4]      [,5]      [,6]       [,7]
##  [1,]   NA 0.727715 0.3858489 0.6284964 0.8925423 0.6015452 0.12869719
##  [2,]   NA       NA 0.5930323 0.8906371 0.6844467 0.4019987 0.26625316
##  [3,]   NA       NA        NA 0.6891895 0.3981144 0.1944349 0.63324225
##  [4,]   NA       NA        NA        NA 0.6049424 0.3368111 0.34175681
##  [5,]   NA       NA        NA        NA        NA 0.8122216 0.17543682
##  [6,]   NA       NA        NA        NA        NA        NA 0.04153545
##  [7,]   NA       NA        NA        NA        NA        NA         NA
##  [8,]   NA       NA        NA        NA        NA        NA         NA
##  [9,]   NA       NA        NA        NA        NA        NA         NA
## [10,]   NA       NA        NA        NA        NA        NA         NA
##            [,8]      [,9]      [,10]
##  [1,] 0.6784417 0.4515762 0.15930878
##  [2,] 0.9533879 0.7190394 0.11609606
##  [3,] 0.6268122 0.8156261 0.06358833
##  [4,] 0.9348083 0.8362441 0.10062615
##  [5,] 0.6460839 0.4690602 0.24019345
##  [6,] 0.3599596 0.1931904 0.23683383
##  [7,] 0.2868756 0.4184691 0.02383300
##  [8,]        NA 0.7623547 0.10643536
##  [9,]        NA        NA 0.06883156
## [10,]        NA        NA         NA
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
##         5% 
## 0.00142458
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
