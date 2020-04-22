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
##       [,1]     [,2]      [,3]      [,4]      [,5]      [,6]      [,7]      [,8]
##  [1,]   NA 0.874407 0.4165237 0.5014703 0.6187876 0.8706437 0.6015672 0.8618122
##  [2,]   NA       NA 0.3297883 0.5938027 0.4838801 0.7466488 0.4618071 0.9717632
##  [3,]   NA       NA        NA 0.1512570 0.6498919 0.5132870 0.6437475 0.3594709
##  [4,]   NA       NA        NA        NA 0.1870497 0.4074507 0.1640572 0.6730343
##  [5,]   NA       NA        NA        NA        NA 0.7613682 0.9949823 0.5230915
##  [6,]   NA       NA        NA        NA        NA        NA 0.7489910 0.7481666
##  [7,]   NA       NA        NA        NA        NA        NA        NA 0.5074382
##  [8,]   NA       NA        NA        NA        NA        NA        NA        NA
##  [9,]   NA       NA        NA        NA        NA        NA        NA        NA
## [10,]   NA       NA        NA        NA        NA        NA        NA        NA
##             [,9]     [,10]
##  [1,] 0.15604515 0.4694588
##  [2,] 0.11114211 0.3441944
##  [3,] 0.56951100 0.7756977
##  [4,] 0.04022353 0.1080605
##  [5,] 0.26024349 0.8052754
##  [6,] 0.20992575 0.6051920
##  [7,] 0.24890015 0.7987505
##  [8,] 0.13900349 0.3999999
##  [9,]         NA 0.3270434
## [10,]         NA        NA
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
## 0.002051603
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
