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
##       [,1]       [,2]      [,3]         [,4]      [,5]       [,6]        [,7]
##  [1,]   NA 0.08513512 0.2799025 0.0002308047 0.1210203 0.27922115 0.621743486
##  [2,]   NA         NA 0.5227032 0.1852636491 0.9227462 0.52880806 0.164877225
##  [3,]   NA         NA        NA 0.0386483878 0.6009998 0.99409483 0.474778077
##  [4,]   NA         NA        NA           NA 0.1664585 0.04033686 0.001130798
##  [5,]   NA         NA        NA           NA        NA 0.60705011 0.218206644
##  [6,]   NA         NA        NA           NA        NA         NA 0.472086674
##  [7,]   NA         NA        NA           NA        NA         NA          NA
##  [8,]   NA         NA        NA           NA        NA         NA          NA
##  [9,]   NA         NA        NA           NA        NA         NA          NA
## [10,]   NA         NA        NA           NA        NA         NA          NA
##             [,8]        [,9]      [,10]
##  [1,] 0.09860435 0.306492120 0.03086683
##  [2,] 0.69911803 0.238280488 0.81956761
##  [3,] 0.74724353 0.649482752 0.35936820
##  [4,] 0.04659676 0.001675762 0.23299369
##  [5,] 0.78964938 0.307377198 0.74571918
##  [6,] 0.75452556 0.645061047 0.36529897
##  [7,] 0.22877116 0.679242640 0.07415674
##  [8,]         NA 0.346479557 0.49812907
##  [9,]         NA          NA 0.11375401
## [10,]         NA          NA         NA
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
## [1] 1
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
##          5% 
## 0.002693127
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 3
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
