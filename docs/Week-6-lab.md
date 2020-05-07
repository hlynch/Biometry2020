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
##       [,1]      [,2]       [,3]      [,4]      [,5]       [,6]       [,7]
##  [1,]   NA 0.0426775 0.05717604 0.1678708 0.0437653 0.05630844 0.03203463
##  [2,]   NA        NA 0.68383182 0.5260805 0.6020798 0.94197129 0.58415308
##  [3,]   NA        NA         NA 0.3979771 0.9316721 0.65431090 0.94421482
##  [4,]   NA        NA         NA        NA 0.3384523 0.58708282 0.31056108
##  [5,]   NA        NA         NA        NA        NA 0.57604414 0.98294232
##  [6,]   NA        NA         NA        NA        NA         NA 0.55788060
##  [7,]   NA        NA         NA        NA        NA         NA         NA
##  [8,]   NA        NA         NA        NA        NA         NA         NA
##  [9,]   NA        NA         NA        NA        NA         NA         NA
## [10,]   NA        NA         NA        NA        NA         NA         NA
##            [,8]      [,9]     [,10]
##  [1,] 0.1249838 0.1940988 0.1947497
##  [2,] 0.8517739 0.9848350 0.6652410
##  [3,] 0.8916447 0.7632937 0.4918548
##  [4,] 0.5605973 0.7127715 0.9250705
##  [5,] 0.8306385 0.7061370 0.4322205
##  [6,] 0.8200622 0.9813290 0.7151528
##  [7,] 0.8372113 0.7073141 0.4147073
##  [8,]        NA 0.8755123 0.6385504
##  [9,]        NA        NA 0.7827069
## [10,]        NA        NA        NA
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
## 0.001837348
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
