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
##  [1,]   NA 0.285911 0.6375798 0.4314452 0.3765037 0.5633350 0.7712001 0.9481587
##  [2,]   NA       NA 0.1895935 0.7879477 0.7570673 0.1121087 0.2732914 0.3121245
##  [3,]   NA       NA        NA 0.2733357 0.2420762 0.9880506 0.8892567 0.7016216
##  [4,]   NA       NA        NA        NA 0.9993298 0.1868842 0.3733102 0.4435718
##  [5,]   NA       NA        NA        NA        NA 0.1444828 0.3443067 0.4027809
##  [6,]   NA       NA        NA        NA        NA        NA 0.8650944 0.6466319
##  [7,]   NA       NA        NA        NA        NA        NA        NA 0.8263064
##  [8,]   NA       NA        NA        NA        NA        NA        NA        NA
##  [9,]   NA       NA        NA        NA        NA        NA        NA        NA
## [10,]   NA       NA        NA        NA        NA        NA        NA        NA
##             [,9]      [,10]
##  [1,] 0.19709959 0.59282017
##  [2,] 0.72607666 0.13311894
##  [3,] 0.13504070 0.99849308
##  [4,] 0.56081663 0.21234843
##  [5,] 0.51730885 0.17220595
##  [6,] 0.07939313 0.98471540
##  [7,] 0.19762788 0.88085882
##  [8,] 0.21947484 0.67053307
##  [9,]         NA 0.09371116
## [10,]         NA         NA
```

Now we can see how many of these p.values are "significant". We know these are false positives, because all the data were generated from the same distribution.


```r
false.positives<-sum(p.values<0.05,na.rm=T)
false.positives
```

```
## [1] 0
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
## 0.001660882
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
