Week 2 Lab
=============

Confidence intervals
-----------------------

Before getting too far, we need to circle back and make sure we understand what is meant by a confidence interval. 

A 95th percentile confidence interval say “If I repeat this procedure 100 times using 100 different datasets, 95% of the time my confidence intervals will capture the true parameter”. It does NOT say that there is a 95% chance that the parameter is in the interval.

**Quiz time! (Don't worry, not a real quiz)**

*Important note*: This is an area where Aho is WRONG. I will not repeat Aho's interpretation here because I think he's just wrong. Aho is correct on only one point. It is true that ONCE THE 95th CI HAS BEEN CONSTRUCTED, it is no longer possible to assign a $%$ to the probability that that CI contains the true value or not. Because that CI, once created, either DOES or DOES NOT contain the true value. However, we often talk about the interval in the abstract. When we say "There is a 95$%$ chance that the interval contains the true value" what we mean is that there is a 95$%$ probability that a CI created using that methodology would contain the true value.

Do not let Week 2 pass by without fundamentally understanding the interpretation of a confidence interval. 

Testing hypotheses through permutation
------------------------------------

We'll start off by working through two examples from Phillip Good's book "Introduction to Statistics Through Resampling Methods and R/S-PLUS":

Example #1: Use permutation methods to test the null hypothesis that the treatment does not increase survival time (in other word: $H_{0}$: No difference in survival between the treated and control groups):

Survival.treated=$\{94,197,16,38,99,141,23 \}$

Survival.control=$\{52,104,146,10,51,30,40,27,46 \}$

(Is this a one-tailed or a two-tailed test?)

Make sure that you understand what is being done here, as this example is very closely related to the problem set.


Example #2: Using the same data, provide a 95% confidence interval for the difference in mean survival days based on 1000 bootstrap samples

Note that these two approaches are very closely related. Do you see why either approach can be used to test the null hypothesis? (What is the null hypothesis here?)

Now we will do one slightly more complicated example from Phillip Good's book "Permutation tests: A practical guide to resampling methods and testing hypotheses":

Holmes and Williams (1954) studied tonsil size in children to verify a possible association with the virus \textit{S. pyrogenes}. Test for an association between \textit{S. pyrogenes} status and tonsil size. (Note that you will need to come up with a reasonable test statistic.)

![](Table2categories.png)

Now lets consider the full dataset, where tonsil size is divided into three categories. How would we do the test now? What is the new test statistic? (There are many options.) What 'labels' do you permute?

![](Table3categories.png)

Basics of bootstrap and jackknife
------------------------------------

To get started with bootstrap and jackknife techniques, we start by working through a very simple example. First we simulate some data


```r
x<-seq(0,9,by=1)
```

This will constutute our "data". Let's print the result of sampling with replacement to get a sense for it...


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 1 2 3 4 8 
## 1 2 1 1 2 3
```

Now we will write a little script to take bootstrap samples and calculate the means of each of these bootstrap samples


```r
xmeans<-vector(length=1000)
for (i in 1:1000)
  {
  xmeans[i]<-mean(sample(x,replace=T))
  }
```

The actual number of bootstrapped samples is arbitrary *at this point* but there are ways of characterizing the precision of the bootstrap (jackknife-after-bootstrap) which might inform the number of bootstrap samples needed. *In practice*, people tend to pick some arbitrary but large number of bootstrap samples because computers are so fast that it is often easy to draw far more samples than are actually needed. When calculation of the statistic is slow (as might be the case if you are using the samples to construct a phylogeny, for example), then you would need to be more concerned with the number of bootstrap samples. 

First, lets just look at a histogram of the bootstrapped means and plot the actual sample mean on the histogram for comparison



```r
hist(xmeans,breaks=30,col="pink")
abline(v=mean(x),lwd=2)
```

<img src="Week-2-lab_files/figure-html/unnamed-chunk-4-1.png" width="672" />

Calculating bias and standard error
-----------------------------------

From these we can calculate the bias and standard deviation for the mean (which is the "statistic"):

$$
\widehat{Bias_{boot}} = \left(\frac{1}{k}\sum^{k}_{i=1}\theta^{*}_{i}\right)-\hat{\theta}
$$


```r
bias.boot<-mean(xmeans)-mean(x)
bias.boot
```

```
## [1] -0.0436
```

```r
hist(xmeans,breaks=30,col="pink")
abline(v=mean(x),lwd=5,col="black")
abline(v=mean(xmeans),lwd=2,col="yellow")
```

<img src="Week-2-lab_files/figure-html/unnamed-chunk-5-1.png" width="672" />

$$
\widehat{s.e._{boot}} = \sqrt{\frac{1}{k-1}\sum^{k}_{i=1}(\theta^{*}_{i}-\bar{\theta^{*}})^{2}}
$$


```r
se.boot<-sd(xmeans)
```

We can find the confidence intervals in two ways:

Method #1: Assume the bootstrap statistics are normally distributed


```r
LL.boot<-mean(xmeans)-1.96*se.boot #where did 1.96 come from?
UL.boot<-mean(xmeans)+1.96*se.boot
LL.boot
```

```
## [1] 2.667862
```

```r
UL.boot
```

```
## [1] 6.244938
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.6   6.3
```

Let's compare this to what we would have gotten if we had used normal distribution theory. First we have to calculate the standard error:


```r
se.normal<-sqrt(var(x)/length(x))
LL.normal<-mean(x)-qt(0.975,length(x)-1)*se.normal
UL.normal<-mean(x)+qt(0.975,length(x)-1)*se.normal
LL.normal
```

```
## [1] 2.334149
```

```r
UL.normal
```

```
## [1] 6.665851
```

In this case, the confidence intervals we got from the normal distribution theory are too wide.

Does it make sense why the normal distribution theory intervals are too wide? Because the original were were uniformly distributed, the data has higher variance than would be expected and therefore the standard error is higher than would be expected.

There are two packages that provide functions for bootstrapping, 'boot' and 'boostrap'. We will start by using the 'bootstrap' package, which was originally designed for Efron and Tibshirani's monograph on the bootstrap. 

To test the main functionality of the 'bootstrap' package, we will use the data we already have. The 'bootstrap' function requires the input of a user-defined function to calculate the statistic of interest. Here I will write a function that calculates the mean of the input values.


```r
library(bootstrap)
```

```
## Warning: package 'bootstrap' was built under R version 3.5.2
```

```r
theta<-function(x)
  {
    mean(x)
  }
results<-bootstrap(x=x,nboot=1000,theta=theta)
results
```

```
## $thetastar
##    [1] 5.2 5.7 5.0 3.6 4.6 3.9 5.8 6.4 3.7 4.8 4.8 4.9 4.4 3.4 2.9 5.0 6.2
##   [18] 3.3 3.9 4.3 5.1 5.9 4.3 4.1 3.7 5.3 4.6 3.4 4.1 4.0 5.9 4.9 4.7 4.6
##   [35] 4.2 4.0 3.4 3.3 4.0 3.4 5.2 4.7 4.6 6.2 5.1 5.1 4.8 4.2 4.1 5.0 6.4
##   [52] 4.6 5.0 5.9 4.4 2.8 4.3 6.2 4.1 5.8 5.2 4.8 4.1 3.8 5.2 5.6 4.8 4.0
##   [69] 5.2 3.9 5.3 3.9 4.1 3.1 5.7 3.7 4.6 4.2 5.1 4.1 5.2 3.8 5.3 4.6 5.9
##   [86] 5.9 2.1 3.5 5.2 5.3 4.5 5.0 4.2 6.2 4.0 4.7 4.7 2.8 4.9 5.0 3.6 5.2
##  [103] 4.1 6.0 6.3 4.6 4.7 4.8 2.4 5.1 5.8 4.0 4.1 4.7 4.6 5.5 5.0 4.4 5.0
##  [120] 3.5 4.5 4.9 5.4 2.6 3.6 5.3 3.9 4.7 4.8 4.6 3.4 4.2 5.1 3.0 4.1 4.4
##  [137] 3.3 5.0 3.7 5.3 5.6 4.6 4.4 3.8 4.8 4.9 4.3 4.7 3.8 3.8 3.9 4.8 3.1
##  [154] 4.3 5.2 4.3 3.7 3.2 4.2 4.7 4.1 3.3 2.9 6.1 4.4 4.3 6.1 5.1 4.0 4.2
##  [171] 3.1 4.9 5.5 4.4 5.5 4.4 4.6 5.1 2.9 3.4 3.0 4.2 4.3 4.6 3.6 4.1 5.6
##  [188] 4.5 4.2 3.8 3.2 4.1 3.6 5.1 5.6 5.7 5.4 4.7 4.3 4.6 4.8 4.6 4.7 3.9
##  [205] 4.7 6.4 5.6 3.8 2.7 4.8 4.5 5.2 4.4 5.1 5.5 4.3 5.2 5.6 5.5 6.2 2.2
##  [222] 4.5 4.3 4.9 5.1 4.9 3.0 4.0 5.4 5.8 4.4 3.8 4.7 5.4 4.4 4.8 5.3 2.5
##  [239] 4.4 5.1 5.1 3.5 3.2 3.2 5.2 5.3 4.2 5.0 5.2 4.7 4.9 4.1 4.0 4.3 5.7
##  [256] 3.4 4.6 4.0 2.9 4.0 2.5 2.6 4.2 3.2 4.1 6.0 4.4 4.9 4.8 3.8 3.1 4.5
##  [273] 3.3 4.7 3.8 6.6 5.5 6.1 4.1 4.7 4.0 4.3 5.1 3.1 5.9 5.2 5.4 3.2 4.2
##  [290] 3.9 6.6 5.0 4.2 3.9 7.3 3.2 4.8 5.2 4.0 4.6 5.0 3.2 5.9 5.8 4.3 5.3
##  [307] 5.4 3.9 4.4 4.2 4.0 3.1 5.7 4.0 2.5 3.8 5.0 5.4 3.0 5.1 4.2 5.8 3.9
##  [324] 4.9 4.2 4.2 4.7 5.0 4.4 5.0 3.8 3.5 4.8 3.0 3.7 5.5 4.4 3.5 4.3 4.3
##  [341] 4.6 5.2 5.0 4.3 4.0 5.3 3.8 5.1 3.4 3.6 3.6 4.7 3.3 4.3 5.4 4.0 4.2
##  [358] 4.9 5.9 5.3 5.1 4.6 4.1 4.0 4.6 4.0 5.0 4.6 5.9 5.2 4.6 4.8 4.1 3.4
##  [375] 4.3 4.6 6.3 4.6 4.0 3.7 4.1 5.1 3.3 4.6 4.4 3.5 4.8 4.6 3.6 4.7 5.5
##  [392] 4.6 5.4 4.7 5.5 4.7 5.2 3.4 5.7 3.7 4.8 5.6 3.9 3.9 4.2 4.0 4.2 5.0
##  [409] 4.3 5.5 4.6 4.2 4.0 3.9 5.0 5.9 3.0 3.3 3.9 4.0 4.9 3.6 3.4 6.3 3.8
##  [426] 3.7 3.6 3.4 4.5 6.4 3.0 4.6 4.0 4.8 3.5 5.4 4.5 4.6 4.8 4.8 3.7 5.1
##  [443] 4.9 3.7 3.5 5.4 4.4 4.3 5.1 3.4 2.8 3.6 5.0 4.6 5.0 5.5 5.1 3.8 4.5
##  [460] 4.4 3.8 4.1 5.2 4.2 4.7 3.7 4.4 3.8 4.6 6.1 5.5 4.6 2.7 5.0 3.0 4.4
##  [477] 4.4 3.9 1.9 5.2 5.8 3.5 5.5 4.6 6.0 5.3 4.7 5.0 4.8 5.9 5.6 4.5 6.4
##  [494] 4.9 5.2 4.8 4.7 3.6 4.1 3.4 5.0 5.8 6.0 5.4 4.2 4.4 5.2 3.8 4.5 5.3
##  [511] 4.5 4.8 4.6 3.1 3.6 3.7 3.2 6.1 4.8 6.2 4.7 5.0 4.3 4.2 6.7 3.7 3.7
##  [528] 3.8 4.7 4.2 4.2 5.0 4.1 4.9 4.3 4.0 5.4 4.4 5.0 5.0 4.1 4.4 5.5 5.2
##  [545] 5.0 3.2 3.9 6.1 3.0 3.0 4.9 3.9 2.6 5.7 4.0 4.0 4.0 4.2 4.8 4.1 4.9
##  [562] 4.2 4.4 3.4 3.2 4.7 4.8 4.1 6.2 5.1 6.0 4.6 5.6 3.2 4.8 4.0 2.0 3.5
##  [579] 5.6 5.3 3.8 6.7 6.0 4.9 2.8 4.3 3.8 3.2 5.5 4.7 5.1 3.2 4.5 3.7 4.8
##  [596] 5.0 3.4 2.9 3.7 3.4 4.4 5.4 4.0 3.4 4.6 4.8 3.9 2.9 3.7 4.5 3.6 4.0
##  [613] 4.1 2.8 5.2 5.2 4.3 4.2 3.8 5.0 3.7 3.8 5.7 6.3 4.2 6.5 5.8 4.9 4.3
##  [630] 4.2 4.0 4.8 5.2 4.9 6.1 4.8 5.4 4.9 4.7 3.0 4.3 4.8 5.1 5.2 4.6 3.3
##  [647] 4.8 4.3 4.6 4.8 4.4 6.1 5.3 4.0 4.5 4.6 4.4 3.8 4.0 4.1 5.6 3.9 3.6
##  [664] 5.7 4.7 4.3 4.2 5.2 4.1 5.1 3.6 4.1 5.3 4.4 6.0 4.5 4.4 5.5 4.6 4.3
##  [681] 5.5 3.9 4.8 5.1 4.5 5.6 4.2 4.7 4.7 3.7 6.9 4.5 3.9 4.3 4.4 4.3 5.5
##  [698] 2.4 5.0 3.8 4.8 3.5 4.5 2.6 4.9 5.0 5.5 3.9 4.0 3.7 5.1 4.5 4.3 4.4
##  [715] 3.3 6.3 4.9 2.8 5.1 6.1 5.1 2.8 4.9 6.1 5.8 4.2 4.8 6.3 6.0 4.8 3.2
##  [732] 4.0 4.9 5.2 5.6 4.6 5.5 6.5 5.9 4.8 6.0 5.2 5.2 4.6 5.0 3.6 3.9 4.5
##  [749] 4.0 4.6 5.5 4.7 3.7 3.2 5.5 5.1 4.0 5.0 6.5 3.4 5.7 4.9 4.8 4.7 5.3
##  [766] 4.4 2.9 3.8 4.2 4.7 5.0 5.0 4.9 3.7 3.5 5.6 3.5 2.8 5.0 5.8 5.5 4.8
##  [783] 4.7 5.3 4.4 5.1 3.2 5.5 3.9 5.0 4.8 4.2 5.1 6.3 5.2 4.2 5.1 3.8 5.3
##  [800] 3.6 3.8 4.1 5.9 4.7 4.1 3.6 2.9 3.2 5.1 5.9 4.1 4.8 5.0 6.6 4.0 5.4
##  [817] 3.8 3.6 4.1 4.4 3.7 3.9 4.5 4.4 3.6 6.3 4.4 4.6 4.9 4.2 2.5 5.0 5.2
##  [834] 5.4 4.6 5.1 4.1 4.2 4.1 3.0 4.9 2.8 4.6 4.5 5.2 4.5 4.9 4.7 4.3 4.6
##  [851] 4.3 4.3 3.7 4.3 4.0 3.1 2.8 5.7 5.1 3.9 3.4 3.2 5.6 5.2 3.4 3.7 4.8
##  [868] 5.4 4.3 5.5 3.2 3.9 4.9 4.9 5.1 3.2 3.7 4.8 5.3 4.9 4.5 3.6 4.7 5.9
##  [885] 7.0 4.6 3.3 3.8 5.1 5.0 3.9 5.5 4.4 3.6 4.2 4.5 4.5 4.2 5.1 4.7 4.8
##  [902] 5.1 4.8 3.8 4.8 3.3 3.0 5.0 5.7 4.4 5.2 5.9 4.4 5.0 3.7 3.8 4.3 4.2
##  [919] 4.9 4.1 6.1 4.6 4.6 4.4 4.1 5.1 6.6 4.6 6.3 5.1 5.8 5.8 5.6 4.7 4.3
##  [936] 4.7 5.8 3.2 4.8 4.3 5.8 4.4 6.4 3.5 5.4 5.4 3.1 5.4 5.0 3.7 3.8 4.1
##  [953] 4.1 5.8 4.1 4.9 4.8 5.6 5.0 4.1 3.5 2.9 4.2 6.5 4.6 3.5 3.3 3.5 4.1
##  [970] 5.4 3.3 5.6 5.7 3.5 5.3 4.1 5.0 4.1 5.1 3.8 3.4 4.2 4.5 4.4 5.0 4.5
##  [987] 4.3 5.8 3.9 5.7 5.4 7.0 5.8 3.7 4.9 4.7 4.6 4.7 3.7 3.3
## 
## $func.thetastar
## NULL
## 
## $jack.boot.val
## NULL
## 
## $jack.boot.se
## NULL
## 
## $call
## bootstrap(x = x, nboot = 1000, theta = theta)
```

```r
quantile(results$thetastar,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.8   6.3
```

Notice that we get exactly what we got last time. This illustrates an important point, which is that the bootstrap functions are often no easier to use than something you could write yourself.

You can also define a function of the bootstrapped statistics (we have been calling this theta) to pull out immediately any summary statistics you are interested in from the bootstrapped thetas.

Here I will write a function that calculates the bias of my estimate of the mean (which is 4.5 [i.e. the mean of the number 0,1,2,3,4,5,6,7,8,9])


```r
bias<-function(x)
  {
  mean(x)-4.5
  }
results<-bootstrap(x=x,nboot=1000,theta=theta,func=bias)
results
```

```
## $thetastar
##    [1] 4.9 3.1 3.8 4.4 6.4 3.8 4.2 4.3 5.0 4.5 5.4 3.5 4.8 4.5 3.8 3.5 4.9
##   [18] 5.7 3.9 4.1 4.4 3.4 4.4 4.9 4.4 3.7 3.0 4.2 3.3 5.0 5.8 4.3 3.4 5.6
##   [35] 4.1 3.9 3.5 6.4 2.3 4.3 1.9 3.7 4.4 2.8 5.0 4.4 4.5 4.0 3.5 5.2 4.7
##   [52] 3.2 3.4 3.2 4.5 5.3 4.5 4.9 4.3 6.2 4.8 5.1 5.4 3.8 5.0 3.3 5.5 6.3
##   [69] 4.0 3.3 3.7 7.2 4.3 2.9 5.1 4.5 5.6 3.9 5.1 3.8 4.2 5.3 4.2 5.3 5.5
##   [86] 4.8 3.4 3.4 3.6 3.4 4.9 3.2 3.5 4.0 4.8 4.4 3.9 4.7 4.6 4.2 4.4 3.3
##  [103] 3.9 4.6 3.7 4.5 3.0 3.5 5.7 3.0 2.7 3.3 6.1 5.5 5.3 4.2 3.8 3.9 4.4
##  [120] 4.9 3.7 5.0 4.9 4.5 4.9 5.0 4.8 2.6 4.1 5.5 4.1 3.8 4.3 6.2 4.2 5.6
##  [137] 4.7 3.8 3.9 5.3 4.4 4.5 6.8 5.8 3.9 4.4 4.1 5.4 4.1 4.9 4.8 4.2 3.8
##  [154] 3.1 5.3 5.2 3.1 5.1 4.0 3.8 3.4 5.2 5.3 4.0 4.2 4.6 3.9 4.5 3.7 4.8
##  [171] 4.5 4.7 5.0 3.0 5.2 4.0 5.9 4.2 5.2 3.4 3.3 4.1 6.1 3.9 4.6 4.4 3.8
##  [188] 5.0 3.9 5.2 4.9 3.6 5.2 5.3 5.8 5.7 4.5 5.1 5.0 4.6 3.4 3.4 3.5 4.5
##  [205] 4.2 4.0 4.7 4.4 5.8 4.1 4.0 5.0 3.6 3.9 5.8 2.1 3.8 4.7 3.6 4.0 4.9
##  [222] 3.8 5.3 4.5 5.4 4.2 3.8 4.5 3.6 2.5 3.6 5.7 4.2 5.1 5.4 4.7 1.8 5.1
##  [239] 3.6 4.9 3.0 5.1 4.2 4.7 5.1 2.5 5.0 5.2 4.0 4.1 3.9 4.9 3.0 5.4 4.7
##  [256] 5.3 2.8 3.9 4.0 6.6 5.1 4.2 5.3 5.8 5.5 3.3 5.2 4.6 4.3 4.9 4.4 4.5
##  [273] 4.8 3.6 5.4 4.3 5.3 4.4 5.1 2.0 5.2 4.0 5.2 4.7 4.5 6.1 3.9 4.1 3.9
##  [290] 4.6 4.5 5.1 5.4 3.8 4.4 4.7 4.4 4.0 3.8 3.6 3.8 3.9 3.4 5.4 3.5 3.7
##  [307] 5.7 3.9 4.5 5.2 4.4 4.6 4.0 4.4 5.9 4.1 3.9 4.6 3.6 3.0 4.5 5.3 4.6
##  [324] 4.4 5.3 5.1 4.2 5.8 2.9 5.1 6.9 5.2 4.7 4.1 3.1 3.5 4.3 4.9 4.5 4.1
##  [341] 4.2 5.6 4.9 6.1 5.0 4.9 6.2 4.5 3.9 4.7 4.1 4.2 6.0 3.4 5.5 4.9 5.4
##  [358] 4.8 4.3 3.8 4.4 3.2 3.6 4.4 4.5 5.0 3.9 5.5 4.2 5.0 2.9 4.3 5.1 3.5
##  [375] 5.3 5.3 4.2 3.1 3.5 4.1 5.3 4.5 4.8 3.1 4.9 3.9 4.0 4.4 4.6 4.6 4.6
##  [392] 5.0 4.2 3.9 5.4 4.6 4.9 4.6 4.9 6.1 4.3 5.2 4.5 3.9 4.5 3.0 3.9 5.2
##  [409] 4.9 4.6 4.7 4.6 4.4 4.9 3.2 4.6 5.1 4.4 3.1 4.2 2.9 3.8 4.5 3.2 4.2
##  [426] 3.3 4.3 5.0 3.9 3.9 5.6 4.2 2.7 3.2 5.7 6.5 5.0 4.4 5.1 4.0 4.4 4.6
##  [443] 5.0 4.7 4.2 4.7 3.9 2.8 4.4 4.3 4.9 5.1 4.5 3.3 5.0 3.9 5.4 4.6 5.5
##  [460] 5.2 3.1 4.6 2.8 4.4 4.6 5.0 4.8 4.2 3.8 2.7 4.8 4.4 4.9 4.7 5.7 5.1
##  [477] 4.4 4.6 3.6 3.2 3.9 3.7 3.5 4.4 4.9 4.2 4.1 5.6 5.1 3.7 4.0 4.2 5.1
##  [494] 3.5 4.0 4.5 4.2 3.5 5.8 5.4 6.1 4.7 5.0 4.8 4.3 4.4 5.5 5.2 4.8 4.7
##  [511] 5.9 3.3 5.8 4.5 3.4 5.4 4.4 5.5 3.3 3.8 5.7 5.1 3.0 4.0 4.4 2.6 4.8
##  [528] 3.6 4.3 3.6 2.7 5.2 3.1 3.8 3.7 3.8 3.4 4.8 4.5 4.8 5.7 3.6 4.7 2.8
##  [545] 4.5 5.5 4.3 4.9 2.4 6.0 5.0 4.9 4.5 3.8 5.1 5.5 3.3 5.8 5.3 3.7 5.1
##  [562] 4.5 3.0 4.4 5.0 4.3 3.8 5.3 6.6 5.7 4.0 4.3 4.5 3.6 4.6 5.4 6.5 4.6
##  [579] 5.9 4.5 5.5 3.9 5.5 5.5 3.5 5.7 4.6 3.8 5.6 5.5 4.8 4.8 4.3 3.4 4.4
##  [596] 3.5 5.9 4.6 6.3 3.5 3.3 3.5 5.3 4.5 5.6 4.6 6.9 4.0 3.3 2.6 5.1 3.1
##  [613] 4.6 4.1 2.9 5.3 2.9 4.1 4.1 4.3 5.2 3.5 3.9 4.3 5.2 5.5 3.2 5.7 4.3
##  [630] 4.5 3.5 4.7 4.8 3.5 5.5 4.7 5.0 3.1 4.6 4.5 4.0 5.7 4.9 5.1 5.0 3.8
##  [647] 5.9 4.0 3.6 3.3 5.2 5.3 3.1 4.5 4.2 4.5 5.5 3.2 4.8 4.3 4.7 3.9 4.2
##  [664] 3.1 4.0 4.9 5.2 5.7 3.4 4.2 5.1 4.6 4.5 4.6 3.3 4.8 5.2 4.5 4.6 4.2
##  [681] 5.2 4.7 5.0 6.1 4.1 5.3 4.2 4.7 4.2 4.3 4.7 4.3 2.9 4.8 4.3 4.3 2.8
##  [698] 3.9 4.2 5.3 4.9 4.1 5.8 4.1 5.0 3.8 4.0 3.7 2.6 3.7 3.4 4.0 5.1 4.7
##  [715] 3.3 4.0 4.4 2.8 3.1 5.5 4.7 5.7 3.0 4.8 3.0 4.1 3.1 5.2 4.6 5.7 5.0
##  [732] 5.4 2.6 4.2 6.1 6.4 4.0 4.9 4.3 5.0 3.9 4.7 5.2 3.4 4.3 4.5 5.3 4.0
##  [749] 2.5 4.4 4.0 4.0 3.8 5.4 3.3 6.1 5.6 3.8 3.7 3.4 4.4 4.9 5.6 4.0 6.0
##  [766] 3.9 4.3 4.7 4.6 3.9 5.8 5.3 4.8 5.2 4.9 3.7 4.3 5.8 3.9 3.8 5.5 4.0
##  [783] 2.7 4.8 4.1 5.5 4.4 6.0 5.2 4.0 5.3 5.1 4.1 3.8 4.8 3.3 3.4 3.9 3.7
##  [800] 1.8 4.8 5.3 5.1 5.9 1.3 4.7 4.7 5.1 6.1 4.3 3.9 5.9 5.7 4.7 3.9 4.5
##  [817] 5.5 5.2 3.6 3.3 7.7 5.8 3.0 5.1 3.0 4.1 3.3 4.0 3.9 4.4 4.1 4.2 3.8
##  [834] 3.3 2.7 5.2 3.8 4.5 3.7 5.1 4.3 4.3 4.2 3.4 4.9 4.5 4.0 6.0 6.0 3.3
##  [851] 4.0 6.4 4.8 5.5 3.1 4.7 4.1 3.6 4.5 4.5 3.1 3.5 3.0 3.9 3.1 5.4 2.6
##  [868] 4.9 5.1 4.1 5.5 4.2 3.3 2.7 3.7 4.1 4.1 4.2 4.7 4.1 5.4 3.5 4.8 5.5
##  [885] 3.3 6.1 5.2 4.4 2.2 4.2 6.3 5.8 5.7 4.7 5.4 5.0 3.4 3.4 3.3 3.3 4.7
##  [902] 4.1 3.3 3.2 4.2 4.7 3.8 5.2 5.6 4.1 4.0 4.3 4.4 4.2 4.2 4.5 4.4 2.2
##  [919] 6.2 4.0 3.9 5.8 4.9 5.0 6.1 3.6 5.3 4.7 4.9 4.4 5.3 4.9 3.1 3.3 4.6
##  [936] 5.4 5.1 4.2 3.4 6.7 3.0 3.9 6.4 4.5 4.8 5.5 3.5 3.8 4.4 3.5 5.6 3.6
##  [953] 5.9 3.9 4.6 4.1 3.6 3.5 4.3 4.3 3.3 5.0 4.7 4.8 4.2 6.4 4.3 4.7 5.6
##  [970] 5.2 5.5 4.2 5.6 5.5 4.3 5.2 5.6 4.1 4.3 5.2 3.7 3.4 4.2 4.2 4.6 3.0
##  [987] 4.0 3.3 3.6 3.8 3.6 4.6 4.4 4.4 4.8 3.8 5.2 4.0 5.1 4.1
## 
## $func.thetastar
## [1] -0.0706
## 
## $jack.boot.val
##  [1]  0.44575758  0.30144928  0.15280899  0.07665706 -0.02155689
##  [6] -0.10143678 -0.25187032 -0.31164773 -0.44460227 -0.52026316
## 
## $jack.boot.se
## [1] 0.9064458
## 
## $call
## bootstrap(x = x, nboot = 1000, theta = theta, func = bias)
```

Compare this to 'bias.boot' (our result from above). Why might it not be the same? Try running the same section of code several times. See how the value of the bias ($func.thetastar) jumps around? We should not be surprised by this because we can look at the jackknife-after-bootstrap estimate of the standard error of the function (in this case, that function is the bias) and we can see that it is not so small that we wouldn't expect some variation in these values.

Remember, everything we have discussed today are estimates. The statistic as applied to your data will change with new data, as will the standard error, the confidence intervals - everything! All of these values have sampling distributions and are subject to change if you repeated the procedure with new data.

Note that we can calculate any function of $\theta^{*}$. A simple example would be the 72nd percentile:


```r
perc72<-function(x)
  {
  quantile(x,probs=c(0.72))
  }
results<-bootstrap(x=x,nboot=1000,theta=theta,func=perc72)
results
```

```
## $thetastar
##    [1] 4.6 2.9 5.4 4.6 4.4 4.6 3.2 4.5 4.0 4.6 5.4 3.6 5.4 4.7 6.1 4.8 5.2
##   [18] 4.7 3.7 5.6 5.6 6.6 5.4 3.7 4.8 5.9 4.1 5.1 4.2 4.1 5.9 4.9 5.0 3.8
##   [35] 4.5 2.8 5.9 6.2 3.8 3.2 2.9 2.8 5.5 3.3 2.8 5.2 3.5 2.8 3.7 3.5 7.0
##   [52] 3.7 6.1 4.8 5.1 3.1 4.8 6.0 4.4 5.5 3.6 4.6 4.1 3.7 5.5 5.6 5.3 4.2
##   [69] 3.9 3.3 4.9 5.1 5.1 5.7 4.2 4.0 3.6 5.2 5.4 4.5 5.2 4.9 5.2 5.5 3.8
##   [86] 6.1 3.7 4.9 5.9 4.3 3.2 3.5 4.3 4.0 5.8 5.1 4.9 5.0 4.5 5.0 4.4 3.5
##  [103] 4.2 4.7 5.8 4.9 3.8 3.4 3.5 2.6 2.7 3.9 3.6 5.7 4.0 3.9 4.1 4.9 5.4
##  [120] 4.8 7.0 5.1 3.4 5.0 4.6 4.2 5.0 6.2 4.5 5.8 3.5 3.5 4.1 4.3 6.5 5.6
##  [137] 5.3 5.0 5.7 4.4 3.3 4.6 3.0 5.8 4.6 4.7 5.1 3.0 3.5 4.4 4.1 6.0 3.6
##  [154] 4.6 4.6 4.5 5.8 5.2 5.7 4.9 4.4 3.8 5.4 5.2 4.4 4.3 4.1 5.0 3.7 5.7
##  [171] 4.8 5.2 3.2 3.9 5.1 4.3 5.3 5.9 2.7 3.1 4.9 4.5 4.3 3.9 5.7 4.9 3.9
##  [188] 6.6 5.1 4.6 6.0 4.1 3.2 5.1 4.6 4.2 4.2 5.1 4.9 2.7 6.6 3.3 4.9 4.2
##  [205] 3.2 4.6 5.7 4.7 2.3 3.7 4.0 3.3 4.1 5.0 7.0 4.3 2.7 4.2 2.6 4.3 3.8
##  [222] 5.4 3.9 3.7 4.0 3.1 3.9 4.8 3.5 4.2 4.1 4.0 4.7 3.9 3.5 6.3 4.7 4.5
##  [239] 4.6 5.8 4.2 5.3 3.9 5.0 5.9 5.7 4.2 3.7 4.4 5.5 3.9 3.9 4.7 4.3 5.8
##  [256] 4.6 5.4 3.4 4.1 4.9 5.6 7.0 3.9 4.5 4.3 3.7 3.8 4.5 6.4 4.3 3.4 2.6
##  [273] 3.4 3.9 3.8 5.2 6.3 4.6 4.3 5.6 4.9 4.7 4.0 4.7 3.9 4.8 4.5 4.2 4.2
##  [290] 3.6 4.6 6.4 4.7 3.3 3.7 5.6 4.7 4.9 5.3 5.7 3.5 4.9 3.9 3.9 3.6 4.6
##  [307] 3.5 4.7 6.0 3.9 3.6 4.0 4.3 5.0 5.7 3.6 5.9 3.0 3.0 2.7 4.8 2.9 5.6
##  [324] 5.0 7.3 3.8 5.7 5.4 3.1 4.5 4.3 5.2 5.4 4.1 3.5 4.6 3.8 3.5 4.8 5.4
##  [341] 4.5 3.7 3.7 4.2 5.2 4.1 3.5 4.9 3.8 4.5 5.2 4.6 4.2 5.4 4.4 6.1 5.1
##  [358] 4.0 3.2 3.7 3.7 3.7 4.5 4.7 4.2 6.5 3.3 4.1 4.2 4.1 3.5 3.8 5.0 4.2
##  [375] 5.2 3.0 3.8 4.5 6.1 6.2 4.2 5.4 5.1 4.3 1.6 4.5 5.3 6.1 5.0 3.9 3.7
##  [392] 4.3 3.9 4.9 3.7 4.0 4.4 4.3 4.8 2.1 3.5 5.2 3.4 3.7 6.1 5.5 5.2 3.5
##  [409] 4.7 3.8 4.9 5.6 2.8 5.1 5.9 5.5 5.5 4.7 4.4 5.8 5.2 4.6 4.9 4.5 5.2
##  [426] 3.7 5.3 4.5 5.3 5.0 4.0 4.7 4.7 5.7 4.4 4.2 6.4 5.7 6.0 3.9 3.5 4.9
##  [443] 5.2 3.7 4.6 4.7 5.3 5.1 3.4 4.5 3.4 5.1 4.5 3.4 3.7 5.5 4.8 4.4 4.6
##  [460] 4.6 4.4 5.6 5.1 4.3 6.5 6.0 5.2 4.5 4.1 4.7 3.5 3.9 5.7 4.0 4.7 3.3
##  [477] 5.2 4.5 3.4 5.2 4.1 5.3 4.5 5.3 4.8 6.5 2.8 4.8 4.3 4.9 4.3 3.3 4.0
##  [494] 3.3 3.3 4.2 6.0 5.9 5.1 4.8 4.6 3.4 5.5 5.4 4.4 6.1 3.1 4.6 6.1 4.1
##  [511] 3.9 7.2 5.5 4.7 4.4 5.0 4.2 3.1 4.4 5.5 6.3 5.1 4.7 4.2 7.2 6.5 3.4
##  [528] 4.4 3.4 4.6 4.4 4.1 4.6 4.2 3.7 4.9 3.5 4.1 5.7 3.6 4.6 3.7 5.4 5.2
##  [545] 3.5 4.7 4.8 3.2 3.3 3.6 2.6 3.7 4.1 4.7 6.5 4.8 4.7 3.3 4.5 5.6 6.1
##  [562] 4.4 4.9 4.8 4.5 4.7 3.9 4.9 4.2 5.3 4.9 4.9 3.1 4.3 4.6 5.3 6.3 4.7
##  [579] 4.5 5.3 3.7 6.8 5.5 4.3 4.6 5.1 5.0 3.8 5.1 5.1 2.4 4.4 4.8 5.7 4.3
##  [596] 4.7 4.1 6.7 5.0 2.5 4.2 3.0 3.5 4.4 4.9 5.1 4.2 5.1 3.4 4.5 4.7 4.1
##  [613] 2.8 5.0 5.5 4.8 5.4 3.7 4.7 5.7 6.7 5.0 4.5 4.2 5.7 4.6 4.6 4.6 4.4
##  [630] 4.1 3.1 4.7 4.6 3.7 6.3 5.8 3.3 4.6 5.1 3.1 4.6 4.6 5.5 6.2 4.0 4.6
##  [647] 5.4 4.2 4.4 4.2 5.1 5.2 4.6 4.4 4.2 5.2 4.2 3.9 5.0 3.5 4.0 5.3 4.7
##  [664] 5.0 5.8 4.3 5.6 4.1 3.7 4.6 5.3 6.0 3.6 4.8 2.8 2.4 4.0 5.5 3.9 3.9
##  [681] 4.3 3.8 3.9 3.5 5.3 4.1 4.5 5.0 5.7 3.6 4.9 4.4 3.9 3.4 5.4 5.7 5.3
##  [698] 5.6 5.4 5.2 4.0 3.3 4.6 6.7 4.9 6.5 4.8 4.8 4.4 5.9 4.2 5.3 4.1 4.2
##  [715] 2.8 4.6 5.7 4.4 4.0 5.6 3.9 4.7 6.1 3.6 4.2 4.3 3.5 3.7 3.0 4.6 3.1
##  [732] 5.7 3.8 2.8 4.1 4.8 4.7 3.6 4.8 4.6 5.6 3.1 5.5 6.1 4.8 4.1 5.2 4.6
##  [749] 4.1 5.6 5.8 3.7 4.2 4.8 5.2 3.8 3.0 5.2 5.7 3.8 3.6 6.3 3.8 3.6 3.7
##  [766] 5.3 4.0 3.4 3.8 4.9 5.0 4.8 5.8 2.9 5.6 3.1 3.0 4.9 5.8 5.1 4.6 3.5
##  [783] 3.1 3.2 4.9 3.4 4.7 4.3 5.8 5.5 4.4 5.5 4.2 6.4 5.5 2.5 5.2 3.7 4.4
##  [800] 4.3 5.7 4.9 4.7 5.0 4.5 5.1 3.7 5.8 3.9 4.4 5.0 3.5 6.0 5.0 4.4 4.3
##  [817] 3.0 3.7 3.8 4.7 3.6 3.6 4.0 4.8 3.7 5.0 4.7 4.4 4.0 4.0 3.6 3.9 6.0
##  [834] 3.8 5.4 3.4 4.7 4.4 6.2 3.6 4.8 4.8 4.5 4.8 5.5 3.5 3.8 3.8 4.0 2.4
##  [851] 4.5 3.6 4.2 6.3 4.7 4.4 4.5 3.0 2.6 3.0 3.6 4.1 4.7 2.8 4.1 4.9 5.4
##  [868] 4.3 4.7 4.5 4.0 4.8 4.2 4.2 4.8 5.6 4.5 3.1 2.5 4.7 4.4 5.5 5.4 4.4
##  [885] 4.3 4.4 3.0 3.6 4.8 5.3 2.9 4.8 5.8 5.8 5.5 4.0 4.5 4.5 3.5 5.0 4.5
##  [902] 4.1 3.9 4.6 2.6 4.1 3.1 5.4 5.2 4.7 4.6 4.5 5.0 4.8 3.8 5.3 5.0 5.5
##  [919] 5.2 5.9 4.8 5.3 3.7 6.0 5.1 3.9 3.4 4.6 4.0 4.5 5.0 3.2 4.8 5.5 2.3
##  [936] 4.3 3.2 3.9 4.3 3.8 3.9 5.4 3.8 4.9 4.1 5.2 4.2 4.8 4.5 5.1 4.0 4.7
##  [953] 5.6 4.0 3.2 4.6 5.5 4.0 2.9 4.7 4.3 3.7 6.0 5.4 3.1 4.0 4.7 5.3 3.9
##  [970] 5.1 3.8 4.0 5.2 4.7 4.9 2.7 4.9 4.8 6.3 4.4 3.8 4.7 3.5 5.1 5.8 6.8
##  [987] 4.7 5.1 5.5 4.5 6.1 4.4 4.3 3.9 5.2 4.1 5.0 6.1 4.0 7.2
## 
## $func.thetastar
## 72% 
## 5.1 
## 
## $jack.boot.val
##  [1] 5.5 5.5 5.4 5.3 5.2 5.1 4.9 4.7 4.7 4.5
## 
## $jack.boot.se
## [1] 1.028786
## 
## $call
## bootstrap(x = x, nboot = 1000, theta = theta, func = perc72)
```

On Tuesday we went over an example in which we bootstrapped the correlation coefficient between LSAT scores and GPA. To do that, we sampled pairs of (LSAT,GPA) data with replacement. Here is a little script that would do something like that using (X,Y) data that are independently drawn from the normal distribution


```r
xdata<-matrix(rnorm(30),ncol=2)
```

Everyone's data is going to be different. With such a small sample size, it would be easy to get a positive or negative correlation by random change, but on average across everyone's datasets, there should be zero correlation because the two columns are drawn independently.


```r
n<-15
theta<-function(x,xdata)
  {
  cor(xdata[x,1],xdata[x,2])
  }
results<-bootstrap(x=1:n,nboot=50,theta=theta,xdata=xdata) 
#NB: xdata is passed to the theta function, not needed for bootstrap function itself
```

Notice the parameters that get passed to the 'bootstrap' function are: (1) the indexes which will be sampled with replacement. This is different that the raw data but the end result is the same because both the indices and the raw data get passed to the function 'theta' (2) the number of bootrapped samples (in this case 50) (3) the function to calculate the statistic (4) the raw data.

Lets look at a histogram of the bootstrapped statistics $\theta^{*}$ and draw a vertical line for the statistic as applied to the original data.


```r
hist(results$thetastar,breaks=30,col="pink")
abline(v=cor(xdata[,1],xdata[,2]),lwd=2)
```

<img src="Week-2-lab_files/figure-html/unnamed-chunk-15-1.png" width="672" />

Parametric bootstrap
---------------------

Let's do one quick example of a parametric bootstrap. We haven't introduced distributions yet (except for the Gaussian, or Normal, distribution, which is the most familiar), so lets spend a few minutes exploring the Gamma distribution, just so we have it to work with for testing out parametric bootstrap. All we need to know is that the Gamma distribution is a continuous, non-negative distribution that takes two parameters, which we call "shape" and "rate". Lets plot a few examples just to see what a Gamma distribution looks like. (Note that the Gamma distribution can be parameterized by "shape" and "rate" OR by "shape" and "scale", where "scale" is just 1/"rate". R will allow you to use either (shape,rate) or (shape,scale) as long as you specify which you are providing.

<img src="Week-2-lab_files/figure-html/unnamed-chunk-16-1.png" width="672" />


Let's generate some fairly sparse data from a Gamma distribution


```r
original.data<-rgamma(10,3,5)
```

and calculate the skew of the data using the R function 'skewness' from the 'moments' package. 


```r
library(moments)
theta<-skewness(original.data)
head(theta)
```

```
## [1] 0.03687528
```

What is skew? Skew describes how assymetric a distribution is. A distribution with a positive skew is a distribution that is "slumped over" to the right, with a right tail that is longer than the left tail. Alternatively, a distribution with negative skew has a longer left tail. Here we are just using it for illustration, as a property of a distribution that you may want to estimate using your data.

Lets use 'fitdistr' to fit a gamma distribution to these data. This function is an extremely handy function that takes in your data, the name of the distribution you are fitting, and some starting values (for the estimation optimizer under the hood), and it will return the parameter values (and their standard errors). We will learn in a couple weeks how R is doing this, but for now we will just use it out of the box. (Because we generated the data, we happen to know that the data are gamma distributed. In general we wouldn't know that, and we will see in a second that our assumption about the shape of the data really does make a difference.)


```r
library(MASS)
fit<-fitdistr(original.data,dgamma,list(shape=1,rate=1))
```

```
## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced
```

```r
# fit<-fitdistr(original.data,"gamma")
# The second version would also work.
fit
```

```
##      shape       rate   
##   1.7456725   2.6363070 
##  (0.7185209) (1.2552254)
```

Now lets sample with replacement from this new distribution and calculate the skewness at each step:


```r
results<-c()
for (i in 1:1000)
  {
  x.star<-rgamma(length(original.data),shape=fit$estimate[1],rate=fit$estimate[2])
  results<-c(results,skewness(x.star))
  }
head(results)
```

```
## [1] 0.6869693 1.2565702 0.5520709 0.3008216 0.1920771 0.3859516
```

```r
hist(results,breaks=30,col="pink",ylim=c(0,1),freq=F)
```

<img src="Week-2-lab_files/figure-html/unnamed-chunk-20-1.png" width="672" />

Now we have the bootstrap distribution for skewness (the $\theta^{*}$ s), we can compare that to the equivalent non-parametric bootstrap:


```r
results2<-bootstrap(x=original.data,nboot=1000,theta=skewness)
results2
```

```
## $thetastar
##    [1] -0.2985073233  0.2464441624 -0.4765161824 -0.5984916971
##    [5]  0.0124193796  0.3890992365  0.0628915809  0.1561811491
##    [9]  0.7014815291  0.0182763042 -0.0423447548 -0.5978992279
##   [13]  0.3639380642 -0.7422112659  0.3877632798 -0.2706662359
##   [17] -0.4804477981  0.1254434119 -0.3425444962 -0.1631575261
##   [21]  0.0448391020 -0.2709227846 -0.4811533248  0.0215045210
##   [25] -0.0383677867  0.3212035127  0.3448641019 -0.0041380293
##   [29]  0.1247535179 -0.2581440976  0.1359001578  0.0200513255
##   [33] -0.2263504617  0.7477637886 -0.0750344043  0.4207512676
##   [37] -0.3824400995  0.6972052431  0.0328378768  0.3601263030
##   [41] -0.2527958345  1.3342697971  0.0697921008 -0.6821359391
##   [45] -0.4813242056  0.3994446004 -0.4723073930 -0.0500306873
##   [49] -0.3905444491 -1.1738173574 -0.2429346863 -0.2595676842
##   [53]  0.1165925567  0.0270020827 -0.0457578616  0.5423250840
##   [57]  0.2736189952  0.0024869052  0.6259427583  0.5445808389
##   [61] -0.9363163511  0.7554142833 -0.4639111675 -0.4413301193
##   [65] -0.1678132594  0.7817642343  0.9697986412 -0.2576040664
##   [69] -0.2411890195  0.0882535509  0.4485862249  0.6970351413
##   [73] -1.0551090136  0.0444994865  0.1582613440  0.0060373564
##   [77]  0.7081994847  0.2149976503  0.0817267982 -0.8504605911
##   [81]  0.3409701104 -0.4976176705 -0.9112572937  0.2197363810
##   [85]  0.1864971428  1.2403202332 -0.0966809412 -0.1915403227
##   [89]  0.2133825071  0.0211429142  0.2086501659 -0.0924251605
##   [93]  0.4706039322  0.6331704789 -0.3039037115 -0.3387798017
##   [97] -0.7399163279  0.3233670005 -0.0364739797  0.7335927132
##  [101]  0.6233837800  0.0194587164 -0.4396094652 -0.1857222260
##  [105]  0.5801303526  0.1792207098  0.4763991408 -0.5093277007
##  [109] -0.6371878358  0.0367991254  0.0958537237  0.0648831280
##  [113]  0.1150629899  0.0854082989  0.0201975984  0.1744625225
##  [117]  0.0249155363  0.7483065621  0.0358279784 -0.0401060940
##  [121]  0.3201678478 -0.7479979259 -0.1982790839 -0.9668976366
##  [125]  0.3470389429  0.0664068358  0.0560247360  0.5662330433
##  [129] -0.2673128206 -0.5581490287  0.0005477915  0.5506466527
##  [133] -0.4467672412  0.6233837800  0.0731774563  0.2174335803
##  [137]  0.0227391236 -0.7734928842 -0.2450586585  0.4556003954
##  [141] -0.9662208976  0.4411799417  0.2719176630 -0.1487190593
##  [145] -0.0411965470 -0.0521620098 -0.0383679986 -0.4506708173
##  [149] -0.2598434525  0.7212122833  0.1083780047  0.1846331437
##  [153]  0.7700626881  0.1941462005 -0.4039740106  0.0732907877
##  [157] -0.1373689941 -0.3886902670 -0.6327017258 -0.5219786946
##  [161] -0.1003502583  0.4047121023 -0.9687851134  0.2550481308
##  [165] -0.4037043200 -0.5986534527  0.6061880094  0.8775529414
##  [169] -0.9961235212  0.6250373896  0.3160085067 -0.2923881551
##  [173]  0.2714448377  0.0764196133  0.9242764501  0.4802625229
##  [177]  0.2766734887  0.1650258730 -0.0598459419  0.7584719955
##  [181]  0.0732473986 -0.2669182423  0.1924986228 -0.2316311275
##  [185] -0.0559624976  0.0417927325 -0.8486656053  1.0548712724
##  [189] -0.6106285630  0.0927647002  0.5877188892 -0.2767885457
##  [193]  0.5916254825  0.1235967859  0.9430602073 -0.5900466617
##  [197] -0.4706310568  0.2373370410  0.1482320046 -0.1996991400
##  [201]  0.4959812376  0.4682163361  0.1995396591  0.1355149006
##  [205]  0.0001454452  0.0552044990 -0.3189455611  0.4026477817
##  [209] -0.2850626771  0.2603191021  0.2235982208  1.2389670268
##  [213]  0.4738787240  0.2317860931 -0.3709680447  0.0245024266
##  [217] -0.2179720487  0.1235616470  0.5986912323 -0.4439977767
##  [221]  0.2618779589 -1.2210949257  0.7945174265  0.2910320449
##  [225]  0.0212088198 -0.7649149648 -0.8638878374  0.9332026899
##  [229] -0.6286216963  0.1012140088 -0.5167669033  0.8611902650
##  [233]  1.1208595689 -0.2157968735 -0.4479334721 -0.0452679438
##  [237]  0.2861427154 -0.1566134179  0.0345516217 -0.3661612354
##  [241]  0.3432799883 -0.9798572672  0.0054314744 -0.0976817232
##  [245] -0.0233361467  0.1075295890  0.2201789866  0.3108417422
##  [249]  0.4293971988 -0.3972452827  0.4375128822  0.2873944861
##  [253] -0.3344088649  0.8126806812  0.2364579678  0.2382817734
##  [257] -0.9491995889 -0.3059263875 -0.1779185101 -1.0079839877
##  [261] -0.2387648314  0.3667906497 -0.1900355935  0.0954566608
##  [265]  0.1824055392  1.6730239828  0.2695262947  0.0992034236
##  [269] -0.2056400287  0.0711427286  0.5594435478 -0.2446106204
##  [273]  0.4299096645 -0.1140613473 -0.2431245333  0.0287646957
##  [277] -0.2188833832 -0.4608969955  0.1284614449  0.0295386731
##  [281]  0.4939943495 -0.0846423952 -0.3629517275  0.6408816231
##  [285]  0.3757530093  0.6044952341 -0.2466479292  0.2695013802
##  [289]  0.2580202323 -0.2476350066  0.5571893236 -0.2157968735
##  [293] -0.2205333177  1.0792276747  0.7778814646  0.1417436770
##  [297] -0.1426088236  0.1772395373  0.1637473638 -0.0154871574
##  [301]  0.1345746135  0.2373727406 -0.2306253897  0.4116097081
##  [305]  0.4496974813  0.6774209737  0.0339545936 -0.4652174325
##  [309] -0.1099034242 -0.4258548240  0.2419977030  0.8503600890
##  [313]  0.3532436599 -0.5548672328 -0.3284988180  0.0532947011
##  [317]  0.1670886700 -0.1157341424 -0.1866648360  0.0390290441
##  [321]  0.3981793474 -0.3105100425  1.3601981280  0.4045507024
##  [325] -0.3687219857  0.0408241641 -0.7248065239 -0.0317615759
##  [329] -0.2942287050  0.1808599013  0.0717910838  0.0343411807
##  [333] -0.2232589527  0.9983713425 -0.3304947371 -0.0066794075
##  [337] -0.3777253071  0.5331550166 -0.2166085819 -0.4517790618
##  [341] -0.0585210814 -0.0138678244  0.0174857274  0.2059129875
##  [345]  0.4391484670  0.4818988921 -0.6161812212 -0.5921299213
##  [349]  0.0981839505  0.0854133076 -0.2156933016  0.2360929555
##  [353]  0.2190170413 -0.7285229776  0.2268965576 -0.6517533347
##  [357] -0.0843197984 -0.4880902041 -0.0326174226 -0.5476807069
##  [361]  1.0190954215  0.8691551716  0.0318041487  0.1669281255
##  [365]  0.3012852247  0.4517750512  0.4153330421 -0.2748700990
##  [369]  0.1294561665 -0.0231303679  0.1385687699 -0.1772757287
##  [373] -0.8437134217 -0.4909904551  0.8634522604  0.2810624841
##  [377]  0.3472619303 -0.2131034352 -0.4436237750  0.0530450316
##  [381] -0.2093289070  0.2153731750  0.0464515565 -0.1511805772
##  [385] -0.7134534077  0.3566168536 -0.6599589814  0.8235842976
##  [389]  0.2018084943 -0.4027791540  0.0713718835 -0.5398537163
##  [393] -0.0322906640  0.5390179349  0.4266570053 -0.0883644632
##  [397] -0.3813699886  0.9981925636  0.6970311097  0.0556293276
##  [401] -0.4437600764 -0.3445731536 -0.0824720300 -1.0105536586
##  [405] -0.3584281858  0.2955871717  0.1814199417  0.4585633931
##  [409] -0.1191078924 -0.0087300145  0.5219187343 -0.3888905110
##  [413] -0.6488333113  0.1074819689  0.7664429585  0.0041777054
##  [417]  1.2057909704  0.2718787153  0.3741137982  0.3684833569
##  [421]  0.3476068105 -0.1143970571  0.5404691672  0.8227454962
##  [425] -0.4007156186  0.2873272308  0.3944775059 -0.3503981967
##  [429] -0.1415589821  0.1232286527  0.3337266026 -0.5894965694
##  [433] -0.4302517566  0.8824520081 -0.3474081534  0.4932076586
##  [437] -0.3826769546  0.9664384366 -0.8107674092 -0.6256217912
##  [441] -0.9343569240  0.2909572434  0.0186208743 -0.0402943150
##  [445]  0.0487643110  0.3981793474  0.4408035908  0.3852277741
##  [449] -0.0677535655  0.2086786733 -0.4827411992 -0.2825289328
##  [453] -0.1544111522 -0.6317210315 -0.3126403652 -0.3395581298
##  [457] -0.0910275435  0.5175468077  0.3117537478  0.1538712966
##  [461]  0.3597037096  0.4845120896 -0.0714476101  0.3712162232
##  [465]  0.7686712886  0.2506229738  0.2047312196  0.3094690348
##  [469]  0.3214668744  0.1029158970 -0.0205089293  0.4962476357
##  [473]  0.6848453299  0.0080018518  0.2349954930  0.2092596929
##  [477]  0.4834299377  0.4128156792 -0.4523163146 -0.5400397169
##  [481]  0.3982655609  0.0133611461 -1.0541674448 -0.0442556604
##  [485]  0.3267755160 -0.0416000898 -0.5019659088  0.0800026832
##  [489]  0.7839245678  0.4239808854  0.6530132181  0.6256233530
##  [493]  0.0884179364 -0.1495110305  0.0198426496 -0.3330931940
##  [497]  0.7167966785  0.7526652947 -1.2067610085  0.0287139675
##  [501]  0.2446967335  0.8527198930  0.9746011808  0.0766431113
##  [505]  0.6824013132 -0.5647308309 -0.1157341424  0.0977590087
##  [509]  0.4843151885 -0.3649232554  0.3867812923  0.4028917467
##  [513] -0.6923950968  0.6337520833  0.6989321257  0.0080082423
##  [517]  0.3952242183  0.5533095168 -0.2739500914  0.3416363185
##  [521]  0.0109772506  0.8624697167 -0.3846491844 -0.5119217848
##  [525]  0.3292562157 -0.2330854687 -0.0867784255  0.4223154890
##  [529] -0.2121690037 -0.3315642586 -0.6992717899  0.2123417080
##  [533] -1.7009024979  0.0702190106 -0.3177971366 -0.1991157063
##  [537] -0.0735531594 -1.1739378307  0.0710382056 -0.0364562154
##  [541] -0.0017167762 -0.4041127159 -0.4510442082 -0.3687757328
##  [545] -0.7397282947  0.6088486103 -0.0297072223  0.3264831879
##  [549]  1.7455572450 -0.1347422490 -0.3442428836 -0.2261020372
##  [553]  0.6820150043 -0.4629470318  0.4532918096 -0.2207289778
##  [557]  0.2985794598 -0.3486188483 -0.2905826626  0.0318041487
##  [561]  0.5022629146  0.1081843862  0.1388689729 -0.0053783322
##  [565] -0.3497055465  0.6710946670  1.1247561183 -0.2035271260
##  [569]  0.1840576597 -0.6243261518 -0.2485868914  0.4893826806
##  [573]  0.0838751982 -0.4102421562 -0.3686316138 -0.3673059587
##  [577]  0.3237085406 -0.0317615759  1.0558849018 -0.1756332139
##  [581]  0.1868283581  0.6642074716 -0.5578765021 -0.3534067383
##  [585] -0.6236582187  0.5659788624 -0.4025222973 -0.2878461722
##  [589] -0.6508963141  0.6624343248 -0.0949131931  0.3105348542
##  [593]  0.5889581663 -0.1015608443 -0.1645176571 -0.3443350275
##  [597]  0.6437469830  0.5841604198 -0.2573979206  0.0819118223
##  [601] -0.4876008256  0.6322977545 -0.6391137507 -0.2707722839
##  [605] -0.2570058717 -0.0802837143  0.0162417675  0.1008185058
##  [609] -0.0952717802  0.5017141412 -0.2425529683 -0.3005941047
##  [613]  0.1970731085 -0.2467233924 -0.1041023449  0.8246584109
##  [617] -0.5178485818  0.3123661636 -0.2016006255 -0.0882239563
##  [621]  0.7776889065  0.2003379427 -0.3974552370  0.0481090694
##  [625]  0.0343772599  0.4696377109  0.0029492217 -1.2188788168
##  [629] -0.5985589835 -0.1715601503 -0.5659362447  0.4304445978
##  [633]  1.1147415654  0.2901751220  0.7304186065 -0.4560523246
##  [637] -0.1731051111  0.5644311661  0.8621113006  0.5542184973
##  [641] -0.0741278621  0.4866189732 -0.7754237734  0.1332857540
##  [645]  0.6073314348  0.1071137319  0.9194266921  0.2210118153
##  [649]  0.6639575424 -0.8715881169 -0.5180452434  0.3632940943
##  [653] -0.2768844090  1.4388894763  0.0487449648  1.6874240810
##  [657] -0.4055133821 -0.3730269864  0.0110003745  0.4328652970
##  [661] -0.4949805694 -0.4811533248 -0.0376417602 -0.1931458784
##  [665] -0.1192639341  0.9113665102 -0.2736961475  0.1889905051
##  [669]  0.3710628205  0.3046612427 -0.3648379305  0.3017182886
##  [673]  0.0310170500 -0.0557758420 -0.2161533214 -0.8148937357
##  [677]  0.4797531779 -0.2974194164 -0.9562465263  0.4439627553
##  [681]  0.0465317663  0.1594443607  0.5659852529  0.6260719897
##  [685] -0.6273586889  0.3233559670 -0.9368327066  0.1887941572
##  [689] -0.2453865954 -0.1872500201  0.6212458169  0.3895705575
##  [693] -0.4270461258 -1.1240449027  0.7451020784  0.4815448051
##  [697]  0.3123328386 -0.3846491844 -0.1481872973  0.0343772599
##  [701] -0.0308438127 -0.1437254469  0.5698860225 -0.0149653377
##  [705] -0.3962111559 -0.2633578458  0.3719462201 -0.3841942563
##  [709] -0.7804367080  0.4889516043 -0.1049604076  0.5892847454
##  [713]  0.4955956747 -0.2560756655 -0.4434262701  0.0405640174
##  [717]  0.5100423018  0.1213998683  0.3222764652 -0.4297912316
##  [721] -0.0678186671  0.6612618297 -0.0671339621 -0.3998089316
##  [725]  0.2564157735 -0.3871034389 -0.1434582891 -0.5870962232
##  [729]  0.8632957741  0.0119605093  0.4339542450  0.0261725698
##  [733] -0.2798235126  0.6725112952 -1.9750362865  0.4935407718
##  [737]  0.6346921540 -0.7166634456  0.7011370833  0.0492666281
##  [741]  0.4538977283  0.2609081377  0.5396961014 -0.4877584049
##  [745] -0.7276876210 -0.1513914718  0.5150995628  0.4682101811
##  [749] -0.1318909293  0.5651254226 -0.2196536918  0.6068304085
##  [753] -0.0099083178  0.4221835585 -0.3508318462 -0.1086526811
##  [757] -0.3387658952  0.3006213708 -0.0852829820 -0.5600081902
##  [761] -0.7648285885  0.7494130444 -0.5815162166  0.4839957739
##  [765] -0.1416077251  0.0212239031  0.4389748847  1.2024517828
##  [769]  0.2421525622 -0.5374443036 -0.4357098890  0.2288246968
##  [773]  0.7642384199  0.0872274186 -0.0543453504 -0.0515808984
##  [777]  0.4391226553  0.2892073951  0.0250171362 -0.5658695919
##  [781]  1.0968380161 -0.0568885863  1.2847495339  0.2019849592
##  [785] -0.5532708396 -0.2581784518 -0.0350326833 -0.2820478641
##  [789] -0.6519042211  0.7865635865  0.4788153376 -0.3230305955
##  [793]  1.1480261742  0.2565792671  0.6488789820 -0.2845223395
##  [797] -0.0636210713 -0.3292894066  0.4566618488 -0.2175912082
##  [801] -0.3037517069  0.9361980678  0.1330976609  0.3366502198
##  [805] -0.2308779109 -0.0281302568 -0.3195169988  0.9688232919
##  [809]  0.3595831644 -0.4798527209  1.2308762207 -0.0633338555
##  [813]  0.0390290441 -0.0735706893  0.6664340170  0.8091893977
##  [817] -0.0842374706  0.2308913565 -0.7009788150 -1.1902851012
##  [821]  0.1556513940  0.5799654571 -0.1227450704 -0.1305695735
##  [825]  0.7807591847  0.6912025199  0.2138990603 -1.2247273391
##  [829]  0.0635412494  0.6036208647  0.3983256525 -0.3047201739
##  [833] -0.0701216034  1.1539079096 -0.1570152223  0.0529390205
##  [837]  0.7174455266  0.6561665863  0.7825217006  0.2578281515
##  [841] -0.0755679332  0.2165827730  0.6576093804  1.1592200163
##  [845]  0.2166447345 -0.4886388561 -0.2149331355 -0.1132871809
##  [849]  0.1224029397 -0.3076102315 -0.5393776224 -0.9433536177
##  [853]  0.0479829329 -0.0425592639  0.2010134838 -0.5048253009
##  [857] -0.4767760720  0.5609834249 -0.2578005647  0.4600317022
##  [861]  0.5423250840  0.7334877621 -0.3869987579  0.3562060511
##  [865]  0.3071069009  0.2758213361  0.6513642343 -0.3629517275
##  [869] -0.3966496109  0.0147477544 -0.6497672690  0.2745968292
##  [873] -0.5210390555 -0.1476293530  1.5281286042  0.3240594110
##  [877] -0.3092779993  0.2795706759 -0.0455851806  0.2341025080
##  [881]  0.1412514104  0.3838230309 -0.7820370885  0.5291655082
##  [885] -0.6406852243  0.5416175287 -0.0591807393  0.8387578238
##  [889] -0.1010053785 -0.9226836470  0.0432827230  0.2923985929
##  [893]  0.5243725822  0.8336705162 -0.0635203994  0.6058001721
##  [897]  0.6088566491 -0.4398720145 -0.3062876921  0.4883815748
##  [901]  0.3651158933 -0.4023834435  0.9700388228 -0.5590766744
##  [905] -0.0458369954 -0.5368672273 -0.5515277516 -0.0526580692
##  [909] -0.0048037430  0.0560406585  0.3377874447 -0.2830215676
##  [913] -0.9326376700 -1.1272211305  0.8333004663  0.3618588359
##  [917] -1.1239288415 -0.1045345655 -0.6791221930 -0.2849347266
##  [921]  0.3308596138  0.5626578091 -0.3055415510  0.3805769686
##  [925] -0.3028076860  0.2852027412 -0.6842002657  1.0642508535
##  [929] -0.1165538725 -0.0158410112  0.0804550439 -0.2476178412
##  [933]  0.2070707255 -0.3340649522  0.1548645211  0.3168951586
##  [937] -0.1468194149  0.5933476200  0.5545024326  0.2631166748
##  [941] -0.4399568636  0.3502473842  0.1478579139 -0.0754817101
##  [945]  0.4784157760 -0.2598434525 -0.5398763986  0.5494921345
##  [949] -0.9823145911  0.0159009579 -0.1284936212  0.7821894596
##  [953]  0.2212580447 -0.1547420403  0.3091298664 -0.1713647711
##  [957] -0.5701832843 -0.2459568667 -0.3292262165 -0.8291510747
##  [961]  0.2966646263 -0.5813252984 -0.8045527356 -0.0032807587
##  [965] -0.9678526440 -1.1419101147 -0.1932530042  0.4536972002
##  [969]  0.9744984981  0.2868108025 -0.3039037115  0.2946260162
##  [973]  0.0929223672  0.5506466527  0.4190963079  0.7746087126
##  [977]  0.1483270162  0.2608715203  0.2065557099 -0.0027695943
##  [981]  1.0518648616  0.7142962515  0.6681795410  0.5763063262
##  [985] -0.5511759163  0.3169824505  1.8473642373  0.2297263010
##  [989] -0.4676714042  0.3100462864  0.0961655884 -0.7054857677
##  [993]  0.0560406585  0.5669289855 -0.2760023054  0.1995396591
##  [997] -0.5265975370  0.3869910273  0.3169441753 -0.6151743943
## 
## $func.thetastar
## NULL
## 
## $jack.boot.val
## NULL
## 
## $jack.boot.se
## NULL
## 
## $call
## bootstrap(x = original.data, nboot = 1000, theta = skewness)
```

```r
hist(results,breaks=30,col="pink",ylim=c(0,1),freq=F)
hist(results2$thetastar,breaks=30,border="purple",add=T,density=20,col="purple",freq=F)
```

<img src="Week-2-lab_files/figure-html/unnamed-chunk-21-1.png" width="672" />

What would have happened if we would have fit a normal distribution instead of a gamma distribution?


```r
fit2<-fitdistr(original.data,dnorm,start=list(mean=1,sd=1))
```

```
## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced
```

```r
fit2
```

```
##      mean         sd    
##   0.6621666   0.4032788 
##  (0.1275280) (0.0901740)
```

```r
results.norm<-c()
for (i in 1:1000)
  {
  x.star<-rnorm(length(original.data),mean=fit2$estimate[1],sd=fit2$estimate[2])
  results.norm<-c(results.norm,skewness(x.star))
  }
head(results.norm)
```

```
## [1]  1.0249052  1.3473262 -0.1693246 -0.1950661  0.5518892  0.3147267
```

```r
hist(results,breaks=30,col="pink",ylim=c(0,1),freq=F)
hist(results.norm,breaks=30,col="lightgreen",freq=F,add=T)
hist(results2$thetastar,breaks=30,border="purple",add=T,density=20,col="purple",freq=F)
```

<img src="Week-2-lab_files/figure-html/unnamed-chunk-22-1.png" width="672" />

All three methods (two parametric and one non-parametric) really do give different distributions for the bootstrapped statistic, so the choice of which method is best depends a lot on the situation, how much data you have, and what you might already know about the underlying distribution.

Jackknifing is just as easy at bootstrapping. Here we will do a trivial example for illustration. We will write a little function for the mean even though you could put the function in directly with 'jackknife(x,mean)'


```r
theta<-function(x)
  {
  mean(x)
  }
x<-seq(0,9,by=1)
results<-jackknife(x=x,theta=theta)
results
```

```
## $jack.se
## [1] 0.9574271
## 
## $jack.bias
## [1] 0
## 
## $jack.values
##  [1] 5.000000 4.888889 4.777778 4.666667 4.555556 4.444444 4.333333
##  [8] 4.222222 4.111111 4.000000
## 
## $call
## jackknife(x = x, theta = theta)
```

Why do we not have to tell the 'jackknife' function how many replicates to do?

Let's compare this with what we would have obtained from bootstrapping


```r
results2<-bootstrap(x,1000,theta)
mean(results2$thetastar)-mean(x)  #this is the bias
```

```
## [1] 0.0371
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8890917
```


Everything we have done to this point used the R package 'bootstrap' - now lets compare that with the R package 'boot'. To avoid any confusion (a.k.a. masking) between the two packages, I recommend detaching the bootstrap package from the workspace with


```r
detach("package:bootstrap")
```


The 'boot' package is now recommended over the 'bootstrap' package, but they give the same answers and to some extent it is personal preference which one prefers to use.

We will still use the mean as the statistic of interest, but we will have to write a new function for it because the syntax of the 'boot' package is slightly different:


```r
library(boot)
theta<-function(x,index)
  {
  mean(x[index])
  }
boot(x,theta,R=999)
```

```
## 
## ORDINARY NONPARAMETRIC BOOTSTRAP
## 
## 
## Call:
## boot(data = x, statistic = theta, R = 999)
## 
## 
## Bootstrap Statistics :
##     original     bias    std. error
## t1*      4.5 0.01161161    0.910188
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 3 4 6 8 9 
## 1 1 2 4 1 1
```

```r
xmeans<-vector(length=1000)
for (i in 1:1000)
  {
  xmeans[i]<-mean(sample(x,replace=T))
  }
mean(x)
```

```
## [1] 4.5
```

```r
bias<-mean(xmeans)-mean(x)
se.boot<-sd(xmeans)
bias
```

```
## [1] 0.0217
```

```r
se.boot
```

```
## [1] 0.9063225
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

