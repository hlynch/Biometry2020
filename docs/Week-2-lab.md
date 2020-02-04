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
## 0 2 3 6 7 8 
## 2 1 1 2 1 3
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
## [1] 0.0193
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
## [1] 2.718956
```

```r
UL.boot
```

```
## [1] 6.319644
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.7   6.3
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
theta<-function(x)
  {
    mean(x)
  }
results<-bootstrap(x=x,nboot=1000,theta=theta)
results
```

```
## $thetastar
##    [1] 5.1 4.8 5.2 5.6 5.2 4.3 4.5 4.0 3.7 6.8 4.0 3.8 5.6 4.0 4.5 4.5 5.0 5.1
##   [19] 4.4 5.4 3.5 3.4 3.3 5.5 3.7 4.8 2.3 3.7 4.9 4.6 3.0 4.4 6.0 5.4 5.1 3.9
##   [37] 4.2 3.5 5.0 3.8 5.1 5.2 4.9 4.6 3.3 5.8 3.0 5.1 4.8 3.2 4.2 5.0 5.0 4.8
##   [55] 4.7 5.6 4.9 3.4 5.2 6.7 4.6 3.8 5.7 3.6 4.2 3.5 3.4 7.2 4.3 5.4 4.5 4.5
##   [73] 4.9 5.8 5.6 4.8 3.8 4.8 5.1 3.6 4.7 4.5 3.1 4.5 3.0 2.9 5.1 3.5 5.4 5.0
##   [91] 4.7 5.2 4.9 3.8 4.2 4.9 4.2 6.7 4.0 5.5 4.2 4.5 4.5 3.5 4.0 5.3 3.7 4.4
##  [109] 5.4 4.5 5.6 4.0 4.6 3.8 4.1 3.8 3.1 3.7 5.6 3.4 4.5 4.8 5.4 4.8 5.2 4.6
##  [127] 4.7 4.5 6.1 4.9 3.5 3.7 2.6 5.8 3.5 4.7 6.6 4.2 3.0 5.0 5.1 4.8 4.5 3.7
##  [145] 5.7 4.1 5.0 4.1 5.3 5.9 5.8 4.8 4.1 2.9 4.1 3.7 4.4 3.7 4.7 4.6 3.0 2.2
##  [163] 3.9 4.4 4.4 6.0 6.1 5.1 4.6 3.9 5.1 4.8 4.5 4.7 5.4 4.1 4.8 4.4 4.8 3.2
##  [181] 5.5 4.2 4.7 4.6 4.6 2.8 4.3 5.0 3.8 5.0 5.6 4.5 3.9 5.5 3.5 4.1 4.2 5.1
##  [199] 5.4 2.9 3.5 5.2 5.0 5.1 4.5 4.2 4.6 5.4 3.0 4.7 4.6 2.4 4.4 5.6 6.3 5.0
##  [217] 4.8 3.9 5.0 6.5 4.3 6.3 4.3 4.2 5.2 3.2 4.9 4.2 4.6 3.9 4.8 2.7 4.0 4.1
##  [235] 5.7 4.8 3.1 5.1 4.7 6.2 3.9 5.5 3.3 4.3 3.1 5.0 4.3 5.2 4.2 5.2 4.7 5.1
##  [253] 5.3 6.9 4.2 4.9 4.3 3.7 5.3 4.2 3.3 4.4 4.0 3.6 4.3 4.9 4.7 3.7 3.1 4.9
##  [271] 4.8 2.7 5.0 4.7 3.8 2.7 4.3 6.5 4.9 4.9 5.4 3.1 3.8 5.0 4.0 5.7 3.8 5.4
##  [289] 3.7 4.5 4.9 1.7 5.2 5.3 4.8 3.4 5.1 4.8 4.8 3.5 5.1 5.3 5.8 5.1 4.2 4.2
##  [307] 4.6 5.4 5.6 4.4 5.1 5.2 6.6 5.5 4.6 6.4 5.7 5.3 4.0 4.1 4.1 4.5 3.9 5.0
##  [325] 5.2 5.0 4.7 4.5 6.0 3.4 4.9 5.8 3.8 3.9 5.3 5.3 4.6 4.5 4.1 5.0 4.7 4.7
##  [343] 5.0 3.8 4.9 4.5 3.1 4.4 3.3 5.1 3.8 3.0 4.4 4.5 5.1 5.1 3.6 4.4 4.6 3.9
##  [361] 5.3 3.3 4.4 5.4 4.1 4.1 5.0 4.0 4.8 5.5 4.2 3.8 4.6 3.8 5.1 5.4 3.8 4.0
##  [379] 5.5 3.5 5.2 5.5 4.5 3.6 5.1 4.7 4.5 4.1 4.1 3.0 5.8 4.7 5.2 4.9 6.4 4.7
##  [397] 3.7 4.2 3.1 4.5 3.1 2.3 5.1 5.9 5.5 4.3 5.5 5.1 5.9 4.7 4.6 4.3 4.5 3.3
##  [415] 4.7 4.9 3.2 5.9 5.0 4.2 4.3 3.3 4.4 6.2 3.1 4.1 3.7 4.9 6.0 5.8 5.7 4.8
##  [433] 5.2 4.0 5.6 3.5 4.7 4.8 5.5 4.2 5.4 4.7 4.4 4.7 4.6 4.9 4.9 5.4 3.9 5.7
##  [451] 5.6 6.5 4.2 4.5 3.8 3.9 4.8 6.6 4.0 3.0 4.2 4.7 5.6 4.6 4.7 3.6 4.4 3.3
##  [469] 5.4 5.1 5.0 5.7 6.7 4.7 7.0 2.6 3.8 2.6 4.5 4.9 4.4 3.0 5.2 4.6 4.0 4.5
##  [487] 4.3 3.6 6.5 5.6 5.3 5.9 6.1 5.4 4.3 5.0 5.9 3.9 2.4 3.9 4.9 5.4 4.6 5.1
##  [505] 4.6 5.1 4.0 4.8 5.3 5.2 5.3 4.8 4.3 3.9 4.1 4.6 3.2 4.7 4.1 4.0 5.8 3.9
##  [523] 4.0 3.7 4.1 3.0 6.1 2.6 5.3 5.5 6.0 6.4 4.1 5.2 4.4 2.9 3.8 5.3 4.2 5.0
##  [541] 5.4 5.4 6.4 5.6 5.7 5.1 4.8 4.5 5.1 5.1 4.3 4.0 5.1 2.8 5.6 5.4 5.8 4.6
##  [559] 5.2 4.4 4.6 4.3 4.8 4.3 6.9 5.3 3.9 4.7 3.0 2.9 6.0 5.4 4.4 5.6 3.4 3.6
##  [577] 4.5 5.0 4.4 4.4 3.6 4.2 2.7 5.0 3.9 4.3 3.5 4.5 4.5 5.2 3.5 6.3 5.0 4.6
##  [595] 5.0 5.4 3.7 4.5 4.4 4.9 3.6 5.8 5.2 5.3 4.9 3.3 4.0 3.4 3.1 5.4 5.8 5.1
##  [613] 6.4 5.6 4.8 5.1 4.5 4.0 3.8 6.3 4.3 4.3 3.7 3.7 4.5 5.7 3.7 4.2 4.6 4.2
##  [631] 5.2 5.9 4.2 5.3 4.5 3.2 4.6 4.5 3.8 5.3 5.3 3.7 4.0 2.7 5.6 4.2 4.6 4.3
##  [649] 4.6 4.3 4.7 5.3 4.2 3.7 4.8 4.6 3.8 5.6 5.1 4.4 4.4 3.3 4.2 4.6 4.3 5.1
##  [667] 4.9 5.9 5.7 4.6 3.0 4.0 4.2 4.8 5.1 4.6 4.5 5.7 4.1 4.8 4.2 5.5 5.4 4.4
##  [685] 5.0 5.4 4.1 4.1 4.8 3.2 5.7 5.8 5.0 5.3 5.0 3.4 4.5 4.5 4.7 2.8 3.6 5.2
##  [703] 4.9 4.6 4.0 3.6 4.7 5.5 3.5 3.8 5.6 5.3 4.1 4.9 5.9 4.9 3.5 5.3 5.8 4.2
##  [721] 3.7 4.8 4.5 5.5 4.5 4.5 3.9 4.8 4.7 3.5 4.1 3.9 6.6 5.2 5.1 4.1 4.6 3.7
##  [739] 5.9 5.4 5.1 3.6 4.4 4.0 5.3 4.8 5.2 5.2 5.8 5.0 4.8 4.6 2.6 4.7 4.5 3.7
##  [757] 5.6 2.8 5.0 4.6 5.3 6.1 4.1 5.9 4.3 2.8 4.9 3.0 3.3 4.6 4.7 5.0 6.4 6.3
##  [775] 2.4 6.1 4.3 4.9 3.6 6.2 3.5 5.9 4.6 1.8 3.9 5.0 4.8 6.4 2.9 4.4 3.5 3.3
##  [793] 3.7 3.8 4.9 4.0 3.3 4.3 3.9 4.5 4.1 3.8 5.4 3.4 4.7 3.8 6.0 4.4 4.8 4.9
##  [811] 4.5 4.1 5.1 5.0 3.9 4.4 2.9 3.1 3.2 5.1 6.0 3.9 3.1 3.9 4.8 4.7 3.9 5.8
##  [829] 4.2 5.0 3.5 5.5 5.3 4.3 5.6 5.9 3.4 4.0 3.0 5.1 6.3 4.9 4.1 4.9 5.7 5.3
##  [847] 4.7 6.0 5.5 3.8 5.6 4.9 5.0 2.3 5.7 4.8 5.5 4.1 3.5 4.6 5.8 4.3 6.2 2.7
##  [865] 5.1 3.9 4.2 2.8 4.0 4.6 3.6 3.8 4.0 5.4 4.1 3.8 4.6 4.1 4.9 3.3 3.1 6.4
##  [883] 4.0 4.3 3.9 6.0 4.5 3.5 4.9 5.1 5.4 5.2 5.2 3.8 5.0 4.4 3.7 4.1 5.1 4.2
##  [901] 5.0 5.4 4.6 5.4 4.5 4.9 5.0 4.5 5.3 2.9 4.9 4.7 4.6 4.7 3.6 4.3 3.4 4.2
##  [919] 5.2 4.2 5.0 5.1 4.0 5.8 4.7 4.9 4.0 4.5 4.9 3.1 4.8 4.6 3.2 3.4 3.8 2.8
##  [937] 3.4 4.1 4.5 4.9 4.0 3.0 4.5 4.6 4.0 5.7 4.8 3.0 2.8 5.2 3.1 7.1 4.7 4.6
##  [955] 2.9 4.2 5.0 5.5 4.3 3.3 5.0 2.9 2.9 2.2 5.0 4.6 6.0 3.9 4.8 4.0 3.2 5.3
##  [973] 5.6 4.5 3.1 4.1 4.7 3.5 6.7 3.8 4.8 3.7 5.5 4.4 5.3 4.0 4.9 5.3 4.9 5.3
##  [991] 4.0 5.7 4.5 6.2 3.8 3.3 3.7 5.1 3.3 5.9
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
##   2.8   6.4
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
##    [1] 5.0 3.0 3.8 4.7 2.6 4.2 5.0 4.0 5.2 5.3 3.3 3.9 4.4 5.0 4.2 3.3 4.3 3.7
##   [19] 3.2 4.6 3.7 6.0 4.5 4.6 4.2 4.0 4.7 4.2 5.1 4.6 4.1 5.0 4.6 2.3 5.4 6.4
##   [37] 4.2 5.3 2.5 4.4 2.7 2.5 4.0 4.1 6.9 4.2 5.2 6.4 4.7 5.0 4.9 4.1 2.5 5.5
##   [55] 4.5 3.8 3.2 5.6 3.8 6.1 5.1 4.1 5.2 6.6 5.0 5.5 4.5 5.5 4.9 2.5 4.8 5.2
##   [73] 4.8 4.4 3.3 4.1 5.1 2.8 4.0 5.5 4.7 4.0 5.3 4.4 5.0 2.7 4.7 5.7 4.4 6.2
##   [91] 1.6 4.6 5.2 3.7 4.2 4.2 6.2 4.8 5.5 5.0 3.7 4.3 2.4 3.8 5.4 2.7 5.9 5.0
##  [109] 4.1 4.1 4.2 4.6 4.5 5.9 5.1 5.9 4.5 4.9 4.8 4.8 2.7 4.7 3.8 4.2 4.8 4.7
##  [127] 5.6 4.3 5.1 4.3 3.3 3.9 2.9 4.6 5.6 5.6 3.8 3.4 3.5 3.1 4.5 4.8 5.6 4.1
##  [145] 4.3 5.2 4.3 4.9 5.3 3.7 4.1 4.7 3.6 5.3 3.5 4.1 4.0 3.9 4.9 5.4 6.1 2.4
##  [163] 4.3 3.8 4.8 4.0 5.6 4.2 3.3 5.2 4.3 3.8 3.5 4.6 3.9 4.0 4.4 4.6 3.1 5.4
##  [181] 4.6 3.6 4.4 6.2 4.7 4.9 4.2 6.1 4.7 5.5 3.2 3.8 3.1 5.0 4.9 5.8 6.6 3.2
##  [199] 3.8 5.4 5.5 4.0 4.9 5.5 4.9 5.2 3.2 4.8 4.8 4.7 4.0 4.3 4.1 4.9 4.5 5.2
##  [217] 5.0 3.4 3.9 5.2 4.8 4.6 3.3 5.7 4.5 4.2 3.0 4.9 4.7 5.2 2.7 4.7 4.1 3.4
##  [235] 4.0 3.4 6.3 5.9 3.4 6.1 4.0 5.4 5.9 5.5 2.4 3.9 4.9 5.9 4.3 5.1 5.1 6.1
##  [253] 5.1 4.9 4.2 4.6 4.8 5.6 5.2 4.3 3.5 5.4 4.3 4.9 4.0 4.5 5.2 3.7 4.4 4.8
##  [271] 2.6 4.5 6.1 3.9 5.4 3.2 5.4 6.3 6.1 4.4 4.6 5.4 5.1 6.5 4.3 4.1 4.8 5.4
##  [289] 5.5 2.8 5.0 3.5 5.2 5.4 4.4 5.7 5.6 3.9 5.5 4.5 3.2 3.3 3.0 4.1 5.0 4.0
##  [307] 5.8 5.3 4.2 4.1 3.6 5.1 4.4 2.8 4.6 3.7 4.0 5.1 5.6 4.3 4.8 4.4 5.3 3.0
##  [325] 4.5 4.1 3.9 5.1 5.3 4.8 3.6 5.3 5.2 4.8 5.1 4.7 4.7 3.6 5.1 4.2 4.1 5.0
##  [343] 4.8 3.1 3.0 5.0 3.7 4.3 3.4 3.7 5.3 3.2 4.7 3.9 4.5 4.4 1.8 1.5 4.3 4.8
##  [361] 5.3 4.3 5.3 3.3 3.7 4.4 3.7 6.1 4.3 4.7 4.9 4.7 4.0 4.4 4.2 5.2 5.5 5.2
##  [379] 4.7 4.8 3.4 3.6 4.1 3.5 4.2 4.7 4.5 4.2 4.5 4.1 2.8 4.0 4.8 6.0 4.7 4.6
##  [397] 4.7 6.0 4.5 5.1 3.9 3.7 3.4 5.1 4.2 4.8 4.2 4.5 6.0 3.2 6.3 4.7 5.3 3.4
##  [415] 4.1 3.8 3.4 3.5 4.3 4.0 4.8 4.6 6.1 3.4 4.0 5.2 6.0 6.1 5.1 5.2 3.0 4.0
##  [433] 4.2 6.1 5.4 3.9 5.2 5.7 5.3 3.9 4.3 4.4 5.9 4.4 2.9 4.3 3.4 4.3 4.2 5.1
##  [451] 5.3 3.3 4.5 3.4 5.5 5.7 4.0 4.3 2.7 4.2 4.9 3.4 4.4 3.0 4.7 4.0 4.6 4.7
##  [469] 4.9 4.3 3.8 4.0 4.0 3.5 5.4 3.3 4.6 5.2 5.6 4.6 4.1 5.4 6.1 3.9 4.7 3.8
##  [487] 3.6 3.9 5.6 4.6 4.6 4.1 5.6 4.5 5.3 3.2 5.5 4.9 4.5 5.2 3.8 5.4 6.0 2.8
##  [505] 4.3 5.0 5.3 5.3 5.2 4.8 4.4 4.3 4.1 4.8 5.2 4.7 6.4 4.1 4.3 5.5 4.8 4.6
##  [523] 5.2 3.9 4.5 5.4 5.6 3.9 4.0 4.3 5.0 3.3 4.7 4.2 5.5 3.8 4.4 4.8 4.5 4.8
##  [541] 4.0 3.8 5.6 4.5 5.4 3.3 2.6 6.6 5.9 4.4 5.9 4.8 3.9 2.2 4.4 3.8 3.6 4.1
##  [559] 3.3 5.3 5.3 4.4 2.9 6.2 5.2 5.1 4.9 4.7 5.7 5.4 4.3 6.8 3.3 3.3 5.4 4.8
##  [577] 4.3 5.1 4.3 4.7 4.3 3.7 4.7 4.4 3.2 5.0 6.0 4.1 4.5 5.5 4.8 5.3 3.4 4.1
##  [595] 5.1 4.4 3.9 5.3 3.9 4.6 5.2 4.2 3.5 3.9 5.9 4.9 5.7 5.0 3.3 4.1 4.6 5.4
##  [613] 4.3 4.0 5.4 3.8 5.0 4.1 4.0 5.9 5.0 4.2 4.2 3.6 3.4 3.7 4.3 5.3 3.8 2.5
##  [631] 6.1 4.4 4.5 5.2 4.3 4.7 5.3 4.4 3.4 3.0 4.7 4.9 4.7 3.8 4.7 3.4 4.3 4.8
##  [649] 5.1 3.2 4.3 4.7 3.5 3.1 4.4 5.5 4.9 4.1 4.1 4.7 4.2 3.6 5.2 4.3 5.8 4.8
##  [667] 4.5 3.9 4.7 4.0 4.9 5.0 6.4 4.3 4.1 3.0 5.4 5.0 4.0 5.4 5.6 3.8 4.2 3.9
##  [685] 6.5 2.2 4.4 3.0 4.0 6.5 4.2 4.9 5.8 3.4 3.2 3.8 2.3 3.7 5.2 4.7 5.7 4.4
##  [703] 4.8 4.3 3.0 4.8 4.3 4.7 4.8 4.8 4.7 4.7 5.8 4.8 4.4 4.0 5.0 4.7 6.0 4.5
##  [721] 4.1 4.2 5.6 3.3 5.4 4.3 5.1 4.1 3.9 4.9 4.8 6.1 2.0 5.3 3.4 3.3 4.2 3.4
##  [739] 5.5 5.1 3.0 3.8 3.6 5.3 5.5 5.4 3.6 4.5 3.8 3.9 5.2 4.1 4.2 4.3 6.2 3.4
##  [757] 3.7 4.9 3.3 3.4 3.8 4.3 5.2 4.7 4.7 3.8 3.1 3.1 6.1 3.3 4.8 2.6 4.7 4.7
##  [775] 3.3 4.0 4.4 3.5 5.5 3.2 4.8 3.4 3.8 4.9 3.7 4.7 4.7 4.4 5.2 3.0 5.3 3.1
##  [793] 4.5 4.0 4.7 5.0 4.8 4.6 4.8 6.1 4.3 3.8 5.8 5.1 4.3 4.1 5.5 4.2 4.9 5.3
##  [811] 3.6 5.1 6.0 6.0 4.3 5.5 4.3 4.0 5.9 5.2 5.5 4.6 5.2 5.3 4.2 3.9 4.5 5.0
##  [829] 3.7 3.7 4.6 4.3 6.9 3.9 4.8 4.0 4.2 4.5 4.9 5.3 4.9 3.6 4.1 5.3 5.0 4.2
##  [847] 5.8 3.6 4.7 4.5 5.7 4.8 5.2 5.4 2.9 4.8 4.9 5.2 3.3 4.7 4.3 4.6 6.2 4.0
##  [865] 5.1 3.4 3.7 2.7 3.8 4.7 3.2 3.7 5.6 4.5 5.1 5.5 3.5 5.2 4.5 4.5 3.6 4.3
##  [883] 5.2 4.4 5.5 3.9 5.2 4.8 4.0 4.2 5.3 4.8 4.8 4.6 3.1 4.5 4.7 2.6 3.5 2.9
##  [901] 4.3 3.5 5.1 5.8 3.4 3.6 3.2 4.0 4.9 5.2 4.6 4.6 4.4 4.7 3.2 5.0 3.9 5.0
##  [919] 4.8 4.8 2.8 5.7 5.3 4.9 5.8 3.3 4.8 4.9 2.9 4.6 3.8 5.7 5.0 3.4 3.7 5.8
##  [937] 3.8 3.7 3.7 4.6 3.5 5.8 3.7 3.9 3.5 3.5 3.8 4.8 5.3 3.4 4.5 2.4 5.4 5.6
##  [955] 3.3 4.1 5.5 5.2 3.7 3.8 6.3 4.3 4.0 4.2 4.0 3.7 4.5 6.7 4.0 3.9 4.1 5.2
##  [973] 4.2 4.6 4.8 5.0 5.9 3.6 4.9 6.2 4.1 4.0 4.2 5.4 4.6 3.0 5.0 4.0 5.7 3.8
##  [991] 4.9 4.8 4.7 3.8 4.6 4.6 4.7 3.9 5.8 4.7
## 
## $func.thetastar
## [1] -0.012
## 
## $jack.boot.val
##  [1]  0.446685083  0.351704545  0.231104651  0.171978022  0.002236422
##  [6] -0.117824773 -0.132773109 -0.197105263 -0.419075145 -0.492415730
## 
## $jack.boot.se
## [1] 0.8961196
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
##    [1] 5.5 4.4 4.9 4.7 3.0 5.2 3.9 3.6 3.8 4.3 4.4 2.5 5.3 5.3 3.2 3.9 4.1 5.0
##   [19] 4.4 5.3 4.8 4.1 6.0 4.7 4.8 4.3 4.5 5.8 4.6 5.4 4.8 3.8 5.6 4.5 5.6 5.6
##   [37] 5.3 4.2 5.9 4.9 3.2 3.0 4.4 5.6 4.2 3.3 4.2 3.7 3.8 4.1 5.1 4.9 4.4 5.4
##   [55] 3.2 3.8 4.8 5.8 3.0 4.5 3.3 5.3 4.9 5.1 4.8 4.7 4.0 4.6 4.8 3.9 5.0 5.0
##   [73] 4.2 4.1 4.7 4.3 5.3 4.2 5.7 3.7 4.6 3.5 5.1 4.4 5.0 4.6 5.1 3.5 2.7 6.2
##   [91] 3.7 3.0 5.5 6.4 4.2 4.8 4.8 5.0 3.8 4.5 4.6 4.4 5.2 4.0 5.3 4.1 3.7 4.0
##  [109] 4.2 4.1 5.7 4.4 3.3 3.6 4.0 2.8 3.1 4.9 4.5 5.8 5.0 1.9 4.6 4.3 4.3 4.7
##  [127] 3.3 3.9 3.5 6.4 3.1 4.2 3.0 4.7 3.3 3.9 5.5 4.1 4.9 6.6 4.6 6.7 5.8 3.2
##  [145] 4.2 4.3 5.7 5.5 4.1 4.9 4.6 3.9 3.7 4.8 4.0 4.4 3.6 4.6 3.5 3.9 3.6 3.6
##  [163] 3.7 5.2 3.8 4.8 5.6 4.2 4.8 3.4 4.1 5.6 3.9 3.8 4.7 5.2 6.1 5.2 4.5 4.5
##  [181] 5.5 4.4 4.4 2.5 3.8 5.2 4.0 4.2 2.7 4.1 5.4 4.2 4.4 3.1 4.2 3.7 4.9 5.0
##  [199] 3.6 5.5 5.1 5.2 3.2 4.1 4.2 5.9 4.9 4.6 4.2 4.4 3.3 4.1 4.4 4.1 3.7 4.2
##  [217] 4.9 3.1 3.6 3.2 3.0 4.7 4.9 2.3 4.5 4.3 5.3 5.3 5.3 3.9 4.8 5.3 3.6 2.2
##  [235] 5.2 3.5 3.6 5.0 2.4 5.1 3.1 3.6 4.2 4.3 4.0 6.4 1.3 4.3 3.7 5.7 5.0 3.9
##  [253] 4.9 3.5 5.5 3.9 3.5 4.4 4.2 4.0 4.8 3.7 6.6 4.4 6.3 5.3 3.2 4.1 3.9 4.5
##  [271] 4.2 6.0 5.5 4.5 5.5 5.4 4.7 4.5 4.2 5.5 4.4 5.2 5.1 5.9 5.5 4.8 3.5 3.3
##  [289] 5.1 3.9 4.7 4.8 3.7 5.5 4.2 6.0 5.2 5.0 5.4 3.7 5.1 4.2 5.5 4.1 5.0 4.9
##  [307] 3.9 5.3 6.6 4.7 5.8 5.5 3.9 5.2 5.7 5.8 4.0 5.2 4.1 4.5 6.0 3.0 5.2 4.8
##  [325] 5.6 3.0 3.5 4.5 5.4 3.0 3.4 6.3 3.3 5.0 4.8 5.7 7.1 5.2 3.9 4.0 6.1 4.1
##  [343] 4.1 3.8 4.3 4.6 4.1 6.0 4.0 5.0 4.9 4.6 4.0 4.7 4.3 3.1 4.8 4.3 4.6 3.9
##  [361] 5.4 6.0 4.8 4.5 4.5 5.3 3.7 3.6 5.2 5.9 5.1 2.4 4.4 3.3 4.9 5.2 4.1 5.3
##  [379] 4.2 4.7 5.0 4.9 3.5 5.7 3.7 2.7 4.2 4.8 2.7 4.1 4.3 3.8 3.0 2.8 5.3 3.8
##  [397] 4.7 3.9 4.5 5.1 3.4 5.3 2.9 5.7 4.4 5.1 4.1 4.2 3.5 5.2 3.2 3.6 5.1 4.4
##  [415] 4.3 5.2 6.7 5.4 4.5 4.8 4.3 6.3 3.5 4.0 4.3 3.0 4.4 4.5 5.6 3.8 4.1 3.7
##  [433] 2.9 4.8 4.9 3.9 4.8 3.4 6.0 5.2 4.7 5.3 5.4 6.9 4.8 3.5 3.6 5.1 4.8 4.9
##  [451] 4.4 5.3 6.3 3.9 4.2 3.8 3.3 3.2 3.3 5.5 5.6 5.3 4.9 4.5 4.9 5.0 3.9 3.4
##  [469] 5.3 4.6 4.5 5.9 2.5 4.1 4.4 4.5 4.3 4.8 4.6 4.9 5.6 5.6 2.6 6.5 3.2 4.4
##  [487] 4.3 3.7 5.4 5.0 4.8 2.5 4.6 3.1 3.9 4.0 3.8 4.0 7.0 3.5 4.2 3.9 4.4 4.0
##  [505] 4.7 4.0 4.6 3.9 3.7 3.7 4.9 5.2 3.8 4.6 5.2 3.9 4.7 4.5 4.2 4.5 3.1 5.4
##  [523] 3.6 3.8 6.1 4.4 5.0 4.2 3.9 4.8 4.2 5.2 4.2 3.1 4.6 3.7 2.5 4.9 4.5 4.5
##  [541] 4.2 4.6 2.8 3.3 2.1 5.5 2.3 5.1 5.5 5.5 4.6 5.2 3.7 5.1 4.1 6.1 4.0 4.3
##  [559] 4.8 3.9 4.8 5.1 3.4 6.7 3.7 3.4 5.1 3.8 3.7 5.8 5.1 4.2 4.5 5.1 4.8 3.9
##  [577] 4.8 4.9 5.1 3.7 5.2 3.1 5.9 4.7 4.0 5.3 4.7 3.9 4.4 5.5 3.6 4.8 3.6 4.0
##  [595] 4.6 4.7 3.3 4.6 4.2 4.6 5.4 5.3 4.4 5.8 5.9 3.7 3.7 4.4 4.3 6.2 4.0 5.4
##  [613] 4.6 4.8 3.0 5.4 4.9 4.0 6.0 4.1 3.9 3.4 5.2 4.2 3.4 6.0 4.9 6.1 4.6 3.6
##  [631] 3.9 3.5 4.3 5.7 4.5 4.9 4.4 3.8 5.5 4.6 3.6 3.5 5.1 2.7 3.4 5.5 2.7 5.0
##  [649] 5.0 4.1 4.8 5.0 4.7 4.7 5.2 5.3 5.0 2.8 5.2 4.9 4.3 5.4 5.7 5.0 3.7 3.7
##  [667] 5.1 4.3 2.7 3.8 4.6 5.2 3.9 5.0 3.9 3.9 4.0 4.0 4.5 6.7 3.8 4.0 3.1 5.8
##  [685] 4.8 5.0 3.9 5.1 5.5 4.6 6.2 4.8 5.6 4.9 4.4 4.8 5.7 5.2 3.2 4.1 3.6 5.2
##  [703] 4.9 5.0 5.1 4.8 5.5 5.0 4.6 5.6 1.7 5.4 5.5 4.1 4.2 3.8 4.3 4.9 3.2 4.5
##  [721] 5.4 4.8 4.3 2.8 6.1 4.2 5.8 3.4 4.7 5.5 5.3 3.1 5.3 5.7 5.7 3.7 3.7 5.8
##  [739] 4.8 4.1 5.2 4.9 3.8 4.8 3.7 5.5 5.4 4.5 5.1 4.9 2.9 3.5 3.6 4.9 5.0 3.9
##  [757] 5.3 3.2 4.6 6.3 4.3 5.5 4.4 4.6 3.8 4.5 3.2 4.0 3.8 4.5 6.3 4.7 5.3 5.8
##  [775] 4.8 5.4 4.8 4.8 2.5 5.3 4.1 4.0 6.2 3.5 3.7 4.5 6.7 5.5 3.6 4.0 5.5 4.2
##  [793] 4.4 4.7 4.1 5.3 5.4 3.9 6.4 4.1 4.7 3.6 4.6 4.3 3.7 4.2 3.1 3.2 3.8 3.7
##  [811] 4.3 4.1 4.9 3.9 4.7 5.6 5.9 5.9 2.5 3.8 3.5 3.2 4.8 4.7 3.4 2.1 5.0 4.9
##  [829] 4.6 2.8 4.2 4.2 3.4 4.0 4.3 3.2 3.9 6.3 6.1 3.7 5.9 4.2 4.4 4.5 5.3 4.1
##  [847] 4.2 3.0 5.9 4.3 6.7 4.0 4.9 3.5 6.6 6.3 4.0 4.1 4.4 3.4 3.8 5.4 5.6 4.4
##  [865] 4.2 5.2 5.3 4.2 5.7 4.9 4.7 5.1 5.0 5.0 4.8 4.3 4.0 5.5 3.5 3.8 4.1 4.5
##  [883] 4.2 5.2 4.9 4.7 4.1 5.3 3.4 5.7 4.9 5.6 5.7 3.5 4.4 4.4 5.4 5.2 5.1 5.0
##  [901] 3.0 4.1 5.8 4.9 4.1 3.7 3.6 4.3 6.0 4.4 3.8 4.8 3.8 5.2 6.4 5.7 6.6 5.4
##  [919] 4.3 3.3 2.4 4.7 3.3 3.7 5.8 3.6 4.3 4.0 5.7 5.9 4.1 4.0 4.2 5.5 4.8 4.3
##  [937] 4.3 4.5 3.5 5.3 3.5 3.8 4.5 4.8 5.2 4.1 6.1 5.2 5.2 6.2 3.8 4.8 5.5 4.7
##  [955] 5.1 3.0 5.3 3.5 4.2 5.1 4.9 2.7 3.0 3.8 4.5 3.0 4.9 5.6 4.6 2.7 5.7 5.3
##  [973] 5.5 3.2 4.3 5.6 4.4 5.1 4.8 4.3 4.7 5.2 3.6 4.5 4.0 4.8 2.9 4.5 3.8 4.5
##  [991] 3.6 5.0 5.2 5.5 5.0 3.1 4.0 3.4 5.2 4.0
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.5 5.5 5.3 5.3 5.1 4.9 4.9 4.8 4.7 4.5
## 
## $jack.boot.se
## [1] 0.9790301
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
## [1] 2.000592
```

What is skew? Skew describes how assymetric a distribution is. A distribution with a positive skew is a distribution that is "slumped over" to the right, with a right tail that is longer than the left tail. Alternatively, a distribution with negative skew has a longer left tail. Here we are just using it for illustration, as a property of a distribution that you may want to estimate using your data.

Lets use 'fitdistr' to fit a gamma distribution to these data. This function is an extremely handy function that takes in your data, the name of the distribution you are fitting, and some starting values (for the estimation optimizer under the hood), and it will return the parameter values (and their standard errors). We will learn in a couple weeks how R is doing this, but for now we will just use it out of the box. (Because we generated the data, we happen to know that the data are gamma distributed. In general we wouldn't know that, and we will see in a second that our assumption about the shape of the data really does make a difference.)


```r
library(MASS)
fit<-fitdistr(original.data,dgamma,list(shape=1,rate=1))
# fit<-fitdistr(original.data,"gamma")
# The second version would also work.
fit
```

```
##     shape       rate  
##   6.008633   9.620888 
##  (2.615944) (4.368612)
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
## [1]  0.7807901 -0.3935708  0.5609180  0.9594277  0.3451040  0.3501046
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
##    [1]  1.2450749199  2.1408259849  1.9234257833  0.6204991372  0.9784909628
##    [6]  1.7469658018  2.1020222650  2.0606040760  1.2888460966  1.2159721681
##   [11]  1.2832578844  1.9644279786  0.9776508766  1.8349556843  2.0287666276
##   [16]  1.2338115220  2.2361852440  2.3540183059  0.7438754393  1.1418874464
##   [21] -1.0531034364  2.0194064025  1.8439942045  1.9932480424  1.2851148478
##   [26]  1.2904458938  1.9510325038 -0.0707135701  2.1506868451  1.2996665667
##   [31]  1.1937548878 -0.0272837620  2.0987051120  2.2834820531  0.6379880495
##   [36]  1.2632024991  0.4664189690  2.0550421793  2.1707541280  1.3629790737
##   [41]  1.8991913040  2.3803107055  2.2811870545  1.2286553966  0.7403299469
##   [46] -0.2466206333  1.9154039329  2.0844580319  1.4434007075  0.5509806276
##   [51]  0.0253317317  1.3226770616  2.2996970790  1.7762462980 -0.3202297721
##   [56] -0.5146614229  1.2636521965  2.1692015476  0.7343927103  2.1596542011
##   [61]  1.2884263941  0.1831809394  2.0759483534  0.6884487451  2.1591326749
##   [66]  1.2156034194  1.8645497091  1.2858651968  0.1864338684  1.5041164911
##   [71]  1.3032559554 -0.2737038088 -0.4472725326  0.0846015601 -1.6144389358
##   [76]  0.4069452671  1.8468318643  0.4501708170  0.7087343488  1.4156628445
##   [81]  2.0096748537  0.0758575237  2.0756855558  0.0295972804  2.1009584488
##   [86]  2.2205108415  1.7057061511  2.0663994918  0.3001293962  0.8631030244
##   [91]  0.1912724294  1.3842675795  2.2613571736  1.3195302734  0.3161299916
##   [96]  1.2229222026  1.8462303936  0.8660692346  0.0470347871 -0.4017067839
##  [101]  2.1499187364  1.9668662429  2.1239476344  1.8073484408  1.8673893273
##  [106]  0.6841892201  1.9768569112  1.9437062149 -0.3881448414  2.1655402652
##  [111]  1.4340179317 -0.1710552476  2.0598030499  0.0807492950  1.2659783492
##  [116]  0.1612197687  1.3869034662  1.9386582437  2.0178499813  2.0863694530
##  [121]  2.0882756962  1.9516191659 -0.1404968136  1.1397499761  1.8765523506
##  [126]  0.7673488143  1.9728412254  0.7406349328  0.6930608859  1.9297246498
##  [131]  0.1093782702  2.0005915954 -0.9598212601  2.1043737627  1.9469931115
##  [136]  1.2844529569  0.6971496755 -0.7928339409  1.4151276172  1.9835433953
##  [141]  0.1150390091  1.9621450102  2.3208747506 -0.0233724237  1.8708597843
##  [146]  0.1665153811  1.2367688397  2.5624377581  1.0629099592  1.0199275205
##  [151]  2.2314774714 -0.1340202842  2.0311710466 -2.2239342648  0.1515888606
##  [156]  2.0160134074  0.2034785226  1.3933623078  0.2047595791  1.1930729648
##  [161]  1.3015853268  2.0279555077  1.3013388474  0.2991244066  0.6217325485
##  [166]  1.2388351328 -0.0546020395  1.7154332601  0.4382381434  0.5353015720
##  [171] -0.0101791827  2.0005915954  2.2732856299  1.9888682246 -0.4591184613
##  [176] -0.0453495282  2.0957672671  1.9808955795  2.0376131554  0.3040043562
##  [181]  1.9655795356  0.7327507521  2.2250542141  0.4958113020  0.6819530495
##  [186]  0.8579993871  1.6850636303  2.1771671538 -0.4373505458 -0.0006212691
##  [191]  1.2457886019  1.9891779676  1.1705266438  2.2116103105  1.3626694822
##  [196]  1.8579092357  2.3297887769  2.3811047297  1.8447883023  0.0505511973
##  [201]  1.9596664276  1.2117039568  0.7851670303  2.1549177509  1.1278497321
##  [206]  2.0038637388  2.0455840907  1.1611670499  2.2631349817  2.1325161613
##  [211]  1.3342341048  2.0347256731  1.2046676910  0.1848733149  1.3283796144
##  [216] -0.3796952634  1.3000386166  0.6409642925  2.0029465036 -0.5222917854
##  [221]  0.0591175058  1.9675806947  1.7626459909  1.3361623830  1.9836016541
##  [226]  0.9669397106  0.7458505669 -0.0262917131  0.6154010421  0.7232441331
##  [231] -0.8351219687  0.7778029388  1.4264776066  1.7450836067  2.2861068390
##  [236] -0.9131014324 -0.5855833665  0.0119672764 -0.3260558626 -0.3952937285
##  [241]  1.6947710028 -0.8396839974  2.1509134239  2.0937721069  2.2014109058
##  [246]  0.6731809138  1.2486595260  2.0931683784  1.1957488371  2.0066190656
##  [251]  1.8071168898  2.2323119902  2.1325421887 -0.2545434975  2.0625125586
##  [256]  1.3900267813 -0.1678484192 -0.6477670715  1.3419250815  0.7387306570
##  [261]  1.3164180723  2.1487093582  2.1731459709  0.5039113874  0.1401209256
##  [266]  0.2529929687  2.0757596065  1.9922490630  0.6362563023 -0.0376265806
##  [271] -0.3497988606  1.2442831753  2.1263765457  1.2002552965  1.7368992971
##  [276]  2.1569294621  2.1441516708  0.2327590679  1.9053979645  1.3405519037
##  [281]  1.3469955883  0.2504330218  1.3181568461 -0.4058590698  1.1919832101
##  [286]  1.1877895965  0.8476025985  1.2053583385  2.1607799479  2.2949776532
##  [291]  1.9429751795 -0.6350388758  2.0076847483  2.2164080870  2.1441997886
##  [296]  1.3574957290  2.0657068035  1.1572748621  2.0572863925  1.2386854273
##  [301] -0.4449557724  0.7462837649  2.1104688012  1.8113759131  1.3038637648
##  [306]  2.3977122962  0.5949759309  1.3026517262  0.6525801121  0.9071502607
##  [311]  1.9199889789  1.8820101314  1.2332902460  0.8110253944  1.2766881954
##  [316]  1.3119127454  1.9352828614  1.8356551477  2.1960128182  2.1983243336
##  [321]  0.2524785865 -0.1741344298  0.3579570595  2.3030434235  0.7417098430
##  [326]  2.1233397614  0.5372339254  0.2744420264  0.6575691733  0.2947229861
##  [331]  1.7099079409  2.5850699487  0.2999504481 -0.5602448039  2.0565070352
##  [336]  0.3254566524  2.1947832220 -0.2898719024  1.2147335169  1.2722078242
##  [341]  0.4692250741  1.3805025720  1.9009748090  2.0330063722  1.3288954588
##  [346] -0.0233724237  0.6388985017 -0.4420511800  1.9669528025  1.3181568461
##  [351]  0.0090159031  2.1628062300  1.1757024611  0.7751636013  1.7629882794
##  [356]  0.7446390901  2.3409436088  0.3111604504  1.9082355475  0.6879038473
##  [361]  0.6666731309  1.3001824094  1.2470287120  2.3780576177  1.2067331955
##  [366]  1.4179229574  1.6633965048  1.4297876348  1.0234401108  1.4777290244
##  [371]  0.6507019548 -0.0247925754  2.0315671893  1.5207103925 -0.8183805009
##  [376]  0.0514851338  1.3705007037 -0.2226125808  1.9063847088  1.2977888987
##  [381]  0.5524541134  1.2897460441  1.9398292797  2.1849556707  0.7206129139
##  [386]  1.8477168261 -1.0419819843 -0.4627116169  2.2421717072 -0.0053584085
##  [391]  2.3165810203 -0.1687420523  2.0473466490  2.1298247473  0.0899812990
##  [396]  0.5799971213  1.9669426918  2.2118808920  2.2786485651  0.1212073698
##  [401]  0.1245705982  0.6845324851 -0.0429495945  2.1704770499  0.2358129756
##  [406]  0.8971180244  1.1248493848  0.3759227602  1.1532365486  1.1275947862
##  [411] -0.1158893749  0.8714908965 -0.1979185338  2.1988830733  1.3082632108
##  [416]  2.2490398212  1.1504758315  0.5498593822  0.7934838125  2.0546806900
##  [421]  2.0234050792  2.0580855667  1.7030323816  0.7349157884  0.7251947245
##  [426]  2.1925919097 -0.0302001357 -0.3610717558  1.2258713879  2.6326493246
##  [431]  0.7480446753  2.2335412256  1.7446716746  1.8455778990 -0.2783168336
##  [436]  1.2510228080  2.0265206324  1.3507053378  2.1638861632 -0.2080500969
##  [441]  1.3102174402  1.9109253158  1.9707144827  0.8037404451  1.9499403700
##  [446]  0.1556154136  2.3459141584  1.3581038937  0.2429549209  0.1350831882
##  [451]  1.8176154332  0.6955243641  2.2275578255 -0.4060525741  1.3361309095
##  [456]  2.1185934439  0.1873042359  2.0533801231  0.6621329756  1.1670205606
##  [461]  0.6787861956  2.0631997790  1.9195723778  1.0408717945  1.8564413287
##  [466]  2.2096363912  2.1705708233 -0.7243462934  1.1471990057  2.3407839222
##  [471]  0.8144260319  2.0773781507  0.7429437384  0.5982799217  1.9688896142
##  [476]  2.0490973474  2.1843035244  1.9181554149  2.0867426924  2.1964262637
##  [481]  1.1818344789 -0.0946687468  1.7882146554  1.8147273814  1.8398196205
##  [486]  1.2536067429  1.1290172752  0.7023855962  0.6150118773  2.0630195277
##  [491]  0.6783902860  0.1759358532  2.1677873062  0.6859006447  2.1803827419
##  [496]  2.1263765457  0.6222731436  1.9984655721  1.2314978496  2.1807858767
##  [501]  1.9941334820  1.3763564116  0.6736911781  1.9757981810  1.3504332973
##  [506]  1.7616953127  2.2334140155  2.0352247260  2.0306134146  0.1890809687
##  [511]  1.2526006089  2.3160775707  2.1538481331 -0.1128936059  0.6002249362
##  [516]  2.0437419405  0.1539653122  0.8357215018  1.2172867926  0.1808578172
##  [521]  2.1029902110  1.2851508517  1.3248925054  0.8725973516  1.9671063712
##  [526]  1.9870860514  0.0237615340  0.6816062068  0.0530965979  2.0022693137
##  [531]  0.7426684078  2.4530196145  0.7382632832  1.6795559179  0.6369848972
##  [536]  1.9401995273  1.2997666989  2.3472891599  1.8634806773  1.2731008143
##  [541]  1.8625677986  1.0655206501  1.9255997260 -1.3423316010  1.1807511403
##  [546]  0.9634471344  0.6609603016  1.2395746254 -0.3170120203  1.3334549534
##  [551] -0.7651044537  2.2586493684  0.8558370052  0.6777852104  0.7635215389
##  [556] -0.0200727712  0.7186004448  1.2584015270  1.2932259544  0.6061670729
##  [561]  2.0560939175 -0.1102478509  2.3679810034  0.4492089543  1.3335992950
##  [566]  2.1186003650  1.1984154239  1.9699626637  1.3049088477  1.3010458384
##  [571]  2.0747061896  0.7057420099  2.1817104014  2.2001512850  1.0908734747
##  [576]  0.2312370867 -0.2082739983  2.2495811029  2.0911596459  0.5620741299
##  [581]  1.8301861694  0.3251686896  0.7355739546  1.1816865962  1.2270485620
##  [586]  0.6660415552 -0.5670442721  0.4909117295  2.0319975608 -0.3559091092
##  [591]  1.2650949728  1.9891779676  2.3028127192 -1.3649849141 -0.1429242098
##  [596]  1.1972203590  1.1581121518  0.1615226089  2.0564315909  1.9672959503
##  [601] -0.7115847654  1.1723052615  1.9782732924  1.9967150733  2.1384636018
##  [606]  1.3236964856  1.9915440261  1.6524055044  1.1756754298  1.8499641414
##  [611]  1.2501339035  2.1844832825  2.2289961934  1.9211973745  0.7511537734
##  [616]  0.0211119437 -0.1305833619  2.1743218173  1.4179347248  0.2113244085
##  [621]  0.3849064952  0.2135568595 -0.0155143636  0.1659437704  1.1567384489
##  [626] -0.0329618218  2.1010285820  2.0533573090  1.3252394736  2.0061158718
##  [631]  1.9022968862  0.5382180453  1.9339396563  2.0906004688 -0.1146341605
##  [636]  2.1536975094  0.7715918252  1.1191929544  2.1192967240 -0.0331048460
##  [641]  2.1067592938  1.8636022307  2.1668973896  2.0380362774  2.0828410349
##  [646]  1.2147070029  0.0551818024  2.3241087983 -0.0986334633 -0.2834322089
##  [651]  2.1134585340  0.0018529676  2.0701853278  0.4277250445  1.9568922354
##  [656]  2.2542877282  1.7127326740  1.9294226108  0.0012378050  0.8542067019
##  [661]  1.7988814502  2.3172625220  1.2312189027  2.0529835990  0.2862274996
##  [666]  0.8147039018  2.0159374031  2.0565070352  0.6157216492  1.2483568218
##  [671]  1.8553624128 -0.1037346546 -0.1682148425  1.9497403525 -0.3968965764
##  [676]  2.0983904193  1.2621880872  1.7989088134  1.1529219669  0.3074159956
##  [681]  0.9721506202  2.0341881850  2.0456309725 -0.2463809367  0.1499431280
##  [686]  1.7486717188  1.9984432513  0.6654032488  0.4829514326 -0.2585316355
##  [691]  1.2256412648  1.9776946020  1.2850809869  0.6215117771  2.0652472134
##  [696]  1.2650949728  2.1612707987  1.9077591121  1.2766881954  2.0393853119
##  [701]  2.0793886946  1.2848552715  1.9538103144  1.3019013436  0.3134424672
##  [706]  2.2291039578  0.2977564491  0.2977564491  1.9890576534  0.1192447704
##  [711]  1.3421565213 -0.3304228313  2.4092516680  1.1705261759  1.2546016471
##  [716]  2.1191088683  1.9970874040  2.1425401069  0.0282511466  1.8726963468
##  [721]  2.0771495635 -0.0379495130  1.3067292445  1.9984325769  1.9913177507
##  [726]  1.2214868623 -0.9744151015 -0.0087507335  0.7565896650  2.4063172623
##  [731]  0.8339439223  1.8482736573  1.1352365249  2.4213755026  2.3423264826
##  [736]  0.4981127981  1.9000878489  0.2902979515  1.2446055455  2.2453517453
##  [741]  0.0413960906 -1.1909215078  1.2433290554  1.9528508334  2.0004863160
##  [746]  1.3661686223  0.2950612549 -0.0821909867  0.6435164841 -0.2061039429
##  [751]  1.2100002882  1.0601705091  0.2483726792  0.3981801259  1.1798222317
##  [756]  2.1841782482  1.2210976264  2.2916741874  2.3876292614  0.8689365210
##  [761]  1.2305190710 -0.2154357228  1.2271354418  0.3458527532  0.1609344329
##  [766] -0.7258515499 -0.7839246412  1.1809832396  2.0484258713 -0.6141116120
##  [771]  0.5030892074 -0.4444216667  1.3134588063  0.2296027793  1.8815772923
##  [776]  0.0846603525  1.2377456320  0.0628364635 -0.2543310685 -0.8387708617
##  [781]  2.0864605980  0.6814078633  1.0917656058  1.2588063191 -0.1508399561
##  [786]  1.9611791106 -0.2055115547  0.6647517793  0.2903940223 -0.1981785916
##  [791]  2.0005915954  0.2197082714  0.7376762926 -0.4642843259  0.3887548408
##  [796]  0.7618473936  1.9635673793  1.2920049301 -0.6058551281  0.1730133222
##  [801]  0.4097768103  0.6448749263  0.7618987554  1.8955085974  1.2800484889
##  [806]  0.3652870096  0.8144845972  1.2358831419  0.3147590058  1.8547954818
##  [811]  1.8243187020  2.1933548198  1.2847241132 -0.0325850523  1.8528474365
##  [816] -0.1849671237  0.5825025245  2.0122453270  1.9299998393  1.1420795744
##  [821]  1.9912037836 -0.8406335201  0.3044251029  0.5178783883  1.7993651405
##  [826] -0.5673270696  0.9391947749  1.2969008077  1.9954007464  2.0130567093
##  [831]  2.2057730650  1.3305994986  2.0312605986  2.5394232935  0.6541041624
##  [836]  1.9893825818  1.1909604233  1.9392160522  1.8716397783  1.2699504618
##  [841]  0.3270405237 -0.6773970391  0.2998107079 -0.9670664600  0.7624388904
##  [846] -0.0253049028  2.2972628398  1.9685345667  2.0748841725  0.9887422723
##  [851]  2.2522814912  2.4037972387  2.0830382300 -0.3770450667  2.3982994008
##  [856]  0.6541041624  1.2405217249 -0.3491155388  2.0673491168  2.3462160826
##  [861]  0.1477544763  0.3654435621  1.1422322992  1.2470960315  0.7209916330
##  [866]  0.5504803285  2.0565070352  1.1779342042  0.3845410024  1.3833660473
##  [871]  2.3338706486  1.2322511263  1.8410054086  0.3425614817  2.1882124851
##  [876] -0.9131014324  2.0527930477 -0.3542241570 -1.2072561615  0.6973088670
##  [881]  2.1926664606 -0.0977380599  2.2374835895  1.3195302734  0.7230308754
##  [886] -0.6767002289  2.1435145788  2.1936048211  0.7172653055  1.0393331126
##  [891]  0.7599797142 -0.4155074612  1.0642381526  2.3301864461  2.2396649905
##  [896]  1.9865695772  1.1947559311  1.8066970701 -0.4790175815 -0.1013775790
##  [901]  2.2727370799  1.2795036394  2.2386048799  2.0417851383  2.3177469778
##  [906]  0.4047097279  2.2987459096 -0.8317919261  2.1953154294  2.0358337574
##  [911]  2.2963857746  1.1507185764 -0.0074676030  1.8736423912  1.1314585765
##  [916]  0.8374614605 -0.3558900982  1.8584022392 -0.0020216386  2.2936198689
##  [921]  0.3262307582  1.9723570719  1.9150690012 -0.1522406975 -0.3869210116
##  [926]  2.1495298075  2.1107435299  0.8639659825  1.8129597353  1.9326022935
##  [931]  0.7540386485  2.0419012654  2.1758145167  0.5730658970  2.1085087263
##  [936]  0.7149424882  1.1063029736  1.1723052615  0.7633174093  2.0294025268
##  [941]  0.3306421538  1.1985552184  0.7332826965  1.9342787404  0.6750386008
##  [946]  1.3720940472  1.1740859000 -0.1522631180  1.5360738033  2.1353064536
##  [951]  2.2071680071  2.0855467548  2.0751264276  2.0606362120 -0.0465879282
##  [956]  2.3184921109  0.6451113680  1.9986702520  2.1572917682  1.2428816839
##  [961]  1.7768197168  0.4136358400 -0.5563031545 -0.6215109358  1.3203466766
##  [966] -0.3347082881 -0.1063632998  1.7341255394  0.7996771814  2.2261915425
##  [971]  1.2347757221  0.1536808041  1.1759097351  0.7022419184  2.1422744206
##  [976]  1.8643565838  1.9585998551  2.2159395560 -0.6304807208  1.1661663222
##  [981]  1.2716727969  1.8133235776  0.2152144362  1.2531835112  1.9553994405
##  [986]  0.3141875949  1.2696148997 -0.4521766087  0.5201596266  2.2916248692
##  [991]  0.0283733681  0.7600915493  1.9577658698  2.0890497887 -0.3077903123
##  [996]  0.6968762978  1.2451010664  1.3858702461  1.0774924707  2.2532205673
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
##       mean          sd    
##   0.62454389   0.30583364 
##  (0.09671309) (0.06838873)
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
## [1] -0.2333330 -0.1349436 -0.4423278  0.3395940 -0.4811803 -1.0898689
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
##  [1] 5.000000 4.888889 4.777778 4.666667 4.555556 4.444444 4.333333 4.222222
##  [9] 4.111111 4.000000
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
## [1] -0.0142
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9097615
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
##     original      bias    std. error
## t1*      4.5 0.008908909   0.9032621
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 1 3 4 8 9 
## 2 1 1 2 3 1
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
## [1] 0.0204
```

```r
se.boot
```

```
## [1] 0.9115453
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

