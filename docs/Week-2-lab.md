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
## 0 1 3 4 5 7 8 9 
## 2 1 1 1 1 2 1 1
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
## [1] -0.0298
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
## [1] 2.638249
```

```r
UL.boot
```

```
## [1] 6.302151
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.7   6.4
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
##    [1] 4.6 7.4 4.0 6.5 3.8 5.0 4.4 4.7 5.0 5.4 5.0 3.6 5.4 4.5 5.6 4.1 3.9 4.4
##   [19] 3.8 4.6 5.9 4.8 4.7 3.4 4.9 3.4 6.3 5.2 4.2 3.2 4.0 3.7 4.7 3.6 4.8 4.0
##   [37] 3.9 4.1 2.6 4.4 4.8 5.3 3.6 3.4 4.5 3.6 4.2 6.0 3.5 3.4 5.1 5.1 4.6 4.5
##   [55] 4.2 4.4 4.2 3.4 1.8 4.8 4.9 3.8 5.4 4.8 4.5 4.8 5.9 4.5 4.7 4.7 3.4 5.0
##   [73] 4.6 4.6 4.5 4.2 3.9 4.4 5.2 5.7 5.6 3.5 4.7 4.1 3.9 6.2 5.0 5.0 2.8 5.1
##   [91] 3.2 3.3 4.6 3.6 3.5 5.1 5.3 5.1 5.6 3.8 3.5 3.8 5.3 4.5 4.4 5.9 4.3 4.2
##  [109] 3.6 4.9 2.9 4.8 4.4 3.8 2.6 4.1 6.4 4.8 4.1 4.4 3.9 6.0 5.3 2.7 5.6 4.7
##  [127] 4.3 4.5 3.9 4.7 5.1 4.4 3.4 4.8 4.0 5.2 3.5 3.7 3.1 4.1 4.8 5.1 5.5 4.9
##  [145] 3.0 5.3 5.6 5.3 5.3 5.0 4.9 3.8 4.9 4.5 5.2 5.4 6.0 4.4 4.9 3.9 6.3 3.7
##  [163] 5.3 3.9 4.0 4.2 4.1 2.5 4.3 3.4 3.9 5.3 5.4 3.9 4.5 4.5 3.9 4.5 4.7 3.1
##  [181] 6.6 4.0 3.8 3.3 4.1 4.7 3.9 6.1 5.6 6.5 3.7 4.2 3.1 4.1 3.9 3.8 5.8 4.4
##  [199] 6.0 3.5 5.2 3.0 5.2 2.9 4.6 4.0 4.6 4.1 4.7 3.0 5.0 3.4 5.0 4.7 2.9 3.4
##  [217] 4.8 5.0 4.1 5.6 3.7 3.9 6.1 5.2 6.9 3.4 3.9 4.1 5.4 4.4 4.1 3.0 3.8 4.4
##  [235] 4.3 5.1 4.9 5.1 5.2 3.8 4.5 3.3 4.6 4.8 3.8 4.4 5.0 5.9 6.2 4.4 6.0 4.3
##  [253] 4.6 4.2 3.5 4.1 4.2 1.9 2.6 5.6 4.3 4.1 5.0 3.6 4.6 4.0 4.9 5.1 3.9 4.0
##  [271] 6.1 4.7 6.2 5.3 3.5 3.9 4.7 6.0 4.7 5.3 2.0 4.6 3.8 5.5 3.9 4.4 5.9 4.9
##  [289] 6.1 4.3 5.2 4.9 5.5 5.4 3.2 5.3 4.4 5.4 3.9 3.7 5.2 5.0 3.5 2.8 4.1 3.4
##  [307] 5.0 4.9 6.3 4.5 5.3 3.7 4.7 5.0 3.0 3.1 3.8 2.9 4.1 5.0 5.3 4.2 2.6 5.6
##  [325] 4.7 3.8 3.0 5.5 4.6 5.4 5.0 6.0 5.5 5.1 4.3 4.7 5.5 5.8 3.4 4.6 4.0 5.5
##  [343] 4.1 4.9 5.7 5.1 4.1 3.2 4.6 4.7 3.6 5.5 4.5 3.8 3.8 5.0 3.8 5.2 3.2 3.0
##  [361] 4.4 5.0 2.7 5.9 5.4 2.9 4.5 4.7 4.3 3.6 5.4 6.7 4.1 3.4 3.9 2.7 4.6 3.3
##  [379] 2.4 2.1 3.3 5.7 4.2 4.5 4.6 5.9 5.7 2.5 4.5 5.8 5.5 4.3 4.3 5.5 4.4 5.9
##  [397] 5.0 6.4 4.0 5.4 3.4 4.4 2.7 5.3 2.6 4.9 4.7 6.1 4.7 4.0 3.7 3.6 4.4 5.2
##  [415] 4.4 5.0 4.8 2.8 3.0 3.0 4.7 5.5 4.0 4.4 3.8 5.3 3.0 3.3 3.8 6.0 4.9 5.4
##  [433] 4.1 4.4 4.0 5.2 3.1 5.2 3.0 5.3 4.6 5.4 2.9 4.5 4.1 5.0 3.9 6.1 5.0 4.3
##  [451] 4.1 6.0 3.5 4.9 5.5 3.8 3.8 5.0 3.1 4.9 2.5 1.0 2.2 4.2 3.8 6.5 4.2 2.5
##  [469] 3.4 4.5 4.7 3.6 3.7 4.9 4.2 3.6 6.4 5.2 4.0 2.9 4.5 4.9 4.3 4.7 5.2 4.8
##  [487] 4.8 5.3 5.0 6.0 5.8 4.2 4.3 2.7 3.1 3.9 5.5 3.8 4.0 5.3 4.6 4.8 5.3 4.3
##  [505] 4.4 3.7 4.6 5.0 6.2 4.1 5.4 5.4 5.2 5.8 6.8 5.3 5.5 3.0 4.3 4.1 3.8 4.8
##  [523] 2.8 4.0 3.4 5.4 5.2 4.1 3.5 3.4 4.5 4.0 3.4 3.6 5.4 5.6 5.3 4.0 3.6 4.1
##  [541] 3.2 3.3 3.8 4.3 3.5 4.4 5.5 5.1 5.3 3.2 6.2 3.4 4.3 5.2 3.9 4.7 5.9 4.3
##  [559] 3.8 4.0 4.8 5.5 5.0 5.6 3.2 4.2 5.0 5.4 3.8 5.3 4.3 4.2 4.6 3.9 5.2 5.5
##  [577] 6.2 3.6 2.8 5.9 4.1 4.7 4.4 3.9 5.3 3.6 4.4 4.6 2.6 4.0 4.4 6.6 2.9 3.4
##  [595] 4.7 4.1 4.6 6.0 5.2 3.8 5.3 5.0 2.2 6.3 4.2 4.7 3.9 4.9 4.3 3.8 4.0 4.6
##  [613] 5.2 4.6 3.6 4.3 4.6 6.3 4.1 3.3 4.3 4.5 3.5 3.9 3.9 6.5 3.0 3.5 4.4 4.5
##  [631] 5.0 6.4 5.9 3.4 4.6 3.1 4.2 3.1 5.6 5.0 3.7 6.1 3.9 3.1 4.2 4.1 3.6 4.6
##  [649] 4.4 5.7 4.4 4.2 3.6 5.2 3.7 4.8 5.1 4.7 3.2 4.7 3.9 4.1 5.8 4.5 4.8 4.9
##  [667] 4.3 4.4 3.8 4.0 3.1 4.4 4.6 4.3 5.0 4.2 3.2 2.6 6.2 3.4 4.9 4.5 6.5 4.1
##  [685] 4.2 3.8 4.0 3.9 4.1 6.3 4.1 5.1 4.4 4.2 4.4 5.6 5.4 6.3 4.6 5.4 6.1 5.7
##  [703] 3.8 3.6 3.7 5.3 4.7 5.0 3.9 3.9 4.2 2.9 3.9 4.8 5.5 3.4 5.1 5.5 5.6 5.6
##  [721] 4.8 4.9 4.9 4.7 3.8 4.5 5.3 5.1 5.2 2.3 4.7 4.3 5.6 5.0 4.9 5.3 5.6 4.3
##  [739] 4.3 5.9 4.5 3.2 4.3 3.4 5.4 2.9 6.1 3.8 5.0 5.9 3.8 5.3 5.1 4.5 5.1 4.1
##  [757] 4.9 4.4 3.7 4.8 5.8 5.3 4.0 5.7 3.9 4.7 4.0 7.2 4.8 4.2 4.6 3.5 3.6 4.8
##  [775] 4.3 6.1 6.1 4.4 3.0 3.8 5.6 4.0 4.3 3.6 3.2 3.1 6.7 3.8 4.2 3.4 4.2 4.7
##  [793] 4.1 5.0 4.1 4.8 3.1 4.2 4.5 4.6 2.4 4.3 4.8 3.8 4.3 4.7 4.0 4.3 1.8 4.1
##  [811] 5.4 5.4 1.8 4.6 5.3 2.8 3.9 5.4 3.7 4.1 3.2 5.6 5.8 5.1 3.8 5.5 4.9 3.9
##  [829] 4.9 4.9 4.9 3.5 4.2 3.7 4.4 4.9 2.9 3.7 5.6 4.0 3.9 4.3 5.9 5.4 3.8 4.6
##  [847] 5.1 6.1 4.2 3.4 3.7 2.1 4.2 4.1 2.9 6.2 4.4 4.2 3.7 3.4 3.6 5.4 2.7 5.2
##  [865] 4.7 3.0 4.7 5.0 4.9 5.7 6.1 5.7 4.2 4.8 3.8 3.7 2.3 4.1 5.3 5.2 3.8 2.9
##  [883] 3.9 4.2 3.0 4.1 5.7 6.4 3.7 5.7 5.7 6.9 2.8 3.9 4.1 4.8 6.0 4.3 2.7 4.0
##  [901] 3.8 4.6 5.8 5.2 5.2 5.2 4.0 5.0 5.6 4.0 4.5 4.6 4.4 4.4 4.7 3.2 4.2 5.4
##  [919] 4.8 4.5 5.4 6.2 4.2 3.7 4.7 4.3 3.5 2.5 4.5 3.8 4.9 5.4 5.3 3.4 3.5 4.6
##  [937] 4.1 4.6 5.1 5.4 5.3 4.8 4.2 2.5 5.9 4.2 4.3 3.3 4.9 4.8 4.6 4.4 4.4 2.4
##  [955] 4.2 4.7 3.9 4.3 6.0 4.8 5.7 4.1 4.8 5.5 4.9 5.8 4.7 5.8 3.4 3.3 3.9 5.0
##  [973] 4.6 5.6 3.4 4.5 3.2 4.7 4.9 3.7 4.5 3.5 2.7 3.4 3.9 5.4 5.3 6.2 4.4 4.3
##  [991] 4.9 3.4 3.6 3.4 4.5 4.5 6.1 3.9 3.7 5.5
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
##   2.6   6.3
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
##    [1] 4.7 3.1 5.3 5.3 6.1 3.2 4.1 6.1 4.7 4.9 5.0 4.7 3.3 4.8 4.2 4.4 4.6 6.6
##   [19] 6.0 6.5 3.6 6.0 6.0 4.2 5.7 5.9 4.5 6.3 3.5 4.9 2.8 4.6 4.4 5.9 5.6 4.6
##   [37] 3.1 3.0 4.6 4.6 3.8 5.5 5.8 4.8 4.6 6.3 4.4 5.2 6.0 5.7 4.2 2.6 5.7 3.0
##   [55] 5.8 4.9 3.9 3.8 3.2 3.7 4.6 4.6 3.5 4.1 4.5 3.3 5.0 5.0 4.9 4.4 5.5 4.7
##   [73] 2.7 5.5 3.9 3.9 4.1 6.0 2.3 3.3 4.0 4.9 6.3 3.7 4.6 4.5 3.6 5.1 4.3 5.0
##   [91] 3.4 5.5 4.4 5.0 5.9 3.9 4.9 3.2 5.1 4.4 3.7 3.3 4.4 4.9 4.5 5.6 4.8 5.1
##  [109] 2.6 4.4 3.7 2.4 5.2 4.0 4.9 2.5 4.1 2.9 4.6 5.3 5.0 4.5 3.5 6.2 4.1 5.6
##  [127] 2.3 5.1 4.3 5.7 4.7 5.0 4.8 4.6 5.5 4.2 5.4 5.1 3.3 6.0 5.6 5.2 4.7 5.4
##  [145] 5.1 4.8 6.5 3.7 4.0 4.0 4.3 4.9 4.4 4.0 4.5 3.6 4.0 4.3 4.9 5.2 5.5 4.9
##  [163] 4.9 3.9 4.8 4.5 3.3 4.3 2.9 6.0 3.3 4.3 5.5 5.4 3.8 4.4 5.6 4.8 5.5 4.3
##  [181] 5.2 6.0 4.1 4.8 4.1 3.4 3.1 3.9 5.2 5.3 4.8 4.9 6.0 5.4 4.3 3.6 4.4 5.0
##  [199] 5.8 3.2 4.8 3.8 3.7 4.5 4.2 3.9 2.9 4.7 6.0 3.6 5.5 4.9 6.0 3.5 5.0 3.5
##  [217] 5.7 2.7 4.3 3.9 5.2 4.2 3.7 3.8 3.9 5.1 4.6 3.9 5.5 4.8 4.8 3.6 4.2 4.3
##  [235] 6.4 5.3 4.7 5.3 4.7 4.7 4.4 6.5 3.1 5.3 4.0 3.5 5.4 5.9 4.7 2.9 5.4 5.7
##  [253] 5.5 3.7 4.4 6.1 3.9 3.9 2.7 4.8 5.8 3.7 3.3 4.8 5.4 3.6 5.1 3.8 6.3 4.4
##  [271] 3.5 4.7 3.2 5.1 3.5 4.0 4.8 5.1 3.5 5.7 4.9 5.1 5.0 4.8 4.3 5.3 5.0 3.6
##  [289] 5.2 3.6 3.9 4.8 4.4 5.4 2.0 4.6 2.6 6.6 4.2 3.6 3.2 6.5 4.7 4.5 4.9 3.9
##  [307] 3.0 3.4 3.8 4.6 4.4 4.2 6.9 5.8 3.5 5.3 4.4 3.9 5.1 5.3 6.1 4.5 4.0 4.8
##  [325] 4.2 4.5 3.9 4.4 3.6 3.2 4.4 4.0 5.2 3.6 4.9 4.5 5.4 3.7 3.8 4.3 5.9 5.7
##  [343] 4.8 5.2 3.2 4.6 3.8 6.1 3.8 4.1 3.7 4.3 5.1 4.5 7.0 5.2 4.8 5.6 4.3 4.7
##  [361] 3.9 3.3 4.3 5.1 5.5 5.8 3.5 4.4 4.6 5.0 4.8 4.8 4.3 4.3 5.2 4.1 4.4 5.0
##  [379] 4.8 5.0 5.4 3.5 3.9 2.3 4.3 5.7 3.9 3.6 2.9 3.5 5.6 3.4 5.6 3.7 5.1 2.9
##  [397] 5.2 4.8 4.4 5.9 4.0 4.5 2.8 4.6 3.5 5.6 4.5 4.6 4.1 5.5 4.5 2.0 4.2 4.5
##  [415] 4.5 4.4 5.2 4.6 5.6 3.0 5.2 1.4 3.1 4.2 3.2 5.7 4.2 4.9 3.8 3.8 3.5 3.7
##  [433] 4.6 2.6 2.9 5.8 3.7 5.3 4.7 4.6 5.1 4.7 5.5 4.0 5.8 6.9 4.4 5.3 4.4 4.8
##  [451] 4.0 4.3 4.3 3.4 4.5 6.7 4.2 5.2 3.7 3.9 5.7 4.4 3.4 3.4 3.4 3.9 4.2 5.3
##  [469] 3.7 5.0 3.7 4.2 4.0 4.2 4.7 5.0 3.6 4.2 3.6 5.4 4.6 4.1 4.4 4.7 4.8 3.2
##  [487] 3.5 4.2 5.9 3.4 5.2 3.6 3.6 4.0 5.8 2.7 3.3 4.5 4.6 6.9 6.1 4.0 5.1 4.7
##  [505] 4.3 4.1 4.5 5.2 3.6 5.5 3.9 4.8 5.3 4.6 3.9 4.3 3.7 4.3 5.4 4.7 3.7 4.4
##  [523] 4.0 4.2 4.6 3.8 5.0 5.3 5.7 4.8 3.8 3.9 4.4 3.5 3.9 5.4 6.0 4.5 4.5 4.6
##  [541] 3.4 4.2 4.1 3.4 4.1 5.0 6.5 4.6 4.4 3.1 3.7 4.4 4.5 5.4 4.0 5.0 4.7 5.9
##  [559] 4.7 3.8 4.6 4.3 5.2 4.5 4.7 4.0 2.7 5.4 4.3 4.9 5.5 3.8 3.9 4.7 4.5 3.7
##  [577] 3.7 3.4 4.1 5.0 5.4 2.8 6.1 3.2 3.9 3.3 3.9 6.0 6.0 4.4 4.4 3.2 5.2 4.7
##  [595] 4.9 3.6 3.3 4.3 4.2 4.6 5.8 5.6 4.6 5.1 4.3 5.0 4.6 3.6 5.5 5.0 3.4 4.9
##  [613] 6.3 4.6 4.9 4.8 5.0 4.1 4.1 5.4 5.4 4.4 3.1 3.7 3.7 4.4 4.1 3.7 4.5 5.4
##  [631] 5.9 3.9 4.1 6.5 5.3 6.2 4.6 4.1 3.7 5.1 4.4 5.2 4.8 5.9 3.4 3.5 2.8 5.5
##  [649] 4.6 5.9 5.5 4.7 3.7 4.9 3.1 3.4 3.4 5.6 4.1 4.2 3.6 4.1 5.1 3.9 4.0 4.2
##  [667] 1.8 3.4 6.4 3.2 3.8 4.9 4.8 3.4 2.2 6.2 5.1 4.6 5.4 4.9 3.8 3.9 6.4 5.1
##  [685] 3.9 3.8 3.6 3.2 3.9 4.6 3.8 5.2 5.1 4.2 3.2 4.0 4.4 5.2 4.3 4.7 4.5 4.8
##  [703] 4.9 4.4 4.6 4.2 5.0 5.4 4.6 3.9 4.0 4.1 5.5 4.8 4.9 4.2 3.8 5.9 4.9 3.7
##  [721] 3.3 3.6 5.0 4.2 3.9 3.6 3.8 4.4 3.2 4.3 4.9 3.8 5.4 4.0 3.4 5.6 4.8 3.6
##  [739] 6.0 4.7 4.0 3.7 5.5 4.8 3.7 5.2 4.9 3.6 4.3 5.5 4.4 2.8 5.1 5.2 5.4 4.7
##  [757] 3.9 2.1 4.5 4.0 7.2 4.0 4.7 5.1 3.6 4.4 5.3 4.2 5.1 3.8 4.8 3.6 4.9 3.7
##  [775] 3.6 4.3 2.9 4.8 5.2 3.2 4.0 3.2 4.0 5.6 4.6 5.9 5.0 3.9 4.5 5.2 3.9 3.9
##  [793] 6.4 5.1 4.1 2.2 4.5 4.6 4.6 3.6 4.3 4.6 1.8 3.2 3.4 4.6 4.6 3.9 5.1 4.4
##  [811] 4.4 3.5 4.9 5.3 3.0 3.8 6.7 5.9 5.3 5.8 3.9 4.6 3.5 4.1 3.9 4.6 3.6 4.0
##  [829] 4.6 1.8 4.9 4.7 2.6 3.4 5.7 4.5 3.7 4.4 4.5 3.4 4.1 4.1 3.5 2.8 5.0 6.2
##  [847] 5.5 3.4 3.7 5.7 2.9 3.2 3.2 3.2 4.5 7.5 3.5 5.1 5.0 4.7 3.0 4.1 4.0 3.9
##  [865] 5.0 3.4 4.0 3.4 3.9 2.5 3.9 4.8 6.6 3.6 4.5 4.4 4.2 4.9 4.9 3.2 4.7 4.2
##  [883] 5.2 4.1 3.1 5.0 4.6 3.7 5.5 4.1 3.5 3.8 4.4 5.7 4.6 5.0 3.9 4.7 4.4 2.8
##  [901] 5.8 3.2 4.6 5.2 4.4 4.7 3.3 5.7 6.3 3.8 4.2 5.3 4.5 5.5 5.8 4.4 2.9 3.5
##  [919] 4.5 4.0 3.2 4.3 5.9 3.9 4.0 4.2 4.6 4.4 4.8 5.2 4.4 4.0 4.3 3.0 4.5 3.7
##  [937] 6.6 4.3 4.0 3.1 4.4 4.1 2.2 3.8 3.9 4.7 4.1 4.1 3.6 4.7 4.3 4.4 2.9 5.0
##  [955] 5.6 5.8 4.9 4.3 3.0 4.9 4.3 4.0 5.3 2.0 3.2 5.2 3.3 4.0 4.0 4.5 2.8 4.1
##  [973] 3.9 4.5 4.8 5.4 5.5 4.2 5.8 6.3 5.2 4.4 5.6 3.6 3.7 2.9 4.4 3.7 3.5 4.9
##  [991] 3.6 5.1 4.3 4.4 5.6 5.5 4.4 4.5 3.2 3.4
## 
## $func.thetastar
## [1] -0.0503
## 
## $jack.boot.val
##  [1]  0.45701220  0.28873239  0.28095238  0.12507205  0.01742857 -0.08275862
##  [7] -0.32822086 -0.33250689 -0.47325905 -0.49943343
## 
## $jack.boot.se
## [1] 0.9747624
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
##    [1] 5.5 3.9 4.9 5.2 4.8 4.9 4.2 6.8 3.7 4.8 4.9 3.4 4.8 4.7 4.1 4.7 4.6 6.1
##   [19] 5.7 4.3 6.4 4.1 5.1 5.4 3.5 4.7 5.1 4.1 4.4 3.3 4.1 3.6 5.5 4.0 2.4 4.9
##   [37] 3.5 3.9 5.4 4.8 4.2 5.0 4.7 5.3 3.9 3.9 3.5 4.8 4.1 2.7 5.0 5.2 6.6 5.7
##   [55] 4.2 5.2 4.5 6.6 4.2 6.6 3.4 4.8 5.5 4.4 3.9 3.3 6.8 4.8 4.6 4.4 4.6 4.9
##   [73] 4.7 4.7 4.3 5.4 4.2 3.7 3.2 5.2 3.6 4.9 3.2 4.6 3.1 3.6 4.7 3.4 5.0 5.3
##   [91] 3.4 5.5 5.0 5.3 3.6 5.8 4.7 5.0 2.9 4.6 3.7 5.4 4.9 5.2 4.7 2.6 5.7 4.7
##  [109] 4.0 5.5 3.6 6.1 4.7 5.0 2.7 4.4 4.9 3.7 5.5 4.4 4.0 4.7 5.3 3.2 6.0 4.5
##  [127] 3.2 3.9 3.3 4.7 4.7 5.7 3.5 3.5 5.7 4.9 3.6 5.0 4.9 3.5 5.9 4.6 4.9 4.9
##  [145] 6.0 4.2 5.1 3.6 5.2 4.8 3.4 4.9 4.2 3.3 5.7 6.3 4.1 4.4 5.3 3.9 3.4 6.0
##  [163] 3.6 3.6 4.3 3.5 4.7 5.3 4.9 3.8 3.8 4.7 3.1 3.9 4.6 5.1 3.6 2.6 4.8 5.3
##  [181] 4.9 2.6 4.8 5.9 3.9 4.0 3.5 3.8 3.1 5.0 2.5 5.3 3.1 3.2 4.8 5.8 4.0 4.5
##  [199] 6.4 4.2 4.5 4.0 5.2 5.2 5.1 4.9 4.1 5.4 5.6 4.9 5.0 4.3 4.2 4.9 4.3 5.4
##  [217] 4.3 4.5 3.3 5.0 3.8 2.1 4.7 5.0 4.3 4.8 4.2 4.5 3.2 4.0 4.2 3.5 4.1 5.2
##  [235] 4.6 6.0 5.1 6.4 4.0 4.4 4.1 4.2 4.0 3.4 4.7 3.7 3.9 5.3 6.2 4.9 7.5 5.5
##  [253] 5.4 3.8 4.0 3.6 2.7 3.8 4.1 3.6 4.6 3.9 3.0 4.5 4.8 4.8 5.3 4.2 4.3 2.8
##  [271] 4.0 4.9 5.5 6.4 6.3 3.3 4.5 4.6 4.6 4.2 3.7 5.7 3.9 5.2 2.1 2.5 5.9 5.9
##  [289] 4.8 4.2 3.7 4.1 4.0 5.4 3.8 3.0 4.0 5.1 3.9 4.6 5.2 3.9 5.5 4.7 3.9 3.6
##  [307] 4.2 4.7 4.0 5.6 5.9 4.3 2.8 3.3 3.7 4.6 3.8 4.5 5.2 4.3 6.0 4.1 4.6 6.0
##  [325] 3.9 5.0 5.1 4.8 3.5 5.0 4.3 4.3 3.5 5.1 4.1 4.7 3.4 4.1 3.3 5.7 4.5 5.7
##  [343] 4.4 3.7 3.1 4.2 3.5 5.2 5.7 4.8 4.5 4.4 4.3 4.0 6.0 3.9 4.1 4.9 5.6 3.5
##  [361] 4.4 5.0 3.1 4.9 5.9 4.8 2.8 5.6 5.5 6.5 4.6 3.9 3.5 4.2 4.4 3.7 6.0 4.7
##  [379] 4.2 4.1 6.3 5.6 4.3 5.9 3.9 5.1 4.6 4.6 5.1 4.5 3.9 5.5 3.8 5.1 3.2 4.4
##  [397] 4.4 4.7 4.3 3.5 5.3 5.2 3.7 4.7 5.5 4.6 4.3 5.1 2.8 3.7 3.6 4.9 4.8 5.1
##  [415] 3.8 4.8 3.2 5.0 5.4 4.0 4.2 6.3 5.0 5.2 3.3 4.9 3.4 5.0 4.0 4.0 3.9 4.7
##  [433] 5.1 4.1 3.9 5.8 5.5 6.0 5.5 5.0 5.1 5.4 3.1 1.3 4.4 4.9 3.0 3.7 3.0 3.7
##  [451] 5.6 3.4 4.8 3.4 5.3 4.7 4.5 4.8 3.8 4.9 6.5 3.8 5.5 4.3 5.0 4.3 4.4 5.0
##  [469] 6.2 2.9 5.3 5.5 6.0 3.4 4.0 5.0 5.0 4.1 5.2 4.6 4.6 3.8 5.4 4.1 4.1 4.0
##  [487] 3.5 4.1 3.4 5.8 5.5 4.3 4.5 5.5 3.8 5.8 3.7 5.7 3.7 3.3 5.8 5.2 4.1 4.2
##  [505] 3.9 4.2 5.1 3.7 3.4 5.2 4.6 3.4 4.8 5.5 5.6 5.3 6.4 2.6 4.5 4.5 5.8 5.5
##  [523] 3.6 4.8 6.7 4.8 5.2 4.3 3.0 4.3 3.6 3.3 3.9 5.4 4.7 5.1 3.5 5.5 3.9 4.7
##  [541] 3.8 5.6 4.9 2.8 3.8 6.0 5.1 5.9 3.8 4.5 4.7 4.0 6.4 4.5 4.8 3.0 4.0 6.0
##  [559] 5.1 2.7 4.1 3.5 4.9 4.4 4.0 4.5 4.1 2.5 4.5 3.1 2.9 3.8 3.3 4.4 5.4 4.6
##  [577] 4.8 4.4 4.7 4.1 3.6 6.5 3.5 4.1 6.7 5.9 4.5 5.7 6.1 4.5 4.8 5.5 3.9 4.1
##  [595] 4.5 4.3 5.4 4.3 5.7 4.3 6.2 6.3 4.8 6.1 2.6 3.7 2.8 4.3 4.3 1.7 3.6 5.0
##  [613] 4.9 4.3 5.5 4.3 5.5 3.7 4.9 4.7 6.0 5.0 4.9 4.8 3.3 5.7 3.8 3.9 5.7 5.4
##  [631] 4.4 6.2 5.0 4.5 4.8 5.8 4.5 4.9 5.2 3.7 3.2 5.5 3.7 4.0 6.0 4.4 4.3 4.8
##  [649] 4.9 3.5 3.8 5.6 4.0 3.7 5.1 3.2 4.2 3.2 5.9 3.4 3.5 4.1 3.8 4.9 4.7 4.5
##  [667] 4.4 3.3 4.6 4.5 5.4 5.1 4.5 4.2 4.0 4.0 4.6 5.8 2.2 4.3 6.1 4.1 4.4 3.5
##  [685] 5.1 4.3 4.3 4.1 4.7 4.6 3.9 3.9 5.2 5.7 3.7 3.5 4.5 5.3 4.6 3.1 3.5 4.7
##  [703] 6.1 5.4 5.1 4.5 5.2 4.3 5.8 4.1 5.0 3.2 3.9 4.8 4.7 4.1 5.4 5.1 5.6 4.0
##  [721] 5.9 4.2 4.9 4.4 5.7 5.3 3.8 5.2 4.6 4.2 4.3 3.5 4.5 5.7 5.3 2.7 5.2 4.6
##  [739] 4.5 3.9 3.6 5.2 4.7 4.6 3.2 4.8 4.9 5.4 5.4 4.0 4.8 5.9 5.4 5.5 2.7 4.6
##  [757] 4.1 3.1 4.7 4.9 4.6 3.4 4.9 4.7 3.7 5.3 4.5 2.9 3.5 4.1 5.5 4.0 4.2 4.8
##  [775] 4.3 4.8 3.7 3.6 4.0 4.3 3.2 3.6 3.7 3.6 5.9 4.1 4.6 5.2 3.7 5.2 3.6 5.1
##  [793] 5.1 5.8 4.7 4.6 4.4 3.1 5.3 5.0 3.9 5.2 4.3 5.9 2.6 6.2 5.4 5.5 3.4 5.1
##  [811] 4.3 4.9 5.5 3.6 4.3 4.3 4.7 4.1 4.5 3.3 3.7 5.4 7.0 5.4 5.4 4.0 3.7 5.5
##  [829] 4.1 2.7 5.3 4.6 3.2 6.3 4.9 5.2 4.2 3.7 6.3 5.0 4.6 4.2 4.1 5.2 3.4 4.7
##  [847] 3.3 3.9 6.1 4.8 3.9 4.9 4.2 3.4 5.1 2.6 4.5 5.0 5.7 4.4 4.5 4.8 4.8 4.9
##  [865] 5.1 4.0 5.1 4.1 4.2 5.5 3.5 4.7 4.1 3.6 4.3 3.3 5.1 4.0 5.5 5.2 4.0 4.9
##  [883] 5.5 4.3 4.4 3.2 4.3 3.3 3.4 5.6 6.4 3.0 6.6 4.2 4.9 4.2 3.3 3.6 3.4 4.2
##  [901] 6.3 5.0 3.9 5.0 5.6 3.7 5.0 5.0 5.5 5.5 3.1 5.7 3.4 3.6 4.5 4.4 5.1 4.0
##  [919] 4.5 5.7 4.5 4.0 3.7 4.7 3.9 3.9 4.4 4.0 6.0 3.1 5.7 4.6 4.2 4.7 4.8 5.3
##  [937] 2.8 3.6 4.2 5.7 5.6 4.2 4.2 5.6 4.2 4.4 2.9 3.0 4.2 3.3 4.1 5.2 4.0 3.3
##  [955] 5.1 3.1 5.3 6.0 3.7 4.6 4.6 4.6 3.3 3.9 5.5 4.3 5.9 3.4 3.0 2.4 5.6 5.3
##  [973] 4.0 4.1 4.7 4.4 4.6 6.4 5.3 3.1 4.0 2.7 4.0 4.0 2.4 4.8 3.8 4.2 5.3 4.2
##  [991] 4.5 4.2 4.4 6.6 5.0 4.8 3.6 4.6 5.1 4.7
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.5 5.5 5.3 5.3 5.1 5.0 4.9 4.7 4.6 4.5
## 
## $jack.boot.se
## [1] 1.032279
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
## [1] 0.5254647
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
##   2.507848   3.291666 
##  (1.055485) (1.533395)
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
## [1]  1.9216041  1.1734862  0.8551375 -0.1333007  1.0379693  1.4464812
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
##    [1]  0.1724851873  1.5148718294  0.4494139675  0.4274335877  0.3582663865
##    [6]  0.9197173680  0.5005712965  0.6529565179  0.6657371928  0.0140865527
##   [11]  0.9539932985  0.9055051090  0.4426493795  1.0201904817  0.2645349911
##   [16]  0.9498816921 -0.8469341272  0.1971075240  1.4737828369  0.5328164845
##   [21]  0.7793720242  0.8238266983  0.2859993284 -0.1160494303 -0.6164977999
##   [26] -0.1521716798 -0.1654611455 -0.1816271544  0.2108903885  0.3368995197
##   [31]  0.1509251225  0.6418486669 -0.6028038712 -0.1732672747  0.3725974626
##   [36]  1.0166461520  0.0769299929  0.4572351743  0.8647450805  0.4561582471
##   [41]  0.7210641934  0.6023279347 -0.3009171250 -0.6294330297  0.6123566889
##   [46] -0.4451354040  1.0436528443 -0.0315194792 -0.2998008377  0.0708124192
##   [51]  0.3087577962  0.9879596065  0.8258937023  0.4085103462  0.0135920105
##   [56]  0.0830065660  0.7887255054 -0.0144438984  0.3399970113  0.1653974286
##   [61]  1.3087975619  0.1276301023 -0.0845032547  0.5003888609  0.6445680619
##   [66]  0.2310389816  1.1898500480  1.2094828345  0.9033617854  1.1295949361
##   [71]  0.2637507206  0.2791380771  0.2131798328  1.0323690784 -0.6441865762
##   [76]  0.0242735664  0.4833092390  0.0432667193  0.0613391195  0.0517928233
##   [81]  0.9759060108  1.1444493376  0.3354827728  0.8673325115  0.9366335099
##   [86]  0.4070706810  0.1465141490  0.7599561060 -0.2298403168  0.5645097113
##   [91]  0.6070829027  0.1638150337 -0.1667363726  0.1285983310 -0.0469119042
##   [96] -0.2686443636  0.8009862704  0.1842796673  1.7296561022 -0.1863859355
##  [101]  1.2708394065  0.7707214808  0.6052013714  0.8199859225  0.6745922007
##  [106]  0.5453412598  0.7590446090  0.7653729020  0.6468858873  0.8064436192
##  [111]  1.3088282372  0.2845443812  0.5255281319  0.9891965331 -0.2507914394
##  [116]  0.2106287492 -0.1304183576  1.2922916072  0.4329244022  0.0435779118
##  [121]  0.3791759828  1.0378040999  0.4224623219 -0.2596598905  0.4706790543
##  [126]  0.8515031194  1.3283725172  0.3029971851  0.3137763782  0.9248441779
##  [131]  0.4957209876 -0.0161648071  0.7855926337  1.4604416220  0.1107643516
##  [136]  0.7828553882  0.1993725463  0.8922224567 -0.3889907288  0.5770955447
##  [141]  0.6969254781 -0.0938023625  0.1943014506  0.0049095103  0.8464053782
##  [146]  0.5279669870  0.6442473219  0.3403134339 -0.0463867236  0.5335327977
##  [151]  0.5622790236  0.4038782000  0.4708520833  0.5738735199  1.4130427015
##  [156]  0.1265762602 -0.3965006537  0.8156577345  0.8337553717  0.2514023928
##  [161]  0.0771481013  0.9831418731  0.5521827902  0.6751886880  0.4711958522
##  [166]  0.4009165386  0.5743742013 -0.0751550510  0.0373520606  0.6883913253
##  [171]  0.8741527531  1.1301772001  0.7782945428  0.6744760499  0.8777814158
##  [176] -0.0419951805  0.9913516736  0.3456519802  0.3643950426 -0.3569481747
##  [181]  0.2517892473 -0.0938023625  1.0517938829  0.2632757549  0.2415545463
##  [186]  0.2452239319  1.2413415934  1.1301772001 -0.3675230472 -0.0884015571
##  [191]  0.6296564104  0.2496771673  0.2127989417  1.0623239966  0.1358791499
##  [196]  0.9579574140  1.2758660636  0.2692287528  0.2996663518  0.5073309368
##  [201]  0.4296694741  0.5428611464  0.3087693495  0.4311663182  0.8899773350
##  [206]  0.2297093967 -0.3266663893  0.1818299559 -0.3863944891  0.0440898051
##  [211]  0.4968285438  0.8568532901  1.1655934452  1.2566218766  0.0063734099
##  [216] -0.2242076804  0.1380247035 -0.0173151201  0.2285204373  0.4890481199
##  [221]  0.0723873725  0.1429430878  0.4612338813 -0.1003279423  0.5987165307
##  [226]  1.4074100609  0.5434836894 -0.1030387399  0.1178633527 -0.5842833950
##  [231] -0.0247697103 -0.3561008464  0.5785474621 -0.3382780858  0.9889093417
##  [236]  0.5468736256  0.1367515095 -0.1884913777  0.5100543868  0.5486663122
##  [241] -0.1622374369 -0.0126001686  0.5954994792  1.0472426224  0.5825789046
##  [246]  0.8977636644 -0.3045264116  0.4102386876 -0.0856264634  0.4602955599
##  [251]  0.7496177061  0.7998910788  0.3520330118  0.4138836734  0.8921635557
##  [256] -0.1943305475  0.4194893038 -0.0149258163  0.2360962476 -0.0334719419
##  [261] -0.0238528597 -0.3394442221  1.0280798055  0.7062556910  0.3527825870
##  [266]  0.2485884411  0.4383073511 -0.1168757905  0.4682741525  0.2175921455
##  [271]  0.1858944660  0.1750279133  0.4274335877 -0.3563883999  0.5539233975
##  [276]  0.7358001494  0.0496303126  0.7533782042  0.9357838560  0.2152906835
##  [281]  0.9145245208  0.3221365621  0.7356140327  0.0002247954 -0.1279567115
##  [286]  0.0810371215  1.6189294100  0.6420492182 -0.3148637318  0.9218587917
##  [291]  0.4551195874 -0.1105272643  0.3206264238 -0.2110531603 -0.1374442257
##  [296]  1.4817685636 -0.2882431602  0.1668594140  0.1409764011  0.1013245918
##  [301]  0.7397194779  0.7793443104  0.5442772731  1.2342854619  1.2700662276
##  [306]  0.4268405724  0.6692923350 -0.0373732117  1.1029416809  0.1489625379
##  [311] -0.0812589665  1.0622090195  0.8795262770  0.1006542840  0.2640557826
##  [316]  0.7204517719  0.4224623219  0.4913777502  1.0141295322  0.6358150238
##  [321]  0.3125786177  0.9690226727  0.2259999464  0.0836893722 -0.0037843556
##  [326]  0.3259211419  0.8695242484  0.2327859447 -0.0324100319  0.1684332667
##  [331]  0.7410804978  0.4899633970  0.4476244685  0.4891640727  0.3075323323
##  [336] -0.2202583423  0.4174344240  0.9108528734  0.7205374531  0.6408124474
##  [341]  0.6116445286 -0.7462413575  0.2153960802  0.6228577921  0.9842890773
##  [346]  1.1493110263  0.5093524031  0.3246926019 -0.0222050279  0.7980903922
##  [351] -0.1446492289  0.4968985045  0.8135074908  0.6027625915  0.1746592551
##  [356]  0.4450907734  0.0428042515 -0.3992936208  0.7069037987  0.3895094600
##  [361] -0.0663359798  0.3415881502  0.5573647164  0.5940040557  0.4240388133
##  [366]  0.1875320540  0.4744201280  0.3120789972  0.5422687615  0.7021542247
##  [371]  0.0468300417  0.8509286405 -0.1680968878  0.2457099096  1.1475151935
##  [376]  0.8740987091  1.1724337454  0.3928683312  0.2083979399 -0.1948455983
##  [381]  0.6503860453  0.6524175458  0.5720864443  0.3204566528  0.6856342735
##  [386] -0.1045282224  0.0609945907  0.1174121009  1.1864986624  0.8675732956
##  [391]  0.4340363930  0.0277484187 -0.1237454552  0.9437018736 -0.3603395507
##  [396] -0.1365935684  0.8040946503  0.1557176354  0.2969037874  0.7312419188
##  [401]  1.0661250324  0.1678838966 -0.2198874058  0.9571634017  0.0158935049
##  [406]  1.2824131438  0.4690903092  0.3916545314  0.2505304084  0.7526140532
##  [411]  0.5828979062 -0.2935329281  0.7215536398  0.7532867075  0.5237231982
##  [416] -0.4914741735  0.5991069688  0.6030077869  0.8997826983  1.0094521791
##  [421]  1.2525458302  0.9069065455 -0.0155887279  0.7082342945  0.5451445829
##  [426]  1.1134750594  0.1699973447  1.5132858158  0.7064875378 -0.0748625035
##  [431]  0.1714244838  0.4554897586  1.2680040307  0.6520327348 -0.3045073268
##  [436]  0.3742266983 -0.0995035513  0.7352589650 -0.4615557799 -0.5374187444
##  [441] -0.0997987740  0.2212311221  0.9250360378  0.4166067970  0.1872574869
##  [446]  0.2131599187  0.3159486537 -0.4778841414 -0.2103739625  0.0559927169
##  [451]  0.9283869462 -0.0970867906  0.3056896710  0.3741709745  0.4469586626
##  [456]  0.6729708014  0.8424424765  0.1535411945  1.0272993089  0.3100897059
##  [461]  0.7305661760 -0.4053569217  0.1341048976  0.0150417295  0.1220028157
##  [466]  0.4551793019  0.8837103624  0.8258374947  0.7072782102 -0.1732672747
##  [471]  0.5019066662  0.3183024355  0.2794290090  0.3832461950  0.5072341444
##  [476] -0.1508599104 -0.2516661623  1.3000866733 -0.0209266030  0.4511047280
##  [481]  0.7816462950  0.7139492155  1.0610195682  0.7188594542  0.0463725997
##  [486] -0.1025455804  1.0312635303  0.3795158638  0.3418041742  0.3779919331
##  [491]  0.6638194090  0.5477400302 -0.1131589468 -0.1396870699  0.4827768289
##  [496]  0.1570498147  0.7311014830  1.3470188952  0.3494066613  0.9528134532
##  [501]  0.6313285642  0.2200583806 -0.7458306439  0.5832740346  1.1103308344
##  [506]  0.8101453970  0.6404170458  0.8293102890  0.7396326144  0.1837128558
##  [511] -0.1676922926  0.7240055653  0.6461612733 -0.1530556150  0.7815801742
##  [516] -0.7124850238  0.0366892146 -0.5602734550  0.4720719101  0.7065704121
##  [521]  0.1338167457  0.2634598267  1.0869213124  1.0375785726  0.8454144941
##  [526]  0.6473550149  0.1281930114  1.1215152858 -0.4994701752  1.3297053480
##  [531]  0.1253172882  0.2358467230  0.2330173725 -0.1147528704  0.3048725765
##  [536]  1.0338971395  0.2433975874  0.5090078503 -0.0902763695  0.2004979746
##  [541]  0.7248126099  0.4824166467  0.2578463005 -0.0182489011  0.8642761912
##  [546]  1.0793225291  0.2575740352 -0.1386616254  0.2498495028  0.3357929265
##  [551]  0.9498539174  0.5243076989  0.1795688970  0.5435131937  0.2354772878
##  [556]  0.0963557263  0.3108145793  0.6246985379  0.0822641272 -0.5438758107
##  [561]  0.8262713374  0.0889634193  0.6125157202  0.6428072279 -0.4223470131
##  [566]  0.6370218556  0.0674858239  0.8487897153  0.7506484581  0.2186024357
##  [571]  0.1935051873  0.4342307363  0.8656871482  0.9836164963  0.8698487274
##  [576]  0.3491972901  0.3884786955  0.5476216573  0.1440584509 -1.1349696988
##  [581]  0.3425394373  0.0402832156  0.9670666624  0.1195427752  0.1833421239
##  [586]  0.2865319898 -0.5680940188  0.5075924556  0.8902888979  0.6410008186
##  [591]  0.6661851465  0.9672688759  0.1406629078 -0.3173194345  1.1608915064
##  [596]  1.4944085729  0.8712181574 -0.3554897074  0.2772805996  1.2027577635
##  [601]  0.4323063424  0.7627026869  0.7294696330  0.1807301317  0.8513721020
##  [606]  0.7042238884  0.6263790103  1.0640248514 -0.0601542194  0.9371482492
##  [611] -0.2495581125  0.5206902401  0.3814126831  0.9069918904  0.4545625661
##  [616]  0.8088300820  0.3514986632  0.6790063595  0.2888189604  0.6686346288
##  [621]  0.2244658135  1.0046761683  0.3950366159 -0.1798836518  0.0378107698
##  [626]  0.2173668772  0.7046453568  0.4518430732  0.9498765124  0.8977993316
##  [631]  0.4410576960  0.1972327959  0.9896163800  0.7030311495  0.7283002436
##  [636]  0.4729009311  0.0505430942 -0.6280341862  1.2265637359  0.5850979916
##  [641]  1.0194245368  0.4597521802 -0.3331387165  0.0104959779  0.7085145766
##  [646] -0.0653098223  0.5579177904  0.8640587777  0.1691351876  0.0118039959
##  [651] -0.2418613943  1.2481054324  0.9876111818  0.3636951811  0.4013637994
##  [656] -0.1496464865  0.1522404102  0.0033144138  0.5118042690  0.9160008434
##  [661]  0.3590016332  1.2508818305  0.8774905948  0.8306429604  0.3494066613
##  [666]  0.6470723700  0.4192307369 -1.2364732658 -0.4609893010 -0.0316589537
##  [671]  0.3864317355  0.0995088384  0.0942184908  0.5418659893  0.9313214541
##  [676]  0.7269653575  0.7649697292  0.3297851182  0.4987467314  0.8792309179
##  [681]  0.2495791224  0.0460480539  0.3935616769  0.1604385758  1.2446112146
##  [686]  1.3982044184  0.1885301499  0.2527630344  0.5426909626 -0.2595733578
##  [691]  0.5470928034  0.2534850945  0.3297860452 -0.2972194399  0.5207616099
##  [696]  0.7925978496  1.2152526572 -0.2097897609  0.3493802794  0.1366123515
##  [701]  0.6856802253  0.6569095194  0.8241238775 -0.2753431469  0.0153981169
##  [706]  0.3216662332  0.2337086156  0.3339662458  0.8316343500  0.8091135197
##  [711]  0.4191958711  1.1904441089  0.4558720907  0.9537116276 -0.3173194345
##  [716]  1.3226505380  0.1888075873  0.0613825757  0.0779937057  0.0468802003
##  [721]  0.5878189734  0.9267142904  0.0733057713  0.7162782714  0.4627509796
##  [726]  0.5683921628 -0.0001527970  0.9732217832  0.4347825076  0.2956250885
##  [731]  0.1295548058  0.3638061121  0.3699646715 -0.0802350185  0.2516088710
##  [736]  0.1227755521  0.2971228438  0.8550108982 -0.2838491931  0.5913702005
##  [741]  1.5523047635  0.1973934713  1.0045346021  0.8499558936  0.4703521974
##  [746]  1.2798523976  1.3113349518  0.4563320110  0.4435276274  0.6985197364
##  [751]  0.4837238513 -0.2864208532 -0.7047045646  0.5686376828  0.5303391438
##  [756]  0.0878971314  0.1753396599  0.6631870779  0.6384021568  0.3991677700
##  [761]  0.2366674270  0.9314664295  0.7343948711  0.1719865892  0.4551793019
##  [766]  0.3701021383  1.3154077629  0.2906371009  0.6224930471  0.5485280540
##  [771]  0.9923247597 -0.2805473417 -0.0450979948  0.0583790816  1.4207225203
##  [776]  1.2152693901  0.0347563658  0.7918404654  0.4984072768  0.3559823040
##  [781]  0.1432782561 -0.4719327122  0.8879210511  0.5853979435  0.7629987504
##  [786]  0.4601063182  0.7610188392  1.1423836838  1.0010961369  1.5274403792
##  [791]  0.4629741833 -0.0284510509  0.2960200950  0.9456837873  0.2316042874
##  [796]  1.3209268665  0.7556474290  0.5863692321  1.0865296774 -0.0147537996
##  [801]  0.0048223304  0.5003929933  0.3442303390  0.0804379385  0.3602201686
##  [806]  0.2687971141 -0.0006205369 -0.5608132192  0.7269653575  0.1310165375
##  [811]  0.3015458323  0.7207517219  0.1776464304  0.2795806748  0.8067081899
##  [816]  0.1303915862  0.4853685461  0.8119391964 -0.0183370887  0.1590545694
##  [821]  0.3742545906 -0.0592579959  0.8393238601  1.1551496179 -0.0024175836
##  [826]  0.6813019584 -0.0431733897  0.2093818220 -0.3517636029  0.6295193239
##  [831]  0.4123857428  0.1033498697  0.5592446264 -0.0221449046  0.1517006132
##  [836]  0.5623293716  0.1965679370  0.1167537157  1.7022265647  0.2814421685
##  [841]  0.7990207193  0.3253109361 -0.1859251869  1.2117261970  1.1758770254
##  [846]  0.6605838081 -0.0245471156  1.3670361313  0.3145725581 -0.0036782127
##  [851]  0.3502282897  1.0633020036  1.3589440705  0.2770124907  0.4900634227
##  [856]  0.3586620970  1.4950387214  0.7615790995 -0.2123231603  0.6229212777
##  [861]  0.0432082151 -0.8265821510  0.0995088384  0.2483048358  0.6550089938
##  [866] -0.4833653650  0.7568266018 -0.5723495808 -0.0128481601  0.6635656132
##  [871]  0.2727121456  0.3938550893  0.7177705976  0.2561682750  0.5065380936
##  [876]  1.2393335481  0.5519502018  0.9226399504  0.2213879665  0.3576126566
##  [881]  0.6210498252  0.2613363782 -0.1084442107 -0.1977596281 -0.2365172209
##  [886]  1.5359013350  1.7257836080  0.5676048779 -0.0057596120  0.2850348950
##  [891]  0.2251745924  0.6597432558  1.1939747201  0.5052012465 -0.6908150343
##  [896]  0.2480399328  0.2693382014  0.3663391855 -0.1872574633  0.4311588482
##  [901] -0.1414572090  0.8580488861  0.7011513466  0.2186228126  0.3068614943
##  [906]  1.1010476035  0.6184300219  1.1078646874  0.1251884729  0.9513098919
##  [911]  0.3655869227  1.3708905299  0.7601443294  0.0882929835  0.1443844312
##  [916]  0.6936703491  1.0395730685  1.0418681350  0.8429227712  1.6127640997
##  [921]  0.3685853896  0.2479255535  0.7604074831  0.8707161094  0.8435945430
##  [926]  0.3299018047  0.5216598788  0.4623182590  0.0022847586  0.7816823190
##  [931]  0.7091391051  0.5998538170  0.0802077838  0.2113589823  0.3907022530
##  [936]  0.0857536060  0.7327266101 -0.1273701813  0.5933454673  0.8126939612
##  [941] -0.3772445647  0.0894231616  0.5423100757  0.1763503957  0.9564024045
##  [946]  0.6028774061  1.1124167645  0.7498568814 -0.0764298855  0.5310015103
##  [951]  0.7946522855  0.2741888670  0.5277277838  0.5769458688  0.2297093967
##  [956]  0.5370742449 -0.4992613389  0.3981118947  0.2653084612  0.1074265772
##  [961]  0.8294879623  0.3881040715  0.6888486105  1.5952625098  0.0306195581
##  [966]  0.6487346050  0.1021410851  0.4668981165  0.5569525006  0.5582557511
##  [971]  0.3523255159  0.2794676356  0.4322930736  1.9582003689  0.7292927603
##  [976]  0.1878219538  0.6197640161  0.4731715823 -0.2046057329  0.2357743331
##  [981]  0.1058480539  0.5282244605  0.4508534751  1.0683798339  0.5657162185
##  [986]  0.3663391855  0.1632630902  0.6089889470  1.0880339768  0.3708223882
##  [991]  1.0550550645  0.4047869838 -0.4818171866  0.7953566453 -0.3543304128
##  [996]  0.2894680420  0.5128078985  0.4844840556 -0.5588091900  0.1985815518
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
##   0.7618783   0.4521457 
##  (0.1429810) (0.1011010)
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
## [1] -0.3596893  0.4158531  0.6754132  1.4081010  0.4029788 -0.3029577
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
## [1] 0.0205
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8904845
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
## t1*      4.5 0.06136136   0.9020843
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 1 2 3 6 8 9 
## 2 1 1 1 3 1 1
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
## [1] -0.093
```

```r
se.boot
```

```
## [1] 0.8986203
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

