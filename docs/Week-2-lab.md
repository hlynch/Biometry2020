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
## 1 2 5 6 7 9 
## 2 1 2 1 1 3
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
## [1] -0.0168
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
## [1] 2.693739
```

```r
UL.boot
```

```
## [1] 6.272661
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.8   6.3
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
##    [1] 3.6 4.9 5.0 4.2 4.8 5.3 4.7 5.2 4.2 4.0 4.9 5.4 6.1 5.7 5.5 5.2 4.4 3.3
##   [19] 3.2 3.3 5.1 3.5 5.7 4.0 3.7 5.8 3.6 5.3 4.3 5.4 5.9 2.0 3.9 4.5 4.2 3.3
##   [37] 4.3 3.7 6.1 4.9 4.4 3.4 4.5 3.4 4.3 4.1 5.1 5.0 4.4 3.4 5.3 2.7 4.7 4.0
##   [55] 4.1 3.6 5.1 2.8 5.7 4.8 5.5 5.1 3.2 3.8 5.4 4.8 5.1 3.9 4.4 4.4 3.5 4.5
##   [73] 3.6 4.6 5.1 4.7 4.9 3.4 5.6 4.6 4.6 4.9 4.1 5.6 5.2 4.7 3.6 5.3 5.2 5.1
##   [91] 4.3 6.4 3.7 4.1 3.4 4.6 4.2 4.9 3.6 3.0 3.3 4.8 4.2 4.3 4.6 3.3 3.4 4.0
##  [109] 4.2 3.2 5.0 5.2 5.0 3.2 3.7 4.0 2.9 4.5 4.8 2.9 3.8 4.0 4.2 3.4 3.3 4.0
##  [127] 6.0 4.0 4.1 3.1 2.8 4.4 5.4 5.2 3.2 4.8 4.1 5.6 3.3 3.4 4.4 3.6 5.0 4.4
##  [145] 5.9 6.0 5.1 4.3 4.3 3.3 3.9 3.7 4.6 4.0 4.5 4.3 6.3 3.6 4.1 3.9 4.7 6.0
##  [163] 3.6 4.7 3.1 5.5 4.0 4.9 5.8 5.4 4.4 6.2 5.8 2.9 4.9 5.1 4.5 4.0 4.4 4.5
##  [181] 3.9 4.9 3.0 5.0 6.0 6.4 4.6 3.4 5.1 2.2 4.9 5.2 3.3 4.7 4.2 4.8 2.7 3.6
##  [199] 4.4 5.0 5.7 3.5 3.8 4.5 5.5 6.2 4.3 3.2 4.4 4.5 4.2 3.8 4.7 4.2 3.9 4.5
##  [217] 3.3 4.4 4.9 4.5 3.4 3.9 4.8 2.9 4.5 6.4 4.2 4.9 5.5 3.7 4.8 4.5 3.3 4.3
##  [235] 4.6 5.1 5.4 4.6 4.1 4.8 5.1 4.5 3.4 4.7 4.1 5.1 4.7 5.5 4.6 4.9 5.1 2.9
##  [253] 2.9 4.8 4.3 5.6 5.4 4.0 4.1 4.7 4.5 5.5 5.5 3.7 4.0 2.8 4.5 3.3 5.3 5.5
##  [271] 4.3 3.7 5.8 5.1 3.6 2.5 4.4 6.3 4.7 5.4 3.7 4.4 3.8 3.5 6.4 3.5 4.9 2.7
##  [289] 5.4 3.1 3.9 4.4 5.5 3.9 6.1 5.4 4.5 5.2 5.8 4.3 4.2 4.1 4.2 4.4 3.6 4.9
##  [307] 5.3 5.4 4.9 4.2 5.5 5.1 5.0 5.7 4.3 5.4 5.4 2.4 5.4 3.7 3.2 2.6 2.9 4.0
##  [325] 3.4 4.5 6.0 3.8 5.5 4.0 4.2 4.4 6.6 4.5 5.4 3.5 4.7 4.3 4.2 5.2 4.3 5.0
##  [343] 4.7 4.7 5.6 4.4 4.6 4.8 3.0 4.3 6.1 5.8 3.8 4.9 5.5 4.3 2.8 2.3 5.9 4.3
##  [361] 3.5 5.1 4.3 5.6 4.1 4.3 3.4 4.7 5.9 4.6 4.1 5.0 4.4 3.9 3.9 4.1 4.5 5.7
##  [379] 4.3 5.6 6.2 4.7 3.0 4.7 4.8 4.4 4.1 4.1 5.2 4.8 7.0 5.3 3.2 3.5 5.2 4.8
##  [397] 4.2 4.5 4.2 3.7 5.6 3.7 3.9 6.1 4.4 4.5 4.0 6.3 5.2 4.8 4.6 6.3 5.2 3.6
##  [415] 4.4 5.3 5.9 4.1 3.5 4.5 6.6 5.8 5.7 4.9 3.2 5.4 2.7 3.3 3.6 5.4 5.6 3.5
##  [433] 2.6 5.4 4.8 5.6 4.1 4.0 4.6 4.5 4.1 4.0 5.4 3.1 4.6 4.4 3.3 5.5 4.7 5.7
##  [451] 5.1 3.7 5.0 4.1 6.7 2.3 3.7 4.7 4.0 6.1 5.7 4.0 2.6 3.9 1.9 3.1 3.2 4.8
##  [469] 4.2 4.6 4.3 4.2 3.9 4.7 4.0 4.7 4.3 5.1 4.7 5.1 4.8 4.1 6.3 4.1 5.9 4.9
##  [487] 6.1 4.7 4.3 5.1 3.6 3.6 5.3 3.2 6.2 4.9 4.4 4.5 3.7 3.9 4.7 5.0 3.7 4.9
##  [505] 5.2 4.0 5.5 5.5 5.3 4.2 3.7 4.8 3.4 3.7 3.4 5.7 4.5 4.5 5.3 4.9 3.9 3.5
##  [523] 3.3 3.9 5.3 3.9 4.6 4.5 5.9 6.5 3.3 4.9 4.0 5.5 3.4 3.7 3.0 5.3 4.2 4.5
##  [541] 4.9 5.1 4.9 4.5 3.0 4.7 4.5 4.3 5.8 4.0 5.1 5.5 3.7 4.8 5.6 4.8 4.6 4.1
##  [559] 4.6 5.2 4.2 6.2 4.6 3.7 5.3 5.1 3.5 4.0 4.7 5.1 4.7 4.6 4.1 6.5 4.0 3.9
##  [577] 5.0 5.3 4.3 2.9 4.1 4.7 5.1 4.2 4.1 3.3 4.4 4.2 3.7 4.3 5.8 4.5 5.3 4.5
##  [595] 5.7 5.8 3.8 4.7 4.0 5.3 3.6 4.5 5.0 4.8 4.5 4.3 6.3 5.1 4.4 4.4 4.6 4.4
##  [613] 5.9 5.8 3.7 4.6 3.9 4.8 3.7 5.4 3.6 3.7 5.1 4.8 3.8 3.0 4.6 7.3 5.3 4.4
##  [631] 3.8 7.1 4.9 4.4 1.9 4.0 3.6 3.1 3.3 4.6 5.1 4.9 4.1 5.5 5.9 5.1 3.7 3.0
##  [649] 5.8 4.2 4.1 6.0 3.1 5.0 4.6 3.3 4.0 3.8 5.0 5.2 5.0 5.4 3.8 4.1 4.2 4.5
##  [667] 3.8 5.2 3.2 4.6 6.5 5.1 3.4 5.5 3.2 3.3 4.3 4.5 4.3 4.0 5.0 4.3 4.6 5.4
##  [685] 4.5 6.8 3.8 4.6 5.6 3.5 4.2 4.8 5.3 4.1 2.7 5.1 3.8 4.8 6.2 4.5 6.0 5.4
##  [703] 4.6 4.6 5.5 4.2 5.1 1.8 4.1 4.9 3.1 3.7 2.9 2.9 4.6 3.6 5.2 5.3 3.7 5.1
##  [721] 3.5 4.9 4.0 4.2 7.1 5.2 4.3 6.4 3.9 3.7 4.9 5.1 4.0 4.9 4.4 5.1 2.8 5.2
##  [739] 3.5 4.0 3.9 5.3 3.9 5.0 4.5 4.5 2.8 4.9 4.3 5.2 4.2 3.8 5.0 3.6 5.0 3.7
##  [757] 4.7 4.6 4.2 4.4 4.8 5.6 4.3 5.5 3.6 4.1 3.4 3.1 5.4 6.2 5.1 5.3 4.8 3.9
##  [775] 6.5 4.3 5.1 4.6 3.8 5.0 3.7 2.4 3.2 4.3 3.5 5.6 4.7 4.0 5.4 3.9 4.0 5.2
##  [793] 3.1 3.0 4.7 1.8 6.2 5.6 5.0 4.1 3.7 2.4 4.2 4.2 5.3 3.9 5.1 3.8 4.4 4.7
##  [811] 4.6 3.4 4.4 3.0 5.6 4.9 3.5 5.4 4.1 4.0 5.8 5.3 6.0 4.8 4.7 5.0 4.0 3.4
##  [829] 4.6 3.3 3.8 4.1 5.5 2.7 4.8 4.7 3.8 5.4 4.4 4.5 4.9 4.3 3.7 3.8 5.4 5.6
##  [847] 3.9 4.8 2.5 4.0 4.1 3.8 5.5 3.6 3.3 5.0 4.2 5.1 4.7 4.2 3.1 5.2 5.1 5.0
##  [865] 2.8 4.3 4.5 5.0 4.8 3.9 5.0 2.6 7.1 3.7 3.3 4.7 5.5 5.1 3.5 3.7 3.6 5.0
##  [883] 4.3 4.7 6.5 3.7 4.4 5.3 3.3 4.4 4.3 4.3 3.4 3.8 4.8 2.4 3.2 5.0 5.5 4.2
##  [901] 4.5 5.9 5.4 4.3 3.5 3.7 4.7 5.5 4.1 5.4 4.1 5.5 4.1 4.2 4.5 3.5 4.9 3.7
##  [919] 4.6 5.8 4.8 5.0 5.1 4.5 3.5 4.0 3.9 5.6 3.3 4.7 5.2 6.0 4.2 3.1 5.0 3.8
##  [937] 2.7 4.7 6.9 4.9 2.3 5.1 3.9 3.3 4.3 4.4 4.7 2.6 2.6 4.8 4.3 4.8 6.2 3.6
##  [955] 5.5 4.8 4.0 4.4 5.8 3.1 6.4 3.9 5.4 4.6 4.5 3.5 3.8 5.8 4.6 5.3 4.7 4.2
##  [973] 4.8 4.3 4.9 3.8 5.6 4.4 6.0 5.2 3.9 4.3 4.8 5.8 3.8 2.4 5.0 4.0 4.0 4.3
##  [991] 4.0 3.9 3.5 5.3 5.4 4.4 3.8 3.8 5.2 4.0
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
##   2.7   6.3
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
##    [1] 5.0 3.4 4.2 4.7 5.1 5.3 4.3 3.1 5.6 5.0 5.6 5.0 5.7 4.1 5.4 3.1 4.8 5.4
##   [19] 4.6 2.7 4.0 3.6 3.7 6.1 5.3 5.0 4.8 3.5 5.7 5.4 5.3 4.6 6.2 4.8 3.4 3.3
##   [37] 4.5 3.8 3.7 4.8 3.7 4.5 2.8 4.0 4.4 4.6 4.6 6.3 4.1 4.6 3.8 4.4 3.0 4.9
##   [55] 6.0 3.9 3.7 4.1 3.9 6.0 3.6 4.5 5.7 5.8 5.2 3.7 5.6 2.1 4.2 4.3 4.4 5.5
##   [73] 5.1 5.5 6.0 5.1 5.1 4.6 3.5 4.1 5.0 4.1 5.5 4.2 3.2 3.3 4.0 3.7 3.4 4.0
##   [91] 5.0 4.1 4.4 3.9 4.6 4.4 4.3 4.2 2.9 4.9 5.0 4.5 4.7 4.5 4.1 4.4 3.3 5.5
##  [109] 4.3 3.4 4.2 3.6 4.6 5.2 3.6 4.4 3.1 4.5 5.2 4.6 4.2 3.6 3.6 5.4 4.3 3.4
##  [127] 4.2 4.8 4.2 6.5 5.1 5.2 4.9 4.4 4.7 4.3 5.1 5.3 4.4 4.5 4.8 4.5 5.3 3.9
##  [145] 3.6 4.8 4.8 2.7 2.4 5.3 3.8 3.6 3.5 3.8 2.7 5.1 4.1 5.0 5.2 3.7 4.9 5.5
##  [163] 3.7 4.4 3.6 4.3 2.7 5.0 4.8 5.5 5.3 5.0 4.9 3.8 4.7 5.3 3.4 4.4 3.4 3.8
##  [181] 2.5 4.8 3.2 4.2 2.8 4.9 4.3 5.1 6.1 4.2 3.8 5.3 5.6 6.6 4.7 4.9 4.3 3.8
##  [199] 4.5 4.4 3.6 4.4 3.6 3.9 4.7 5.9 4.0 5.4 4.8 4.1 5.4 5.9 3.5 4.6 3.9 4.0
##  [217] 4.3 3.9 6.2 6.3 3.8 3.2 5.6 2.3 4.2 4.7 5.0 5.3 4.1 4.4 4.3 4.2 5.1 4.7
##  [235] 5.2 5.3 5.8 4.5 5.6 4.8 5.1 5.4 4.3 5.6 3.3 5.3 4.8 2.9 3.4 6.8 5.2 4.3
##  [253] 5.2 4.9 3.7 3.4 4.8 3.3 3.5 5.4 6.3 3.7 2.9 5.5 4.7 4.3 5.3 4.9 4.6 4.5
##  [271] 4.1 6.2 4.0 6.0 5.7 5.9 3.8 4.6 5.2 3.9 4.4 5.2 4.2 5.6 5.4 4.3 5.8 3.7
##  [289] 4.1 4.7 5.7 5.2 3.0 5.0 5.1 5.3 4.2 3.9 5.2 5.9 5.6 5.3 4.4 4.5 5.9 3.5
##  [307] 4.8 4.9 5.2 4.4 4.0 5.2 6.0 4.7 3.8 4.8 3.2 5.8 4.8 3.3 3.3 3.9 4.9 6.8
##  [325] 4.0 4.5 5.1 3.7 5.2 4.1 3.4 6.5 2.9 2.8 4.5 5.6 3.7 4.2 3.6 4.3 4.1 3.6
##  [343] 4.3 3.6 3.9 5.2 5.1 4.4 4.9 4.6 3.8 5.2 5.7 4.0 4.5 4.2 3.1 6.0 4.8 3.9
##  [361] 4.4 4.5 6.0 6.5 2.8 3.8 5.9 3.3 5.1 5.4 5.6 6.0 3.4 5.2 5.7 4.5 4.1 3.7
##  [379] 5.7 5.3 4.4 3.7 5.8 5.0 3.8 3.4 5.4 4.0 4.2 4.5 4.6 4.0 5.1 5.1 4.7 3.7
##  [397] 5.6 5.5 6.4 4.2 5.2 5.7 4.2 5.4 5.4 3.5 3.3 5.6 5.5 4.4 4.3 5.6 5.2 6.7
##  [415] 4.3 4.9 5.5 4.6 5.7 3.7 5.1 5.1 4.2 4.6 3.0 3.6 4.9 3.2 2.7 2.0 4.9 4.1
##  [433] 4.9 3.9 3.9 6.4 5.3 3.7 4.3 5.4 3.8 4.8 3.5 3.8 5.0 6.3 4.8 2.6 3.2 5.0
##  [451] 3.9 3.7 5.2 3.6 4.1 2.6 3.8 4.6 5.8 6.3 3.5 2.9 4.7 6.2 4.3 4.2 4.5 4.7
##  [469] 3.1 4.2 3.1 3.4 4.7 4.4 5.0 2.5 4.0 3.8 4.4 4.6 5.0 5.9 4.3 5.2 5.0 4.5
##  [487] 4.5 4.9 5.1 5.1 2.5 3.9 3.4 5.5 5.1 3.6 5.1 2.8 3.0 5.6 4.9 4.0 3.7 3.0
##  [505] 4.6 4.5 5.5 4.0 4.1 3.5 5.8 4.0 5.1 3.7 3.6 6.3 4.1 4.8 5.5 5.3 5.0 5.1
##  [523] 5.4 4.6 5.4 4.9 4.5 4.9 4.2 4.9 3.6 4.4 4.6 3.5 4.4 6.3 4.4 5.2 4.7 5.6
##  [541] 4.7 3.2 5.2 5.2 2.6 4.9 5.7 5.5 4.0 3.7 3.5 5.0 6.1 4.9 4.2 2.5 4.4 5.0
##  [559] 4.5 3.5 3.4 4.4 3.5 6.9 4.6 2.6 3.5 5.9 4.8 3.7 5.6 5.8 5.6 4.9 3.4 5.3
##  [577] 3.4 2.2 4.3 6.6 4.2 5.0 3.8 4.1 5.1 4.8 5.6 4.6 5.4 4.9 4.2 3.5 4.9 3.9
##  [595] 4.6 5.2 4.5 4.2 6.2 4.0 4.4 3.8 5.8 3.8 4.6 5.0 5.3 5.7 5.1 5.4 3.3 4.6
##  [613] 4.1 3.8 4.7 4.3 3.3 4.1 3.5 5.3 4.0 5.2 4.2 3.2 3.5 3.7 5.8 5.3 3.8 3.6
##  [631] 4.7 3.7 5.3 4.9 6.0 4.6 6.3 5.7 5.3 3.5 5.1 5.2 4.8 4.5 4.7 4.3 5.4 4.3
##  [649] 3.6 4.6 4.3 3.9 3.3 3.6 4.4 5.0 2.9 4.7 3.8 5.5 3.2 4.2 3.7 3.8 5.6 3.7
##  [667] 5.3 3.1 4.2 4.5 4.8 5.0 4.7 3.6 2.8 3.6 4.2 4.7 5.3 4.8 5.7 5.7 6.8 4.7
##  [685] 5.5 5.4 3.2 5.2 4.6 3.2 6.2 5.4 4.9 5.5 3.2 5.1 6.3 4.8 3.2 5.7 4.1 5.4
##  [703] 5.5 4.3 4.3 5.1 3.4 4.4 6.2 4.7 6.5 5.9 3.3 6.2 5.2 3.1 5.2 5.8 4.0 5.7
##  [721] 4.5 4.9 5.8 3.1 4.3 4.9 6.4 4.4 4.5 4.9 5.7 4.5 4.5 4.7 5.7 4.5 5.0 4.9
##  [739] 4.2 3.5 5.1 5.3 5.9 4.6 6.2 4.2 5.4 4.2 4.6 4.1 4.6 3.1 3.6 4.7 4.4 3.5
##  [757] 3.7 4.0 3.7 4.1 3.9 5.9 3.4 2.8 5.6 3.9 4.8 4.4 3.4 4.0 2.9 3.1 3.5 5.1
##  [775] 3.1 4.0 5.6 3.5 4.5 5.5 4.9 5.3 3.6 5.2 3.8 2.7 3.4 5.5 3.4 4.8 3.1 4.9
##  [793] 4.9 6.3 4.3 3.2 5.3 4.8 6.0 4.9 5.0 5.1 5.1 5.0 4.5 6.2 2.8 3.3 3.9 4.3
##  [811] 3.9 4.7 4.6 4.0 3.9 6.2 5.6 3.1 5.9 5.4 4.9 4.2 4.4 3.4 4.0 4.4 3.6 5.3
##  [829] 4.2 5.9 4.8 3.9 5.3 4.8 5.0 3.5 4.4 4.5 6.0 4.8 5.6 3.9 4.0 3.4 5.9 4.0
##  [847] 5.5 4.4 2.8 3.1 4.1 3.9 4.5 4.0 5.1 3.3 2.9 5.6 5.6 4.4 5.2 4.3 4.3 3.5
##  [865] 3.6 4.6 3.3 6.5 4.8 5.7 5.3 3.2 4.6 4.4 6.4 3.7 6.0 3.9 3.2 5.6 3.3 5.4
##  [883] 3.1 4.3 3.8 5.6 2.9 3.6 4.4 4.8 4.1 3.1 3.8 4.0 3.7 4.7 4.7 4.9 4.6 3.6
##  [901] 3.4 4.2 4.7 3.9 6.0 5.2 3.4 5.3 4.2 5.3 4.3 4.1 3.2 5.7 6.0 6.0 4.2 5.6
##  [919] 4.2 5.6 4.6 4.2 4.1 6.1 3.7 3.2 5.5 5.2 5.5 3.8 4.8 3.8 3.6 4.1 5.3 4.5
##  [937] 3.7 4.1 5.0 4.3 5.6 2.2 3.8 3.2 5.2 4.6 3.6 5.4 3.3 5.1 5.9 5.7 4.0 3.3
##  [955] 5.7 5.6 5.4 5.3 4.6 5.8 3.4 4.1 4.4 5.2 4.9 4.4 5.0 4.0 3.8 6.7 4.7 4.4
##  [973] 4.6 5.0 4.7 5.2 3.9 5.0 6.4 4.5 3.8 4.0 4.0 4.1 2.5 5.7 6.2 3.8 5.4 3.4
##  [991] 5.3 4.9 5.5 3.2 4.2 4.3 4.4 6.0 5.5 3.8
## 
## $func.thetastar
## [1] 0.0301
## 
## $jack.boot.val
##  [1]  0.575522388  0.462162162  0.322560976  0.182142857  0.008994709
##  [6]  0.010215054 -0.100835655 -0.290000000 -0.388414634 -0.466961652
## 
## $jack.boot.se
## [1] 1.010901
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
##    [1] 5.4 3.1 4.6 4.2 6.0 3.9 4.3 3.3 5.6 4.0 2.8 5.3 3.7 4.5 6.2 4.1 4.3 6.0
##   [19] 4.3 5.2 5.1 4.8 5.1 5.1 1.6 5.9 4.8 5.5 4.4 4.6 3.4 4.2 3.1 4.0 4.9 4.2
##   [37] 4.6 3.6 4.8 4.5 5.3 3.3 4.2 5.3 5.1 5.1 5.6 4.0 3.5 4.1 4.7 3.3 4.3 4.2
##   [55] 4.8 7.1 3.3 4.8 4.4 6.0 5.6 6.0 6.0 3.1 4.2 4.8 3.6 4.6 4.9 3.4 3.9 4.4
##   [73] 3.3 3.7 4.5 4.2 3.8 6.0 4.5 4.2 5.5 3.8 3.7 5.8 4.4 4.2 4.9 4.9 5.0 5.0
##   [91] 6.0 4.7 5.4 6.1 5.9 4.8 3.4 4.2 3.8 4.3 5.1 5.5 6.4 4.1 4.2 5.3 5.3 5.6
##  [109] 3.9 5.6 3.5 4.3 5.2 4.0 5.3 3.2 4.4 4.5 4.4 4.7 3.8 4.5 2.4 4.1 4.8 5.0
##  [127] 5.1 3.8 5.1 4.8 5.8 4.3 5.7 5.3 4.8 3.8 4.0 4.7 3.7 6.5 3.7 4.8 5.0 5.3
##  [145] 3.5 4.1 5.8 4.6 4.3 4.2 4.5 3.3 3.8 5.0 3.3 4.5 3.6 4.4 3.1 3.3 4.6 5.3
##  [163] 6.3 2.4 4.8 4.4 5.1 4.5 5.1 5.0 4.7 4.5 5.3 4.4 5.3 4.2 5.0 5.0 4.0 1.7
##  [181] 5.0 4.8 5.5 3.8 5.2 6.5 4.2 4.8 3.9 5.1 3.2 4.0 4.4 6.1 4.5 4.0 3.3 4.5
##  [199] 5.8 4.6 4.2 2.8 4.9 3.6 3.9 3.9 4.2 5.1 4.5 5.5 4.8 4.2 4.1 3.0 4.9 3.9
##  [217] 4.1 4.1 5.1 4.9 3.3 5.1 4.9 3.5 4.9 3.7 5.2 3.4 4.6 3.6 2.3 3.4 5.3 4.5
##  [235] 6.6 5.1 5.3 5.7 5.0 5.8 3.8 3.4 4.5 3.0 3.7 3.0 3.8 3.7 5.3 4.3 4.9 4.4
##  [253] 4.8 3.6 6.9 5.6 3.3 5.1 3.6 3.4 3.7 4.0 5.0 5.4 4.8 5.4 5.1 3.3 6.1 4.3
##  [271] 5.0 4.7 3.9 4.3 4.8 4.3 4.7 3.3 5.4 3.5 5.1 4.1 5.4 4.6 4.0 3.0 4.7 5.8
##  [289] 4.6 4.1 2.6 3.6 5.6 3.8 4.5 5.8 4.5 3.8 3.6 4.3 4.1 5.7 5.6 3.4 4.2 4.4
##  [307] 4.3 4.4 5.3 5.0 5.4 4.0 5.4 4.3 3.5 4.3 4.2 3.6 4.7 4.0 4.6 3.2 4.2 4.7
##  [325] 3.7 3.0 5.0 3.4 5.7 4.8 5.0 5.4 4.5 5.6 5.3 4.7 2.8 4.3 4.0 4.8 3.3 4.3
##  [343] 3.8 3.5 4.3 3.5 4.9 5.4 5.6 4.1 4.0 5.6 5.5 3.5 4.1 4.0 5.9 4.7 5.8 3.8
##  [361] 3.0 4.9 4.2 4.8 4.8 3.4 5.6 5.5 5.6 5.1 6.4 3.8 5.4 3.6 4.7 5.5 3.4 5.0
##  [379] 3.9 5.0 3.2 4.0 2.5 6.2 6.0 3.4 5.6 5.3 4.6 4.4 4.8 5.0 4.9 5.1 4.3 3.0
##  [397] 5.4 3.2 5.9 2.6 3.3 4.1 5.9 5.0 3.4 4.6 5.2 3.3 4.7 5.1 4.0 5.6 5.5 4.2
##  [415] 4.6 2.9 4.2 3.4 5.1 4.3 5.2 4.6 4.2 4.9 5.5 3.4 4.7 5.0 5.1 5.1 5.1 5.2
##  [433] 5.0 4.3 3.3 5.8 5.4 4.3 4.3 3.9 4.4 5.5 7.1 4.1 5.2 6.9 4.6 4.1 5.0 5.6
##  [451] 4.6 4.6 4.3 5.2 3.0 5.7 5.9 6.2 4.9 2.6 2.8 4.8 5.0 3.9 4.5 5.8 3.1 3.8
##  [469] 3.5 3.8 3.6 2.4 3.1 4.8 6.0 3.5 4.9 5.1 5.9 3.4 3.9 3.2 3.2 4.9 3.8 5.0
##  [487] 4.2 3.0 4.7 3.5 4.0 4.0 4.1 5.5 4.3 6.1 4.5 3.9 3.0 5.9 3.5 5.4 3.9 3.9
##  [505] 5.4 4.2 4.5 5.2 6.0 4.1 4.2 5.2 4.9 5.2 4.8 5.6 4.1 5.0 4.4 3.4 5.9 3.6
##  [523] 3.8 5.5 4.4 2.9 4.5 4.6 5.2 6.3 3.7 5.1 3.4 4.3 5.6 4.8 5.5 4.5 3.6 4.5
##  [541] 6.8 4.4 2.2 5.6 3.9 5.3 5.4 6.0 3.3 3.1 4.7 3.6 3.5 4.0 5.0 4.4 3.5 4.7
##  [559] 4.1 4.4 2.5 5.1 6.1 2.7 5.0 7.0 4.1 4.0 4.1 4.2 2.8 3.8 4.3 3.1 5.8 4.0
##  [577] 4.5 2.7 4.1 4.9 3.1 6.2 3.5 4.1 5.1 3.5 2.7 3.5 4.7 4.3 4.0 5.0 5.1 5.0
##  [595] 4.6 3.9 4.6 4.6 3.8 4.5 3.7 3.8 4.3 4.6 4.1 4.7 5.5 5.8 3.6 4.7 4.4 4.0
##  [613] 5.3 3.8 5.1 4.1 3.8 6.4 4.1 4.2 7.1 3.6 5.1 4.7 3.2 4.2 4.2 5.4 4.0 5.1
##  [631] 4.4 4.3 5.6 4.9 3.6 4.7 4.5 3.4 5.8 4.8 5.0 5.6 5.5 3.6 4.4 3.9 3.5 5.5
##  [649] 4.1 3.6 4.0 4.3 3.7 3.6 5.2 4.3 4.4 5.5 6.6 3.5 3.7 3.4 6.2 4.1 4.5 4.6
##  [667] 3.0 3.5 4.5 3.7 4.5 5.2 5.3 3.3 5.9 4.7 4.5 4.0 3.6 3.3 4.6 3.3 3.2 5.3
##  [685] 4.2 3.9 4.6 3.6 4.3 4.0 3.5 4.5 5.0 3.5 4.2 4.0 2.9 4.0 5.2 4.5 3.0 4.3
##  [703] 6.3 2.5 3.0 4.3 5.2 4.8 4.9 3.8 5.4 5.7 4.5 5.1 4.2 3.9 3.7 4.0 3.8 3.6
##  [721] 3.9 5.5 5.4 5.9 5.0 6.0 2.9 4.0 4.2 5.6 5.0 4.5 4.3 4.2 5.0 3.9 4.9 4.5
##  [739] 5.1 6.0 6.4 5.3 3.2 3.8 4.2 4.0 4.4 6.6 5.5 5.0 3.8 5.5 5.3 4.4 3.8 5.7
##  [757] 4.1 5.2 5.5 4.6 4.4 5.2 5.6 3.6 4.1 5.4 4.3 4.7 7.2 4.2 4.7 3.3 5.9 4.6
##  [775] 3.7 4.2 6.4 4.5 4.4 2.3 4.4 3.9 4.9 4.6 4.3 5.6 5.8 3.3 4.6 4.6 5.8 4.8
##  [793] 5.2 4.4 4.4 4.9 3.7 3.1 4.9 4.3 5.9 4.8 5.0 5.3 3.8 3.6 4.6 5.8 4.8 6.1
##  [811] 5.4 5.9 4.6 3.5 5.1 3.3 4.0 5.2 5.2 4.1 3.9 5.4 4.5 4.5 4.3 4.7 5.1 5.1
##  [829] 4.2 4.8 4.1 5.7 4.9 3.5 3.4 4.5 5.2 4.7 6.6 3.9 5.4 4.9 5.0 4.3 3.3 5.2
##  [847] 3.4 3.8 1.9 4.2 3.6 3.1 2.9 3.7 6.0 5.0 3.9 4.1 5.2 5.6 3.5 5.3 3.9 3.9
##  [865] 6.0 3.3 4.8 5.1 3.7 3.9 3.4 5.3 4.9 4.5 5.0 4.0 5.7 4.7 3.8 5.9 4.9 5.0
##  [883] 4.2 3.3 4.3 3.8 2.6 5.0 4.6 6.5 5.7 4.6 5.3 5.8 3.2 4.1 5.1 5.2 5.0 4.5
##  [901] 4.5 6.1 5.8 5.9 4.8 5.2 3.8 4.8 4.1 4.6 3.8 5.7 5.5 4.6 4.7 5.4 6.6 3.9
##  [919] 5.4 3.3 4.1 3.1 4.4 5.0 5.1 4.6 4.2 4.6 4.0 4.9 4.1 4.1 6.1 3.9 5.2 4.5
##  [937] 4.8 5.3 5.2 5.9 5.5 4.4 4.0 5.3 3.6 5.7 5.1 4.8 3.6 3.8 4.4 4.3 4.5 3.2
##  [955] 4.1 5.9 3.4 3.5 5.5 5.9 3.3 3.4 5.3 5.2 3.0 5.5 3.6 4.6 5.1 3.9 5.7 3.4
##  [973] 5.1 3.9 4.6 4.3 2.8 5.7 4.7 5.6 3.8 3.4 5.8 5.2 3.6 2.9 6.5 3.2 5.3 5.1
##  [991] 5.0 4.6 4.4 6.2 5.2 2.6 5.6 4.3 4.4 4.6
## 
## $func.thetastar
## 72% 
## 5.1 
## 
## $jack.boot.val
##  [1] 5.500 5.468 5.300 5.400 5.100 5.000 4.800 4.800 4.604 4.500
## 
## $jack.boot.se
## [1] 1.036439
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
## [1] 0.20955
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
##   21.70207   38.67205 
##  ( 9.63180) (17.36300)
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
## [1]  0.5617233  0.4855255  0.8723190 -0.2510972 -0.7146691  0.5216789
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
##    [1] -0.584310600  0.496896534  0.360743890 -0.409449965 -0.614063986
##    [6] -0.487025126  0.102215426  0.682931648  0.583178560  0.168765647
##   [11] -0.038600117 -0.047610576 -0.233337509  1.168159526  0.729571753
##   [16]  0.392609575  0.231470575 -0.551957585 -0.431932245 -0.490279467
##   [21]  0.280976965  0.168174355  0.958698258  0.340126168  0.628050093
##   [26] -0.751154124 -0.042393035  0.500929916 -0.113512851  0.034554907
##   [31] -0.233494712 -0.012537256  0.420371488  0.459471783  0.176069245
##   [36]  0.429807383  0.394528980  0.960410092 -0.108918101  1.234054312
##   [41]  0.329619120  0.210749132 -1.141579977 -0.331045309 -0.070419662
##   [46]  0.425411498  0.160378828  0.708278338  0.702435659  0.694060357
##   [51]  0.379772311  0.082639477  0.505414161 -0.455046603  0.575463806
##   [56]  0.089446392  0.014051567  0.465907827 -0.378037720 -0.225311149
##   [61]  0.459987102 -0.166755292  0.566558030  0.847762940  0.632081769
##   [66] -0.156597265  0.012212813  0.908575968  0.322973416 -0.714711137
##   [71]  0.492834864 -0.534706567  0.086827339  0.892490874  0.129222463
##   [76] -0.091413763  0.906138305  0.670666844  0.495910461  0.351149571
##   [81]  0.095509744 -0.316152314 -0.132196791  0.116214239  1.188918398
##   [86] -0.312451771  0.553155659  0.520690111  0.151564450  1.095602906
##   [91]  0.480083133 -0.073186104 -0.073680952 -0.638840770  1.122738812
##   [96] -0.250161398 -0.034794401  0.159882281  0.480727538 -0.024462878
##  [101]  0.014546051 -0.601132674  0.673192875 -0.423677962  0.931878224
##  [106]  0.061900005 -0.127196668  0.629438510  1.000501634  0.308899927
##  [111]  0.056596946  0.386252865  0.008538002  0.179818088 -0.087634306
##  [116]  0.283190616 -0.789666765  0.257409231  0.111048494  0.180462367
##  [121]  0.119239327  0.398798141 -0.196816128  0.329312649  0.070613377
##  [126] -1.088659712  0.231283022  0.116415863  0.454745691  0.236218138
##  [131]  0.587989845  0.314426339 -0.548385525 -0.800044350  0.400047165
##  [136] -0.609688219  0.294859529  0.087417097  0.220204242 -0.084025003
##  [141]  0.403715591  0.501245604 -0.120718185  0.597709548  0.247123109
##  [146]  0.956451546  0.251847936  0.062661240  0.014336361  0.679328729
##  [151]  0.182880049  1.293270738 -0.010928548  0.567869665  0.491021072
##  [156] -0.199479862  0.447141556  0.379331018  0.506607130 -0.069035492
##  [161]  0.187468175 -0.567386832  0.154574812 -0.118255045  0.599263595
##  [166]  0.164692735 -0.388753658  0.541346079  0.735175175 -0.515535211
##  [171]  0.020220446  0.079694637  0.728567540  0.553460990  0.317365027
##  [176]  0.006642652  0.221730303  0.862080479 -0.184647966 -0.304079692
##  [181] -0.364612945  0.534427667  0.637733677  0.394543885  0.131634805
##  [186] -0.916052343 -0.172127298 -0.852765851  0.123729109  0.785946618
##  [191] -0.031232585 -0.529727216  0.220815420  0.209550037  0.967788526
##  [196]  0.185784465  0.326733414 -0.012374122 -0.075491841 -0.200952629
##  [201] -0.379499032  0.488295304  0.168588372  0.220348913  0.607752967
##  [206] -0.206046328  0.210430136  0.170397999  0.295256650  0.785345124
##  [211]  0.525663171  0.589585377 -0.819764910 -0.414334051  0.258681280
##  [216]  0.072029629  0.288966699  0.307293737 -0.255800931  0.248842681
##  [221] -0.366192160  0.372480084  0.268027883 -0.529469463 -0.032213554
##  [226] -0.356470470  0.163248252 -0.325158481  1.216960659 -0.154885364
##  [231]  0.297355539  0.789729383  0.351223015  0.167185744  0.675507053
##  [236]  0.583198802  0.368808943 -0.224078347  1.132255246 -0.116509648
##  [241]  0.014051567  0.489860945 -0.273945924  0.510347213  0.192216102
##  [246]  0.285668834  0.131660225  0.916628634  0.976484235  1.462398579
##  [251]  0.640220741  0.881074454 -0.015158034  0.334124762  0.702312597
##  [256]  0.396738594 -0.094016775  0.478729094  0.126746265  0.111876136
##  [261]  1.295366191  0.619177857 -0.274606480 -0.457086517 -0.967794865
##  [266]  0.307146138  0.155217164 -0.065832302  0.820410203  1.083398417
##  [271]  0.370941497  0.391106823 -0.025287096 -0.623769911 -0.710788444
##  [276]  0.910509338  0.478685651  0.572908073  0.501361725 -0.201092148
##  [281] -0.756141250  0.404992301 -0.599236581  0.277290136  0.148250127
##  [286] -0.752173391  0.794147368  0.282788444 -0.091195253  0.622353297
##  [291]  0.574942130  0.273473754  0.613914700 -0.667299543 -0.153706170
##  [296] -0.316834573  0.623968767 -0.097609129  0.719378529  0.379894823
##  [301] -0.077117585 -0.097366320 -0.452669976 -0.056883846  0.206798655
##  [306]  0.140904658  0.222744726  0.286859885  0.151365326  0.559898096
##  [311]  0.579241813  0.399199756  0.747174070  0.101213630  0.623797164
##  [316] -0.389656994 -0.136468760  0.030998335  0.824248710 -0.145612705
##  [321]  0.461653886  0.201647124  0.172170927  0.052815604  0.502300599
##  [326]  0.317840007 -0.228588932 -0.440749255  0.390339871  0.622087227
##  [331] -0.482813084  0.102934327  0.262786614 -0.380175206  0.729227060
##  [336]  0.100727436  0.232975441  0.902574378 -0.145258582 -0.059984860
##  [341] -0.273898009  0.230005085  0.461959970 -0.564897508  0.195157221
##  [346]  0.360050458  0.226321409  0.316680952  0.640220741 -0.375856546
##  [351]  0.229838019 -0.396464222  0.242766945  0.500571636  1.237682573
##  [356]  0.552365492  0.707605396  0.383454690  2.295770801  0.225935367
##  [361]  0.747610013  0.537632298 -0.007585074  0.206480235  0.100059307
##  [366]  0.077359225  0.208342158  0.010727805  0.207744811  0.157689212
##  [371] -0.559802013  0.515856499 -0.319284881 -0.029777485 -0.436815346
##  [376]  0.544705014  0.629518272  1.178252212  0.335021670 -0.303178490
##  [381]  0.708782701 -0.714791224 -0.150621040  0.125297407  0.772149300
##  [386] -0.120389919  0.185863188  0.091499774 -0.117489092  0.732360825
##  [391] -0.227294078  0.053586966  0.115130896  0.147746291  0.044723230
##  [396] -0.279235289 -0.010290203 -0.133254168  0.250416211 -0.460814543
##  [401]  0.075642643 -0.031646680 -0.360083699  0.584557963  0.531146310
##  [406]  0.912558422 -0.446502141  0.415461459 -0.378150696 -0.312537604
##  [411]  0.431368454 -0.070826929  0.591226571  0.521240828  0.129599534
##  [416]  0.475252409  0.143426586 -0.403122355  0.425784415  0.919117126
##  [421]  0.385850599 -0.073963186 -0.176346570  0.005534130  0.170536817
##  [426]  0.052890292 -0.748276063  0.098285225  0.955568673  0.403410378
##  [431]  0.092352136  1.029951221 -0.234155535  0.219477579  0.061500329
##  [436]  0.661598896  0.021834749  0.409723726  0.530897010  0.623797164
##  [441] -0.233335712 -0.154792053  0.513574558  0.232458919 -0.633798210
##  [446]  0.812649529  1.099714652  0.195626816  1.343341475 -0.005161731
##  [451]  0.177757541  0.474198118 -0.409124948  0.894227670  0.269044220
##  [456]  0.268825219 -0.006367433  0.196728526 -0.658065537 -0.815201253
##  [461]  0.379337850  0.247000236  0.206215772  0.027738876 -0.265944818
##  [466] -0.396464222 -0.557278731 -0.037527818  1.310776776 -0.422118425
##  [471]  0.569815680  0.761087386  0.198276289 -0.488075807  0.075646017
##  [476]  0.610336205  0.481016159  0.034820323  0.338077010  0.547924228
##  [481]  0.108726511 -0.342771066  0.294206995  0.140315711  0.117110714
##  [486]  0.194054416  0.244669778  0.269047572  0.183736555 -0.311071706
##  [491]  0.799961945 -0.296444212  0.792101716  0.016724721  0.559898096
##  [496]  0.307293737  0.253645084  0.501125769  0.364282760  0.047648182
##  [501]  0.044531960  0.750406935  0.415843530  0.596683962  0.340566474
##  [506] -0.525847144  0.810146115  0.040163601  0.996694228  0.482200254
##  [511]  0.625439702  0.187468175  0.032318483 -0.394653928 -0.229477342
##  [516]  0.264401977 -0.468847037 -0.010828083  0.436044516  0.485191823
##  [521]  0.165395932  0.591007834  0.083608793  0.247924733  0.413157160
##  [526]  0.074474088 -0.552771474 -0.061495030  0.905629700  0.719165291
##  [531]  0.569815680  0.281855047 -0.037238118  0.802047365 -0.196653039
##  [536]  0.436588983  0.303858367  0.478392472 -0.024760076  0.062752628
##  [541]  0.469420219  0.577870945 -0.340169418  0.135089263 -0.293456517
##  [546] -0.222421775  0.417143358  0.012561062  0.115769189  0.068190864
##  [551] -0.075855235 -0.516899832  0.323347892  0.234836389  1.166401061
##  [556] -0.101362321 -0.176268803 -0.041275596 -0.163411672 -0.761392840
##  [561] -0.991370396  0.356769626 -0.066310366 -0.522194650  0.612788177
##  [566]  0.183353026 -0.546843854  0.383568096  0.239233880 -0.072957541
##  [571]  0.073646000  0.049688570  0.421718625  0.288237951  0.248409414
##  [576]  1.259803351  0.373115584 -0.702517860 -0.216477658  1.023755713
##  [581]  0.462934393  0.395805791  0.853887442  0.024231784  0.362007623
##  [586] -0.067557788 -0.066354578 -0.258866092 -0.205270130  0.784391286
##  [591]  1.303500408  0.231283022 -0.595044989  0.058835962  0.277942940
##  [596] -0.651813801 -0.560174624  0.407295252 -0.410782971  0.317205775
##  [601]  0.476691096 -0.043892508  0.217078215  0.495218430 -0.229013762
##  [606]  1.071364315  0.421718625  0.017196806  0.013634270  0.531144903
##  [611]  0.230759317  0.880699581 -0.026868550  0.593150933  0.428384266
##  [616] -0.375494497 -0.143137534 -0.287645253  0.615621169  0.520733044
##  [621] -0.159063819 -0.522416999 -0.038986422  0.527485705  0.011449617
##  [626] -0.518773699  0.170983549 -0.137358913  0.479627129  1.153119078
##  [631] -0.253408686 -0.160021748 -0.039296317 -0.299509016 -0.688068823
##  [636]  0.377381506 -0.556388340  0.319891972 -0.645159116 -0.165766036
##  [641]  0.489389839 -0.270528735  0.006034065  0.045757135  0.380287844
##  [646] -0.341451427 -0.003854327 -0.483144907 -0.028913873 -0.343924954
##  [651] -0.236590525  0.279706562  0.584718744 -0.686853366  0.008629052
##  [656] -0.001423016 -0.563674575  1.241694645 -0.537183181  0.126665207
##  [661]  0.157675631 -0.037142860  0.799166622  0.303413590 -0.377927098
##  [666] -0.008060420 -0.006038365  0.068044265  0.492477012 -0.360704668
##  [671] -0.171424187 -0.111013333  0.170348513 -0.031554970 -0.552656876
##  [676] -0.149182361  0.038232695 -0.739794677 -0.338679337  0.266336095
##  [681] -0.989832920  0.133610945 -0.449388303  0.842223358  0.001590991
##  [686] -0.592141287 -0.015823674  0.674864060  0.239429873 -0.253063688
##  [691]  0.098409806  0.200846168 -0.285741643  1.552687053  0.044531960
##  [696]  0.093148682  0.586752945  0.952405059  1.344961813  0.048637908
##  [701]  0.568648908  0.343601522  0.828503957  0.312418501  0.932598350
##  [706]  0.168326053 -0.070312498 -0.099383773 -0.450459191 -0.150049330
##  [711]  0.234544751  0.869191629  0.322351855  0.399121915 -0.050237393
##  [716] -0.601768700  0.615239262  0.248280633  0.076051106  0.714734200
##  [721]  0.504540682  0.468997014  0.389476750  0.555420866  0.080872712
##  [726] -0.610593802  0.666130615  0.017916060  0.058256992  0.031322244
##  [731]  1.048021032 -0.297598960  0.046783883 -0.459150055 -0.607933115
##  [736]  1.143059157  0.312846373  0.539462237  0.629990701 -0.919211194
##  [741]  0.555201882  0.209216403 -0.326340251  1.557731306  0.144522419
##  [746] -0.431305636  0.399121915  0.289465628 -0.082346632 -0.221315915
##  [751] -0.128092094  1.020481548  0.490589016 -0.060648263  0.643124622
##  [756] -0.207942443 -0.261514553  1.498760160 -0.390542283  0.526979477
##  [761]  1.004319506 -0.072644850  0.362217354 -0.031646680  0.429255993
##  [766]  1.418645914  0.361318804  0.188353173  0.961778356 -0.609620260
##  [771]  0.375357657  0.094748223  0.273676218 -0.510599706  0.443031409
##  [776] -0.530400837  0.273676218  1.425026496 -0.392569363 -0.130091166
##  [781]  0.199319669  0.846080168 -0.216151363  0.039255112 -0.170190497
##  [786]  0.190803051  0.187553231  0.220152695 -0.023599780  0.574069155
##  [791]  0.646414770  1.613099924  0.845654110  0.374858769 -0.217374280
##  [796] -0.195297772 -0.425238010 -0.558362132  0.299825280  0.560316426
##  [801]  0.325874231  0.125743076 -0.152339962  0.145912806  0.944082849
##  [806] -0.207735166  0.820910604  0.670078840 -0.381952046  0.035960381
##  [811]  0.343881602  0.446092578 -0.264302416  1.019463865 -1.034653334
##  [816]  0.229200577 -0.130140820  0.196859791 -0.176932587 -0.088402620
##  [821]  0.995868456  0.544705014  0.159285047 -0.160021748 -0.500741023
##  [826] -0.430664378  1.059335485  0.541482916  0.708231512 -0.225316271
##  [831]  1.179459963  0.430286442 -0.223154820 -0.585943877  0.286815796
##  [836]  0.096979054  0.947482621  0.184257757 -0.404594932 -0.268348152
##  [841]  0.271579491  0.032069814 -0.139870284  0.584752357 -0.029720186
##  [846] -0.190461117  0.855631623 -0.130374604 -0.580818533 -0.305811467
##  [851] -0.018864267 -0.465983682  0.262892350 -0.249258048  0.379060415
##  [856] -0.636426236  0.375860718 -0.746616641  0.030679800 -0.098743112
##  [861]  0.199698752  0.070613377  0.275278232 -0.199664573  0.801319563
##  [866] -0.519428943  0.029095166 -0.010178079 -0.285294146 -0.060338247
##  [871] -0.299934837  0.632824051  0.807587958 -0.734862093 -0.062885886
##  [876] -0.396369963 -0.121987779  0.086159574 -0.739191512  0.465525004
##  [881]  0.547924228 -0.639304727  0.977273348 -0.106414022  0.892490874
##  [886]  0.284135572  0.280279984  0.078341234  0.001815884  0.116265851
##  [891]  0.189932864  0.430689914  0.467597245  0.220863042 -0.098601450
##  [896]  0.297097582 -0.621308098 -0.739191512  0.777876422 -0.170407052
##  [901]  0.159090649  0.343040953  1.079978737 -0.435152841 -0.907679258
##  [906]  0.273283490  0.829065510  0.586710446  0.520593505 -0.046986534
##  [911]  0.897778694 -0.366533687 -0.601487115 -0.111975352  0.493740079
##  [916] -0.047194747  0.221711644  1.142604361  0.084190692  0.963076163
##  [921]  0.900472221 -0.164304440  0.775921000  0.016015640 -0.119565208
##  [926]  0.353390363  1.101948369  0.469464163  1.766180266  0.977757288
##  [931] -0.534040640  0.804766660  0.266192508  0.458712645  0.006236778
##  [936]  1.156491940 -0.420233078  0.434695022  0.181419907 -0.214247642
##  [941]  0.137241414  0.069004694  0.150943077  0.821827470  0.378084779
##  [946] -0.311825732  0.437769531 -0.081834009  0.340547572 -0.586242036
##  [951] -0.258180502  0.249790754  0.033438439  0.330152574 -0.391765519
##  [956]  0.160239970 -0.016463789 -0.362931893  0.433322967 -0.341326708
##  [961] -0.242204546  0.531610617  0.115279778  0.433240540  0.743386505
##  [966]  0.703761986  1.089143896  0.287648200 -0.424089448  0.300944793
##  [971]  0.692772060 -0.735639651  0.535308122  0.961747594  0.607634541
##  [976]  0.281094073 -0.259969905  0.941221554  0.195528396  0.378900263
##  [981]  1.093171392  0.871986029  0.144009296  0.490531250  0.151555378
##  [986]  1.104600776  0.311178273  0.821622274  0.582177536 -0.441494841
##  [991]  0.064191225  0.376203518  0.897993395  0.591575775 -0.304599184
##  [996]  0.085241289  0.236725787  0.556333396 -0.234532513 -0.248026970
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

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

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
##   0.56114990   0.12020121 
##  (0.03801096) (0.02686921)
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
## [1]  0.2075966 -0.1054708  0.6885978  0.5295741 -0.4246298 -0.2727695
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
## [1] -0.0232
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9158715
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
##     original       bias    std. error
## t1*      4.5 0.0002002002   0.9047304
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 3 4 5 7 9 
## 1 1 3 2 1 2
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
## [1] -0.0481
```

```r
se.boot
```

```
## [1] 0.9175174
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

