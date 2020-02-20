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

**Example #1**: Use permutation methods to test the null hypothesis that the treatment does not increase survival time (in other word: $H_{0}$: No difference in survival between the treated and control groups):

Survival.treated=$\{94,197,16,38,99,141,23 \}$

Survival.control=$\{52,104,146,10,51,30,40,27,46 \}$

(Is this a one-tailed or a two-tailed test?)

Make sure that you understand what is being done here, as this example is very closely related to the problem set.


**Example #2**: Using the same data, provide a 95% confidence interval for the difference in mean survival days based on 1000 bootstrap samples

Note that these two approaches are very closely related. Do you see why either approach can be used to test the null hypothesis? (What is the null hypothesis here?)

Now we will do one slightly more complicated example from Phillip Good's book "Permutation tests: A practical guide to resampling methods and testing hypotheses":

Holmes and Williams (1954) studied tonsil size in children to verify a possible association with the virus \textit{S. pyrogenes}. Test for an association between \textit{S. pyrogenes} status and tonsil size. (Note that you will need to come up with a reasonable test statistic.)

<div class="figure" style="text-align: center">
<img src="Table2categories.png" alt="Data on tonsil size and S. pyrogenes status. Source: Good (1994)" width="40%" />
<p class="caption">(\#fig:unnamed-chunk-1)Data on tonsil size and S. pyrogenes status. Source: Good (1994)</p>
</div>

Now lets consider the full dataset, where tonsil size is divided into three categories. How would we do the test now? What is the new test statistic? (There are many options.) What 'labels' do you permute?

<div class="figure" style="text-align: center">
<img src="Table3categories.png" alt="Fill dataset on tonsil size and S. pyrogenes status. Source: Good (1994)" width="50%" />
<p class="caption">(\#fig:unnamed-chunk-2)Fill dataset on tonsil size and S. pyrogenes status. Source: Good (1994)</p>
</div>

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
## 2 5 6 7 8 
## 2 1 2 3 2
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

<img src="Week-2-lab_files/figure-html/unnamed-chunk-6-1.png" width="672" />

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
## [1] -0.0015
```

```r
hist(xmeans,breaks=30,col="pink")
abline(v=mean(x),lwd=5,col="black")
abline(v=mean(xmeans),lwd=2,col="yellow")
```

<img src="Week-2-lab_files/figure-html/unnamed-chunk-7-1.png" width="672" />

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
## [1] 2.661334
```

```r
UL.boot
```

```
## [1] 6.335666
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
##    [1] 3.0 5.0 5.8 4.5 4.8 5.5 5.0 4.9 5.9 2.5 2.9 4.7 5.2 5.2 5.0 3.7 4.6 4.4
##   [19] 4.2 5.2 5.4 4.9 3.8 3.7 3.7 5.9 4.2 6.3 3.8 3.0 3.4 3.9 5.3 4.2 5.1 5.0
##   [37] 5.2 4.6 3.2 3.2 3.7 4.8 4.3 5.7 4.6 5.1 5.2 4.4 4.7 5.1 5.0 5.1 4.5 3.3
##   [55] 4.6 4.9 6.0 4.5 5.9 3.9 6.3 4.9 4.0 4.6 4.7 5.6 2.6 4.4 6.8 5.1 5.6 3.7
##   [73] 4.8 3.6 2.6 5.2 3.8 4.5 3.9 4.7 4.4 5.5 5.7 4.3 5.1 6.2 5.1 5.0 3.4 4.2
##   [91] 3.8 4.9 5.6 4.2 4.7 4.5 4.7 4.6 2.4 4.3 4.4 4.7 4.1 5.4 3.8 3.2 4.9 3.7
##  [109] 3.8 4.8 3.8 5.4 5.2 3.5 4.4 4.3 4.2 3.9 5.1 4.5 3.6 3.7 4.5 3.1 3.7 3.4
##  [127] 4.7 4.3 3.2 5.6 3.2 3.6 3.9 4.6 5.8 5.0 5.5 5.6 4.3 5.5 4.5 4.0 5.8 3.8
##  [145] 4.6 4.9 4.3 4.2 3.8 4.0 3.9 3.3 3.3 3.9 5.7 3.8 5.6 4.5 2.9 3.7 4.9 3.8
##  [163] 4.5 5.4 6.5 4.2 5.0 4.0 5.7 6.0 5.3 3.6 5.3 5.4 4.1 3.6 5.4 3.3 6.4 5.0
##  [181] 3.8 5.3 4.7 4.1 4.3 4.7 4.9 5.2 3.7 3.7 3.9 4.2 4.0 4.4 4.5 4.8 4.5 4.8
##  [199] 3.7 3.9 4.4 4.7 5.5 5.4 3.9 5.7 3.0 4.4 4.9 5.2 6.1 4.5 5.0 4.7 5.7 3.5
##  [217] 3.2 5.5 4.2 4.7 5.3 5.7 4.9 5.1 5.4 4.5 5.2 6.7 4.1 2.9 4.4 5.0 2.9 4.5
##  [235] 3.3 4.0 4.8 3.6 3.2 5.0 4.8 3.2 5.3 6.0 4.6 5.7 3.2 5.6 4.4 4.8 3.1 4.7
##  [253] 4.9 6.2 5.6 3.3 4.3 5.7 3.6 4.1 4.0 3.5 6.2 4.4 5.3 4.4 2.8 4.7 3.5 3.9
##  [271] 5.8 4.8 4.5 5.6 4.7 4.3 4.3 4.6 5.1 5.2 5.8 4.0 4.4 3.4 4.7 4.0 2.5 5.8
##  [289] 4.9 3.9 3.2 3.5 5.4 3.9 4.6 4.6 4.1 4.4 3.3 5.9 4.6 6.5 5.8 4.7 5.4 5.2
##  [307] 3.4 4.9 3.9 5.9 4.8 3.5 6.9 2.1 3.9 5.0 4.3 4.1 5.0 3.1 4.4 5.9 4.5 4.7
##  [325] 5.4 5.9 4.2 4.3 3.5 4.6 4.5 5.4 5.0 4.1 4.1 4.7 4.9 4.2 4.4 4.2 4.7 5.5
##  [343] 5.0 5.3 5.2 3.7 3.1 3.3 5.0 3.6 2.8 4.6 3.3 4.7 4.8 4.5 3.5 3.2 5.0 2.8
##  [361] 5.4 4.8 3.7 3.2 4.3 5.2 5.2 5.2 4.6 3.8 4.7 6.0 5.7 4.4 3.1 5.8 4.5 5.7
##  [379] 4.7 3.7 4.9 5.5 4.3 6.4 4.6 4.3 4.1 4.3 4.1 6.0 3.8 4.1 4.3 2.9 4.4 5.2
##  [397] 5.2 6.4 4.9 4.7 5.1 3.7 3.4 5.4 5.3 4.3 3.1 4.8 4.2 5.8 4.7 3.5 5.3 6.1
##  [415] 3.2 3.6 4.7 6.2 4.0 4.8 5.3 5.5 3.4 4.3 6.1 4.7 3.9 4.3 4.0 6.7 4.8 4.0
##  [433] 4.9 5.2 2.9 6.1 5.3 5.1 4.4 3.6 5.8 4.6 4.2 4.3 4.3 4.5 3.5 3.9 6.0 3.6
##  [451] 5.2 4.9 4.3 3.5 5.5 4.0 4.0 4.6 3.7 4.9 6.6 5.0 4.0 3.9 3.8 4.9 4.2 6.5
##  [469] 4.7 4.3 4.9 4.4 5.3 3.7 5.0 2.2 4.8 3.8 4.3 4.5 4.2 5.2 4.7 3.7 4.6 4.6
##  [487] 3.2 3.8 4.9 5.0 4.8 3.3 6.7 3.5 3.9 4.4 5.3 4.1 3.3 3.4 4.9 5.3 4.2 4.0
##  [505] 3.7 5.8 3.7 3.6 4.6 4.7 2.4 5.0 2.3 3.0 4.2 5.5 3.7 4.9 5.5 3.1 5.8 5.1
##  [523] 6.7 4.5 4.0 4.8 4.1 5.4 4.8 3.7 3.3 5.6 1.7 4.7 4.8 4.6 4.6 5.4 4.4 4.3
##  [541] 4.6 4.0 4.4 5.4 4.1 5.7 5.2 4.9 4.3 5.4 4.2 5.9 4.9 3.7 3.7 4.0 4.0 4.6
##  [559] 4.2 4.1 2.9 3.9 2.9 4.4 3.3 4.6 4.0 4.9 3.0 4.6 4.2 4.8 5.3 5.5 7.0 3.7
##  [577] 4.9 4.4 3.1 4.4 5.1 3.9 4.3 4.7 3.5 5.4 3.0 4.1 3.2 3.7 5.1 3.6 3.8 5.5
##  [595] 4.9 5.0 4.8 4.7 4.3 4.4 1.1 5.1 4.7 4.1 5.1 4.7 4.3 4.9 2.9 3.9 4.0 5.3
##  [613] 5.3 3.3 5.8 3.6 2.6 2.8 5.1 3.3 3.6 5.4 4.4 5.4 4.2 4.1 5.4 4.2 5.5 3.7
##  [631] 4.5 2.9 2.4 6.1 4.0 4.1 3.8 4.8 5.5 4.5 3.0 5.3 5.3 2.9 5.9 4.0 6.0 6.0
##  [649] 3.9 4.1 5.1 4.8 5.9 5.6 4.6 5.2 5.3 3.8 5.3 6.1 3.4 4.9 5.0 4.5 4.5 4.8
##  [667] 5.4 5.3 4.7 3.9 4.4 4.0 3.2 4.9 5.6 3.8 4.8 6.2 3.9 5.1 4.6 4.1 4.5 4.6
##  [685] 6.3 3.1 3.8 4.5 2.2 2.5 4.0 4.0 6.3 3.6 4.1 3.6 3.5 2.4 3.5 5.4 3.0 5.0
##  [703] 3.6 4.6 5.5 4.2 3.7 2.7 4.1 3.5 3.7 3.2 4.5 3.8 5.0 5.9 3.8 5.9 4.9 4.9
##  [721] 3.3 3.1 6.3 5.5 5.5 4.8 3.3 5.7 4.2 5.9 3.9 4.3 5.2 4.3 4.8 4.6 5.0 4.3
##  [739] 2.8 4.6 4.5 4.4 3.4 4.9 6.3 4.0 3.7 4.0 3.9 5.4 5.8 4.2 5.7 4.7 4.8 4.5
##  [757] 3.7 4.2 3.7 5.6 4.0 5.5 3.8 4.1 4.2 5.7 4.5 5.8 4.3 3.2 2.8 3.8 5.2 3.5
##  [775] 7.0 4.0 2.9 4.8 5.4 4.7 4.0 4.9 5.6 3.8 4.9 4.1 4.7 4.4 4.2 4.4 4.0 5.8
##  [793] 5.4 2.2 5.2 5.1 4.9 5.1 5.4 3.3 3.3 4.6 4.9 6.0 4.3 2.7 6.1 4.8 3.4 3.5
##  [811] 4.2 4.4 5.2 4.3 6.4 5.2 3.8 5.1 4.6 4.4 4.6 4.7 3.2 5.5 2.0 4.5 3.3 2.6
##  [829] 3.2 5.7 3.3 6.0 4.7 3.9 5.5 5.2 4.3 4.9 3.6 3.0 4.7 4.5 4.7 4.7 3.4 4.2
##  [847] 4.6 3.9 4.0 5.8 3.3 5.6 6.1 3.0 4.6 5.3 4.9 5.4 6.6 3.4 3.9 5.2 2.9 3.2
##  [865] 4.7 4.2 4.2 5.2 3.9 3.1 4.6 6.0 4.5 3.3 4.7 4.6 4.8 5.4 3.4 4.5 5.1 4.8
##  [883] 4.7 5.8 2.6 4.4 4.3 4.8 3.8 4.9 4.1 4.8 4.1 4.6 4.5 3.3 4.1 4.7 3.6 5.5
##  [901] 2.2 3.3 4.1 3.0 4.3 4.9 4.6 4.2 5.6 5.4 3.0 4.5 3.7 4.5 3.8 5.2 2.9 4.3
##  [919] 6.1 5.6 5.5 5.3 3.9 4.9 5.0 4.6 4.7 5.1 6.3 2.8 5.3 3.9 4.7 3.6 5.1 4.4
##  [937] 3.7 4.8 4.1 5.1 4.0 4.9 5.3 4.2 3.4 4.2 4.4 4.7 5.4 6.7 4.0 2.8 6.2 4.1
##  [955] 4.0 4.0 4.6 3.0 5.8 4.5 4.2 4.6 3.7 4.0 5.5 3.7 3.4 5.6 3.1 4.4 3.6 5.7
##  [973] 5.1 4.3 5.4 3.8 5.1 5.4 4.1 3.3 2.6 4.1 5.8 4.1 2.8 3.4 4.3 4.1 5.4 4.8
##  [991] 3.1 2.7 5.0 4.5 4.8 5.7 4.4 4.5 2.6 5.2
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
##   2.5%  97.5% 
## 2.7000 6.2025
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
##    [1] 4.2 5.1 4.8 5.0 5.4 5.0 3.7 4.2 5.2 3.8 3.9 3.5 4.1 3.4 5.9 4.2 4.2 5.0
##   [19] 4.5 6.5 5.1 5.8 3.7 3.7 4.5 2.3 4.0 3.9 4.1 5.1 4.2 5.5 4.0 2.0 2.9 4.2
##   [37] 5.0 5.2 4.3 5.5 4.5 4.8 4.9 5.5 4.2 4.7 2.6 4.8 6.7 4.5 4.5 4.4 3.7 5.8
##   [55] 4.7 4.8 4.1 3.7 3.1 4.7 4.9 4.1 5.3 5.1 5.2 5.4 5.2 3.4 4.8 4.8 3.5 4.8
##   [73] 4.8 4.5 3.2 4.0 4.7 4.7 5.6 5.3 5.2 5.5 3.2 3.0 5.2 5.8 5.9 5.1 4.8 5.0
##   [91] 4.2 5.2 5.1 3.2 4.7 6.2 4.4 4.8 4.1 5.0 3.5 4.9 3.3 4.6 4.1 5.6 5.0 4.3
##  [109] 5.3 5.2 4.3 4.7 3.7 5.2 5.1 4.1 3.4 4.9 5.5 4.3 4.1 3.3 4.2 4.2 6.3 4.2
##  [127] 5.9 6.4 4.8 3.9 4.9 5.6 5.8 3.4 3.0 4.5 4.7 5.0 6.1 4.2 4.3 3.8 4.4 4.1
##  [145] 4.8 4.5 4.6 4.1 5.4 3.9 4.5 5.1 4.3 5.1 4.8 3.9 4.3 4.8 5.2 4.7 6.2 4.4
##  [163] 5.2 3.5 2.0 4.4 5.7 5.6 6.2 3.0 4.5 3.9 4.2 5.1 3.9 4.1 5.2 5.7 3.8 5.0
##  [181] 5.6 5.5 2.1 4.6 3.8 5.4 3.9 4.8 3.2 3.8 4.7 3.6 3.6 6.3 4.9 4.3 4.5 3.2
##  [199] 5.7 3.8 4.8 4.2 4.6 5.6 6.2 5.3 4.7 4.8 4.5 3.5 4.2 5.6 3.9 4.1 3.1 5.7
##  [217] 6.4 4.7 3.6 2.7 3.2 3.3 4.5 3.9 2.4 3.4 3.6 5.0 4.6 4.6 2.6 4.8 4.9 5.0
##  [235] 4.1 4.5 4.3 4.5 6.1 4.3 4.0 4.4 2.9 5.3 4.7 3.6 4.4 5.3 6.1 4.2 3.9 5.4
##  [253] 2.7 4.2 4.7 5.2 6.3 5.6 3.4 6.4 3.5 2.7 4.5 5.6 4.3 4.7 5.7 4.2 3.8 4.5
##  [271] 4.7 2.9 3.9 4.8 3.9 3.6 4.7 3.3 4.8 5.3 3.9 5.0 6.2 4.3 3.1 5.5 3.4 5.8
##  [289] 3.4 4.8 4.8 3.2 4.8 4.6 3.1 4.7 4.8 3.8 3.8 4.6 3.3 6.2 4.9 4.1 4.5 5.1
##  [307] 3.4 4.3 5.2 4.6 3.8 5.2 4.9 4.1 2.8 2.4 4.4 4.6 4.6 5.7 3.7 3.3 4.6 3.0
##  [325] 5.1 5.1 5.7 5.8 4.1 4.6 5.7 2.8 4.1 3.9 4.5 3.4 4.4 4.2 5.6 5.6 6.1 5.8
##  [343] 4.4 4.3 5.1 2.4 5.9 3.8 5.1 4.8 5.7 4.4 3.7 3.8 3.2 4.4 4.4 3.2 3.7 5.4
##  [361] 4.0 5.5 3.7 3.1 4.7 4.6 5.7 4.7 2.9 5.2 5.2 4.2 3.3 5.9 4.9 5.1 4.6 4.4
##  [379] 5.6 5.3 4.7 5.4 5.1 3.9 4.2 3.4 4.4 5.0 5.7 3.4 7.9 3.7 4.5 4.7 6.1 3.1
##  [397] 3.5 5.1 4.1 3.9 5.2 5.6 6.0 6.1 3.8 6.2 4.9 4.9 3.1 2.7 5.6 3.8 4.4 4.5
##  [415] 4.3 2.5 5.1 4.2 3.5 4.6 3.7 2.7 5.0 4.2 5.4 5.7 5.4 4.2 5.9 4.5 4.1 4.0
##  [433] 4.2 2.4 4.7 4.7 4.9 5.8 2.9 4.8 6.1 5.1 4.3 5.5 6.4 4.8 4.1 5.0 3.5 5.2
##  [451] 4.0 5.3 4.3 5.3 4.0 5.3 4.7 4.2 3.4 3.8 6.0 3.5 5.0 5.5 2.9 5.5 5.0 4.9
##  [469] 4.6 3.9 5.0 4.3 4.5 4.8 5.1 4.6 2.6 3.6 5.3 5.1 4.1 3.9 5.5 4.7 5.4 6.1
##  [487] 3.6 3.7 3.6 6.4 4.3 4.0 3.3 3.6 5.2 5.7 4.3 5.1 3.8 4.2 3.2 4.3 5.4 3.6
##  [505] 3.9 4.3 5.3 5.6 4.2 5.4 4.8 3.9 2.8 4.5 3.3 4.4 2.0 5.0 4.1 4.5 5.2 2.4
##  [523] 5.1 5.1 5.0 3.8 5.6 5.5 4.8 3.8 4.8 5.0 4.3 4.7 4.9 6.2 4.0 3.6 3.7 4.7
##  [541] 3.6 4.2 4.4 5.2 4.4 6.5 3.5 5.4 3.8 5.0 1.9 3.9 3.2 3.5 3.6 4.6 5.1 4.4
##  [559] 5.5 5.3 4.8 5.6 3.8 6.9 4.0 4.7 3.9 3.7 5.8 4.4 3.5 5.6 3.6 6.2 3.0 3.0
##  [577] 2.8 3.6 4.9 5.9 4.7 4.9 5.0 5.4 5.3 4.2 4.2 3.6 4.8 3.2 4.5 4.0 5.1 5.9
##  [595] 3.5 4.3 3.6 4.5 3.5 3.8 3.2 4.7 4.8 4.5 5.1 4.7 4.4 3.9 4.4 5.6 4.4 5.2
##  [613] 3.3 5.8 4.3 4.2 4.9 3.5 5.1 5.1 5.6 4.9 5.1 3.7 5.4 6.0 4.5 4.5 5.1 4.9
##  [631] 2.8 6.8 5.9 3.0 3.2 3.8 4.8 4.6 4.5 5.7 5.3 5.2 5.3 2.6 4.7 5.4 4.9 4.9
##  [649] 5.2 3.0 5.7 4.5 4.7 4.5 3.6 4.0 5.9 4.5 4.6 5.8 4.5 6.3 5.6 4.8 3.9 5.6
##  [667] 3.7 6.6 5.8 5.3 3.7 3.0 4.6 5.4 4.9 5.0 4.1 6.4 4.5 4.6 4.4 5.5 6.1 3.2
##  [685] 5.0 4.6 6.5 4.7 5.0 5.4 3.7 4.5 6.3 4.2 4.5 3.4 3.9 4.8 4.7 6.2 5.3 5.8
##  [703] 3.9 5.3 5.4 4.2 3.6 4.7 4.6 5.5 4.3 4.3 4.6 4.5 2.9 4.9 4.8 4.3 6.8 4.6
##  [721] 4.4 2.7 5.6 4.9 3.8 4.0 4.5 5.2 4.7 4.0 5.1 3.2 4.9 4.0 6.4 5.3 4.1 4.5
##  [739] 4.8 6.1 4.0 4.7 4.2 5.0 5.8 4.5 4.5 5.1 5.3 4.7 6.2 3.8 5.2 3.7 5.0 3.2
##  [757] 5.4 3.8 5.5 3.4 7.0 5.1 3.3 5.1 3.8 5.5 4.4 4.7 3.5 4.6 3.8 6.2 4.9 4.7
##  [775] 5.9 3.7 5.1 5.7 4.4 5.0 5.6 4.8 3.9 4.1 3.9 4.1 4.6 4.7 5.5 3.7 3.7 3.8
##  [793] 4.2 4.3 5.3 3.5 5.2 6.6 3.6 5.2 4.8 5.6 4.7 3.9 2.9 3.7 4.5 3.9 6.1 3.1
##  [811] 5.5 6.4 4.1 5.1 4.4 4.5 3.4 4.3 5.3 4.0 3.1 5.5 4.9 4.3 5.1 4.5 4.8 4.1
##  [829] 4.0 6.0 4.4 4.9 3.8 3.9 4.0 4.9 4.0 3.5 4.3 3.0 5.0 4.9 4.7 4.4 5.0 4.8
##  [847] 3.7 5.8 5.7 5.0 5.9 3.6 4.0 4.3 6.5 4.0 5.9 4.8 4.8 5.0 6.3 3.8 4.3 3.5
##  [865] 5.0 5.3 3.5 5.2 5.2 4.3 4.3 3.6 6.6 5.0 4.1 4.5 4.9 5.5 4.7 4.3 3.0 5.2
##  [883] 5.1 4.1 3.8 4.8 4.7 5.4 3.3 5.4 4.9 4.1 4.0 4.6 3.3 4.0 2.7 4.7 4.2 6.9
##  [901] 2.9 3.2 4.7 4.8 5.5 3.2 5.9 3.3 5.4 5.3 5.4 3.6 4.4 4.0 5.4 5.3 3.7 3.6
##  [919] 3.4 3.4 3.4 4.7 5.2 7.2 4.0 5.4 3.6 5.6 3.0 4.4 4.4 3.7 3.9 4.9 4.8 6.5
##  [937] 5.5 3.7 4.7 3.6 5.8 4.5 4.4 4.5 4.7 5.5 3.5 5.5 5.1 6.6 3.5 3.9 4.6 5.5
##  [955] 4.3 4.0 4.0 4.1 4.3 3.2 5.9 3.1 4.2 5.0 5.4 4.6 5.0 3.9 4.4 4.5 3.8 5.2
##  [973] 3.0 4.5 4.7 5.7 3.5 5.0 4.9 3.5 4.7 5.9 5.4 5.5 2.2 5.8 5.3 2.0 3.6 4.1
##  [991] 4.1 3.9 5.7 4.5 4.1 4.3 3.1 4.2 4.8 3.1
## 
## $func.thetastar
## [1] 0.0467
## 
## $jack.boot.val
##  [1]  0.57000000  0.42240437  0.32634561  0.19349112  0.10408163  0.02112676
##  [7] -0.17791045 -0.25434174 -0.39107692 -0.42463343
## 
## $jack.boot.se
## [1] 0.9826874
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
##    [1] 4.3 3.7 4.2 4.9 4.7 5.0 5.7 4.6 5.2 5.2 4.6 4.4 4.6 5.2 4.1 4.1 5.4 5.5
##   [19] 4.5 4.6 5.4 6.7 4.1 5.4 3.8 4.0 6.8 4.1 4.8 3.6 4.5 3.9 6.0 4.0 5.4 5.0
##   [37] 4.5 4.3 3.9 3.6 5.0 5.8 3.2 3.6 2.6 3.1 4.7 3.7 2.6 4.3 4.9 4.6 4.2 2.5
##   [55] 4.4 4.8 4.3 3.5 6.1 2.0 3.7 3.3 3.5 2.4 2.8 5.9 3.7 4.1 4.6 4.7 3.8 3.7
##   [73] 6.2 6.2 4.0 3.7 4.5 4.2 4.0 6.4 3.7 4.6 4.8 4.6 5.3 3.8 3.8 4.0 3.9 5.3
##   [91] 6.1 5.9 4.2 4.2 4.9 5.6 3.2 5.6 6.1 5.1 5.5 4.7 4.5 3.7 3.0 3.5 3.6 3.3
##  [109] 3.1 5.0 2.1 4.7 4.6 4.5 4.8 4.6 3.3 2.8 5.7 5.4 4.0 3.7 4.0 4.2 5.2 3.7
##  [127] 2.9 4.3 3.8 4.3 2.7 4.5 3.4 4.6 7.2 5.0 3.3 4.6 3.8 5.3 3.1 4.3 4.3 5.5
##  [145] 4.4 4.2 4.8 5.6 4.1 4.8 3.5 4.8 5.1 6.4 3.5 4.1 4.5 4.9 5.4 3.2 5.1 5.7
##  [163] 3.1 6.0 4.8 4.0 3.7 3.4 3.7 4.9 4.5 5.0 3.1 4.4 5.1 3.4 4.5 4.0 5.4 3.0
##  [181] 5.1 3.6 4.1 4.8 3.6 4.8 5.7 6.7 5.2 4.7 4.5 4.0 3.3 5.3 5.2 2.3 6.4 3.1
##  [199] 3.1 5.9 6.1 4.8 6.3 5.0 5.3 4.2 3.8 4.7 3.2 3.2 4.9 4.6 3.5 3.3 2.6 5.2
##  [217] 3.6 5.0 3.2 3.9 4.7 5.9 5.9 4.8 5.8 4.2 4.6 4.5 4.8 3.2 4.9 4.4 4.3 4.9
##  [235] 4.5 3.7 4.3 4.5 4.6 4.6 5.0 4.0 6.7 4.0 4.9 6.1 4.2 4.7 4.1 6.0 2.3 2.7
##  [253] 3.7 4.6 4.6 4.7 4.5 4.0 6.2 4.3 5.4 4.0 3.2 4.5 2.9 2.6 4.5 4.5 4.4 4.5
##  [271] 5.2 3.9 3.5 4.0 2.6 5.9 4.0 4.3 4.4 3.9 4.4 5.8 4.8 4.8 5.4 4.9 4.1 4.0
##  [289] 3.9 2.7 3.7 4.3 5.8 3.6 3.9 4.4 5.6 5.6 5.1 2.6 4.3 3.8 5.0 5.5 5.2 5.7
##  [307] 3.2 5.9 3.3 4.5 5.6 4.3 1.5 4.2 6.2 4.4 6.5 4.2 3.6 3.5 4.3 4.7 6.6 4.4
##  [325] 3.8 4.6 5.6 5.6 5.5 5.7 4.9 4.9 5.2 2.9 5.8 3.9 4.5 3.8 4.3 5.0 5.6 3.9
##  [343] 5.3 3.9 4.9 3.9 3.6 4.4 4.6 3.5 3.5 3.6 4.3 5.3 5.8 5.0 3.4 4.2 3.9 2.7
##  [361] 3.3 3.8 4.3 3.6 4.2 5.9 5.5 3.2 5.7 5.1 2.8 3.5 4.6 3.6 3.2 3.6 4.4 3.1
##  [379] 4.8 3.3 4.2 4.1 3.7 5.5 4.5 3.8 3.7 3.2 5.1 4.3 5.2 5.8 2.1 5.3 5.2 3.6
##  [397] 4.2 5.0 5.7 4.8 5.2 3.7 3.9 6.2 4.9 3.9 3.4 4.6 2.7 4.1 4.1 4.6 4.6 4.5
##  [415] 6.0 4.2 4.1 4.6 5.2 3.6 4.4 5.4 4.4 4.8 2.6 4.7 4.2 5.1 5.5 5.9 4.2 3.7
##  [433] 3.6 4.1 4.4 3.7 5.4 3.7 4.8 5.0 4.8 6.5 3.7 3.8 5.2 4.4 1.2 5.1 3.8 4.2
##  [451] 4.1 4.7 5.2 3.5 5.7 4.2 4.6 4.7 4.8 5.4 4.9 3.1 4.0 4.2 3.7 5.1 2.7 5.4
##  [469] 2.9 4.5 5.8 4.4 4.9 3.8 4.4 4.5 4.1 5.7 2.6 2.2 5.0 4.0 5.6 5.0 5.5 5.2
##  [487] 3.9 5.2 3.8 4.6 3.7 5.6 5.5 4.6 4.4 6.4 4.9 4.6 4.4 3.9 5.1 5.7 6.0 5.0
##  [505] 4.6 4.7 5.4 4.4 3.9 3.8 5.4 5.7 3.1 5.7 5.8 5.0 5.3 5.6 4.1 3.0 3.9 4.8
##  [523] 4.2 2.6 4.8 3.6 6.1 4.1 5.5 3.2 5.7 2.9 3.5 4.4 4.2 5.3 3.2 5.2 5.6 4.2
##  [541] 5.3 5.8 4.8 4.6 3.5 4.8 6.6 4.5 4.4 5.5 3.1 4.6 2.8 5.8 6.6 4.4 4.3 3.4
##  [559] 3.8 3.4 5.8 3.6 3.9 3.0 4.1 4.6 5.3 3.5 4.9 5.5 4.8 4.6 4.6 5.5 4.2 4.6
##  [577] 4.4 5.2 4.5 4.1 4.8 3.8 5.2 6.0 4.7 4.1 4.4 3.3 4.7 4.0 4.9 5.4 4.7 3.3
##  [595] 3.5 4.3 4.2 6.2 3.3 4.4 5.9 4.4 3.1 4.2 4.8 5.3 5.2 3.9 4.1 3.9 6.8 5.1
##  [613] 4.3 3.5 3.2 4.0 4.9 5.5 4.2 4.1 4.9 2.6 3.4 4.3 2.7 3.4 2.9 4.6 5.3 5.0
##  [631] 3.8 3.9 4.6 4.7 3.5 5.1 4.9 4.8 3.9 3.3 4.4 4.1 5.9 3.1 4.6 3.2 4.7 4.5
##  [649] 4.8 3.2 5.2 5.7 5.5 5.4 4.0 6.6 3.9 3.9 4.4 3.6 4.7 5.1 3.9 4.4 4.2 3.6
##  [667] 3.8 4.2 4.5 2.9 3.7 4.3 5.0 2.4 4.6 2.9 2.9 4.2 4.9 5.5 5.5 2.2 6.1 3.6
##  [685] 4.3 3.9 5.6 3.8 5.7 4.5 2.7 4.5 3.8 4.1 4.8 4.3 3.0 3.9 5.9 4.7 4.2 4.4
##  [703] 4.3 3.7 4.6 4.8 4.6 5.0 4.1 5.8 2.3 4.2 5.0 3.7 4.4 4.2 4.4 6.7 5.9 6.0
##  [721] 5.5 5.1 4.1 6.6 6.4 3.5 4.7 3.8 5.0 4.7 3.7 5.0 4.5 4.0 6.0 3.3 3.6 3.6
##  [739] 5.2 4.5 3.3 4.2 5.5 5.1 3.1 4.4 5.0 5.1 4.5 4.8 3.6 4.5 4.3 5.7 4.6 3.1
##  [757] 4.9 5.2 3.6 4.9 4.1 4.1 5.4 4.5 4.3 4.9 6.3 4.3 5.3 3.4 5.4 4.8 5.6 3.7
##  [775] 3.8 4.8 5.2 4.7 4.4 4.8 4.9 7.4 6.6 3.5 5.0 3.6 3.7 4.2 3.1 4.4 3.8 4.4
##  [793] 3.9 4.3 5.5 4.8 4.9 3.8 3.2 4.5 6.7 3.8 3.1 3.5 5.2 4.5 5.8 4.4 3.7 3.7
##  [811] 3.6 5.8 5.0 5.6 5.9 4.3 5.6 4.2 3.9 3.3 3.3 4.3 5.6 5.2 6.1 4.1 5.1 5.5
##  [829] 5.4 4.3 6.1 3.8 5.4 4.1 5.7 4.1 5.4 5.1 3.9 3.5 6.8 6.4 4.2 6.0 4.5 3.6
##  [847] 4.9 4.4 5.1 5.4 3.9 5.6 3.4 5.3 3.9 3.2 3.9 3.4 4.5 3.4 3.6 5.6 4.9 2.2
##  [865] 4.3 4.4 3.7 4.6 4.9 4.3 2.1 5.3 4.0 4.4 4.8 4.1 3.9 3.6 4.2 6.0 4.8 5.8
##  [883] 5.7 5.0 5.2 5.1 4.1 6.1 3.1 5.4 4.1 5.2 4.5 4.3 4.2 2.8 4.2 4.4 4.1 3.9
##  [901] 3.9 5.0 4.0 5.8 4.2 5.5 6.0 3.9 3.8 3.2 2.6 3.7 5.9 5.5 5.0 5.4 5.7 5.8
##  [919] 4.1 4.4 4.6 5.0 4.8 4.7 3.7 4.3 3.1 3.9 4.3 5.5 4.9 4.2 2.7 5.3 5.3 5.8
##  [937] 4.3 4.4 4.9 5.4 3.5 5.4 4.6 4.7 5.8 3.1 3.7 4.6 4.9 5.1 3.4 4.9 3.7 4.5
##  [955] 4.9 2.7 3.3 5.1 6.1 3.7 2.5 4.5 4.2 3.7 5.9 4.4 5.5 5.1 2.8 4.6 4.6 4.0
##  [973] 4.3 4.1 3.3 3.9 3.6 4.7 4.9 4.5 3.2 4.4 4.7 5.3 4.0 3.5 3.9 5.2 4.5 4.1
##  [991] 5.2 4.1 5.9 4.7 5.0 3.9 5.1 3.9 3.5 4.0
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.500 5.500 5.256 5.300 5.200 5.000 4.920 4.800 4.600 4.400
## 
## $jack.boot.se
## [1] 1.056102
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

<img src="Week-2-lab_files/figure-html/unnamed-chunk-17-1.png" width="672" />

Parametric bootstrap
---------------------

Let's do one quick example of a parametric bootstrap. We haven't introduced distributions yet (except for the Gaussian, or Normal, distribution, which is the most familiar), so lets spend a few minutes exploring the Gamma distribution, just so we have it to work with for testing out parametric bootstrap. All we need to know is that the Gamma distribution is a continuous, non-negative distribution that takes two parameters, which we call "shape" and "rate". Lets plot a few examples just to see what a Gamma distribution looks like. (Note that the Gamma distribution can be parameterized by "shape" and "rate" OR by "shape" and "scale", where "scale" is just 1/"rate". R will allow you to use either (shape,rate) or (shape,scale) as long as you specify which you are providing.

<img src="Week-2-lab_files/figure-html/unnamed-chunk-18-1.png" width="672" />


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
## [1] 0.8864666
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
##   1.9500252   3.8650113 
##  (0.8085797) (1.8260561)
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
## [1] 0.3590103 0.7239850 2.1095803 1.7624179 1.9803335 1.5759712
```

```r
hist(results,breaks=30,col="pink",ylim=c(0,1),freq=F)
```

<img src="Week-2-lab_files/figure-html/unnamed-chunk-22-1.png" width="672" />

Now we have the bootstrap distribution for skewness (the $\theta^{*}$ s), we can compare that to the equivalent non-parametric bootstrap:


```r
results2<-bootstrap(x=original.data,nboot=1000,theta=skewness)
results2
```

```
## $thetastar
##    [1]  1.0945064054  0.3162993281  0.9088319394  0.3387727637  0.6482379481
##    [6]  1.0914877418  1.1671603757  0.1664445562  0.4881364608  0.9750884519
##   [11]  1.0369624897  1.2309660585  0.7007297624 -0.2840808169  1.3626144063
##   [16]  0.1551648983  0.8096161821  1.1146245932  0.5960230466 -0.9282628011
##   [21]  0.9207310697  1.5593727253  1.5126154103  1.4053135887  1.3905879838
##   [26]  0.1476832444  1.2999590989  1.2313546828  1.6835583224  0.6575052919
##   [31]  1.2137534695  0.7469437241  1.0102493537  0.5673555822  1.0144508800
##   [36]  1.1567371753  0.0544417567  0.2459536531  1.0921510574 -0.0005671645
##   [41]  0.2675132055  0.7471139930  0.9817076499  0.9568062732  0.9211299739
##   [46]  0.0605523953  1.0527596755  1.2641232464  0.5622090699  1.1476016161
##   [51]  0.2574896981 -0.4353297700  0.9903489550  1.0844599617  1.0604993496
##   [56]  0.3113737648  0.3574811576  1.0357498539  0.5610407169  1.0945525574
##   [61]  0.5050114512  0.7566249733  0.8472045908  1.3278685861  0.0230222110
##   [66]  0.7528737581  0.1119943303  0.1349781561  1.0163289602  1.1937099112
##   [71]  1.1057308594  0.8871769045  0.3294964101  0.3308593549  1.4140762686
##   [76]  1.0385237638  0.1882019858  1.7651113540 -0.0780312805 -0.3866736934
##   [81]  0.6905090461  0.9621603137  0.0640612477  1.3934949890  0.5651377571
##   [86]  0.9790658991  1.0100280718  0.7017356497  0.6481526647 -0.3500640195
##   [91]  0.7036963444  0.8617605809  0.0696773244  1.1516520057  1.1361141476
##   [96]  0.5744298736  0.8708456346  1.6393119562 -0.4858047836  1.4220244476
##  [101]  0.9500163339  0.4800802078  0.6518094206  0.9500672894  1.2200166582
##  [106]  0.5208606848  1.5525246113 -0.0990644922  1.4711597732  0.9913417863
##  [111]  0.7286581225  0.9819779002  0.6393091357 -0.4886884908  0.3432851488
##  [116] -0.1328802672  0.6783055079  0.6482366569  1.1134577018  1.6340696888
##  [121]  0.2505297149  0.5960230466  0.4612082535  0.5303773522  0.6588323755
##  [126]  0.2231911509  1.4474062370 -0.8936927806  0.6918595306  0.0960315238
##  [131]  1.5746888231  0.3074646312  1.4909637087  0.6669789690  0.5785275798
##  [136]  0.5673055420  0.7213697511  1.1728358111 -0.1613169029  0.0980816400
##  [141]  1.2716861105  0.8376620368  0.8473219066  0.5595630762  0.5885747844
##  [146]  0.8034692578  1.0100981120  1.1506054671  1.1109563667  0.6450636817
##  [151]  0.4925763553  0.1508541066  0.8336384738  1.3227575159  0.4981529482
##  [156]  0.3093650582  0.0459672749  0.4956392910  0.5670295579 -0.1023891382
##  [161]  0.8175919621  0.8365029355  0.9962216518  0.5207428516  1.1976256891
##  [166]  0.2296895246  0.6027603415  0.6743328806  0.6401980517  0.4291063972
##  [171]  1.1451614359  0.9550495413  0.9329299023  0.2408484258  1.0535183154
##  [176]  0.4681651587  0.8304358769 -0.0589900441  0.7651138418  1.1351758258
##  [181]  0.4605988157  0.6751964705  0.3974606641  0.4811671930  0.5998525487
##  [186]  0.5761388989  1.0791463628  0.8303387805  1.6587490830  0.6975098553
##  [191] -0.1114418754  0.0497109845  0.9582928148  1.8197913169  1.0501210586
##  [196] -1.7293070578  0.8175227968  0.7999409032  0.8252477048  0.8023313419
##  [201]  0.2048430661  0.8772109553  1.6219630302  1.3870961906  1.0525517147
##  [206]  0.1702462284  1.3486908486  0.0205715996  0.5087349213  0.2610085648
##  [211]  1.2511754918  0.7068754485  1.6163386041 -0.3381753861  0.6276809286
##  [216]  1.6440404404  0.6390878786  1.0162322435  1.0628122044  1.1255675841
##  [221]  0.7734676948  1.7699102059  0.4965427296  0.7615338194  0.6571403017
##  [226]  0.4814078863  0.1680199013  0.7767570744  1.6441129189  0.8075369425
##  [231]  1.5362062614  0.9490148116  0.5896339821  0.5704004263  0.6138533713
##  [236]  0.8177799152  0.6407408621  1.0527309421 -0.8121307577 -0.6231304800
##  [241]  0.4252995153  1.4677352732  0.2781159854  0.9693127049  0.2913988107
##  [246]  0.7528675873  1.7999978785  1.6851011083  1.1793611794  0.8077390745
##  [251]  1.0262830428  0.3899205767  1.2198439020 -0.1776754278  0.8104068372
##  [256]  0.1290207227  1.1530880540  0.1790880815  0.6643936486  0.7994278155
##  [261]  0.1269852532  0.7232271149  0.7829280761  0.8456704861  0.3056474732
##  [266] -0.0282444223  0.7737700074  1.4330537458  0.7360318314  2.0948156404
##  [271]  0.8746244807  0.7737449228  0.5428314325  0.5504765761  0.9700628056
##  [276]  0.4027424080  0.5679771046 -0.1640146980  0.6462837724  0.6585113686
##  [281] -0.7615865907  0.6571235473  0.6805237718  0.6322535136  1.5206463704
##  [286]  0.9700628056  1.4151162424  0.6858983941  0.9561284388 -0.2616321006
##  [291] -0.3757153293  0.9512118367  1.9651601333  1.2846625035  0.8629876243
##  [296]  1.4191290746  0.1941138228  0.6011528900  1.2841179405  0.1690433977
##  [301]  0.3327678210  0.4407449148  1.0098407673  1.0577754230  1.6314903959
##  [306]  0.6610349373  0.7520433178  0.9412221780  0.6003014627  0.5120114895
##  [311]  0.3672577674  0.0924936147  1.0468535857  0.5819818759  1.2167216191
##  [316]  0.7112997395  0.0971077439  0.8709432257  1.3400220350  1.4414665196
##  [321]  1.2843861546  0.8792676102  0.1294894611  0.6645836686  1.0928611747
##  [326]  0.8953314253 -0.5854759390  1.0236006894  0.4810579394  0.6552375184
##  [331]  1.0663912415  1.0967738385  0.2696058488  0.7903782024  1.1867425068
##  [336]  1.0794809762  1.3865661157  0.7720854154  0.5306439061 -0.1571671654
##  [341]  0.8639251451  0.7005826936  0.8187878998  0.8408585530  1.0577754230
##  [346] -0.1190028176  1.7399194176  0.1181821696  0.5064863072  0.6038274327
##  [351]  0.1925847845  0.0485365674  0.9271595618  0.4072130279  0.0243838198
##  [356]  1.1664484928  1.1723791304  0.7175119884  1.5569521981  1.4606054645
##  [361]  0.8171681969  0.9207807930  0.2666751746  1.9245732167  0.5931808904
##  [366]  1.1140982326  0.4510065042  0.7462636902 -0.2993508021  0.9354819548
##  [371]  0.2260078721  0.8328417646  1.5926801876  1.9174312802  0.7902994994
##  [376]  0.4662836490  0.6643936486  0.4217893518  0.8076887375  0.4372280909
##  [381] -1.5436342458  0.6132412442 -0.6342827031  0.7530340507  0.0346087837
##  [386] -0.0246358738  1.0798288057  0.7856386854  1.2177493098  0.3110571237
##  [391]  1.1175144999  0.9137142883 -0.3030288449  0.8202908959  0.7279452092
##  [396] -0.6364367717 -0.3486784795  0.5381554426  0.9720562906 -0.0765229805
##  [401]  1.2379722111 -0.0670025147  1.6593969742  0.7650688490  0.6290297471
##  [406]  1.6252937043  1.8395355336  1.4922434030 -0.1470181518 -0.4149644274
##  [411] -0.4134030523  0.3308648482 -0.0415239279  0.1505490147  0.2363448372
##  [416]  0.9574447773  0.2354304299  1.1183000342  1.7228365348  0.6394571141
##  [421]  0.9486525297  0.4743169508  1.9804792810  2.2046009192  1.2228015552
##  [426]  0.4440186923  1.3052548515  0.2094967819  0.2019548257 -1.0108816012
##  [431]  1.1506054671  0.1660825220  0.5460978634  1.3693547854  1.3700309272
##  [436]  0.4531630526  0.5429172243  0.8672504002  0.0938839187  0.7996744877
##  [441]  0.4692299804  0.6608623697  0.8847063973  0.5673863821  0.6545262203
##  [446] -0.1972024860  1.6975353115  0.1341880966  1.9634450718  0.4102944546
##  [451] -0.1838245199  0.7589488841  0.2632113895  0.7426285178  0.8805649421
##  [456]  0.3383450808  1.0134635906 -0.3009612596  1.4535304356 -0.6319513980
##  [461]  0.3105626520  0.5479193593  0.4163337536  0.2405967307  1.0754315584
##  [466]  1.1938887151  0.9863171457  0.7689462510  1.7546437272  0.0179632545
##  [471]  0.0391443489  1.6504172826  0.7770129235  0.5103307550 -0.0765006467
##  [476]  1.1725183169  0.6859529816  0.3425579626  0.3099866062 -0.9800639877
##  [481] -0.1089697502  0.6207473278  0.2909639183  0.6191438988 -1.1796144275
##  [486]  0.9512118367  1.0333002228  0.6533744938  0.6228233154  1.1739675604
##  [491]  1.1390032383  2.3664650137  0.4045923517  0.0149703720  0.4886046533
##  [496]  1.1517498707  0.7947667776  1.0385546736  1.5641376564  1.1075157932
##  [501]  1.7220942282  0.9228696667  1.4214940170  0.1310680155 -0.3124167350
##  [506]  1.1208183544  1.0194878410 -0.0976810834  0.6075789708  0.7528737581
##  [511]  0.3311493333  0.4746309821  1.3989983266  0.5617828834 -0.1696541786
##  [516]  0.8696495150  0.3998054919 -0.2035399127  1.0513309711  0.1300608700
##  [521]  0.9876894621  0.7414152175  0.8914530961  1.5950585118  1.0615910062
##  [526]  0.7826633222 -0.0804602632  0.5150918220  0.0340878972  1.4803060688
##  [531]  0.8371696382  1.4419700729  1.6081146851  1.2134522989  0.8569558984
##  [536]  0.6565906093  0.7142722826  1.0364877347  1.2775030864  1.2797808619
##  [541]  0.3235514708 -0.7933218502 -0.0370359192  0.6815213326  1.4802836848
##  [546]  0.6281846046 -0.8274756199  0.8098346387  0.1021003785  1.7444482451
##  [551]  0.4272150773  0.7277748673  0.5532928817  1.2616795867  0.1037403459
##  [556] -0.2469273869  0.6416716772  0.7605163195  0.9619711725  1.2846115771
##  [561]  0.2649734827  0.9691967994  1.4702215910  0.1874382806  0.9190656102
##  [566]  0.6350302502  0.3501768035  0.7650203997  1.2595850853  0.5329938487
##  [571]  1.0608947149  0.7413579487  1.6487969439  0.1192335930  0.5133294139
##  [576]  1.4447284422  0.8983643651  0.7246478805  0.2872149796 -0.2010571064
##  [581]  0.9348817494  0.1269011976  0.3729485466  1.1154710375  1.2650883978
##  [586]  1.3839996052  1.5387230659  0.5438717205  0.9221062708  0.7680985817
##  [591]  1.5207848667 -0.0665189484  0.8505240206  0.8183147509 -0.5985363263
##  [596]  1.4855121921  0.6773467464 -0.1355674250  0.4624833886  1.7012802277
##  [601]  0.9696011679 -0.0235223533  0.7462675575  0.4539583154  0.5485056104
##  [606]  0.5734332553  0.4097726683  0.6898577289 -0.2280308250  0.5203858585
##  [611]  0.5020211727  0.1779265629  0.8253778267  0.1224073888  1.1603640381
##  [616]  0.1541548414  1.2511754918  1.3139554012  0.5960068578  0.7945678829
##  [621]  0.2390597884  0.8056584123  0.5969592449  1.4999372676  1.3702824735
##  [626]  0.4982470205  1.7957437590  0.5363477461  0.4800196121 -0.2289282076
##  [631]  1.7424237941  1.0158329885  1.0509144845  0.7455599698  1.5804296385
##  [636]  0.5753133017  1.3176603935  0.9997957141  0.6337269182 -0.1473862290
##  [641]  0.4321180125  0.2690565039 -0.1755634269  1.1397932872  1.0615077501
##  [646]  1.1380110103 -0.3546556132  0.5656769333 -0.9260403362  0.6488138593
##  [651]  0.2765917563  1.5350959207  1.3630939355  0.9241617766  0.8367707091
##  [656]  0.2968464352  0.3482769621  0.6919481620  0.7251904851  0.2202939227
##  [661]  1.3827045924  0.7942374565  1.3882883829  0.4540376445 -0.4964856597
##  [666]  1.6827587180  0.1145203032  0.8769240360  1.1486219818  0.7401091605
##  [671]  0.5823688819  1.1218288699  0.6487633210  1.6189017361  0.0369831969
##  [676]  0.6886799508  0.5743454538  0.2307568506  0.2539257616  0.9929756803
##  [681]  0.5964349029  1.5470506875  0.8281290841  1.0962796637  1.2812292836
##  [686]  1.4443817026  1.7067164174  0.8122473182  0.7186314931  0.7726472356
##  [691]  0.1161701165  1.1187354082 -0.1837692693  0.8846585635  1.4121850063
##  [696]  1.7001119133  1.4265647930  0.5542508764 -0.0407423566  1.0642117545
##  [701]  1.1836592512  0.6211322853  0.6269385323  0.5504636978  1.1245226511
##  [706] -0.3127169198  1.0772758306  1.2127835112  1.0273144936  0.8032733394
##  [711]  0.7676847902  1.0790216770  0.3728576076  1.0791966269  1.4876917957
##  [716]  0.8625427627  1.4907520011 -0.4693036925  0.4714194202  0.9228339985
##  [721]  0.5349625916 -1.0473977146 -0.1414327475 -0.1329821048  0.4492343864
##  [726]  0.4253140200 -0.4646303934  0.8188730156  0.5573529195  0.4800156120
##  [731]  0.7312779142 -0.0183091063  1.1799666644  0.9833390330  0.7519025572
##  [736]  0.8855319053  0.5014790474  1.2300702002 -0.0299924576  0.8410628347
##  [741]  1.0429053852  1.5083101980  1.1410863712  1.7614395508  1.9813702899
##  [746]  0.7652953359 -0.1076808337  0.6803082070  0.8517757933  0.1346645794
##  [751]  1.0466914843  1.0908456830  0.5677095364  0.7988898476  0.3607940804
##  [756]  0.2010999470  0.4434144186  0.0411584989  0.5043670897  0.3695247918
##  [761]  2.2761393631  0.0656060690  0.9583606293  0.0179369188  0.3117275987
##  [766]  2.1429333879  0.7244424031  1.0279220823  0.5708205889  0.1781726758
##  [771]  1.8216210123  0.2340258286  0.5895919603  1.2355025627  0.9848142141
##  [776]  0.4893863176  0.0600452528  1.2196460705  1.0930608788  1.1279156182
##  [781]  0.2539257616  0.6179138490 -0.7418568614  1.8255436954  1.0617363235
##  [786]  0.5195789397  1.4786662512  0.3798075743  0.9717562275  0.8635860860
##  [791]  1.6572169788  1.6731746086  0.0582137797  1.1534152849  0.9973535963
##  [796]  0.6311551550  1.5011765071  0.2276355541  0.6320663592  0.5161889443
##  [801]  1.7761559553  1.6184595700  0.7588271448  0.7377684001  2.0891114098
##  [806]  1.2221277705  0.4021465117  0.7869823791  0.2744553532  0.7364570114
##  [811]  0.6548055936  0.5764292270  1.9905726246  0.6169035231  0.1196897524
##  [816]  1.5185490664  0.4846030476  1.0461386300  1.4347211351  1.0302751872
##  [821]  0.6401536355  0.9081129853  0.5303500634 -0.1092267410  0.1569093937
##  [826]  1.9098790448  0.6574167626  0.6137897063  1.1098353272  0.4886046533
##  [831]  0.1949394750 -0.5011744280  1.7086379145  0.6090731282  1.5731970254
##  [836]  1.2567395236 -0.3485830794  0.6092817147  0.4927024109  0.5794462958
##  [841]  0.7415788551  1.5891007742  0.0335531442  0.7907106033  1.3835977853
##  [846]  1.1686560108  0.6784582697  0.6979466514  1.2277976184  1.1068507134
##  [851]  0.7200619311  1.4227520071  0.6640875629  2.1530009895  0.4397093309
##  [856]  0.5458311108  0.1832482005  0.6565738745  2.0797437018  0.9138673913
##  [861]  0.4074590545  0.2153912431  0.6641538992 -0.1817055060  1.1094777368
##  [866] -0.2830760572  0.9012405224  0.9385369827  0.6952165843  0.4732734983
##  [871]  1.1134577018  1.6544295955  0.6219317650  1.0159486515  0.7358249735
##  [876]  1.1649925101 -0.4936479681  0.7090660737  0.8086302280  0.4886046533
##  [881]  0.9471730356  0.6842045569  1.4514468855  0.2588411558  0.5673272631
##  [886]  0.7557202432  1.0440873138  1.6633349280  0.2877308172  0.6265508670
##  [891]  0.2018215014  1.7498261657  1.1523340071  0.9883111924  1.3262421681
##  [896]  1.0459386147  1.9762672915  0.2034494547  0.1627590512  1.4316047822
##  [901]  0.5261797329  0.4668688256  1.9751810797  0.2831861022  0.9968429158
##  [906] -0.0165630206  1.7214272022  0.9451284852  0.5639446679  0.5721289569
##  [911]  0.5296685644 -0.1945612092 -0.4493385448  1.0532858754  0.8593391671
##  [916]  0.2889147688 -0.3547554008  0.3832563625  0.9829465248  0.7999318368
##  [921]  0.1505490147  1.0033857861  0.1717411384  0.2743813214  0.9590085585
##  [926]  0.4958582052  1.2283556944  0.7236140259  1.0241934570  0.0585000726
##  [931]  1.6097736054  0.0770296263  1.1587738148  0.6090075175  0.5425291129
##  [936]  0.4436010423  0.2382624792  0.6099766208  1.6078658281  0.7007297624
##  [941] -0.1091166683  0.5562272825  0.7238892243  1.9631468249  0.1016388570
##  [946]  2.5003781912  0.7971397914  0.6806354296  0.3576345379 -0.7918890860
##  [951]  1.1754695190  0.3687193998 -0.0814322894  1.1495286318  1.0676194651
##  [956]  0.3720289667  1.3826900780  1.1711840442  0.9918604725  1.2454589224
##  [961]  1.3926628682  0.7082208938  1.7413682814  0.6613479603  1.2192967024
##  [966]  1.0640564502  0.1478135924  0.3543753027 -0.6638980005  0.6934521800
##  [971]  0.4049034662  1.0777682086  1.9278828405  0.3185369518  0.7656188003
##  [976] -0.4425368669  0.9504819228 -0.4575490916  0.6443085015  0.7200619311
##  [981]  1.1083816812  0.6577323489  1.1762995664  0.7181498386 -0.0319910412
##  [986] -0.9408271617  1.4794446292  0.2639965653  0.0009671164  0.0716414614
##  [991]  0.3035245488  1.0480075997  0.7942828525  1.0031281401 -0.2300566175
##  [996]  0.8525458842 -0.6479938372  0.3010171207  0.3210696707  1.0712979623
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

<img src="Week-2-lab_files/figure-html/unnamed-chunk-23-1.png" width="672" />

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
```

```r
fit2
```

```
##       mean          sd    
##   0.50453283   0.35599832 
##  (0.11257655) (0.07960104)
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
## [1] -0.21313461 -0.97462537  0.63142948 -0.10785476  0.16303171 -0.01877426
```

```r
hist(results,breaks=30,col="pink",ylim=c(0,1),freq=F)
hist(results.norm,breaks=30,col="lightgreen",freq=F,add=T)
hist(results2$thetastar,breaks=30,border="purple",add=T,density=20,col="purple",freq=F)
```

<img src="Week-2-lab_files/figure-html/unnamed-chunk-24-1.png" width="672" />

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
## [1] -0.0283
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8778039
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
## t1*      4.5 0.01321321   0.9177964
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 2 3 4 5 7 8 9 
## 1 2 2 1 2 1 1
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
## [1] -0.011
```

```r
se.boot
```

```
## [1] 0.946845
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

