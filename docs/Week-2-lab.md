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
## 0 1 3 4 5 7 
## 3 1 3 1 1 1
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
## [1] -0.0125
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
## [1] 2.612743
```

```r
UL.boot
```

```
## [1] 6.362257
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##   2.5%  97.5% 
## 2.6975 6.4025
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
##    [1] 3.9 4.2 5.1 4.3 3.8 3.8 3.0 3.7 5.9 5.5 4.1 3.9 6.1 4.9 5.2 3.5 5.5 3.0
##   [19] 5.9 5.0 3.7 2.6 5.1 4.1 6.7 3.6 4.5 5.6 4.2 4.3 3.3 4.3 6.3 5.5 4.0 3.9
##   [37] 5.0 4.8 4.9 4.3 4.0 4.0 6.1 4.2 6.7 4.7 3.1 4.8 4.3 3.6 4.4 3.0 4.0 5.1
##   [55] 4.3 2.9 3.8 5.6 2.6 4.0 3.9 5.0 5.7 5.1 5.8 4.7 4.6 4.8 3.7 4.7 3.3 2.6
##   [73] 4.4 4.8 5.6 3.9 2.5 3.7 5.5 4.1 5.0 4.9 5.4 4.1 5.3 4.1 4.9 3.9 4.0 4.2
##   [91] 5.3 4.4 5.3 4.4 2.9 4.0 3.7 3.0 4.8 4.1 3.8 5.9 4.7 4.0 4.2 6.5 2.4 4.0
##  [109] 4.0 4.7 4.6 3.4 3.1 4.9 6.1 4.8 4.1 5.5 5.7 5.2 5.0 5.0 4.6 3.8 4.1 5.3
##  [127] 5.2 4.1 4.6 3.5 3.5 3.5 5.1 4.4 5.9 5.1 4.3 5.8 4.6 5.2 4.8 2.2 5.3 5.0
##  [145] 4.0 5.1 3.6 3.8 2.2 4.9 4.7 5.1 5.7 5.2 4.4 4.6 3.3 5.2 4.8 4.0 6.2 4.0
##  [163] 3.4 4.9 2.5 3.0 3.5 4.3 4.3 6.2 5.3 5.6 3.8 3.1 5.0 4.1 4.5 3.0 4.9 6.0
##  [181] 3.0 5.5 7.4 4.3 5.4 4.6 4.4 5.3 5.6 3.5 3.8 5.3 4.3 4.9 4.0 4.0 2.8 3.0
##  [199] 4.6 5.4 4.4 5.7 2.8 3.9 5.8 4.2 4.2 6.1 2.5 3.1 4.8 4.2 3.8 4.1 5.1 3.1
##  [217] 4.0 3.2 3.8 4.8 4.3 3.3 4.7 5.1 6.4 5.9 4.2 4.6 5.6 3.3 5.6 4.5 3.6 5.0
##  [235] 4.2 5.3 4.1 5.6 4.1 5.5 4.0 4.0 3.1 4.3 5.4 4.2 2.7 5.3 4.1 4.1 6.5 5.4
##  [253] 4.2 5.4 5.7 5.3 3.8 2.4 3.8 4.4 3.7 4.5 4.1 4.9 3.3 4.7 3.3 4.1 5.0 5.5
##  [271] 4.9 3.4 3.7 3.0 3.9 5.9 3.4 4.4 3.9 5.2 4.0 3.4 3.9 3.9 4.9 6.4 5.3 3.6
##  [289] 4.4 5.0 3.2 4.6 5.3 4.7 3.7 3.2 3.9 4.4 5.5 3.9 4.1 4.8 4.9 4.1 6.4 4.6
##  [307] 6.0 4.4 4.2 3.4 5.9 3.7 6.0 4.8 4.2 5.2 5.1 3.4 3.7 3.2 5.7 3.2 4.3 6.2
##  [325] 5.4 5.6 2.8 5.7 3.5 3.3 4.1 3.9 4.3 4.2 3.6 3.8 3.6 5.7 3.9 4.6 4.8 4.2
##  [343] 4.7 5.5 5.2 3.8 4.4 4.8 3.4 3.5 5.1 2.7 4.1 3.9 3.9 3.6 4.4 3.7 4.5 4.3
##  [361] 3.7 5.8 4.6 5.2 5.6 4.4 5.1 4.5 6.4 4.9 4.4 4.6 3.7 3.4 4.1 4.7 4.3 3.8
##  [379] 5.4 5.7 4.2 4.9 6.1 5.8 5.4 5.8 4.6 4.4 6.6 4.1 6.5 3.8 7.0 3.8 5.3 4.3
##  [397] 5.6 4.1 5.5 3.9 5.3 6.3 4.9 4.1 4.6 6.0 4.5 4.8 4.7 5.4 5.7 5.0 4.5 5.8
##  [415] 3.9 4.9 5.2 3.1 4.8 5.4 5.7 3.5 4.5 4.6 4.6 5.8 3.9 4.4 4.7 4.5 3.7 5.1
##  [433] 4.1 5.0 5.6 4.7 5.0 5.1 3.3 4.5 5.9 5.2 6.1 3.8 6.7 6.2 4.0 4.3 6.4 4.5
##  [451] 5.0 3.4 4.4 3.3 3.5 2.9 3.3 5.0 4.3 5.5 2.4 5.6 4.8 4.3 4.5 3.7 4.2 4.0
##  [469] 3.8 5.1 4.0 4.5 4.5 3.7 3.5 5.2 4.5 5.2 5.2 5.4 4.0 5.0 4.3 5.2 5.8 3.5
##  [487] 5.7 5.6 5.1 4.2 4.8 3.8 5.2 6.1 3.4 3.5 5.4 6.0 4.0 3.7 3.2 4.2 4.4 4.0
##  [505] 4.4 2.9 3.7 3.1 5.9 4.8 6.8 4.7 3.8 3.0 3.6 5.9 2.1 4.9 4.6 3.5 3.0 3.4
##  [523] 3.3 4.9 4.3 4.2 4.8 5.3 4.6 4.9 4.3 2.7 2.9 6.3 4.6 4.5 5.4 3.7 4.3 3.4
##  [541] 4.6 3.8 4.2 3.5 4.4 3.9 3.0 4.2 3.8 3.6 5.0 4.3 5.3 4.3 4.3 3.9 4.7 6.0
##  [559] 1.8 3.5 4.4 4.9 4.3 4.6 5.5 5.0 5.0 5.0 4.2 4.7 5.5 3.2 3.1 2.7 3.4 3.9
##  [577] 3.0 4.0 3.0 4.3 5.5 3.6 5.9 6.8 5.1 2.7 5.9 4.5 5.6 4.4 3.9 2.7 4.7 2.2
##  [595] 4.2 4.8 4.2 4.1 4.9 3.7 4.0 4.2 5.2 6.0 4.0 5.2 4.2 4.8 4.4 5.1 4.6 5.0
##  [613] 4.1 4.4 5.2 3.2 5.7 2.7 5.2 4.4 5.0 5.0 3.6 4.6 4.6 4.7 3.6 4.5 4.4 3.6
##  [631] 2.7 2.9 2.4 5.4 3.3 4.4 3.4 3.5 5.4 3.9 3.9 3.4 4.3 5.1 6.5 2.8 3.7 5.1
##  [649] 5.1 6.2 4.8 4.0 5.1 4.8 3.9 5.2 3.7 3.1 5.7 2.8 5.0 5.7 3.8 5.4 3.3 6.0
##  [667] 3.8 4.9 2.5 4.1 3.6 6.5 5.0 5.7 2.9 3.4 4.1 5.7 6.2 4.5 4.4 5.0 3.3 5.5
##  [685] 5.8 2.7 3.9 4.9 4.6 5.7 4.5 4.3 4.8 4.4 5.0 5.1 5.5 4.5 3.9 5.9 3.6 4.2
##  [703] 4.5 3.8 3.7 4.0 4.4 4.4 3.9 3.9 4.8 5.2 5.8 4.3 4.0 4.3 5.5 4.5 3.5 4.3
##  [721] 5.2 4.6 4.6 4.1 6.3 4.4 2.2 4.6 4.4 4.6 4.4 4.3 4.7 5.1 4.7 6.2 4.8 4.6
##  [739] 3.6 4.1 4.7 5.1 5.8 4.0 3.8 5.0 4.0 4.7 5.3 5.7 4.3 3.9 2.8 3.3 5.9 4.3
##  [757] 1.9 6.0 3.1 5.2 3.9 4.4 4.5 5.4 4.7 4.5 7.6 4.5 4.1 4.1 4.7 4.0 5.2 3.6
##  [775] 5.8 3.2 3.7 4.5 4.9 2.5 5.2 5.5 4.1 4.4 3.0 5.0 3.6 5.8 4.7 6.0 4.4 3.7
##  [793] 4.5 5.1 2.2 5.3 6.5 4.6 4.3 3.7 5.5 4.3 4.1 2.7 6.0 4.0 3.3 5.9 3.8 4.0
##  [811] 4.9 2.6 5.5 6.2 4.8 5.1 5.6 4.6 6.1 4.8 4.9 4.5 5.2 5.4 5.2 5.8 4.6 3.7
##  [829] 5.7 4.9 5.3 4.5 5.7 4.3 4.3 6.1 3.5 5.4 4.0 4.5 5.0 4.3 4.6 3.7 3.7 5.3
##  [847] 4.1 7.0 5.5 4.3 2.0 5.1 4.6 3.7 6.4 3.9 4.3 2.1 4.4 6.5 3.4 2.6 1.7 4.6
##  [865] 5.1 4.3 3.2 5.2 3.2 4.7 5.2 5.9 2.7 5.4 5.7 5.7 4.8 5.4 3.4 4.1 5.2 3.8
##  [883] 5.2 4.3 5.0 5.0 3.8 3.5 5.3 3.7 3.5 4.8 4.3 4.6 5.2 3.9 3.9 4.1 3.5 5.6
##  [901] 4.8 4.3 3.7 5.1 3.2 4.7 6.0 5.2 4.8 5.1 4.1 5.5 4.1 4.2 2.0 4.9 4.8 4.6
##  [919] 4.5 2.7 4.5 5.6 4.6 3.1 4.1 4.3 3.8 5.0 4.6 4.0 4.3 4.3 5.3 6.2 4.6 4.0
##  [937] 4.3 4.6 4.6 4.6 5.8 5.4 6.4 3.8 4.3 4.2 3.5 6.4 5.0 3.4 4.2 4.0 4.9 4.6
##  [955] 3.3 4.1 4.2 2.8 4.4 5.2 5.8 5.8 4.8 6.0 3.6 3.0 4.6 4.9 5.9 4.4 3.2 5.3
##  [973] 3.3 4.2 4.3 4.1 5.3 3.7 2.9 4.7 4.7 5.2 5.5 6.0 3.7 4.6 5.2 4.1 3.0 4.9
##  [991] 3.0 3.1 4.9 4.5 6.5 3.1 4.9 3.3 3.1 3.5
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
##   2.6   6.4
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
##    [1] 4.8 5.4 3.6 5.1 5.2 3.5 4.2 4.9 3.2 5.8 4.7 2.1 4.8 5.6 3.8 5.4 3.7 5.0
##   [19] 6.0 4.1 4.3 4.0 4.6 4.6 3.0 5.1 6.2 3.4 6.0 6.3 5.2 4.2 3.9 4.7 5.8 4.6
##   [37] 3.9 5.4 4.1 3.7 5.1 5.7 4.6 4.1 4.0 4.4 5.7 4.4 4.6 6.6 2.9 5.3 5.1 3.1
##   [55] 3.5 4.7 3.8 4.3 5.2 5.5 6.7 2.9 4.4 6.0 5.7 6.4 4.1 3.6 3.3 7.2 5.7 2.6
##   [73] 5.4 3.9 5.0 5.2 4.1 5.3 3.6 5.7 5.0 4.2 5.7 5.0 4.6 5.2 3.6 4.3 4.8 6.0
##   [91] 4.6 3.4 5.1 5.9 4.2 5.4 3.8 4.2 3.4 4.5 3.3 4.1 4.7 5.4 4.6 2.5 4.6 3.8
##  [109] 4.2 5.3 3.0 5.3 5.7 3.9 4.4 4.0 5.1 3.5 4.3 5.3 5.0 5.5 6.2 3.2 4.4 3.3
##  [127] 3.6 4.7 4.0 4.3 6.0 4.7 6.0 5.3 5.1 4.0 5.6 5.9 4.9 4.2 5.6 4.3 5.4 4.5
##  [145] 3.3 3.5 5.3 4.9 5.1 6.4 4.3 4.7 5.6 5.9 4.7 4.7 5.4 4.2 6.0 5.0 3.1 5.2
##  [163] 6.1 4.7 5.7 3.0 4.3 5.8 4.2 5.5 3.2 3.8 4.6 5.0 6.0 5.2 4.1 4.7 3.2 4.1
##  [181] 4.6 3.9 4.0 4.0 3.0 5.2 4.7 2.8 5.1 4.1 6.0 4.4 4.1 4.2 6.0 3.6 5.0 2.7
##  [199] 4.7 4.6 5.3 4.4 4.4 3.2 3.3 5.0 4.4 4.8 4.7 4.1 5.0 5.8 5.4 6.1 3.9 5.1
##  [217] 4.0 5.2 3.5 5.7 4.6 4.6 5.5 2.7 4.0 5.3 3.6 2.9 3.3 5.0 6.5 4.3 3.1 4.7
##  [235] 2.9 6.4 5.6 3.5 5.0 4.8 4.7 4.9 5.3 3.4 3.5 3.3 4.8 5.9 3.4 4.1 4.5 4.2
##  [253] 3.0 5.5 5.1 4.4 4.4 4.0 5.3 3.8 6.6 4.2 3.1 3.9 5.0 3.5 3.8 4.3 6.5 3.9
##  [271] 4.0 4.6 5.1 3.6 5.4 4.2 5.1 4.9 3.8 4.2 3.4 3.0 3.8 5.1 5.2 5.6 4.7 4.9
##  [289] 4.4 5.9 5.1 4.3 4.2 3.4 4.2 4.9 4.1 4.5 5.2 3.5 5.4 5.0 4.1 4.9 4.7 4.9
##  [307] 4.8 5.2 2.7 2.7 3.9 5.9 5.0 5.3 4.1 3.3 2.7 5.8 6.3 4.9 3.8 4.7 4.2 5.5
##  [325] 5.1 5.7 5.7 3.6 4.9 5.6 3.5 2.8 3.9 3.3 5.8 4.8 5.0 4.5 5.0 2.9 4.6 4.0
##  [343] 4.9 3.7 2.8 3.1 6.0 5.4 4.2 4.4 3.1 4.8 5.8 3.9 6.0 3.2 5.2 3.1 5.1 4.8
##  [361] 6.6 5.1 4.4 4.5 4.3 4.9 4.9 6.0 5.0 4.3 5.3 5.8 6.3 5.0 4.3 5.5 4.8 5.0
##  [379] 3.0 2.7 5.5 4.0 4.7 4.3 4.4 6.1 3.7 5.4 4.7 2.5 4.2 5.8 5.2 5.4 5.0 4.6
##  [397] 4.6 3.8 4.0 4.4 5.5 3.2 2.6 4.5 3.4 5.5 4.9 4.2 5.4 3.8 4.9 6.0 4.6 2.8
##  [415] 3.2 4.0 4.9 5.1 5.0 4.6 4.0 5.5 2.5 3.5 3.5 4.0 4.9 5.1 4.9 6.2 5.9 3.1
##  [433] 4.0 6.1 4.7 4.3 5.2 3.7 5.3 5.0 5.3 4.3 4.8 4.5 3.9 4.6 4.8 5.0 4.1 4.4
##  [451] 3.2 4.7 4.2 3.6 4.7 4.9 5.4 4.9 4.5 3.7 3.8 3.7 4.6 5.2 5.1 4.8 4.7 3.2
##  [469] 2.8 2.9 3.6 3.3 5.0 4.9 5.3 3.2 5.8 4.7 4.5 4.3 5.5 4.6 4.1 2.4 4.2 3.8
##  [487] 4.5 3.8 3.8 4.2 4.8 5.5 3.6 5.9 5.8 7.6 3.9 2.9 4.7 6.0 4.7 5.0 5.0 5.4
##  [505] 3.8 5.2 4.3 5.3 3.9 3.1 5.3 5.6 4.3 3.7 3.0 5.4 5.1 4.5 5.2 4.1 5.2 4.3
##  [523] 5.5 3.4 6.5 2.1 4.3 6.3 5.6 5.3 5.5 4.7 3.5 5.7 3.6 4.5 3.9 3.6 3.6 3.0
##  [541] 6.2 5.3 4.1 5.0 2.8 4.4 3.7 3.7 4.8 3.9 4.7 5.5 5.0 3.7 4.5 3.5 4.9 5.4
##  [559] 4.7 5.6 3.9 5.2 4.6 4.5 4.6 4.9 4.6 4.7 4.9 5.3 3.7 3.8 4.7 4.7 5.7 5.0
##  [577] 4.6 4.3 3.7 4.1 3.9 3.5 4.5 3.6 3.9 3.9 3.9 6.2 2.8 3.7 3.6 5.9 4.5 5.6
##  [595] 5.3 3.9 5.8 5.9 5.8 4.0 4.6 5.0 4.4 3.1 5.2 4.8 4.5 3.8 5.7 3.9 3.8 6.8
##  [613] 4.3 5.4 4.8 5.4 4.0 4.6 3.9 5.2 5.0 4.3 5.6 3.8 5.7 4.1 4.8 4.6 4.2 3.8
##  [631] 5.2 5.6 4.8 5.1 4.4 2.3 5.5 5.9 7.2 5.2 4.0 5.0 3.8 5.3 2.8 6.4 5.8 3.5
##  [649] 5.6 3.6 6.1 3.0 6.0 4.2 2.7 5.1 3.7 3.7 4.5 3.1 3.6 3.9 4.9 3.4 4.6 5.2
##  [667] 4.8 5.2 3.5 4.7 4.4 3.4 4.9 4.8 3.5 4.8 6.7 5.2 5.9 5.2 5.4 3.8 3.6 3.7
##  [685] 4.4 4.9 4.5 4.6 4.3 4.6 4.7 5.1 3.6 2.6 4.7 5.4 3.3 4.5 3.6 3.4 5.5 5.6
##  [703] 4.3 3.5 4.8 5.8 3.4 5.4 4.9 5.0 3.0 3.9 5.4 5.0 4.3 6.4 6.0 5.8 5.6 5.7
##  [721] 3.1 4.8 4.5 4.7 2.8 3.0 5.5 3.4 3.9 5.8 4.0 4.2 5.4 3.6 4.9 3.1 3.7 4.7
##  [739] 3.9 3.0 3.5 4.8 3.1 4.8 4.3 5.0 3.7 4.7 5.3 3.1 4.9 5.0 4.8 5.2 5.0 3.3
##  [757] 4.2 3.9 3.2 5.4 5.2 5.5 3.7 3.2 3.8 4.8 2.5 4.4 5.1 5.2 5.2 5.5 5.3 2.9
##  [775] 3.7 5.8 4.8 5.7 3.1 3.9 6.4 5.1 4.2 5.3 4.3 5.4 4.3 4.7 2.7 5.8 4.8 5.0
##  [793] 5.7 3.5 4.1 6.6 3.0 4.6 4.5 4.2 4.6 4.2 4.8 5.1 5.5 6.1 5.4 3.9 6.4 4.0
##  [811] 6.1 5.3 3.3 3.7 5.2 3.5 5.2 4.9 3.6 5.6 4.6 4.7 5.2 3.4 4.4 6.6 3.8 6.2
##  [829] 5.1 5.6 5.7 4.9 5.8 5.9 5.4 5.0 5.2 4.6 3.7 5.1 4.5 5.2 4.2 4.9 6.3 3.3
##  [847] 3.6 4.5 3.4 4.1 4.3 4.9 3.6 3.3 5.5 4.4 4.7 5.3 5.5 3.4 4.2 3.1 4.0 3.1
##  [865] 4.5 3.8 5.1 4.2 4.3 5.6 5.6 5.4 4.2 3.9 4.1 5.4 5.4 4.0 2.5 3.6 4.4 4.8
##  [883] 4.6 3.6 4.6 3.9 5.2 6.2 4.3 4.0 4.6 5.2 4.2 3.7 4.8 5.2 4.2 5.3 4.8 5.0
##  [901] 3.7 4.4 4.5 4.7 4.7 3.8 3.6 5.3 2.0 4.6 6.4 3.7 3.5 4.5 5.7 5.9 4.1 5.4
##  [919] 5.4 5.4 4.9 4.0 3.8 4.9 3.9 3.8 5.0 3.0 3.6 4.2 4.8 3.9 3.0 3.8 5.6 4.0
##  [937] 5.2 2.9 5.4 4.8 3.2 2.4 5.0 2.6 4.1 3.7 4.6 3.7 4.0 4.0 6.0 4.1 4.2 4.6
##  [955] 4.8 6.4 4.6 6.0 4.6 6.3 5.1 6.2 5.2 4.7 4.2 5.2 5.2 3.9 4.3 4.9 5.3 5.3
##  [973] 5.5 5.1 4.4 5.3 3.2 5.1 2.9 3.6 4.8 5.0 5.6 4.2 4.9 4.2 3.3 5.0 4.5 5.6
##  [991] 4.0 3.3 6.2 4.7 3.4 5.0 4.4 5.2 4.5 3.4
## 
## $func.thetastar
## [1] 0.0582
## 
## $jack.boot.val
##  [1]  0.62338028  0.49570201  0.31325967  0.24722222  0.18347107 -0.03422619
##  [7] -0.15312500 -0.28357771 -0.35347432 -0.42148997
## 
## $jack.boot.se
## [1] 1.037836
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
##    [1] 4.2 5.3 4.3 4.9 3.3 3.6 4.2 4.8 4.0 4.5 4.2 3.9 4.2 3.5 4.8 6.0 5.8 4.1
##   [19] 3.6 3.3 3.9 4.2 4.2 4.1 5.5 6.0 2.7 4.4 5.6 4.8 4.1 3.0 2.4 5.0 3.9 4.5
##   [37] 5.0 5.7 3.2 3.8 3.7 3.8 4.8 3.5 3.9 4.7 4.5 4.4 5.7 3.0 4.9 5.1 4.8 4.8
##   [55] 6.2 5.3 5.8 3.9 4.6 3.7 3.5 5.0 4.9 2.9 5.0 5.7 5.1 4.2 3.1 3.7 4.2 4.5
##   [73] 4.6 3.1 4.6 4.7 4.8 4.3 5.0 4.7 4.3 5.0 6.0 5.6 4.5 4.8 3.3 5.6 3.7 6.5
##   [91] 4.2 3.5 3.6 4.8 5.2 3.9 3.9 4.8 5.5 4.5 5.2 3.9 3.9 4.9 3.7 5.2 5.1 4.1
##  [109] 5.4 5.9 5.3 4.0 4.6 5.3 4.0 4.8 3.2 4.5 3.7 5.1 4.8 5.6 3.3 5.8 6.0 6.4
##  [127] 3.8 5.9 3.3 5.9 5.3 5.0 5.3 3.2 5.6 4.5 4.9 3.9 4.6 5.2 4.4 3.6 5.1 4.6
##  [145] 5.8 5.8 5.1 5.2 5.5 5.0 5.0 3.5 4.9 5.4 3.4 3.1 4.1 4.7 5.0 3.8 4.1 5.8
##  [163] 5.4 4.6 4.8 4.0 3.7 4.0 4.6 3.2 4.4 2.7 4.4 3.6 4.9 4.0 5.2 4.0 3.7 4.7
##  [181] 4.8 4.0 4.4 5.3 3.5 3.9 5.4 4.3 3.9 4.5 3.7 5.7 3.2 4.1 3.6 4.9 5.3 4.7
##  [199] 4.7 4.2 3.0 3.3 4.2 4.0 4.3 5.5 5.4 6.4 4.5 3.7 4.2 5.9 4.7 4.5 4.3 3.4
##  [217] 5.9 5.2 4.3 5.4 7.3 3.4 4.9 4.8 4.3 3.3 5.7 4.2 5.4 3.8 4.1 5.9 4.9 3.7
##  [235] 4.9 4.1 4.2 3.6 4.9 4.8 4.8 3.8 6.0 3.6 6.1 4.2 3.6 5.6 4.6 5.5 5.3 4.5
##  [253] 5.5 2.8 6.1 6.5 4.4 3.1 3.8 2.7 4.7 4.9 4.7 3.8 4.0 3.4 4.2 5.2 5.0 3.2
##  [271] 3.2 5.3 5.0 5.6 4.7 4.6 4.3 4.2 4.1 3.4 3.8 4.6 5.0 5.3 5.8 5.0 4.7 4.8
##  [289] 5.4 2.6 5.3 5.2 3.7 5.5 3.4 3.9 3.7 4.3 4.9 6.1 5.8 4.4 3.9 3.8 3.9 5.5
##  [307] 2.6 5.3 5.7 4.9 5.2 4.5 5.5 4.0 4.9 6.0 4.6 4.5 4.9 4.1 4.6 6.2 2.7 3.7
##  [325] 3.9 4.6 4.2 4.6 6.9 4.1 4.4 4.5 5.0 5.1 5.9 4.7 4.0 4.6 5.8 4.9 4.0 4.7
##  [343] 3.3 4.1 5.9 4.1 3.2 3.7 4.5 4.0 2.9 5.4 4.7 3.2 5.6 5.8 4.9 3.6 3.8 5.7
##  [361] 6.0 4.9 5.7 5.5 4.6 5.4 5.8 3.9 3.8 4.3 4.6 4.8 5.9 4.4 4.5 5.4 3.0 5.9
##  [379] 3.9 4.6 6.3 5.4 4.3 3.8 5.3 4.5 5.1 5.7 4.7 4.8 3.6 3.8 4.4 4.4 3.3 3.2
##  [397] 3.6 4.1 4.2 3.9 5.5 4.8 4.5 4.1 4.0 4.4 6.1 4.4 4.5 4.1 3.4 6.1 3.9 4.2
##  [415] 4.3 4.2 2.7 3.4 5.4 4.2 4.8 3.3 3.6 5.2 5.4 3.9 5.0 5.5 4.6 5.4 3.2 3.4
##  [433] 4.7 3.4 3.3 3.9 6.3 4.2 5.8 4.1 3.7 4.3 2.9 4.9 4.9 4.2 3.7 3.5 4.4 4.8
##  [451] 5.0 4.6 4.9 3.9 4.4 2.9 4.3 4.3 6.1 3.6 7.0 4.6 4.5 5.2 6.0 2.6 5.2 5.8
##  [469] 3.9 2.3 4.8 3.5 4.2 4.2 5.1 4.5 4.8 5.2 4.5 3.5 3.4 4.7 4.6 4.2 3.7 5.6
##  [487] 4.5 3.3 4.2 5.3 3.4 3.5 4.1 3.7 6.2 4.5 4.7 4.5 5.4 3.6 3.7 3.7 4.3 4.4
##  [505] 4.9 3.5 3.8 5.0 4.9 4.9 3.0 4.8 4.5 4.2 2.1 4.0 4.5 4.2 4.8 5.1 4.0 4.1
##  [523] 2.2 5.0 5.1 3.6 4.0 3.9 4.7 6.2 6.3 5.7 5.2 4.2 3.9 5.4 2.7 5.4 3.8 4.2
##  [541] 5.2 3.3 4.8 5.5 4.0 4.6 3.8 4.6 5.9 6.6 5.1 4.4 4.2 3.7 4.9 2.9 4.3 5.0
##  [559] 3.8 3.7 2.9 3.9 4.2 5.4 4.1 4.4 4.5 5.8 3.3 5.8 4.2 5.6 4.3 6.6 4.9 4.3
##  [577] 3.8 4.0 4.9 4.9 5.8 6.0 4.9 4.9 3.4 5.3 2.0 5.8 4.7 4.3 4.9 4.4 3.8 5.5
##  [595] 5.4 3.9 4.9 3.8 4.4 4.9 5.7 3.8 4.4 4.9 5.9 5.1 6.1 5.2 2.5 3.5 4.0 4.6
##  [613] 3.0 4.1 3.8 4.9 5.2 4.2 2.2 6.3 4.5 4.4 5.9 4.0 5.1 4.2 4.2 4.5 3.9 4.8
##  [631] 4.1 4.7 3.2 4.2 4.4 5.4 4.0 4.9 4.4 5.7 4.5 6.0 6.4 3.7 5.5 5.0 3.0 3.8
##  [649] 3.8 4.1 6.0 4.7 4.9 3.7 5.1 3.5 5.5 6.7 3.1 4.4 4.9 6.2 4.3 4.2 4.3 5.3
##  [667] 3.8 3.1 4.5 6.3 4.1 4.2 4.2 5.2 4.2 3.6 4.4 5.2 4.5 6.0 5.0 5.4 3.3 4.8
##  [685] 3.3 5.7 5.3 5.8 3.7 3.4 4.8 3.2 3.6 4.6 3.6 3.9 5.6 3.3 4.2 4.4 5.8 4.7
##  [703] 3.7 5.5 4.7 4.7 4.7 5.0 2.4 3.2 4.8 4.9 4.9 5.5 3.5 3.0 4.0 5.5 5.6 5.5
##  [721] 3.9 5.1 5.2 6.2 3.7 4.2 4.2 3.8 3.5 3.9 3.9 4.7 4.6 5.8 5.3 5.1 5.0 4.5
##  [739] 5.9 4.3 5.8 4.5 4.6 3.4 3.5 5.0 4.2 5.0 4.0 5.6 5.8 6.4 4.7 5.5 4.8 4.7
##  [757] 4.9 5.2 3.9 4.0 4.1 4.6 4.9 3.8 4.0 5.0 5.0 4.6 5.3 3.7 3.9 4.3 4.4 5.1
##  [775] 4.8 4.7 4.3 3.1 3.8 4.8 4.6 4.1 6.0 5.2 4.3 4.0 4.0 3.4 4.8 2.2 5.2 4.6
##  [793] 5.4 5.1 2.8 4.0 4.2 4.1 5.9 5.1 3.8 4.0 5.6 4.6 4.2 4.3 5.2 5.4 4.6 6.2
##  [811] 3.7 5.9 4.8 6.9 2.4 4.1 4.3 5.4 5.0 4.3 4.0 6.0 4.1 3.6 3.6 4.0 4.6 5.2
##  [829] 4.6 4.5 6.7 5.1 4.1 5.4 4.3 3.5 6.5 5.2 4.9 2.8 4.0 5.6 4.0 3.5 3.9 2.4
##  [847] 5.3 4.6 5.3 4.4 5.7 4.7 6.1 4.4 5.1 3.2 3.8 4.4 4.8 2.6 2.7 5.4 4.7 4.7
##  [865] 3.9 3.5 5.3 6.4 3.7 2.9 4.5 3.9 4.4 3.1 4.9 5.0 3.0 5.2 4.5 5.2 5.2 3.8
##  [883] 4.6 5.0 3.1 2.9 4.7 6.5 4.1 6.2 4.2 4.9 4.8 3.2 4.0 3.7 5.8 4.5 4.1 2.7
##  [901] 4.4 4.7 5.0 4.9 4.0 4.6 5.4 4.1 4.5 2.6 3.9 4.4 3.7 4.9 4.9 4.5 4.8 3.7
##  [919] 5.9 4.8 7.1 5.0 5.1 4.1 5.6 5.5 5.1 3.3 3.9 3.1 3.6 2.8 4.7 3.6 4.5 4.7
##  [937] 3.7 4.9 2.2 4.3 3.6 3.3 5.2 2.4 3.0 4.2 4.1 5.8 3.9 5.3 2.7 4.3 6.3 5.7
##  [955] 4.5 5.5 5.2 4.6 4.6 5.3 3.9 3.5 4.8 4.0 3.9 4.2 5.2 4.4 5.8 5.5 4.8 3.5
##  [973] 6.3 5.5 5.7 6.4 4.0 3.1 5.0 4.5 3.8 4.1 5.7 4.2 3.1 3.5 4.0 4.1 5.5 4.0
##  [991] 5.3 4.8 4.2 2.9 3.6 3.6 3.0 4.2 3.8 3.3
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.5 5.4 5.4 5.2 5.2 5.0 4.9 4.7 4.6 4.5
## 
## $jack.boot.se
## [1] 1.014692
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
## [1] 0.6864247
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
##   3.746750   5.121765 
##  (1.606584) (2.350182)
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
## [1] -0.791516750  0.003366605  0.340187606 -0.126075003  1.103621580
## [6]  0.542936105
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
##    [1]  0.8173727310 -0.1218063703 -0.2183310142  1.2071100200 -0.0083868488
##    [6]  0.4748322890  0.1953249194  0.8000264812  0.7906628657  0.8033131313
##   [11]  0.4697232627  0.0255249424 -0.0989853938  0.4390106619  0.9225653727
##   [16]  0.3865872163  0.2902061874 -0.1038715645  1.4015669443  0.3142406754
##   [21]  0.4517173355  0.0586229436  0.1578307924  0.2392126605  0.4476097947
##   [26] -0.1084062870  0.6841378120  0.7266403819  0.3099392991  0.3968988312
##   [31]  1.1408928953 -0.0412662406  0.4527437923  0.3882685776 -0.9802287510
##   [36]  0.2712248344 -0.1188604349  0.0179936585  0.1513633364  0.6406437459
##   [41]  0.0365168652  0.0383988730  1.0341768555  0.5618657447 -0.2744132450
##   [46]  0.3163605697  0.3422153832  0.8766216453 -0.0946120171  0.6740362371
##   [51]  0.7901385557 -0.0923728713  0.1931095270  0.5646944212  0.5541841706
##   [56]  0.2338936957  0.6573795111  0.9868984782  0.4957690491  0.8105008159
##   [61]  0.1268435270  1.0757763446  0.6270794223  0.6220930212  0.3322071322
##   [66]  0.1960848980  0.3454757840  0.3772467565  0.8297422565  0.2537233926
##   [71]  0.7749591732 -0.4789943255  0.1247552049 -0.3147155027 -0.5037935956
##   [76] -0.0597517474  0.9386581536  1.5472973268  0.6817840350  0.5576673732
##   [81]  0.2082048345  0.5051771120  0.6353045813  0.3344962962  0.7502344008
##   [86]  0.5893365891  0.0475511886  0.2857517179  0.1059164570  0.8411986634
##   [91]  0.0329439999  0.7592815006  0.6635271691  0.4377120802  0.8330463518
##   [96]  0.7186418436  0.5748062751  0.7322823721  0.1726708035 -0.2667680956
##  [101]  0.4288474632  0.8974734349  0.5261871870 -0.0878668836  1.0682768755
##  [106]  0.3556377231  0.6061978372  0.3801863304  0.1687141182 -0.1014172857
##  [111]  0.4498025769  0.0658165821  0.4337317691 -0.1663222491  1.1171949646
##  [116]  1.0202089502  0.9965211278  0.2527887557  0.2364657441  0.1602849565
##  [121]  1.0406719590  0.6699828065 -0.0609970437  0.4625894043 -0.3987053957
##  [126]  0.1326387005  0.0991761053  1.0857260234  0.6735851734  0.8894683739
##  [131]  0.6716698374  1.0253107835  0.5527354384  0.1461515875 -0.0639184439
##  [136]  0.3121016545  1.3451189868  1.0616456575  0.1390847467  0.5253125571
##  [141]  1.4525126594  0.5247126018  1.3808906956 -0.2020988706  0.7756209204
##  [146]  0.8164433951 -0.2055196319  0.9732187972  0.4943040435  0.5467986131
##  [151]  1.2027033827  0.8970794174 -0.5869086894  0.5180951197  0.2837587710
##  [156]  1.3398355946  0.1177517614  0.4218608963  0.3279027203  0.9432667777
##  [161]  0.4879791966  0.6525465844 -0.1835429320  0.4787474509  0.1747584641
##  [166]  1.0718424773  0.8500403432  0.9802599048  0.6250707842  0.5452717350
##  [171]  1.3307491067  1.5712441285  0.2832561397  0.3019071628  0.0893748850
##  [176]  0.1990587755  0.7728911566  0.8887505026 -0.6297769461  0.2105392958
##  [181]  0.3737024850  0.2442447921  0.3420414985  0.2019666291  0.0880950821
##  [186]  0.3610964147  0.1560118666  0.0010584035  0.3834625871  0.1870660833
##  [191]  1.2156773453  0.9299223639  0.1040489167  0.8329177109  0.4261062418
##  [196] -0.2775158052  1.1515857297  0.4007390745  0.1688146440  0.3011505170
##  [201]  0.3745232830 -0.1428480357 -0.1954879546 -0.3281134846  0.4078856213
##  [206] -0.3723766582 -0.2873551465  1.0279685731  0.6977760299  0.3155209034
##  [211]  0.6822573653  0.8770533161  0.1268435270  0.4932678962  1.0286625853
##  [216]  0.2873750282  1.1491548037  1.3101028936  0.8523168997 -0.8807756052
##  [221]  0.8559505058  0.6989363905  0.5931916295  0.6507974450  0.1464129561
##  [226]  0.4409559849  1.1994427879  1.2516918837  0.0847872184  0.2187220462
##  [231]  0.8201170143  0.8006537761  0.7206297852  0.8720290016  0.7708869261
##  [236] -0.1058341831  1.2139504284  0.1212493410  0.3352295490  0.8722336241
##  [241] -0.0715901857  0.7084891682  0.6913628138  0.5095172168  0.4765107544
##  [246]  1.1910774710  1.6695297087  0.4403049615  0.8001925774 -1.0845599127
##  [251]  1.1103790560  0.0511737996  0.6943447491  0.7580527424  1.4487149873
##  [256]  0.1170444262  0.6824132818 -0.2289065422  0.9380434816  1.1026161842
##  [261]  0.5783693850  0.4082173424  0.3369642953  0.6618672148  0.7442201005
##  [266]  1.0723608186 -0.0391701165  0.3399469661  0.0124395542  0.3055758084
##  [271]  0.5923252390  0.1197129817  0.6111550921  0.7487505452  1.1688522945
##  [276]  0.4389282049  0.0877736801  1.2854384547  0.7187714835 -0.0110689833
##  [281]  0.0904981041  1.3827575905  0.2424103927  0.1687141182 -0.2382807538
##  [286]  0.0214172456 -0.2827728002 -0.0173653830  1.0914728727  0.5595364749
##  [291]  0.2988656490  0.6661986080  0.7361524409  0.2886366779  0.6435108182
##  [296]  0.0658165821  1.1115262376  0.9461694972  0.6162317553  0.1265508904
##  [301] -0.6801646737 -0.7603731386  1.5605329964  0.5506946187  0.9747037271
##  [306]  0.4323026731  0.2396154386  0.8611181213  0.6193371524 -0.0373678538
##  [311]  0.5634178820  0.8694757914  1.0736947053  0.3682213461  0.8287297305
##  [316]  0.0105991315  0.1243225306  0.3890606711  0.8092243570  0.8346082567
##  [321]  1.1256725598 -0.5503712112  1.5196922964  2.0166880135 -0.0550454117
##  [326]  0.7227861377  0.6480562741  0.0940202671  0.8703807227 -0.0616439265
##  [331]  0.3165105331  0.2670982640  1.5471527814 -0.0250369881  0.2421043590
##  [336]  0.2537233926  0.2023540626  0.6824819128 -0.2134857281 -0.0293617543
##  [341] -0.0713588651  0.3898750936  0.3339265376  0.1751614864  0.3707616380
##  [346]  0.8849968830  0.3356238314  0.2886366779 -0.1610568154 -0.1461988137
##  [351]  0.9255623826  0.7880274633 -0.0023351462  0.2598483967 -0.1461988137
##  [356]  0.3103640752  0.7388750789  0.3767424750  0.7005887600 -0.1430299420
##  [361]  1.3731431047  0.2420398035  0.5691508261  0.7997853327 -0.0755053263
##  [366] -0.1683775717 -0.0390576203 -0.2270459691  0.8404183846  0.4921616324
##  [371]  0.2691158472 -0.1778497121  0.3209716358  1.6033908584  0.4142440301
##  [376]  0.8988879005  1.3043047412  0.8297475973  1.5768316027  0.4403707695
##  [381]  0.9587677133  0.1279562528  0.7005175627 -0.0263596513  0.2260660807
##  [386]  0.3557246635  0.7837015094  1.1327483360  0.0007117679  0.3343458679
##  [391]  0.0790039141  1.1570525188  0.6538449725  0.6793653643  0.5103133986
##  [396]  0.4907380001  0.7720214227 -0.4167662364  0.7958800912  1.5178798266
##  [401]  0.7000641316  0.1049348358  1.1515409215  0.1957395603  0.4577543674
##  [406]  1.3057488830  0.8919540241  1.0750431921  0.0350115402 -0.1342991375
##  [411]  0.5657121922  0.5926205773  0.1090536913 -0.5433174799  0.8612363012
##  [416]  0.3285257771  1.3956908175  0.7792716626  0.3283687056  0.3627354088
##  [421]  0.1258366569  0.8571878154  0.4737514004  0.4822838594 -0.1216651353
##  [426]  0.2636010760 -0.0492017569 -0.1609309868  0.4971053727 -0.4192264990
##  [431]  1.6101687746  0.5665181292 -0.0484264193  0.2743352182  0.3372193838
##  [436] -0.0843124517  0.2970463277  0.6987884756  1.1938324561  0.6814671476
##  [441] -0.3860207610  0.5052746586 -0.7158291737  0.6138411481  0.6184101789
##  [446]  0.0158243329  1.0658886353  0.8876438121  0.3898627942  0.7889893309
##  [451]  0.0680770186  0.4266980433  0.3520209481  0.1503269310  0.2611089143
##  [456]  0.4932950295 -0.0806994740  0.3028711039  0.4690222600 -0.0597517474
##  [461] -0.2074409122  0.7266403819 -0.0536866477  0.5388784353  1.1267731688
##  [466]  0.5580161153 -0.2329223678  0.5156200832  0.2307097273  0.5653231678
##  [471]  0.6549050929 -0.0056272859 -0.0514407393 -0.2732383765  0.1983569846
##  [476]  0.0511737996  0.7083132605  1.3245972798  0.3909055795  0.2337299363
##  [481]  0.2594780280 -0.1207834810  0.7825000026  0.1021820554  0.3033469801
##  [486]  0.3037756988 -0.4282419659  0.6394918089  0.5382510626  0.0853432381
##  [491] -0.3342417863 -0.2194238719  0.2085531464 -0.0775615928  1.5004772230
##  [496] -0.0101175126  0.3760862313  1.9783995032  1.2350294159  0.2434590823
##  [501]  0.3888893858  0.7226335913  0.5426651196  0.3784882406  1.9083168525
##  [506]  0.9972609137  0.8159768914  1.7997475798  0.3619027108  0.0874365611
##  [511]  0.1278779369  1.0845706751  0.2357649152  0.6577327785  0.0755031553
##  [516]  0.7402814193 -0.2953925119  0.5089566197 -0.1277015799  1.1200982161
##  [521]  1.0033262193  1.0628011893  0.9298266039  0.6732170860  0.4956137176
##  [526]  0.8898587865  0.5234807436  0.3570068375  0.8428103628 -0.1006980602
##  [531]  0.4234429482  0.6582383401  1.0532533325  0.2302819679  0.2738840122
##  [536]  0.1451050552  0.4379823858  0.9966280278  0.3468044675  1.6802641492
##  [541]  0.5616945864  0.9383403471  0.2389280654  0.5219376019  0.8124123075
##  [546]  0.7963431508  0.7444442669  0.6168077421 -0.2600721122  0.5192955432
##  [551]  0.5006863360  0.8970420451  0.8798636241  0.0665004563  0.6106021349
##  [556]  0.3635518978 -0.0291982780 -0.1796794351 -0.2264119901  1.3935520719
##  [561]  0.7640475570  1.2163637967  0.3319246083  1.0667949071 -0.3405590881
##  [566]  0.8837338603  0.1490755660  0.7118532017  0.6201968538  0.2223772416
##  [571]  0.1181898229 -0.0309844363  0.1354721399 -0.4224529028  0.1069255583
##  [576]  0.2122598516 -0.1749487367  0.9090412520  0.3723203608  0.1523472187
##  [581] -0.2338356299  0.6552368145  0.5679636695  0.2932163497  0.3965303395
##  [586]  0.4792863310  1.0593038239  1.1403458261  1.1050553418  0.6432907261
##  [591]  0.5306496201  0.9614195636  0.7178676321  1.0578623027  0.3783397565
##  [596]  0.2857547814  0.2162472825  0.2128574658  0.7767523888  0.7272623545
##  [601]  0.3739303756  0.3067903649  0.6009703988  1.5731383708  0.5777284565
##  [606]  1.1643194936  1.6223738163  0.5125686509  0.8783009367  0.0847604043
##  [611]  0.0059360918 -0.0578190085  0.6356598472  0.4143958533 -0.0136088139
##  [616]  0.5812755776  0.9250971490  0.8913823412 -0.1557028847  0.0698938093
##  [621] -0.2307929915  0.9729517382  1.3129278552  0.9089570803 -0.0469175146
##  [626]  0.8550272497  0.1337374006  0.7771795694  1.0078643280  0.5855426030
##  [631]  0.8628878969  0.0147784583  0.4663663701  0.5102932445  0.6885249762
##  [636]  0.9435872229  0.7740095113  1.3058010420  1.2690724465  1.0080825237
##  [641]  1.4464396714  0.6391489949  0.1950839084  0.1345237488  1.8080705348
##  [646]  0.8313532089  0.0832214177  0.6146829454  0.8926428218  0.7016861217
##  [651]  0.4069975931  0.4587372703  0.2627903347  1.0312614356  0.2680517125
##  [656]  0.5179167755  0.4595978951  0.2717897564  0.5593628506  0.5864829264
##  [661]  0.7578412990 -0.2160048313  0.5352282038  0.0088055920  0.8466073068
##  [666]  0.2771174443 -0.2314960126 -0.5737147796  0.7119765575  1.0821062049
##  [671]  1.3307491067  0.6159646138  0.7961930254  0.1838563959 -0.2832480011
##  [676]  1.3509473802  0.7206282377  0.5346919390  0.2296396823  0.5542051099
##  [681]  0.4565682962  0.4954681299  0.6212247592  0.4627106351  0.0645786631
##  [686]  0.5049818329 -0.0608903491  0.7686527203  0.5985613500  1.1750557153
##  [691]  0.4233436516  0.8798636241  0.4751400302  0.4694373808  0.1935569224
##  [696]  0.6585197165  0.0594779540  0.4497633905  0.4120794215  1.0220242537
##  [701]  0.8480519343  1.1899188441  1.2186500178  0.9181582946 -0.4931970661
##  [706]  0.2783494522 -0.5488128290  1.7614929847  0.6344092361  0.5999933412
##  [711]  0.0305048588  1.0311137284  0.4505706787  0.7235010182  1.2323958468
##  [716]  0.5207967885  0.8647398725  0.6454151653  0.8574739903  1.2325088672
##  [721]  0.5553037267  0.0243846709  0.5371610025  0.2801173216  0.7542438189
##  [726]  1.2522786693  0.8003506202  1.0578476059  1.1326640777  1.1115262376
##  [731]  1.3667823903  1.0294597304  0.4807225565  1.1491548037  1.0489478592
##  [736]  1.7942643153  0.3628430318  0.9084608763 -0.2677368983 -0.0184280397
##  [741] -0.1137598141  1.2699558897 -0.0874181004 -0.1721736668  0.0642890915
##  [746]  0.7532769355  0.0922986370  0.6277173258  0.3212478450  0.4322595668
##  [751] -0.2654778293  0.7672961174  0.0865653422  0.3947420554  0.5245488235
##  [756]  0.3571713268 -1.1494681634  0.7930617696  1.2258750688 -0.0320072094
##  [761] -0.1099857997  0.1630748403  0.5211302884  1.0201751589  0.9435308609
##  [766]  0.2154394431 -0.0670016912  0.7985583243 -0.0979970172  0.7443681391
##  [771]  1.1251555653  0.3520365028  0.5667990579  1.1265563692  0.8429573885
##  [776]  0.1216479545  1.1468396588  0.7345311693  0.0119675071  0.1221358252
##  [781]  0.1563503022  0.7491628714  0.6486940765  0.8757286965  0.4980710722
##  [786]  1.7938417322 -0.0463083572  1.0365204404  0.7053218955  0.2235043354
##  [791]  1.1730360074  1.3827575905  0.7817598483  0.1861284696  1.1579583125
##  [796]  0.9426575968  0.8989868557 -0.5184549581  0.7736957517  1.4675649016
##  [801]  0.7344282233 -0.0434933980 -0.0004052952 -0.0194383541  1.5401788810
##  [806]  0.6072195508  0.9802824677  0.7431321460  0.1345237488  0.1393549385
##  [811]  1.1393921535 -0.2364251074  0.1944609017  0.7569710570  1.2516918837
##  [816]  0.0292544200  1.4034907648  0.4311979960  0.1717826442  0.3372193838
##  [821]  0.0784209175  0.2257314099  0.2884596014  0.2326385363  0.4020009922
##  [826]  1.0952683049  0.8857652536  0.5922782663 -0.4854392957  1.0369624604
##  [831]  0.0537774874 -0.0302642718  0.7571552367  0.7758574234 -0.2508221000
##  [836]  0.7617488073  0.9686164218  0.4849912886 -0.1739898643  0.8586932142
##  [841]  0.5772460375  1.3117327651  0.8023000306  0.7682735361  2.2016285963
##  [846]  0.7218100674  0.6738721009  1.1172471084  0.4556055316  0.8581546214
##  [851] -0.1016655830  0.4823852179  0.2280775967  0.9734309246  0.4097152020
##  [856]  0.3258701237  1.6590040707  0.4992934050  0.1598352464 -0.1040727690
##  [861]  0.2203310292  0.6338510276  0.6417769913  0.4847057962  0.9948985983
##  [866]  0.6264297999  0.1265007046  0.1364075457  0.2653500978  1.0668180846
##  [871]  0.9715621624  0.1943141841  0.7916773012  1.1132228672  0.4425463112
##  [876] -0.1470285178  0.6660011050 -0.2437139064  0.4319335415  1.1917337617
##  [881]  0.3492241522  0.2608796157  0.7672390678  0.3610784999  0.6891182602
##  [886]  0.6333108459  0.7559554580  1.3670606722  0.3989208947  0.6954866070
##  [891]  0.6786337912  0.1846674011  0.3903386776  0.0023028235  0.2992772920
##  [896]  0.3651875187  1.5687161149  1.5981718935  0.7303920334 -0.2249474259
##  [901] -0.0203880401  0.7416647913  1.4486364790  0.5201927576  0.3799319451
##  [906]  0.4250115426 -0.4087841456  0.5871848149  0.5360572844  0.4400998628
##  [911]  0.7871970433  0.1340669346 -0.0527170768  0.6941562761  0.9530341156
##  [916]  1.3270857832  1.0376151711  0.7956429958  1.2992200069  0.2880662161
##  [921]  0.1409789569  0.0903267441  0.3088076640  0.3951696906  0.0736163751
##  [926]  0.3320627112  0.4910086577 -0.0553991950  0.3224858860  0.3012460301
##  [931]  0.2283520128  0.7577000775  0.2459768136 -0.3748035327  0.6620338298
##  [936]  0.1334375344  0.3380205778  0.4376599583  1.3724882024  0.1355765047
##  [941] -0.0333947128  0.9585317137  1.0306276262 -0.0061037050 -0.1154133195
##  [946]  0.5918392430  0.1402047814  0.7886761453  0.9582751565  1.0696243454
##  [951] -0.0953293741  0.2914059309  0.2398475112  0.6345678325 -0.2896834967
##  [956]  1.0110698196  0.2999869561 -0.1329191028  0.8968751613  1.3927628116
##  [961]  0.5385641800  1.1016680595 -0.3628980126  0.9813183978  0.6528594267
##  [966]  0.4814444437 -0.0770744715  1.2747170313  0.9869504251 -0.2270382410
##  [971]  0.7931582704  1.0970374099 -0.7756234288  0.2547363254  1.0512293152
##  [976]  0.5201927576 -0.1928861023  0.4920096787  1.3161079071  1.1605901351
##  [981]  0.1247552049  0.5889567838  0.5414230425 -0.2659778652  1.1013666384
##  [986]  0.8186513058  0.8598308553  0.5809126990  0.9776774117  0.2837857124
##  [991]  0.7206297852  1.0699472094  0.1741854363  1.0803126998  0.6779086346
##  [996]  1.7497368745  1.4268969694  0.4385637372  0.2781390916  0.6785026972
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
```

```r
fit2
```

```
##       mean          sd    
##   0.73153779   0.37552655 
##  (0.11875192) (0.08396856)
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
## [1] -0.2062151  0.1592828 -0.3238363 -0.7366258 -0.1249817 -0.1602956
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
## [1] -0.0103
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9009082
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
## t1*      4.5 -0.03003003   0.9165037
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 2 4 5 6 9 
## 1 3 1 1 3 1
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
## [1] -0.0343
```

```r
se.boot
```

```
## [1] 0.8663857
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

