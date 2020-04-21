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
## 0 1 3 4 6 8 9 
## 1 3 1 1 2 1 1
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
## [1] -0.0192
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
## [1] 2.655137
```

```r
UL.boot
```

```
## [1] 6.306463
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##   2.5%  97.5% 
## 2.7000 6.3025
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
##    [1] 4.4 4.6 5.9 4.6 5.3 5.3 6.0 3.8 4.5 5.2 4.1 6.0 3.7 3.4 5.7 5.8 5.7 3.3
##   [19] 4.1 4.6 5.0 4.8 3.0 4.3 4.7 5.6 6.0 4.9 4.7 4.0 3.0 5.3 4.6 4.9 4.1 3.0
##   [37] 3.0 5.6 3.6 6.2 5.7 4.0 3.4 4.8 5.0 5.3 5.1 4.6 4.9 3.5 5.0 2.9 4.8 4.4
##   [55] 4.3 3.9 5.9 5.0 4.5 2.2 5.3 4.6 4.6 6.1 4.7 4.1 3.3 5.1 3.4 4.5 4.2 3.5
##   [73] 5.7 5.1 5.0 3.1 4.2 5.7 5.9 5.8 5.7 4.5 3.0 4.1 4.1 2.6 7.0 4.5 5.4 3.1
##   [91] 4.4 3.9 6.4 4.7 3.4 5.4 4.9 5.1 4.4 4.0 4.6 4.9 3.3 3.1 5.5 4.7 6.0 4.7
##  [109] 4.7 4.7 6.2 3.6 4.4 3.6 4.6 3.8 7.6 4.6 5.7 4.4 2.7 3.4 5.0 3.6 3.1 3.5
##  [127] 6.1 4.8 5.5 3.4 5.5 4.9 4.7 4.7 4.3 5.7 2.3 4.1 5.6 4.7 4.7 5.0 4.6 5.9
##  [145] 4.5 5.2 3.6 2.8 4.4 3.8 3.0 6.7 4.4 4.7 4.7 4.5 4.1 4.3 5.5 3.9 5.1 4.3
##  [163] 4.8 6.2 4.9 5.6 3.7 5.6 4.4 5.5 4.6 5.6 5.6 2.5 3.7 4.0 4.9 6.6 4.6 4.0
##  [181] 3.0 4.4 6.1 4.1 4.6 5.2 5.5 4.6 6.7 5.4 3.9 3.3 3.5 4.7 4.8 5.4 4.3 3.1
##  [199] 5.7 3.7 4.2 3.4 2.8 2.8 5.1 5.1 3.9 4.7 4.6 6.0 4.4 2.6 2.3 6.4 5.2 2.4
##  [217] 4.5 4.8 6.2 4.7 6.5 5.1 4.3 3.4 4.2 5.6 4.7 4.0 3.8 4.5 4.3 5.3 5.7 4.3
##  [235] 4.7 5.1 6.1 3.2 6.0 2.9 5.4 3.9 5.8 4.1 4.0 5.0 2.8 4.4 4.5 3.6 5.6 5.2
##  [253] 4.1 4.2 2.8 4.1 4.2 4.2 5.1 4.7 3.9 4.2 4.1 5.1 3.3 4.1 4.8 3.2 2.6 5.0
##  [271] 3.3 6.2 4.3 4.2 6.3 4.1 2.9 5.9 3.1 4.3 3.0 5.0 3.4 4.1 5.3 3.7 5.3 4.5
##  [289] 3.8 3.9 4.2 4.0 5.5 4.1 5.1 5.6 6.3 4.1 6.7 2.2 3.8 4.8 5.3 3.1 5.0 4.1
##  [307] 3.4 3.4 5.4 4.0 4.3 3.3 2.7 4.1 6.1 5.1 4.3 5.0 4.1 3.7 3.5 5.1 4.2 3.5
##  [325] 4.4 4.3 5.0 4.9 3.7 5.3 5.8 4.2 5.4 5.3 4.4 5.4 4.5 4.9 5.3 5.1 5.2 6.0
##  [343] 4.1 3.9 4.7 3.8 4.0 4.7 3.8 5.3 4.1 3.6 5.1 5.6 4.9 4.8 4.5 5.0 3.7 4.5
##  [361] 5.1 5.4 6.2 4.4 5.8 3.7 5.0 4.4 4.5 3.7 4.8 5.7 5.9 5.7 3.5 5.0 4.2 4.7
##  [379] 4.2 3.7 4.3 4.1 4.9 3.1 5.5 4.5 5.5 4.9 4.4 5.7 3.1 3.6 5.3 6.4 5.1 5.7
##  [397] 7.1 4.7 4.3 4.6 2.6 5.3 3.8 6.4 5.2 4.3 5.0 2.7 4.6 6.0 3.5 3.2 4.9 6.0
##  [415] 4.1 4.1 5.9 4.8 6.2 4.0 5.6 5.5 3.8 5.8 4.2 4.9 4.0 5.7 3.6 3.5 6.7 4.2
##  [433] 4.5 4.2 5.3 4.5 3.7 4.4 5.3 4.6 3.4 4.9 5.9 4.5 3.7 5.3 4.6 5.7 3.9 5.7
##  [451] 4.1 4.5 3.2 4.8 4.4 6.4 4.6 4.3 5.2 4.3 5.5 6.3 4.7 6.8 4.3 3.7 4.6 5.4
##  [469] 5.4 4.2 4.2 5.1 3.8 6.0 3.6 4.7 4.8 4.3 4.1 6.8 5.0 3.9 3.6 5.0 3.9 5.3
##  [487] 3.9 5.1 3.8 5.2 3.1 5.3 4.3 3.3 4.6 3.9 4.6 4.9 3.5 4.2 5.4 4.1 5.0 5.4
##  [505] 5.6 4.8 2.8 4.4 4.7 5.9 3.1 3.7 4.1 4.3 3.3 5.6 3.6 4.7 4.5 5.3 2.7 4.1
##  [523] 6.5 4.8 4.0 4.5 4.6 4.5 4.1 2.8 4.9 4.5 3.7 5.6 3.8 3.7 4.3 4.2 4.8 4.8
##  [541] 6.6 3.4 4.3 4.0 3.6 4.9 3.8 5.8 3.8 3.5 3.7 5.3 3.9 4.9 5.3 4.8 4.6 4.0
##  [559] 5.6 4.5 4.2 2.6 5.1 4.3 4.0 4.1 5.5 4.6 3.8 4.0 3.9 3.5 5.5 5.2 4.0 3.8
##  [577] 5.3 4.3 3.0 3.7 4.7 4.2 4.6 4.9 2.8 5.0 5.4 4.1 3.8 6.3 4.1 5.0 4.2 4.3
##  [595] 5.9 5.0 4.7 5.3 6.3 4.1 4.3 4.0 5.7 4.5 5.4 5.1 3.7 5.7 4.8 4.4 5.3 3.7
##  [613] 4.3 5.1 6.4 3.8 4.0 3.8 4.9 4.6 4.9 4.2 4.1 4.0 3.5 4.4 5.7 4.2 3.9 6.5
##  [631] 3.2 4.9 5.3 4.7 5.5 4.6 4.3 4.1 4.8 5.2 4.0 4.1 4.9 3.0 4.3 3.4 5.5 3.8
##  [649] 3.8 6.2 3.7 4.2 4.6 5.7 4.0 4.5 4.6 3.7 3.2 3.9 3.5 5.3 5.3 5.5 5.0 5.2
##  [667] 4.7 4.5 3.0 2.3 6.6 4.3 5.0 4.6 5.1 4.9 3.7 3.0 3.5 5.4 4.4 6.0 4.6 4.0
##  [685] 5.9 3.3 3.7 3.1 5.5 2.4 4.3 5.2 3.7 5.5 5.5 5.5 5.3 4.0 3.6 5.0 5.8 5.1
##  [703] 5.3 3.3 5.8 6.8 5.0 4.1 4.1 4.1 4.2 3.8 5.6 4.8 2.6 4.9 4.2 3.4 3.0 3.4
##  [721] 4.3 3.3 4.6 4.0 5.3 5.2 5.1 4.2 4.8 3.7 4.0 5.2 4.7 3.9 4.6 5.4 3.2 4.1
##  [739] 3.4 2.1 4.9 3.7 4.6 5.5 3.4 5.8 3.7 4.1 5.4 5.4 5.2 5.5 2.4 5.4 3.3 3.6
##  [757] 5.0 2.6 5.0 5.0 5.5 5.2 4.8 4.3 4.7 4.3 4.3 4.7 4.2 3.7 4.7 4.6 5.6 6.2
##  [775] 4.9 4.2 5.9 3.8 4.9 3.8 4.1 5.8 3.4 4.7 4.1 5.9 3.4 5.2 4.4 4.4 4.0 4.3
##  [793] 5.6 5.0 5.4 4.7 4.9 4.5 4.8 4.6 4.5 6.0 4.7 5.2 2.7 5.4 6.9 4.6 4.3 3.9
##  [811] 6.4 3.0 3.0 5.3 6.3 4.4 3.5 6.1 3.4 4.8 3.7 2.0 5.0 5.0 4.5 5.5 3.9 4.6
##  [829] 5.3 3.9 4.2 3.4 4.1 4.8 3.7 6.4 4.1 4.2 3.1 4.7 5.2 2.7 5.8 5.6 6.0 4.6
##  [847] 3.7 4.0 3.0 2.9 6.2 3.7 3.6 3.7 4.5 4.9 4.1 3.3 4.8 3.5 4.8 5.1 3.8 3.3
##  [865] 4.2 4.3 4.8 6.1 4.6 4.4 3.8 4.3 5.1 3.8 3.9 3.1 4.1 4.5 6.4 5.1 3.8 5.1
##  [883] 3.4 4.7 3.1 4.2 5.3 5.8 6.0 5.1 5.9 4.0 5.1 5.1 4.4 3.5 3.6 5.2 4.2 5.6
##  [901] 4.0 4.3 4.9 4.7 3.2 2.8 3.4 4.4 5.7 2.3 3.9 4.7 2.3 5.7 5.3 4.7 4.6 5.2
##  [919] 5.5 4.6 4.4 2.7 5.5 3.8 5.5 4.2 2.5 5.2 4.0 4.7 4.2 4.1 3.9 5.7 3.3 5.1
##  [937] 4.9 4.4 4.8 6.0 3.8 3.9 4.4 4.3 5.2 5.4 2.9 5.5 5.4 3.1 3.5 4.2 3.7 4.4
##  [955] 6.4 5.3 4.8 4.8 3.7 3.8 5.7 4.9 3.7 2.8 4.9 3.7 5.0 3.9 5.4 3.4 3.0 5.1
##  [973] 2.0 4.4 3.9 5.8 4.9 4.5 4.3 4.8 5.0 5.3 4.8 5.0 3.9 5.9 5.7 4.1 3.1 5.6
##  [991] 5.8 4.6 4.2 4.4 5.3 4.9 5.9 6.3 4.3 4.9
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
##   2.7   6.4
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
##    [1] 4.9 4.3 4.2 5.7 3.2 5.4 4.5 4.7 5.6 5.5 4.6 6.1 4.6 4.1 5.3 4.4 5.1 3.5
##   [19] 4.1 5.3 4.4 4.3 4.2 3.5 4.7 5.5 4.3 3.7 6.5 4.0 2.9 5.1 4.1 5.1 3.1 4.9
##   [37] 4.0 5.0 4.5 3.9 5.0 5.7 2.8 4.0 4.4 5.4 6.6 6.0 4.0 3.0 4.5 4.5 5.6 4.8
##   [55] 4.9 5.4 5.1 3.8 3.8 5.8 5.5 5.2 4.0 3.7 5.4 4.6 3.9 5.6 3.6 4.5 4.4 6.2
##   [73] 3.8 3.7 4.2 3.3 5.2 3.1 3.5 5.1 3.2 6.1 5.2 4.8 2.8 5.1 6.0 4.2 4.2 3.6
##   [91] 5.1 5.7 3.2 4.5 5.4 3.8 3.7 5.4 3.3 6.0 4.4 4.7 4.8 4.7 3.3 4.0 4.3 5.1
##  [109] 4.8 5.2 5.9 4.2 4.0 3.7 4.7 5.5 4.1 5.9 5.8 5.6 3.0 4.0 4.9 4.9 2.6 4.6
##  [127] 4.5 4.8 5.8 5.0 5.0 3.9 4.2 2.7 3.3 3.8 4.9 4.7 5.6 3.8 3.7 4.8 5.1 2.7
##  [145] 5.0 3.4 4.7 3.5 3.6 2.7 4.0 3.8 5.1 4.9 3.8 5.7 4.6 3.6 4.3 3.3 3.6 6.0
##  [163] 3.1 4.6 3.0 4.6 4.3 5.3 6.0 3.9 5.1 2.5 3.3 5.6 4.7 3.3 5.9 3.9 5.1 4.7
##  [181] 4.9 5.3 3.7 3.8 4.2 5.6 4.7 5.9 3.6 5.1 3.5 4.9 3.8 5.3 4.0 5.4 3.8 4.0
##  [199] 5.5 4.0 4.3 4.0 4.5 4.2 5.4 5.4 3.5 3.7 4.5 6.2 3.6 5.7 5.3 3.1 5.2 6.0
##  [217] 4.9 4.6 4.8 4.0 2.9 3.2 4.9 3.7 5.2 3.3 4.1 4.6 5.4 4.9 5.1 4.6 6.5 5.8
##  [235] 5.7 5.0 4.1 4.7 5.5 4.6 4.5 4.7 3.2 5.0 5.8 5.1 4.4 5.0 5.0 4.0 4.2 4.6
##  [253] 3.0 5.0 3.6 6.0 4.8 4.4 4.4 5.0 4.5 4.3 4.6 4.9 3.5 5.8 3.9 3.6 4.6 6.0
##  [271] 4.8 4.1 3.7 4.6 4.7 3.1 3.5 5.4 3.9 6.7 4.5 3.9 3.0 3.5 4.1 3.5 4.7 6.0
##  [289] 4.3 3.6 4.6 4.9 3.2 4.9 4.3 4.5 3.9 4.7 3.9 5.4 4.5 4.7 1.4 4.9 5.1 5.4
##  [307] 4.7 5.1 3.9 6.2 5.4 3.9 3.4 5.0 4.5 5.2 5.2 5.6 5.2 3.0 4.1 4.7 6.1 3.6
##  [325] 3.8 5.1 4.7 5.6 4.0 2.9 5.8 5.3 3.5 5.2 4.9 4.8 5.0 4.0 5.3 5.0 3.4 6.2
##  [343] 5.6 5.2 3.9 4.5 5.2 4.0 5.9 5.1 4.7 5.4 4.2 6.2 4.0 3.9 5.7 5.6 5.3 6.0
##  [361] 4.7 4.9 5.9 4.4 4.1 3.4 2.8 2.9 4.4 5.4 4.2 3.1 3.1 5.4 3.2 4.7 3.6 4.4
##  [379] 4.4 3.7 4.7 4.8 5.2 3.9 3.5 5.2 2.1 4.1 3.8 4.2 5.8 4.6 4.5 6.1 5.0 5.1
##  [397] 3.9 3.4 4.1 3.2 5.1 4.9 4.6 4.1 4.1 3.5 3.5 6.1 6.2 5.1 6.7 4.7 4.6 3.7
##  [415] 5.3 5.5 5.3 4.4 3.1 4.6 5.0 2.6 4.8 4.6 3.6 4.2 3.7 3.9 5.4 5.2 3.0 4.3
##  [433] 4.1 5.1 4.7 4.1 6.1 4.3 3.7 3.8 4.7 4.3 3.8 4.9 5.3 3.0 4.5 4.0 3.6 4.1
##  [451] 6.0 3.8 3.7 6.1 5.6 3.8 2.6 3.8 4.4 4.9 4.8 3.7 2.8 6.0 4.2 4.9 4.4 5.6
##  [469] 4.6 3.0 2.1 5.3 3.9 6.0 5.3 2.8 3.8 5.3 2.8 4.4 4.1 3.9 4.2 3.7 5.7 4.3
##  [487] 5.1 5.8 2.8 5.0 4.2 5.8 4.8 4.8 4.7 5.3 5.4 4.3 3.6 4.8 3.4 5.4 5.1 4.6
##  [505] 5.7 4.9 5.6 3.9 5.1 5.7 4.0 6.2 5.5 3.6 4.7 5.0 5.0 6.5 4.8 4.0 4.5 4.2
##  [523] 3.1 5.5 4.5 3.8 3.5 4.6 3.8 4.0 2.7 4.6 3.7 4.7 4.6 5.8 5.6 3.8 4.0 5.3
##  [541] 3.0 4.5 4.6 5.3 3.4 3.6 3.8 6.3 5.2 4.6 5.9 4.3 4.8 5.0 5.0 4.2 3.8 4.6
##  [559] 4.7 4.1 6.0 6.1 3.8 4.8 3.2 4.0 5.4 3.7 3.2 2.9 4.1 4.9 3.4 5.7 4.9 5.8
##  [577] 2.7 5.6 4.9 4.4 5.1 5.0 4.7 3.6 4.3 5.4 5.0 4.6 4.6 5.1 2.3 4.7 3.5 5.1
##  [595] 3.7 4.4 3.5 2.9 5.2 3.8 3.8 5.8 3.2 3.1 4.3 4.1 4.4 4.3 3.8 5.0 3.3 5.7
##  [613] 4.2 3.2 4.1 4.8 5.8 4.7 4.6 4.4 4.4 5.5 5.1 4.2 5.2 4.0 4.8 4.8 4.1 4.4
##  [631] 3.3 4.5 4.2 4.3 4.8 3.7 3.0 3.1 4.2 4.3 4.2 3.7 3.6 4.1 4.8 4.4 4.1 4.3
##  [649] 5.3 3.4 4.1 4.4 4.9 3.7 4.5 4.9 5.6 4.7 3.9 3.0 4.5 4.3 4.6 5.6 3.7 4.7
##  [667] 6.7 5.6 5.3 3.2 4.8 5.6 4.2 3.9 2.4 2.9 4.8 4.2 5.3 5.5 3.0 3.2 3.9 2.5
##  [685] 4.7 5.2 4.7 4.1 5.1 4.9 4.6 3.4 3.8 5.3 5.7 5.4 4.5 3.8 3.9 3.9 5.3 5.4
##  [703] 3.7 3.6 4.7 5.9 4.8 4.9 4.7 5.4 5.3 3.4 3.9 3.4 3.8 3.7 5.4 2.9 5.1 3.8
##  [721] 5.5 3.5 3.1 4.3 4.2 4.5 4.9 2.6 5.7 4.3 6.1 2.6 4.5 4.4 3.8 5.5 5.0 4.3
##  [739] 3.8 2.2 5.2 4.5 4.9 5.0 4.0 4.8 4.7 3.6 2.5 5.9 4.3 5.4 4.1 5.6 4.5 2.5
##  [757] 2.6 3.9 3.3 5.2 4.8 4.1 4.5 6.6 4.2 5.8 7.3 5.3 6.4 5.1 3.6 4.1 4.0 3.5
##  [775] 3.0 2.6 6.0 6.1 4.6 3.9 4.2 3.0 5.9 4.9 6.0 3.7 3.8 4.4 5.5 5.1 4.0 3.3
##  [793] 4.4 3.3 2.7 3.8 4.0 3.8 4.6 6.0 6.4 4.9 3.7 3.7 4.3 5.0 4.9 3.1 3.9 3.7
##  [811] 3.5 5.2 4.9 4.9 3.9 5.4 4.4 5.4 5.3 5.2 3.9 3.3 6.3 3.8 6.4 5.4 4.4 4.0
##  [829] 4.8 4.5 3.4 6.2 4.7 6.3 5.3 4.4 5.5 4.2 5.1 4.3 5.9 5.0 5.4 3.3 3.8 4.9
##  [847] 4.2 5.8 3.8 3.8 6.0 4.0 5.5 4.0 5.8 5.2 5.2 6.9 4.6 5.4 4.9 4.3 5.0 4.7
##  [865] 3.6 5.4 4.7 4.5 5.5 3.4 5.6 3.0 3.9 4.8 2.1 4.1 5.9 4.6 5.5 4.1 3.6 3.9
##  [883] 3.8 4.0 4.9 5.4 4.3 3.6 5.5 4.8 5.4 4.2 4.8 2.9 5.3 5.8 4.3 3.0 4.4 4.0
##  [901] 4.6 3.8 2.7 5.5 4.6 4.3 3.7 4.3 5.1 3.6 4.2 5.0 6.5 4.8 4.5 4.8 3.2 4.7
##  [919] 3.3 4.4 3.1 4.8 4.6 3.3 6.2 4.8 4.0 5.0 4.6 4.3 5.4 3.3 4.1 3.2 5.3 4.9
##  [937] 3.5 2.4 5.7 4.7 4.3 3.4 2.9 4.4 4.6 5.4 3.9 2.9 4.1 3.3 5.5 4.8 4.6 4.6
##  [955] 6.1 4.3 5.6 3.8 4.4 4.8 6.5 2.8 5.0 4.9 4.6 4.3 4.1 3.5 6.3 5.3 3.7 4.3
##  [973] 4.7 4.7 4.7 4.2 5.8 5.0 4.7 4.5 4.5 3.5 5.6 4.7 4.0 5.6 3.3 5.6 5.3 5.3
##  [991] 4.4 4.9 3.1 3.7 4.5 3.5 4.6 4.5 5.6 4.6
## 
## $func.thetastar
## [1] -0.0037
## 
## $jack.boot.val
##  [1]  0.514005602  0.344507042  0.294117647  0.186227545  0.001369863
##  [6] -0.013649025 -0.264564565 -0.282300885 -0.398011364 -0.514124294
## 
## $jack.boot.se
## [1] 0.9835848
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
##    [1] 4.8 3.0 5.9 4.6 4.6 3.9 4.7 2.7 4.0 4.3 4.2 4.6 3.5 4.7 4.7 5.8 4.4 4.2
##   [19] 3.8 5.1 5.1 3.9 4.7 4.8 3.2 2.8 4.9 5.4 5.2 4.6 6.3 3.2 4.2 3.9 5.4 3.0
##   [37] 3.7 5.8 4.3 5.9 4.9 3.3 3.9 4.5 4.5 6.1 4.2 4.8 5.1 6.9 2.4 4.7 5.5 2.6
##   [55] 4.4 2.9 4.4 4.9 5.3 5.2 4.2 5.8 4.9 4.0 4.3 3.3 5.4 5.1 4.1 3.7 5.3 4.1
##   [73] 3.7 5.2 4.4 3.1 3.9 5.3 5.5 5.0 6.2 3.9 4.0 6.7 2.9 2.9 4.3 5.7 4.9 4.5
##   [91] 3.0 2.2 5.0 4.3 4.1 3.3 4.0 2.0 4.9 5.1 5.1 5.7 4.2 4.5 5.9 4.5 4.8 6.3
##  [109] 3.7 4.3 3.3 5.0 4.4 5.4 4.0 3.7 5.9 3.0 5.4 2.1 4.9 3.0 3.2 3.4 4.7 5.1
##  [127] 5.4 5.4 5.1 5.0 4.6 4.6 4.4 5.6 5.6 2.9 4.3 4.3 4.0 3.9 3.7 5.0 5.9 5.2
##  [145] 4.5 3.6 4.9 4.6 5.1 6.0 4.1 3.4 4.0 5.9 3.5 4.3 5.9 5.4 4.9 4.6 3.6 4.4
##  [163] 4.4 4.5 4.6 5.8 4.3 4.7 3.5 3.3 2.6 4.6 6.3 4.0 3.9 4.2 4.3 4.8 5.5 4.8
##  [181] 5.5 3.5 3.9 5.0 5.6 4.8 4.7 3.1 4.9 4.4 3.4 4.3 4.2 3.4 4.6 3.7 3.1 3.5
##  [199] 5.9 3.1 4.0 2.2 5.1 4.8 4.5 3.9 5.1 4.1 5.7 4.9 3.5 3.8 5.2 5.3 4.3 4.7
##  [217] 3.4 4.7 4.1 4.8 4.3 5.3 4.6 6.2 4.6 3.6 4.0 4.8 3.7 5.9 6.9 5.9 4.8 5.1
##  [235] 4.8 3.4 6.4 4.8 3.6 4.0 4.2 4.9 5.5 4.7 3.4 2.8 4.8 3.5 5.3 4.0 4.0 4.1
##  [253] 5.1 5.0 6.5 4.8 4.1 7.3 5.5 4.2 2.7 3.5 3.3 4.9 4.4 4.3 4.8 4.1 3.8 4.8
##  [271] 3.8 5.1 4.8 3.8 4.8 4.9 6.0 4.4 5.0 3.7 4.9 4.4 4.4 4.0 3.9 5.3 4.2 5.6
##  [289] 5.0 3.5 3.7 3.5 3.8 6.7 4.4 4.6 5.3 3.1 3.2 2.7 5.1 4.8 3.6 4.1 5.5 3.5
##  [307] 4.9 5.9 5.8 4.9 4.2 4.4 4.4 4.5 4.4 5.0 4.1 4.8 3.7 4.3 4.4 3.7 2.9 4.3
##  [325] 5.6 5.8 3.9 3.0 6.0 5.6 5.8 5.5 5.4 6.0 2.3 4.1 3.8 5.4 3.7 3.1 5.4 5.4
##  [343] 4.8 5.7 3.9 3.9 5.3 3.4 5.0 3.7 4.8 5.6 4.6 2.7 4.6 5.5 4.2 4.0 3.8 4.5
##  [361] 5.4 5.7 3.9 2.1 3.7 3.0 4.7 5.4 2.2 3.6 4.6 5.2 5.0 5.7 3.2 5.6 5.7 4.7
##  [379] 4.8 5.8 4.9 3.4 4.4 3.1 5.9 4.4 3.2 4.3 5.5 5.7 3.6 4.0 4.7 3.3 5.5 3.4
##  [397] 4.4 4.3 4.1 5.9 2.7 4.0 5.7 5.1 4.4 4.0 4.7 6.2 3.4 5.4 4.8 4.8 5.6 4.5
##  [415] 5.1 4.5 4.5 4.6 5.4 4.2 5.3 5.7 5.4 4.5 4.7 3.0 3.6 4.1 5.6 5.5 3.9 4.5
##  [433] 4.4 4.3 6.5 3.3 5.6 4.5 5.1 4.2 5.9 4.1 2.8 2.6 6.4 5.3 6.0 4.2 4.2 5.5
##  [451] 4.7 5.3 2.5 6.0 4.3 5.4 5.5 5.4 3.3 5.5 5.2 3.8 2.9 2.8 5.1 3.5 3.4 2.9
##  [469] 4.9 5.2 5.4 5.4 4.8 5.6 1.9 4.7 6.3 3.5 4.2 4.4 3.7 5.0 4.6 4.0 4.6 5.2
##  [487] 5.2 4.6 4.9 6.5 5.3 6.3 2.8 4.8 3.3 4.8 5.1 4.5 3.5 4.0 4.8 2.8 5.1 3.8
##  [505] 3.6 6.6 2.8 3.3 4.5 3.3 4.4 3.4 4.9 4.7 3.7 5.2 5.0 4.0 3.0 3.5 4.2 4.4
##  [523] 4.2 2.8 5.3 4.4 3.9 5.1 3.8 4.8 2.4 4.2 5.5 4.1 4.9 4.3 5.5 5.7 4.3 4.5
##  [541] 3.8 4.2 4.5 5.7 4.6 3.3 4.6 5.0 5.1 1.8 5.8 5.2 5.0 3.1 5.1 3.9 5.7 4.4
##  [559] 5.5 4.0 3.8 5.1 3.3 3.6 5.9 3.3 5.1 3.9 3.7 4.7 3.8 3.8 3.6 5.6 4.4 5.6
##  [577] 3.3 5.5 3.0 4.4 4.0 4.4 4.8 4.1 4.8 5.6 5.6 4.1 5.6 4.8 4.7 4.1 4.3 5.0
##  [595] 3.5 5.0 3.9 4.6 4.3 5.3 5.8 5.6 5.2 5.2 4.9 5.4 3.9 5.3 4.7 3.6 3.8 5.4
##  [613] 4.0 4.3 2.8 4.7 4.3 4.8 4.2 4.5 5.7 3.8 6.1 4.1 7.0 5.3 5.4 4.2 3.9 4.2
##  [631] 3.2 4.0 4.9 4.8 4.1 4.6 3.2 5.7 5.2 5.2 3.6 4.3 4.4 5.0 4.8 4.7 5.4 4.8
##  [649] 5.3 4.3 4.2 4.6 3.9 4.8 4.0 4.4 5.2 4.2 4.8 3.3 4.0 3.5 5.4 3.7 2.6 5.3
##  [667] 5.3 6.3 3.3 4.8 3.7 5.2 4.6 4.2 5.0 6.4 5.9 5.7 4.5 3.6 6.2 3.8 4.2 3.2
##  [685] 5.1 4.7 5.0 4.6 4.6 4.6 4.1 4.9 4.6 3.4 4.0 6.8 6.1 3.9 5.8 4.2 5.1 3.2
##  [703] 4.4 3.0 4.1 3.1 4.6 5.3 5.9 4.7 4.6 4.0 3.5 4.9 5.2 3.9 5.0 6.0 3.9 4.9
##  [721] 4.1 4.7 5.9 4.4 2.0 4.4 5.0 3.5 3.8 5.2 6.2 2.5 4.0 4.7 4.6 4.4 4.4 5.6
##  [739] 5.2 4.5 5.4 5.1 5.2 4.1 5.7 4.1 2.5 4.5 4.3 5.2 4.7 6.4 3.8 5.6 4.5 2.5
##  [757] 4.3 3.7 5.9 5.1 5.5 4.6 4.3 3.7 3.8 6.4 5.2 4.8 4.3 4.6 4.8 3.6 4.2 3.1
##  [775] 5.7 4.4 4.3 4.7 4.4 3.2 5.4 5.8 6.1 4.5 4.3 5.1 3.7 4.0 5.1 5.6 2.3 3.2
##  [793] 4.2 4.1 4.5 5.0 4.0 4.6 4.9 3.8 5.0 3.3 4.7 5.0 3.9 3.7 4.7 1.9 3.8 4.7
##  [811] 4.1 4.4 5.0 4.9 3.0 3.6 4.3 5.7 3.3 4.8 3.2 3.2 3.9 4.6 3.3 4.0 5.1 4.0
##  [829] 4.0 3.5 4.1 4.6 4.2 4.6 4.4 3.3 3.7 5.1 4.5 3.2 4.5 4.8 4.5 6.4 4.5 4.7
##  [847] 4.3 2.8 2.7 4.6 4.5 5.9 4.4 5.2 3.3 5.3 3.5 3.9 5.0 3.2 4.6 4.9 5.1 5.4
##  [865] 2.7 4.3 3.1 4.7 4.7 3.5 4.8 4.7 3.6 4.7 3.4 5.2 3.7 4.4 3.3 4.0 6.5 5.0
##  [883] 4.5 2.7 3.7 3.4 4.4 4.3 5.3 5.1 2.5 3.4 4.2 3.2 4.2 5.4 4.9 5.2 5.3 5.7
##  [901] 2.5 4.7 5.6 3.2 4.3 4.0 6.2 5.3 3.7 5.1 2.7 4.8 4.6 6.2 4.0 5.2 5.3 3.7
##  [919] 5.5 5.0 4.6 4.3 4.2 3.0 4.8 3.3 3.5 6.4 5.3 3.9 4.4 4.1 5.7 3.3 3.2 5.3
##  [937] 3.7 3.2 4.1 3.0 6.3 4.0 3.3 5.4 6.5 4.4 5.0 6.3 7.3 3.9 4.4 5.7 4.8 3.5
##  [955] 3.4 5.0 4.4 4.6 5.1 4.9 4.6 4.2 4.1 4.7 3.4 4.1 4.1 5.1 6.6 3.1 3.2 3.6
##  [973] 4.6 4.8 4.7 3.9 5.3 3.7 5.6 3.4 3.5 3.5 5.1 4.8 5.7 4.9 4.3 3.5 5.0 3.5
##  [991] 4.3 5.2 4.0 4.8 6.4 5.4 5.9 5.5 5.8 5.7
## 
## $func.thetastar
## 72% 
## 5.1 
## 
## $jack.boot.val
##  [1] 5.500 5.432 5.300 5.300 5.100 5.056 5.000 4.800 4.600 4.500
## 
## $jack.boot.se
## [1] 0.9698799
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
## [1] -0.3603941
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
##     shape       rate  
##   2.033970   3.419282 
##  (0.845648) (1.611143)
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
## [1]  1.07203558  0.68785841  1.46745753  1.55939789  0.33144917 -0.07215129
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
##    [1] -0.841101235  0.479457737 -0.549252328  0.333256276 -0.490182201
##    [6] -0.262689139 -0.454980551 -0.535381687 -0.625267385 -0.024152457
##   [11] -1.179335694 -0.680148731 -0.847735813  0.348741340 -0.795613984
##   [16] -0.073305125  0.713003126 -0.007746502 -0.007829859 -1.569724795
##   [21] -0.931839741  0.178487078 -0.517462223 -0.047966725 -0.168931867
##   [26]  0.505118347 -0.535849085 -0.014427209 -0.312700000 -0.920765604
##   [31] -0.392377591 -0.011635242  0.014970060 -0.910922026 -0.095961929
##   [36]  0.168394055  0.196756066 -0.882726853 -0.604891149 -0.295777539
##   [41] -0.275377401  0.120142644  0.580790222 -0.807946003  0.276280522
##   [46] -0.117860339  0.278611412 -1.037002912 -0.626641392 -0.838621973
##   [51] -0.163206515  0.171270863 -0.689154356 -0.081260337 -0.634734464
##   [56]  0.081600681  1.172459550 -0.687612885 -0.259612481 -0.122815007
##   [61] -0.194621298 -0.118617765  0.175728854 -0.422031649  0.669510288
##   [66] -0.529737214 -0.073177476  0.175149544 -0.022040316 -0.139684674
##   [71] -0.919433859 -0.853687012 -1.036305305 -0.912910159 -0.989287689
##   [76]  0.864801555 -0.595889338 -0.222704402 -0.503311126 -0.937157675
##   [81] -0.284252512 -0.361986199 -0.562196895 -0.316235828 -0.906172603
##   [86] -0.281063131 -0.152690182 -0.672339476 -1.145910321 -0.466943525
##   [91]  0.199774869  0.159284501 -1.140678500 -0.248918433 -0.734267932
##   [96]  0.553706852 -0.610915634 -0.496699932  0.040612948  0.661857398
##  [101] -0.466826454 -0.182123046  0.031338066 -1.593785174  0.795095490
##  [106] -1.146889073 -0.604343823 -0.398498643  0.188173613  0.075589618
##  [111] -0.348398607 -0.008027782 -0.044061096 -0.729638516 -0.090363590
##  [116] -0.547330169  0.145898148 -0.647480014  0.424996539  0.048086539
##  [121] -1.145755997 -0.364615615 -0.571590581  0.763049003 -0.114127237
##  [126] -0.798806391 -0.654994940 -0.125452806 -0.849927370 -0.030512808
##  [131] -0.252644184 -0.506622161 -0.489722974 -0.452623376  0.372446847
##  [136] -0.418690640 -0.233232988 -0.827825877 -0.505246601 -0.983858644
##  [141] -0.564162045  0.134023775  0.057985396 -1.157470308  0.300084406
##  [146] -0.350932377 -0.592395758 -0.824308374 -0.034483120  0.078781269
##  [151] -1.659766039 -0.553842274 -0.064510105 -0.522768112  0.135121848
##  [156] -0.269857737 -0.447966238 -0.284729482 -1.258593944 -0.256898195
##  [161] -1.540089617 -0.494306209 -0.462368071 -0.171960952 -0.183978709
##  [166] -0.697270694 -0.222487589 -0.150172540 -1.123714906 -0.561934951
##  [171] -0.097415741 -0.458973219 -0.659640628 -0.410857191 -1.116057603
##  [176]  0.174515988 -0.208413591 -0.911734043  0.524320785  0.409045560
##  [181] -0.490522636  0.308091287 -0.534268276 -0.434435724 -0.428777251
##  [186]  0.407166206 -0.166360796 -0.769066622  0.085596744 -0.287002017
##  [191] -0.256891255 -0.234261223 -0.427332617 -0.663666997 -0.259783444
##  [196] -0.302007479  0.580189526 -0.913125789  0.012545102 -0.753352533
##  [201]  0.847589250 -0.041964248 -0.233251832 -1.105449880  0.005193517
##  [206] -1.196345534 -0.595659301 -0.308061875 -0.357932793 -0.108103721
##  [211] -0.420033880 -0.321491224 -0.424073758 -0.968529657 -0.321983087
##  [216]  0.714789690 -0.639427310 -0.127683803 -0.024284846 -0.865519376
##  [221] -0.057496077  0.486822917 -0.917103011 -0.288851445 -0.309579501
##  [226] -0.211479493  0.941710490 -0.542735876 -0.385386796 -0.435051968
##  [231] -0.486056596 -1.048419487 -0.948852851 -0.118287026  0.196707570
##  [236] -0.255106046  0.806319423 -1.009932223 -0.099041277 -0.536513947
##  [241] -0.396406595 -0.634592972 -0.143946409  0.487632260  0.341423114
##  [246] -0.561623000 -0.379722279 -0.138139008  0.184194555 -0.240885558
##  [251] -0.430543921 -1.291634062 -0.553234892  0.198482137  0.401341693
##  [256] -1.086102590 -0.674705263 -0.363283135 -0.715129227  0.243276315
##  [261] -0.729560573 -0.457236899 -0.215568270  0.166341279  0.341200720
##  [266] -0.034750110 -1.421746174 -0.268771400  0.136773775  0.045573177
##  [271] -0.596836475  0.141755146  0.423745051 -0.411176082 -0.587150547
##  [276] -0.262412969  0.117033705 -0.079024615 -0.562867494 -0.793033779
##  [281] -0.849363841 -1.860531806 -0.081673533 -0.235070048 -0.632955104
##  [286] -0.575399322 -0.978231437 -1.100816382 -0.067537425 -0.153698149
##  [291] -0.531734778 -0.467766623  0.107017527 -0.888754769 -0.582559109
##  [296]  0.244681130 -0.594015223 -0.474676204 -0.613463839 -0.497191026
##  [301] -0.360394102  0.955399370 -0.963918885  0.029686322 -1.173012592
##  [306] -0.614422027 -0.326977149 -0.406918248  0.346057013 -0.200286358
##  [311]  0.005603147 -0.119812014 -0.898749415 -1.040134088 -0.444552826
##  [316] -1.172515532 -0.986953649 -0.905965673  0.135974081 -0.448402957
##  [321] -0.538975655 -0.806110568 -1.106584760 -0.264723487 -0.421803404
##  [326]  0.010323010 -0.385370182  0.078577116 -1.309182273 -1.019534543
##  [331] -0.208170880 -0.694978713  0.582280802 -0.384105385 -0.417932166
##  [336] -0.146631684  0.016401256 -0.227792036  0.285686850  0.138779047
##  [341]  0.072556015 -0.539937315 -0.777645899 -0.139199614 -0.114557343
##  [346]  0.442743688 -1.225452635 -1.186706177  0.176919996 -0.269904894
##  [351] -1.259240630 -0.593387478 -0.007935749  0.255233547 -0.019162963
##  [356]  0.295102968  0.083702205 -0.290369016 -0.739849417 -0.254116527
##  [361] -0.133079272 -0.168043745 -0.385831902 -1.129566930  0.433782987
##  [366] -0.305426607  0.023725242 -1.863714323  0.210101297  0.348006108
##  [371]  0.127156223 -0.366504861 -1.051938929  0.260189453 -0.660296073
##  [376] -0.769733499 -0.186491516 -0.204226260 -0.046681071 -0.233079639
##  [381]  0.771026043  0.057985396 -0.112181495 -0.535915619 -0.276669363
##  [386] -1.181283473 -0.306315329 -0.733051310 -0.273344684 -0.770794919
##  [391]  0.169314446 -0.136576516 -0.807740382 -0.109180496 -0.229127468
##  [396] -0.170577151  0.256166745 -0.672993548  0.111874095 -0.319897683
##  [401]  0.281892273 -0.726348337 -0.434674658 -0.126436981 -0.623235920
##  [406] -1.144075147  0.029978168  0.134552804 -0.435914018 -1.157802937
##  [411] -1.103770976 -0.687838719  0.314584111  0.251019885 -0.998280172
##  [416] -0.086316103 -0.102345140  0.241966657 -0.464264073 -0.289054352
##  [421] -0.630054168  0.014391174 -0.277639420 -1.482230990  0.365000503
##  [426] -1.116109379 -0.375799311 -0.604328181  0.056315185 -1.017030509
##  [431] -0.510089170 -0.401693736  0.246304385 -0.105058496 -0.344122252
##  [436]  0.053045908 -0.751219900 -0.046681071 -0.487178014 -0.439672388
##  [441] -0.191317300 -0.305090232 -0.633962661 -0.464938348 -0.779706287
##  [446] -0.333419877 -0.078791298 -1.043929287 -0.031557650  0.863599176
##  [451]  0.319223881 -0.376806360 -0.279758188 -0.197643568 -0.076929117
##  [456] -0.212305928 -0.832597047 -0.628972767 -0.494547620 -0.014410429
##  [461] -0.257772559 -0.751656171 -0.151480966  0.133917054 -0.586004549
##  [466] -0.253094133 -0.280395090 -0.462124783 -0.394542928 -0.060671171
##  [471]  0.535843345 -0.210013316  0.290747388 -0.047779836 -1.324910852
##  [476]  0.115742606 -0.662403321 -0.336762897 -1.333314003 -0.074084088
##  [481] -1.288071487 -0.121554040 -0.349290827 -1.300136420  0.273288448
##  [486]  0.380510987 -0.588158610 -0.136912389 -0.160780611  0.002969360
##  [491]  0.142820401 -0.864048300 -0.904143453  0.460366667 -0.644647842
##  [496]  0.318246968 -0.803963819  0.427988729 -0.357871416 -0.597383592
##  [501]  0.270542457 -0.935303868 -0.327253658 -0.154744620 -0.297369070
##  [506]  0.234434045 -0.390103455 -0.365552130 -0.337904362 -0.141177327
##  [511] -1.239849597 -0.472557761 -0.578381177  0.221792022 -0.205289026
##  [516] -0.142161737  0.173773865 -0.663831378 -0.396998197 -0.764235352
##  [521]  0.054999655 -0.783343580 -0.078449987 -0.728685628  0.032932375
##  [526] -1.292314908 -0.609305939 -0.848129602 -0.315058358 -0.666121310
##  [531] -0.329069977 -0.604343823 -0.063376269 -1.153464385  0.313418900
##  [536]  0.195641971 -0.082970920  0.224458279 -0.688604616 -0.111449600
##  [541] -0.237534756 -0.798332065 -0.349373699 -0.385394789  0.070820823
##  [546] -0.580966577 -0.330184190 -0.514792198  0.368515462  0.417333933
##  [551]  0.049195586 -0.145141867 -0.749206264 -0.388078138 -0.271865693
##  [556]  0.223447142 -0.915601329 -0.210978149 -0.987246772 -0.175520840
##  [561] -1.624072583  0.412328744 -0.276293412  0.424036286 -1.362599271
##  [566]  0.606394426 -0.588961400 -0.521492996 -0.184689574  0.295992559
##  [571] -0.377348076 -0.272757746 -1.529453524  0.119082353 -0.961161325
##  [576] -0.037972680 -0.340551059 -1.083969427  0.419599891  0.262731486
##  [581] -0.425702335 -0.204968907 -0.782817077 -0.451501037  0.282942027
##  [586] -0.730948572 -1.517131582 -0.229604645 -1.624072583 -0.514606480
##  [591] -0.256231746 -1.058754672 -0.589231605 -0.537309059 -0.797602387
##  [596] -0.421297868 -1.501511810 -0.385394510 -0.973873400  0.294157656
##  [601] -0.689482365 -0.483388831 -0.218515329  0.033539963  0.342095550
##  [606]  0.033758254  0.935813093  0.035452628 -0.647446106  0.164901237
##  [611]  0.734485654  0.235093018  0.143550493  0.181292094 -0.190162019
##  [616]  0.169314446  0.007594781 -0.183658054 -0.401315397  0.315809392
##  [621] -0.808613400 -0.064142474  0.184672952 -0.291944306 -0.766452556
##  [626]  0.835516529 -0.739548691  0.030373650  0.165791034 -0.122258704
##  [631] -0.077364040 -0.464283814  0.010977273  0.265892187 -0.601402206
##  [636] -0.225261167  0.063476246 -0.968485888 -1.094224774 -0.662391739
##  [641] -0.393930078 -0.284128320 -0.466395977 -0.469974896 -0.944907025
##  [646] -0.420886443  0.014476750  0.058442415  0.133917054 -0.799553595
##  [651] -0.167545964  0.120953793 -1.762868556 -0.115153172 -1.131703490
##  [656] -0.664839837 -0.533130334 -0.232079199 -0.593387478 -0.340582423
##  [661] -0.388465317 -0.607516005 -0.105142046 -0.838621973 -1.152191776
##  [666] -0.691993256 -0.670786237 -0.666788033 -0.054703764 -0.704882433
##  [671]  0.119434219 -0.245159131  0.186548706  0.300259113 -0.375215084
##  [676] -0.218647366 -0.387101976 -0.803306139 -0.845308403 -0.611159250
##  [681] -0.187065258 -0.370001075 -1.326946742 -0.703407022 -0.449546015
##  [686] -1.580827289 -0.774194593 -0.102848463  0.183100876 -0.048234047
##  [691] -0.286222692 -0.596931519 -1.115570572 -0.084520370 -0.391267867
##  [696]  0.464314039 -0.334270068 -0.712513236 -0.526364354 -1.129798137
##  [701] -0.823738935  0.030503008 -0.295111879 -0.127683803  0.167367096
##  [706] -0.419369169 -0.426821283 -0.744096471 -0.837933647 -0.729106232
##  [711] -1.168189267 -0.290556664 -1.195349965 -0.148827356 -0.729691238
##  [716] -0.474297674 -0.166393762 -0.709044447  0.186610446 -0.649810867
##  [721] -0.402594984 -0.230184324 -0.565157627 -0.590795520 -0.107082295
##  [726] -0.659570874 -0.630934394 -0.842258431 -0.599801683 -0.207985153
##  [731] -1.065286428 -1.262562720 -0.245432333 -0.355423154 -0.385414551
##  [736] -0.314034526 -0.322269775 -0.593119251  0.145081891  0.226289621
##  [741]  0.549816130 -0.267251257 -0.763315641 -0.354434563 -0.355602194
##  [746] -0.765884048  0.471435777  0.611141960 -1.216968395 -0.024383274
##  [751] -0.466826454  0.401997638 -0.707327518  0.390774621  0.096843189
##  [756] -0.244996844 -1.066183218 -0.504356446 -1.124555060 -0.656586749
##  [761] -0.387221467 -0.620517366  0.030988247 -0.617962200  0.050093942
##  [766] -0.718226553 -0.222289953  0.070796451 -0.207097240 -0.415198641
##  [771] -0.554133443 -0.947099316  0.181690280  0.251816179 -0.598386329
##  [776] -0.521424209 -0.134965566 -0.237812251 -0.262067948  0.265167936
##  [781] -0.119577504  0.466161900  0.155347837 -1.132705108 -0.714839878
##  [786] -1.771937735  0.132704019  0.238805385 -0.138372783 -0.144536956
##  [791] -0.887605464 -0.472963965 -0.530083798 -0.295041432  0.162910265
##  [796] -0.549962747 -1.157692501 -0.988654073 -0.382802277 -0.051633331
##  [801]  0.014911604  0.043537728  0.332283934 -0.696126006  0.215028473
##  [806] -0.495015489 -0.307366882  0.143273760  0.417209902 -0.333329636
##  [811]  0.063118652 -0.561802912  0.078507764 -0.947487794  0.017745996
##  [816] -1.049493283  0.424435946 -0.302081215 -0.239812579 -0.872755674
##  [821] -0.239668045  0.018832622 -0.409999961 -1.091892060  0.062915726
##  [826] -0.631763898  0.332894056  0.122860948  0.192900429  0.115221976
##  [831]  0.175383715 -0.562823691 -0.313238417 -0.857559757 -0.593026047
##  [836] -0.038892283 -0.589998905 -0.701999695 -0.421803404  0.404454367
##  [841] -0.296788813 -0.047938897 -0.270548333 -0.379708675 -0.205830880
##  [846]  0.022373016 -0.668385439 -0.484520611 -0.808613400  0.094081839
##  [851] -0.170851745 -0.337687292 -0.115000663 -0.193961761  0.202737103
##  [856]  0.231920286  0.285837546  0.468297160 -0.223118845 -0.137294586
##  [861] -0.989858907 -0.222759074 -1.655164306 -0.037979208 -0.190028637
##  [866]  0.385594166 -0.398525459 -0.355265829 -0.145559926 -0.648085165
##  [871] -0.686950003  0.295592803 -1.527763988 -1.133493712 -0.115966282
##  [876] -0.851746610 -0.136183547  1.976166392 -0.338946070 -0.596838134
##  [881] -0.095943786  0.017984832 -1.283681527 -0.446599188 -1.027250628
##  [886] -0.348591015  0.201526061 -0.231206637 -0.580017436  0.817762139
##  [891] -0.149157595 -0.222891867  0.370189383  0.093038265 -0.033210735
##  [896] -0.166673881 -0.264029736 -1.103440535 -0.321664972 -0.290914888
##  [901] -0.510107908 -0.509535415 -0.026558645 -0.115360954 -0.544296592
##  [906] -1.442315822 -0.250404302 -0.383005710 -0.589603684 -0.234541025
##  [911] -0.648205589 -0.443368398  0.346148697 -1.334736454 -0.391818894
##  [916]  0.138779047 -0.121567648 -1.103223018 -0.101886061 -0.022832036
##  [921]  0.045850523 -0.841954274  0.623686977  0.233673313 -0.591647459
##  [926] -0.422620739 -0.098560837  0.183939503 -0.160943448 -0.507418371
##  [931] -0.577946227 -0.360887322 -0.274658913 -0.225932649  0.135317018
##  [936] -0.074502837  0.569561194 -0.136952258 -0.273326153 -1.157470308
##  [941] -0.430948645 -1.090899787 -1.366859322 -0.262528809  0.113936564
##  [946]  0.121853757 -0.404958361  0.488143447  0.178680592 -0.706771542
##  [951] -0.676115698 -1.049751793  0.119453657 -0.331843507  0.940482295
##  [956] -0.382915517 -0.855333919 -0.881859510 -0.428189212 -0.295458987
##  [961] -0.753908805 -0.765200137 -0.346802247 -0.495115945  0.111383052
##  [966] -1.919204209  0.047568323 -0.670680584 -0.192772625 -0.085471984
##  [971] -0.273657666 -1.040984442  0.476713365 -0.769786859  0.108529163
##  [976]  0.225462425  0.085343078  0.194793169 -0.859747059 -0.521619137
##  [981] -0.969692088 -0.052581617 -1.184493014 -0.112766255 -0.134273844
##  [986] -0.033019143  0.312340994 -0.095961929 -0.334766119  0.088215724
##  [991] -0.627639639 -0.192888885 -0.755278915 -0.502629702 -0.680148731
##  [996] -1.027369480  0.276780384 -0.136229817 -0.029827481 -0.128871513
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
##   0.59486216   0.31342660 
##  (0.09911419) (0.07008172)
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
## [1]  0.28057078  0.79887473 -0.55421457  0.09970380 -0.02768275  1.11663515
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
## [1] 0.0233
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9155137
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
## t1*      4.5 -0.01111111   0.8893529
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 1 3 4 6 8 9 
## 1 2 2 2 1 1 1
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
## [1] -0.0098
```

```r
se.boot
```

```
## [1] 0.9313492
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

