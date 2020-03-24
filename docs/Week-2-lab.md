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
## 0 1 3 6 7 8 
## 2 3 2 1 1 1
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
## [1] -0.0077
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
## [1] 2.659406
```

```r
UL.boot
```

```
## [1] 6.325194
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##   2.5%  97.5% 
## 2.7975 6.4000
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
##    [1] 5.3 4.7 4.7 3.0 4.6 4.3 4.4 5.3 5.3 6.2 5.1 5.8 3.3 6.0 5.7 4.8 5.9 2.7
##   [19] 4.4 3.0 4.2 3.8 4.5 5.1 3.7 4.1 4.3 4.3 3.3 4.4 5.0 5.2 4.6 4.2 4.6 3.7
##   [37] 6.0 4.3 6.6 5.2 5.4 2.9 4.1 2.9 4.7 4.4 4.1 6.0 4.4 3.8 4.6 3.2 5.0 3.0
##   [55] 3.2 4.4 3.8 4.5 5.4 3.2 4.0 4.8 4.7 5.0 5.2 6.7 3.9 4.3 3.8 4.2 4.1 4.9
##   [73] 5.5 3.2 3.2 3.8 3.6 3.4 3.9 3.2 4.8 4.0 4.1 3.8 3.8 3.4 4.3 3.8 4.6 4.4
##   [91] 6.2 4.1 5.0 4.2 6.6 4.4 4.8 4.2 5.1 3.0 3.2 5.0 4.0 5.1 5.2 3.1 3.7 5.0
##  [109] 5.1 5.1 4.2 6.4 4.8 4.0 4.2 4.6 4.4 5.2 3.4 4.5 4.0 3.3 5.1 3.8 5.7 5.6
##  [127] 5.4 5.3 2.7 3.8 3.8 4.4 3.3 4.0 3.6 3.1 4.6 5.0 5.6 5.0 4.6 3.1 3.8 3.2
##  [145] 4.0 3.6 5.0 5.0 6.4 4.9 4.8 5.4 4.7 4.0 5.4 3.1 5.5 3.2 5.3 5.2 5.8 4.0
##  [163] 5.4 6.6 3.9 4.7 5.1 3.1 3.9 5.0 3.9 4.1 5.0 5.1 5.9 5.7 4.3 4.5 5.9 4.6
##  [181] 3.9 4.8 3.9 4.7 3.4 5.0 4.8 4.9 4.2 4.5 5.4 4.9 4.8 3.1 4.8 4.3 4.9 4.1
##  [199] 6.0 4.4 4.4 4.9 4.5 5.4 5.4 4.3 4.2 3.8 4.7 5.1 5.2 4.5 3.2 5.1 4.1 4.6
##  [217] 4.5 2.9 4.0 4.6 5.2 5.7 3.7 4.2 5.8 4.3 2.7 5.8 3.7 4.2 4.5 5.3 5.5 5.8
##  [235] 3.8 3.9 5.2 3.9 2.7 5.9 4.2 3.8 6.0 4.2 3.8 4.3 4.2 4.0 4.9 3.6 6.6 3.6
##  [253] 4.8 3.8 3.2 6.0 4.5 3.0 4.0 5.0 4.4 5.5 5.7 2.9 4.6 3.2 5.1 4.6 3.5 4.7
##  [271] 4.8 4.2 3.7 1.9 3.9 2.3 5.7 5.8 4.2 4.5 5.0 3.8 4.3 4.4 3.7 5.3 5.2 4.4
##  [289] 4.4 5.1 5.6 4.6 3.3 5.1 5.2 5.4 6.4 4.7 3.3 5.6 3.7 5.1 2.6 4.7 6.1 5.8
##  [307] 5.5 5.1 5.4 3.4 5.2 4.4 5.8 5.2 5.6 4.7 4.2 3.7 4.7 4.3 4.9 4.4 4.3 3.0
##  [325] 5.0 3.5 3.6 3.5 4.0 5.8 3.0 5.8 5.5 3.6 3.4 3.9 5.1 4.6 4.6 5.7 4.7 4.6
##  [343] 4.5 4.7 3.6 4.4 4.6 4.7 5.0 5.7 2.6 3.3 5.1 3.5 4.5 4.5 4.4 3.4 5.0 3.5
##  [361] 4.4 5.9 3.6 5.9 3.3 3.9 3.9 4.2 3.0 5.9 5.0 4.5 3.5 6.6 4.6 4.0 4.2 5.6
##  [379] 4.0 4.1 5.3 4.4 5.2 4.4 4.8 5.6 4.7 5.0 5.6 3.1 5.9 3.4 5.4 4.0 5.1 4.2
##  [397] 4.7 3.9 4.4 4.8 5.5 5.3 4.7 4.3 4.3 5.3 4.0 4.9 4.9 3.5 3.7 4.1 4.8 4.2
##  [415] 5.6 4.8 3.4 4.5 6.5 2.4 6.1 5.5 5.0 3.9 4.8 5.1 6.5 2.4 4.1 3.1 2.5 4.0
##  [433] 4.8 5.5 4.4 3.6 4.2 4.2 3.8 5.5 1.6 5.1 2.5 5.0 3.8 4.3 4.2 5.0 5.6 3.2
##  [451] 4.5 3.9 6.2 3.0 4.9 3.7 6.2 3.5 4.2 5.4 5.2 5.0 3.9 3.6 6.3 4.7 3.7 4.6
##  [469] 3.9 3.5 5.0 4.2 4.0 4.9 4.2 4.5 2.4 3.8 5.6 3.2 2.5 3.6 5.0 5.7 5.1 4.8
##  [487] 4.0 4.7 5.7 5.6 4.5 5.5 5.6 3.9 4.7 3.5 5.3 4.2 5.8 4.7 4.0 4.5 5.8 3.4
##  [505] 5.3 4.3 5.0 3.3 3.7 3.4 4.0 3.6 4.2 3.1 3.4 5.0 4.3 5.5 6.4 3.2 6.2 5.2
##  [523] 3.7 3.4 4.7 4.6 5.0 4.7 3.1 4.3 4.6 4.3 4.5 5.5 6.0 4.7 4.5 3.4 4.0 4.5
##  [541] 6.1 5.1 4.9 5.3 3.9 5.2 4.2 4.8 3.6 4.6 3.8 2.7 6.4 3.6 3.3 4.2 4.7 5.3
##  [559] 3.9 4.0 3.7 3.5 5.5 3.2 6.5 4.9 2.9 3.4 2.2 5.5 4.4 3.7 4.5 3.6 4.7 4.2
##  [577] 5.3 4.8 4.9 4.0 5.8 4.9 4.7 6.0 2.4 4.0 3.4 4.4 4.5 4.1 3.9 3.8 5.4 4.2
##  [595] 5.0 5.7 4.2 2.6 5.3 4.3 5.7 4.2 3.2 4.6 6.2 3.7 5.6 5.7 5.1 4.0 5.0 4.0
##  [613] 3.9 4.2 4.8 4.7 4.8 3.7 5.6 4.5 5.2 3.4 5.2 6.2 5.2 4.8 3.4 3.6 4.2 5.0
##  [631] 3.8 4.7 4.1 5.4 4.3 5.1 1.5 5.4 3.7 4.3 4.6 5.0 4.1 4.9 3.8 5.3 4.2 4.5
##  [649] 3.9 4.8 4.1 4.0 3.6 4.3 5.4 4.5 4.0 3.6 4.2 5.3 5.3 5.7 4.8 3.6 3.7 4.8
##  [667] 4.6 3.2 5.0 4.3 3.0 6.0 5.0 3.1 5.6 3.9 6.5 3.3 4.2 4.4 3.7 3.5 5.0 4.4
##  [685] 3.7 4.8 4.6 5.7 4.1 4.4 5.3 5.2 5.8 3.5 4.5 4.9 5.3 5.9 4.4 6.1 4.1 4.9
##  [703] 3.1 3.4 4.1 4.7 3.0 4.8 2.6 3.4 4.0 5.0 4.8 3.9 4.2 4.6 5.0 4.4 2.9 5.0
##  [721] 5.6 6.0 4.4 3.6 5.8 5.1 3.9 4.3 6.6 2.5 5.9 4.0 5.0 5.0 5.1 3.4 4.6 4.6
##  [739] 4.3 4.3 4.2 3.3 2.6 4.9 4.1 4.0 4.0 4.4 5.2 3.8 3.8 4.8 5.1 4.5 2.5 4.6
##  [757] 3.4 4.2 4.8 3.0 3.4 4.7 4.1 3.5 4.1 5.7 5.9 3.7 5.4 5.1 4.1 3.9 3.1 4.1
##  [775] 4.6 5.0 6.3 3.8 4.7 4.1 6.1 5.9 3.5 4.4 3.6 4.7 3.8 5.6 4.8 3.0 6.2 5.4
##  [793] 4.5 4.4 4.8 3.6 4.4 4.5 4.4 5.6 4.4 4.2 4.9 4.9 3.5 5.3 4.0 4.3 6.3 3.7
##  [811] 5.0 5.1 4.6 3.8 5.7 3.5 4.0 3.7 3.2 4.3 3.4 3.5 4.2 3.0 2.8 4.6 3.9 3.5
##  [829] 5.0 4.7 4.1 5.6 4.5 5.2 4.4 5.1 5.3 4.9 4.5 4.4 5.5 5.4 4.2 5.5 4.7 4.7
##  [847] 5.2 4.7 3.3 4.9 5.6 4.5 4.7 4.4 5.0 5.6 3.4 3.8 6.2 6.8 4.5 3.8 5.3 3.6
##  [865] 5.2 4.4 4.2 3.6 5.1 5.6 5.0 4.3 4.0 4.4 4.3 6.7 5.4 4.2 5.1 7.0 4.6 6.2
##  [883] 4.4 1.9 5.4 4.7 4.4 3.9 2.6 3.6 3.9 5.6 2.6 4.3 5.5 5.7 4.3 6.0 3.3 3.3
##  [901] 5.2 5.0 6.1 4.1 4.7 3.2 5.8 4.6 3.0 4.5 4.2 4.4 3.4 5.6 3.9 3.8 3.4 4.1
##  [919] 2.8 4.9 3.1 5.1 6.0 4.0 4.3 4.0 3.5 6.1 7.2 3.5 4.0 4.2 5.1 7.0 3.7 4.1
##  [937] 4.6 4.8 5.5 5.3 4.4 4.4 4.8 5.3 4.3 4.5 5.0 5.2 4.1 5.0 4.9 5.9 3.7 4.9
##  [955] 3.0 5.6 4.5 4.6 3.9 5.6 5.4 4.4 4.6 4.6 4.9 6.2 5.1 5.6 2.9 3.8 2.7 5.1
##  [973] 5.9 5.0 3.8 4.4 4.5 3.5 5.4 5.3 4.6 3.4 4.2 4.3 4.7 5.7 5.7 3.1 3.5 2.7
##  [991] 2.9 4.1 4.6 5.0 4.5 4.5 5.1 5.5 3.7 3.6
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
##   2.7   6.2
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
##    [1] 3.0 5.1 4.0 5.0 3.5 5.4 4.9 4.1 3.6 4.9 4.1 4.2 4.0 5.0 2.8 4.6 5.0 3.4
##   [19] 5.9 5.0 3.9 5.5 3.7 5.3 4.7 4.2 3.8 5.7 4.2 4.3 2.8 4.7 2.8 5.3 3.9 4.8
##   [37] 2.8 5.2 6.1 4.0 4.7 4.4 4.2 3.4 5.2 4.5 4.6 6.1 3.6 4.3 4.8 4.4 5.1 4.2
##   [55] 2.3 2.9 5.7 4.8 3.8 3.4 3.8 5.2 2.9 4.2 3.1 4.1 3.6 3.2 3.9 5.8 5.2 1.8
##   [73] 4.6 4.9 4.4 1.9 4.2 2.9 3.9 4.6 5.4 4.3 3.2 3.4 4.6 2.6 3.3 4.6 4.2 4.5
##   [91] 4.9 3.1 2.4 3.8 4.5 4.9 3.6 5.6 4.6 3.9 3.7 4.9 3.8 3.1 4.1 4.0 6.1 3.6
##  [109] 3.3 6.4 5.7 5.1 5.0 4.9 3.6 5.2 6.1 3.3 5.1 4.4 5.4 5.3 4.8 4.2 5.1 3.8
##  [127] 3.6 4.8 5.1 4.3 3.7 5.5 2.7 3.2 5.0 3.2 3.3 3.8 4.8 5.1 4.1 2.8 5.8 5.1
##  [145] 5.2 4.6 6.3 4.9 4.7 3.7 5.4 4.4 5.9 5.1 3.1 4.8 4.4 3.9 5.0 4.0 5.1 6.4
##  [163] 4.2 3.8 3.9 5.4 2.7 0.9 2.7 5.1 3.6 5.0 4.6 5.6 5.0 3.7 5.0 6.3 3.9 4.2
##  [181] 3.5 5.9 4.2 4.1 3.8 4.2 4.6 2.5 4.5 3.4 5.0 4.7 6.6 3.8 3.8 2.0 4.5 6.0
##  [199] 3.9 4.0 5.6 4.7 5.0 3.4 4.0 3.8 4.3 5.5 4.7 3.0 5.6 4.9 3.0 3.8 3.6 5.0
##  [217] 5.2 2.3 3.4 6.3 5.0 6.2 4.0 3.6 4.6 6.3 3.8 3.7 3.1 5.1 5.0 4.7 4.0 5.4
##  [235] 4.6 5.9 5.4 3.4 4.9 3.9 4.2 5.1 2.8 4.2 4.2 5.2 4.3 3.9 4.0 4.6 4.0 4.7
##  [253] 4.7 3.7 4.8 5.4 4.2 5.1 4.3 4.1 5.4 6.7 2.8 4.5 3.8 5.0 5.4 3.4 7.1 3.0
##  [271] 4.9 4.5 6.5 5.4 4.8 5.3 4.8 4.7 4.6 4.8 3.8 5.0 3.3 4.6 4.7 4.0 4.9 5.7
##  [289] 3.0 4.4 3.9 5.3 5.5 5.1 2.8 3.1 4.5 4.3 3.8 5.1 2.5 4.0 5.7 5.4 2.9 4.2
##  [307] 4.0 4.2 4.7 3.2 3.1 5.7 4.3 4.6 3.8 4.2 5.2 4.0 6.5 5.2 5.1 4.2 4.8 4.0
##  [325] 4.6 5.4 4.0 4.8 4.6 4.9 5.3 6.6 4.9 4.7 3.8 4.7 5.1 3.9 5.2 7.0 4.8 5.7
##  [343] 3.8 3.4 5.3 4.3 4.9 4.1 3.6 3.5 5.4 4.2 4.2 5.2 4.5 4.9 5.4 3.8 3.9 3.3
##  [361] 5.5 3.6 4.4 4.9 4.8 6.6 2.5 4.3 5.7 5.9 5.7 4.8 4.8 5.6 3.8 2.9 4.8 4.7
##  [379] 2.4 3.6 5.4 4.8 4.0 4.7 2.9 5.2 4.2 4.0 4.9 4.7 3.8 4.2 4.5 4.9 6.1 3.0
##  [397] 4.8 4.7 3.6 2.8 3.6 4.1 4.9 4.6 5.8 4.9 4.2 3.2 4.3 3.6 4.1 3.8 4.2 3.9
##  [415] 3.6 5.1 4.9 5.4 5.0 4.7 5.7 4.2 4.4 4.1 3.7 5.4 4.7 5.4 4.9 5.4 4.8 3.9
##  [433] 5.7 5.7 3.5 5.8 5.0 3.5 3.6 4.4 4.8 4.4 6.4 4.3 5.6 4.6 3.5 3.2 3.3 5.0
##  [451] 5.4 4.0 5.8 3.4 4.9 6.8 4.9 4.9 4.2 6.0 3.7 4.1 3.5 4.5 5.1 6.1 3.1 4.4
##  [469] 4.8 4.4 3.2 4.6 4.4 4.4 4.2 4.5 3.4 5.5 4.2 3.5 4.7 3.5 4.1 5.2 3.8 2.5
##  [487] 4.7 3.6 6.4 3.1 6.1 2.8 3.0 3.5 5.2 4.9 5.5 4.1 4.9 2.7 5.0 5.4 4.9 2.4
##  [505] 4.2 4.6 6.6 3.9 6.7 4.4 5.7 4.5 4.6 5.9 4.0 6.5 3.2 4.3 4.9 5.2 5.3 6.0
##  [523] 5.1 5.3 5.6 3.3 4.3 3.6 4.5 4.9 4.3 5.8 4.7 3.1 3.7 4.3 4.5 2.0 5.6 3.5
##  [541] 4.7 4.3 4.8 4.7 4.2 4.5 4.5 4.9 4.3 6.1 4.1 6.0 4.9 3.3 3.6 5.2 5.3 4.7
##  [559] 5.4 3.8 5.2 4.8 4.9 3.8 5.0 4.6 4.3 5.5 5.2 4.1 5.7 4.0 4.0 4.0 3.3 4.2
##  [577] 3.9 4.8 5.6 3.9 5.0 3.5 4.3 4.1 4.4 5.8 4.5 5.1 3.3 5.0 4.1 4.6 4.1 4.1
##  [595] 5.6 3.1 4.4 4.9 3.7 5.4 6.5 5.1 4.5 4.7 5.2 6.4 3.2 3.6 4.6 5.3 3.0 5.0
##  [613] 4.3 3.8 4.2 3.7 4.8 4.5 4.1 3.7 4.8 5.0 4.9 3.4 4.3 4.8 5.0 3.8 5.5 4.3
##  [631] 5.3 5.1 3.6 5.6 4.4 5.2 3.7 2.9 7.4 4.7 6.4 4.5 6.2 4.5 3.6 5.6 4.3 4.4
##  [649] 4.7 5.2 3.5 4.4 3.5 2.3 2.7 3.0 6.0 3.2 3.8 3.6 5.0 4.7 4.2 3.4 4.4 5.4
##  [667] 4.5 5.2 4.9 5.3 4.7 5.0 3.5 3.7 3.3 4.6 5.4 3.7 6.4 5.1 5.7 4.3 2.7 3.2
##  [685] 3.4 4.6 3.2 3.9 3.6 4.2 4.0 5.3 3.5 5.0 4.5 4.2 5.3 4.9 4.5 3.5 3.7 6.4
##  [703] 6.2 5.4 4.3 4.5 4.1 4.8 5.0 5.2 4.5 3.5 3.2 4.7 5.1 3.7 4.7 3.7 2.9 4.7
##  [721] 6.1 3.9 5.0 3.3 4.4 5.6 4.9 4.3 6.0 5.0 3.6 5.0 5.1 7.0 5.4 4.9 4.1 3.1
##  [739] 4.5 4.3 3.7 3.3 5.0 4.8 2.6 5.5 2.8 6.7 2.7 4.1 4.6 4.0 4.3 4.5 5.0 4.6
##  [757] 5.3 3.7 4.6 3.3 5.0 2.1 6.4 2.5 4.3 4.4 4.8 3.3 6.0 4.6 4.8 4.0 3.8 4.5
##  [775] 3.9 4.9 3.8 5.9 5.3 4.0 4.7 4.3 3.9 5.7 4.4 4.0 5.8 4.6 5.8 4.6 5.1 4.2
##  [793] 4.5 3.9 5.7 4.3 5.6 5.6 5.5 5.5 4.7 3.9 4.0 5.3 4.0 4.2 3.3 5.5 4.3 4.2
##  [811] 4.4 4.4 4.5 4.9 4.4 5.7 4.9 4.1 3.2 4.0 4.0 4.9 4.7 3.3 6.1 7.4 4.3 4.5
##  [829] 4.2 3.5 3.9 4.5 4.7 5.0 5.2 5.6 4.9 2.8 3.9 2.2 4.6 4.8 4.5 3.8 5.8 5.9
##  [847] 6.9 4.8 4.0 3.9 4.6 4.0 4.8 3.7 6.7 4.7 5.0 3.1 4.4 4.3 4.6 4.9 3.8 5.0
##  [865] 4.9 3.4 5.1 3.8 4.2 4.3 2.4 5.6 5.6 4.3 5.6 4.2 3.9 4.3 5.0 4.6 5.9 1.5
##  [883] 4.5 3.6 4.4 4.2 4.5 3.8 3.8 4.9 5.0 4.1 5.2 3.6 4.5 3.6 3.8 3.4 4.9 4.4
##  [901] 4.0 5.8 4.9 5.2 3.9 5.3 3.8 5.2 4.1 4.9 4.2 3.4 4.7 4.8 3.7 4.2 3.8 4.3
##  [919] 4.9 4.1 3.0 4.5 5.3 3.7 4.5 4.0 4.4 5.5 5.2 3.5 5.3 4.6 4.3 3.5 5.0 5.2
##  [937] 6.2 4.5 4.9 4.2 4.9 4.7 5.3 3.5 5.2 6.1 4.3 4.0 4.4 5.5 3.5 5.3 5.5 4.2
##  [955] 5.4 3.9 3.9 4.2 4.3 5.1 4.3 5.0 4.2 3.5 5.4 3.6 4.1 6.0 3.9 4.7 2.6 4.8
##  [973] 4.5 4.3 3.8 4.0 3.4 5.0 4.2 6.0 3.5 5.1 5.7 5.0 4.6 4.6 5.5 5.0 4.0 3.6
##  [991] 5.2 4.9 5.3 4.9 5.0 3.7 4.5 3.7 4.5 4.2
## 
## $func.thetastar
## [1] -0.0254
## 
## $jack.boot.val
##  [1]  0.48885449  0.35480226  0.32735562  0.11842105  0.05902965 -0.05102041
##  [7] -0.27535817 -0.31312849 -0.41761364 -0.55470588
## 
## $jack.boot.se
## [1] 1.013437
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
##    [1] 4.2 4.5 4.4 5.6 4.2 5.1 4.8 4.2 4.4 4.5 5.4 3.9 4.3 4.8 3.5 3.0 4.3 4.2
##   [19] 5.3 4.5 4.3 3.2 5.0 4.1 4.5 6.6 3.6 2.9 4.5 4.6 3.9 4.4 4.9 4.7 5.4 5.2
##   [37] 5.0 6.0 4.8 3.7 3.7 5.6 4.4 4.2 5.7 5.1 4.7 5.0 3.5 5.6 4.7 3.6 2.7 4.8
##   [55] 3.8 3.1 3.6 4.4 4.9 3.4 3.0 2.2 4.1 3.5 2.9 5.3 3.9 5.8 3.0 4.8 6.0 3.6
##   [73] 3.8 3.8 6.1 5.0 4.4 5.1 4.8 4.5 5.9 4.7 6.0 3.7 5.9 4.5 3.3 4.1 5.0 4.8
##   [91] 4.6 3.0 5.0 3.8 5.7 4.3 5.1 5.0 6.3 4.9 3.8 3.6 5.5 5.8 3.9 4.3 5.5 5.8
##  [109] 3.9 3.2 4.8 3.1 4.7 3.9 4.4 5.8 5.9 4.5 4.6 3.3 6.1 4.3 5.9 5.3 4.6 3.1
##  [127] 6.2 4.0 5.7 4.5 4.1 4.4 4.4 4.1 4.4 3.2 5.5 4.3 4.9 5.7 3.9 2.8 5.1 5.8
##  [145] 5.1 5.0 4.1 3.0 4.3 4.0 4.7 4.6 4.2 4.9 4.5 4.5 5.5 3.4 4.6 4.5 4.6 4.9
##  [163] 4.9 4.6 4.8 3.8 3.2 5.0 5.3 3.7 4.9 6.4 5.0 3.6 5.7 5.6 4.1 3.6 4.9 5.3
##  [181] 4.4 4.8 4.7 4.5 4.2 6.1 5.4 3.6 4.7 4.8 3.7 4.9 5.1 5.9 4.3 3.7 4.8 4.0
##  [199] 3.8 3.9 2.5 4.6 5.3 1.1 4.2 3.7 4.6 5.7 5.6 6.6 4.9 5.0 5.2 5.0 4.6 4.0
##  [217] 3.2 5.3 5.9 3.2 5.2 3.3 3.0 2.4 4.6 2.7 6.1 6.2 5.2 4.3 4.9 3.6 4.4 4.7
##  [235] 5.2 4.5 5.6 3.6 4.0 4.4 4.5 4.6 4.2 3.5 4.1 5.1 4.3 4.5 4.6 4.2 5.7 4.4
##  [253] 5.0 3.4 5.2 5.6 5.0 3.5 3.7 5.0 4.0 4.9 5.5 3.9 3.3 2.5 4.9 3.7 5.6 5.0
##  [271] 4.7 3.5 5.2 5.6 4.4 4.8 3.7 4.9 4.1 5.0 4.8 4.0 4.0 5.0 5.1 5.6 3.0 3.9
##  [289] 6.8 4.6 3.9 4.1 4.8 4.1 4.3 5.3 6.1 3.4 5.5 5.0 4.6 3.7 4.7 3.4 4.5 3.5
##  [307] 5.3 4.5 4.6 4.3 4.6 4.0 4.6 5.1 5.4 4.3 4.2 5.0 4.7 3.8 4.1 4.1 3.9 4.7
##  [325] 7.0 3.8 5.4 5.0 5.3 4.5 5.2 5.2 3.9 4.3 5.7 2.7 5.0 6.1 5.2 4.6 3.3 4.4
##  [343] 3.9 4.8 3.6 4.5 4.4 4.7 3.7 4.3 4.7 3.4 4.7 4.5 3.0 3.7 6.2 4.5 5.4 4.8
##  [361] 5.5 5.2 5.2 4.8 6.5 4.8 4.8 3.2 5.7 3.8 4.3 4.9 3.6 2.5 3.9 4.5 4.1 3.6
##  [379] 5.9 4.1 4.5 4.3 5.3 3.9 3.5 5.1 4.4 3.6 4.7 3.6 5.1 4.1 7.4 5.4 3.7 5.2
##  [397] 4.8 4.4 5.3 3.7 5.0 3.6 5.6 5.8 4.3 3.0 5.9 4.9 4.1 3.4 4.7 5.3 3.7 5.5
##  [415] 4.4 4.7 5.6 4.0 4.3 2.7 4.8 4.5 4.6 3.6 4.4 4.4 3.8 5.4 3.6 4.7 4.8 3.0
##  [433] 3.9 5.8 3.6 4.5 3.7 5.0 4.6 5.5 6.2 4.7 4.4 3.3 5.9 4.4 4.2 3.4 5.4 4.1
##  [451] 5.8 4.2 6.1 5.1 5.4 5.6 4.5 5.3 4.7 5.5 5.4 5.1 3.5 3.9 4.2 3.4 6.0 3.5
##  [469] 4.3 3.6 3.1 3.6 5.1 3.9 4.4 4.4 5.1 5.7 3.8 5.9 5.5 4.3 3.2 4.2 4.1 3.8
##  [487] 5.0 3.2 4.7 4.5 3.3 4.5 5.1 4.7 4.7 4.7 4.2 4.5 6.1 5.0 4.8 3.5 4.0 3.8
##  [505] 5.2 3.9 4.9 5.3 4.0 4.4 4.9 4.8 4.0 4.1 4.5 5.0 5.0 3.6 4.7 5.1 2.8 5.5
##  [523] 5.3 4.1 5.1 4.8 5.0 4.5 5.1 4.0 6.2 5.3 3.5 3.9 3.5 4.8 4.1 4.9 3.8 4.0
##  [541] 4.2 4.1 4.5 5.3 4.1 4.0 3.2 5.1 4.9 4.2 7.3 4.2 4.8 3.2 5.6 4.9 5.1 4.0
##  [559] 4.7 4.8 3.7 4.9 5.3 5.6 4.9 4.5 5.0 4.3 3.8 3.3 6.3 4.2 4.7 5.5 3.6 4.3
##  [577] 5.2 6.0 2.6 3.4 3.7 2.6 5.6 5.1 5.9 4.8 5.2 6.5 4.3 4.8 3.6 3.5 4.3 4.3
##  [595] 5.9 5.1 5.6 5.5 3.2 4.9 3.2 4.3 4.5 5.0 5.9 3.9 5.3 3.4 2.9 5.2 4.5 4.2
##  [613] 3.1 5.7 4.1 4.1 3.7 4.1 4.1 5.1 3.5 4.7 5.4 4.2 4.6 4.6 5.7 4.4 4.8 4.7
##  [631] 3.3 4.5 4.6 2.2 5.6 5.6 4.9 4.7 4.1 4.1 4.3 4.5 5.0 4.5 5.4 3.6 3.4 4.2
##  [649] 4.7 4.7 4.4 4.3 5.0 5.7 4.0 2.6 5.9 6.0 5.4 5.4 4.2 5.4 6.0 5.0 5.8 3.3
##  [667] 4.4 3.5 3.1 3.9 4.0 5.2 4.3 5.5 4.6 3.1 4.3 5.2 4.1 4.3 4.7 5.1 5.3 5.6
##  [685] 3.2 5.0 3.6 4.1 2.5 5.5 5.1 4.7 4.9 3.9 5.7 4.3 3.6 3.4 3.9 4.4 4.4 4.5
##  [703] 5.3 4.9 5.1 5.4 5.7 4.6 3.4 3.7 3.9 4.2 4.1 5.5 4.7 5.0 4.5 4.6 5.3 4.4
##  [721] 5.6 3.6 4.2 4.1 4.2 5.5 5.0 6.3 4.0 5.4 3.8 4.9 5.5 4.5 3.0 3.5 3.8 4.8
##  [739] 5.8 4.6 6.0 4.9 4.2 5.3 3.8 4.6 4.9 5.5 3.3 5.1 4.9 4.8 5.9 4.3 4.0 4.9
##  [757] 4.9 4.8 3.0 4.5 5.1 3.7 3.9 5.0 2.7 4.3 4.7 4.9 5.3 4.3 4.6 5.1 6.0 5.2
##  [775] 2.9 2.5 4.6 5.4 5.0 3.9 4.6 5.2 4.8 4.0 4.1 5.9 5.3 3.3 5.1 4.4 3.4 4.6
##  [793] 4.5 5.0 3.2 4.3 4.7 4.3 4.6 2.8 5.1 4.3 4.5 3.0 4.9 4.0 3.9 5.6 6.2 4.3
##  [811] 5.0 5.4 4.8 4.4 5.5 5.0 5.8 3.2 5.3 4.7 2.6 4.2 4.7 5.3 5.0 4.1 3.3 5.7
##  [829] 4.5 4.3 1.6 5.8 4.8 3.9 4.5 4.7 3.6 4.7 5.5 4.8 5.1 4.9 3.3 4.3 3.3 5.7
##  [847] 3.0 5.9 3.0 3.7 4.1 6.2 5.2 4.3 3.6 4.9 5.9 3.9 6.2 5.0 2.8 4.7 5.1 5.9
##  [865] 3.2 4.3 4.8 5.6 3.8 3.3 3.9 6.4 5.6 3.8 5.1 3.0 3.6 4.2 5.8 5.9 4.3 6.1
##  [883] 4.6 4.7 4.9 5.1 4.6 3.9 4.3 4.1 3.5 4.7 5.8 4.6 4.6 4.0 4.3 4.6 6.3 4.8
##  [901] 3.7 5.7 4.3 6.0 3.1 3.3 5.0 5.9 4.6 3.9 3.6 3.4 4.2 3.9 4.9 4.8 4.0 4.8
##  [919] 5.6 4.8 5.2 4.6 4.2 3.8 3.3 4.8 3.9 4.1 4.3 4.3 4.2 5.0 2.9 4.3 5.5 3.8
##  [937] 3.8 4.8 2.5 5.0 5.9 5.2 5.0 4.9 4.3 5.4 3.9 4.9 3.8 4.2 3.5 5.0 4.2 5.1
##  [955] 5.1 4.6 5.0 4.6 3.3 5.4 5.0 3.4 4.3 4.7 5.2 4.6 4.1 4.4 4.3 4.4 5.3 3.5
##  [973] 2.6 5.0 3.5 3.5 4.8 5.6 3.9 3.7 4.9 3.4 3.9 4.9 4.1 5.4 4.0 6.7 4.1 5.5
##  [991] 4.9 3.4 4.1 5.0 6.2 4.1 6.2 5.7 4.7 5.1
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.4 5.4 5.3 5.2 5.1 5.0 5.0 4.8 4.7 4.5
## 
## $jack.boot.se
## [1] 0.861162
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
## [1] 1.151162
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
##   3.226127   3.782111 
##  (1.374713) (1.743771)
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
## [1]  0.04467142 -0.06038342  0.68033845  0.27125790  1.18053674  0.02503591
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
##    [1]  1.550294883  1.645479788  1.389873320  0.429817090  0.389678087
##    [6]  1.086398453  1.263311899  0.564816775 -0.100639275  1.523294750
##   [11]  1.000101391  1.580429166  0.672143511  0.911008180  1.043010196
##   [16]  1.877812312  1.354849510  1.455345928  1.041440779  1.274627134
##   [21]  0.735920556  0.442108022  0.936599358  1.854913825  0.719214798
##   [26]  1.162784506  0.788207077  1.303341300  0.221380758  1.481086850
##   [31]  0.922827188  1.088287796  1.174767408  0.659890585  0.889441940
##   [36]  0.452117411  1.047287177  0.514090382  0.665257365  0.683744316
##   [41] -0.620709321  0.831908943  0.496601949  1.458760781  0.520246859
##   [46]  1.205476981  0.951470173  0.746164733  0.678835737  0.192503148
##   [51]  1.113836700  1.078609781  0.912697113  0.719388466  0.195977232
##   [56]  0.380997330  1.985099338  1.576777237 -0.052294928  0.838076068
##   [61]  0.813488829  1.403591216  0.553303118  0.580703174  0.985715919
##   [66]  0.990985024  0.477556939 -0.205976662  0.262979163  0.348698693
##   [71]  1.450784247  0.464349857  1.760233701  0.411398440  0.484340376
##   [76]  1.408734156  0.860545911  0.758944253 -0.064862045 -0.099303662
##   [81]  0.782484032  0.320444019 -0.210911351  0.868242909  0.138026590
##   [86] -0.882725422  0.909839186  0.985428825  1.940243324  1.199817052
##   [91] -0.150097579  2.367093217  1.484194652  1.991869734  0.211813860
##   [96]  1.255249473  0.926798492  0.663670084  1.574696147  0.882039476
##  [101]  0.524625349  0.322270849  0.787596509  0.430921984  0.290741466
##  [106]  0.644164605  1.747439316  0.520269227  1.725090650  1.671128386
##  [111]  1.481752136  0.884425841  0.233655817  0.946601404  1.193268670
##  [116] -0.279176635  1.604759156  0.893640737  1.047575874  0.937752538
##  [121]  0.560295252  0.376891780  0.414811075  0.613650171  1.227798490
##  [126]  1.101083976  0.665221803  0.389936910  0.836014447  0.051535697
##  [131]  0.157346529  1.775226237  0.807594230  0.961316245  1.507108207
##  [136] -0.036531154  1.159780246  0.292891887  1.363092580  0.848267036
##  [141]  1.001857268  0.463577893  1.651348133  0.206121348  0.654317481
##  [146]  0.910693903  2.031201431  0.902255089  2.097059386  1.082435146
##  [151]  1.356060677  0.589220910  0.539525867  1.145038944  0.845882565
##  [156]  0.822239582  0.356449314  0.246939652  1.215608975  0.229172092
##  [161] -0.062164940  0.779332298  0.864898255  0.174615558  0.433557691
##  [166]  1.292096998  2.129495712  0.542132294  0.509946110  0.488923365
##  [171]  0.520344316 -0.021109853  1.510025399  0.229730036  0.668369151
##  [176]  0.259386100  0.223655490  1.634540511  2.116162068  0.309420521
##  [181]  0.919542033 -0.380800366  1.801538563  0.986774693  0.967929999
##  [186]  1.524294985  0.079678668  0.834736605  1.512186416  0.125537237
##  [191]  2.161987909 -0.177797482  0.342441195  0.152548800  0.654169948
##  [196]  0.779867879  1.054346342  0.024846992  0.525495883  0.586016001
##  [201] -0.530918934  0.700463037  0.055776021  0.646471409  1.560746841
##  [206]  1.121725884  1.054345290  1.069921845  1.057842739  0.380816765
##  [211]  1.165019746  1.153791574  2.276621177  0.854269368  0.558234210
##  [216]  1.018249899  0.706557163  0.300048461  0.870757514  0.286861043
##  [221]  0.066126452 -0.834989968  0.584494194  2.150663850 -0.289041923
##  [226]  0.542511294  0.356895081  0.198186783  0.429259018  1.054501304
##  [231]  0.627777370  0.755164460  0.566488591  0.963052196  0.502059473
##  [236]  0.700497343  0.510641746  0.456616042  0.234828003  1.413142507
##  [241]  1.262496238  0.498550007  0.623445971  1.869180686  0.573078726
##  [246]  1.025680347  1.285756484  0.594775454  1.649185408  1.309589451
##  [251]  1.080769956  1.249423482  0.443109927  1.065733182  0.253661712
##  [256]  0.279804341  0.826007641  1.979446760  0.824782630  1.475065861
##  [261]  0.778685264  0.105837444  0.267984545  0.740122568  0.587076308
##  [266]  0.062217981  0.196446029  0.587028816  1.201711534  1.203654684
##  [271]  0.346475127 -1.031184025  0.706801828  1.373418264  0.932335549
##  [276] -0.004887053  0.179451350  0.723972112  1.042395674  0.614647351
##  [281]  0.168323131  1.234179995  0.138538672  0.193037591  0.223453062
##  [286]  1.611097442  0.924652015 -0.269770265  1.670776994  1.037510859
##  [291]  1.083660361  0.818250834  0.783245694  1.256373311  0.776764735
##  [296]  0.113213720  0.776388995  1.034371304  1.501755741  0.627314826
##  [301]  1.545840748 -0.279687092  1.751440496  1.557375579  1.848884933
##  [306]  0.546440809  0.746736986  0.284997818 -0.433790167  1.313789870
##  [311]  0.797617689  1.179115597 -0.055424610 -0.321375104  0.815474251
##  [316]  1.395062575  0.085819369  0.430184377  0.577696345  1.076073804
##  [321]  0.241484728  1.159296706  0.443025486  1.220689182  1.440607420
##  [326]  1.136298098  0.794327230  2.138134614  0.743568670  0.471367436
##  [331]  0.400143939  1.170003644  0.103208894  0.540122204  1.349555552
##  [336]  0.611961139  1.151816220  0.744186836  0.174640509  0.736709590
##  [341]  0.529718885  0.250767000  1.077565034  0.647228456  1.207630998
##  [346]  0.689671684  0.226842540  1.530452838  0.218778668  0.337422233
##  [351]  1.147727180  0.342797013  0.837240350  1.470381466  0.742280674
##  [356]  2.408462899  0.633265931  1.225709928  0.337422233 -0.435424108
##  [361]  1.186521875  0.742280674  0.284153014  0.683191332  0.958132087
##  [366]  1.181407037  1.095110160  0.469757717  1.105688424  0.294489804
##  [371]  0.348044530  1.112578452  1.267456001  0.057763635  0.932222837
##  [376]  0.535313090  1.391908233 -0.019848560  1.185663326  0.569971423
##  [381]  0.301833587  0.244943838  1.520105680  1.124065153  1.083813162
##  [386]  1.224736014  1.330874083  0.143815764  0.786802602  1.609956834
##  [391]  0.408961048  1.596350803  0.353325404  0.631199090  1.651560645
##  [396]  0.238800886  0.209996001  0.887370624  1.042314968  0.729194473
##  [401]  1.386108045  1.101082411  1.005073711  1.187922561  0.872585646
##  [406]  0.627872947  2.297134004  0.593613002  0.381501397 -1.059987425
##  [411]  0.737582500  1.199692995  0.872172132 -0.032081616  1.701129195
##  [416]  0.448621773  0.907014675  0.798435313  0.996899106  1.002216653
##  [421]  0.107244151  0.362341155  1.556901242  0.507094712  0.788377959
##  [426]  1.473647964  0.702363679  0.626152007  1.599014146 -0.758062408
##  [431]  0.778679381  2.022457290 -0.107142307  1.201119343  0.675977488
##  [436] -1.075702906  0.212604500  0.672296289  0.790099466  1.186676721
##  [441]  0.905549398  0.909881623  0.222168069  1.055139657  0.126758345
##  [446]  1.028537237  0.259018533  1.133904472 -0.362891257  0.664149928
##  [451]  0.813959809  0.300705835  0.675212956  1.207328219  0.102046911
##  [456]  1.199432558  0.714856109  0.520912343  0.815474251  0.257369137
##  [461]  0.292641533  0.867552745  0.599643279 -0.041569711  0.899925503
##  [466]  0.714122226  1.903358746  0.756197922  1.316161605 -1.485985058
##  [471]  0.722378552  1.445216406  0.434017237  0.177476970  1.003883006
##  [476]  1.503107431  0.550784367 -0.582414226  0.552388236 -0.005885250
##  [481]  1.223775355  0.498019156  1.176980597  1.053599207 -0.013852254
##  [486] -0.175803442  0.366780615  1.548275941  0.203614524  0.965214551
##  [491]  1.284496925  0.474633526  0.497020527  0.738379265  0.113706835
##  [496]  0.351069921  1.273993579  0.319358763  0.863028101  1.143538033
##  [501]  0.807087401  0.643954230  1.937517123  0.665135179  0.365592817
##  [506]  0.625365433  0.706801828  0.493409893  0.710152107  1.337974834
##  [511]  0.559026811  1.086951626  2.233841167  0.463557569  0.572758539
##  [516] -0.047120773 -0.021964544  0.553554594 -0.644145756  0.977029088
##  [521]  0.158794176  0.493736383  1.086463173  0.927115873  2.328431389
##  [526]  0.673419440  0.105903775  0.013326865  0.579266812  0.515592109
##  [531]  1.230353116  0.835935103  1.294476999  0.858677030  0.184958409
##  [536]  0.740612667  0.769737915  1.131678009  0.505492967  0.969652580
##  [541]  0.186382609  0.488971093  0.049448407  0.774934632  0.496642918
##  [546]  0.573078726  0.918186363  0.203492682  1.857680875  1.234458964
##  [551]  1.001102598 -0.171850793  1.042180492  0.971852773  0.884632183
##  [556]  0.529363477  0.324745815  0.802339924  0.938996142  2.027818538
##  [561]  0.434629568  0.466258002 -0.324283168  1.633383545  0.305899895
##  [566]  0.597939024 -0.201561148  0.205757832  1.277403596  1.409950002
##  [571]  0.747725563  1.451557961  1.451603260  0.637956629  1.090443993
##  [576]  0.742861278  0.395335663  1.412273018  1.568133018 -0.723626342
##  [581]  0.589976019  0.594438582  0.167615950  0.531876602  0.758734146
##  [586]  0.638227420  1.511365593  1.199287724  1.828411235  0.426145879
##  [591]  1.114685467  1.207328219  0.488971093  0.415940075  2.121088715
##  [596]  0.400730265  0.949364048  1.899792353  1.273381476  1.186414134
##  [601]  0.586016001 -0.306909350  0.024066215  0.334782896  1.024770903
##  [606]  0.889419851  1.186696702  0.623317118  0.909287490  1.432169140
##  [611]  0.549263316 -0.141542435  0.442452051  1.243840366  0.289814656
##  [616]  0.452153662  0.529000588 -0.259470493 -0.210528578  1.185886445
##  [621]  1.027295718  1.056009462  0.077011707  0.820811248  0.889114250
##  [626]  1.574983077  0.870397189  0.368310605  0.908142090  0.224145063
##  [631]  1.090720139  1.442623548  0.926928015  0.353073370  0.396606942
##  [636]  0.406466559  0.793624449  1.476101995  0.216370024  0.416270070
##  [641]  0.540399168  0.441440019  0.646190694  0.475071829  1.237471760
##  [646]  1.239886141  0.164642274  1.273454170  1.369893065  1.117595493
##  [651]  1.096724036 -0.607340011  1.383483514  1.193745273  1.217895857
##  [656]  0.508815019  1.596024409  1.395227651  0.222314339 -0.532896253
##  [661]  0.552081443  0.350723011  1.040676732  0.903752262  1.201119343
##  [666] -0.416180375  0.890650139  0.607587827  0.811240671  0.662694766
##  [671]  0.072099913  1.333504745  0.402640441  1.706300471  0.444852578
##  [676]  1.207497467  0.709659327 -0.144181486  1.343293268 -0.693293895
##  [681]  0.185972331  1.692692341  0.198850710  0.901230136  1.089965430
##  [686]  0.507557662  0.640620371  0.918258685  0.364069792  1.012306511
##  [691]  0.961843844  0.190351411 -1.171781954  1.124185631 -0.170223748
##  [696]  1.506992520  1.574696147  2.008325113  0.817915927  2.219041158
##  [701] -0.560907431  1.480027439  0.365457275  0.242837574  0.361138477
##  [706]  2.114049137  0.996240683  1.565197352  0.821843684 -0.279887404
##  [711]  0.818524645  0.634198165  1.820849289  0.117266097  0.650408239
##  [716]  1.191242356  1.407186614  1.143538033  1.024517543  0.664260126
##  [721]  0.452633798  1.334016157  1.897679645  1.015942305  0.902298928
##  [726]  1.506213016  0.806776021  0.893150926  1.634477712  1.095577542
##  [731]  0.402128132  1.249126373  1.405386105  1.241066588  1.594577505
##  [736]  0.780722772 -0.012572855  0.329790749  0.742505405  1.233033296
##  [741]  0.723566954  0.788513416  1.323968135  0.150116313  0.275264560
##  [746]  0.780863149  1.326816657  0.081398714  1.539769779  1.965447867
##  [751]  1.230383503  0.578790597  1.641089409  1.329220417  0.532710745
##  [756]  0.659677973  1.307103086 -0.019914018  0.265052193  1.192333950
##  [761]  0.715180456  0.856517846  1.320329119  0.034317968  1.008127019
##  [766]  1.026956743 -0.313922561  1.127911538  0.728777067 -0.335611577
##  [771]  0.553276910  0.150104383  1.620738080  0.305055214  0.769058768
##  [776]  0.924400042  0.255419013  0.455342741  1.940238414  0.569895191
##  [781]  1.902584894  0.749130188  1.589865181  0.963052196  1.139005356
##  [786]  0.794490835  0.500723183  1.740657869  0.668621262  0.741913181
##  [791]  1.678235516  1.444816183  0.535600099  0.236203302  1.579017811
##  [796]  0.233310777  1.342037298  1.008345566  0.588751934  0.569310602
##  [801]  1.534733154  0.410162499  0.376808372  1.205337429  0.707265701
##  [806]  0.564334336 -0.288375239  0.554292271  0.029504380 -0.046049755
##  [811]  1.555281421  1.315475602  0.494977125 -0.447915789  1.315868550
##  [816]  0.331608600  0.657745852  0.819599047  0.677862662  0.207778621
##  [821]  1.668503415  0.203745888  0.746097239  0.958283262  0.405561903
##  [826]  0.872585646  1.049121279  1.368355963  1.571457628  0.604045305
##  [831]  0.369137075  0.650977526  0.433365843  0.419650404  0.899269578
##  [836] -0.024600600  1.216079786 -0.173602787  1.187479245  0.616118392
##  [841]  0.028191650  0.958887464  0.443444435  0.524675691  2.084726023
##  [846]  0.783900864  0.822795059  1.270584571  1.988520296 -0.247954557
##  [851]  0.484765449  0.953588623  0.634731623  0.066116597  1.362916146
##  [856]  0.993784943  0.722355483  0.654913046  0.375200474  1.584987338
##  [861]  1.377318046  1.374280354  0.737809980  0.736080771  0.275008815
##  [866]  1.307368474  1.606318795  0.981698222  0.743410361  1.501029893
##  [871]  0.982150169  1.325995437  0.646751311  0.244993716  0.913171143
##  [876] -0.073104535  1.279608869  1.764610142  0.010490099  0.631426619
##  [881]  0.882369536  1.029040639 -0.223820293  0.937841212  0.404628581
##  [886]  0.097385882  1.495822570 -0.188116670  0.952675364  0.993456359
##  [891]  1.098466400  1.265754017  0.588656544 -0.061102689  0.814701546
##  [896]  1.151162363  0.181058501  0.615787821  0.302055844  1.539053316
##  [901]  0.822128834  0.256202344  1.590804871  0.916617438 -1.371565545
##  [906]  0.689384523  0.122988803  1.269959409  0.853627065  0.497290105
##  [911]  0.870785863  0.889934353  2.061276201  1.457615483  0.694585636
##  [916]  0.372263661  1.805633125  0.543248982  0.245245138  1.795762094
##  [921]  1.771115762  0.094918245  0.412049307  1.284404622  0.504036703
##  [926]  0.123856768  1.453328196  1.734871907  1.051761429  1.960238388
##  [931]  0.496422290  0.821347448  1.252991244  0.773526176 -0.085828157
##  [936]  0.029196155  0.638369858  0.698721043  0.817233539  1.379736433
##  [941]  0.560850047  1.179799675  1.221551845  0.749900189  1.234833560
##  [946]  0.594669357  1.246188088  1.030384607 -0.312087331  2.108889456
##  [951]  0.363765173  0.229522994  1.728357196  0.945814921  0.404855770
##  [956]  0.041546748  1.349720474 -0.105602469  1.002203865  0.468069030
##  [961]  0.769901506  1.105378493  0.858677030  0.907716051  0.523370720
##  [966]  0.765434261  0.047381442  0.740663669  0.305976321  1.255965849
##  [971]  1.289804020  0.973645129  0.516047208  0.588041675  0.942814224
##  [976] -0.062600630  0.662875213  1.441119005 -0.299627572  0.661247473
##  [981]  0.728289350  0.551951742  0.832784751  1.026098546  1.546991000
##  [986]  0.342115969  0.839164073  1.037426512  0.356867607  1.943551528
##  [991]  0.605675030  1.686177922  1.599014146  0.715158360  0.738329388
##  [996]  0.709299275  1.177575100  1.020452975  0.744175134  1.925618119
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
##      mean         sd    
##   0.8529965   0.5061640 
##  (0.1600631) (0.1131799)
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
## [1] -0.2501182 -0.2376035 -0.0081845  0.9798824 -0.3392879  0.4289510
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
## [1] -0.0075
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9237896
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
## t1*      4.5 -0.03953954   0.9474405
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 2 3 4 5 6 7 8 
## 1 1 1 2 1 1 2 1
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
## [1] 0.0464
```

```r
se.boot
```

```
## [1] 0.9329617
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

