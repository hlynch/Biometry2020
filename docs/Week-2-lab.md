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
## 0 1 4 5 6 7 8 
## 2 1 1 1 1 3 1
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
## [1] -0.0173
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
## [1] 2.686833
```

```r
UL.boot
```

```
## [1] 6.278567
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.6   6.2
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
##    [1] 4.3 4.5 3.8 5.3 4.2 4.3 5.8 3.0 6.0 4.3 4.4 3.9 4.1 6.6 5.6 3.1 5.6 5.0
##   [19] 4.5 5.2 3.5 6.8 4.3 4.5 3.6 3.7 4.3 4.1 4.2 4.2 6.0 3.8 3.4 3.1 3.9 4.3
##   [37] 4.3 6.8 5.3 5.4 4.7 3.4 4.5 2.8 4.7 5.5 4.6 5.9 3.4 5.9 3.5 3.8 5.0 3.7
##   [55] 4.6 4.9 4.2 3.8 4.4 4.5 3.3 4.7 5.3 4.3 5.5 3.4 2.7 3.8 5.4 4.2 3.8 4.9
##   [73] 4.9 4.3 4.3 4.0 5.0 3.7 5.6 4.6 5.6 4.5 4.0 5.3 5.1 5.0 5.6 3.5 4.6 4.8
##   [91] 4.9 5.9 4.5 4.8 4.1 3.0 3.7 3.9 4.6 4.2 4.3 4.0 3.0 4.0 3.7 4.8 4.7 4.3
##  [109] 3.3 6.7 3.9 5.3 4.3 3.5 3.7 4.3 4.4 5.4 5.8 3.8 3.7 3.5 4.1 4.7 4.8 4.0
##  [127] 3.7 6.3 5.1 4.0 4.8 4.1 3.3 5.7 4.8 5.7 5.6 3.3 4.9 3.7 3.5 4.2 3.8 5.2
##  [145] 3.5 3.8 6.1 4.3 3.9 4.3 3.5 2.8 4.9 4.0 4.2 3.9 4.7 3.7 2.4 6.0 5.1 3.3
##  [163] 5.4 4.0 4.5 3.5 4.6 5.6 4.3 4.9 4.5 4.4 4.1 3.2 4.6 4.6 6.0 4.0 4.9 5.6
##  [181] 6.1 3.1 3.7 5.0 7.0 5.7 3.6 4.6 3.4 6.0 4.7 4.4 6.7 5.4 4.3 3.8 5.2 5.2
##  [199] 6.0 4.5 4.1 5.2 4.0 6.4 3.0 3.6 5.6 5.4 4.8 6.3 5.2 5.0 3.8 4.1 5.8 6.5
##  [217] 5.0 3.2 4.8 4.6 5.3 4.6 3.7 3.1 4.1 5.5 4.4 1.9 5.6 2.5 3.2 3.6 5.2 5.0
##  [235] 4.0 5.2 3.9 4.2 5.6 4.9 5.6 5.2 3.7 3.1 4.5 4.8 5.3 5.2 3.8 2.8 4.6 4.6
##  [253] 4.2 4.1 4.4 5.3 3.4 4.8 5.6 3.8 5.6 6.6 3.7 4.6 5.4 4.9 7.0 3.2 2.4 3.9
##  [271] 4.8 3.7 4.9 5.5 3.6 3.0 4.5 5.1 5.0 4.6 5.3 4.5 3.5 3.7 6.1 4.1 3.4 5.5
##  [289] 4.9 3.8 4.4 3.7 4.8 4.3 3.5 5.5 3.0 1.8 4.2 4.7 5.4 3.7 5.4 5.1 4.2 5.2
##  [307] 2.1 5.2 4.3 4.6 3.9 3.9 3.7 4.7 4.0 4.8 4.2 3.3 3.6 4.4 3.7 4.7 4.2 5.8
##  [325] 4.4 5.3 4.9 5.9 4.7 3.1 4.1 2.8 4.2 5.1 5.6 5.6 5.2 5.1 2.5 5.1 3.2 6.2
##  [343] 3.2 3.5 5.0 4.9 5.4 5.7 3.9 5.3 5.2 4.1 5.9 3.9 4.0 4.5 4.7 5.1 4.9 4.1
##  [361] 4.4 4.6 3.3 4.6 4.4 4.5 5.9 5.9 4.3 4.1 4.1 4.3 5.5 4.1 4.8 5.0 4.8 5.4
##  [379] 4.3 4.3 4.3 4.4 3.6 4.7 4.4 5.3 5.9 4.5 3.7 6.9 6.3 4.3 4.1 2.5 2.6 5.2
##  [397] 3.3 4.4 5.0 4.8 3.2 4.0 4.9 4.6 5.6 5.1 4.0 3.6 4.4 3.6 5.2 4.4 2.2 3.4
##  [415] 5.3 3.4 4.7 3.2 4.7 4.1 4.1 4.9 4.9 4.5 3.9 3.9 4.2 4.6 4.4 5.7 5.3 5.7
##  [433] 5.1 6.6 3.3 6.0 4.2 5.1 4.7 5.0 4.5 5.4 3.2 3.2 4.8 5.4 5.0 4.0 3.9 4.3
##  [451] 3.9 6.1 4.9 4.2 6.1 2.6 3.8 4.9 3.5 4.9 3.7 4.7 4.4 5.1 4.6 2.2 4.1 4.0
##  [469] 4.5 4.5 4.8 3.0 5.6 4.2 6.3 4.8 5.4 4.8 5.3 3.7 4.6 4.7 5.0 2.9 4.7 3.9
##  [487] 5.8 4.5 2.4 5.1 5.2 4.4 4.4 2.6 4.8 3.7 5.0 3.8 3.7 4.3 4.7 6.3 4.4 3.4
##  [505] 5.4 5.0 2.4 3.6 4.1 3.9 4.0 4.3 4.2 4.3 3.1 4.6 6.0 4.5 4.0 3.5 4.5 3.3
##  [523] 5.3 4.5 4.5 4.2 5.9 5.0 5.4 3.0 5.8 5.8 3.5 2.9 4.1 3.0 5.9 4.3 3.8 4.8
##  [541] 2.9 4.0 3.5 4.6 4.0 4.6 5.3 3.4 5.0 5.0 3.4 3.5 4.8 5.3 4.0 2.9 4.9 4.9
##  [559] 5.1 3.7 5.2 5.6 5.6 4.6 5.1 4.7 3.5 4.7 4.4 4.7 3.9 3.1 3.1 4.7 2.9 3.3
##  [577] 4.6 3.0 4.6 4.1 4.4 4.9 5.7 4.1 4.6 6.3 5.9 3.6 5.9 4.3 4.9 4.5 4.7 3.6
##  [595] 5.3 4.5 3.6 4.9 3.7 4.7 3.8 3.7 3.2 4.5 3.8 5.7 3.8 4.1 4.5 5.0 5.7 4.5
##  [613] 2.7 5.1 3.7 4.9 4.4 6.6 4.2 4.9 4.5 4.4 4.5 4.7 3.8 4.7 4.7 4.4 3.9 3.2
##  [631] 4.4 5.7 5.0 2.9 5.3 3.1 4.1 5.1 5.8 5.4 3.8 3.8 4.1 5.8 4.9 5.5 3.0 5.9
##  [649] 3.5 3.9 6.4 4.2 4.9 3.4 4.6 4.6 4.5 4.7 4.6 6.3 2.0 5.2 5.2 4.0 3.5 7.5
##  [667] 3.6 4.3 6.1 4.6 3.3 5.1 3.7 4.4 3.6 4.1 3.0 4.5 3.7 5.3 4.8 5.7 4.9 5.9
##  [685] 4.8 4.4 4.5 3.5 4.3 5.4 4.9 4.4 4.4 5.3 3.3 4.7 4.9 4.6 4.4 4.4 4.1 4.9
##  [703] 5.2 3.1 2.2 4.6 4.2 4.2 4.5 5.9 4.2 3.5 4.3 3.1 5.4 4.6 4.6 2.5 5.2 5.7
##  [721] 6.3 4.3 5.3 4.2 4.3 4.1 3.7 3.9 4.3 3.7 5.1 5.8 4.2 4.6 5.4 3.7 2.7 5.3
##  [739] 2.8 4.3 2.9 5.1 4.9 4.2 4.0 4.3 5.0 5.0 4.9 3.8 5.0 3.6 5.7 5.2 3.7 4.4
##  [757] 4.3 4.1 3.5 3.9 4.5 4.4 4.7 4.3 3.3 4.8 6.3 4.7 6.3 4.6 4.6 4.7 4.8 3.6
##  [775] 4.5 4.5 4.6 3.9 4.0 6.0 5.6 4.4 3.5 5.0 4.4 2.0 5.0 4.8 5.6 5.4 3.5 4.5
##  [793] 4.0 4.4 4.6 4.9 3.7 4.4 5.4 3.6 5.1 2.9 2.9 5.0 4.0 4.1 6.1 4.8 3.9 4.2
##  [811] 5.6 4.4 5.4 4.0 4.7 3.1 3.9 4.6 2.7 5.3 4.5 4.7 4.3 5.2 5.1 3.9 5.0 4.6
##  [829] 4.7 6.4 4.9 4.3 5.0 5.2 4.2 4.8 4.2 4.1 2.7 4.6 4.2 5.0 4.3 4.3 4.3 6.3
##  [847] 5.6 5.2 5.3 4.2 4.4 5.8 4.3 3.8 4.5 4.8 4.7 4.5 5.4 3.8 3.4 3.6 4.4 4.0
##  [865] 4.9 5.7 3.9 3.9 5.2 4.9 4.2 6.1 5.1 6.3 3.1 4.5 3.6 4.4 4.1 5.3 4.6 2.2
##  [883] 3.2 4.4 4.7 5.1 4.2 5.9 4.2 3.5 5.6 4.6 3.9 3.7 4.0 5.1 4.3 4.0 5.8 4.2
##  [901] 5.0 3.7 5.8 4.8 4.8 3.1 3.9 4.7 3.5 5.8 5.4 3.6 4.5 5.0 5.7 5.3 4.6 4.4
##  [919] 3.0 4.9 5.2 5.3 6.7 4.8 4.4 5.1 4.6 4.1 5.2 5.8 4.7 5.7 4.5 5.1 4.9 3.7
##  [937] 3.9 5.0 4.8 4.3 3.1 3.6 3.8 4.4 4.4 5.6 5.2 3.1 3.8 4.5 4.6 5.7 3.8 4.3
##  [955] 3.6 3.8 3.9 4.6 4.6 4.3 4.5 3.9 4.1 4.1 5.1 3.4 3.7 4.6 4.9 4.8 3.7 4.7
##  [973] 5.1 5.5 4.6 2.9 5.4 5.8 3.3 4.2 5.4 3.6 4.3 5.3 5.5 2.6 4.4 3.5 3.4 4.3
##  [991] 5.7 6.3 4.2 4.3 4.5 3.5 4.5 4.6 4.6 5.0
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
##    [1] 5.0 4.5 4.6 5.7 6.0 5.3 5.3 4.5 5.1 3.8 4.1 5.8 4.2 4.7 5.6 4.7 4.4 4.0
##   [19] 5.1 3.9 4.6 5.1 5.1 4.1 5.7 5.3 4.8 5.0 4.7 5.2 4.3 4.0 5.7 5.0 6.2 4.2
##   [37] 6.1 5.4 4.3 4.1 4.3 3.9 4.6 4.6 5.3 4.5 5.7 4.9 5.9 4.5 3.3 4.3 3.6 5.6
##   [55] 4.5 5.0 4.9 5.3 3.1 3.0 5.0 4.1 5.0 3.8 5.2 4.0 5.1 4.4 4.1 4.3 4.3 5.2
##   [73] 5.3 4.6 3.7 5.0 3.9 3.5 5.1 3.3 5.1 5.3 4.8 3.4 5.6 6.4 4.2 5.2 4.1 4.0
##   [91] 5.5 4.4 3.6 6.1 3.8 3.7 4.7 6.3 4.9 4.4 4.2 4.5 2.9 5.7 3.9 5.5 3.9 4.7
##  [109] 3.9 6.1 5.8 5.6 5.0 4.9 4.8 6.5 4.1 3.8 5.6 5.4 4.2 4.2 2.8 4.7 4.6 3.4
##  [127] 2.1 5.9 5.1 4.8 4.3 3.3 4.6 5.9 4.0 4.1 4.5 4.5 4.5 3.9 4.3 3.9 4.7 5.1
##  [145] 3.6 3.9 4.3 4.1 3.7 4.4 4.7 4.4 3.1 5.2 4.9 5.0 4.9 4.2 4.5 4.0 3.6 5.6
##  [163] 4.2 4.0 3.2 5.0 2.6 5.0 2.2 4.6 6.3 4.2 2.8 3.6 4.8 4.7 3.3 5.0 4.7 5.3
##  [181] 6.2 5.4 6.1 4.2 4.6 4.5 3.0 5.4 5.0 3.6 5.0 5.0 5.7 4.0 4.9 4.3 3.3 3.8
##  [199] 4.7 3.3 4.2 6.3 3.5 3.8 3.3 3.2 3.7 5.5 2.6 4.3 3.8 3.4 4.5 3.5 2.8 5.0
##  [217] 4.6 5.7 4.4 4.2 5.1 4.9 4.2 3.4 3.9 3.5 5.2 2.9 5.9 6.1 4.3 4.6 6.1 6.9
##  [235] 6.8 2.6 3.4 3.9 2.7 7.0 4.8 2.6 3.7 4.8 2.7 3.5 5.1 5.2 2.6 4.6 3.8 3.8
##  [253] 3.5 3.8 4.6 4.8 5.3 6.4 4.7 5.2 5.8 4.4 3.9 4.0 4.9 5.1 3.1 3.9 5.9 5.8
##  [271] 5.7 6.3 4.1 4.2 5.1 5.3 4.9 4.4 4.0 3.2 5.0 4.1 3.6 5.7 4.3 4.6 2.6 4.2
##  [289] 3.6 4.6 5.1 4.7 4.8 4.7 5.5 5.4 3.4 3.5 5.6 4.0 5.1 4.4 6.4 4.7 5.4 4.8
##  [307] 5.0 5.7 5.3 3.0 4.4 5.5 5.5 5.5 3.1 3.9 4.1 4.4 6.7 4.0 4.1 4.7 4.5 5.1
##  [325] 3.8 4.1 3.5 4.1 5.1 3.9 4.4 5.2 4.0 3.6 6.1 4.3 3.7 3.4 4.8 5.9 5.6 5.6
##  [343] 3.0 4.6 4.4 4.9 3.7 3.1 4.0 3.6 6.1 2.8 4.1 5.4 4.7 5.0 5.6 6.6 3.5 3.3
##  [361] 3.6 5.7 4.3 4.1 5.9 6.0 5.0 5.1 5.0 4.4 3.9 4.0 5.8 6.2 5.2 5.2 4.5 4.6
##  [379] 3.0 3.8 4.3 5.2 4.2 2.6 5.4 3.7 4.7 3.1 4.4 4.1 4.6 5.7 5.7 5.1 3.0 4.0
##  [397] 4.0 4.3 4.2 6.1 5.9 4.2 5.4 5.7 4.7 5.4 4.1 4.9 5.1 4.1 3.2 3.8 6.0 4.6
##  [415] 4.1 4.3 2.7 5.6 4.1 5.2 5.4 4.4 3.3 3.7 4.2 4.1 5.1 5.7 4.1 4.2 5.1 3.8
##  [433] 5.7 3.5 3.6 6.4 4.8 5.3 4.3 3.7 5.3 6.2 5.3 5.5 4.7 5.3 5.1 5.7 5.7 5.8
##  [451] 5.0 4.6 4.7 4.6 4.7 4.8 2.9 3.9 4.9 2.4 4.5 3.3 4.9 4.7 4.6 3.6 5.5 5.0
##  [469] 5.2 4.2 5.5 4.8 4.3 3.8 6.3 3.7 3.9 5.4 4.1 4.6 5.0 5.0 3.9 4.4 5.6 5.3
##  [487] 5.2 3.3 5.3 3.5 3.9 4.4 6.1 5.2 4.9 3.4 3.9 4.8 4.0 2.9 6.3 4.6 3.0 4.4
##  [505] 4.6 5.9 3.5 4.2 5.8 3.9 6.0 4.1 4.7 4.4 3.8 5.5 5.3 3.6 4.1 4.1 4.5 3.5
##  [523] 2.9 4.2 4.7 4.6 4.5 3.6 4.1 4.7 4.3 5.7 5.4 4.4 4.6 3.6 5.3 5.2 5.5 3.2
##  [541] 4.2 4.4 5.7 4.2 5.8 4.7 4.0 3.9 4.1 2.0 5.2 4.5 4.3 4.9 4.9 5.6 4.9 4.6
##  [559] 4.5 3.4 3.4 5.4 5.4 5.1 4.3 4.6 2.9 5.8 4.0 5.3 5.2 3.7 3.3 1.7 4.9 3.9
##  [577] 5.2 4.5 3.7 5.3 5.1 3.6 4.9 4.7 3.2 5.0 4.4 3.6 5.4 5.3 4.5 4.6 2.7 4.1
##  [595] 6.8 3.9 3.6 4.4 3.7 4.2 2.9 5.0 3.4 4.2 5.5 3.1 3.9 5.7 4.9 5.5 5.5 3.4
##  [613] 3.8 4.7 5.1 4.1 5.6 4.3 5.3 4.9 4.8 4.6 3.9 4.0 3.6 4.1 5.6 2.9 4.6 4.1
##  [631] 4.8 3.4 3.5 5.0 3.3 4.6 3.9 5.0 2.8 3.2 4.8 4.1 4.7 5.8 5.8 5.7 5.8 4.4
##  [649] 4.4 5.0 4.4 4.8 4.2 4.1 5.1 2.8 4.3 4.5 5.7 4.3 4.1 4.5 4.0 4.3 3.8 4.5
##  [667] 4.4 5.8 5.6 4.4 2.8 4.4 3.0 4.3 4.3 4.5 4.1 4.7 3.6 4.6 7.1 5.3 5.0 4.7
##  [685] 4.9 4.2 4.7 5.7 4.8 4.2 4.6 4.4 4.0 2.9 4.8 4.0 4.2 4.9 5.0 5.2 4.0 3.5
##  [703] 3.9 3.1 3.9 4.5 3.4 3.7 5.3 4.6 4.9 5.8 4.1 3.7 4.2 4.9 5.4 5.2 4.3 5.2
##  [721] 4.6 2.1 5.0 3.6 5.5 4.5 3.1 7.3 5.4 5.3 5.1 4.2 4.3 3.6 4.9 3.5 4.3 4.2
##  [739] 4.6 6.4 6.2 3.9 4.9 5.0 5.1 3.9 4.7 4.2 5.5 5.8 3.3 3.3 5.5 4.2 5.9 5.2
##  [757] 4.2 5.4 4.3 4.5 6.3 4.8 3.5 4.4 5.0 5.0 4.9 5.0 3.7 5.3 3.9 3.8 5.9 3.3
##  [775] 4.5 4.1 3.7 3.8 4.3 4.6 5.4 4.4 4.3 3.8 5.4 4.9 5.8 4.9 3.5 6.1 3.7 3.9
##  [793] 4.7 5.1 4.5 3.8 2.8 3.9 4.3 5.2 4.9 4.5 4.2 3.8 4.2 3.5 4.2 3.7 4.8 2.9
##  [811] 3.2 3.1 4.0 5.6 5.8 4.4 4.5 5.0 3.0 4.6 3.2 4.0 3.9 4.3 5.2 5.0 4.6 4.5
##  [829] 7.1 5.9 5.2 3.1 5.1 4.7 4.7 3.9 5.1 4.0 4.0 3.7 4.6 4.4 4.7 5.4 3.5 5.7
##  [847] 4.4 4.1 4.4 3.5 4.2 4.1 3.9 3.0 5.1 5.1 3.3 6.2 4.1 4.5 4.2 4.1 4.2 3.9
##  [865] 5.3 5.6 4.1 5.9 4.9 5.3 3.5 5.7 4.7 4.1 5.1 5.4 3.3 1.6 5.2 5.3 4.0 4.9
##  [883] 2.6 5.7 2.9 4.5 3.6 6.1 4.2 4.0 4.8 3.0 3.9 4.4 2.8 4.5 5.0 3.5 5.0 3.9
##  [901] 3.9 3.7 5.6 4.8 5.4 3.5 5.7 5.7 6.6 4.9 4.7 3.1 5.0 4.7 3.9 3.4 4.9 3.8
##  [919] 3.0 4.4 5.1 5.5 5.1 4.8 3.4 3.6 5.5 4.8 4.0 4.5 5.7 4.7 4.8 3.8 4.8 4.5
##  [937] 3.9 4.2 4.7 4.7 2.9 4.0 3.6 4.0 2.5 3.2 5.7 4.5 3.2 5.1 3.3 4.0 4.8 6.4
##  [955] 5.8 5.2 5.5 5.5 2.9 5.2 3.2 2.8 4.6 4.7 5.5 3.0 4.9 7.1 5.0 2.9 4.9 4.6
##  [973] 3.3 4.3 5.1 4.0 6.0 3.9 5.8 4.8 3.8 4.9 4.3 2.5 4.8 4.7 3.6 4.9 4.3 3.9
##  [991] 3.2 3.8 5.9 3.9 6.0 3.6 3.7 5.0 5.4 5.3
## 
## $func.thetastar
## [1] 0.0146
## 
## $jack.boot.val
##  [1]  0.515988372  0.344444444  0.264130435  0.229737609  0.053581662
##  [6]  0.000591716 -0.140109890 -0.270434783 -0.359718310 -0.498005698
## 
## $jack.boot.se
## [1] 0.9384336
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
##    [1] 4.2 5.3 4.5 4.8 2.4 4.9 3.4 3.8 5.5 4.3 2.2 6.1 4.7 4.0 2.4 5.2 3.7 5.3
##   [19] 4.4 3.3 3.6 4.4 3.4 3.3 4.3 3.4 3.3 4.3 3.7 5.6 4.6 4.0 5.3 3.9 5.1 3.3
##   [37] 4.7 2.2 3.8 2.8 4.4 4.5 4.6 5.0 5.3 3.8 6.6 5.8 5.1 4.4 4.8 3.2 5.7 2.9
##   [55] 6.5 4.4 4.2 6.7 4.6 4.3 4.9 3.0 5.7 3.3 4.0 4.7 4.9 4.2 4.9 4.0 5.6 5.4
##   [73] 4.6 4.8 5.5 3.7 4.2 3.2 4.9 4.6 3.5 3.4 3.8 3.1 6.1 6.0 4.5 3.7 4.3 5.9
##   [91] 5.7 4.3 4.4 3.6 5.3 3.2 4.4 3.8 4.0 3.5 4.2 5.3 5.0 4.2 3.6 5.3 4.9 5.1
##  [109] 4.1 4.3 4.0 5.5 5.5 4.1 3.5 4.2 5.7 5.0 3.5 3.4 4.5 4.8 4.3 4.9 3.4 3.2
##  [127] 4.9 3.8 5.9 6.0 2.1 5.0 5.2 4.0 5.0 4.0 3.3 3.6 4.2 4.1 5.3 2.7 3.8 4.7
##  [145] 4.2 3.5 4.5 4.4 4.9 6.8 4.6 5.6 4.7 3.7 6.0 5.1 6.6 4.2 3.6 4.8 6.5 4.9
##  [163] 5.1 4.9 6.1 5.2 4.9 4.3 5.2 4.8 6.5 4.8 5.7 4.8 2.8 5.5 4.6 4.6 2.8 4.1
##  [181] 4.6 5.4 4.3 4.6 5.0 5.3 4.1 4.4 5.4 3.4 3.1 6.1 3.8 5.0 5.7 3.9 5.6 4.5
##  [199] 4.9 5.4 5.0 1.7 4.8 4.0 5.4 4.6 4.3 3.6 5.9 5.3 3.2 4.9 4.1 4.5 3.8 3.5
##  [217] 3.2 3.9 4.8 4.3 4.7 3.7 2.2 3.6 4.9 3.8 4.6 3.7 3.1 4.6 3.2 5.8 3.6 5.2
##  [235] 4.4 5.6 4.8 5.6 4.0 4.6 3.9 4.6 3.2 5.3 4.0 6.1 4.2 4.3 4.7 4.1 4.6 5.8
##  [253] 6.2 5.2 4.1 5.3 6.1 3.8 4.8 3.1 5.8 4.8 6.0 4.7 5.5 5.4 4.7 5.8 3.7 4.1
##  [271] 5.4 3.4 4.5 4.1 4.5 5.2 2.7 5.3 4.8 5.0 4.7 5.0 5.7 5.5 5.3 4.4 4.5 4.9
##  [289] 5.0 5.1 4.6 3.7 2.7 6.8 5.6 3.4 5.0 4.6 4.5 5.4 4.6 5.1 5.5 3.9 4.0 4.3
##  [307] 5.7 5.2 4.5 3.8 4.1 3.9 4.2 3.9 3.6 4.5 5.1 3.6 4.5 5.4 4.5 3.9 4.3 5.1
##  [325] 4.1 6.4 4.1 4.1 5.5 3.9 4.8 3.7 5.4 4.7 4.4 4.9 3.8 5.3 2.4 3.5 3.5 6.0
##  [343] 5.6 3.3 5.8 4.9 3.0 7.1 4.2 4.1 4.2 5.0 5.2 5.1 4.5 4.5 4.9 4.9 4.1 4.4
##  [361] 5.1 5.0 4.4 4.1 3.0 5.0 4.4 3.5 3.7 2.6 5.0 4.6 4.2 5.2 3.8 4.0 3.1 5.4
##  [379] 4.7 4.7 4.2 5.1 4.0 3.9 4.9 5.3 4.8 5.9 4.7 4.5 3.6 3.3 4.8 4.3 6.7 6.2
##  [397] 5.2 5.6 4.4 4.3 4.5 5.8 3.6 4.9 4.8 4.5 5.0 4.0 3.3 3.4 4.0 3.8 3.2 3.9
##  [415] 5.1 5.7 5.5 3.5 3.8 3.7 5.0 4.1 4.5 4.4 3.6 6.1 4.4 4.6 3.8 4.1 4.0 3.2
##  [433] 5.1 3.0 4.2 4.4 4.2 4.1 4.2 4.4 5.1 3.6 3.6 3.4 4.4 4.1 5.4 4.2 5.2 5.8
##  [451] 4.4 4.9 5.1 5.0 5.6 3.8 4.1 3.9 5.5 3.8 6.1 5.5 4.4 6.5 4.6 3.8 1.9 3.5
##  [469] 4.7 3.9 4.7 5.6 4.1 3.6 4.1 5.7 4.1 4.4 5.5 5.5 5.0 3.9 5.2 4.5 4.3 3.0
##  [487] 6.8 5.1 3.6 4.2 4.8 3.1 5.8 4.9 6.4 4.9 6.1 4.9 6.0 5.5 5.1 5.7 3.3 5.3
##  [505] 6.7 4.6 3.7 2.7 6.1 4.0 3.5 4.5 4.1 5.1 4.5 3.2 4.1 5.4 5.7 4.3 4.2 3.8
##  [523] 5.0 4.4 3.6 4.1 3.4 4.5 5.3 3.8 5.5 4.1 3.6 4.4 4.6 6.0 4.4 5.2 3.2 3.4
##  [541] 3.6 3.9 4.8 4.6 3.1 3.1 3.2 5.1 3.7 4.1 6.1 4.3 4.9 4.6 5.0 4.3 5.6 4.5
##  [559] 5.5 4.0 3.6 2.7 5.9 4.3 3.7 4.6 6.0 5.4 5.4 4.2 3.3 3.5 4.3 5.1 4.3 3.0
##  [577] 3.9 5.2 4.7 3.6 4.3 5.6 4.0 5.3 3.9 4.3 4.8 6.6 4.1 5.0 3.8 3.8 4.4 5.7
##  [595] 4.1 4.6 5.4 5.3 3.9 3.3 5.4 2.3 3.2 4.0 4.4 6.4 4.4 5.3 4.1 4.9 5.9 5.2
##  [613] 4.6 4.0 4.8 2.7 6.0 5.0 4.8 4.1 4.4 3.9 3.6 4.8 3.3 4.2 2.5 2.9 5.8 4.8
##  [631] 5.9 2.5 3.8 4.3 3.2 4.9 4.6 3.9 6.0 3.9 5.4 4.6 3.8 3.4 4.7 3.4 4.9 6.0
##  [649] 4.8 5.4 3.5 3.6 5.6 4.8 4.5 4.6 3.3 3.8 5.7 4.3 4.5 3.8 3.4 4.9 4.3 4.6
##  [667] 4.9 6.0 3.9 4.1 4.1 3.8 3.6 3.7 3.8 4.5 4.8 5.4 5.0 4.4 5.6 3.2 5.7 4.7
##  [685] 5.0 3.1 3.9 3.6 4.4 4.6 4.6 4.6 3.9 4.1 3.8 4.5 3.3 3.4 5.7 4.7 4.1 5.0
##  [703] 4.9 3.6 5.7 6.1 4.4 7.1 5.1 5.0 3.4 5.0 4.4 3.4 4.5 3.4 3.0 5.2 3.8 5.8
##  [721] 4.3 4.3 3.6 4.6 3.1 5.9 5.1 3.0 4.0 6.2 4.6 4.9 4.0 5.1 6.6 4.5 3.6 3.3
##  [739] 3.9 4.6 4.5 6.5 3.8 4.5 4.4 5.4 5.1 4.8 4.4 3.8 4.8 5.5 4.1 6.3 4.4 4.1
##  [757] 4.6 5.2 4.5 5.3 5.3 3.7 3.1 5.4 3.2 4.3 5.3 6.2 3.2 5.6 5.0 5.6 2.2 4.0
##  [775] 4.9 5.3 5.6 4.1 3.2 3.5 5.3 3.0 4.3 4.9 4.2 4.3 5.5 3.9 4.3 3.7 2.9 4.9
##  [793] 4.4 2.6 5.1 3.7 6.0 6.6 3.8 3.8 6.8 4.8 4.1 4.7 1.8 5.7 7.0 6.2 3.8 4.1
##  [811] 5.0 5.3 4.6 2.7 4.8 5.7 5.6 5.6 5.3 4.1 5.6 4.9 4.5 5.1 4.4 3.4 3.9 3.9
##  [829] 3.5 4.9 5.1 4.8 5.1 4.0 3.7 4.8 4.4 4.4 3.9 4.1 5.5 4.7 4.3 4.1 4.9 4.0
##  [847] 3.7 4.0 7.2 4.5 4.1 3.9 4.4 4.2 3.9 4.7 5.5 4.1 4.2 4.0 3.9 3.8 3.3 4.9
##  [865] 5.0 3.3 4.2 5.0 3.5 2.8 3.1 5.5 3.5 4.3 2.6 5.2 4.2 4.9 3.4 4.5 3.6 3.6
##  [883] 5.2 4.6 4.5 3.1 4.0 5.6 5.7 4.5 3.9 4.5 4.9 4.7 3.9 6.1 4.8 3.5 4.0 4.2
##  [901] 5.2 2.3 5.2 4.7 5.6 5.9 3.6 4.3 4.0 5.9 4.9 5.1 3.3 4.4 4.4 5.5 4.4 4.2
##  [919] 5.2 4.8 4.4 4.6 4.8 4.8 3.3 5.4 6.3 4.1 4.8 6.3 5.7 4.1 3.2 4.5 4.5 4.1
##  [937] 4.5 5.9 4.6 3.8 5.8 4.5 3.9 4.6 3.8 6.2 4.1 4.0 4.9 4.6 4.9 4.2 4.9 4.0
##  [955] 3.8 5.0 5.1 4.5 3.8 3.2 4.3 5.2 3.0 4.5 5.7 5.4 4.3 5.6 3.8 3.8 5.0 2.8
##  [973] 3.9 5.2 5.2 3.6 2.9 4.6 3.9 3.9 4.8 3.7 5.8 4.1 4.2 2.6 2.5 3.9 6.2 5.5
##  [991] 5.1 4.8 4.3 4.4 4.0 4.4 4.1 3.2 3.8 4.1
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.492 5.400 5.300 5.300 5.100 5.000 4.900 4.800 4.600 4.400
## 
## $jack.boot.se
## [1] 1.018913
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
## [1] 0.7406021
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
##   2.865052   6.199635 
##  (1.214112) (2.871181)
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
## [1]  1.88094573  1.21155971  0.82876218 -0.12284953 -0.46897578  0.04447859
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
##    [1]  5.219082e-01  3.498726e-02  3.770124e-01  1.087683e+00  1.081588e-01
##    [6]  4.285963e-01  1.204074e+00 -5.352483e-02 -8.727567e-02  1.910828e+00
##   [11]  5.861146e-01  7.070715e-01  7.771419e-01  9.596764e-01  5.019032e-01
##   [16] -4.545154e-01  1.238973e-01 -4.822768e-01  2.870321e-01  1.000164e+00
##   [21]  1.155264e+00  8.839838e-01  1.192488e+00  8.826359e-01  8.242738e-01
##   [26]  8.110235e-01  7.169033e-01  6.537136e-01  1.218102e-01  2.503083e-01
##   [31]  7.304845e-01  5.088551e-01  3.985572e-01  8.807726e-01  7.120198e-01
##   [36]  1.166847e+00  1.099304e-01  1.970399e+00  6.747472e-01  1.383255e+00
##   [41]  1.895162e+00  4.113013e-01 -3.908842e-01  2.079646e+00  5.909712e-02
##   [46] -3.953200e-01  2.001388e+00  3.035236e-01  3.057449e-01  1.322705e+00
##   [51]  8.114809e-01  1.151483e+00  3.933225e-01  6.520006e-01  1.116362e+00
##   [56]  4.469913e-01  3.004184e-01  1.329923e+00  1.253133e+00  3.820642e-01
##   [61]  9.989295e-01  1.854299e+00  7.546925e-01  8.270159e-01  2.630000e-01
##   [66]  7.466804e-01  1.819577e+00  1.318244e+00  7.931390e-01  3.550346e-01
##   [71]  4.787068e-01  1.621439e+00  7.157202e-01  6.526625e-01  1.747690e+00
##   [76]  7.818902e-01  1.182215e+00  8.420632e-01  7.125115e-01  1.323346e+00
##   [81]  6.634952e-01  5.276654e-01  1.406941e+00  6.786231e-01  2.617208e-01
##   [86] -3.409282e-01  3.389177e-01 -3.645625e-03  4.137244e-01  5.848711e-01
##   [91] -6.838878e-01  1.220791e+00  8.466911e-02  4.866264e-01  8.063747e-01
##   [96]  2.216361e+00  2.153341e-02 -3.084506e-02  8.236557e-01  3.772812e-01
##  [101]  4.691963e-02  1.157011e+00  7.928838e-01  1.100700e+00  1.417643e+00
##  [106]  7.786178e-01  7.211279e-01  1.045622e+00  6.677557e-01  3.421123e-01
##  [111] -4.511807e-02  1.434872e-01  1.173602e+00  1.745640e+00  3.089161e-01
##  [116]  1.017905e+00  1.225187e+00  1.829845e+00  6.584816e-01  2.339528e+00
##  [121]  1.193151e+00  2.919502e-01  8.515676e-01  1.095067e+00  4.902836e-01
##  [126]  3.128158e-01  3.163594e-01  2.947956e-01  6.933541e-01  1.234418e+00
##  [131]  1.108262e-01  7.719747e-01  6.937635e-01  7.313975e-01  1.589758e+00
##  [136]  1.121386e+00  4.609533e-01  7.049554e-01  1.249409e+00 -2.144893e-01
##  [141]  9.936815e-01  1.928961e+00  1.383057e+00  6.666161e-01  1.416563e+00
##  [146]  1.877799e+00  3.791480e-01  1.754482e+00  7.803219e-01  7.093635e-01
##  [151] -4.706134e-01  6.367420e-01  1.360615e+00  8.255175e-01  7.486358e-01
##  [156]  3.692285e-01  8.316951e-01  1.117208e+00  6.720215e-01  6.592394e-01
##  [161] -2.734794e-01  9.830508e-01  5.253003e-01  1.732870e+00  3.874068e-01
##  [166]  9.866188e-01  5.291226e-01  9.410190e-01  1.224383e+00  6.643705e-01
##  [171]  4.607106e-01  1.092309e+00  2.522841e-01 -6.581926e-01  7.434419e-01
##  [176]  3.930426e-01  1.302700e+00  4.152380e-01  4.458932e-01  9.982943e-01
##  [181]  1.480018e+00  1.747221e+00  9.373927e-01  4.271608e-02  3.827118e-01
##  [186]  6.136098e-01  2.553006e-01  6.602595e-01 -4.286955e-05  2.050467e+00
##  [191]  7.522367e-01  7.155022e-01  3.086897e-01  3.825019e-01  1.084019e+00
##  [196] -8.935475e-02  2.231649e+00  6.656479e-01  1.465317e-01  5.241458e-01
##  [201]  1.564944e+00  7.398818e-01  2.840597e-01  1.140824e+00  1.271031e+00
##  [206]  1.118667e+00  3.729657e-01  1.337775e+00  4.767045e-01  3.458530e-01
##  [211]  7.742666e-01  5.387893e-01  1.051024e+00  9.536118e-01 -4.074652e-01
##  [216]  4.097567e-01  1.251655e+00  1.957549e+00 -1.190403e-01  1.160863e+00
##  [221]  3.699698e-01  1.110982e+00  2.073483e-01 -5.141353e-01 -3.720725e-02
##  [226]  6.770947e-01  6.270410e-01  6.664175e-01 -2.479248e-01  5.131637e-01
##  [231]  7.009398e-01  1.461265e+00  1.107903e+00  7.910425e-01  3.707031e-01
##  [236]  1.223348e+00  1.540730e+00  1.499401e+00  1.231543e+00  1.841293e+00
##  [241]  7.071464e-01  1.984175e+00  1.295298e+00  3.165360e-01  1.062659e+00
##  [246]  7.536474e-01  1.058312e+00 -4.417174e-01  7.742666e-01  6.364408e-01
##  [251] -3.997194e-01  7.495583e-01  6.523460e-01  8.181029e-01  5.706239e-01
##  [256]  1.198186e+00  9.277351e-02  8.885091e-01  7.878965e-01  7.823936e-01
##  [261]  1.569748e-01  1.963045e-01 -3.726553e-02 -6.838878e-01  5.042231e-01
##  [266]  8.244277e-01  4.411267e-01  9.757408e-01  2.056011e+00 -1.189096e-01
##  [271]  6.467777e-02  3.411952e-01  5.639481e-01  3.951618e-01  1.105149e+00
##  [276] -3.969889e-01  4.293836e-03  1.136592e+00  9.602799e-01  1.264615e+00
##  [281]  5.661393e-02  1.456922e+00  8.001857e-01  8.539107e-02  7.506446e-01
##  [286]  1.243124e+00  7.228447e-01  1.851409e+00  1.546926e+00  4.762746e-01
##  [291]  1.036127e+00  1.265414e+00  7.313115e-01  1.028030e+00 -8.749794e-03
##  [296]  8.375057e-01  1.671988e+00  2.080348e+00  7.769973e-01  5.033550e-01
##  [301] -3.126356e-02  1.262186e+00  1.040059e+00  2.016918e+00  4.297421e-01
##  [306] -4.026198e-01  8.993350e-03  8.175670e-01  2.129944e+00  2.108822e-01
##  [311]  1.252657e+00  7.774279e-01  8.533362e-01  5.975728e-01  1.445229e+00
##  [316]  2.349294e-01  7.073657e-01  1.872275e+00  7.304845e-01  3.886889e-01
##  [321]  3.909181e-01  2.669204e-01  2.827169e-01 -4.083774e-01  3.045636e-01
##  [326]  7.627228e-01  3.727393e-01  5.216364e-02  6.999938e-01  7.214275e-01
##  [331]  2.905834e-02  9.249547e-01  2.144237e+00  3.337610e-01  1.353515e+00
##  [336]  6.608050e-01  2.432759e-01  7.738006e-01  4.369662e-01  2.686288e-01
##  [341]  9.558418e-01 -2.073912e-02  1.016537e+00  1.416246e+00  6.680354e-01
##  [346]  5.673926e-01  1.038098e+00  4.456933e-01  7.101188e-01  4.166966e-02
##  [351]  1.277850e+00  4.107287e-01  1.182711e-01  1.200474e+00  7.943408e-01
##  [356]  8.157128e-01  1.182455e+00  1.256547e+00  5.688601e-01  1.898580e-04
##  [361]  3.794408e-01  4.180692e-01  2.349309e+00  3.727393e-01  4.721627e-01
##  [366]  4.416767e-01  7.935608e-01  4.239828e-01  7.800583e-01  6.064207e-01
##  [371]  2.648828e-01  3.365894e-01  1.803115e+00  1.370620e+00 -8.830387e-02
##  [376]  8.462301e-01  6.315902e-01 -6.002776e-02  7.991083e-01  1.637161e+00
##  [381]  1.156557e+00  1.045853e+00  7.884076e-01  1.262928e+00  7.895938e-01
##  [386]  1.287768e+00  7.751298e-01  3.373091e-01  1.009015e+00  4.478547e-01
##  [391]  1.238035e+00  1.578925e+00  1.498067e+00  4.535571e-02  6.537251e-01
##  [396]  3.537451e-01  1.715848e+00  3.904665e-01  5.019685e-01  1.194336e+00
##  [401]  1.818143e+00  7.108372e-01 -1.826516e-02 -1.006290e-01 -2.383950e-01
##  [406]  7.987097e-01 -3.317771e-01  7.888755e-01  9.734421e-01  1.258200e+00
##  [411]  1.798875e+00  3.048382e-01 -2.882807e-01  4.741132e-01  1.477451e+00
##  [416]  3.249852e-01  8.725190e-01  2.160430e-02  6.803596e-03  4.184111e-01
##  [421]  3.443342e-01  1.124304e-01  1.172016e+00  4.096790e-01  8.054816e-01
##  [426] -2.676340e-01 -1.457319e-02  7.798189e-01  1.670619e+00  8.313614e-01
##  [431]  1.154687e+00  1.170245e+00  1.230842e+00  7.486424e-01 -1.632370e-01
##  [436]  8.270091e-01  1.703050e+00  1.219309e+00  7.613539e-01  9.951941e-01
##  [441]  7.316819e-01  7.849692e-01  3.736794e-01  1.259928e+00  1.414189e+00
##  [446] -3.346725e-03  1.997880e+00  1.104520e+00  8.096038e-01 -1.099542e-01
##  [451]  8.410059e-01 -1.260332e-01  4.558147e-01  4.221363e-01  1.220969e+00
##  [456]  3.949056e-01  8.645891e-01  7.424618e-01  3.603806e-01  7.336468e-01
##  [461]  7.110985e-01  3.762574e-01  5.770799e-03  5.006716e-01  1.400327e-02
##  [466]  1.816422e+00 -3.794412e-02 -8.890030e-01  3.059504e-01  3.740945e-01
##  [471]  2.027819e+00  1.224118e+00 -1.043645e-01  1.223267e+00  7.155850e-01
##  [476]  8.050038e-01  6.662536e-01  9.803521e-04  6.206075e-01  2.835860e-01
##  [481]  5.873008e-01  1.209358e+00  3.766985e-01  7.575034e-01  7.787805e-01
##  [486]  1.206705e+00 -7.787864e-02  7.161799e-01  1.227008e+00  3.527724e-01
##  [491]  2.591561e-01  6.959517e-01  1.967646e+00  2.277786e+00  3.697004e-01
##  [496]  6.827838e-01  8.007719e-01  4.840799e-02  3.702505e-01  1.307032e+00
##  [501]  1.759334e-02  1.058589e+00  7.741956e-01  9.201112e-02  6.934387e-02
##  [506]  1.260290e+00  1.327145e+00  7.792781e-01  4.578609e-01  8.243086e-02
##  [511]  8.408240e-01  1.141820e+00  9.041072e-01  1.111756e+00  7.142491e-01
##  [516]  1.736526e+00  1.822508e+00  9.162873e-01  3.638805e-01  4.231917e-01
##  [521]  4.408706e-01  3.526965e-01  8.310684e-01  6.373085e-01  6.488235e-01
##  [526]  1.379005e+00  3.521294e-01  1.376497e+00  1.254824e+00  1.338948e+00
##  [531]  9.000609e-01 -3.707379e-01  1.802020e+00 -3.946548e-01  1.331284e+00
##  [536] -4.589773e-01  2.277505e-01  6.630810e-01  1.203780e+00  7.082700e-01
##  [541]  8.403644e-01  1.186842e+00  1.401152e+00  7.443817e-01  6.281438e-01
##  [546]  9.571344e-01  1.230874e+00  5.390188e-03  5.861146e-01  6.842040e-02
##  [551]  3.952117e-01  1.158262e+00  4.854400e-01 -1.794990e-01  3.338293e-01
##  [556]  9.454875e-01  6.942750e-01  7.749362e-01  1.164743e+00  3.780879e-02
##  [561]  2.175613e+00 -6.126292e-02  5.379968e-01  1.290259e+00  3.735707e-01
##  [566]  1.250448e+00  2.423476e-01  5.050124e-01  7.572264e-01  1.131487e+00
##  [571]  7.361950e-01  1.276652e+00  7.726783e-01  6.648449e-01 -8.498034e-01
##  [576]  7.518204e-01  7.893253e-01  4.676346e-01  1.503537e+00  1.227612e+00
##  [581]  1.265121e+00  8.138931e-02  5.199525e-01  8.610499e-01  9.749591e-01
##  [586] -1.599816e-01  6.101996e-01  8.372999e-01  1.115090e+00  1.289567e+00
##  [591]  1.811878e+00  5.281550e-01 -5.518651e-02  1.980578e+00  4.897541e-01
##  [596] -2.478728e-02  1.103941e+00  1.910515e+00  1.188478e+00  7.894254e-01
##  [601]  8.989950e-02  1.305981e+00  3.961116e-01  6.494659e-01  9.174859e-01
##  [606]  1.016795e+00 -2.736291e-01  7.771419e-01  8.459112e-01  1.083536e+00
##  [611]  8.656425e-01  6.618010e-01  1.342077e+00  7.843903e-01  7.171140e-01
##  [616]  8.274892e-01 -1.098871e+00  1.556663e+00  2.146434e+00  7.652674e-01
##  [621]  2.722943e-01  1.259524e+00  8.647681e-01  3.181496e-01  6.863655e-01
##  [626]  1.372648e+00  8.070765e-01  7.954059e-01  8.838923e-01  2.142174e-01
##  [631] -2.366824e-01  7.786132e-01  3.852931e-01  1.165232e+00  7.728768e-01
##  [636]  3.967594e-01 -1.988962e-02 -5.141578e-01  8.890773e-01  1.797452e+00
##  [641]  4.263776e-01  1.233401e+00  1.384441e+00  1.140997e+00  1.245917e+00
##  [646]  6.622660e-01  1.232610e+00  7.030833e-01  2.167012e+00  6.958429e-01
##  [651] -1.767140e-01  1.279653e+00  5.271278e-02  4.162132e-01  2.874577e-01
##  [656]  1.252674e+00 -2.216975e-02  1.286623e+00  2.567670e-01  1.308829e+00
##  [661] -7.466142e-02  3.858401e-01  6.662536e-01  3.731062e-01  7.293393e-01
##  [666]  8.875615e-01 -1.186200e-02  1.422074e+00  5.224139e-01  1.025605e+00
##  [671]  1.078169e+00  6.988911e-01  7.575353e-01  1.221574e+00  1.262807e+00
##  [676]  1.988908e+00  3.453752e-01  8.247568e-01  1.142740e+00  4.200204e-01
##  [681]  8.533362e-01  1.964361e+00  1.699155e+00  2.303760e+00  6.317892e-01
##  [686]  6.485058e-01  1.442976e+00  3.811597e-01  3.704782e-01  1.065156e+00
##  [691] -9.265432e-02  8.581737e-01  4.536997e-01  4.776607e-01  1.036413e+00
##  [696]  3.739176e-01 -3.003410e-01  2.809225e-01  6.989029e-01 -2.057975e-01
##  [701]  1.159473e+00  7.480935e-01  6.512874e-02  1.187373e+00  6.750172e-01
##  [706]  8.934245e-01  7.845661e-01  9.006992e-01  8.448617e-01 -3.861269e-01
##  [711]  2.458032e-01  1.260163e+00  4.305853e-03  1.097494e+00  1.745330e+00
##  [716]  4.324522e-01  3.719177e-01  1.146838e-01  1.283935e+00  8.874822e-01
##  [721]  8.061516e-01  1.826875e-01  6.114883e-01  9.129347e-01  7.166606e-01
##  [726]  1.252595e+00  7.906519e-02 -9.191026e-01  7.409524e-01  4.035973e-01
##  [731]  1.449176e+00  1.334338e+00  7.376280e-01  5.213536e-01  8.379765e-01
##  [736]  9.010043e-01  1.265779e+00  1.967194e+00  7.560049e-01  9.835502e-01
##  [741]  1.229029e+00  8.432661e-01  3.207983e-01  1.033507e+00  2.248786e+00
##  [746]  1.320894e+00  1.040700e+00  7.406021e-01  1.046032e+00  6.099293e-01
##  [751]  8.471926e-01  1.118965e+00  4.759331e-02  1.458778e+00  9.952114e-01
##  [756]  1.018996e+00  1.308562e+00  1.154902e+00  4.496872e-01  7.534547e-01
##  [761]  2.013661e+00  1.165185e+00  7.675421e-01 -1.449928e-02 -3.424080e-01
##  [766]  7.106064e-01  1.872275e+00  1.974697e+00  8.243611e-01  1.223829e+00
##  [771]  2.453511e-01  6.954069e-01  1.607812e+00  2.031944e+00  6.153008e-01
##  [776]  9.014385e-01  1.338026e+00  8.279105e-01  8.651026e-01  1.163608e+00
##  [781]  1.125739e+00  6.604747e-01  8.421473e-01  1.587246e+00  7.359289e-02
##  [786]  1.369296e-01  6.876045e-01  3.175154e-01  7.435734e-01  8.905750e-01
##  [791] -1.185279e-01  6.536000e-01  6.799590e-01 -1.200035e-01  1.064639e-01
##  [796]  1.144182e+00  3.848960e-01  7.769973e-01  4.022093e-01  7.152589e-01
##  [801]  6.464387e-01  3.693813e-01  1.257703e+00  1.304042e+00  4.005575e-01
##  [806]  1.257479e+00  7.520418e-01  1.800646e+00  8.853148e-01  4.189321e-01
##  [811]  1.154919e+00  1.290668e+00  1.671494e-02  6.338685e-01 -6.821338e-02
##  [816]  1.211964e+00  1.149995e+00  4.305361e-01  3.946858e-01  3.662729e-01
##  [821]  9.364116e-01 -4.009245e-01  1.463810e+00  1.074739e-01  1.235011e+00
##  [826]  1.037293e+00 -5.553693e-01  7.949020e-01  7.263831e-01  7.829503e-01
##  [831]  3.047630e-01  7.105043e-01 -3.565317e-01  6.963332e-01  1.956657e-01
##  [836]  7.252130e-01  8.073220e-01  4.836394e-01  7.953880e-01  1.330474e+00
##  [841]  7.890998e-01  1.673392e+00  3.626937e-01  7.074271e-01  1.715906e+00
##  [846]  6.816871e-01  7.286663e-01 -8.699058e-02  4.277585e-01  4.010902e-01
##  [851]  1.418407e+00  6.542176e-01  1.678535e+00  1.399713e+00  3.847876e-01
##  [856]  1.025549e+00  7.455198e-01  1.547474e+00  6.372357e-01  4.114879e-01
##  [861]  1.275338e+00 -2.033865e-01 -6.328067e-01  4.015214e-01  1.800646e+00
##  [866]  3.586242e-01  3.762574e-01  2.627145e-01  8.073673e-01  1.689731e+00
##  [871]  1.297975e-01  5.957103e-01 -1.862491e-02  2.134216e+00  2.364151e+00
##  [876]  1.954989e+00  4.505234e-01  1.015688e+00  1.445919e+00  7.536981e-01
##  [881]  9.032838e-01  1.124591e+00  7.145254e-02  1.404804e+00 -3.079049e-01
##  [886]  7.578363e-01  7.131296e-01  2.455882e-02  9.734019e-01  7.023563e-01
##  [891]  1.200037e+00 -9.021567e-03  7.898766e-01  1.863191e-02  1.474233e-02
##  [896]  1.938404e-01  7.139318e-02  1.264157e+00  2.324653e-01  7.653987e-01
##  [901]  1.189126e+00 -2.708191e-01  8.068470e-01  5.110752e-01  1.259318e+00
##  [906]  2.880065e-01  4.021625e-01  7.643624e-01  3.248277e-01  8.386254e-02
##  [911]  2.527636e-01  7.318374e-01  7.272003e-01  7.440975e-01  3.662193e-01
##  [916]  6.952411e-02  1.289117e+00 -7.427131e-02  3.501682e-01  9.338025e-01
##  [921]  4.481956e-01  1.059930e+00  5.142229e-02  6.863655e-01  6.736806e-01
##  [926]  6.756197e-01  3.270383e-01  9.445729e-05  1.210796e+00  1.830451e+00
##  [931]  1.562897e+00  6.963858e-01  8.483242e-01 -2.149553e-01  1.237877e+00
##  [936]  4.831385e-01 -6.654274e-02  4.028825e-01  6.460799e-01  1.565524e+00
##  [941]  1.374935e+00  4.009318e-01  2.055140e+00  1.355217e+00  1.069231e+00
##  [946]  1.161057e+00  1.157397e+00  6.903774e-01 -3.696015e-01  1.235884e+00
##  [951]  1.557728e-01  3.225363e-01  7.506899e-01 -6.074654e-02  6.047211e-01
##  [956]  7.862777e-01  6.979008e-01  6.804212e-01  8.054816e-01  7.158044e-01
##  [961]  4.144304e-01  4.557619e-01  1.100372e+00 -3.692974e-01  5.609960e-01
##  [966]  1.090489e+00  1.871822e+00  6.360225e-01  1.280600e+00  1.882084e+00
##  [971]  2.222066e+00  6.436324e-01 -8.599112e-02  7.576719e-01  7.076053e-01
##  [976] -3.982521e-01  6.633751e-01  4.781478e-01  7.096312e-01  8.116909e-01
##  [981]  7.522367e-01  6.646488e-01  2.143222e+00  2.499883e-01  7.101586e-01
##  [986]  1.223829e+00  4.406675e-01  1.193312e+00  1.553834e+00  1.746715e+00
##  [991]  7.764223e-01  7.943430e-01  1.102003e+00  7.774279e-01  4.180818e-01
##  [996]  1.432724e-01  7.536981e-01  7.548815e-02  1.108977e+00  7.609720e-01
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

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced
```

```r
fit2
```

```
##       mean          sd    
##   0.46213041   0.28241988 
##  (0.08930901) (0.06314820)
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
## [1]  0.7163145 -0.4426811  0.1125541 -0.7551891 -0.9823004  0.2696337
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
## [1] -0.017
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8934928
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
## t1*      4.5 0.006806807   0.8998295
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 2 3 4 7 8 
## 3 1 2 1 1 2
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
## [1] -0.0352
```

```r
se.boot
```

```
## [1] 0.9129592
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

