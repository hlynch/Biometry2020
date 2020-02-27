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
## 0 3 4 5 8 
## 3 2 1 3 1
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
## [1] -0.0259
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
## [1] 2.717255
```

```r
UL.boot
```

```
## [1] 6.230945
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##   2.5%  97.5% 
## 2.7000 6.2025
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
##    [1] 3.7 4.8 5.5 4.4 2.1 5.3 5.3 5.9 3.3 4.1 4.3 5.2 4.0 5.3 3.9 5.9 4.8 4.9
##   [19] 5.4 4.5 6.5 5.7 4.2 5.1 5.6 5.2 2.6 5.4 4.9 4.6 3.6 3.9 4.7 4.9 5.6 5.6
##   [37] 2.8 3.5 4.1 2.4 3.5 4.6 5.0 5.1 6.0 4.2 4.7 5.6 2.0 4.8 3.3 4.6 4.2 4.7
##   [55] 5.1 2.5 3.9 4.4 3.6 5.8 3.8 4.6 3.1 3.0 6.0 5.6 5.3 4.6 4.4 4.7 4.0 3.9
##   [73] 5.6 5.6 4.2 5.5 4.2 4.3 6.0 5.5 4.0 4.2 5.3 4.9 4.9 4.5 4.5 5.9 4.2 4.9
##   [91] 5.0 4.0 4.3 5.9 6.0 4.5 4.0 2.9 3.1 5.1 4.0 5.6 5.1 4.2 5.9 3.8 5.7 5.4
##  [109] 5.1 6.0 3.7 4.2 4.2 4.7 5.8 5.0 4.5 5.8 4.8 3.5 5.9 6.1 3.7 4.4 4.4 4.7
##  [127] 3.7 5.0 4.5 6.1 7.0 3.4 3.6 5.0 4.6 2.8 3.5 4.6 4.6 4.7 6.6 5.0 4.6 4.4
##  [145] 3.6 3.9 5.4 5.7 4.9 3.7 3.5 6.7 5.3 4.4 4.3 4.7 4.5 4.0 4.5 5.0 5.1 7.6
##  [163] 4.8 3.0 3.7 5.5 6.3 4.6 4.0 5.6 3.6 3.3 4.8 3.9 5.9 3.3 4.9 3.3 4.1 4.2
##  [181] 3.5 4.9 3.4 4.7 4.3 3.6 6.0 5.3 5.0 3.4 4.6 2.9 5.4 5.5 5.3 4.2 3.7 4.5
##  [199] 4.1 3.1 4.7 4.7 4.5 2.8 3.1 4.7 4.4 6.2 4.5 5.3 4.4 5.0 4.0 3.6 2.1 2.9
##  [217] 3.2 3.8 3.8 5.3 6.1 3.6 3.9 3.8 3.8 3.7 2.8 4.2 5.0 3.9 4.0 4.4 4.5 4.1
##  [235] 4.9 4.9 3.9 5.8 5.9 6.0 4.3 5.7 4.3 4.5 7.2 5.2 5.4 4.0 4.6 5.5 4.6 3.0
##  [253] 5.1 3.5 3.8 4.5 3.4 4.3 5.2 4.7 4.3 3.1 5.3 5.7 3.7 5.0 4.8 4.1 4.4 5.8
##  [271] 2.7 5.5 4.8 4.8 3.6 3.6 5.2 5.4 4.9 4.8 4.5 2.6 5.1 4.7 4.4 5.3 4.0 4.4
##  [289] 4.6 2.5 2.9 5.1 5.1 4.4 5.9 5.5 5.5 4.0 4.1 6.1 4.6 3.5 5.2 4.7 4.7 4.9
##  [307] 4.5 4.4 4.0 5.9 4.7 6.2 5.4 4.9 3.2 3.8 4.7 3.9 5.3 3.1 4.6 5.6 4.0 4.0
##  [325] 4.7 3.9 4.8 4.7 5.4 4.0 4.1 4.0 3.9 5.1 5.8 5.1 5.1 4.5 4.0 4.9 4.8 4.7
##  [343] 4.4 4.4 4.1 3.8 5.5 4.2 4.7 3.8 4.1 5.7 5.2 3.8 4.1 3.8 4.6 3.7 4.1 4.5
##  [361] 5.0 6.7 4.3 4.4 4.2 4.1 4.8 3.3 4.0 4.9 4.6 4.8 5.7 6.3 5.1 3.6 5.8 4.7
##  [379] 4.4 4.1 4.3 4.5 5.0 4.4 3.5 4.6 2.3 4.3 3.5 5.7 3.6 3.8 5.1 5.6 4.2 4.2
##  [397] 4.5 5.3 4.4 5.2 3.6 3.3 4.4 4.8 4.5 3.6 4.0 5.1 3.9 4.4 4.9 3.9 4.0 4.8
##  [415] 3.5 3.6 4.0 5.7 4.4 5.7 3.8 5.2 5.4 5.3 5.5 5.9 3.7 4.3 3.7 3.4 4.1 4.8
##  [433] 4.2 3.9 3.1 5.1 5.9 5.1 5.5 4.0 4.3 5.2 4.0 3.1 4.2 3.6 5.2 4.3 4.8 3.6
##  [451] 4.9 3.5 4.1 3.7 3.6 4.1 6.5 4.4 3.4 4.6 5.4 3.9 4.0 4.9 3.9 3.2 2.8 4.0
##  [469] 6.4 4.4 5.3 5.3 3.3 3.5 4.4 4.3 5.2 4.7 2.8 4.4 3.8 5.6 4.9 3.6 4.3 4.2
##  [487] 3.7 3.6 2.5 5.6 5.7 3.8 5.1 2.4 5.4 3.9 3.7 3.0 5.7 2.8 3.5 2.7 3.4 2.7
##  [505] 6.5 5.3 5.3 4.7 4.8 4.1 4.8 5.4 5.7 2.5 4.8 3.6 5.1 3.9 5.5 5.0 3.6 5.9
##  [523] 4.9 3.9 5.0 4.9 3.5 5.7 4.7 4.4 5.5 5.4 4.7 2.7 4.7 5.0 4.1 4.6 4.7 3.4
##  [541] 3.2 4.4 3.9 2.8 5.8 4.8 4.5 4.6 3.4 5.3 6.2 3.9 5.1 3.4 3.9 5.5 4.7 5.6
##  [559] 4.2 3.1 5.4 4.3 4.9 4.8 2.0 4.9 5.0 3.9 4.6 3.8 3.5 4.5 6.5 5.9 4.1 6.0
##  [577] 5.6 4.7 6.2 4.1 5.4 6.4 5.9 4.0 4.1 3.2 6.2 3.8 4.1 2.8 5.4 5.1 4.0 5.2
##  [595] 4.0 5.9 5.7 6.0 4.6 4.1 3.8 4.4 5.1 2.7 4.3 5.7 5.1 4.7 4.3 3.1 3.6 4.8
##  [613] 5.1 3.3 4.8 4.2 5.8 4.6 3.3 4.5 5.2 5.4 3.9 3.7 4.2 3.0 4.1 4.3 4.7 4.4
##  [631] 5.2 4.7 4.9 3.9 5.4 4.0 5.0 3.0 5.3 5.3 5.2 4.4 3.7 5.4 5.2 4.1 4.7 4.8
##  [649] 6.4 4.4 3.5 4.4 5.1 5.3 4.7 5.0 4.1 4.2 5.2 4.6 5.4 3.3 4.0 4.4 4.4 5.6
##  [667] 6.2 3.2 4.6 5.3 3.2 3.2 5.3 4.2 4.2 5.4 4.6 6.3 4.1 4.6 4.0 4.3 4.2 5.9
##  [685] 4.4 4.6 4.2 4.1 3.4 3.5 4.4 4.5 5.3 4.1 4.3 5.0 3.7 5.6 4.7 4.9 5.2 6.1
##  [703] 4.2 4.7 4.8 4.1 4.2 5.3 4.8 3.7 3.2 5.1 4.0 3.9 4.1 4.7 4.9 5.5 5.0 4.2
##  [721] 4.5 3.4 5.1 6.6 3.4 3.9 5.2 5.2 3.6 5.2 4.3 5.5 2.9 3.4 4.5 4.8 5.3 4.1
##  [739] 4.4 5.5 5.4 4.6 4.1 4.9 5.5 4.3 4.4 4.2 4.4 4.9 3.1 4.7 4.2 4.0 3.8 4.9
##  [757] 4.9 3.7 5.4 4.2 3.6 7.4 3.3 4.6 4.6 5.6 4.9 4.4 2.6 5.3 3.8 5.5 4.2 3.7
##  [775] 4.7 4.7 3.0 4.2 5.6 3.3 4.5 4.3 4.4 5.2 4.3 5.2 3.0 3.9 5.4 4.5 5.1 6.1
##  [793] 5.1 4.3 4.6 4.4 3.7 4.5 4.7 4.7 4.8 4.7 2.8 3.5 6.6 5.8 3.9 4.6 4.5 4.1
##  [811] 4.9 3.8 3.6 4.3 5.2 3.5 3.8 3.9 4.9 4.9 5.9 3.9 5.3 4.6 4.1 5.5 4.8 5.3
##  [829] 4.8 4.0 4.7 7.8 4.9 4.3 4.1 4.2 3.6 3.5 3.8 4.5 2.9 2.7 6.3 3.3 3.5 4.0
##  [847] 4.1 4.0 4.5 5.4 3.9 3.8 5.2 3.9 5.2 3.8 5.1 5.5 3.7 5.1 5.3 5.6 4.8 4.8
##  [865] 5.1 3.7 4.6 4.8 4.6 4.1 3.9 5.0 4.6 5.2 3.5 3.9 4.9 5.2 5.2 4.2 5.8 3.1
##  [883] 4.1 5.9 4.1 4.1 4.6 5.3 5.5 6.2 4.4 5.8 4.3 4.7 4.5 4.4 6.1 3.5 5.9 5.4
##  [901] 3.0 4.9 5.6 4.5 6.2 2.8 4.7 3.9 4.2 3.4 4.7 4.6 4.9 4.7 4.1 3.3 3.4 4.4
##  [919] 5.4 4.2 5.8 5.1 5.4 5.0 3.7 5.0 5.7 4.0 3.6 2.9 4.3 3.0 5.9 5.5 4.3 3.2
##  [937] 2.6 5.5 4.9 4.4 5.2 3.4 4.4 4.8 6.0 5.7 4.7 6.4 3.2 2.6 4.3 4.7 5.0 3.3
##  [955] 3.2 5.2 5.2 4.9 4.0 3.6 3.9 5.0 5.5 4.0 3.5 5.1 4.7 6.1 4.6 5.6 5.1 3.7
##  [973] 4.8 4.0 5.6 4.3 6.0 4.9 5.1 4.3 4.9 4.3 4.1 4.2 5.1 5.2 5.0 5.5 4.7 2.9
##  [991] 3.5 4.6 4.6 4.0 4.3 3.3 3.4 4.2 4.7 2.7
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
##   2.8   6.2
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
##    [1] 4.0 3.6 5.8 4.2 4.0 5.1 6.1 5.8 4.5 4.9 4.9 5.0 4.3 4.4 3.1 2.9 3.8 4.6
##   [19] 3.0 4.1 2.3 3.7 4.9 5.6 4.9 3.7 4.4 4.6 5.3 4.8 3.6 4.1 2.9 5.9 3.5 4.1
##   [37] 5.5 4.8 3.7 4.2 3.9 4.4 3.5 4.6 3.9 4.3 3.5 3.2 4.4 4.3 4.3 3.7 5.8 4.5
##   [55] 4.5 4.1 6.3 4.2 3.6 3.9 3.9 4.7 3.9 4.6 4.3 5.1 4.6 3.6 5.5 4.6 4.9 4.5
##   [73] 4.0 5.1 5.5 5.3 3.4 4.2 4.7 4.2 3.8 2.8 4.3 3.7 3.5 4.1 4.5 3.6 3.9 3.3
##   [91] 4.3 3.8 3.4 4.2 5.1 4.6 5.5 3.8 4.1 5.1 4.4 4.1 4.5 4.6 4.1 4.5 6.1 5.5
##  [109] 2.6 2.7 5.0 6.2 2.2 5.6 2.8 5.4 4.8 3.6 5.2 4.7 3.9 4.0 3.9 4.5 4.9 4.3
##  [127] 4.8 3.4 5.3 4.1 4.6 3.5 4.5 5.9 4.8 4.5 4.5 4.3 3.4 6.2 4.0 4.4 4.8 4.4
##  [145] 2.9 4.5 4.3 4.7 5.2 5.2 3.1 4.9 4.4 4.4 4.9 4.9 4.1 6.1 4.1 4.9 3.7 4.3
##  [163] 4.2 5.5 4.8 3.0 3.5 4.0 4.2 5.3 4.4 4.5 5.7 4.6 4.2 4.2 2.9 3.6 5.7 4.9
##  [181] 2.8 2.4 3.2 5.1 4.5 4.2 5.0 2.7 4.1 4.1 3.2 5.0 4.0 6.4 5.8 4.6 5.0 5.0
##  [199] 3.2 5.1 3.5 4.3 3.9 5.4 4.2 6.1 3.7 4.0 6.6 4.6 4.9 2.6 4.5 5.5 3.8 5.2
##  [217] 3.7 3.4 3.9 5.2 4.3 3.8 3.9 1.8 5.4 4.2 5.1 5.6 4.4 6.1 4.1 4.1 6.1 5.2
##  [235] 3.7 4.9 4.5 5.5 3.3 4.5 5.6 3.7 6.2 5.4 3.9 3.2 6.4 5.0 4.5 5.5 5.3 3.7
##  [253] 4.5 5.4 5.6 3.3 4.2 4.3 3.3 4.4 5.7 2.8 5.4 3.1 6.5 5.8 4.6 4.7 5.6 3.4
##  [271] 5.2 5.7 4.0 6.6 5.0 5.0 4.9 4.5 3.0 5.2 3.6 3.5 4.3 4.3 3.3 3.1 3.6 4.3
##  [289] 4.7 5.8 5.2 5.4 3.0 4.2 4.4 4.8 4.2 3.7 4.3 4.7 3.6 4.0 5.7 5.0 5.6 4.5
##  [307] 5.5 5.7 6.0 3.9 5.1 4.7 5.0 4.8 4.9 5.1 4.8 5.1 5.4 3.0 3.6 4.1 4.7 4.2
##  [325] 3.2 5.4 5.0 3.3 5.1 3.0 3.9 2.2 3.7 3.9 4.0 3.9 3.9 6.5 4.9 4.7 4.8 4.2
##  [343] 5.4 4.4 5.4 4.9 4.5 6.1 5.5 3.4 2.9 4.7 4.5 4.2 4.3 4.5 2.9 5.6 4.6 4.2
##  [361] 5.6 5.1 2.9 2.9 3.3 4.7 3.8 5.2 4.3 5.2 2.9 6.4 5.6 3.5 5.3 3.9 4.8 4.8
##  [379] 3.3 4.5 5.7 4.8 3.2 4.2 4.2 4.2 3.8 6.4 3.6 5.4 5.1 6.2 6.4 6.0 5.1 3.3
##  [397] 5.0 5.0 3.2 6.0 2.5 5.3 4.7 4.9 4.0 4.9 4.6 5.8 3.7 3.1 4.0 4.2 4.5 4.4
##  [415] 3.8 5.2 4.4 3.7 2.7 3.6 2.6 5.7 5.8 3.3 5.0 4.8 2.2 3.3 4.4 5.5 6.7 5.0
##  [433] 4.1 4.3 5.7 5.4 4.5 4.7 3.9 5.8 4.3 5.5 4.7 5.5 6.3 4.0 4.3 4.2 4.1 4.0
##  [451] 4.3 4.3 4.7 4.4 4.2 5.6 4.2 5.4 6.4 4.2 4.9 4.7 3.3 4.5 5.0 4.5 5.0 5.5
##  [469] 4.4 6.3 4.1 3.9 4.1 4.8 4.4 4.8 5.2 4.7 6.1 5.3 6.1 5.6 5.1 4.7 4.6 5.2
##  [487] 5.4 4.5 3.7 3.2 3.5 5.6 6.2 2.6 3.2 5.4 3.3 4.3 4.5 5.7 4.3 6.2 3.6 4.8
##  [505] 5.3 3.4 4.8 5.4 3.8 4.4 4.9 3.2 3.9 4.0 3.0 3.9 3.6 6.0 4.6 5.6 4.8 5.5
##  [523] 2.1 4.1 4.9 4.4 5.1 4.1 3.4 3.7 5.1 4.3 4.7 5.2 3.8 4.7 3.0 5.1 3.7 4.7
##  [541] 5.4 5.4 4.6 4.4 4.0 3.8 5.8 4.4 4.3 4.9 5.7 4.1 2.7 4.5 5.5 6.6 5.3 5.0
##  [559] 4.4 4.0 4.1 5.4 2.9 5.4 4.2 3.9 4.2 5.0 4.5 2.9 3.7 4.5 4.8 4.6 3.3 3.9
##  [577] 3.5 5.3 3.0 2.1 4.2 3.7 5.8 4.4 5.8 2.2 3.8 4.0 3.9 4.3 4.2 4.9 3.0 3.7
##  [595] 5.8 2.8 6.0 5.1 6.6 4.6 4.2 3.9 3.2 3.8 4.7 6.2 5.4 4.5 3.2 4.2 6.0 4.5
##  [613] 3.6 5.8 5.6 4.8 5.1 5.1 3.1 5.5 4.1 5.3 4.5 3.3 5.0 4.6 1.9 3.5 4.7 3.8
##  [631] 3.3 5.2 5.1 5.2 3.8 4.8 4.6 3.4 3.3 4.9 3.3 4.8 5.2 5.0 4.2 4.4 4.5 6.0
##  [649] 4.3 3.9 5.1 4.9 6.1 5.2 4.2 5.7 4.9 4.2 5.3 4.7 4.7 4.2 3.7 4.4 3.8 5.3
##  [667] 5.7 2.2 5.4 5.7 4.8 5.3 4.7 5.9 3.9 3.6 5.8 4.1 4.8 3.8 4.6 4.2 4.9 3.7
##  [685] 2.4 3.6 4.7 3.2 6.5 5.1 2.9 6.2 4.0 4.6 4.3 5.4 4.4 5.5 4.6 4.9 4.6 3.3
##  [703] 5.3 3.6 3.1 4.5 5.9 6.0 4.8 4.8 5.4 5.1 3.4 4.5 3.8 5.6 6.1 5.4 3.1 2.9
##  [721] 3.6 3.7 3.4 4.3 4.9 3.9 3.7 5.4 2.8 4.5 5.1 4.1 3.6 3.4 3.1 2.5 3.9 4.4
##  [739] 4.6 5.1 3.3 3.9 5.4 3.9 4.7 2.6 6.2 5.0 4.7 3.4 4.2 3.4 5.5 4.6 4.5 4.9
##  [757] 4.9 4.4 4.9 5.0 4.6 4.4 3.0 4.4 4.6 3.7 4.7 5.0 5.8 4.4 4.7 6.2 4.5 5.1
##  [775] 3.9 3.3 7.4 3.5 3.4 3.5 5.2 5.2 4.5 6.0 2.2 4.3 4.4 4.1 3.8 4.3 5.3 5.6
##  [793] 4.7 5.3 4.6 3.6 5.0 3.8 3.2 4.1 4.4 5.3 5.1 5.9 4.0 5.0 5.3 5.5 5.3 3.5
##  [811] 4.6 4.8 4.3 4.7 3.4 4.6 4.3 4.9 4.0 3.5 5.7 5.4 5.3 5.2 5.7 3.8 4.6 3.6
##  [829] 2.8 4.4 4.3 4.7 4.2 4.7 5.2 2.8 4.7 4.5 4.7 4.0 4.8 2.2 6.0 4.0 5.2 2.3
##  [847] 4.3 4.5 4.5 2.6 5.6 3.3 3.6 4.9 5.3 4.1 3.3 5.5 4.8 4.8 4.7 4.2 5.2 5.6
##  [865] 6.8 5.7 4.8 5.6 4.2 4.5 5.4 4.0 2.6 4.6 5.8 5.3 5.2 5.0 4.2 4.0 3.7 4.7
##  [883] 3.1 6.1 3.5 4.5 5.3 4.4 4.0 4.6 4.6 3.8 4.2 3.8 4.3 4.8 4.8 3.3 4.6 4.1
##  [901] 4.4 4.1 3.0 3.8 2.9 5.7 3.7 5.1 2.9 3.9 5.1 4.9 5.5 3.5 5.6 3.5 6.0 4.6
##  [919] 3.4 3.2 5.4 4.8 4.4 5.3 5.0 5.7 4.7 5.5 4.4 5.3 3.0 5.0 7.4 5.5 3.2 6.8
##  [937] 6.6 4.3 4.5 4.4 3.8 4.6 6.2 6.2 3.1 3.5 3.7 4.9 4.2 5.8 3.4 4.6 5.2 5.1
##  [955] 7.0 4.7 4.7 3.9 4.4 5.5 6.0 5.0 3.6 2.8 3.7 3.7 4.1 3.2 4.0 3.5 3.8 5.4
##  [973] 5.0 4.6 4.9 4.3 4.0 3.3 4.9 4.1 4.4 2.7 3.4 4.1 2.9 3.7 3.3 6.5 5.7 4.9
##  [991] 3.6 4.7 5.1 4.7 4.0 5.0 3.8 6.1 5.0 3.9
## 
## $func.thetastar
## [1] -0.0205
## 
## $jack.boot.val
##  [1]  0.52882883  0.46060606  0.27636888  0.17404130  0.04573171 -0.13057851
##  [7] -0.22286501 -0.26942857 -0.48600000 -0.55454545
## 
## $jack.boot.se
## [1] 1.073499
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
##    [1] 3.8 5.2 5.4 4.1 6.0 2.9 6.2 3.8 4.8 4.9 5.9 4.9 4.3 4.4 2.7 3.8 5.7 3.6
##   [19] 5.6 5.4 5.0 5.8 4.9 5.1 4.2 5.6 4.1 4.8 5.6 3.1 4.6 3.7 3.6 3.6 3.5 5.1
##   [37] 4.1 4.8 4.5 4.1 4.2 4.1 4.3 4.1 5.0 4.9 4.5 3.7 3.9 3.3 6.5 5.0 4.1 4.6
##   [55] 3.9 4.5 4.8 5.0 4.9 5.1 4.4 5.0 4.7 3.2 5.0 4.7 5.0 3.4 2.5 4.6 4.8 5.1
##   [73] 3.0 5.3 5.8 5.1 4.5 4.5 4.5 4.6 3.5 4.4 6.6 4.0 3.5 3.5 5.9 4.5 4.0 5.0
##   [91] 5.6 5.6 4.6 3.6 5.2 4.9 3.8 4.6 3.2 3.9 4.4 4.0 5.0 4.3 4.0 4.9 3.9 4.1
##  [109] 4.1 4.3 4.6 4.0 4.5 4.0 4.4 4.0 3.7 4.5 3.8 4.7 5.7 4.7 4.6 5.0 3.6 4.8
##  [127] 4.8 6.1 4.4 4.6 3.7 5.7 4.3 3.8 3.8 3.2 5.1 3.2 3.4 4.7 4.2 5.1 5.1 3.5
##  [145] 4.9 5.7 2.6 4.5 5.3 4.8 4.9 4.2 5.6 4.5 5.2 4.8 3.9 4.0 4.5 5.0 4.1 5.8
##  [163] 3.8 4.4 5.4 4.2 4.6 5.0 5.5 4.6 3.4 6.2 4.6 4.3 4.3 4.4 5.4 4.9 3.7 4.5
##  [181] 3.2 4.7 4.5 4.1 4.0 4.7 6.8 3.8 3.5 3.5 6.2 3.9 5.5 4.1 4.9 3.4 3.7 3.2
##  [199] 3.8 5.2 4.8 3.8 6.4 5.0 3.0 6.5 4.7 5.5 4.9 4.1 3.6 4.7 6.2 3.4 5.7 3.8
##  [217] 4.2 5.4 4.0 4.6 4.9 4.7 2.2 4.2 5.2 4.9 3.2 6.7 3.9 5.9 4.4 5.5 4.5 4.7
##  [235] 4.0 4.2 4.4 4.0 5.3 4.0 5.6 6.0 4.3 6.1 5.3 3.4 6.3 5.1 5.0 4.9 4.5 3.3
##  [253] 4.1 5.1 3.1 3.4 6.1 3.3 3.8 4.0 4.2 4.8 5.3 4.4 5.4 4.8 5.2 3.8 4.0 4.7
##  [271] 5.0 4.9 4.4 4.3 4.3 6.0 5.2 4.3 5.2 4.9 4.9 2.4 3.9 5.3 5.7 4.7 4.8 4.5
##  [289] 4.7 4.8 4.7 5.0 4.0 4.8 5.4 4.1 4.8 7.2 3.2 3.7 4.3 3.3 5.5 2.8 4.4 5.2
##  [307] 5.2 4.1 3.3 5.7 3.9 5.7 6.3 4.6 4.0 3.6 3.4 2.8 4.4 4.3 4.7 5.0 5.2 4.0
##  [325] 3.3 3.4 5.3 4.3 4.4 4.4 4.6 4.8 3.6 5.0 5.2 5.0 3.0 6.2 4.3 5.6 3.8 4.1
##  [343] 4.5 3.7 4.7 4.3 3.7 5.3 5.6 4.4 4.6 4.6 4.6 4.6 5.1 4.3 5.1 5.1 4.1 3.8
##  [361] 5.8 4.9 3.7 5.2 4.9 5.3 5.4 4.0 4.0 4.1 5.1 4.4 4.2 4.8 4.9 4.8 5.7 4.1
##  [379] 4.1 3.6 3.9 3.4 2.6 3.2 5.9 5.9 4.0 3.8 3.1 3.3 3.0 4.3 3.7 4.5 3.4 5.8
##  [397] 5.0 4.2 4.9 2.6 4.7 4.8 5.6 4.9 4.6 2.9 6.5 4.0 4.7 4.8 4.0 4.6 5.3 4.1
##  [415] 3.7 6.1 4.6 4.1 4.8 3.6 4.6 4.6 5.1 5.3 5.7 4.9 2.4 4.3 5.5 5.6 4.9 2.8
##  [433] 5.0 3.1 2.8 4.7 4.3 5.9 4.4 5.3 3.0 4.5 2.8 3.8 4.5 3.4 3.7 3.5 4.5 5.4
##  [451] 3.8 6.1 4.7 4.0 5.8 5.3 6.2 6.3 3.3 4.3 5.6 4.8 5.3 5.7 3.4 4.4 4.5 3.9
##  [469] 4.3 3.6 4.4 4.3 4.4 4.8 3.0 3.5 3.0 4.4 7.0 4.8 4.5 4.1 4.5 4.6 4.8 4.0
##  [487] 4.4 4.2 4.7 4.9 4.8 4.8 4.3 4.1 3.1 2.7 3.6 4.9 4.5 4.4 4.2 5.3 2.5 4.9
##  [505] 4.7 5.5 4.7 4.5 4.3 3.3 4.7 4.4 4.3 3.9 4.1 4.4 4.4 4.8 4.7 3.7 5.4 3.5
##  [523] 5.5 5.0 5.4 3.9 4.0 5.7 3.0 4.6 5.0 6.1 4.9 3.5 3.6 4.0 4.2 5.0 4.8 3.8
##  [541] 3.2 4.3 5.5 5.3 4.0 4.6 4.7 4.6 4.4 4.9 5.3 3.7 3.6 4.5 6.1 4.8 4.3 5.7
##  [559] 5.7 4.6 4.5 4.3 4.3 3.9 6.4 5.7 4.7 4.0 4.7 4.1 3.1 6.0 3.1 5.7 3.3 4.9
##  [577] 5.0 3.2 4.0 5.8 6.2 5.9 5.1 5.0 3.2 5.2 4.8 3.5 5.4 5.0 4.9 5.0 3.4 5.4
##  [595] 3.9 4.8 3.0 4.5 6.4 4.0 3.9 4.5 4.1 5.5 3.6 3.4 5.3 3.6 5.6 4.9 4.9 4.5
##  [613] 4.4 4.6 5.9 5.4 2.8 4.1 5.2 4.1 5.9 4.9 4.0 4.8 4.0 3.9 5.6 4.3 4.8 2.6
##  [631] 4.4 4.4 4.8 5.2 5.9 3.9 5.2 4.3 4.1 4.7 3.9 6.1 4.6 4.8 5.3 3.9 4.7 4.4
##  [649] 5.1 3.5 4.8 4.3 5.7 5.6 3.7 4.3 3.1 3.6 3.9 4.1 2.6 4.5 4.8 2.2 6.1 5.5
##  [667] 3.9 2.5 6.3 4.5 5.1 2.4 3.2 5.4 4.5 5.9 5.7 5.4 3.9 3.5 4.0 4.0 5.5 4.2
##  [685] 3.6 3.1 5.0 6.3 4.5 5.7 5.7 4.3 5.2 4.1 4.4 5.1 4.8 4.4 5.0 5.7 4.7 3.2
##  [703] 5.3 5.3 3.4 5.8 3.7 2.4 4.7 5.4 1.9 4.3 3.8 4.3 3.7 3.6 5.7 4.4 4.6 5.9
##  [721] 5.4 4.4 5.2 3.1 4.6 4.9 4.9 3.1 3.8 4.7 5.5 6.3 4.4 3.7 4.9 6.1 4.3 4.4
##  [739] 3.8 3.9 4.0 5.6 5.6 3.2 4.0 4.1 5.4 5.2 4.3 5.4 3.4 5.1 4.5 5.1 4.3 4.3
##  [757] 3.3 4.6 4.5 4.7 3.7 4.0 4.0 3.5 5.4 5.7 2.5 3.9 4.2 4.7 5.3 3.9 5.5 5.4
##  [775] 3.9 2.9 4.6 5.7 4.2 5.0 2.5 5.0 3.5 6.7 4.2 3.4 4.1 4.1 4.8 3.6 4.8 4.5
##  [793] 5.3 5.4 4.1 2.3 4.6 6.4 4.5 5.7 4.0 4.7 4.7 3.5 5.2 4.3 5.3 4.4 3.3 3.9
##  [811] 3.2 3.9 4.0 3.5 5.4 4.8 4.1 3.9 5.3 5.6 5.6 4.8 3.3 4.1 3.6 4.1 3.6 3.8
##  [829] 3.8 5.1 5.5 4.9 5.4 4.4 5.0 6.1 4.2 3.5 4.4 5.3 4.9 3.6 3.8 4.1 5.9 4.1
##  [847] 4.7 3.7 4.9 3.5 3.7 2.6 4.8 4.1 5.6 4.9 2.8 4.1 3.8 4.1 5.1 4.8 5.2 4.9
##  [865] 5.9 3.9 2.4 5.5 3.9 3.4 4.0 6.0 5.3 4.2 4.8 5.5 3.5 6.0 4.4 4.1 5.3 4.5
##  [883] 4.4 7.0 3.5 4.1 6.7 5.0 4.8 4.1 6.9 4.6 5.5 5.4 6.9 3.2 3.5 4.4 4.5 4.2
##  [901] 5.3 5.2 5.0 4.3 4.5 3.2 4.8 3.9 5.2 5.5 3.4 4.3 6.8 3.7 3.8 4.5 4.9 5.4
##  [919] 4.3 4.5 4.7 5.4 3.0 5.7 5.0 3.4 5.8 6.1 4.0 4.9 3.4 4.5 3.0 3.9 5.3 4.5
##  [937] 6.0 5.9 3.6 5.6 4.5 3.1 6.2 4.7 3.9 3.8 4.4 3.7 4.5 5.0 5.3 4.0 4.9 4.3
##  [955] 5.5 5.2 4.9 2.7 5.2 4.9 3.9 4.5 3.4 3.0 4.8 4.4 4.7 3.9 4.0 2.9 5.4 4.3
##  [973] 6.9 2.8 3.5 4.1 4.9 3.6 3.5 4.5 5.2 4.8 4.4 4.0 5.3 3.7 3.2 3.5 3.6 4.3
##  [991] 5.5 5.8 4.7 4.4 5.0 4.7 4.3 3.5 6.1 2.8
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.5 5.4 5.3 5.2 5.1 5.0 4.9 4.8 4.7 4.5
## 
## $jack.boot.se
## [1] 0.911921
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
## [1] 0.806202
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
##   1.937906   2.823411 
##  (0.803231) (1.334488)
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
## [1] 0.4491101 1.6494614 1.3156445 0.3542904 0.6753906 0.8322717
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
##    [1]  0.539918400  0.254643744  0.523515798  0.891798642 -0.475352256
##    [6]  0.717021543  0.516701062  0.348816182  0.161855686  0.494332459
##   [11]  0.239384861  0.404255814  0.119467063 -0.319200011  1.526133758
##   [16]  0.079305981  0.041379116  0.715241242 -0.401698782  0.429131800
##   [21]  1.578199397  0.920747211  0.141031056 -0.842843943  0.908841463
##   [26]  0.394477482  0.748993322 -0.333817659  1.173019594  1.352607709
##   [31]  1.132251941 -0.320185579 -0.154525964  1.186049367  0.846089128
##   [36]  0.744845950 -0.286148558  0.502818821  0.334574654 -0.267065661
##   [41]  0.340671396  0.790497776  0.902817486  0.399599047  0.470642575
##   [46] -0.200573892  0.965293342  0.477862405  0.623090110  1.404214388
##   [51]  0.106612580  0.777937911  0.648239777  0.911914018  1.548956814
##   [56]  0.027551947  1.648581571  0.261743459  0.396196799  0.507272541
##   [61] -0.225186956  0.163181771  0.863655614  1.523024816  0.868392632
##   [66]  0.893725070  1.116246431  0.712184069 -0.199869886  0.278508737
##   [71]  0.549990302  0.566006302  0.916008410  1.000526902  1.150561148
##   [76]  1.203802712  0.055145496 -0.288983304  0.744766961  0.689906616
##   [81] -0.208690778 -0.097681092 -0.200274317  1.042445199  0.900144713
##   [86]  0.474709564 -0.545334851  1.124721322  0.324526599  0.199441705
##   [91]  1.140938645  0.480058562  0.587878178  0.223715116 -0.884355957
##   [96]  0.734850478  0.237002137  0.815033838  0.404255814  0.770483798
##  [101]  1.059838362  0.918585128  0.826598089 -0.359661539  0.828712602
##  [106]  0.194756777  1.141307210  0.665308519  0.434943348  1.483228636
##  [111]  0.666217761  0.688858517  0.261247313  0.819740978  0.321145978
##  [116]  1.196458147 -0.825404976 -0.430725884  0.467724184  1.863902851
##  [121]  1.498471009  0.324753481  0.062237743 -0.538149027  0.630035518
##  [126]  0.297518239  0.078113797  0.672047431  0.614647173  0.965083186
##  [131]  1.123854346  1.173432237  0.091823602  0.702636490  1.531126240
##  [136]  0.361310373  0.934521937 -0.105306790 -0.174009100  1.134013473
##  [141]  1.050020945  0.780801756 -0.823601683  0.699673828  0.310280787
##  [146]  0.243886689  0.435398435  0.219048365  0.731328822  0.189156785
##  [151]  0.485216487  0.795809260  0.348494806  0.614629033  0.720972903
##  [156]  1.427589842  1.046802288 -1.164988211  0.033886052  1.077879242
##  [161]  0.142354548  1.188865414  1.053548841  0.427444323  0.246850839
##  [166]  0.072955848  0.913052458  0.383609217  0.892353455 -0.634696118
##  [171]  0.385748057  1.044503321  0.111790862  1.006933928  0.826034761
##  [176]  0.226043568  0.413819622  1.315357615  0.945305600  0.084064359
##  [181] -0.242752009  0.922199845 -0.276773562  0.525051985  0.003411710
##  [186]  0.878155000  0.686483395  0.936195600 -0.021794119 -0.293709481
##  [191]  0.167180181  0.383550194  0.073871608  0.819077845  0.763604201
##  [196]  0.602916193 -0.964775414  0.399672719  1.086306673  1.329224513
##  [201]  0.792298168 -0.081361852  0.329519098  0.923152586  1.648558820
##  [206] -0.058472706  0.994407501  0.321752580  0.695583401 -0.790513170
##  [211] -0.342547179  0.271889421  0.540733715 -0.204722519  1.259896459
##  [216]  0.108008126  0.134038896  0.999591432  0.978508547 -0.029517013
##  [221]  0.033572044  0.902064705  0.979250769 -0.073910372  0.487999371
##  [226]  0.735672789 -0.737887699  0.392704533  0.422821156  0.176885431
##  [231]  0.593157204  0.339581860  0.022487282  0.883735698  0.830075197
##  [236] -0.298496451  0.769094337  1.336344444 -0.275528681  0.589325660
##  [241]  0.771248477  0.280736526  1.088020223  0.358071127  0.889828675
##  [246]  1.130310453 -0.753525485 -0.042447484  0.571937421 -0.051875077
##  [251]  1.102917954  0.862508199  1.450898347  0.326735842 -0.013879776
##  [256]  1.006699290  0.226233827  1.304387005  1.344339714  0.711919339
##  [261]  0.539139549  0.990594893  0.032181029  0.607854575  0.604326638
##  [266]  1.236518595  0.793399029  0.229928626  0.581107494  1.144906146
##  [271]  0.588826659  0.413308037 -0.215397446  0.209413866  0.828067000
##  [276]  0.537272241  0.630537355  0.366541559  0.385049685  0.519279974
##  [281]  0.698767813  1.152241314  0.424601161  0.227222617  1.188454032
##  [286]  0.442070959  0.739647987  0.925599631  0.142932925  0.725171381
##  [291]  0.428366205  1.191095793  0.464565504  0.486387327  0.646734360
##  [296]  0.212192980 -0.193192325  0.132181935  1.486662877 -0.146519024
##  [301] -0.159464238  0.389813464 -0.554749385  0.541364972  1.598943814
##  [306]  0.562424583  1.782923323  0.934703107 -0.393387526  1.583123044
##  [311]  1.302391460  0.838236622  0.290085519  0.455677373  1.067539238
##  [316]  0.318629488  0.944880852  0.259088649  0.499166381  1.071946979
##  [321]  0.273992437  0.193157438  0.522750801 -0.305690758  0.414683810
##  [326] -0.473287692  0.809659507  0.768377953  0.914429514  0.401483149
##  [331]  0.134300953  0.936418236  0.629167101  2.126538236  0.777298667
##  [336]  1.012323414  0.536479211  0.265228466  0.347764963  0.545089689
##  [341]  0.010236809  0.481197989  0.497467770  0.986361245  0.482699081
##  [346]  0.182914580  0.804603987 -0.570279219  0.970211823  1.113392281
##  [351]  0.783147475  0.265102412  0.851039547  0.136682076  0.327497234
##  [356]  0.646240582  0.457349665  0.792142907  1.012282363  0.157743603
##  [361]  0.863492205  0.665583043  0.548756673  0.664996060  0.130028577
##  [366]  0.834176527  0.154219767  0.032131018  0.638006009 -0.052479612
##  [371]  0.706238549  1.007294984  1.005140811  0.546713445 -0.387735488
##  [376]  0.758720820  0.893770526  1.530537122  0.153094845  0.676613445
##  [381]  0.427480236 -0.463416341  0.300670767  0.920406406  0.624495318
##  [386]  0.796264350  0.450538358  0.994222057  0.155227646 -0.094969314
##  [391]  0.034051138  0.520688911 -0.109374900  0.755688924  0.926937525
##  [396] -0.870422916  1.274011554  1.346247247 -0.528572565  1.563702970
##  [401]  1.287026710 -0.224912295  0.060016876 -0.043379259  0.236322198
##  [406]  0.858058416  0.730741817 -0.377744143  1.681395333  0.056140832
##  [411]  0.678793284  0.733497784  0.136682076  0.472030951  0.738397600
##  [416]  0.272920206  1.087337107  0.470914079  0.446801547  1.229838964
##  [421]  0.163819955 -0.646727865  1.428569021 -0.166826355  0.865174686
##  [426]  0.826017895  1.479299183  0.949242317  1.697967678  0.983609611
##  [431]  0.164914403  0.285303498  0.610155814  0.478800615 -0.253240754
##  [436] -0.511621600  0.413011789 -0.471831339  0.768501961  0.083886128
##  [441]  0.645994630 -0.616493765  1.340231641  0.495723021  0.880110135
##  [446]  0.308501821  1.027325660  1.381592590  0.783013477  0.735651161
##  [451]  0.984797456  0.570176446  0.020367256  0.813742815 -0.837064260
##  [456] -0.466772042  0.650601304  1.499147010  1.419572566  0.112102837
##  [461]  0.972892819  0.104327344 -0.141050706 -0.718711375  1.359129378
##  [466]  0.465164772  0.493226476  0.753172950  0.346964995 -0.439071498
##  [471] -0.808852935  0.607677714  0.235642881  1.090259063 -0.471252857
##  [476]  1.052615289  0.841796649  0.501462667  0.902531310  0.430424113
##  [481]  0.660336468  0.209443865  0.358994813  0.459603839  1.010726617
##  [486]  0.944963663  0.712336270  1.279100453  1.522231924 -0.260411816
##  [491]  1.133574604  1.293319644 -0.157954062  0.491568625  0.609934185
##  [496]  0.749690323  0.863472185  0.202156741  1.034612727  0.714735026
##  [501]  1.319045909  0.542044613  1.330669861  0.424111643  0.412002339
##  [506] -0.074224553  1.809746536 -0.109060442  0.806079919  0.907448290
##  [511]  0.453787264  0.898387769 -0.103538702  0.819259969  0.396875667
##  [516]  0.273692023 -1.216985083  1.390535683  0.824196192 -0.780504383
##  [521]  0.700905505  0.904656077  1.168784034  0.636819137  0.736392811
##  [526]  0.504087612  0.307268979  1.187920171  0.205093645 -0.515185240
##  [531]  0.189214046  0.084294331  1.493454278 -0.232828326  0.594658059
##  [536] -0.297936046  0.153172484  0.190948140  0.550131764  0.400319564
##  [541] -0.767861648  1.757704125 -1.853006994  1.146169993  1.209277001
##  [546]  0.632826643  0.154299417  1.041090813  1.006699290  1.023986143
##  [551]  0.526766062 -1.630399966 -0.647485750  0.250516605  0.629167101
##  [556]  1.098169489  0.613318330  1.785284978 -0.485526456  0.476544080
##  [561]  1.028213177  0.204838833  0.072361875  0.465033895  0.870659044
##  [566]  1.428714384  1.113641921 -0.197501095  0.756459554  0.217427993
##  [571] -0.014878970  0.354290676  0.441089903  0.628136608  0.807208756
##  [576] -0.831990926  0.883536697  1.025341666  1.290959266 -0.109935949
##  [581]  0.327644088  0.226544375  0.233938809 -0.022069367  0.803810533
##  [586]  0.914429514  0.824734117  1.010726617 -0.136347638  1.349908480
##  [591]  0.784975192  0.283493869  0.952867665  0.260179406  0.030570910
##  [596] -0.546881277 -0.446246352  0.343556918  0.431670126  1.277205444
##  [601] -0.171258936  0.609776134  0.031174280 -0.180250645  0.901801493
##  [606] -0.218995114  0.754619597  0.393868441  0.641986244  0.438655068
##  [611]  0.979971527  2.288116854 -0.215409322  0.182188016 -0.043845830
##  [616] -0.257795499  1.771909732 -0.489155829 -0.137487563  0.622052549
##  [621]  0.839256758  0.628127965 -0.301625983  0.508466026  1.096244676
##  [626]  1.647148472  0.843349098  0.830652621 -0.866104059 -0.811972733
##  [631]  0.563197415  0.896401704  0.545487130  0.704149046  0.456226152
##  [636]  0.771913601  0.995381748  0.332024892  1.379755005  0.507403608
##  [641]  1.086250913  0.394107829  0.951496993  0.754445580  0.759743859
##  [646]  0.640468559 -0.089499410  1.008796252  0.433827128  0.414636631
##  [651]  0.674150973 -0.514042946  0.563926238  1.881247913  0.913753900
##  [656] -0.567583617  0.432524970  0.367268004  0.310353516  0.753637128
##  [661] -0.399104926  0.318680167  0.146903868 -1.302743823  0.929362281
##  [666]  0.098675794  0.634935875  0.410401862  0.718520784  0.566967735
##  [671] -0.092255802 -0.387133011  1.544119872  1.480326288  1.260568401
##  [676]  0.780753033  0.367584910  0.325141099  0.720858503  0.469784429
##  [681] -0.021148389  0.639712690  0.450879036 -0.525653433  0.912504912
##  [686] -0.782713704  0.209413866  1.358109845  0.285618724 -0.008227882
##  [691] -0.395262360 -0.002770746 -0.639193613  0.854900269  0.623296013
##  [696]  1.130491491  0.363513848  1.445845109  0.442221731  0.775364705
##  [701]  0.378344298 -0.033574284  1.200007075  0.668830391  0.235115259
##  [706]  0.332347751  1.185910928 -0.610457261  0.506119900  0.630502634
##  [711] -0.506753205  1.344339714 -0.506753205  0.893383066  0.597733542
##  [716]  0.764631001  0.543831954  1.140023639  0.824196192 -0.110651641
##  [721]  0.882543517  0.576315723 -0.231594095  1.214633634  0.359824886
##  [726]  0.856073544 -0.113471823 -0.789783136  0.835779255  0.461952458
##  [731] -0.491904099  0.483618505  1.104369392  1.177847187  0.479728546
##  [736]  0.413658563  1.205673930  0.454411587 -0.894638173  0.971225568
##  [741]  0.464371630  0.421968574  0.847988314  0.545143055  0.924725532
##  [746] -0.143836966 -0.317563421  0.749571252  0.455318267  0.716395213
##  [751]  0.507357325  1.157812087 -0.405930480  1.314848248  1.296846068
##  [756]  0.561865137  0.253633200  1.032686121  0.901548574  0.506119900
##  [761] -0.415637940  0.642039870 -0.293009317  0.818789094  0.894487196
##  [766]  0.673428946  1.308695004 -0.050468563 -1.302823365  0.971701756
##  [771]  0.526229660  0.995201961  0.201563240  0.437859174  0.421523618
##  [776]  0.396468186  1.047221109  0.089427854  0.582787974  0.396161876
##  [781]  1.116783906 -0.179720526  0.393821079  0.229783822  1.244071894
##  [786]  0.670561929  0.347097244  1.129823327  1.012736729  1.200573783
##  [791]  0.615484304 -0.046900604  0.828330385 -0.159752648  0.385441803
##  [796] -0.391200196  0.274009124  0.322059125  0.926002439 -0.230624059
##  [801]  0.993433324 -0.353524024 -1.100695333  0.599841165  0.777937911
##  [806] -0.363240014  0.603502963  1.229306080  1.206604051  0.790692957
##  [811]  0.475221641 -0.563891358  0.732155583  0.252791363  0.702120537
##  [816] -0.377800750  0.437393117  0.707713668  1.152635582 -0.434149578
##  [821] -0.533937838  0.247033920  0.519559645 -0.322524283  0.862350608
##  [826]  0.012859592  1.734795874  0.601177289  0.926554178  0.597210574
##  [831]  0.393975099  0.519279974  0.394163190  1.262245551 -1.489436065
##  [836]  0.495092788  0.073503853  0.601677988 -0.030914120 -0.391723036
##  [841]  0.351494747  1.150061076  0.474951801  0.695452853 -0.096733459
##  [846] -1.310548115 -0.573678477  0.217940624  0.636163958  0.950023878
##  [851]  0.604622049 -0.447267759 -0.265832756  1.222519103  0.276349795
##  [856]  0.070938095  0.194760284  0.765167382  0.235552812 -0.189518914
##  [861]  0.926338053  0.486466709  0.159348981  0.680821200  0.756352580
##  [866]  0.602595817  1.203776104  0.549764661  0.301257675  1.461092410
##  [871]  0.777328709  0.527146744  0.700954158 -0.267336784  1.366925751
##  [876]  0.536770942  0.631748445  1.434621115  0.596742708  0.921751485
##  [881]  1.263551531  0.535324155 -0.528243764 -0.584696959  0.978938113
##  [886]  1.194671303  0.816751653  0.454292566 -0.402437219  0.893144605
##  [891]  0.728734404  0.338418919  1.270319006  0.606026845  1.350610198
##  [896]  0.864772865  1.143407447  0.219823337  1.062220466  0.402105732
##  [901]  0.192379350  0.695460621  0.552813494 -0.601049618  0.229276425
##  [906]  0.806470773  0.540959675  0.498838110  1.173534420  0.204141943
##  [911]  0.146074048  0.196834847  1.157131982  0.225711429  0.161770218
##  [916]  0.866356666  0.269700454  0.875142705  0.595022553  0.145937138
##  [921]  0.138509269  0.735943631  0.065109738  0.739976745  1.735584137
##  [926] -0.331012820  0.459287227  0.441264682  0.213566497 -0.163557078
##  [931]  0.497025702  0.326735842  0.471691968 -0.495917600  1.025520206
##  [936]  0.245037384  0.697617920  0.889882333  0.830673280  2.015008829
##  [941]  0.198751387  0.242172193  1.572594091  0.538590085  0.571992128
##  [946]  0.581706846  0.666792178 -0.411846787  0.582305106  0.727718414
##  [951]  1.276460832  0.350166831  0.384424907 -0.280411097 -0.057085916
##  [956]  1.384530639  0.629631043  0.786033585  0.514253190  1.271790668
##  [961]  0.739589576  1.023237757  0.971943912  0.200883410  0.444786468
##  [966]  0.599431607  1.413530236  0.368179384  0.386514628  0.609390575
##  [971]  0.870093454  0.953754673  0.882281865  0.575819971  1.360242717
##  [976]  0.170016073  0.666228558 -0.439879133  0.619152820  0.039843109
##  [981] -0.308427436  0.956169530  0.834030471  0.236428395  0.383110575
##  [986]  0.878570447  0.876768597  0.483984608 -1.006608701  0.535483743
##  [991]  0.837877837  0.119761879  0.527306341  0.630456744  0.713715937
##  [996]  0.944137800  0.590259152  0.443441284  0.349565297  0.555725815
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
##      mean         sd    
##   0.6863600   0.4617867 
##  (0.1460298) (0.1032545)
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
## [1] -0.12323976 -0.06839665  0.30371350 -0.24975632 -0.09191863  0.25829765
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
## [1] 0.0018
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8943135
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
## t1*      4.5 -0.01121121   0.9516437
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 2 3 6 7 8 9 
## 1 1 3 2 1 1 1
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
## [1] 0.0192
```

```r
se.boot
```

```
## [1] 0.8994445
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

