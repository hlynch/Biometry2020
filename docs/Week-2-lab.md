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
## 0 2 4 5 6 9 
## 1 2 2 2 1 2
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
## [1] -0.0517
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
## [1] 2.684292
```

```r
UL.boot
```

```
## [1] 6.212308
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
##    [1] 4.5 4.5 4.5 5.9 4.4 4.1 3.4 3.6 4.6 4.5 4.2 4.5 4.9 4.2 3.4 6.8 4.1 3.7
##   [19] 2.9 4.5 5.3 3.6 5.3 4.3 5.1 6.3 4.6 4.3 4.7 5.3 3.4 5.9 4.5 3.6 5.4 4.3
##   [37] 4.4 5.7 5.4 4.6 4.0 4.6 5.2 4.2 3.8 4.8 4.5 4.4 4.7 5.5 4.4 5.8 3.4 5.6
##   [55] 4.8 4.8 3.5 5.7 4.5 3.5 2.9 4.4 3.9 5.1 4.6 4.9 5.0 4.0 6.7 3.9 4.7 4.8
##   [73] 4.2 4.8 4.4 5.4 5.1 5.9 5.7 3.7 5.1 4.7 6.9 3.7 3.4 5.0 4.1 3.3 3.6 3.9
##   [91] 4.3 6.0 4.3 4.6 2.5 4.6 4.4 3.5 5.5 4.4 6.2 5.0 4.0 4.2 4.5 5.6 6.1 5.3
##  [109] 3.4 4.9 3.5 5.4 3.5 2.3 4.2 3.3 3.5 4.9 4.7 4.0 5.1 3.9 5.5 5.4 6.4 4.5
##  [127] 4.3 4.6 4.8 5.6 4.0 4.6 5.5 4.1 5.1 2.7 4.8 3.7 4.0 6.9 4.2 4.0 4.1 3.9
##  [145] 2.9 4.6 3.3 4.6 5.1 4.4 5.1 2.1 5.6 5.3 3.5 4.3 4.0 5.2 4.7 6.5 5.5 5.4
##  [163] 5.2 4.2 4.1 4.8 3.1 5.1 3.5 3.5 5.4 5.3 3.9 4.2 5.8 5.6 3.6 5.7 4.1 4.1
##  [181] 4.9 4.7 5.8 4.3 2.9 1.3 3.7 5.3 5.6 5.2 6.0 2.7 5.6 5.5 5.7 2.6 4.0 5.0
##  [199] 5.1 3.3 4.1 3.4 6.5 3.7 2.9 4.0 5.6 3.7 6.0 3.7 2.8 3.6 5.7 4.8 5.0 4.5
##  [217] 4.5 4.1 5.5 6.0 3.3 4.1 3.6 3.6 4.9 4.3 4.5 5.3 5.1 4.8 4.0 4.3 3.1 3.7
##  [235] 4.6 3.7 4.0 5.9 3.9 6.5 3.8 5.0 4.4 3.9 4.5 3.4 4.7 3.9 4.5 3.7 5.6 5.2
##  [253] 5.2 5.3 3.9 3.7 4.5 4.6 5.6 4.9 3.2 4.9 6.0 5.2 3.5 5.2 3.8 5.7 4.8 6.1
##  [271] 3.5 4.7 5.1 5.0 3.5 4.4 4.1 4.3 4.9 6.2 4.3 6.2 5.0 5.2 4.0 5.4 5.3 3.0
##  [289] 4.9 5.0 3.2 3.7 5.3 4.2 5.0 5.2 3.8 5.1 5.2 4.9 2.5 5.6 5.8 4.8 5.5 5.5
##  [307] 5.0 3.2 3.7 6.1 5.5 5.4 3.3 5.1 4.8 4.6 4.9 4.1 4.8 4.3 2.5 4.3 4.1 5.0
##  [325] 4.2 5.6 3.3 4.6 3.1 5.4 5.4 4.0 3.9 5.3 4.9 3.9 5.2 5.3 4.8 5.7 3.6 3.4
##  [343] 5.4 4.9 3.4 3.7 4.7 3.2 4.7 4.1 4.6 4.1 4.6 4.7 5.5 2.1 5.4 2.7 5.0 4.0
##  [361] 5.0 4.8 4.8 5.2 4.7 3.8 4.0 4.8 5.6 3.8 2.8 3.8 3.3 4.7 3.7 4.1 3.7 4.1
##  [379] 4.4 3.4 4.9 6.1 5.6 5.0 3.1 4.8 5.0 4.5 3.6 5.0 3.0 5.8 3.6 5.0 4.6 3.6
##  [397] 5.3 4.2 4.7 3.0 5.5 3.9 6.0 3.8 3.8 3.5 3.6 6.0 5.6 3.6 3.7 5.3 4.9 3.5
##  [415] 4.9 4.1 5.6 5.1 4.4 7.6 3.5 2.4 3.9 3.3 5.2 4.0 4.9 5.0 4.4 4.3 4.8 5.3
##  [433] 4.9 5.5 4.3 5.2 4.4 4.2 6.3 4.7 3.3 5.6 4.8 3.0 4.3 3.6 5.4 4.8 4.4 4.6
##  [451] 5.1 5.1 5.4 4.5 5.2 4.0 4.3 3.0 3.8 4.0 4.4 4.6 5.4 4.4 3.6 3.3 5.1 4.2
##  [469] 5.5 3.1 3.2 2.8 6.3 4.3 3.5 6.1 4.8 3.0 3.9 5.1 3.2 5.1 4.3 2.6 5.2 4.2
##  [487] 4.1 4.7 4.3 3.8 3.9 3.1 3.8 5.6 4.5 4.4 4.4 4.8 3.9 3.7 4.2 3.2 4.3 3.8
##  [505] 4.0 4.5 4.3 4.5 4.8 4.4 2.2 4.0 3.1 6.3 4.1 2.8 4.2 3.3 4.4 3.5 4.5 4.8
##  [523] 7.0 4.6 4.8 4.6 4.6 5.2 4.6 4.6 5.2 5.5 4.7 5.2 4.6 3.7 3.1 5.1 3.9 3.3
##  [541] 4.8 4.9 5.0 3.9 4.6 5.4 4.4 3.0 3.3 4.8 4.1 4.3 3.8 4.3 4.1 4.7 5.6 3.9
##  [559] 3.7 7.1 6.0 4.6 5.3 3.9 3.7 3.1 5.7 4.9 4.2 5.1 4.0 4.6 7.3 3.8 3.1 4.3
##  [577] 3.6 6.9 4.4 4.0 4.3 3.0 3.7 4.4 6.3 3.7 4.0 3.9 3.6 5.6 3.5 3.2 5.6 5.3
##  [595] 4.3 3.4 3.6 4.8 4.2 3.2 4.3 5.5 5.2 3.5 5.4 3.4 5.2 3.8 4.3 4.1 4.3 4.5
##  [613] 4.7 6.7 4.5 4.7 4.2 3.4 4.8 4.8 4.4 5.2 2.6 4.0 5.9 3.9 4.7 4.3 5.3 3.9
##  [631] 3.7 4.1 3.4 3.6 5.6 3.5 4.9 5.8 3.3 4.2 4.2 4.8 5.0 4.8 4.7 5.3 3.9 3.4
##  [649] 4.6 4.6 3.0 3.8 6.6 4.3 3.7 4.0 4.7 1.9 4.2 5.0 4.3 4.0 5.9 4.8 4.1 3.4
##  [667] 2.8 5.0 4.0 3.0 5.2 5.3 4.0 4.9 4.2 5.2 4.1 5.5 3.0 5.2 6.5 4.6 5.0 4.6
##  [685] 5.2 3.7 5.3 5.6 5.0 5.4 3.2 4.2 5.2 4.8 4.4 4.4 4.8 3.6 5.6 4.5 4.6 5.9
##  [703] 5.2 3.4 4.6 3.8 4.8 4.2 2.8 6.1 3.9 4.1 4.2 4.0 3.5 4.2 5.8 2.9 5.3 4.4
##  [721] 4.1 4.1 5.5 4.7 3.5 4.3 4.3 3.9 5.5 5.7 5.1 5.6 5.3 5.3 5.7 3.9 5.5 4.0
##  [739] 4.3 4.1 4.5 2.7 2.7 3.9 5.6 5.3 3.5 5.8 4.6 6.3 3.2 4.8 6.7 4.2 3.9 6.7
##  [757] 4.9 3.0 4.2 4.6 4.5 5.3 5.6 2.9 5.8 3.6 6.1 4.3 4.4 4.4 4.9 5.0 2.6 4.8
##  [775] 4.8 4.8 4.0 4.9 4.8 4.9 5.6 4.3 5.0 4.2 4.7 5.0 4.1 3.8 4.2 4.8 4.8 2.8
##  [793] 4.9 5.3 3.4 4.5 5.7 5.5 4.4 5.9 4.6 4.4 3.7 4.4 4.2 3.7 4.9 5.0 4.9 5.4
##  [811] 5.7 3.6 5.5 4.8 3.5 4.0 4.2 4.2 4.4 3.6 3.6 4.9 5.9 5.8 3.9 4.8 4.6 5.1
##  [829] 2.7 3.3 4.8 5.2 4.8 3.7 2.7 4.0 5.1 4.4 3.7 3.7 5.4 3.1 3.8 4.7 5.0 4.6
##  [847] 3.1 4.9 2.9 4.9 5.7 5.1 4.7 2.9 4.3 6.8 4.8 3.5 4.4 5.2 5.8 4.3 3.9 3.6
##  [865] 6.1 4.1 5.2 4.6 5.5 3.4 4.9 5.8 4.9 3.7 5.4 5.0 3.5 3.7 4.4 2.9 2.7 4.7
##  [883] 4.1 3.6 5.8 4.1 4.3 5.4 4.2 6.1 4.1 3.4 3.3 4.8 3.8 3.6 4.6 4.8 5.5 5.5
##  [901] 3.2 4.8 5.6 4.8 3.3 3.1 4.0 5.7 3.4 6.1 3.8 1.8 3.5 3.8 4.6 4.8 4.2 3.0
##  [919] 4.4 2.8 4.4 3.7 3.4 6.8 5.6 3.0 5.1 4.5 4.2 5.4 4.0 4.6 5.4 3.7 4.8 4.5
##  [937] 4.2 3.6 3.9 4.3 3.4 3.3 6.0 6.1 3.9 4.2 5.2 3.9 5.8 5.4 4.9 4.5 4.1 4.0
##  [955] 4.6 5.5 3.6 6.1 3.1 4.1 3.3 3.9 5.3 4.5 5.6 3.7 3.7 3.9 1.7 3.6 3.2 2.5
##  [973] 5.4 4.1 4.1 5.0 4.9 3.9 4.9 4.0 5.2 4.6 4.4 4.5 4.6 4.6 3.2 4.1 3.4 4.7
##  [991] 5.1 5.2 4.2 4.1 5.3 4.3 5.9 4.1 5.7 4.9
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
## 2.7975 6.3000
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
##    [1] 4.6 5.1 3.6 4.8 3.2 4.1 3.8 3.9 6.8 6.1 4.0 5.5 5.8 3.3 2.9 4.7 2.9 5.3
##   [19] 3.9 5.6 3.9 4.0 5.6 3.9 4.2 5.0 3.0 5.2 5.5 4.8 4.9 5.6 2.9 5.1 4.9 3.6
##   [37] 4.0 5.4 3.9 4.5 4.9 3.5 5.0 4.9 5.9 4.2 5.8 3.7 5.0 3.8 4.2 5.3 4.6 2.2
##   [55] 5.1 5.8 5.9 5.2 2.5 4.2 5.9 3.3 4.5 4.8 5.9 4.6 3.9 5.3 3.9 4.6 6.0 6.6
##   [73] 4.2 4.0 6.2 3.0 5.4 4.6 5.3 5.6 3.3 5.1 3.1 5.0 4.5 4.2 5.6 5.2 4.2 4.7
##   [91] 3.1 2.9 5.2 3.5 5.2 5.3 4.6 3.5 3.1 4.4 4.3 6.6 3.3 4.6 4.6 4.2 3.9 3.3
##  [109] 4.1 3.0 4.6 4.9 2.6 4.4 3.3 4.7 4.2 5.3 5.6 3.2 6.0 4.9 4.2 3.7 5.2 4.7
##  [127] 5.5 6.1 4.6 3.2 3.8 4.0 4.3 4.0 4.8 4.6 3.8 4.6 5.4 4.4 4.4 5.1 4.0 3.6
##  [145] 4.4 5.6 4.7 2.8 2.1 4.9 3.8 4.9 3.5 5.1 3.9 4.6 3.7 4.1 5.0 6.3 5.2 3.9
##  [163] 4.8 3.5 4.5 5.5 3.5 5.6 4.7 3.1 2.7 2.8 5.4 4.1 4.5 4.2 5.2 3.7 5.5 5.1
##  [181] 4.5 5.2 3.9 5.3 5.1 3.9 4.7 4.9 4.6 5.7 6.0 3.2 5.8 5.0 4.4 5.4 3.0 5.3
##  [199] 2.7 5.8 4.2 4.3 4.6 4.2 5.4 6.4 4.7 3.9 3.6 4.5 5.0 5.5 4.3 4.3 5.2 2.7
##  [217] 4.0 4.5 4.2 3.8 4.3 4.4 5.0 4.7 5.1 5.3 3.9 4.4 5.4 5.0 4.4 4.9 4.4 2.2
##  [235] 5.2 4.0 4.6 3.0 4.9 4.8 3.1 4.3 4.3 4.2 4.6 4.2 3.3 5.5 5.3 4.4 4.5 3.7
##  [253] 3.8 3.7 3.6 4.3 6.4 4.5 4.6 3.1 4.7 5.0 4.1 5.2 5.1 4.9 6.2 5.1 4.8 4.9
##  [271] 5.3 3.0 4.3 6.6 3.6 3.3 5.3 5.1 3.8 4.3 5.1 4.4 4.8 3.3 5.3 5.7 4.2 4.6
##  [289] 4.6 2.2 6.0 3.9 3.7 4.0 3.9 3.5 4.4 4.0 4.7 6.2 1.8 4.2 3.6 4.6 3.7 3.8
##  [307] 4.1 3.4 5.0 5.0 3.9 3.9 4.8 4.2 4.1 3.8 4.4 5.7 5.0 4.0 5.9 5.1 4.5 5.9
##  [325] 4.6 4.3 4.7 3.7 2.4 4.4 3.4 5.4 4.5 4.2 5.0 5.1 3.1 4.0 4.0 4.5 5.0 4.2
##  [343] 3.7 3.8 6.5 4.3 4.4 3.4 4.3 3.1 4.2 5.2 4.7 5.1 5.1 5.2 4.0 5.1 4.7 5.0
##  [361] 3.7 4.5 4.7 3.9 4.2 3.7 6.3 3.2 5.7 6.1 6.0 2.7 4.2 4.1 4.2 5.1 4.1 4.9
##  [379] 4.0 4.7 4.6 5.1 5.2 4.9 5.6 4.3 5.6 4.2 4.3 4.6 4.2 5.9 5.3 4.4 4.8 4.8
##  [397] 5.2 4.0 5.6 4.4 3.7 3.4 4.2 3.6 5.0 4.1 3.5 5.5 4.1 5.3 3.0 5.2 2.5 3.0
##  [415] 4.0 3.6 4.3 5.0 4.5 4.8 6.5 4.2 5.0 4.1 5.5 4.2 4.5 5.4 5.7 5.2 4.7 5.3
##  [433] 3.1 5.0 5.0 3.4 3.5 4.0 4.6 5.4 4.7 4.9 5.8 4.8 6.2 5.8 5.9 4.7 6.1 5.1
##  [451] 4.7 3.4 4.6 5.2 4.5 2.9 4.5 5.2 2.4 4.6 4.4 5.1 3.5 4.1 3.3 6.4 5.9 5.4
##  [469] 3.9 4.8 5.0 3.5 4.6 3.9 5.3 5.6 4.0 3.2 4.5 3.6 4.9 3.6 5.2 4.8 4.1 4.8
##  [487] 3.4 4.2 4.4 4.5 4.4 3.6 3.4 3.7 4.5 3.9 2.7 3.6 5.1 4.3 3.9 4.1 5.8 5.1
##  [505] 5.7 4.8 4.6 3.3 3.9 4.0 4.0 3.0 5.0 6.3 4.1 4.3 3.8 6.3 4.9 4.9 4.1 4.9
##  [523] 4.7 5.0 4.7 4.0 4.3 4.7 4.6 2.8 5.2 4.5 4.3 4.9 3.9 4.6 5.2 4.9 4.4 4.3
##  [541] 3.9 4.4 4.4 5.9 4.3 4.0 4.7 3.1 3.9 5.3 4.2 5.9 5.0 5.4 5.1 6.1 4.4 5.8
##  [559] 5.4 4.4 4.9 4.4 3.7 4.3 5.6 4.0 5.1 3.7 4.8 4.1 5.6 3.6 4.1 5.8 3.5 3.4
##  [577] 5.2 4.1 4.9 2.8 5.0 3.2 4.6 5.1 3.3 3.5 4.2 5.8 4.3 5.6 3.6 4.2 3.9 4.9
##  [595] 4.2 5.9 4.8 4.0 6.8 5.2 4.6 5.2 5.8 4.7 4.7 3.7 5.7 3.2 6.8 4.1 4.3 5.3
##  [613] 3.8 5.8 4.7 4.8 4.2 5.9 4.3 3.9 4.8 4.5 3.5 4.0 2.6 4.2 4.6 4.9 4.6 5.5
##  [631] 3.5 4.9 4.4 1.9 4.6 4.3 4.9 5.8 4.1 2.5 5.1 3.8 3.5 5.3 5.0 5.4 5.0 5.6
##  [649] 4.0 4.1 3.2 4.0 3.0 6.2 4.7 4.6 4.7 5.2 3.8 4.2 3.3 3.2 3.1 5.1 5.3 4.3
##  [667] 4.2 5.5 4.6 4.5 4.1 3.4 6.0 2.1 4.7 5.0 4.5 4.8 5.2 3.4 3.1 4.4 2.8 6.2
##  [685] 5.2 3.0 2.2 5.5 5.4 4.1 4.1 4.6 3.8 5.9 5.4 5.7 4.6 3.3 4.1 5.6 4.2 4.5
##  [703] 5.1 2.5 3.7 5.2 4.7 5.0 4.0 4.2 4.8 2.4 3.8 6.5 4.5 4.2 4.9 4.8 5.2 5.7
##  [721] 4.8 3.7 5.9 4.3 3.3 3.8 6.1 5.1 6.1 4.9 4.1 4.1 4.6 4.0 4.6 4.4 5.0 4.6
##  [739] 4.1 3.7 3.3 5.1 4.9 4.8 6.5 5.1 3.9 4.3 3.7 3.8 4.6 4.5 6.1 5.8 3.4 5.3
##  [757] 5.1 5.2 2.5 4.9 3.7 5.0 4.5 5.5 5.3 3.7 4.9 4.5 5.3 4.8 6.6 4.3 4.9 3.5
##  [775] 5.3 4.9 3.9 4.0 5.2 5.0 4.2 4.7 5.6 5.5 4.7 5.6 4.3 6.5 4.1 5.3 4.0 4.1
##  [793] 5.9 4.3 3.9 4.7 2.4 4.2 5.0 4.3 3.9 3.6 3.3 4.2 4.0 3.6 4.0 4.4 5.2 4.7
##  [811] 4.9 4.1 2.8 5.3 5.2 5.5 4.2 3.2 4.0 4.0 3.4 6.3 5.2 5.0 5.0 4.4 4.7 4.6
##  [829] 5.7 4.0 5.0 3.2 3.4 5.3 4.5 4.5 4.9 5.7 4.4 4.2 4.8 5.3 4.2 4.8 4.1 4.3
##  [847] 4.4 3.7 6.1 4.5 4.1 2.9 3.6 5.3 3.0 4.2 4.6 5.2 4.7 5.1 2.7 3.4 4.4 5.2
##  [865] 3.0 4.1 4.7 6.0 2.7 4.2 4.9 2.3 3.9 3.6 4.2 4.5 5.1 3.1 3.8 3.7 4.5 4.0
##  [883] 5.0 3.5 3.2 5.1 5.0 5.2 6.4 4.9 4.0 5.9 4.5 3.7 3.1 4.2 3.8 6.5 3.6 4.2
##  [901] 4.9 4.4 3.4 1.7 2.8 3.3 4.6 5.5 5.7 4.1 6.9 4.4 4.1 5.0 5.3 5.3 5.7 4.3
##  [919] 4.8 4.4 3.8 4.8 5.0 2.4 4.0 3.7 5.0 3.7 4.9 2.8 4.0 4.8 2.9 4.1 6.7 5.8
##  [937] 5.3 4.2 4.1 4.5 3.6 5.6 6.4 3.8 3.6 2.7 5.1 4.4 2.9 5.6 5.4 4.7 5.5 6.3
##  [955] 4.0 4.4 5.2 3.9 5.1 4.8 4.2 3.5 3.4 4.0 3.2 3.8 1.9 5.8 4.5 4.5 5.6 2.9
##  [973] 4.3 4.5 4.4 2.9 5.5 6.3 3.8 4.0 5.0 3.2 5.8 5.8 4.5 5.5 3.7 4.7 3.8 4.4
##  [991] 4.0 4.4 5.5 4.5 5.3 5.1 3.8 4.8 5.9 3.3
## 
## $func.thetastar
## [1] -0.0105
## 
## $jack.boot.val
##  [1]  0.49751381  0.29740634  0.30550459  0.19350282 -0.02234957 -0.09106145
##  [7] -0.19116022 -0.30792350 -0.38862974 -0.53584337
## 
## $jack.boot.se
## [1] 0.9670764
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
##    [1] 4.6 3.3 4.3 4.2 5.4 3.2 4.3 4.1 4.0 6.1 5.0 5.5 5.7 4.5 5.1 3.3 5.3 5.6
##   [19] 5.2 3.4 3.6 3.1 5.5 4.9 5.0 4.0 5.1 4.7 4.4 5.4 4.3 3.6 6.8 4.3 4.6 2.9
##   [37] 2.6 4.8 3.7 4.4 4.7 5.6 3.6 5.0 3.9 3.9 5.3 3.2 2.6 5.2 3.1 4.6 5.9 5.1
##   [55] 4.6 5.3 4.7 4.8 5.2 4.8 4.2 4.7 3.0 3.5 4.3 4.6 4.5 5.3 5.3 3.8 4.2 3.9
##   [73] 4.5 4.4 3.4 2.4 5.1 5.1 4.2 4.2 4.3 3.8 3.6 3.4 4.2 4.7 2.9 3.5 3.8 4.3
##   [91] 6.1 4.1 3.0 5.3 2.7 3.9 4.0 4.8 3.0 2.7 3.7 5.7 4.2 4.1 3.9 3.7 5.1 5.9
##  [109] 4.6 4.4 6.2 3.2 5.4 4.3 4.0 5.5 5.9 5.4 4.1 6.1 5.0 4.0 4.6 5.1 5.0 4.0
##  [127] 4.8 5.0 3.7 4.3 5.6 4.9 4.7 5.2 3.7 5.2 2.5 5.1 5.3 4.9 3.6 3.0 4.4 4.2
##  [145] 3.0 4.1 4.2 5.7 4.0 4.4 3.3 5.2 4.5 4.0 4.5 2.3 4.0 4.5 5.5 5.0 3.1 3.4
##  [163] 4.0 5.1 4.7 4.8 4.9 5.5 4.4 4.0 2.8 6.0 5.8 4.5 4.3 4.8 4.3 3.8 5.1 4.5
##  [181] 4.9 4.3 4.5 3.2 4.5 2.6 3.8 3.8 4.6 4.1 5.6 6.1 4.5 4.6 5.0 4.8 4.6 3.0
##  [199] 3.8 4.3 5.5 3.8 4.6 5.0 2.8 4.2 3.5 6.5 6.4 3.7 5.0 4.7 5.7 2.9 5.8 4.3
##  [217] 5.5 5.1 4.6 5.3 5.1 5.3 4.4 5.8 4.5 3.7 4.3 5.6 3.7 4.9 6.0 4.8 3.0 4.8
##  [235] 3.6 5.1 2.8 3.5 2.9 3.5 4.4 4.6 5.5 4.0 4.3 5.3 4.1 4.0 5.9 4.8 5.7 5.6
##  [253] 4.6 5.7 3.9 4.4 4.6 4.7 6.6 3.7 4.8 3.1 5.4 3.2 5.1 6.2 5.5 4.6 5.7 5.5
##  [271] 5.0 5.2 3.4 4.4 5.8 3.9 3.0 5.2 5.5 4.6 5.2 3.9 5.9 5.3 6.0 5.2 4.3 4.1
##  [289] 2.3 5.3 4.1 3.0 5.4 3.9 5.3 3.9 4.1 4.7 4.7 4.3 5.8 4.2 6.6 4.6 4.5 4.6
##  [307] 3.6 5.9 4.2 3.1 4.0 3.8 4.6 3.8 4.7 4.9 3.7 4.5 2.8 6.1 5.4 4.2 4.0 4.9
##  [325] 3.3 5.6 4.1 2.7 5.0 5.2 4.5 5.3 3.9 5.3 4.6 5.5 4.9 4.3 3.8 6.1 5.0 4.6
##  [343] 4.0 2.8 5.7 4.0 4.6 4.0 5.5 5.2 4.0 4.1 4.5 4.7 4.3 4.6 2.8 4.3 2.9 4.0
##  [361] 4.3 3.2 3.7 4.2 4.8 2.9 4.6 5.3 6.3 5.8 5.1 3.2 4.1 5.5 4.8 4.6 4.0 5.0
##  [379] 5.9 5.1 3.9 4.5 3.3 4.5 5.6 5.9 4.3 4.2 4.0 4.4 4.8 4.1 3.2 4.2 6.1 5.3
##  [397] 4.4 6.5 2.8 3.7 5.3 4.3 5.3 4.3 2.2 6.0 6.5 3.8 3.1 5.7 4.7 4.5 5.1 4.0
##  [415] 3.0 5.6 4.1 5.0 4.5 4.6 4.1 4.5 4.0 4.0 5.0 4.1 4.0 4.7 6.0 5.3 4.8 4.6
##  [433] 3.9 5.0 3.7 3.3 3.8 5.0 5.3 4.8 5.4 3.5 5.0 4.8 3.6 2.7 4.3 5.5 6.5 4.3
##  [451] 5.7 4.6 3.6 3.8 3.1 5.3 5.9 4.8 5.8 3.7 2.3 5.6 3.3 4.1 3.5 4.6 4.1 5.7
##  [469] 4.4 4.4 3.8 5.8 3.2 5.9 4.8 5.1 2.6 5.0 3.8 4.1 4.3 3.0 5.3 5.2 5.5 4.4
##  [487] 3.2 4.1 4.3 5.9 4.2 6.7 3.9 5.1 3.7 2.5 4.4 3.9 3.1 4.2 3.2 3.2 4.4 3.9
##  [505] 5.8 3.6 4.2 3.8 4.6 3.7 3.7 4.1 3.9 3.1 5.0 3.5 3.3 3.2 4.8 3.8 3.8 3.7
##  [523] 4.8 5.6 4.3 4.2 5.2 5.4 2.0 6.7 2.3 4.5 4.1 4.3 5.1 5.7 5.6 4.9 4.8 4.5
##  [541] 4.9 3.5 5.4 4.0 3.1 3.6 6.4 4.9 3.8 5.7 5.2 4.8 5.6 3.5 3.9 3.9 3.3 4.2
##  [559] 5.0 5.1 5.5 4.9 4.8 3.4 4.7 6.5 5.1 3.5 5.3 4.0 5.8 4.4 4.3 2.9 3.5 5.7
##  [577] 4.7 5.6 4.2 5.4 5.0 4.1 3.5 4.9 4.1 2.2 3.4 5.3 1.3 3.8 4.5 3.1 4.4 5.0
##  [595] 4.4 2.5 5.2 3.5 4.5 4.0 6.6 2.6 4.3 3.7 4.2 4.2 4.4 5.5 6.5 6.9 3.9 5.2
##  [613] 7.1 4.1 6.4 5.7 5.8 4.3 4.7 4.2 4.1 4.9 6.4 5.1 2.9 5.1 6.3 3.6 5.3 4.4
##  [631] 5.3 5.5 5.2 4.9 6.1 3.6 3.7 3.7 4.0 4.6 3.6 4.5 2.5 5.8 5.2 5.4 5.7 4.2
##  [649] 4.8 3.0 3.8 5.9 4.3 5.3 5.1 4.5 4.6 4.9 6.4 5.1 5.8 3.9 5.0 4.5 5.4 3.2
##  [667] 4.3 4.0 4.4 5.5 3.5 4.6 4.4 3.6 4.2 5.5 3.7 5.4 3.2 5.2 3.1 3.3 5.3 4.9
##  [685] 4.4 5.7 4.5 5.3 4.6 4.1 4.7 5.1 5.8 4.2 5.3 3.0 3.9 4.7 4.7 4.8 5.1 3.6
##  [703] 5.4 4.1 4.3 4.2 4.9 4.1 4.7 4.1 4.6 4.8 4.4 3.0 3.8 5.5 3.0 5.8 4.1 4.8
##  [721] 4.9 3.9 3.3 4.0 5.5 4.4 5.7 4.1 3.8 5.2 5.3 4.5 5.4 4.7 5.6 4.8 3.2 4.5
##  [739] 5.6 5.3 2.9 5.0 5.0 4.7 6.2 2.9 4.7 3.8 3.9 3.6 4.7 4.0 2.9 2.8 3.8 3.7
##  [757] 4.8 4.5 3.9 5.2 4.5 5.2 4.4 4.5 4.7 4.9 5.2 5.2 5.3 4.4 5.2 5.2 4.5 3.4
##  [775] 3.8 4.5 5.0 3.4 3.7 4.0 4.4 5.3 5.2 3.4 5.9 3.1 5.4 5.4 6.2 3.0 4.4 5.5
##  [793] 3.6 5.7 3.9 4.1 4.9 6.2 3.1 4.5 5.1 4.5 4.2 5.3 4.5 3.8 5.3 3.7 3.7 4.5
##  [811] 5.3 3.0 3.7 6.4 4.2 5.1 4.5 4.4 5.6 3.8 5.5 3.5 5.8 5.6 4.9 5.0 5.1 5.4
##  [829] 4.4 5.2 3.9 3.2 4.7 4.3 4.6 5.1 4.4 4.2 5.3 4.1 4.2 4.6 4.1 5.9 5.1 5.1
##  [847] 4.9 4.8 2.1 3.9 5.2 4.3 4.4 5.5 4.8 5.6 5.8 5.8 3.8 4.2 5.3 3.9 5.8 4.0
##  [865] 5.4 5.2 2.9 2.9 5.6 4.1 5.1 4.9 5.6 3.3 4.3 4.5 5.1 5.9 6.4 2.4 4.8 3.8
##  [883] 3.9 4.3 4.8 4.7 5.0 5.6 6.0 4.7 4.1 4.6 3.6 4.1 3.3 4.5 5.8 3.6 6.0 4.0
##  [901] 4.8 5.2 5.7 5.1 5.2 4.1 5.2 4.6 2.9 5.2 3.4 4.6 3.9 2.8 4.9 5.2 5.2 3.2
##  [919] 3.1 3.8 4.3 4.7 3.9 4.5 4.0 4.8 3.1 4.5 5.3 4.5 4.7 6.2 3.6 5.1 3.9 4.5
##  [937] 3.5 5.2 4.1 5.9 4.0 3.7 3.8 4.4 3.4 3.0 3.3 4.7 4.0 6.6 5.4 5.3 3.9 5.9
##  [955] 4.5 3.4 3.2 4.5 4.8 5.2 4.9 6.0 4.2 4.3 5.0 3.7 2.7 5.1 2.4 4.4 5.1 4.5
##  [973] 6.2 3.9 3.3 5.9 5.0 5.0 4.5 4.1 3.8 4.4 4.6 4.9 4.8 5.1 3.7 5.8 3.5 4.7
##  [991] 3.8 4.3 5.7 5.0 2.8 4.7 3.8 5.3 4.0 5.4
## 
## $func.thetastar
## 72% 
## 5.1 
## 
## $jack.boot.val
##  [1] 5.500 5.388 5.300 5.300 5.200 5.100 4.900 4.800 4.600 4.500
## 
## $jack.boot.se
## [1] 0.9748706
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
## [1] 0.2879141
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
##   4.764749   8.174204 
##  (2.060596) (3.728146)
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
## [1] 0.2523278 0.8468052 0.3980653 1.5697386 0.2026111 0.9616714
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
##    [1]  0.0204616826 -0.0951355836  0.6093184470 -0.0902277636  0.8992414777
##    [6]  0.1700186940  0.2234458540  0.6832185966  0.0125023555  0.0415735684
##   [11]  0.1016936046  0.0536620390 -0.1290922731  0.4726391388  0.6432034652
##   [16]  0.0992029460  0.7875824918  0.3021228236 -0.0344558903  0.8591316515
##   [21]  0.2536596072 -0.2239487912  0.1933694261  0.4842886682 -0.2794604684
##   [26]  0.3462737546 -0.2072768768 -0.5086520834 -0.4010170678 -0.0065885453
##   [31] -0.2892012439  0.5857222806 -0.4162310487 -0.0169978523  0.0738584756
##   [36] -0.2701582233  0.3982100754  0.4816242842  0.7055093967  0.6200367019
##   [41]  0.2178065415 -0.3206378075  0.0637469510  0.6671585122  0.1445760234
##   [46]  0.2526799466  0.8945988910  0.0119791781 -0.1126314775  0.1833929159
##   [51]  0.8512601715  0.0859120803  0.5091530058  0.1455638329  0.5039255136
##   [56]  0.1773090690  0.4402416267  0.4172765593  0.1828690357  0.5301160306
##   [61]  0.0151668681  0.7017422400  0.6135271545  0.2100281078  0.2131295160
##   [66]  0.0400718285  0.1432947372  0.3862874819 -0.3557341688  1.2554245675
##   [71] -0.6424260830  0.2124378541  0.1281004225 -0.1344412156 -0.0331475023
##   [76]  0.7196200873 -0.5182455356  0.3263813551  0.7364242361 -0.0453071553
##   [81]  0.5595160929  0.0375517186 -0.3189509584  0.6069874009  0.4374707709
##   [86]  0.7131840915  0.8603368295 -0.0607337582 -0.2310410737  0.3020185979
##   [91]  0.1032920839 -0.0378603939 -0.2677846785  0.5838547018  0.6315263932
##   [96] -0.0795101560 -0.0591476486  0.1051385402 -0.0577621270  1.2351037963
##  [101]  0.7356449183  0.8432778552  0.2750471571 -0.0693315243 -0.3138455267
##  [106]  0.0448260261 -0.1210057140  0.9470997848 -0.1235846602  0.3024775088
##  [111] -0.6221900837  0.7131288677  0.2522985510  0.1851286553  0.1138620575
##  [116] -0.3318496322 -0.4911865909 -0.2111643767  0.4418136891  0.6222775744
##  [121]  0.2587078660 -0.0421272238  0.7345562710 -0.6400637662  0.2760795102
##  [126]  0.2247747046 -0.6975050857  0.7616149455  0.3366947617 -0.1658526199
##  [131]  0.9536863090 -0.0069772477  0.3456652298 -0.0660510513  0.4621398195
##  [136] -0.4440669793 -0.2206987532 -0.2446288681  1.1832034002  1.0054350167
##  [141] -0.3481628128 -0.0592810278  0.1317499231  0.3637289807  0.5997638276
##  [146] -0.1527116387  0.3479258434  0.0851981420  1.2539130470  0.2198029370
##  [151]  0.2337359741  0.5209826178  0.2694539233  1.0706698479  0.3713176273
##  [156]  0.2269794274  0.1625013364 -0.0717251209  0.7436344416  0.1202441684
##  [161]  0.3437567053 -0.3995363007  0.5907165145 -0.5914283700  0.1241604767
##  [166] -0.0150487860  0.4340283014  1.0172897325 -0.0403144451  0.1372081986
##  [171]  0.5412599911 -0.1032032552  0.5749454733 -0.2425259544  0.0355101507
##  [176]  0.2947995452 -0.3046147620  1.3973907491  0.1653036263  0.4132362067
##  [181]  0.3630670173 -0.3522982672  1.2914069646  0.6177478132  0.5884830747
##  [186]  0.1536400127 -0.4058202883 -0.1612485226  0.4579771240 -0.0636969028
##  [191] -0.6349802920 -0.6047551881 -0.9526365692  0.7544637469  0.6252387197
##  [196]  0.2217016058  0.2502807904  0.4949306563  0.0594972876 -0.2341641205
##  [201] -0.5173613640  0.1343401118  0.2446170116  0.0768304376  0.3158321917
##  [206]  0.3517154685  0.5630613129  0.5978075247  0.0275668674  0.6148186585
##  [211]  0.7031752106  0.8088918755 -0.3763321218  0.1190146660 -0.5781936737
##  [216] -0.2140568092  0.1499526737 -0.4145727935  0.2681844710 -0.2462998426
##  [221] -0.6159348534 -0.4377658900  0.4633678215  1.0694634717  0.7183029605
##  [226] -0.0568771926 -0.0376260756  1.0383144844  0.4355803116  0.6623600956
##  [231] -0.5926264500  0.6520251611 -0.3691104954  0.1210439628 -0.0401614115
##  [236]  0.6509372108  0.2691292801  0.7031752106  0.5181729229  1.3948123240
##  [241]  0.5786407389 -0.0258658602 -0.4787521172  0.7717254053  0.1191719857
##  [246] -0.4992480315  0.3140760460  0.1301110341  0.4170655106  0.3891932637
##  [251] -0.2047174409  0.5413095751  0.2184540438  0.6910934059  0.6692367475
##  [256]  0.5679746643 -0.4779669550  1.1194119762  0.1809772707  0.3298047506
##  [261] -0.4165728645  0.6026037589  0.7585057810  0.4058955651  0.5164549681
##  [266]  0.4384591356  0.1207249430  0.4306446048 -0.3544540584  0.4407557109
##  [271] -0.1796105137  0.7156664702 -0.5065304249 -0.2672933059 -0.3360987214
##  [276]  0.7063996330  0.9904649779  0.9550666609  0.0658302507 -0.3540015101
##  [281]  0.8393357884 -0.0623316010  0.2303566051  0.3145175572  0.4826282842
##  [286] -0.0920273926  0.0640455208 -1.0976310051  0.8432048472  0.3487407967
##  [291]  0.3996861773  0.4200363804  0.4909989467 -0.1088359073  0.3526853672
##  [296] -0.6948774706  0.3260693685  0.5306920790 -0.1112234534 -0.1618886250
##  [301] -0.4504881608 -0.7720210688  0.3419072837  0.4400625936 -0.0689165472
##  [306]  1.2213520373 -0.1902376188  0.6313921927  0.3512258172  0.0910845344
##  [311]  0.5417071934 -0.5088900861 -0.0912010489  0.9575511522  0.5952289004
##  [316]  0.3862388217 -0.6261565705 -0.7008161587  0.0624006861 -0.2525233655
##  [321] -0.0874712495  0.5979757999  0.7249435654 -0.6742353594  0.2815165359
##  [326]  0.4352152616  0.1012996763  0.6013593097  0.1682812343  0.6520251611
##  [331] -0.3672497520  0.5379230789 -0.0094619264  0.1355619769 -0.6574990277
##  [336] -0.3219900039  0.0205384536 -0.2184984970  0.0429791667  0.7922866750
##  [341]  0.6550145128  0.2470621875  0.4105259944  0.5408031293 -0.1783507025
##  [346] -0.5113775304  0.4314950013  0.0800672152  0.3303664244  0.8839643295
##  [351]  1.0877616932 -0.4307336512  0.4603964759  0.1760682221  0.0549931404
##  [356]  0.1847850262 -0.2362177073  0.2046629183  0.0777289601 -0.5342154118
##  [361]  0.4114186471  0.2835923933  0.1091029063 -0.5447327837  0.0537040919
##  [366] -0.1379119007 -0.7340435355  0.3858158115  0.0607447848  0.1846491585
##  [371]  0.6000128808  1.2971000732 -0.6850637620 -0.3851888070 -0.0116678540
##  [376]  0.0694201600  0.4553755570  0.2445997266  0.5404340865 -0.2485827012
##  [381]  1.1581312614  0.0205373373  0.1169560660  1.4280800730  0.4684478354
##  [386] -0.4416180958  0.1156832562 -0.5716249278 -0.5785247632 -0.6621440038
##  [391]  0.6119600818 -1.2153554446 -0.2189751889  0.4481094217  0.9344651268
##  [396] -0.6865824916  0.4725873088  0.1548452691  0.5585842816 -0.0231009870
##  [401] -0.1194456928 -0.0222384744 -0.0855643259  1.3445429952  0.8825948135
##  [406]  0.7556290865  0.7110562676 -0.3481057334 -0.1449320410  0.0903726956
##  [411] -0.0222162161  0.2302405586  0.0024028893 -0.0856948948  1.0248130814
##  [416] -0.2645640384 -0.1336485166 -0.0838582953  0.0700950468  0.1050993707
##  [421]  0.2740447764 -0.4029201620  0.0333413757 -0.8157655021 -0.1718059042
##  [426]  1.0390158764  0.7063996330  0.4989874117  0.4540556421 -0.0762058740
##  [431]  0.0841229940 -0.2603062815  0.9152354663  0.2122199195  0.2256917472
##  [436]  0.3241025712  0.1751545028  0.5210675598 -0.2511245025 -0.1255878017
##  [441] -0.0114871941 -0.2013139449  0.7002631170  1.2011153625  0.9177945854
##  [446]  0.1073145721  0.5370296178  0.1726635100  0.3752243618 -0.0260036312
##  [451]  1.0085072117 -0.3155110632  0.2446170116  0.5314811304  1.1616905254
##  [456] -0.6879974623  0.7280684278 -0.1449320410  0.3636984489  1.0364800768
##  [461]  0.3206978975  0.2634207875  0.3369236613 -0.1111286838  0.3462190139
##  [466]  0.1279613757  0.9427738259 -0.4158737656  0.3727163195  2.1198214842
##  [471]  0.1453791147  0.4742260751 -0.2465176627  0.9923223528 -0.5904480254
##  [476] -0.6324172033  0.7706082857  0.6484040159  0.4707035594 -0.2050614185
##  [481]  0.7831917920  0.1280315080 -0.1161870196 -0.2396808231  0.2261776031
##  [486] -0.2933112951  0.2383084008 -0.0002051643  0.3176935272 -0.2651428016
##  [491]  0.4123344243  0.2315128344  0.6251426626  0.6789548651  0.8761662647
##  [496]  0.7113399165  0.0988094832  0.9300879120  0.7835926517  0.5380178980
##  [501]  0.7610576373 -0.0234710536  0.6958278378  0.3037192652 -0.2765016596
##  [506]  0.0697493445  1.4148682062  0.2224217045  0.4953254808  1.2574763714
##  [511]  0.2491986598 -0.0300271738  0.6360973009  0.7597535438  0.2884961026
##  [516] -0.0965537326 -0.0443805541 -0.2482127099  0.6095026861 -0.2558491754
##  [521]  0.1954178074  0.7633232405 -0.1122192398  0.0954298871 -0.8447677189
##  [526]  0.1017174226  0.2757365000  0.2314882507  0.0979311761  0.2232115446
##  [531] -0.0678904860  0.1736673577  0.1498238665  0.7668979547  0.0057032913
##  [536] -0.1359512139  0.1537461031 -0.3031924978  0.2827258855 -0.2707483816
##  [541] -0.2351038975  0.7776778283 -0.2506102165 -0.0315772537  0.5027877251
##  [546]  0.6681759327 -0.6933694029  0.6541372322  0.3450073877 -0.3362064566
##  [551]  0.2854635793  0.0525273776  0.1377540111 -0.2041744344 -0.7553510416
##  [556]  0.2436137626 -0.1726693599  0.4852982497 -0.1734356837  0.5104113646
##  [561]  0.0142983625  1.3420592133  0.5454807480  0.9609095327  0.1623637081
##  [566]  0.1139429951 -0.0082175539  0.1279694144  0.0005942766 -0.5802811922
##  [571] -0.2535435927  0.1216984606 -0.0512156077 -0.1110832920  0.7417665086
##  [576]  0.0628880160  0.7104475876 -0.0441751982 -0.0040813066 -0.0010801604
##  [581] -0.8565778211  0.4315337969 -0.1283635436 -0.3216821245 -0.5461731847
##  [586]  1.0143925770 -0.0107447193  0.5673458354 -0.6152495246  0.8533876473
##  [591]  0.3208486901  0.2122026343 -0.2072974212  0.4605271491  0.3146109458
##  [596]  0.1343907641  0.5724778745  0.3822901860  0.6887462536 -0.0426055259
##  [601]  0.2156403430  0.5122755499 -0.2455530822  0.0768304376  0.7400955358
##  [606] -0.2451801773 -0.5184664379 -0.2303172289  0.6573340257  0.9765066244
##  [611]  0.6480685835  0.0637501982 -0.5157697450  0.0955173573 -0.4162091290
##  [616] -0.2760877151  0.5274544448 -0.0115300205 -0.1763349654  0.1915909847
##  [621] -0.0585590284  0.6817128167  0.4871100240 -0.8613116257  0.5740591359
##  [626] -0.6458796381  0.1212178495 -0.6069658792  0.5036155535  1.1843663821
##  [631] -0.1498797095  0.3642178203  0.0052911208  0.3954383938 -0.0091063862
##  [636]  0.3270758700 -0.6942337449  0.9064437629  0.2827560142 -0.1116938311
##  [641]  0.6897880903  0.4056887170  0.5886890347 -0.5151409859  0.6188437145
##  [646]  0.1558775712  0.2710876086  0.3004405028 -0.1468995711  0.5405071086
##  [651]  0.5599410022 -0.2351038975  0.4204970909 -0.2580802114  0.5287653657
##  [656] -0.4569742855 -0.8600460829  0.1129184440  0.0782362957 -0.3373527580
##  [661] -0.1530965476 -0.2479425231  0.6207647456  0.4993975351 -0.3514605640
##  [666]  0.9711760450 -0.0569236637  0.7841120088  0.1359459947 -1.3577211070
##  [671]  0.3460477859  0.2434276612  0.5339727064  0.2461391662  0.1740998502
##  [676]  0.9748395049 -0.6940263635 -0.2796269872  0.3457337359 -0.1348233024
##  [681]  0.2937974209  0.2192944664  1.1458675341 -0.2543772246 -0.4802660213
##  [686]  0.5297972913  0.3466076362 -0.4273845984  0.0570853482  0.1348986322
##  [691]  0.1489612898  0.6112001533  0.0050737439  0.6645990698  0.7411629415
##  [696]  0.3369963392 -0.2376477261  0.2870952402 -0.1053156985 -0.0055570326
##  [701] -0.1056378467  0.8226724384 -0.5780816886  0.3122197249  0.6042581363
##  [706] -0.1235342952 -0.1103685535  0.1338235416 -0.2258545179  0.5900743280
##  [711]  0.3295963871  0.4797070842  0.7216874176  0.6704967644  0.2183979108
##  [716]  0.6660524495  0.7766100862 -0.1457646151  0.4603324962  0.0114170854
##  [721]  0.6079838182 -0.3446850873  0.2553175752  0.1308922492  0.2878683664
##  [726] -0.0258256102  0.3741901193 -0.1690105495  0.3567145013 -0.2703240284
##  [731] -0.2795185550  0.5066185823 -0.0296830964  0.9111487432  0.7604152878
##  [736]  0.6834514582  0.4077980850 -0.0060324599 -0.2153274782  0.8763019834
##  [741]  0.1616805513 -0.5415999327  0.9154703962 -0.3077723463  0.7633083940
##  [746] -0.0167126840  0.2684785383  0.8481402185  0.4939410076  0.2636190445
##  [751] -0.5182714639  1.5903369165  0.0330666128  0.0782714136  0.3513156427
##  [756] -0.1584596114  0.3372801688  0.3155503053 -0.4469901123 -0.1272108272
##  [761]  0.9141261943  0.2359649686  0.3759082208  0.3214574323  0.1347555334
##  [766]  0.1187436718  0.5649422600 -0.5783429667 -0.2610588530  0.7557679142
##  [771] -0.3078126094  0.6871698840 -0.1094701260  0.6436796761  0.0763747565
##  [776]  0.7169562425  0.0844559670  0.1385605029  0.1307390743  0.5846371288
##  [781] -0.3366263877 -1.5092569656  0.5822554167  0.2491986598  0.1114469410
##  [786]  0.9735710725  0.8321976536 -0.3655195885  0.5318305555  0.0450241576
##  [791] -0.3701849124  0.3193773017  0.4094224355  0.2567892912 -0.6621987598
##  [796]  0.0997143057 -0.3218395945  0.7983067751  0.8235531153 -0.3652998828
##  [801]  0.9299755200 -0.1218463561 -0.7734904229  0.3637095984 -0.3743316164
##  [806]  0.1336162452  0.5417454636  1.1920817742 -0.0957200324  0.1810953455
##  [811]  0.1610298254 -0.2349020662  0.8303666835  0.8588886409 -0.5000158870
##  [816]  0.1751895568 -0.6835512876  0.2317685850 -0.1433282898 -0.3331576850
##  [821]  0.3373673283  1.2758259154  0.1341344116 -0.5229512169  0.0542749033
##  [826]  0.0081104729  0.5445041702  0.0742639140 -0.0139691000 -0.2287284279
##  [831]  0.0924187037  0.0797883397  0.7032115180  0.6169072048 -0.2567615637
##  [836] -0.1855358110 -0.5316616430  0.5891954998  0.0828502564  0.5923032567
##  [841]  0.1769289391  0.0527175889 -0.6685737131  0.5667896915  0.2201410708
##  [846]  0.3495152655  0.4928186741  0.3468648171 -0.2242394385  0.0906844328
##  [851] -0.2195039087  0.0471453361  0.9113933558  0.1038553328 -0.2397324783
##  [856]  0.3228166597 -0.1489527791  0.2854361644 -0.8476497553 -0.5182455356
##  [861]  1.0239247852  0.1751895568  1.1391056053  0.6798804955 -0.4513218507
##  [866] -0.5923217529 -0.0866919028  0.0020634503 -0.0258710884  0.0683359225
##  [871]  0.1859187233  0.9323307802  0.3145733901  1.3131661783 -0.1783520227
##  [876] -0.5508063316  0.2009677921 -0.1174186750  0.7636454741  0.8658206283
##  [881]  0.3375115777  1.1037726218 -0.4728379023 -0.0353688441  0.2449779768
##  [886] -0.1602373410  0.5771768543  0.5038870817 -0.0479475452  0.5586123971
##  [891] -0.2101170969 -0.7082192658  0.3527175165 -0.1862419661  0.4908884125
##  [896]  0.2614642055 -0.2620938956  0.3540108370 -0.4039472844  1.7370538298
##  [901]  1.0248596910  0.5522970376  0.0612628381 -0.0053872293  0.8139157918
##  [906] -0.3838371761  0.3780077928  0.3752243618  0.5207568342  0.2600752417
##  [911]  0.0757564135 -0.1099368446 -0.2079791418  0.2236008416 -0.2553085116
##  [916]  0.1951395659 -0.4058202883 -0.0640136587  0.7332664270 -0.1747436811
##  [921] -0.1772547963  0.1451876224  0.5538576928 -0.9071634008  1.1802436794
##  [926]  0.6218612886  0.9411432949  0.7027668335  0.0514318249 -0.0861825653
##  [931] -0.5903437489  0.0908286652 -0.3354193163  0.9899463269 -0.4192105833
##  [936]  0.1942076589 -0.5235259226  0.6284046601 -0.0494591861  0.0765005396
##  [941]  1.3799283188  0.2096921259  0.7738556306  0.0361199054 -0.2670332506
##  [946]  1.0085072117 -0.4396086098  0.1313049027  0.2624841001  0.1876145480
##  [951]  0.1592358118  0.4052979193  0.1521728951  1.0525471759  0.2834149334
##  [956]  0.3324947473 -0.5887638732  0.7027668335  0.3203796408 -0.0438198312
##  [961] -0.2664670564  1.0489768823 -0.3373830711  1.1635390814 -0.1731265175
##  [966] -0.1635783837  0.6614288691  0.0904366020  0.1994612388  0.2270687937
##  [971]  0.6505128963  0.3523097040 -0.0419767062 -0.5594196659  0.0872967746
##  [976]  0.0134748119  1.2310634013 -0.7830625055 -0.0964470704 -0.2842322718
##  [981]  0.9598637413  0.5246258435  0.5557135609 -0.3921889075  0.7058714444
##  [986]  0.5964410512 -0.9555373162 -0.2864730775  1.0679001597 -0.7867333600
##  [991] -0.6397141142  0.3852527181 -0.2567615637  0.1236064582  0.6971436540
##  [996]  1.8247430608  0.2792418011  0.6005846940  0.5399837920 -0.3580143691
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
##   0.58290308   0.25626285 
##  (0.08103743) (0.05729766)
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
## [1]  0.81355193  0.53791712 -0.55707969  0.08067673 -0.27037212 -0.65389044
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
## [1] -0.0033
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9251243
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
## t1*      4.5 -0.005905906   0.9175614
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 3 4 5 6 7 
## 2 3 1 2 1 1
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
## [1] 0.0519
```

```r
se.boot
```

```
## [1] 0.901637
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

