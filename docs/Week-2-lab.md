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
## 0 1 2 4 8 9 
## 1 3 1 2 1 2
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
## [1] -0.0419
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
## [1] 2.614805
```

```r
UL.boot
```

```
## [1] 6.301395
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.7   6.2
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
##    [1] 3.3 5.9 3.7 6.0 5.1 4.8 5.1 4.7 5.2 3.6 5.6 5.6 3.7 3.7 5.9 4.4 4.2 3.6
##   [19] 4.0 4.7 5.4 5.0 4.9 3.8 5.3 4.5 4.6 5.2 4.7 6.4 5.2 2.8 5.0 5.1 3.8 5.4
##   [37] 3.5 3.7 5.9 4.2 5.9 4.1 4.7 4.5 4.8 6.3 4.8 3.7 5.9 4.6 4.1 3.8 3.8 4.3
##   [55] 5.0 2.5 4.2 3.3 4.2 4.7 4.6 4.1 4.6 4.9 3.8 5.1 3.8 4.6 4.8 5.3 5.5 4.9
##   [73] 3.4 4.0 4.7 5.6 4.1 2.9 4.6 2.9 5.1 5.5 5.0 4.6 6.3 4.0 4.7 3.2 5.1 5.5
##   [91] 3.6 5.9 4.5 4.2 4.6 5.1 4.6 3.6 4.9 4.3 4.7 5.4 4.9 6.3 3.7 6.5 4.3 5.0
##  [109] 3.9 4.7 4.7 5.3 3.9 4.3 3.7 4.5 2.9 4.5 4.0 3.6 4.7 6.5 4.8 5.4 3.5 4.1
##  [127] 3.9 4.2 6.2 4.8 5.8 4.5 5.6 4.6 4.2 4.4 3.7 5.8 3.2 4.1 5.1 3.3 5.8 4.2
##  [145] 5.8 5.1 4.4 5.3 5.5 3.8 5.2 4.8 4.1 4.2 6.0 4.6 4.7 3.9 4.3 3.1 3.4 6.3
##  [163] 5.1 5.1 4.8 5.2 3.9 5.1 4.5 3.6 4.5 4.8 5.0 4.7 3.1 4.9 3.9 3.7 4.0 3.8
##  [181] 4.5 2.9 3.8 4.8 3.0 5.4 4.6 4.8 4.6 4.6 4.4 3.2 2.2 5.1 5.5 4.1 4.3 3.7
##  [199] 3.6 5.1 5.2 5.4 3.6 3.1 3.8 4.8 4.1 4.5 4.7 4.4 4.8 4.8 5.6 4.4 3.2 3.5
##  [217] 4.6 4.6 4.3 4.1 4.2 5.2 4.1 4.1 6.6 5.4 6.1 4.9 4.8 4.8 3.3 3.8 3.8 5.0
##  [235] 4.7 5.1 3.4 4.9 3.6 6.1 4.9 4.7 4.4 4.5 3.8 4.3 4.7 4.4 5.8 4.6 5.3 3.9
##  [253] 4.2 4.2 3.7 4.1 5.5 5.5 4.3 5.6 4.3 3.9 3.2 3.4 5.5 5.2 3.5 4.0 4.2 4.3
##  [271] 5.0 4.0 5.4 4.8 4.2 4.0 4.3 3.8 4.5 4.6 2.9 4.5 5.6 4.2 3.8 5.3 5.3 5.3
##  [289] 3.2 4.7 3.9 5.0 5.8 5.4 4.2 5.0 5.8 4.1 4.2 4.8 5.1 5.3 3.5 5.1 4.8 5.3
##  [307] 5.4 5.4 5.6 4.0 4.9 4.6 5.5 3.4 3.2 3.3 3.5 3.9 6.6 5.4 5.0 2.2 5.3 5.0
##  [325] 5.6 5.3 3.6 4.0 4.7 5.1 3.6 4.3 3.0 4.8 5.9 4.5 5.7 6.2 4.1 2.6 4.9 4.5
##  [343] 4.6 3.4 3.7 5.0 5.1 5.1 2.8 6.2 7.0 5.0 2.9 4.1 4.5 3.5 4.3 4.6 5.6 3.7
##  [361] 4.0 4.9 4.3 6.1 4.0 3.7 5.6 4.7 5.8 5.3 5.0 5.2 3.7 5.1 2.6 4.8 2.5 4.4
##  [379] 3.8 5.4 5.7 3.6 4.5 5.3 4.6 3.5 3.9 4.2 4.6 4.3 5.8 4.9 4.0 3.9 5.7 4.4
##  [397] 4.2 4.8 4.9 3.8 6.5 3.3 5.1 2.9 4.3 3.5 5.8 3.9 6.0 3.9 4.9 4.1 6.6 3.5
##  [415] 3.8 5.3 5.2 4.7 4.3 4.9 5.6 3.7 4.4 4.0 5.3 4.5 5.6 5.2 4.7 3.7 5.1 4.7
##  [433] 5.1 4.1 4.4 4.1 4.8 4.9 4.0 4.2 4.8 4.5 3.4 6.0 5.1 6.2 3.9 3.1 5.1 4.4
##  [451] 4.7 5.8 4.7 3.9 4.3 3.9 3.9 4.4 5.1 4.8 4.2 4.7 4.2 4.0 3.9 3.8 5.1 5.8
##  [469] 5.5 4.3 4.2 5.3 4.9 4.5 5.3 5.3 4.6 4.7 4.9 4.0 5.7 4.8 5.3 2.3 4.0 4.2
##  [487] 4.7 3.6 3.6 4.3 4.3 5.6 5.8 2.9 4.2 5.1 4.9 5.5 5.4 4.8 4.0 5.3 5.2 5.3
##  [505] 3.5 3.5 4.1 5.4 4.2 4.5 3.0 3.0 5.4 3.2 6.0 5.2 6.3 3.1 2.7 5.2 4.5 4.6
##  [523] 3.8 4.7 4.0 4.5 3.3 4.2 5.8 5.1 5.0 5.0 4.3 5.3 3.3 3.8 4.6 4.2 2.6 2.9
##  [541] 5.2 3.9 5.0 3.5 3.1 4.9 4.5 3.6 4.8 3.3 4.4 6.8 4.3 4.1 2.9 5.0 3.4 4.0
##  [559] 4.7 5.3 5.9 4.0 4.4 5.7 5.4 4.8 2.4 2.3 4.6 4.6 5.1 4.7 4.9 3.6 4.6 4.0
##  [577] 4.5 5.1 6.0 4.7 4.0 5.7 4.3 4.8 2.1 4.4 4.5 5.2 4.8 2.8 4.2 5.3 4.7 3.9
##  [595] 5.2 4.8 6.0 5.9 6.7 4.5 4.6 5.1 4.1 4.2 4.2 2.8 4.9 4.6 5.4 4.6 4.9 3.6
##  [613] 4.6 4.7 5.5 6.0 3.9 3.8 3.7 3.9 4.1 4.1 3.9 4.1 3.9 3.3 5.4 3.7 3.5 3.2
##  [631] 4.0 4.1 5.0 4.1 3.8 4.7 4.3 5.6 5.2 3.3 5.3 4.8 5.0 3.7 4.5 3.7 3.2 3.9
##  [649] 2.8 6.6 4.7 3.2 4.1 4.3 4.3 3.5 3.2 6.0 4.5 4.7 3.2 4.5 4.4 4.5 5.1 3.2
##  [667] 4.8 3.3 3.9 4.6 4.0 5.5 3.6 5.7 4.0 5.0 4.3 3.9 3.9 5.0 5.3 4.4 3.9 4.8
##  [685] 4.0 4.2 3.3 2.8 4.2 4.6 4.4 4.4 3.7 3.3 4.3 3.9 5.3 3.2 5.5 3.4 5.0 5.3
##  [703] 3.8 5.4 4.0 5.2 5.1 4.5 3.8 3.5 5.2 5.7 4.1 2.8 4.2 6.3 4.1 5.5 3.5 4.7
##  [721] 6.3 4.8 4.4 4.6 6.2 4.9 4.2 4.8 4.2 3.8 2.1 4.8 5.0 6.3 4.0 4.7 4.1 5.0
##  [739] 3.4 2.7 4.3 4.8 4.3 3.6 3.6 5.7 4.5 3.5 5.0 4.0 6.0 5.2 4.5 5.8 3.8 4.4
##  [757] 4.3 5.2 3.9 3.2 3.8 4.1 5.2 3.9 2.6 3.8 4.4 2.8 4.6 3.8 5.6 4.6 3.7 5.1
##  [775] 4.1 5.4 4.5 3.9 4.5 5.4 4.9 3.1 2.4 4.6 2.8 3.9 4.6 5.1 3.0 4.6 5.1 3.3
##  [793] 3.8 5.7 3.6 3.3 4.9 4.7 3.8 4.7 4.4 3.8 3.2 5.8 4.0 5.8 3.2 5.5 4.0 4.1
##  [811] 3.5 4.8 5.0 5.7 4.3 3.8 5.0 4.4 5.9 3.5 3.7 3.8 4.6 5.1 5.9 6.1 5.7 4.9
##  [829] 4.4 3.7 3.2 4.4 6.0 4.9 5.0 5.0 3.9 4.5 4.6 5.3 6.2 4.7 2.0 5.6 3.8 4.7
##  [847] 5.5 4.8 5.3 4.6 3.2 3.1 5.3 3.9 4.0 5.7 3.7 3.6 4.4 5.2 5.1 4.2 3.3 5.1
##  [865] 4.2 3.9 4.4 2.6 4.3 3.5 3.8 3.8 3.3 6.0 3.1 4.7 4.0 2.0 4.2 5.5 6.4 4.8
##  [883] 4.3 4.5 3.9 3.2 4.6 4.4 3.2 2.8 4.1 3.8 4.3 5.4 3.8 4.0 4.2 5.0 3.0 4.8
##  [901] 4.2 4.1 4.9 5.0 5.7 4.2 3.9 4.0 3.7 5.7 4.0 4.7 4.3 4.5 4.4 5.2 4.4 4.0
##  [919] 3.9 5.6 5.4 4.2 4.0 3.7 2.6 5.1 4.2 5.1 4.2 4.3 4.2 4.5 3.4 4.3 5.0 5.8
##  [937] 4.8 3.9 3.6 4.1 3.3 4.8 5.1 5.1 3.9 4.1 2.7 3.9 4.3 5.1 5.1 4.6 2.6 3.8
##  [955] 6.1 3.6 3.3 3.3 3.7 4.4 5.5 3.4 4.9 3.2 5.1 5.3 3.6 5.9 4.4 4.0 3.9 6.0
##  [973] 4.3 4.0 4.4 4.8 4.1 3.3 4.7 4.6 3.7 4.3 5.5 3.8 5.1 5.6 4.6 3.4 6.3 6.3
##  [991] 2.7 4.0 4.6 3.4 4.3 4.2 3.1 4.1 3.9 5.4
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
##    [1] 3.4 2.9 4.6 4.4 5.4 4.6 5.2 5.2 4.9 5.7 2.8 5.3 3.5 4.0 4.5 5.8 4.0 3.0
##   [19] 4.2 5.5 5.0 5.5 2.9 6.3 5.3 3.6 5.7 6.1 3.5 4.3 4.0 4.1 3.7 4.6 5.1 2.8
##   [37] 5.8 3.6 4.5 3.7 5.1 6.3 2.8 3.4 4.1 6.0 4.8 3.9 5.3 3.9 4.8 4.8 3.3 5.2
##   [55] 4.1 4.8 4.4 6.1 2.5 4.4 4.3 4.6 5.7 4.8 4.9 3.1 2.8 3.5 4.3 4.0 5.1 5.8
##   [73] 3.7 5.9 5.1 3.2 6.1 4.6 5.3 5.8 4.1 4.8 4.8 3.2 3.8 5.1 4.1 4.3 4.9 3.4
##   [91] 4.6 5.3 4.3 5.2 5.3 5.0 3.3 6.3 4.4 5.5 4.0 3.2 4.6 5.0 3.8 4.4 4.1 4.2
##  [109] 5.8 5.4 4.5 4.5 3.9 4.6 4.1 4.2 4.0 4.9 3.4 4.7 3.5 5.6 4.5 4.8 3.5 5.4
##  [127] 4.3 3.6 3.1 5.6 4.7 4.8 5.3 3.3 4.5 5.6 5.9 4.9 4.4 2.7 5.6 3.1 3.8 4.6
##  [145] 4.9 3.6 3.7 4.2 4.1 5.5 3.4 3.6 4.6 4.5 4.6 4.1 4.4 3.8 5.2 5.0 5.5 4.7
##  [163] 4.6 3.8 3.3 4.4 5.4 4.8 4.5 4.6 4.4 4.3 5.7 5.2 6.0 5.0 3.4 5.0 4.2 4.1
##  [181] 4.8 4.7 3.8 4.2 4.4 4.1 3.8 6.3 4.3 3.2 4.3 4.6 4.9 4.6 5.4 5.3 3.7 6.1
##  [199] 4.2 4.3 5.4 2.9 4.0 4.7 5.7 4.4 4.5 4.8 4.7 4.8 5.8 6.0 5.4 3.4 5.4 6.7
##  [217] 4.3 4.6 3.4 3.7 4.7 3.0 5.8 4.7 2.9 6.2 4.8 4.8 3.2 4.9 4.2 3.0 4.3 4.0
##  [235] 2.6 4.0 5.7 4.2 3.5 4.5 3.8 3.9 5.3 5.6 5.8 4.9 4.8 5.4 3.7 3.3 5.3 5.2
##  [253] 5.0 5.4 5.8 3.2 6.5 4.3 4.5 4.9 4.6 5.2 4.9 3.2 3.2 3.9 4.3 5.3 5.2 5.5
##  [271] 4.4 2.1 3.9 4.3 6.4 4.6 5.2 3.6 4.2 4.4 4.8 4.9 3.2 4.4 3.2 5.1 4.8 4.1
##  [289] 5.4 2.7 3.0 4.0 5.0 4.7 4.6 5.4 3.5 4.7 3.6 4.9 5.5 5.1 4.4 4.3 3.3 4.6
##  [307] 5.5 5.0 4.7 6.4 5.5 5.5 4.8 3.9 3.6 3.5 5.3 4.5 6.5 4.2 3.2 5.0 2.2 3.6
##  [325] 4.2 3.7 4.4 3.7 5.5 4.6 3.6 4.7 3.8 3.7 3.3 6.9 4.4 4.9 3.0 5.0 5.0 4.5
##  [343] 4.7 3.9 3.5 3.7 4.7 4.6 4.0 3.5 4.0 4.0 4.6 4.8 3.5 3.0 4.3 5.6 5.5 3.1
##  [361] 3.5 4.3 4.1 4.4 5.5 3.3 3.7 4.8 5.2 4.6 2.1 5.2 3.2 4.9 5.7 5.3 3.8 3.1
##  [379] 4.4 6.5 4.5 2.9 4.1 5.5 3.3 3.8 3.7 3.3 4.8 6.4 4.6 5.8 6.3 5.1 4.4 4.6
##  [397] 4.4 6.3 2.1 3.1 5.7 4.0 5.2 3.8 3.5 3.0 4.5 5.3 4.7 4.3 4.4 3.4 4.7 3.1
##  [415] 4.4 4.2 3.7 3.6 3.9 6.0 3.9 5.1 4.7 4.8 3.2 2.5 3.3 3.1 6.5 4.8 3.1 4.5
##  [433] 3.6 5.4 6.8 4.9 3.9 4.7 6.0 4.4 4.9 4.4 5.6 3.9 5.5 4.7 3.1 6.1 5.3 5.4
##  [451] 6.1 3.6 6.0 4.7 4.9 4.1 5.5 3.7 5.3 5.4 4.7 3.9 5.5 4.7 3.5 4.5 5.9 3.5
##  [469] 4.0 4.7 4.1 4.5 3.4 3.9 7.5 3.5 5.2 4.4 5.5 5.6 5.7 3.2 4.5 3.7 3.3 3.2
##  [487] 6.4 3.8 3.7 5.8 3.4 4.1 4.8 4.3 5.4 7.0 5.5 2.3 4.7 4.1 4.5 4.8 5.1 5.3
##  [505] 3.5 5.3 3.9 3.5 5.0 4.1 5.3 6.2 5.6 4.0 5.0 4.6 5.0 4.1 4.0 5.3 4.4 2.8
##  [523] 2.9 4.0 4.9 4.6 4.7 4.5 5.5 4.7 4.4 3.3 5.9 4.7 5.9 4.4 5.0 4.7 3.9 5.8
##  [541] 3.7 3.1 6.3 3.6 4.9 4.9 2.9 3.7 3.6 5.0 5.3 5.2 3.6 5.5 6.0 5.6 4.9 4.2
##  [559] 3.7 4.4 4.8 3.9 4.9 3.7 5.3 3.5 4.6 4.4 4.2 5.4 4.8 3.3 5.5 4.0 5.3 3.4
##  [577] 4.4 5.4 4.7 5.0 4.2 3.2 3.3 4.8 4.5 5.9 4.9 4.1 5.5 3.6 5.3 3.5 3.3 4.5
##  [595] 5.3 4.6 4.2 4.8 4.7 5.2 3.5 4.5 4.3 2.0 5.3 4.7 3.2 5.4 4.4 4.7 5.1 5.4
##  [613] 5.2 5.5 4.6 5.7 5.3 4.0 5.5 5.4 4.2 3.2 5.4 5.7 4.0 4.7 5.0 2.5 3.2 3.9
##  [631] 4.3 5.2 4.5 6.7 3.2 3.0 4.7 5.7 4.9 2.8 5.1 4.4 5.5 4.5 4.5 4.0 3.2 5.6
##  [649] 4.8 4.1 5.3 5.1 4.6 4.9 2.8 6.7 4.8 4.5 5.0 5.0 5.7 3.4 5.2 3.1 5.6 4.4
##  [667] 4.1 3.3 4.2 5.4 4.7 5.3 6.1 3.6 3.9 3.7 3.7 4.0 4.5 5.6 5.8 3.7 5.1 5.3
##  [685] 5.0 3.8 4.5 3.7 3.6 3.7 3.7 4.4 5.6 4.3 5.2 3.6 2.7 3.8 1.8 5.0 4.9 3.7
##  [703] 5.2 6.0 6.0 4.6 3.3 5.3 6.1 4.6 5.5 5.4 4.9 4.9 5.3 4.2 4.9 5.2 3.7 4.1
##  [721] 5.1 2.9 3.6 5.1 5.0 5.5 6.1 4.6 4.3 5.0 4.6 4.1 4.5 5.8 5.0 5.0 4.4 3.1
##  [739] 4.7 4.8 4.3 3.7 4.2 6.3 4.9 2.5 4.9 4.7 4.6 3.8 6.5 6.7 5.2 5.2 2.4 2.8
##  [757] 5.4 4.5 4.8 5.3 3.4 5.2 5.1 5.2 4.3 3.7 5.4 6.2 3.7 4.9 4.7 4.9 3.2 4.8
##  [775] 5.1 3.8 4.7 4.7 3.8 5.2 4.1 5.4 3.0 5.9 4.0 5.1 3.3 3.3 6.0 3.4 4.5 4.5
##  [793] 3.4 5.3 3.9 4.5 3.1 3.6 5.4 4.0 5.4 5.6 4.2 3.6 5.9 6.0 5.0 4.9 2.9 4.6
##  [811] 6.0 4.1 3.7 4.9 4.9 4.4 5.1 4.0 4.9 5.0 4.3 5.8 3.9 4.6 4.7 4.0 3.5 3.4
##  [829] 4.6 5.2 3.5 5.3 4.2 4.4 5.4 5.8 4.3 5.1 4.6 5.4 3.8 5.1 4.8 2.9 4.2 3.2
##  [847] 3.4 5.5 5.6 2.6 6.0 2.4 3.7 3.3 4.8 5.3 4.1 7.1 4.1 5.4 3.6 5.0 2.7 6.5
##  [865] 5.8 5.4 4.9 2.9 4.0 3.7 3.5 4.5 5.3 4.3 5.5 3.4 4.1 5.2 4.5 3.0 3.3 3.3
##  [883] 4.7 4.9 5.4 2.2 6.1 2.6 2.5 5.4 5.3 4.1 4.7 5.3 4.2 4.2 4.0 4.7 4.0 4.7
##  [901] 3.3 3.6 3.9 5.0 4.4 6.2 3.0 3.6 5.5 6.8 3.9 3.2 4.5 5.0 5.1 3.1 5.0 4.3
##  [919] 4.9 5.5 3.9 4.9 4.2 2.9 3.2 3.9 3.9 3.2 5.0 4.5 4.0 5.2 2.8 4.3 2.4 3.9
##  [937] 5.4 4.1 4.9 4.4 5.4 4.2 4.1 5.0 4.3 4.9 4.1 4.6 4.4 3.7 4.0 4.2 5.1 3.6
##  [955] 5.2 5.6 4.2 6.1 4.4 2.6 2.8 5.4 4.6 5.0 4.3 4.0 3.9 4.4 5.7 6.3 3.5 5.5
##  [973] 3.6 6.7 4.8 2.2 3.7 4.9 3.8 5.3 4.1 5.5 5.6 5.5 3.4 5.4 4.2 6.0 4.2 3.6
##  [991] 4.9 5.6 3.4 3.6 4.1 3.8 5.5 2.8 5.1 4.0
## 
## $func.thetastar
## [1] 4e-04
## 
## $jack.boot.val
##  [1]  0.57035928  0.37941176  0.27595308  0.20279188  0.04198895 -0.09212121
##  [7] -0.15119760 -0.24126074 -0.43260274 -0.49353100
## 
## $jack.boot.se
## [1] 0.9986959
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
##    [1] 3.7 4.3 6.1 3.7 5.1 4.8 5.1 4.1 3.4 4.3 4.0 4.1 2.9 4.5 4.1 4.3 5.2 4.7
##   [19] 4.3 4.8 4.6 5.0 3.5 4.7 4.5 4.7 4.8 4.6 3.6 4.6 4.4 5.5 3.7 4.4 3.0 3.8
##   [37] 4.1 4.2 4.9 6.0 4.2 5.4 4.2 5.6 4.7 5.2 6.6 4.5 5.3 3.9 4.1 3.4 4.0 4.4
##   [55] 4.9 5.2 4.2 3.4 3.5 3.9 3.4 3.4 4.9 5.1 3.7 4.3 3.5 4.1 4.9 4.7 3.4 5.0
##   [73] 4.0 5.7 5.0 5.0 2.4 5.9 4.0 4.3 4.8 4.9 2.9 7.7 2.9 4.2 3.6 5.0 4.5 2.0
##   [91] 4.7 4.3 5.7 4.6 5.2 3.7 5.5 3.6 4.2 3.9 4.0 4.8 3.8 4.9 4.2 3.4 4.4 4.3
##  [109] 5.2 5.0 2.2 4.1 4.7 5.0 5.4 4.2 3.7 4.5 5.6 5.2 3.6 2.7 4.1 3.7 3.9 6.3
##  [127] 5.3 2.9 4.8 5.8 3.3 5.3 3.1 4.1 5.8 2.3 3.4 5.0 4.4 4.0 3.9 3.1 4.4 5.5
##  [145] 4.6 4.5 3.5 5.7 5.5 3.1 5.1 5.2 4.8 5.3 3.3 4.1 5.1 4.3 5.4 5.7 4.5 5.4
##  [163] 5.4 4.0 5.3 4.3 3.3 5.2 4.6 5.6 3.8 4.8 4.4 3.2 4.0 6.8 4.0 4.8 4.4 5.3
##  [181] 5.1 4.8 4.6 4.9 4.5 5.0 3.8 4.1 5.9 4.0 6.3 5.1 4.8 5.9 3.8 4.8 5.1 4.0
##  [199] 4.4 4.7 3.9 4.4 3.3 4.6 3.6 4.0 5.0 3.5 3.0 4.3 4.2 4.5 3.6 5.4 4.5 5.3
##  [217] 4.6 4.5 5.5 3.4 4.2 6.0 4.5 4.5 5.5 4.7 3.1 3.1 5.5 5.0 5.7 5.9 6.5 5.1
##  [235] 5.6 5.8 4.7 5.5 4.8 4.6 4.5 3.1 3.5 3.8 4.3 4.0 3.0 5.2 4.2 4.1 4.9 4.5
##  [253] 5.0 4.1 4.6 3.9 4.9 4.1 2.5 4.8 5.1 4.1 5.9 4.9 5.9 3.3 5.1 3.4 4.6 4.8
##  [271] 5.7 4.6 5.4 5.1 5.6 4.7 4.8 5.6 6.3 6.3 3.0 5.1 5.2 5.2 4.8 5.0 3.4 5.4
##  [289] 4.0 4.7 3.9 3.3 4.6 4.8 3.8 5.0 3.0 6.0 6.0 5.8 6.2 5.3 4.8 4.2 4.2 3.8
##  [307] 4.4 4.4 4.0 6.0 4.0 4.8 4.5 4.7 5.0 6.1 5.1 4.6 5.2 5.0 4.1 3.8 4.5 4.9
##  [325] 5.2 4.7 3.6 5.3 4.6 4.9 3.2 4.3 2.0 4.9 4.3 4.0 3.1 4.3 5.8 4.6 5.7 6.1
##  [343] 5.5 5.1 4.2 5.2 5.0 4.4 4.7 4.4 4.4 4.9 5.3 3.8 4.4 3.0 4.9 4.2 4.4 5.6
##  [361] 5.3 5.1 5.2 5.1 3.9 4.0 3.7 6.3 3.3 3.9 3.6 3.9 6.4 5.9 4.9 5.3 5.3 3.8
##  [379] 4.9 4.0 5.1 4.5 4.2 4.5 4.9 4.4 4.7 4.0 3.6 5.7 4.0 5.3 2.5 4.8 5.9 3.9
##  [397] 4.4 5.3 6.3 4.6 4.8 4.3 4.2 6.8 6.0 5.7 3.0 6.0 4.1 5.1 5.4 3.1 4.1 2.3
##  [415] 4.1 4.3 4.6 3.6 4.8 5.5 3.7 3.1 4.1 4.1 4.0 6.6 4.4 4.0 5.9 3.5 3.6 4.1
##  [433] 6.2 4.3 4.9 5.2 6.2 3.7 2.6 3.3 5.0 4.3 3.7 2.7 4.3 4.2 3.7 4.8 5.2 4.6
##  [451] 3.1 5.2 4.0 5.4 3.7 6.6 5.7 3.9 5.9 3.3 5.1 4.9 5.2 4.8 4.3 5.2 4.3 2.9
##  [469] 5.0 4.2 5.9 4.5 5.2 2.6 4.0 4.5 4.9 5.6 5.1 3.3 6.0 3.6 5.4 4.9 4.0 4.3
##  [487] 3.2 4.4 6.4 4.7 4.6 4.9 5.6 4.3 5.7 5.3 5.6 3.2 4.9 3.9 3.7 5.1 3.7 4.3
##  [505] 5.7 4.0 5.1 5.8 5.5 4.7 3.0 4.7 3.8 4.9 3.1 5.4 5.1 6.1 4.6 4.4 4.9 4.5
##  [523] 4.9 5.5 3.9 3.5 3.8 5.6 5.5 3.7 3.9 5.1 5.0 2.4 5.0 4.8 4.4 5.9 5.3 4.8
##  [541] 5.0 3.7 5.7 4.5 5.6 6.1 5.0 5.1 5.0 3.4 4.2 5.8 3.5 2.6 5.2 6.3 3.6 3.7
##  [559] 2.7 4.2 5.5 4.5 3.9 4.0 2.0 4.5 6.1 5.0 5.9 4.3 3.7 2.9 3.0 5.5 5.7 4.0
##  [577] 4.9 5.1 6.4 5.1 4.3 4.3 3.3 4.7 3.3 4.0 5.0 3.6 4.7 5.4 5.3 4.3 5.1 5.7
##  [595] 3.5 3.5 5.2 4.4 4.3 5.9 5.9 6.3 3.9 2.9 5.3 4.2 2.5 5.5 4.5 5.3 6.0 4.9
##  [613] 3.3 5.9 4.2 4.7 4.5 4.3 4.8 3.8 5.3 4.7 4.7 5.6 5.6 5.4 2.8 4.6 4.3 5.8
##  [631] 5.8 3.6 4.4 3.3 5.0 4.5 4.2 4.1 3.6 4.2 4.7 4.5 4.4 5.0 3.4 5.6 5.1 4.6
##  [649] 4.2 4.8 5.1 2.7 3.3 5.3 6.1 2.9 6.2 3.9 3.3 3.5 5.3 3.5 3.6 4.7 3.7 4.3
##  [667] 3.4 3.8 5.2 4.5 4.9 4.5 6.3 3.7 5.7 5.1 5.1 3.1 4.1 4.9 5.7 5.6 4.0 5.1
##  [685] 6.4 5.2 5.4 5.6 4.8 3.6 3.0 4.9 5.5 4.7 5.4 6.2 4.4 4.5 3.7 3.8 4.4 4.6
##  [703] 4.9 3.8 3.7 3.9 4.6 5.2 4.7 3.5 4.3 5.9 5.0 3.1 4.5 3.8 3.9 4.4 3.2 4.2
##  [721] 4.9 3.2 4.7 6.0 3.5 6.1 5.2 4.1 5.6 4.4 4.6 3.9 4.5 3.4 3.6 3.5 4.8 4.6
##  [739] 4.4 4.5 3.9 3.8 4.4 4.0 3.9 4.0 5.4 5.6 4.9 5.3 4.4 5.3 5.8 4.0 3.2 4.6
##  [757] 4.4 4.2 5.2 3.5 4.4 5.3 4.6 4.7 5.0 4.5 5.4 5.4 4.3 4.0 5.3 3.4 4.0 2.8
##  [775] 4.0 4.8 5.3 5.3 3.3 4.6 6.2 4.9 4.0 3.7 5.3 5.9 5.3 4.1 4.5 4.9 4.8 5.3
##  [793] 3.1 5.3 6.0 4.8 5.8 4.7 3.7 4.8 3.9 3.6 3.6 4.3 3.8 5.5 5.6 5.0 4.2 4.6
##  [811] 4.5 4.4 6.2 3.9 4.5 4.4 4.9 4.0 4.9 4.2 4.6 4.5 3.4 5.0 6.3 3.6 3.1 4.4
##  [829] 4.4 4.1 4.7 4.8 4.9 3.9 4.8 3.2 4.7 4.0 5.4 4.2 4.2 4.5 4.3 3.9 5.6 2.4
##  [847] 3.8 4.4 4.2 4.6 3.4 5.5 4.8 4.8 5.4 3.0 4.3 3.7 5.7 4.5 4.9 6.5 3.3 6.0
##  [865] 6.6 5.8 5.0 4.3 4.1 5.4 6.8 5.1 4.9 3.7 5.1 6.3 5.0 2.9 5.6 3.3 3.6 4.6
##  [883] 3.5 5.7 5.5 4.2 4.7 4.0 4.2 7.4 6.1 5.9 6.7 4.5 3.8 5.5 3.1 4.6 3.0 5.7
##  [901] 3.1 5.1 4.5 5.0 2.7 6.1 4.8 2.3 3.6 4.7 3.5 3.3 3.6 4.6 4.8 6.5 5.2 4.0
##  [919] 3.7 3.2 3.9 3.9 6.8 4.9 5.2 5.0 5.1 4.8 4.8 5.5 2.9 4.9 4.0 4.2 4.3 3.7
##  [937] 4.2 4.0 4.7 5.1 5.9 5.8 4.6 3.7 3.6 4.9 5.0 5.5 3.5 3.5 5.0 4.6 3.6 4.2
##  [955] 3.8 3.2 4.4 5.8 5.1 3.7 6.0 2.8 5.1 3.5 4.3 4.3 4.7 4.6 5.7 6.1 4.7 4.8
##  [973] 6.4 3.9 5.5 3.0 3.8 3.9 5.0 5.3 5.6 4.5 4.0 4.9 5.3 3.9 5.6 3.8 6.7 3.6
##  [991] 5.2 5.5 4.6 3.3 5.3 4.6 4.2 3.8 5.0 4.0
## 
## $func.thetastar
## 72% 
## 5.1 
## 
## $jack.boot.val
##  [1] 5.6 5.4 5.4 5.2 5.1 5.1 4.8 4.9 4.8 4.6
## 
## $jack.boot.se
## [1] 0.9044888
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
## [1] 0.6186845
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
##   2.801778   4.651230 
##  (1.185992) (2.156056)
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
## [1] 0.5004849 0.8070365 1.2343394 1.6292360 0.9412622 0.6155377
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
##    [1]  0.6908427233  0.2446373181  1.2955345349 -0.4915262870  1.6131529855
##    [6]  0.2694614502  0.9782722350  0.1983332207  0.1565608003  0.4111577753
##   [11] -0.0492284294  0.1821336479  0.0892948516  0.6978514386  0.2851399066
##   [16]  1.7796477988  0.0139820421  0.1717785294  0.5933978057  0.1135941633
##   [21] -0.0355612619  0.5853406088  0.7218070471  0.6396212011 -0.1730136792
##   [26]  0.5961339946  0.6384559157  0.9736703300  0.7552795677 -0.3753419143
##   [31]  0.2376160642  1.0516098641  1.0521081552 -0.0467414950  0.1842857486
##   [36]  1.1803807245  0.5761502857  0.8673552155  0.8287558736  1.4477519102
##   [41] -0.0087089706  0.0781944844  0.0947468753  0.3953227134  0.7569445582
##   [46] -0.0738512631  0.3448732011  0.2525455518  0.6201668022  1.2499006730
##   [51]  0.1225659526 -0.0214302757  0.6127685623  0.5395696273 -0.2167870880
##   [56]  0.3423037918  0.2119147551 -0.2313176988  1.1707396715  0.9820636986
##   [61]  1.4417720278  1.1751491102  0.8475259926  0.3605827616  2.2145620665
##   [66]  0.4271564842  0.6134207040  0.5914819736  0.3959302994 -0.6051368735
##   [71]  0.6104966657 -0.2877049112 -0.0981185920  0.5062272101  0.5918752485
##   [76]  0.6217169555  0.8402336152  0.6306442142 -0.3859341809  0.3093270682
##   [81]  0.3476482524  0.6179751140  0.8540427175  0.8454732743 -0.5086645495
##   [86] -0.6210145986  0.0251183242  1.1674191389  0.9638412541 -0.3886769464
##   [91]  0.2149014243  0.1133971701  1.0868201334 -0.3772528523 -0.3417280249
##   [96]  0.5824601510  1.4766546966  0.1967826716  0.4138986051  1.1392427173
##  [101]  0.1707899391  0.2381915352  0.7934310664 -0.1997131576  1.4336662098
##  [106]  2.1878432054  0.1796480395  0.1154649082  0.6809464255  0.6517883360
##  [111]  1.5486699232  0.3954525585  0.7402377784  0.9706720137  0.0168698715
##  [116]  0.2358094876  0.8207996431 -0.4014046785  1.0491079313  0.5531653683
##  [121]  0.2200379694  0.1186921052  0.8529532070  0.8276318484  0.5584742536
##  [126]  1.9031025064 -0.3487586318  0.5386344671  0.3471356182  0.0061522092
##  [131] -0.0430748244  0.3280278701  0.5646060222  1.4101545990  0.2789232321
##  [136]  0.0962866578  0.4191089375 -0.0587483221  0.3439868265  0.5995750486
##  [141]  1.4126957371  0.5691053344  0.7914092181  0.6070079782  0.0579660377
##  [146]  1.0556505151  0.1829178292  1.4233445732  1.2638122307  0.1104054060
##  [151]  0.6482112677  0.0580741600  0.9746193937  0.3004328070  0.4758297980
##  [156] -0.0214189538  0.7838427589  0.9157903973  0.9903967852  0.1290081919
##  [161] -0.1076704904  1.0502327166  1.8886042200  0.8520395315  1.5442245628
##  [166]  0.6143540867 -0.3741066544  0.3069281451  0.7292679351  1.4083639768
##  [171]  0.7883454193  0.7623267669  0.6889823722  1.5506038120  0.4005587036
##  [176]  1.0485198989 -0.3662254837  0.5257837300  0.5646060222  0.9711418503
##  [181]  1.5684649206  0.1052063006  1.8106562130  0.1591094356  0.3942847002
##  [186]  0.8466472832 -0.6584780589 -0.6088936682  0.1329794187  1.5859919206
##  [191]  0.4289637864  0.6618401633  0.8472154820  1.1763529271  0.3600789772
##  [196] -0.1934571320  1.4561334293  0.5804854499  0.1663409490 -0.0071835413
##  [201]  0.8646834932  0.5234381221  0.3407857123  0.4089401063  1.0643737252
##  [206]  1.0910574798  0.0277771682  0.8811283705  0.0166147754  0.8381621935
##  [211] -0.1951423514  0.4266299488  0.5257232290  0.8946335960  0.4257607725
##  [216]  1.0461063709 -0.1465911921  0.5502036449 -0.3766525264  0.6997614872
##  [221]  0.7269996819  0.2605257613  1.0769851557  0.9865406768  0.1200073307
##  [226]  1.5261504021  0.5711126342  0.4180934903  0.9290125466  0.2561990883
##  [231]  1.3022543032  0.9324083840  0.2860799334  0.0893838728  0.5518082127
##  [236]  2.5084480059  0.1002736215  0.4971971767 -0.1658935548  0.7090355212
##  [241]  1.2756167882 -0.1653401250  1.0446856502  1.2014531764  0.1034076916
##  [246]  0.0074216040  0.4025326451  0.8654593574  1.0718102597 -0.0617714416
##  [251]  0.5778366374  0.5533692817  1.4937620146  0.7276306864  0.8764015520
##  [256]  0.4158960454 -0.1944598116  0.6002361761  0.7434262615  1.8725913431
##  [261]  0.8482500839  0.8476105978  0.9821801331  0.1033505872  0.8689478333
##  [266]  0.1527600261 -0.2567535947  1.4665046319  0.5167422701  1.2010496340
##  [271]  0.4066416500  0.9372631448  0.3341088721 -0.3743571700  0.3374756245
##  [276]  2.6275419998  0.7064127531  0.2022382798  0.5486458096  0.5814435402
##  [281] -0.0177720516  0.0386236412  0.4671646092  0.5916342596  1.3048815130
##  [286] -0.0638764790  0.5495624281  0.8386792704  0.8404808946  0.4054395667
##  [291]  0.7723011661 -0.2082232759 -1.8734439984 -0.0610832447  0.8819875634
##  [296]  0.5521339644  0.9697411003  0.8723953452  0.7295598412  0.7948972247
##  [301]  0.8542182257 -0.1026587569  0.0383169881  0.5676264953  0.9711418503
##  [306]  1.2737277241  0.5740705788  0.3881576685  0.4147503520  0.1653642577
##  [311]  0.9518698242  0.6199877524  0.2749405577  0.8574418403  0.2849445277
##  [316]  1.4295780660 -0.0499029076 -0.0019964579  0.3017690168  0.5825374999
##  [321]  0.2907186614  0.1816377749  0.0024994363 -0.2691585785  0.8593688007
##  [326]  0.4408390268  0.4415794945  0.9800066288 -0.0470901649  0.3211806888
##  [331]  0.4342501892  1.1078731016 -0.0713184709  0.4930136346  1.7639978856
##  [336]  0.7762679406  0.2559034001  0.8870406242  1.3038920343  0.0558281113
##  [341]  1.1200351990  0.1850240801  0.2531272847  1.1760991153  0.0919201079
##  [346]  0.2624424741  0.1799952139  0.4966430420  0.4695162637  1.3072735846
##  [351]  0.7821186126  0.6049582583  2.4792672173  0.0014778868  0.6347472057
##  [356]  0.7283431915  0.0146292583  0.6969403625  0.8774737197  0.2833131732
##  [361]  1.8389468636  0.5062808533 -0.0607438579  0.5730757857  0.0400242481
##  [366]  1.7475684792  1.4051423806  0.4826109444  0.2577149254  0.3938856643
##  [371]  0.5595683016  1.1000012700  1.6665751779  0.2530661776  1.4182389083
##  [376]  0.2151816933  0.6997614872  0.9733353746  1.6270055690  1.8141230480
##  [381]  0.4820422576 -0.1322241758 -0.1021240254  0.4889388152  0.3983068775
##  [386]  0.2654218609  0.3064574653 -0.2942836183  0.3641134684  0.7427138700
##  [391]  2.6078238941 -0.1354329268  0.1138622274  0.3741709611  1.8072873885
##  [396] -0.4585813881  0.0036011040  0.6334110662  0.6417156503  0.0284795394
##  [401]  0.1556484359 -0.9064068080  0.8687623000  0.2977815782  0.7821186126
##  [406]  0.2229112414  0.8544246142  0.5317929331  0.1091652422 -0.3658965766
##  [411]  0.5399051270 -0.3791470558  2.1968099375  0.4046984323  0.1232925985
##  [416]  0.7786611373  0.1919077971  0.9964860324  0.8453605495  1.0529677388
##  [421]  0.5164460809 -0.1050891501 -0.1359027108  0.2420104418  0.0773739423
##  [426] -0.1981277578  0.9842259315  1.5570491219 -0.1122358269  0.8590582253
##  [431]  0.7343576910 -0.1233463321  1.4663389309 -0.1600975639  0.8234998168
##  [436]  0.3863290261  0.0226611565  0.7120156959 -0.0965117454  0.0906917941
##  [441] -1.1258794204  0.0103157772  0.1140588840  0.5293308866 -0.0462918250
##  [446]  0.7751542042  0.2651963841  0.5078918631  0.0001251289  0.8576659979
##  [451]  0.0010232565  0.7713496466  0.1398342250  0.3676389940  0.4030461867
##  [456]  0.5666402378  0.0263920781  0.5860167638  1.6629638420  0.7445467066
##  [461]  0.8750231328  0.5179801512  0.5467989307  1.4525459684  1.0353823532
##  [466] -0.0455159156  1.0589039094  0.2184130588  0.2892784974 -0.0653015829
##  [471]  0.6288013574  0.5606839690 -0.0843211625  0.2227223892  0.8309067832
##  [476]  0.2695611237  0.3594931795  0.5315454900  0.4929082136  0.3499378236
##  [481]  0.6978881454 -0.0164892389  0.8404808946  0.5030699981  1.1751324822
##  [486]  0.1625146907 -0.3159606779  0.1973515612  1.0691149270 -0.1163413620
##  [491]  0.3665665778 -0.1410226967  0.3837329321 -0.0220603949 -0.1492242812
##  [496]  0.9083967992  0.7133749037  0.4298318797  0.6062160863  0.9077136182
##  [501]  0.1499804803  2.1796734298  0.4493777327  0.4985069931  0.6229672180
##  [506]  0.1333443833  0.2213334651  0.1616370465  0.0219777284  1.4632644858
##  [511] -0.0685662465  0.5827511601  2.4063230696  0.9698373627  1.4274882000
##  [516]  0.1287724999  0.4290397881  0.9753150901  0.8856236155  0.9658180499
##  [521]  1.1238685996  0.5913485463  0.6657682371  0.5679627130  0.2892784974
##  [526]  0.9668826649  0.6544675771  0.6029978730  0.6065932702  1.0672922751
##  [531]  0.6231097551  0.5042287072 -0.0302068021  0.4473388778  0.2866141150
##  [536]  0.9562478363  2.5365778876  0.1841558314  0.1361466668  1.0560837504
##  [541]  0.1862248038  1.1758567488  0.9817186839  0.9972858763  1.3186996594
##  [546]  0.0541471271  0.5296956125 -0.3022943567  0.1444308049  0.3939819351
##  [551]  0.1961923882  1.1727405377  0.7723555888  0.4115807599  0.4772615271
##  [556]  0.8593688007  0.6564299722 -0.3411578376  0.3088645078  0.1591035408
##  [561] -0.0331206310  0.0026899656  0.5697104181  0.5092712157  0.0451790792
##  [566]  0.6369684326  1.4749118128  1.0620538755  0.9778656789  1.5421920408
##  [571]  0.5992147693  0.4044357494  0.8554251250  0.5651767516  0.1298006842
##  [576]  0.9827021051  0.4250652056 -0.0183491888  1.1844910776  1.0724995823
##  [581]  0.4522311804  0.4560778289 -0.4638033717  0.1690011423  0.4945383758
##  [586]  0.8687277143  0.1051167197  0.5679907226  0.5930937229  0.1597838817
##  [591]  0.4539266993  0.3966660586 -0.2176448457  0.0074608857  0.9528336900
##  [596]  1.2676444591  0.2050345012  0.9799128315  0.8488354004 -0.2279873615
##  [601]  0.0111304639  1.3999260368  0.8208539456  0.7523670973 -0.2236102225
##  [606]  0.9705255351  0.5248684890  1.0697945524 -0.0512889973  0.7369764526
##  [611]  0.5536891817  0.9660097252  0.9410459516  1.5951513770 -0.2205657333
##  [616]  0.2178498025  0.5356140322  0.3110208741  1.0694432395  0.8540110968
##  [621] -0.0210732704  0.6835836831  0.8595792674  0.0400242481  0.6224204369
##  [626]  0.3944860811 -0.4860087972  1.4295780660  0.9663398415  0.3495761716
##  [631]  0.3986900792  1.5785480347  0.6959705541  0.7360113377  0.4000767913
##  [636]  0.1557061928  0.6043603253  0.1803701144  0.4147498885 -0.4274344631
##  [641]  0.0800214797  1.2931531494 -0.2820991360  0.6735457434  0.3595695730
##  [646]  0.4046196249  0.4249935328  0.5937490370  0.5079443441  0.3239413440
##  [651]  0.6340230086  0.1621165576  1.3116939210 -0.3864789732  0.5719734916
##  [656]  0.1587075531  0.6134085927  1.2900501354  0.3940700669 -0.2088683990
##  [661]  0.2025279744 -0.1864087973  1.0896097887  1.4588004445  0.9982047508
##  [666]  1.0905257104  0.5877630569  0.9516987303 -0.1771703883  0.9820636986
##  [671]  1.5274364187  0.5915652617  0.5779272910  0.5697104181  1.4492670468
##  [676]  0.7639874306  1.0798138988 -0.8083018177 -0.0208555891  0.3058978890
##  [681]  0.8692446840 -0.4744737730  0.4179975929  0.5206833782  1.0914936401
##  [686]  0.8974002323  0.4507732254  0.5770937920  0.0923025806  0.3543051720
##  [691]  0.0771773540  0.1954986157  0.9683608972  0.1894789545  0.9821801331
##  [696]  0.4263248154  0.1566852175 -0.0499826615  1.4271445451  0.5870345942
##  [701]  0.3550302620  0.5616242742  0.8514374312  0.0054632624  1.8572716402
##  [706]  0.2064176151  0.9950396435 -0.9970204836  0.8508982862  0.6296974232
##  [711]  0.2106084667  1.5939950156  0.4962384991  0.1544124327  0.1483448719
##  [716]  0.7075834349  0.5961339946  0.1352717191 -0.3398372310  0.0918674050
##  [721]  2.2484653149  1.1731754290  1.8208885998  0.1636997805 -1.2053217399
##  [726]  1.6780531050  0.5421625104  0.1606005897  0.5991453987 -0.7700062611
##  [731]  1.0488394397  0.4307311460  0.6085340131 -0.2561059428  0.1735702385
##  [736]  1.4104167644 -0.2304141403 -0.3669665963  0.5058338711 -0.2059694024
##  [741]  0.1406127841  0.6360781368  0.6325719589  0.3323553702  0.3107254555
##  [746]  0.0149476275  0.4983988590  0.2613827514  0.8925612693  0.7804590846
##  [751]  0.4315080028  0.3947296170  1.5990461000  0.6592607973  0.2777384519
##  [756]  1.4642491577  0.7544307200  0.8562434230  0.2010268511  1.5353846791
##  [761]  1.4695105376  0.0706735760  0.7782376942  1.2582879863  0.3053835510
##  [766]  0.8884596173  0.3746445255  0.2866243311  1.2735420966  0.4167006831
##  [771]  0.3522261452  0.3181410307  0.2191619620  1.5538809935 -0.6054838379
##  [776]  0.3178562392  0.3387033241  0.9518698242  0.4772595379  0.7571764941
##  [781]  1.4388592504  0.7260573186  0.8581583179 -0.0496264892 -0.2181648267
##  [786]  0.6224204369  1.2917788390  0.3297742505  0.6824038576 -0.0593022301
##  [791]  0.5752610483  0.0874930030  0.2598588757  0.5857560885  0.5887194630
##  [796]  0.4824334653 -0.3965813771  1.5390411046  0.6142473803  0.5974328485
##  [801]  0.3055822558  0.8043363393  0.7040661136  0.9419176777  0.6705022427
##  [806]  0.5469844567  0.6751579671 -0.3238143882  0.3967901154  0.7301998451
##  [811]  0.3846549119  0.5838062978  0.1086005663  1.4460934530  0.8677946875
##  [816]  0.5310559775 -0.5859441258  0.2892784974  0.1569783529  2.2045765978
##  [821]  0.4750031527  0.2552183239 -0.0379466717 -0.0850762700 -0.6186688034
##  [826]  0.8839298805  0.3058419623  0.2466673727 -0.0642113210  0.9402496669
##  [831]  0.3659315810  1.1694261709 -0.2943309965 -0.0110475385  1.6762165007
##  [836]  0.2803259968 -0.2408036943  0.6118779792  0.1384698861  1.5949711439
##  [841]  0.3837329321  1.0797782325  0.8402336152  0.5702218005  0.2022382798
##  [846] -0.2162407111  0.8474367628  0.4314402106  0.4822260188  1.4163278472
##  [851]  0.6587659826 -0.3159642173  0.6925618633 -0.3936942715  0.7638942975
##  [856]  0.8495177801  0.4094284916  0.3386673296  0.3322566943  1.0691149270
##  [861]  0.6652289394 -0.2012628792  0.7923687836 -0.3057559010  0.0857787306
##  [866]  0.4125749336  0.2795907322  0.9438920617  0.5524850550  0.1960088729
##  [871]  0.5589890733  0.1142275474  0.8551084220  1.0709824043  0.2501865574
##  [876]  1.4172543301  1.3382958445  0.4352510631  1.0561694673 -0.6952593100
##  [881]  0.5963906261  1.1151067236  0.4083415606  0.0325677715  1.4365603322
##  [886]  0.4493777327  0.4106051242  0.2791485936  0.2222358331 -0.0669477859
##  [891] -0.0587483221  1.2630625708  0.6912251508  0.8108660912 -0.2334733572
##  [896]  2.5070247824  0.0156047793  0.2577149254  0.1639372867  0.7000527640
##  [901] -1.2736105595  0.4884429282  1.0067456618  0.4959627906  0.9025281777
##  [906]  0.4657038742  0.4183381761  1.5430441639  0.7015356082  0.1544124327
##  [911]  0.5450533194  0.0314517967  2.6152671220  0.1911762371  0.9425716002
##  [916]  0.7175777209 -0.5171865882  0.2166021500  0.4302455910  0.1678880156
##  [921]  0.2808095821  0.4227633853  0.7687233164  0.5984458223  0.1576665191
##  [926]  1.2825606694  0.3619573230  0.0175068063 -1.3752345024  0.6223314384
##  [931]  0.6287430167  0.5872543432  0.1692819620  0.7235840001  0.7617719577
##  [936] -0.4021445620 -0.6078759061  0.8251033630  0.1540049928  0.6331541563
##  [941]  0.0191607221  0.6020938640  0.7972863735  0.6592022796  0.3719174163
##  [946]  0.6910236951 -0.0374404257  0.3603215874 -1.0896258277  0.1960109610
##  [951]  0.8518401485 -0.4758479702  0.9543114189  1.0784040468  0.8867663896
##  [956]  0.4467504146  0.8649051744  0.1306822994  0.8703721373  0.7349827300
##  [961]  0.7523893735  0.6793515299  0.4684784978  0.3706032497  0.3060897279
##  [966]  0.7913241770 -0.2460994882  0.4094263170  0.3713681345  0.8569014414
##  [971]  0.3996963287  0.8439019256  0.6179367654  0.5246710143  0.1488952534
##  [976]  0.4854189552  0.1015107833  0.0427698899  0.2227223892  0.0175064822
##  [981]  0.3882846557  1.3080395680  0.8668267772  1.2654106099 -0.2718740839
##  [986]  0.1033005777  0.7760009466 -0.2686105912  0.1789669402  1.0620974957
##  [991]  0.3691875793  1.8078777219 -0.2877049112 -0.3251382684  0.4942053312
##  [996] -0.6645695023 -1.1258637868  0.7030181845  0.2158274901  1.1045520109
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
##   0.60237360   0.36295619 
##  (0.11477683) (0.08115747)
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
## [1] -0.18851225 -0.09594645  0.80262225 -0.02420263  0.78630951 -0.34759057
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
## [1] -0.0089
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9113734
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
## t1*      4.5 -0.05565566   0.8968245
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 2 3 8 9 
## 1 1 1 4 3
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
## [1] 0.0354
```

```r
se.boot
```

```
## [1] 0.9082135
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

