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
## 0 2 3 4 6 8 
## 1 1 2 2 2 2
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
## [1] 0.0165
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
## [1] 2.759241
```

```r
UL.boot
```

```
## [1] 6.273759
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##   2.5%  97.5% 
## 2.7975 6.2000
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
##    [1] 4.6 5.5 4.8 4.1 4.6 4.5 5.9 4.5 3.5 6.0 4.9 5.4 4.5 3.1 5.3 4.5 3.9 6.4
##   [19] 5.3 5.6 4.8 4.1 5.5 4.6 6.2 4.3 5.9 6.2 3.5 4.0 4.9 4.8 4.2 4.6 3.4 5.0
##   [37] 4.5 3.8 6.0 3.9 4.6 4.0 5.6 4.8 5.4 5.2 4.2 4.2 4.6 5.3 5.0 4.7 5.0 3.9
##   [55] 7.0 4.3 4.3 3.8 4.4 2.6 5.0 2.9 5.6 4.0 4.9 4.1 2.9 4.2 5.2 5.1 4.8 4.7
##   [73] 5.3 3.5 5.1 4.5 4.4 4.5 5.4 4.4 4.8 5.1 3.5 3.4 4.9 7.3 5.6 4.2 4.2 4.0
##   [91] 4.8 5.8 5.0 4.3 5.5 5.6 5.7 5.0 4.6 5.4 5.7 4.2 4.2 4.8 3.7 4.7 4.5 5.5
##  [109] 3.5 5.1 5.6 5.2 4.0 4.5 3.6 4.4 3.9 4.1 5.1 3.4 6.3 5.1 3.5 5.2 3.1 4.9
##  [127] 5.2 4.8 4.2 4.9 5.6 5.2 4.1 4.6 4.4 3.6 3.5 3.6 4.3 5.4 4.9 4.7 4.9 4.5
##  [145] 4.9 5.5 4.4 3.6 4.3 4.0 4.6 4.8 5.8 4.1 3.1 6.0 4.5 4.4 3.0 4.1 4.6 5.4
##  [163] 6.5 3.5 4.9 5.7 5.0 5.2 4.6 2.8 4.2 4.4 3.2 5.3 2.5 4.1 4.1 2.6 5.1 3.8
##  [181] 4.5 4.1 4.9 4.4 4.9 3.6 5.3 4.0 3.5 3.9 2.0 4.2 6.5 3.7 4.0 4.3 5.8 4.9
##  [199] 4.8 4.7 3.4 4.3 3.6 4.0 4.8 4.7 3.0 4.3 5.4 5.3 6.4 4.3 4.5 5.9 5.6 3.6
##  [217] 4.3 4.7 6.8 5.0 5.2 4.8 4.7 4.6 4.9 3.5 5.5 4.6 2.7 4.4 4.4 3.9 3.5 4.5
##  [235] 3.5 6.3 5.6 5.0 4.2 5.9 4.3 4.2 4.3 5.3 4.6 4.1 6.4 3.0 5.7 4.4 5.7 4.2
##  [253] 5.0 3.8 4.9 4.6 5.6 3.3 5.9 5.9 5.4 6.2 4.0 3.0 5.4 5.0 3.6 5.0 4.3 4.6
##  [271] 5.2 4.6 5.3 5.1 6.3 3.9 4.3 4.5 4.7 4.4 5.3 6.1 3.9 4.1 5.9 2.9 5.8 4.9
##  [289] 4.4 4.4 4.1 4.2 2.7 5.3 3.8 5.1 4.1 4.7 4.5 3.2 3.0 5.1 5.9 3.9 4.0 4.4
##  [307] 4.1 4.9 4.8 6.2 4.6 5.5 5.3 3.4 5.1 4.1 6.5 4.2 4.4 5.3 5.7 4.4 4.8 4.6
##  [325] 4.3 4.8 4.3 4.0 4.8 5.4 4.6 5.2 5.2 4.6 5.1 5.8 5.6 5.7 5.5 3.6 5.7 2.8
##  [343] 6.5 3.4 4.1 5.0 4.1 5.7 4.3 4.7 4.2 5.1 5.2 2.7 3.8 4.6 5.1 4.8 3.9 4.3
##  [361] 4.7 4.4 4.4 3.7 4.7 4.9 4.9 6.1 5.3 3.7 3.1 4.8 3.1 4.6 4.7 5.5 5.2 5.6
##  [379] 5.0 3.5 4.9 4.7 3.8 6.1 5.1 3.5 3.8 3.0 3.7 5.4 5.6 4.4 5.3 3.8 5.4 4.1
##  [397] 4.8 5.6 4.0 3.2 6.0 4.5 3.8 6.4 5.3 4.9 4.6 4.3 5.2 4.8 5.2 4.6 4.2 5.0
##  [415] 2.4 4.0 5.0 4.7 4.4 3.3 5.2 2.2 3.9 4.2 5.9 4.9 5.3 5.5 3.4 3.9 4.2 4.1
##  [433] 3.1 3.7 5.0 5.2 4.0 3.8 3.5 4.3 4.2 4.2 5.5 3.8 5.9 4.6 5.0 3.9 3.1 3.6
##  [451] 3.2 4.7 3.9 6.0 5.6 3.0 5.0 3.1 5.6 4.3 3.9 3.4 3.2 3.8 3.8 4.7 4.0 4.6
##  [469] 3.7 4.5 5.4 3.7 4.6 3.1 4.4 3.7 4.0 5.0 5.2 3.7 4.7 4.1 5.3 3.9 4.3 5.8
##  [487] 4.3 5.7 4.4 4.9 4.2 3.0 4.1 4.1 2.1 3.4 4.2 5.7 4.6 4.7 4.5 4.6 6.8 4.7
##  [505] 4.2 4.0 5.2 5.0 6.3 4.4 3.4 3.9 4.4 3.9 4.9 5.7 6.3 6.4 3.1 5.5 4.8 4.8
##  [523] 4.3 3.1 4.7 4.9 5.1 4.4 3.2 4.6 5.0 4.8 2.6 3.8 4.6 3.0 4.0 5.6 3.7 3.2
##  [541] 5.0 4.9 4.3 3.5 3.9 3.8 3.4 5.4 5.1 6.0 3.9 3.0 3.6 5.5 4.8 5.5 4.1 3.9
##  [559] 4.7 4.2 4.5 3.9 5.4 4.2 4.7 5.2 5.5 5.4 4.6 4.3 4.8 5.4 3.1 4.8 3.6 3.2
##  [577] 4.0 4.4 3.5 5.0 3.3 4.6 5.5 3.6 4.7 3.9 3.7 4.9 4.5 4.1 5.3 4.4 3.7 3.8
##  [595] 5.9 6.2 5.7 4.1 4.8 4.2 4.2 4.0 4.4 3.4 4.1 4.5 4.1 3.9 2.7 4.2 4.5 3.6
##  [613] 4.0 3.9 5.1 4.9 4.0 4.2 5.3 5.2 3.9 3.8 4.3 5.7 6.2 2.4 3.4 5.0 4.9 4.8
##  [631] 3.4 4.5 3.8 4.6 5.5 2.3 5.7 3.0 4.6 5.0 4.2 5.8 5.3 4.4 3.9 3.5 3.6 5.7
##  [649] 4.8 4.8 5.3 6.1 4.5 5.8 5.3 4.5 3.3 4.1 4.3 3.3 4.6 6.8 4.1 3.8 5.7 3.1
##  [667] 4.8 3.3 4.2 4.8 6.1 5.8 4.5 4.2 5.3 4.2 4.5 4.3 5.8 5.3 4.4 5.7 4.7 4.1
##  [685] 3.4 4.8 4.8 4.5 3.3 5.2 3.9 4.3 4.5 5.7 4.4 3.8 3.6 5.0 5.1 5.0 5.9 4.9
##  [703] 4.6 3.7 5.3 4.2 5.0 6.8 4.8 3.4 5.8 3.8 4.9 2.8 6.0 7.0 4.7 4.8 5.8 4.2
##  [721] 4.7 5.0 4.7 4.4 4.5 4.6 4.2 4.6 4.0 4.8 3.2 4.0 3.8 3.2 5.6 4.4 4.8 5.4
##  [739] 5.9 5.0 5.1 4.5 4.8 3.0 3.9 5.2 3.8 4.3 5.4 5.1 5.3 4.1 5.6 4.9 5.0 4.1
##  [757] 6.3 3.8 5.0 6.2 3.5 3.7 5.1 3.1 4.6 4.9 3.9 4.4 4.6 4.0 4.1 3.0 4.2 4.4
##  [775] 5.3 4.4 5.6 5.1 4.0 5.3 4.5 6.1 4.8 4.7 4.7 5.3 5.0 5.1 3.9 4.9 4.2 5.7
##  [793] 3.7 3.5 3.7 4.6 4.4 5.2 3.6 4.3 2.4 4.2 3.3 4.6 4.2 5.0 3.6 4.6 5.5 4.6
##  [811] 3.5 3.5 5.3 4.0 3.9 5.5 3.7 4.0 3.6 3.7 5.5 3.7 3.0 5.5 3.9 2.5 5.0 4.3
##  [829] 3.8 3.5 3.8 4.2 3.8 3.6 3.1 4.6 3.5 3.9 4.1 3.8 3.6 3.3 5.9 5.1 2.8 4.0
##  [847] 4.9 5.0 3.5 5.5 4.4 4.9 5.8 3.2 4.9 4.0 4.1 4.7 4.3 3.7 3.8 6.1 5.5 4.7
##  [865] 4.0 4.8 4.2 3.5 3.9 4.2 5.0 3.2 4.7 5.1 4.3 4.9 5.3 5.5 4.7 3.9 2.7 5.5
##  [883] 4.4 4.1 4.8 4.9 3.3 3.0 4.0 5.0 3.9 3.8 5.3 5.5 4.3 4.7 5.3 5.0 2.7 6.0
##  [901] 5.2 5.3 4.4 3.7 5.0 4.1 4.0 3.3 5.2 4.4 5.2 4.2 5.4 3.8 5.0 4.5 4.2 3.8
##  [919] 4.6 5.0 4.4 3.1 2.5 4.5 3.6 4.1 3.8 4.1 5.4 3.4 4.9 6.1 5.6 5.3 5.0 3.5
##  [937] 5.4 4.2 3.7 4.6 6.4 3.5 4.2 3.6 4.8 3.7 4.7 6.1 5.3 5.5 2.7 5.2 5.1 5.5
##  [955] 5.0 5.5 5.3 4.5 5.3 4.2 4.9 4.7 4.6 3.1 4.7 4.2 4.1 5.7 5.3 4.7 4.9 4.3
##  [973] 3.9 5.2 4.9 5.3 4.0 4.2 5.0 3.8 3.7 5.2 3.3 3.1 4.3 7.0 2.5 6.4 3.6 4.4
##  [991] 5.0 4.0 4.3 4.2 4.9 5.1 3.7 4.9 4.8 3.9
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
## 2.8975 6.2025
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
##    [1] 4.2 5.1 5.7 4.3 5.5 3.8 5.6 4.2 3.9 2.3 3.3 6.1 5.1 3.9 3.2 4.2 5.7 5.3
##   [19] 5.3 4.7 4.4 4.3 5.0 3.8 4.9 4.6 5.6 3.6 4.5 4.6 6.0 4.8 4.7 5.0 4.6 5.7
##   [37] 3.1 3.1 6.3 3.0 4.8 5.4 5.7 4.6 5.0 3.8 3.8 4.8 5.1 4.1 3.5 5.7 3.7 5.4
##   [55] 4.4 5.3 6.5 5.7 4.7 4.9 6.1 4.3 4.5 6.4 5.0 5.9 4.7 3.2 4.6 3.5 5.8 4.5
##   [73] 5.3 6.9 4.8 4.4 5.4 4.8 5.6 4.7 3.3 3.2 3.4 2.1 5.0 4.3 3.0 4.5 3.7 4.1
##   [91] 5.5 4.3 4.4 4.5 4.5 3.7 5.3 3.8 5.3 3.6 6.3 3.9 3.7 5.0 4.7 5.0 4.2 4.0
##  [109] 3.3 4.4 4.9 4.8 5.3 5.4 4.5 4.4 6.5 3.7 3.4 4.3 3.4 4.2 3.6 4.4 4.2 4.9
##  [127] 3.9 4.4 4.5 2.6 4.6 4.0 5.6 3.3 3.8 5.2 4.6 5.9 5.5 3.9 4.6 6.5 4.3 3.9
##  [145] 4.9 5.5 4.0 3.6 4.5 4.9 5.3 2.5 3.4 4.5 5.7 4.5 5.1 6.2 2.7 4.8 5.8 5.6
##  [163] 6.0 4.3 4.3 3.4 5.5 6.3 5.9 4.5 4.9 5.1 4.1 4.1 4.6 3.3 6.2 5.8 3.5 5.9
##  [181] 5.7 3.5 3.3 4.6 5.5 3.0 5.3 5.7 6.0 4.9 5.2 3.8 5.1 3.2 3.6 5.0 4.0 6.0
##  [199] 6.4 4.5 4.9 4.7 3.2 4.1 5.4 4.2 3.5 4.5 3.6 6.3 4.6 3.7 5.7 6.2 3.6 4.2
##  [217] 5.1 2.9 2.4 3.5 4.1 5.7 5.2 4.1 5.5 5.0 5.7 4.6 6.1 3.0 5.6 2.8 5.0 3.8
##  [235] 5.0 4.4 4.7 3.8 5.3 4.4 4.5 5.1 4.2 4.4 4.5 6.0 4.9 4.3 4.0 4.8 4.6 6.2
##  [253] 4.0 4.7 4.1 5.9 3.7 5.4 5.9 4.4 4.1 6.3 4.5 5.4 5.0 4.5 3.8 5.6 3.7 5.0
##  [271] 4.4 3.9 4.8 4.4 4.3 4.6 4.6 5.1 4.5 3.8 4.8 2.7 4.3 3.5 3.9 4.5 4.6 4.6
##  [289] 3.9 6.0 3.5 4.8 3.1 3.7 4.5 3.4 4.8 4.5 3.9 4.6 3.9 4.3 5.5 2.6 3.5 4.7
##  [307] 3.6 4.5 4.7 4.3 5.3 5.4 5.2 4.1 4.4 6.0 4.4 3.2 6.6 3.4 5.3 3.4 5.7 5.1
##  [325] 3.8 5.6 4.4 6.3 5.2 2.4 3.7 6.4 5.2 3.8 3.2 3.9 3.2 3.2 2.7 5.6 6.0 5.2
##  [343] 6.2 4.0 5.3 4.8 4.8 3.2 2.8 4.3 3.0 4.5 4.8 6.7 3.9 4.3 4.9 3.7 4.7 5.8
##  [361] 3.5 4.0 5.0 4.4 3.7 6.8 5.9 5.1 6.4 4.2 4.6 6.4 3.4 4.3 4.6 6.0 6.1 5.8
##  [379] 3.9 4.6 5.2 4.7 4.0 5.2 4.3 4.8 4.6 4.1 2.9 5.0 2.8 5.1 3.8 5.4 4.2 2.4
##  [397] 2.9 4.9 3.4 3.2 3.3 4.0 5.4 3.9 4.8 3.8 4.6 6.1 4.0 5.1 3.8 5.8 3.6 4.7
##  [415] 5.0 5.2 4.7 5.0 4.2 3.2 2.3 5.1 4.5 3.7 5.1 4.3 3.2 5.1 4.5 4.5 4.8 4.8
##  [433] 6.4 3.4 4.7 5.2 3.7 5.2 4.8 5.2 5.0 3.5 5.0 3.9 4.8 5.6 5.3 5.5 5.0 5.0
##  [451] 3.6 5.5 4.3 4.2 5.3 4.6 4.2 3.9 4.3 3.6 6.6 3.3 4.5 5.4 4.6 4.4 3.5 4.3
##  [469] 4.5 4.6 3.7 3.1 5.6 5.5 5.0 3.1 3.4 3.3 4.2 4.9 4.8 4.2 3.5 3.2 4.7 3.5
##  [487] 3.1 3.5 5.3 4.3 4.8 4.1 3.0 3.9 4.6 5.1 4.6 5.4 3.6 5.5 4.6 3.5 6.2 5.3
##  [505] 4.1 4.7 3.3 4.2 5.1 2.8 7.0 4.5 5.1 3.0 5.8 6.5 5.5 5.0 4.4 5.7 4.6 4.4
##  [523] 4.9 4.9 4.2 4.8 4.3 4.5 4.9 3.8 5.0 3.7 3.2 5.1 4.7 6.0 3.5 2.2 5.6 6.2
##  [541] 4.7 3.3 4.0 4.8 4.9 4.3 4.8 2.9 4.6 4.5 3.8 3.8 3.0 5.3 4.5 4.1 4.9 6.2
##  [559] 5.8 6.0 4.0 5.0 3.8 5.0 5.2 4.3 3.9 3.1 6.1 5.8 4.5 4.9 5.5 3.7 4.9 2.5
##  [577] 4.3 5.7 3.7 4.9 5.0 3.9 3.8 5.1 4.3 3.9 5.0 2.6 4.0 4.7 4.8 5.6 5.6 3.4
##  [595] 5.7 3.5 4.0 3.9 5.4 2.9 5.1 3.7 4.8 4.2 4.9 4.5 5.3 4.3 5.6 4.1 5.7 5.3
##  [613] 3.8 5.7 3.5 3.6 4.3 4.7 3.0 4.3 3.4 6.1 4.3 4.3 6.4 3.6 3.8 6.3 5.7 5.3
##  [631] 4.4 5.2 5.1 4.7 5.5 4.1 3.9 4.8 4.9 4.4 3.7 4.9 3.4 5.7 5.1 4.3 3.9 6.6
##  [649] 3.7 4.5 4.8 3.7 4.0 4.8 5.0 3.6 3.9 4.9 4.7 4.4 4.2 4.9 4.5 6.3 4.0 5.6
##  [667] 4.9 3.6 4.5 4.4 4.2 5.6 4.0 4.4 6.2 4.6 5.2 5.3 4.9 3.9 5.9 6.0 4.9 5.2
##  [685] 4.4 5.4 4.0 5.7 5.3 3.2 6.1 4.5 4.9 5.7 4.4 3.9 5.1 3.6 4.9 6.5 4.1 4.0
##  [703] 5.2 4.3 4.1 4.1 5.0 3.3 3.8 6.1 7.0 5.0 3.5 4.8 3.3 5.3 4.6 4.3 5.3 3.9
##  [721] 4.6 4.8 3.2 5.8 4.6 3.8 4.1 4.9 5.0 4.9 4.3 4.3 5.6 4.1 5.3 4.8 5.6 4.3
##  [739] 3.0 4.0 5.0 4.0 5.4 5.1 4.9 3.7 4.5 3.8 4.6 4.1 5.5 3.9 4.6 3.9 5.9 4.1
##  [757] 3.4 5.2 3.1 4.1 4.2 6.1 3.9 4.0 6.5 5.0 4.8 4.6 5.0 5.1 2.7 4.8 4.1 3.1
##  [775] 5.9 4.5 6.0 3.8 5.9 5.4 4.3 3.5 5.3 3.4 4.0 3.3 3.8 3.6 3.3 4.3 2.6 3.7
##  [793] 4.0 3.8 3.6 3.6 4.7 2.8 5.2 4.4 3.8 4.5 4.8 5.3 4.2 4.8 5.1 4.1 4.3 4.8
##  [811] 4.0 4.6 2.8 3.6 4.7 4.6 5.3 3.8 4.2 3.8 6.4 4.7 4.6 4.5 5.2 2.5 4.6 2.9
##  [829] 6.0 5.9 4.6 3.8 5.1 3.5 3.3 4.3 4.2 4.3 3.2 5.0 3.3 5.0 3.1 3.8 4.5 5.0
##  [847] 2.6 5.1 4.7 5.0 3.7 5.2 3.9 4.3 5.9 4.3 4.8 4.9 4.7 4.1 3.9 3.3 6.4 4.2
##  [865] 4.0 4.4 4.8 4.4 5.6 5.2 4.9 3.9 4.3 6.8 4.0 5.8 4.6 4.6 4.2 5.3 4.2 4.4
##  [883] 2.4 4.2 3.8 3.1 4.5 4.8 3.7 3.7 5.3 3.2 5.8 4.7 5.7 4.4 3.7 5.1 4.6 6.0
##  [901] 5.1 5.7 5.6 4.8 3.7 5.6 4.5 5.4 4.6 4.2 4.3 5.7 4.3 4.6 5.4 5.3 6.4 4.0
##  [919] 5.3 5.1 5.2 3.7 5.7 4.1 5.2 4.5 3.5 4.2 5.1 6.4 3.5 4.2 3.0 5.0 5.0 4.6
##  [937] 4.7 3.5 3.2 3.6 3.4 4.8 5.9 4.3 4.3 5.0 5.0 6.4 5.4 3.0 5.0 3.6 5.7 3.9
##  [955] 4.5 4.3 3.6 3.7 6.1 4.1 4.2 4.2 5.1 5.1 5.9 4.2 4.2 4.6 4.5 4.0 5.1 5.7
##  [973] 5.4 4.8 5.8 5.3 4.6 3.6 3.7 4.7 5.6 5.2 5.4 4.5 4.5 4.7 3.9 5.0 3.7 4.6
##  [991] 4.7 5.0 4.7 5.3 5.0 5.5 4.6 3.0 5.8 4.4
## 
## $func.thetastar
## [1] 0.0577
## 
## $jack.boot.val
##  [1]  0.57060440  0.40948509  0.33473054  0.21383812  0.11488764  0.04073034
##  [7] -0.16242236 -0.19128065 -0.35820896 -0.44068323
## 
## $jack.boot.se
## [1] 0.9595841
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
##    [1] 5.2 5.8 3.8 5.5 5.4 3.2 4.2 3.1 3.7 3.0 6.0 4.2 4.0 4.0 4.9 5.8 4.9 3.5
##   [19] 4.2 4.0 5.7 4.2 4.8 4.6 3.7 7.0 3.2 3.9 4.0 4.6 3.8 5.1 4.3 4.6 5.3 5.6
##   [37] 5.3 4.4 6.4 5.0 4.8 4.2 5.3 3.8 3.1 4.1 4.5 4.5 6.8 3.9 3.6 2.1 4.7 4.7
##   [55] 4.0 3.9 5.3 4.3 5.2 3.8 3.5 2.9 4.7 5.3 3.9 5.2 5.1 4.9 4.4 3.7 3.9 4.1
##   [73] 3.2 5.1 4.3 2.6 5.1 4.0 4.2 4.1 4.8 4.6 3.8 4.3 3.7 3.9 4.1 3.0 5.1 3.2
##   [91] 2.8 6.0 4.7 2.7 4.0 4.8 3.0 6.1 5.0 5.0 5.1 4.8 4.9 5.9 5.2 5.4 2.7 4.6
##  [109] 5.5 5.4 3.8 5.1 5.3 4.7 3.2 4.5 6.0 5.0 3.5 5.2 5.7 4.9 3.4 2.5 4.8 3.4
##  [127] 4.6 6.2 4.4 4.3 4.5 3.9 4.0 4.2 4.5 3.4 4.0 5.2 4.9 3.5 3.6 4.3 3.2 4.2
##  [145] 5.0 3.3 4.5 5.2 4.6 6.4 3.5 4.6 5.3 3.4 3.2 4.4 3.5 3.0 3.0 4.7 3.3 5.2
##  [163] 5.3 4.5 5.5 3.8 4.9 4.6 4.9 4.7 3.1 3.5 4.4 4.5 3.3 4.7 3.4 2.1 5.4 4.7
##  [181] 4.8 3.7 5.4 4.8 4.6 4.8 5.9 5.9 5.0 4.7 6.3 4.6 4.3 4.2 2.2 5.4 4.0 5.8
##  [199] 4.6 3.9 3.1 5.5 5.8 4.3 3.9 5.1 4.4 5.4 3.3 3.7 4.7 3.4 3.1 3.8 5.9 4.8
##  [217] 4.1 5.7 4.7 5.1 4.3 6.2 4.9 4.8 5.0 5.6 3.0 3.6 2.8 4.2 4.5 4.8 3.3 4.7
##  [235] 5.0 4.2 3.1 5.2 5.1 4.2 4.7 4.8 3.4 4.5 5.9 3.4 4.2 3.1 3.0 4.7 3.5 4.6
##  [253] 5.3 5.4 4.9 3.5 4.5 5.7 3.9 5.5 5.1 3.8 4.0 3.4 3.6 4.9 5.1 4.6 5.2 5.0
##  [271] 5.1 4.2 4.3 3.2 5.1 3.4 3.1 3.5 4.1 4.4 2.9 3.5 3.9 4.9 4.6 5.0 4.4 5.3
##  [289] 3.2 4.0 5.5 5.3 4.4 4.7 5.9 3.9 5.7 5.5 5.1 4.7 6.5 5.4 4.4 4.8 4.4 6.5
##  [307] 4.9 5.1 4.4 3.8 6.9 3.0 3.7 4.7 3.1 4.7 3.8 3.5 3.8 4.2 4.8 3.6 4.9 5.8
##  [325] 6.0 4.6 3.3 5.0 3.7 4.2 3.9 3.6 4.2 4.3 3.9 5.8 4.0 4.2 6.5 4.9 4.6 3.9
##  [343] 5.9 6.2 3.2 4.4 4.3 3.9 3.6 4.8 4.3 4.7 5.2 4.8 4.9 3.7 4.8 4.4 4.5 5.3
##  [361] 6.2 3.8 5.5 4.4 5.2 5.5 6.9 5.0 3.4 2.8 5.2 4.2 5.4 2.8 2.3 4.8 5.2 4.2
##  [379] 5.6 4.5 5.5 7.0 3.8 7.3 2.4 3.7 4.4 4.1 3.7 3.7 4.6 5.2 3.8 4.1 4.2 2.7
##  [397] 5.3 3.8 4.9 3.5 4.9 4.1 5.3 4.5 4.5 4.8 4.2 4.8 3.9 5.2 4.4 3.6 4.6 4.4
##  [415] 4.1 4.8 5.0 3.7 5.7 4.5 2.7 3.7 3.8 3.8 5.0 3.5 6.0 3.6 6.1 4.9 3.4 3.5
##  [433] 5.0 4.8 3.9 6.3 4.0 3.9 4.3 5.0 4.1 3.2 5.5 3.5 4.2 5.0 4.1 3.0 4.6 5.3
##  [451] 3.9 4.5 3.7 3.6 3.6 4.9 5.4 4.1 4.8 4.1 5.8 4.9 3.2 3.9 3.6 3.4 6.0 6.2
##  [469] 5.4 7.0 4.5 2.8 3.7 5.8 5.3 6.8 5.1 5.2 4.6 3.5 3.0 5.9 4.7 5.5 4.2 5.0
##  [487] 4.4 5.9 5.8 6.1 3.8 5.0 4.3 3.9 6.4 4.8 4.4 3.6 3.9 3.9 3.4 3.8 4.5 4.9
##  [505] 5.6 4.8 5.8 5.4 3.4 4.8 2.7 4.6 5.4 4.5 4.7 4.3 5.4 3.4 4.5 4.0 2.8 4.9
##  [523] 2.5 3.2 4.4 5.4 3.4 4.2 4.8 5.6 5.4 4.5 4.9 3.5 6.3 4.9 4.7 4.3 5.2 5.6
##  [541] 4.9 4.5 4.6 3.5 5.1 4.7 4.6 5.1 3.5 3.6 3.5 2.9 5.8 4.6 5.0 4.7 4.4 3.7
##  [559] 6.3 3.7 3.9 3.7 2.7 5.0 4.0 4.0 2.7 5.5 3.4 4.7 3.7 4.3 4.2 4.2 3.2 3.6
##  [577] 3.0 3.3 4.3 4.0 5.4 4.5 3.4 6.2 4.1 4.7 4.5 4.6 4.5 4.6 4.0 4.7 4.1 4.2
##  [595] 4.4 4.6 5.5 4.5 3.8 5.9 3.9 4.4 4.3 5.9 4.7 6.4 4.0 6.4 4.2 5.0 5.3 5.1
##  [613] 3.8 3.9 5.5 4.2 3.6 5.7 4.9 5.0 5.3 3.9 4.4 4.1 4.6 4.9 4.2 3.3 4.8 5.4
##  [631] 5.7 3.4 5.5 4.6 5.2 4.8 4.9 4.6 3.9 3.6 3.5 4.2 4.1 2.7 3.3 5.4 5.1 4.7
##  [649] 5.3 4.6 4.1 3.3 4.0 3.2 5.6 5.7 5.1 5.2 3.7 3.6 5.7 4.6 4.3 4.2 4.4 5.1
##  [667] 4.5 4.6 5.3 4.3 5.3 4.9 4.2 4.4 5.5 5.0 5.6 3.1 3.5 4.7 4.1 3.3 5.1 3.3
##  [685] 5.0 4.7 5.1 4.2 4.5 3.5 4.2 3.7 5.5 4.1 3.6 5.2 3.2 4.6 4.3 5.0 4.4 4.7
##  [703] 2.8 4.0 4.7 5.2 4.7 3.7 3.6 4.1 3.3 5.5 4.1 5.4 4.3 3.5 4.9 4.5 2.6 4.0
##  [721] 4.3 3.8 5.1 5.7 5.0 3.9 4.3 5.0 4.0 4.5 3.9 3.5 3.9 3.5 4.8 4.3 6.0 4.3
##  [739] 4.9 3.6 4.1 4.7 4.6 4.9 4.5 5.4 4.5 4.7 3.6 3.1 3.6 4.6 4.5 4.2 4.5 4.1
##  [757] 6.7 4.8 4.0 5.2 4.2 2.8 4.6 4.4 3.6 4.4 2.6 5.5 6.2 5.9 4.5 5.9 3.1 3.3
##  [775] 3.3 5.6 4.3 4.8 5.5 4.8 5.5 3.9 4.9 5.1 4.5 4.6 4.5 3.8 5.0 4.9 5.3 5.1
##  [793] 6.0 4.1 4.3 3.6 6.2 4.8 4.4 5.0 5.0 5.3 4.3 5.3 4.5 3.3 3.5 5.9 3.8 5.9
##  [811] 4.7 5.5 4.3 2.9 4.1 5.4 4.4 5.7 3.3 4.0 3.9 4.2 5.5 5.6 4.7 5.6 3.9 5.1
##  [829] 4.7 2.6 5.5 4.5 2.6 4.9 4.2 4.9 3.9 4.4 3.5 4.8 4.8 4.2 4.6 5.2 5.3 3.6
##  [847] 4.5 5.1 6.4 2.7 5.0 5.3 3.2 4.9 5.0 4.6 2.7 4.8 3.7 3.9 3.9 5.0 5.7 5.7
##  [865] 5.6 2.9 5.7 3.6 4.3 4.5 4.1 4.2 4.2 3.2 3.8 3.4 4.7 4.1 4.7 4.8 4.4 5.1
##  [883] 2.7 4.9 5.7 3.6 3.3 4.9 3.1 3.9 3.4 3.8 3.7 4.4 4.9 2.4 4.1 4.2 5.9 4.2
##  [901] 3.9 3.9 3.3 5.6 5.5 4.7 5.9 5.5 3.3 5.3 6.4 4.9 4.1 5.7 3.0 6.1 6.0 4.0
##  [919] 3.9 4.9 4.6 5.0 4.8 5.0 4.3 5.3 4.4 4.3 2.5 3.4 4.6 4.5 3.8 4.6 4.2 4.4
##  [937] 5.2 3.3 3.4 5.1 4.4 4.0 6.0 5.2 3.9 3.4 4.2 6.2 4.2 2.1 7.1 4.3 2.4 3.4
##  [955] 4.2 4.7 5.7 3.4 6.2 4.8 3.3 4.3 4.1 4.6 3.3 3.7 5.0 4.4 4.9 5.3 4.3 4.7
##  [973] 5.0 3.8 4.2 4.6 4.1 5.7 4.5 5.5 4.2 6.0 3.5 5.4 4.8 2.6 4.1 3.8 4.4 5.1
##  [991] 4.6 5.1 2.9 2.7 4.6 5.0 3.6 4.6 5.8 4.0
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.46 5.30 5.20 5.10 5.00 4.90 4.90 4.70 4.60 4.40
## 
## $jack.boot.se
## [1] 0.9305998
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
## [1] -0.8206239
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
##   1.9300488   3.3473373 
##  (0.7997642) (1.5825468)
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
## [1]  1.21870601  0.64246826  0.52538683 -0.45731204  1.44340945  0.09656041
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
##    [1]  0.619269375 -0.749185665  0.424592358 -0.961460830 -0.828729721
##    [6] -1.459529758 -0.918070909 -0.933703551 -0.431052881 -0.981741400
##   [11] -1.295167907  0.493470372 -1.421429764 -1.140228006 -0.954336666
##   [16] -0.803604344 -1.059778028  0.114921247 -1.668977542 -0.477199223
##   [21] -0.534161263 -0.024005407 -1.460215021 -1.141716948 -0.635615890
##   [26]  1.013620495 -0.744382812 -0.602907503  0.137726603 -0.793201408
##   [31] -0.643695409 -0.062371474 -0.986261398  0.026754474 -0.020774416
##   [36]  0.495012890  0.430120662 -1.184363885 -0.481098898 -0.238290179
##   [41] -1.684655609 -0.587198076 -0.888081974 -0.836607316  0.653926364
##   [46] -0.953885494 -0.243918839 -0.810320529 -0.104954133 -0.757529448
##   [51] -0.961460830 -0.411177853 -0.694676395 -1.400179720  0.365152956
##   [56]  0.084877087 -0.838023224 -0.273395477  0.264212604 -0.933145728
##   [61] -0.391783455 -1.091044143 -1.838911403 -0.421623563 -0.639680822
##   [66] -0.757833884  0.073694643  0.062699271 -1.308998424 -1.370361314
##   [71] -0.976304792  0.024857157 -0.673819240 -0.522588444 -0.732941562
##   [76]  0.606273115 -0.245092083 -0.403501580 -0.103421390 -0.353231861
##   [81] -1.191245495  0.382671138  0.137746229 -0.318451613  0.303230230
##   [86] -0.244277156 -0.913679641 -0.213987032 -1.150850645 -0.454306260
##   [91] -1.284025265 -0.864491866 -0.411964692 -0.744024651 -0.569662057
##   [96] -1.581413862  0.155809276 -1.541559532 -1.387954862 -0.192280793
##  [101] -1.149986679 -1.286850642 -0.405792478 -1.817940643 -0.710410999
##  [106]  0.298941944 -1.161839295 -1.105388146 -0.701535169 -0.783340807
##  [111] -0.509352392 -0.629750805 -0.259335207  0.764873881 -1.505242541
##  [116] -0.169006146 -0.663433616 -1.293437444 -0.597973999 -0.989623957
##  [121]  0.141437560  0.522138836 -0.987611873 -1.909423737 -0.627476835
##  [126] -1.153800416 -1.267342808 -0.890683227 -0.506470943 -0.701641879
##  [131] -0.465015837 -0.237151118 -0.447026018 -0.935262178 -0.227286541
##  [136] -1.340392349 -0.529349379 -0.312521218 -0.229220495 -1.165017247
##  [141] -0.534975373 -1.043422044 -0.398319792 -1.202008351  0.036067191
##  [146] -0.947621200 -0.064864856 -0.665037471  0.672935025 -0.873977347
##  [151] -0.692242184  1.628543433 -0.285277382 -0.098074539 -0.994949554
##  [156] -0.283962102 -0.784219117 -0.662888104 -0.739304274  0.039823924
##  [161] -0.011726268 -1.170800119 -0.429570975 -1.488355400 -1.211579330
##  [166] -0.948616315  0.310500870 -0.443843336 -1.046437636 -1.277689793
##  [171] -1.854326985 -1.251604505 -0.825866142  0.305571600  0.755879217
##  [176] -1.018350472 -0.644134528 -1.277230296 -0.842252666 -0.766269070
##  [181] -1.363294598 -0.863354537  0.263752680 -0.575663412  0.149945012
##  [186]  0.484910848 -2.472584288 -0.504806791  0.220935599 -0.649043103
##  [191] -1.165534859 -0.634140730 -0.937039378 -1.281773200 -0.201152968
##  [196]  0.561636614 -0.409311490  0.121551180  1.195061400 -1.041722540
##  [201] -1.405175525  0.490461159 -0.727265326  0.420153301 -0.398053716
##  [206]  0.503846703 -0.639470158 -1.375063635 -1.432176658 -1.227544992
##  [211] -0.647275605 -0.504987962 -0.486244527 -0.251712993 -1.536624229
##  [216] -0.763651850 -0.506274420  0.421054106 -0.451712606  1.585289019
##  [221] -0.392087082 -0.557691997 -0.046006221 -2.076683123 -0.527918590
##  [226] -0.143211548 -0.795847931 -1.160513593  0.324844092 -0.780889752
##  [231]  0.129744095 -0.294010497 -0.942464287 -1.029150330 -0.711542599
##  [236] -0.172804644  0.356472455 -0.172852110 -0.783440035 -0.114769225
##  [241] -1.077992920 -0.559199521 -0.500965781 -0.831815741 -1.523884257
##  [246] -0.651732818 -0.490687449 -0.205793545 -1.141716948 -0.538531796
##  [251] -0.001724240 -0.523072599 -0.704311174 -0.368116634 -0.237773896
##  [256] -0.545764130 -0.545853435 -0.040028357  0.981593225 -0.684600623
##  [261]  0.417740384 -0.107049922 -0.358782611 -0.776161333 -1.166086307
##  [266] -0.506370820 -1.048826351 -1.439831142 -0.014233332  0.063009347
##  [271] -0.969152005 -0.460577475 -0.590117680 -1.188412161 -1.305614171
##  [276] -0.622466460 -0.619479410 -0.169702495 -0.794141632 -0.736902106
##  [281] -0.742300503 -0.740104206 -0.756397499 -0.444124355 -0.366900348
##  [286] -0.649043103  0.411662931 -0.690873061 -0.460196320 -0.203582419
##  [291] -0.995955028 -0.507548209 -0.752943143 -0.634052963  1.047286454
##  [296] -0.984932228 -1.923102630 -0.679554376 -0.566624824 -0.761522981
##  [301] -0.543798718 -0.743410285 -1.002295185  0.036726189 -1.116530171
##  [306] -0.946179989 -0.871821046  0.181224591 -0.730633861 -0.512919923
##  [311] -1.167081124 -1.320402432 -1.037270628 -1.234835502 -1.420016909
##  [316] -1.471006886 -0.804932328 -0.538604200 -0.558986766 -0.584728924
##  [321] -1.143593362 -1.101933821  0.008444280 -0.049057114 -0.972897469
##  [326] -0.511516078  0.219065458 -0.280918381 -0.316165286  0.395410683
##  [331] -1.601955314 -0.959013488  0.017250964 -0.367892353  0.613111452
##  [336]  0.617040838 -0.534379941 -0.157596424 -1.490317136 -0.781666050
##  [341]  0.279556834  0.408810585  0.017495939 -0.387235935 -0.901152148
##  [346] -0.979759800 -0.226616477 -1.674628481 -0.588563894  0.384447448
##  [351] -0.249444782 -0.915049539  0.743557845 -0.568729281  0.311088427
##  [356]  0.020596081 -0.411650425 -0.853178507 -0.096214184 -0.892434375
##  [361] -0.903531103 -0.010162760 -1.092707221 -1.060338480 -0.096628981
##  [366] -0.409024079 -0.607082351 -0.461436670 -0.779266120 -0.917711505
##  [371] -0.720270775 -1.706774609 -0.947396554  0.469034535 -1.157641335
##  [376] -1.634081919  0.012150779 -0.720326137 -0.379551796 -0.220791555
##  [381] -1.020010487 -1.065110521 -1.305336882 -1.343149386  0.356823027
##  [386] -0.629393194 -0.493282497  0.424415125 -0.562021869  0.029486617
##  [391]  0.006796312 -0.834510205 -0.984851038 -0.749003256 -1.169627494
##  [396]  0.266236389  0.107349316  0.013702572 -0.442292937 -0.758089661
##  [401] -0.345420546 -1.237007758 -0.867684164 -1.110595867 -0.430630107
##  [406] -0.623751507  0.258726054 -0.706226508 -0.093197776  0.048526114
##  [411] -0.940300447  0.025213559 -0.702417073 -0.293724521 -1.118809921
##  [416] -0.299281531  0.307397837 -0.764927242 -1.120508508 -0.394919034
##  [421] -0.301770791 -0.865659060 -0.460577475 -0.245008717 -0.214256693
##  [426] -1.603099929 -0.783824331 -0.639525237 -0.248674029 -0.812239380
##  [431] -0.227900602  0.677415371 -0.389981470  0.432516961 -0.480906033
##  [436] -1.060659993 -1.255643552 -0.494715656 -0.979352098 -1.275403257
##  [441] -0.710497219 -0.147837233 -0.949323427 -0.872970112 -0.502831969
##  [446]  0.221787131 -0.520732122 -0.831387912 -1.216078118 -0.603911221
##  [451] -0.345873179 -0.516594062  0.054668712  0.473809394 -0.861829638
##  [456]  0.179477910 -0.602907503  0.096044393 -0.520462238 -0.162938068
##  [461] -0.700199662 -0.839163748 -1.061974234 -0.500205869 -1.094247986
##  [466] -0.214702046 -0.955664139 -0.742301149  0.701770251 -1.306313659
##  [471] -0.156930276 -0.313270044 -0.493157499  0.087646193 -0.266343109
##  [476] -0.329056292 -0.067676475  0.355764500 -1.183652531 -0.781174531
##  [481] -0.225921250 -0.622160399 -1.598892290 -0.840093620 -0.376273932
##  [486] -1.032777517 -0.890561183 -0.094293630 -0.606778550 -0.951194223
##  [491] -0.759701340  0.193347948  0.030281447 -0.947411039  0.686460419
##  [496] -0.980083301 -1.182902987 -0.484935316 -0.972451007 -0.399566030
##  [501] -0.122950601 -1.424549698 -1.340007339 -0.720152846 -0.058763967
##  [506] -1.887023237  0.201950025 -0.700479208  0.223909144 -0.538087050
##  [511] -1.196044604  0.099592477 -0.753490561 -0.241796803 -0.992006356
##  [516] -0.471565247  0.347170093 -0.231648108 -0.836282440 -0.643089868
##  [521] -0.734866454 -0.331596226 -0.555191325 -0.220576459 -0.258135829
##  [526] -0.682349934 -0.585328011  0.360501922  0.044367497 -0.953869481
##  [531] -1.631689133  0.159004025 -0.312384753 -0.960180617  0.817902981
##  [536] -0.475507777 -0.891875367 -0.508774605 -1.025168146  0.552820320
##  [541] -1.002201052 -2.414893490 -0.381027733 -0.606609960 -0.564571019
##  [546] -0.095512664  0.691795693  0.208713986  0.102165888 -0.775887596
##  [551] -0.540042919 -0.688695733  0.482629490  0.348597585 -1.228180288
##  [556] -0.447582557  0.859024477 -0.248095902 -1.605071401 -0.277450120
##  [561] -0.756266911 -0.748762514 -0.517884046 -0.665417775 -0.156547852
##  [566] -1.237420032 -0.560639857 -1.280460979 -1.184694131  0.801235746
##  [571] -0.688548040 -0.505208161  0.100555903 -0.141696567 -1.253301874
##  [576]  0.826233244  0.523250977 -0.334730961 -0.122298313 -0.761972567
##  [581] -0.941053313 -1.000778952  0.859787850 -1.081759197 -1.242927136
##  [586] -2.081092812  0.260139047 -0.407599648 -0.695812371  0.303789669
##  [591]  0.263069576 -0.317896777 -0.452050669 -0.513333262 -0.745327223
##  [596]  0.550325599 -0.313936880 -0.365494058  0.011485033 -0.504964953
##  [601] -1.338849108  0.039986763 -0.779183942 -0.711863285 -0.096435860
##  [606] -1.312727888 -0.613090722 -1.962058247 -1.281242246 -1.047874175
##  [611] -0.774226851 -0.399592759  0.033420824 -0.489519439 -1.307938390
##  [616] -0.777237632 -0.656944252 -1.207384858 -1.808644612 -0.764606769
##  [621] -0.337378769 -0.618677565 -1.471944308 -0.738109774 -1.046054926
##  [626]  0.559480381 -0.912387457 -0.677911363 -1.536466985 -0.905491260
##  [631] -1.178531590 -0.516610168 -1.548988863 -0.256367263 -0.588916787
##  [636] -0.390054289 -1.131851768 -1.141261363 -0.868111586  0.191036370
##  [641] -1.562526944 -0.883018527 -0.914233416 -0.426010739  0.420986410
##  [646] -0.481536959 -1.006503286 -1.479233755 -0.071631985  0.248183071
##  [651] -0.087622048 -1.104857965 -1.266321467 -0.169066511 -1.283138930
##  [656]  0.002911777 -1.445384676 -0.491846517 -0.008996967 -0.806288736
##  [661] -0.406541231 -0.794537861 -1.218760555 -0.284361235 -0.529241202
##  [666]  0.265632568 -0.332385254 -0.732816921 -0.716801612 -0.619143429
##  [671] -0.610105185 -0.025929620 -0.649857617 -0.802499759  0.373726797
##  [676] -1.862973023 -0.785105028 -1.781970242  0.861557124  0.189784611
##  [681] -0.128005254  0.385692458 -2.496214232 -0.168558148  0.236344170
##  [686] -1.132534622 -0.315008958 -1.105435306 -1.458588341 -1.196185687
##  [691] -0.824970053 -0.341860012 -0.956800005 -1.081618773 -0.928651324
##  [696] -1.286535691 -0.232687634  0.038605026 -0.039172705 -0.953076061
##  [701] -1.322853156 -0.646821154 -0.703121672 -0.287431579 -0.768890082
##  [706] -0.428026659 -1.982089710 -1.262667075 -0.543747544 -0.820623922
##  [711] -1.301788819 -0.777320750 -0.249285043  0.644899244 -0.182721320
##  [716] -0.998482008 -0.637524756 -0.821455236 -0.916402479 -1.833105442
##  [721] -0.910149215 -0.169452812  0.509983582  0.018294521  0.246829030
##  [726] -0.944069718 -0.614289253 -0.281816703  0.416326650 -0.056375552
##  [731] -0.130143690 -1.325549341 -0.026825416 -1.019093952 -0.216015394
##  [736]  0.078757500 -2.539887857 -0.731800336 -0.861094275 -0.639103200
##  [741] -0.930208364 -1.458484521 -1.162659490  0.040519435  1.058663957
##  [746] -0.885071634 -0.218608953 -0.089908476 -0.883209087 -1.572541115
##  [751]  0.601982826 -0.563701672 -0.052479783  0.005776841 -0.218360001
##  [756]  0.812060935  0.178167836 -1.010781656 -0.689003955 -0.547135370
##  [761]  0.003121750 -0.346614205 -0.476162068 -0.719364541 -0.786591558
##  [766] -1.277689793 -0.647779385 -0.542434674 -1.808232621 -0.724715024
##  [771] -0.527461109 -0.383548857 -0.561264608 -0.931525336 -1.163099982
##  [776]  0.622388199 -1.207293366 -1.025683183 -0.190247786 -0.590857863
##  [781]  0.225703701 -1.354377822 -0.876208829 -0.156862449 -0.959878011
##  [786] -0.649803888 -1.482645674  0.270894083 -0.379861439 -0.350243214
##  [791] -0.922883935 -0.070721566 -1.080300019 -0.036350355 -0.931951136
##  [796]  0.133020067  0.437515702 -0.700757918  0.322975764 -0.755466050
##  [801] -0.314595419 -0.648828784 -0.496221395 -0.197354073 -0.288648166
##  [806] -0.686440011 -0.209020760 -0.467467131 -0.041457879 -0.189527857
##  [811]  1.090548762 -0.327084192 -1.999872903 -0.293255315 -0.559645843
##  [816] -0.842040512 -1.547885412 -0.264505901  0.042885515 -1.181073270
##  [821] -0.460196320 -0.744146922 -0.221144654 -1.263443913 -0.647448267
##  [826] -0.770596746 -1.522898717 -1.297373305 -0.752182638  0.022199929
##  [831] -1.999908645 -0.320771282  0.164104653 -0.156311000 -0.324384403
##  [836] -1.014803141 -1.515743118 -0.650113668 -1.153800416 -0.528864700
##  [841] -0.653082795 -0.042077956 -0.364934364 -0.229205803 -0.545662017
##  [846] -2.341491484 -1.574013803 -0.695434967 -0.485354489 -0.368246719
##  [851] -1.111512270 -1.008534696 -0.659973849 -0.459024880 -0.070765831
##  [856] -0.802235547 -0.693535635 -0.561600206 -0.218892146 -0.643268732
##  [861] -0.446807934 -0.558221402  0.194288707 -0.914277774 -0.390513216
##  [866]  0.083265569 -0.381872253 -0.899354354 -2.042012821 -0.534032971
##  [871] -0.734766663 -0.809277508 -1.138793157 -1.335284018  0.440205964
##  [876] -0.898129522 -0.332551239  0.027728550 -0.782273722 -0.126622673
##  [881] -0.577832632 -0.539018695 -0.599352624  0.190517553 -0.089424750
##  [886]  0.707406209 -0.697199754  0.803796783  1.262369663 -0.283634258
##  [891] -0.794537861  0.168627531  0.833293374 -0.167654541 -0.992877708
##  [896] -0.324806493 -1.235513750 -0.572002915 -1.277071366  0.181619743
##  [901] -0.632522669 -1.588284108 -1.501398872 -1.008541572 -0.670736924
##  [906]  0.246620884 -0.694525459 -1.550211175  0.174093950 -0.523113370
##  [911] -0.485017639 -0.689467035 -0.787580424 -0.451574519  0.331150553
##  [916] -1.253292478 -0.214670876 -0.043052379 -2.012482423 -1.569060669
##  [921] -0.166891200 -0.185107784 -0.336933317 -0.282039913 -1.244094714
##  [926]  0.019932753 -0.678058337 -0.448439338  0.282972255 -0.738189196
##  [931]  0.015765490  0.039282337 -0.614827126 -0.234373756 -0.960086620
##  [936] -0.096263808 -0.829933667 -0.366965989 -0.209252826  0.079606245
##  [941] -1.279710395  0.249402677 -0.855657637 -0.917919501 -0.236914496
##  [946] -0.023124165 -0.160411662 -0.973334383 -0.175534934 -0.926803884
##  [951] -0.453387660  0.328229354 -0.822597417 -0.200259798 -0.155287738
##  [956] -1.703084288 -1.518950586  0.010450507 -0.783052596 -0.200039732
##  [961] -0.476616639 -0.141109572 -0.132819875 -1.015067123 -0.764688099
##  [966] -1.058652106 -0.391781173 -0.665652417 -1.006729628 -0.368715542
##  [971] -0.624071031 -0.951353353 -0.110377688  0.215165566 -0.397634104
##  [976]  0.557414745 -0.321314362 -0.211078007 -0.767000226 -0.108665133
##  [981] -1.238125672 -0.499645905 -0.651260137 -1.047820253 -0.961434138
##  [986]  0.041410445 -0.715308591 -0.473069482 -0.085288401  0.462689041
##  [991] -0.675108666  0.307767382 -1.162376445 -0.197996835 -0.962224218
##  [996] -0.586133926 -0.523498540 -0.839339566  0.232400617  0.004948432
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
##   0.57659551   0.24501040 
##  (0.07747909) (0.05478227)
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
## [1] -0.02166705  1.39385026  0.10811268  0.51960425 -0.99796997 -0.32633691
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
## [1] -0.0377
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9487406
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
## t1*      4.5 -0.01471471   0.9134893
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 2 7 8 9 
## 4 2 1 1 2
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
## [1] -0.0079
```

```r
se.boot
```

```
## [1] 0.9290698
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

