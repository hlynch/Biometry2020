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
## 1 2 3 4 6 7 
## 1 2 2 1 2 2
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
## [1] -0.0191
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
## [1] 2.750336
```

```r
UL.boot
```

```
## [1] 6.211464
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##   2.5%  97.5% 
## 2.6975 6.2000
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
##    [1] 6.2 3.7 2.9 2.7 5.1 3.9 5.1 4.8 4.7 5.1 5.4 4.5 5.3 2.9 4.8 5.3 5.0 5.2
##   [19] 3.8 3.7 3.0 5.2 4.6 2.5 4.3 6.0 4.9 4.4 6.4 3.7 4.3 4.7 4.1 4.4 3.8 4.8
##   [37] 4.3 3.9 4.6 5.3 4.5 4.1 4.7 4.8 4.9 5.5 4.3 4.6 4.5 4.8 5.2 5.0 5.1 5.6
##   [55] 4.8 3.6 4.0 3.7 4.2 5.2 4.8 4.7 6.1 5.5 4.7 2.9 5.2 2.5 3.5 4.1 4.1 7.0
##   [73] 4.6 4.0 4.6 3.9 5.6 5.8 5.4 5.2 4.4 5.3 4.0 4.8 4.3 5.3 4.9 4.2 4.3 4.1
##   [91] 4.0 3.6 3.7 4.5 4.4 6.7 3.2 5.1 4.6 4.1 4.6 4.6 5.1 3.8 3.5 5.2 4.8 5.8
##  [109] 4.1 5.5 3.5 4.9 5.2 5.1 4.7 3.9 6.1 5.6 5.4 4.5 5.4 3.9 5.0 4.1 3.6 4.9
##  [127] 5.2 4.7 4.3 5.0 3.6 3.1 4.1 5.5 3.9 5.5 3.9 5.3 3.7 3.4 3.9 4.3 3.1 2.7
##  [145] 3.6 4.4 3.2 4.8 6.1 5.6 4.6 5.2 4.4 4.6 3.8 5.6 5.0 4.0 4.4 4.5 4.4 4.2
##  [163] 4.5 3.6 3.1 4.6 5.8 5.3 3.1 4.4 4.3 3.5 5.9 4.3 3.3 2.5 3.9 5.7 4.9 5.7
##  [181] 2.7 5.3 4.7 5.1 4.2 4.8 5.0 4.6 4.8 4.9 3.7 5.0 4.1 5.0 3.7 5.4 4.2 3.4
##  [199] 3.7 3.8 5.5 3.8 5.5 3.7 4.4 4.5 5.4 4.6 5.6 5.7 3.9 5.2 5.5 4.7 4.7 5.1
##  [217] 4.5 3.3 5.2 3.2 3.8 3.6 5.5 4.9 4.8 6.4 4.5 5.8 4.6 3.7 2.8 5.0 3.5 4.6
##  [235] 5.1 3.7 5.0 4.4 4.4 2.4 4.9 4.4 3.2 4.6 3.4 4.1 3.7 5.5 3.8 5.8 4.7 4.7
##  [253] 5.3 3.3 5.2 4.3 2.6 3.5 3.4 5.4 4.7 3.2 4.6 4.7 3.9 4.5 2.5 3.3 4.1 3.8
##  [271] 4.8 3.3 5.6 2.3 4.6 4.4 3.3 5.7 4.9 4.2 3.5 4.8 5.6 5.0 5.3 5.0 3.6 4.6
##  [289] 2.0 3.9 4.0 4.1 3.4 2.9 4.0 4.9 4.3 5.3 5.5 3.1 4.9 4.4 4.8 5.7 3.9 5.8
##  [307] 5.6 4.0 5.1 5.9 4.2 4.9 5.3 4.6 5.4 6.0 4.5 4.0 5.5 4.6 4.0 4.6 4.3 5.5
##  [325] 4.3 6.2 3.6 4.5 4.5 4.5 5.6 4.1 5.8 3.2 5.9 2.4 6.0 6.2 4.6 4.4 4.0 4.9
##  [343] 4.6 4.9 5.1 3.5 6.0 4.3 4.7 6.6 5.1 5.3 4.3 3.9 4.4 5.2 5.7 4.5 3.5 5.6
##  [361] 3.3 5.0 4.1 5.5 4.0 4.9 5.1 5.3 6.6 5.5 4.3 4.3 3.7 5.6 2.4 4.3 3.2 3.5
##  [379] 6.0 4.1 5.1 3.7 5.2 5.4 3.9 3.2 5.6 4.3 5.1 6.4 5.3 4.8 3.5 5.3 4.2 3.7
##  [397] 3.7 3.4 6.3 2.8 2.2 4.1 4.1 3.4 4.3 3.4 3.3 4.2 5.1 6.3 4.0 6.1 5.0 4.1
##  [415] 4.7 4.7 4.6 4.8 4.9 3.7 3.6 4.7 3.1 4.3 3.3 5.1 2.8 4.7 3.4 3.6 3.2 6.0
##  [433] 6.0 4.8 5.6 4.8 3.8 3.3 3.8 5.4 5.9 4.3 5.1 4.0 4.6 5.4 4.7 5.2 4.3 5.1
##  [451] 3.5 3.6 5.1 3.3 4.5 3.0 4.8 4.8 4.7 4.3 4.5 4.1 4.7 4.4 4.5 2.6 4.4 5.1
##  [469] 5.7 4.3 4.1 3.3 4.2 5.0 4.3 3.1 4.6 3.7 6.0 4.1 4.7 4.7 5.4 4.7 5.4 4.4
##  [487] 4.5 3.9 4.7 6.6 5.0 4.5 3.9 3.4 5.7 5.4 4.8 4.0 5.9 4.9 3.5 3.3 3.6 4.8
##  [505] 3.1 4.4 3.1 3.7 3.8 4.7 5.6 5.9 3.7 4.4 5.1 6.0 5.7 6.7 4.5 4.0 3.2 4.7
##  [523] 4.4 4.4 3.9 3.7 4.4 3.4 5.0 4.5 3.8 4.0 4.0 6.0 4.1 4.3 4.1 3.9 3.3 4.1
##  [541] 5.9 5.7 4.8 4.0 6.2 5.4 4.8 5.8 5.0 4.5 3.7 5.1 4.9 4.4 3.9 4.0 3.8 4.3
##  [559] 5.8 4.5 5.2 5.7 3.9 6.6 3.3 4.2 4.4 5.7 4.0 4.3 3.5 3.7 4.2 3.9 2.8 5.3
##  [577] 3.9 5.1 4.2 6.3 6.7 7.1 5.6 4.5 3.5 5.1 3.0 5.0 5.3 3.4 5.7 5.2 5.1 6.4
##  [595] 5.1 5.3 6.2 5.7 3.9 4.3 3.8 3.2 4.1 4.7 4.7 4.9 3.3 4.4 4.2 5.2 3.6 5.4
##  [613] 4.5 3.9 5.1 5.3 6.0 4.3 3.9 4.3 4.6 3.4 4.1 5.3 4.2 4.5 3.4 4.7 5.7 4.5
##  [631] 4.2 5.9 4.3 6.3 3.4 3.6 5.1 3.8 5.1 5.9 6.2 4.9 4.4 4.1 4.8 4.7 5.4 4.1
##  [649] 3.9 3.6 5.1 4.2 5.2 3.8 5.0 3.9 5.1 3.9 5.1 4.4 5.4 4.6 3.3 4.5 3.8 5.6
##  [667] 4.5 3.8 5.5 5.3 3.0 6.6 3.9 4.0 4.0 4.6 5.3 4.3 6.5 5.6 5.2 5.4 4.6 5.0
##  [685] 4.5 4.1 3.3 5.6 5.1 3.8 4.5 4.8 4.9 5.3 4.2 4.0 3.4 5.9 5.3 2.8 5.3 4.7
##  [703] 5.8 4.4 5.1 4.1 3.2 5.5 4.5 3.1 4.0 4.8 5.3 4.2 5.2 5.6 4.5 4.8 6.6 5.2
##  [721] 4.1 5.0 4.8 5.4 2.7 4.8 4.6 4.7 2.9 5.8 5.1 4.1 4.0 3.8 4.1 4.7 3.7 3.8
##  [739] 4.7 5.0 5.1 5.3 4.1 3.8 3.9 4.8 4.4 5.5 4.3 5.5 4.5 2.2 4.4 4.3 2.6 4.4
##  [757] 3.9 3.0 3.5 3.2 2.8 5.6 6.0 2.8 3.6 4.4 4.8 4.9 3.9 5.7 2.9 3.7 5.2 5.4
##  [775] 2.7 5.5 4.2 4.8 4.5 4.5 5.1 4.5 3.9 4.4 5.0 4.4 5.5 4.5 5.2 3.9 5.2 5.4
##  [793] 5.1 3.5 4.7 5.3 3.8 3.8 5.1 5.3 2.8 5.1 4.6 5.2 3.6 6.9 3.2 4.8 5.2 3.3
##  [811] 4.8 5.1 5.9 4.4 3.5 3.5 5.9 3.4 4.4 5.0 3.4 5.1 4.9 5.9 5.7 5.7 4.8 4.7
##  [829] 4.4 4.3 4.6 4.9 5.2 3.7 4.6 5.0 3.5 5.6 3.8 4.2 4.6 6.5 4.4 4.1 3.7 5.7
##  [847] 4.9 5.4 3.1 4.4 4.9 2.9 6.0 5.0 6.1 6.4 3.7 4.5 4.3 5.6 4.0 4.8 5.0 5.5
##  [865] 3.1 4.2 6.3 3.7 3.8 4.4 4.5 2.5 4.5 5.1 4.3 5.0 3.5 4.3 4.7 4.1 3.7 3.5
##  [883] 4.9 4.1 4.1 3.9 4.8 3.9 5.3 4.7 4.5 4.6 3.3 3.1 4.8 4.8 4.8 5.5 4.7 2.6
##  [901] 4.7 4.2 4.8 4.6 5.5 4.0 3.8 5.2 5.0 5.2 3.2 4.2 4.5 3.6 4.1 4.1 5.5 4.0
##  [919] 4.3 4.5 3.9 5.0 4.0 4.8 4.7 6.3 4.4 5.8 5.8 3.4 4.9 5.0 5.0 3.6 5.1 3.0
##  [937] 4.9 5.1 5.1 4.1 2.7 6.0 4.9 5.9 2.7 5.5 5.5 4.4 2.6 4.1 4.3 4.7 3.7 3.6
##  [955] 5.3 3.2 4.4 5.0 5.0 5.2 5.6 3.6 4.8 4.7 4.4 4.6 2.3 4.4 4.8 4.7 4.2 4.0
##  [973] 4.9 4.5 3.6 5.6 4.0 4.8 5.5 4.7 4.6 5.3 3.5 5.7 4.2 3.1 4.7 3.6 3.5 3.4
##  [991] 3.5 5.3 2.8 5.0 3.6 4.2 4.4 4.3 5.1 5.6
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
## 2.7975 6.2025
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
##    [1] 3.7 4.0 3.1 3.9 4.7 4.9 5.0 6.3 3.7 6.3 5.3 4.7 4.6 6.6 5.7 3.4 4.5 3.5
##   [19] 5.6 4.3 6.3 5.7 5.4 3.2 4.0 5.8 3.8 3.2 4.9 3.5 5.1 6.4 5.4 2.7 4.1 5.8
##   [37] 5.3 5.5 4.3 4.9 5.1 3.2 4.4 5.4 3.9 3.9 3.0 6.0 4.5 4.9 4.2 4.3 4.9 4.3
##   [55] 5.1 4.4 3.5 4.7 3.4 3.3 4.8 4.1 3.8 3.9 4.0 5.8 4.9 3.3 4.9 4.0 5.1 5.9
##   [73] 3.4 4.2 6.0 5.0 3.6 6.1 4.9 4.0 3.8 4.8 5.5 5.5 5.1 3.5 5.3 5.2 3.6 5.9
##   [91] 4.1 4.1 3.2 3.6 4.4 4.0 5.4 4.6 3.7 5.4 3.8 3.5 4.0 4.6 4.0 3.7 4.8 4.7
##  [109] 4.5 2.9 6.5 4.2 5.0 4.8 3.6 3.8 3.9 4.4 3.3 3.3 3.5 5.4 4.7 4.5 2.0 5.3
##  [127] 3.9 5.0 4.0 3.9 3.6 5.2 5.0 4.2 5.2 3.7 4.3 4.4 2.6 4.9 5.8 5.7 6.9 3.8
##  [145] 3.7 4.3 3.4 5.0 4.1 4.9 3.0 4.0 6.2 5.2 3.9 5.4 5.0 5.4 3.5 3.4 6.0 5.9
##  [163] 5.0 4.5 4.6 3.1 3.3 5.0 4.7 5.1 5.0 3.8 4.2 6.1 5.4 5.2 4.3 5.1 3.7 3.8
##  [181] 5.8 6.0 4.7 5.1 5.0 4.1 3.6 5.1 5.5 3.3 6.7 4.9 3.3 4.7 5.5 3.5 6.2 6.4
##  [199] 5.6 4.8 4.3 3.4 4.6 4.2 4.1 5.6 4.9 3.3 3.9 3.6 5.0 5.0 4.5 2.8 4.6 3.3
##  [217] 4.9 5.0 4.6 6.7 5.7 4.2 3.3 5.3 5.6 4.1 5.5 2.0 3.9 2.2 5.2 4.6 4.2 4.2
##  [235] 5.8 5.8 4.8 3.9 4.8 3.9 4.0 3.8 4.9 5.7 4.7 5.4 6.0 5.7 3.9 3.3 5.8 4.1
##  [253] 5.8 4.9 2.8 5.3 4.1 5.7 4.7 2.7 5.7 3.7 4.6 4.3 3.8 5.6 4.9 3.1 5.9 4.8
##  [271] 4.5 5.0 3.1 4.1 4.4 5.6 4.8 3.8 3.7 3.4 5.1 5.3 5.5 4.1 5.3 6.9 5.6 2.0
##  [289] 6.7 4.9 4.8 4.7 5.1 3.2 3.5 3.6 4.9 2.6 5.9 6.4 3.6 3.4 4.2 4.3 4.4 6.4
##  [307] 4.7 4.8 3.3 3.2 5.4 3.8 4.2 4.2 6.2 4.0 5.2 3.8 6.1 5.5 4.8 4.9 4.2 2.9
##  [325] 5.4 4.1 6.6 5.0 3.2 5.2 4.0 4.3 4.4 4.7 4.7 6.0 3.6 4.5 5.9 4.7 3.2 4.3
##  [343] 5.1 3.9 5.6 5.5 5.1 4.6 5.5 3.4 4.9 4.4 3.2 5.0 4.8 5.1 5.7 5.1 3.9 4.1
##  [361] 4.2 3.2 4.3 3.3 2.9 6.1 3.5 5.8 4.3 3.7 5.5 5.7 5.4 3.2 4.3 3.7 4.6 5.1
##  [379] 5.3 3.0 4.1 3.9 4.5 4.4 4.9 6.8 4.1 3.5 4.6 4.8 4.7 4.1 4.8 3.3 4.1 4.3
##  [397] 6.4 5.3 3.6 5.2 4.7 3.5 3.4 5.4 5.9 5.5 4.6 3.8 4.9 5.1 4.9 4.5 4.5 3.5
##  [415] 5.6 4.7 3.4 4.6 3.7 4.7 5.0 5.5 3.1 4.7 4.9 4.2 5.1 4.4 4.5 4.6 3.9 5.5
##  [433] 4.6 5.7 4.2 4.8 3.2 4.5 2.6 3.3 4.9 3.3 3.1 4.8 4.7 4.5 4.3 4.5 5.1 4.3
##  [451] 5.4 3.4 4.9 6.3 4.3 5.0 5.4 4.5 3.4 5.0 5.4 4.1 3.4 3.9 5.9 5.1 3.5 4.3
##  [469] 3.8 3.6 3.7 3.7 3.8 3.9 4.4 3.1 5.8 5.4 4.6 4.9 3.9 4.0 5.0 4.3 3.7 2.7
##  [487] 6.0 5.5 4.9 3.6 4.7 4.9 5.2 5.1 4.8 5.5 4.9 3.3 3.1 4.9 5.2 4.4 4.1 3.2
##  [505] 4.8 3.5 3.5 5.2 6.1 5.5 4.4 4.5 3.9 6.6 5.2 4.0 3.7 4.6 5.5 5.2 5.1 4.6
##  [523] 4.2 3.6 5.5 5.8 3.4 4.6 5.4 5.2 3.9 4.0 4.3 5.5 3.4 4.4 2.7 4.6 5.2 5.3
##  [541] 4.3 4.5 6.0 4.6 2.9 5.8 5.8 4.5 4.2 5.3 5.3 3.8 5.2 4.6 6.8 4.8 3.3 4.8
##  [559] 4.1 4.0 3.7 4.6 4.1 3.7 4.5 4.3 5.4 4.4 3.1 5.9 5.4 3.5 4.6 3.2 4.5 4.0
##  [577] 5.2 5.0 4.1 3.6 4.8 5.0 4.4 3.4 4.4 4.7 5.3 3.6 3.2 4.4 4.7 3.6 4.1 4.2
##  [595] 3.8 4.4 4.7 4.2 5.1 4.1 4.3 4.9 6.3 4.8 5.7 4.0 2.8 4.1 3.7 4.4 5.8 5.2
##  [613] 4.8 3.9 6.2 4.1 6.3 6.1 3.4 3.5 4.1 4.7 3.5 6.1 5.1 6.3 4.1 3.7 5.6 3.3
##  [631] 4.8 6.4 6.2 5.1 4.6 2.7 4.5 4.0 5.7 5.1 4.6 6.6 4.9 3.5 4.7 5.0 5.0 4.1
##  [649] 3.0 4.8 4.8 4.7 5.2 5.3 3.1 4.0 4.7 4.7 5.4 4.5 6.6 4.3 4.6 4.9 5.4 5.0
##  [667] 5.1 3.5 3.5 4.1 3.2 5.4 3.9 3.4 3.9 3.5 4.9 5.4 4.5 4.7 5.6 4.6 4.0 4.1
##  [685] 4.3 3.2 4.9 5.1 5.0 4.1 5.2 5.5 4.7 3.1 4.3 5.9 6.7 4.9 5.2 4.7 4.9 4.1
##  [703] 4.3 5.1 3.8 5.4 3.8 5.4 4.9 3.3 4.6 4.4 4.9 4.2 3.9 5.8 4.1 4.2 4.5 5.1
##  [721] 6.3 4.2 4.1 3.8 3.6 3.9 5.0 4.8 3.9 3.9 4.0 3.4 3.8 4.8 4.0 5.1 4.1 5.1
##  [739] 3.7 3.3 6.1 6.3 3.0 4.6 4.4 3.2 5.5 6.7 4.9 5.7 5.5 3.0 4.2 3.5 3.7 5.9
##  [757] 5.1 2.4 4.5 3.8 5.4 5.1 6.1 3.5 4.9 2.5 4.7 3.3 4.3 3.7 6.1 5.9 5.0 3.9
##  [775] 4.4 3.6 2.8 4.8 4.9 4.3 5.3 2.5 3.6 3.5 5.1 4.1 5.1 3.7 5.0 3.7 5.3 5.3
##  [793] 4.9 6.3 5.4 3.1 4.7 4.5 4.1 2.8 4.6 6.0 5.7 6.4 4.2 3.3 5.4 4.7 4.5 5.8
##  [811] 4.4 6.4 3.5 5.9 5.5 4.7 4.1 3.8 5.6 4.1 2.9 4.4 5.2 3.7 4.2 4.0 4.1 3.8
##  [829] 3.2 4.6 5.3 4.5 3.4 4.4 3.9 5.0 5.0 5.4 4.9 3.6 3.5 3.4 4.8 1.7 3.8 4.6
##  [847] 7.3 5.4 5.1 3.5 4.0 3.7 5.3 3.6 4.5 4.6 6.2 5.0 4.8 3.8 4.7 4.4 3.5 4.8
##  [865] 4.0 5.4 5.5 5.2 4.1 3.6 3.7 5.8 4.3 5.7 5.7 2.8 2.9 3.3 5.0 5.2 4.7 4.9
##  [883] 5.1 4.7 4.3 4.2 4.8 4.9 4.7 6.3 2.8 5.3 2.9 4.4 3.9 5.0 6.0 3.4 2.8 6.0
##  [901] 4.1 4.4 4.8 3.7 3.9 4.0 4.6 3.5 3.6 5.5 4.5 4.8 4.0 4.1 3.7 5.0 3.6 4.8
##  [919] 4.6 4.2 3.8 4.9 4.3 6.1 3.2 3.8 4.2 3.3 5.3 3.9 4.8 4.7 4.6 5.0 3.5 4.3
##  [937] 5.2 5.0 4.9 5.0 3.6 3.8 4.6 4.3 5.8 5.3 5.0 5.3 4.2 5.8 5.1 3.8 3.8 4.9
##  [955] 5.4 3.5 4.2 5.4 5.2 4.8 3.6 5.2 5.2 3.8 6.0 4.7 4.7 3.7 4.2 4.3 4.6 2.7
##  [973] 4.7 2.5 4.8 3.5 2.8 4.1 2.7 5.1 4.9 5.3 6.3 3.6 3.3 5.0 4.5 5.9 4.9 4.3
##  [991] 3.8 3.2 4.1 5.3 5.4 2.9 5.3 2.2 4.7 5.6
## 
## $func.thetastar
## [1] 0.0287
## 
## $jack.boot.val
##  [1]  0.507303371  0.454069767  0.351744186  0.208226221  0.071823204
##  [6]  0.006069364 -0.133923304 -0.237535014 -0.392238806 -0.539710145
## 
## $jack.boot.se
## [1] 1.015826
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
##    [1] 3.9 4.2 6.0 5.1 5.2 3.7 4.2 4.0 3.8 5.2 5.6 4.2 5.2 6.6 3.8 5.6 3.0 5.2
##   [19] 5.0 4.7 4.1 4.4 5.3 3.9 4.4 4.5 4.7 4.3 4.3 5.4 4.4 2.6 4.4 5.2 4.4 4.4
##   [37] 4.3 5.2 5.1 6.0 4.8 3.8 4.1 4.5 5.9 4.7 5.4 4.3 3.1 4.4 5.2 5.0 5.8 4.8
##   [55] 3.3 4.5 4.8 6.6 5.0 4.5 5.9 5.3 4.2 4.3 3.5 3.4 4.9 5.1 4.1 6.1 5.0 4.1
##   [73] 3.9 4.4 4.6 3.5 3.6 5.7 4.7 4.5 2.9 3.7 5.0 3.1 4.7 3.7 4.5 4.2 5.2 4.9
##   [91] 4.9 3.1 2.7 5.8 5.0 5.4 4.4 5.4 5.0 4.5 4.6 4.5 4.9 4.3 4.9 3.4 5.4 3.6
##  [109] 4.2 5.3 3.4 4.3 5.3 4.4 5.2 5.6 5.5 5.9 4.9 4.4 5.7 3.8 5.8 5.0 3.4 3.8
##  [127] 4.7 4.7 5.3 3.7 4.9 6.8 5.9 4.7 4.2 2.8 4.1 6.4 4.2 5.7 4.3 4.8 3.2 5.3
##  [145] 4.4 6.3 4.8 6.2 5.1 4.3 5.7 5.6 5.8 3.9 1.7 3.3 3.0 4.8 3.2 4.2 5.1 5.5
##  [163] 5.4 5.2 4.7 4.8 2.7 4.9 4.6 3.4 5.5 4.7 4.5 3.7 3.3 4.2 2.7 4.2 6.3 4.9
##  [181] 4.8 5.4 4.0 5.2 5.0 3.1 5.5 5.1 4.1 4.6 5.1 2.9 4.8 5.9 4.5 4.4 6.8 3.9
##  [199] 4.3 4.4 5.4 3.9 3.9 2.0 3.4 5.0 5.2 4.7 5.7 6.4 3.9 3.7 4.8 5.7 2.5 5.4
##  [217] 5.7 4.7 5.0 4.6 5.2 5.6 3.3 5.3 5.3 4.5 5.6 4.4 4.2 2.8 4.7 4.5 3.6 6.4
##  [235] 4.3 5.6 4.1 6.3 4.7 4.2 3.4 4.5 4.9 4.2 4.2 3.1 4.0 3.5 4.7 3.7 5.0 4.7
##  [253] 4.4 3.3 3.1 6.2 2.8 4.8 4.7 4.4 5.8 4.6 4.7 4.3 4.6 3.3 5.0 5.8 5.2 3.4
##  [271] 4.4 4.7 2.8 5.2 4.1 3.6 5.2 4.0 4.9 4.2 4.5 4.7 3.3 4.4 5.4 4.3 5.5 2.6
##  [289] 3.9 5.2 4.8 5.2 4.9 5.2 4.8 5.0 4.2 5.4 4.3 4.8 3.9 5.7 3.3 4.3 4.5 3.6
##  [307] 5.4 4.2 4.5 5.0 4.4 3.5 2.4 4.3 5.0 4.8 4.8 3.4 3.3 4.5 3.6 5.7 3.2 5.3
##  [325] 5.0 5.0 5.3 6.0 4.7 4.8 5.1 5.4 3.8 4.4 5.1 4.8 4.4 4.8 5.2 4.2 5.0 3.3
##  [343] 5.1 4.7 5.9 3.4 4.2 3.6 4.5 3.9 4.5 5.7 3.0 3.4 3.5 4.7 5.2 4.5 4.7 4.5
##  [361] 4.1 5.1 3.6 2.9 4.3 4.5 4.9 3.3 3.4 3.9 4.5 4.2 3.0 4.3 4.8 4.8 4.0 4.2
##  [379] 3.8 4.8 5.2 4.4 5.0 3.4 2.8 3.2 5.0 5.4 3.0 5.1 5.2 5.1 5.8 4.7 2.4 4.1
##  [397] 5.3 4.5 5.7 4.0 5.3 4.5 3.9 4.9 4.5 6.4 4.6 5.5 3.5 4.8 3.5 5.1 5.3 4.3
##  [415] 7.0 2.6 5.8 5.8 3.0 3.9 4.0 4.7 6.0 4.9 4.3 3.0 2.4 6.2 3.8 4.8 5.1 3.5
##  [433] 4.6 5.0 4.7 5.3 4.4 4.6 3.4 3.0 4.6 4.4 4.7 3.4 5.4 3.9 6.3 5.8 5.6 4.6
##  [451] 4.2 2.1 3.4 3.9 3.9 4.5 5.0 6.2 3.9 2.2 4.3 2.8 4.7 3.9 5.8 5.2 4.9 5.3
##  [469] 4.5 3.9 4.5 3.4 3.6 5.4 4.6 5.1 4.1 5.0 4.7 4.5 5.1 6.3 4.3 3.8 5.0 4.3
##  [487] 5.3 4.9 5.4 2.8 4.6 5.8 4.0 4.8 6.0 6.0 5.1 4.0 3.3 4.7 4.7 4.2 5.8 4.5
##  [505] 5.8 3.8 4.5 5.5 5.3 3.9 2.9 5.4 4.9 3.0 5.7 3.8 4.5 5.4 5.2 4.4 5.0 4.5
##  [523] 3.7 4.3 5.4 5.6 5.4 4.0 4.2 1.9 6.0 5.0 3.8 5.2 4.8 4.5 4.8 4.6 2.6 4.8
##  [541] 4.4 4.1 5.1 6.7 5.0 4.3 4.3 5.2 5.3 4.6 4.4 4.1 5.7 3.0 4.2 3.5 4.2 4.7
##  [559] 2.7 4.2 5.8 6.3 5.9 3.9 5.3 4.1 4.7 3.9 6.5 3.4 3.0 3.8 5.5 4.3 4.6 4.2
##  [577] 4.8 4.5 3.9 4.2 3.5 3.7 5.2 3.7 5.5 4.8 3.4 3.4 4.8 2.8 4.3 5.2 3.1 5.2
##  [595] 3.8 2.5 4.5 4.1 3.5 5.5 3.2 3.7 4.6 4.9 3.6 4.7 3.8 3.6 3.4 4.8 4.6 5.9
##  [613] 4.8 3.6 4.7 3.6 5.0 6.3 5.5 3.7 3.6 5.0 5.0 3.1 4.0 4.8 4.9 6.6 3.5 5.0
##  [631] 5.5 4.9 5.3 6.0 3.2 3.9 4.7 5.1 4.0 5.4 4.7 5.2 5.7 4.0 4.0 3.3 1.6 4.9
##  [649] 5.5 4.7 5.0 3.5 5.3 4.6 2.8 3.2 3.6 3.0 4.9 4.1 4.6 4.8 4.1 4.5 3.4 4.3
##  [667] 3.1 3.9 2.2 5.3 6.0 4.6 5.0 5.8 4.2 3.9 3.2 4.9 4.6 4.2 4.3 4.5 3.3 4.9
##  [685] 4.0 4.7 3.7 2.9 4.9 6.2 4.1 4.2 4.3 4.2 2.6 3.0 2.8 5.0 4.3 5.5 3.4 5.1
##  [703] 3.1 4.0 5.2 3.7 4.6 4.4 3.5 4.0 5.3 5.7 4.6 3.5 4.5 4.4 3.8 6.0 3.9 4.3
##  [721] 5.6 5.2 4.1 4.3 5.0 4.7 5.5 5.0 5.3 5.3 4.2 4.1 4.0 5.3 4.1 3.4 4.7 5.1
##  [739] 4.5 5.3 5.0 5.5 5.8 3.7 3.5 4.6 3.6 4.3 4.9 4.5 5.3 2.3 4.4 5.6 5.0 3.3
##  [757] 5.5 4.2 4.8 3.4 5.2 4.0 4.1 3.3 4.0 4.2 4.2 5.1 3.9 4.0 3.1 5.7 2.6 4.9
##  [775] 5.0 6.2 4.6 4.1 2.7 3.8 4.5 4.1 3.9 4.8 4.9 4.8 3.9 4.8 4.1 5.2 4.2 3.4
##  [793] 4.7 6.5 4.8 5.2 5.2 3.9 4.9 4.0 3.4 5.0 4.9 6.6 3.1 3.1 4.0 4.8 4.6 5.9
##  [811] 3.9 4.7 2.9 5.2 4.2 4.6 4.2 3.4 3.2 3.2 5.2 4.5 4.5 4.1 5.0 4.2 3.5 4.4
##  [829] 4.9 4.3 3.9 5.9 4.7 2.6 4.3 6.0 3.2 3.9 5.0 5.1 3.1 4.2 5.2 4.7 3.4 3.3
##  [847] 4.2 3.7 5.6 3.6 4.4 3.6 5.2 4.2 4.5 4.1 5.1 4.4 5.3 3.9 4.2 5.8 4.2 5.5
##  [865] 5.0 4.1 4.9 3.3 5.4 4.4 5.0 3.3 4.7 5.6 3.6 5.0 4.0 4.9 3.9 2.7 4.5 2.9
##  [883] 4.0 3.4 3.9 3.7 4.0 5.0 5.9 3.1 4.7 4.3 5.6 5.4 3.6 4.8 3.3 5.3 5.3 3.6
##  [901] 6.1 6.7 5.8 4.2 4.1 3.1 6.5 4.8 4.6 4.4 4.3 4.1 4.9 4.4 4.9 5.5 4.3 6.8
##  [919] 5.1 4.2 3.5 5.0 6.1 3.3 3.8 4.5 3.9 5.0 2.1 4.0 6.2 2.8 4.4 4.1 5.0 4.5
##  [937] 4.2 5.1 5.1 5.0 3.8 3.7 5.0 5.1 4.9 4.1 4.9 2.2 1.3 3.4 4.7 5.7 4.0 4.9
##  [955] 4.6 5.5 3.4 4.1 6.4 4.4 4.7 4.2 3.9 4.0 5.6 2.9 3.6 4.3 5.2 5.6 3.8 4.4
##  [973] 4.3 5.7 3.5 3.8 3.4 4.1 4.3 4.0 4.7 2.8 4.8 3.8 4.4 5.4 5.3 6.2 5.7 3.3
##  [991] 4.6 4.4 4.0 6.2 4.8 6.0 5.4 4.2 3.5 5.4
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.500 5.300 5.200 5.200 5.000 5.000 5.000 4.700 4.664 4.500
## 
## $jack.boot.se
## [1] 0.8861364
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
## [1] -0.3842826
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
##    7.363413   14.709952 
##  ( 3.221175) ( 6.659387)
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
## [1]  0.9453285 -0.2406668  0.3448688  1.7271260 -0.4968285 -0.4338740
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
##    [1] -0.6553394142 -1.2582199291 -1.1219095850 -0.3675055551  0.1328731390
##    [6] -1.0228498150  0.1926422608 -0.7881833421 -0.1448973127 -0.2810903480
##   [11] -0.6005319454 -0.2658988331 -0.2845519502 -0.0513092802 -1.9178762430
##   [16]  0.0247110481 -0.7499290214 -0.7720608551 -0.2429986810 -0.4950720516
##   [21] -0.5537677768 -1.0967464007 -0.0564008480 -0.0884129427  0.5886749766
##   [26] -0.3638238571 -0.0241126756  0.3128267751 -0.3516949573 -0.5966861185
##   [31] -0.2004305357 -0.4162653875 -0.6238566992 -1.7222525762 -0.0942572313
##   [36] -0.5594884270 -0.3141251647  0.2005955426 -0.5190244332  0.0394322037
##   [41] -0.4454195589 -0.2517821513 -1.4271211521 -0.2347912179 -0.3197410301
##   [46]  0.3969816347 -0.4974111400 -0.1088754828  0.0009661700 -1.2757534125
##   [51] -0.5957984784  0.2513083892 -0.1127849611 -0.2208835755 -0.8253725272
##   [56] -1.2398527799  0.0950910746 -0.2994668983  0.2102583576 -0.3892514777
##   [61] -0.1393813077  0.0961829344  0.3191787520 -0.2098998475 -0.0948145354
##   [66] -0.4018252864  0.0853984174 -0.0945283138 -0.5657490342 -0.3424179093
##   [71] -0.5883693950  0.1353035492 -0.3942915721  0.4442631739 -0.3130811903
##   [76] -0.2983923447 -1.2283545031 -0.6874567720 -0.2658713025 -0.3685348570
##   [81] -0.9171221633 -0.9133866615 -0.2582051376 -0.4063550672 -0.3166270514
##   [86] -0.5144619773 -0.4311699392  0.4963227079 -0.0140456270 -0.5746360405
##   [91] -0.4745302488 -0.3511669940 -0.2519034365 -0.9970717744 -1.0004535077
##   [96]  0.1554705609  0.3583203840 -0.1668210727 -0.2569838882  0.7974500676
##  [101] -0.4283323901 -0.2920314542  0.0610053277 -0.1727311193 -1.4886292270
##  [106]  0.2844939719 -0.7629950755 -0.6194741231  0.2172403897  0.2864961783
##  [111] -0.3481452118 -0.2386331946 -0.3391701317 -0.1284679990 -0.8012131806
##  [116]  0.0689262866 -0.3941171379  0.2561364451 -0.0228123430 -0.6830898253
##  [121] -0.1070382765 -1.1545387712 -0.6948307847 -0.5328034403 -0.4748675584
##  [126] -0.4061532022 -0.3173417785 -0.3767209619 -0.7942340967 -1.4340907222
##  [131] -0.0757293983  0.0705435642  0.2872233256 -0.2356267532 -0.6937312647
##  [136] -0.6389093750 -0.2501581421 -0.4220413559  0.2360188646 -0.4337904555
##  [141]  0.3409100925 -0.1632393367 -0.0306175858 -0.9038972096 -0.8451414442
##  [146] -0.6079439609 -0.4864137311  0.4904737325  0.0378071635 -0.7110771703
##  [151] -0.3426404268 -0.8998358390 -0.2598006250 -0.4500257356 -1.5092771492
##  [156]  0.0204469343  0.0007326236 -0.7451985093 -0.0322344185 -0.1793627356
##  [161] -0.5294770001 -0.3328772521 -0.7680383601 -0.3164304828 -1.4787315758
##  [166] -0.7179560680 -0.3699147941 -0.5100295166 -0.3403682520 -0.4617346449
##  [171]  0.5296851673 -1.7197346469  0.1971032339 -0.0412993064 -0.1796450638
##  [176] -0.2745050854 -0.1589898323 -0.5555038948 -0.2353855584 -0.9628219669
##  [181] -0.7907722362 -0.5523420771 -0.0026109556 -1.0677228670  0.0594619365
##  [186] -0.0521745652 -1.0072320227 -0.3616913627  0.0517988032 -0.2861803885
##  [191]  0.2650221819 -0.1526521496 -0.4705991170  0.1079209944 -0.9504622178
##  [196] -0.0452530322 -0.6970682493 -0.5063026865 -0.7592283139 -0.3113156146
##  [201] -0.8400762968  0.1592739872 -0.0069059025 -0.2056510411 -0.5340576018
##  [206]  0.2855114244  0.0892016449 -0.6604232247 -1.0629936499 -0.7261336658
##  [211] -0.6147048324  0.0068033509  0.2350658949 -0.0922942794 -0.3590535175
##  [216] -0.0675024785 -0.3968352192  0.1123039291 -0.1097340559 -0.2770745353
##  [221] -0.3877510525 -0.5937640750  0.2946010096 -0.3692262259 -0.8216251944
##  [226] -0.6863022951 -0.4808902695 -0.7390567930 -0.9608940875 -0.5908945648
##  [231] -0.3228343347 -0.0312541944  0.0943535692  0.0262936718 -0.3139377745
##  [236] -0.8831243005 -1.0361226291 -0.7911299100  0.0106354873 -0.8653188536
##  [241] -0.2247773735  0.5149659602  0.6404904613 -0.5872553809  0.7720811959
##  [246]  0.6584608639  0.0950269990  0.0394006215 -0.3211476490 -0.5187949775
##  [251] -0.4721776126 -0.0168068299 -0.0559563770 -0.0279476567 -1.3778930006
##  [256] -0.0352258671 -0.4437885795 -0.2108933590  0.4183399246 -0.2037609735
##  [261] -0.4314482199 -0.6932327050 -0.1308508913 -0.2499886745  0.1335016958
##  [266] -0.5286413241 -0.7460617725  0.3794157188  0.0644291252  0.1930050746
##  [271] -0.4187055102 -0.6598180206 -0.6523242760 -0.4629509767 -1.4068027126
##  [276] -0.2143325560 -0.4286015789  0.3148795004 -0.4043811100 -0.2209956560
##  [281] -1.0948918557 -0.6540640263 -0.2334080546 -0.0380003471 -0.7804039163
##  [286] -1.0810897890 -0.5798377343  0.4728565345  0.3115168386  0.2151535203
##  [291] -1.1210219636 -0.7686990695 -0.5522951115  0.1005983259 -0.0893987191
##  [296] -0.1481975709 -0.2851404575 -0.1507992947 -0.8701604434 -0.0044354948
##  [301] -0.7098244707 -0.1179141979 -0.3415927643 -0.3578402476 -0.2656143014
##  [306] -0.4176863710  0.1148924002 -0.7685177248 -0.4294314770  0.2630268979
##  [311] -1.0524650359 -0.0347871324 -0.1676055477  0.3384753735  0.1076554557
##  [316]  0.0144502223 -0.4852543953 -0.1736956112  0.2041727603 -0.5869983420
##  [321]  0.0556377393  0.1813370104 -0.6326388227 -0.3376682556 -0.3174161051
##  [326] -0.2684625707 -0.2407974082 -0.2913549618 -0.7927691760  0.8306297984
##  [331] -0.0840406782 -1.1165506124 -0.8741789989 -0.6740901748 -0.3224540070
##  [336] -0.2319887515 -0.2499886745 -0.7503414846 -0.0123026063  0.0717625691
##  [341] -0.0753122469 -0.5119486962  0.6133976984 -0.6153363140  0.2958004021
##  [346] -0.7637948556  0.0906675480 -0.5940610753  0.0877786933 -0.1461022895
##  [351] -0.4240260381  0.3158230047  0.3465392104 -0.8982397662 -0.1520929298
##  [356]  0.1636699987 -0.0080768914 -0.2976866165 -0.3738676930 -1.2781766516
##  [361] -0.1887139910 -0.1319230318 -0.2553866654 -0.6024739453  0.1855124873
##  [366] -0.2641963899 -0.3182298922 -0.1393704052  0.0330014703  0.0231752997
##  [371] -0.0934447981 -0.3511748139 -0.0295413383 -0.0242832753 -0.1180196793
##  [376] -0.6517207825 -0.6767373190 -0.5951680441 -0.4807941696 -0.3012851830
##  [381] -0.3057523111  0.0282545199 -0.8283176537 -0.2952331584 -0.4814992443
##  [386] -0.0137159388 -0.7168218251 -0.1685358373 -0.9419107664 -0.2769666617
##  [391] -1.1408435291 -1.1095854958  0.1219363618 -0.3271922952  0.2595238743
##  [396]  0.1796911952  0.0264848610 -0.4235299916  0.2866470435 -0.3618970970
##  [401] -0.9526076874 -0.0238249383 -0.0032491122  0.2696197997 -0.0645745983
##  [406] -0.1421976151 -0.0213476300 -0.3032176820  0.0923702499 -0.4984027760
##  [411] -0.6691655762 -0.0953924457  0.0756220933  0.3258218881 -1.4667258582
##  [416]  0.1418701615 -0.3899860301 -0.3175114291 -0.4569623938 -0.2919945068
##  [421] -1.0119080459 -0.2772632312  0.2583131376  0.6819504266 -1.0815376940
##  [426] -0.2627145310 -0.1113645396  0.1095526876 -0.5120642371 -0.8051291799
##  [431]  0.2930926293 -0.1708710191 -0.8051330666 -0.2298951686 -0.1696415973
##  [436] -0.2969312222 -0.5390388256 -0.2705379501 -0.8914941301 -0.2846982661
##  [441] -0.1234130544 -0.2912222469 -0.6175211503 -0.3213992755 -1.0101620092
##  [446] -0.4351712265 -0.3455761692  0.2389151071 -0.4904317279 -0.3175114291
##  [451] -0.2235084732  0.5970429860 -0.5985562483  0.1509727383 -0.7771239595
##  [456] -0.5974271058  0.0595370360 -1.0170097981  0.1037205901  0.4215036704
##  [461]  0.5172623011  0.0349808670  0.0415618918 -0.4796391878 -0.4633818054
##  [466] -0.2740732409  0.3421794249 -0.9956102106  0.2443854623  0.4730289546
##  [471] -0.2316414906 -1.1652014021 -1.4387217198  0.2175452215 -0.0455653020
##  [476]  0.2109612752 -0.5118348265  0.1091962278 -0.7066459738 -1.0990196309
##  [481] -0.5546238961  0.7496904698  0.1613424152 -0.2554161588 -0.2721383449
##  [486] -0.6635144957 -0.2218361575 -0.9282207552 -1.2009000609 -0.1931314680
##  [491] -0.1550821244  0.2403387986 -0.6217785858 -0.3059512225  0.9195627385
##  [496] -0.5739849207 -0.0763655872 -0.9268521817 -1.0497247450 -0.2199244742
##  [501] -0.5447047413 -1.4397639502 -0.6747101984  0.9864892666 -1.0517299253
##  [506]  0.1110129964 -0.0193329996 -0.6815160300  0.0324840856 -0.3901250742
##  [511] -0.6299923536 -0.5079683362 -0.1358514920  0.2573119625 -1.0281230126
##  [516]  0.2841734694 -0.4811198130 -1.1554957746  0.0962852602 -0.4751454835
##  [521] -0.1272651468 -0.8909170726 -0.0980755394  0.9256799065 -0.6912035932
##  [526] -1.3172065059 -0.5497428299 -0.8646539057 -0.1112816637 -0.1763485799
##  [531] -0.5456941989 -0.5701793432 -2.0556336128 -0.5126841996 -0.2763728662
##  [536] -0.1857670583 -0.0841241408 -0.1278378971 -0.1023435396 -1.0775232263
##  [541] -0.3333392302  0.1368884010 -0.0152569190 -0.4450947255 -0.7678999301
##  [546] -0.0260495358  0.2682015461 -0.1729593567 -0.4729125968 -0.0756415544
##  [551] -0.2174118067 -0.7832059241 -0.0226851918 -0.0060362262  0.5220010399
##  [556]  0.1978926799 -0.2351058790 -0.8398536542  0.1154278202  0.2147511145
##  [561] -0.3337169772 -0.3014464762 -0.4028143346 -0.5191175906 -1.0170097981
##  [566] -0.4805558292 -1.3583069292  0.1239576153  0.1379351654 -0.2086259885
##  [571]  0.4611203893 -0.5311714376 -0.0607942306 -0.2755381635 -0.3879668403
##  [576]  0.0402717941 -0.8613651367 -0.0614192556 -0.8198172234  0.8973044614
##  [581] -0.3631026642  0.6084757461  0.1866781130 -0.9226595250  0.0537437669
##  [586] -0.3318999775  0.1484274860 -0.7014504288 -0.3996069566 -0.5409030302
##  [591] -0.5854239039  0.4865703840  0.2149051751  0.0995084445 -0.2911591644
##  [596]  0.5077345695 -0.6416251613 -0.2034286678 -0.2468831713  0.0702701226
##  [601]  0.3622665439 -0.6530113789 -0.9522945744  0.1492261828 -0.2780506159
##  [606] -0.9057785055 -0.8591945975  0.2006390697  0.0792497457 -0.0690890550
##  [611] -1.2457715077  0.5865669275  0.0569023886 -0.7628847934 -1.0557299154
##  [616] -0.2783735671  0.7263020285 -0.0414739999  0.0321046676 -0.1143929235
##  [621] -0.5657770445 -0.3346667303 -0.9972987046 -0.5546380716 -0.0513077272
##  [626] -0.1579442414 -1.0717776247 -0.8272108376 -0.9372157732  0.3238182849
##  [631] -0.9167526512  0.0092614835 -0.0009357756  0.0264848610 -0.9293475367
##  [636] -0.5801933900 -0.5412902206 -0.0471278301 -0.3207319190 -0.4674824320
##  [641]  1.0132694617 -0.5380167576 -0.5563562946 -0.4030574447  0.3212719980
##  [646] -0.1373546278 -2.5742240467 -0.4606164329 -0.5941086764 -0.0151386435
##  [651] -0.5738068644  0.7904015764 -1.4345826288 -0.6644715864 -0.2123977140
##  [656] -0.2286642062 -0.7929493808 -0.4828753251  0.1363061644  0.0160592127
##  [661] -0.5067821611 -0.4346359742 -0.7123437326 -1.0115628561 -0.9264724846
##  [666] -0.4306314320 -0.5507195732 -0.6893903753 -1.1820175820  0.0351475638
##  [671] -0.7914962533 -0.8616265205 -0.1312423110  1.3679364081 -0.9357778691
##  [676] -0.5781344909  0.0086301585 -0.6220585025  0.0567221923  0.6298884504
##  [681]  0.6406001173 -0.1673345633 -0.0112826911 -0.1295773878  0.4327440987
##  [686] -0.2566675008  0.1425587973 -0.1457039689  0.0988742860  0.2751232205
##  [691] -0.2475146957 -0.1428438492 -0.7116108749 -0.5559716740 -0.2098682525
##  [696]  0.1973432282  0.8018529462 -0.2856850820  0.1868098810  0.4782733772
##  [701] -0.6647176390 -0.7883428923 -0.0502941939 -0.5373582421  0.2707263619
##  [706] -0.3554625676  0.0314736219  0.1502737317 -0.4237946688  0.2874157308
##  [711] -0.3664658367 -0.5354072131 -0.7136654294  0.0297178147 -0.4912786313
##  [716] -0.8520081591 -0.5558185964  0.0230052676 -0.2072236831 -0.1031423110
##  [721] -0.9555413894 -0.7025145115 -0.3166388480  0.0035875334 -0.2366791028
##  [726] -0.1369939182  0.3427316704 -0.7814445964 -0.0957015465  0.1730756353
##  [731] -0.8377396516  0.0834733793  0.4434851364  0.1671466471 -0.5232145638
##  [736] -0.0457853704 -0.5680560250 -0.9528876752 -0.0609166229 -0.0502349536
##  [741] -0.0454776887  0.1329022482 -0.5478228043 -0.5595237430 -0.8034305602
##  [746]  0.4084874569  0.3132187938 -0.8220703346 -1.4452326180 -0.6805755406
##  [751] -0.5801434787 -0.4354394331 -0.2680247879 -0.6143479671 -0.4907489177
##  [756] -0.2766267510 -0.7350870266 -0.5748886779 -0.7625897578 -0.4605173961
##  [761] -1.1320994533 -0.5564616655 -0.4243715840 -1.0078819190 -0.8604294433
##  [766]  0.0326071375  0.7224254278  0.5642990232 -0.4203209191  0.1713339064
##  [771] -0.8721227308 -0.3327820008 -0.5907684382 -0.2424795986 -0.3053548849
##  [776] -0.2675312090 -0.2386723260 -0.5910208062 -0.2598873210 -0.5687011928
##  [781] -0.9092585254 -0.1207076836 -0.6189682806 -0.2048315503  0.2710142360
##  [786] -0.7744625043  0.2043929439  0.2958574323 -0.2529625640  0.2527644121
##  [791] -0.1827895608 -0.1730769442 -0.8653147808 -0.2645429990 -0.6952368194
##  [796] -0.2208835755 -0.4826242381 -0.4950720516 -0.6146734005 -0.6325737357
##  [801] -0.4827884534 -1.3031500503  0.1013541856  0.2669835609 -0.8539774614
##  [806] -1.9858827799 -0.7663465083 -0.7928075283 -0.7483566347 -0.1245641470
##  [811] -0.8085291338 -0.5116847989 -0.0844342050 -0.3296226413 -0.4851934470
##  [816] -0.0806115918  0.0300706642 -0.2887355180 -0.9370320623 -0.7332578689
##  [821] -1.2386778344 -0.9180781637 -1.1322481247  0.4255320908 -0.3400580664
##  [826] -0.8910737624 -0.4863399593 -0.2906294823  0.1287557356 -0.8419005837
##  [831]  0.1058319112 -0.2261125122 -0.4449391031 -0.2073144428 -0.5985200229
##  [836] -0.7216554520 -0.4950888779 -0.7298040512  0.1897191239  0.0182733125
##  [841] -0.3897376249  1.3670367263 -0.5185515633  0.1505777814 -1.0896642796
##  [846] -0.7422153826  0.1169306344 -1.1104303500 -0.1527692395  0.1316143907
##  [851] -1.0859778876 -0.9389914705 -0.4666247501 -0.7977647769 -1.5298923443
##  [856] -0.4467514935  0.1002193606 -0.2571100966 -0.2316414906 -0.7930825661
##  [861] -0.3599500631  0.1395414709 -0.8838614634 -0.4838516908 -0.5452001026
##  [866] -0.1716876238  0.1624414924 -0.3133528927 -0.3232168860  0.3382175122
##  [871] -0.7187720583 -1.2458762784 -0.5296255461  0.2546039921 -0.1236746978
##  [876] -0.4913622207  0.6041552091 -0.7392228942 -0.5613671122 -0.9772021486
##  [881] -0.2015697374  0.0015348941 -0.1880727354 -0.2328975866  0.2842570348
##  [886] -0.3596244877 -0.3506337883 -0.4375212086 -0.7466141784 -0.3111388435
##  [891] -0.3429192327 -0.8404680859 -0.6912035932 -0.2323703528 -0.3299975106
##  [896] -0.3649696710 -0.2463624250  0.1840440479  0.3186925152 -0.2939131509
##  [901] -0.8682647577 -0.6838598590 -0.2283820652 -0.7262195568 -0.5869149463
##  [906] -0.2439779415 -0.2160889027 -0.5149790903 -1.0379953836 -0.3605673215
##  [911] -0.1660810252 -0.0506112556 -0.0614408605 -0.3538019281  0.4272941614
##  [916] -0.4035621461 -0.7749383688  0.6017437508 -0.6073796305 -0.1619872994
##  [921]  0.0740504154 -0.6669886270 -0.1875025306 -0.5951680441 -0.1575680598
##  [926] -0.1021216622 -0.2464430194 -0.2941782244  0.3551649854  0.5163677368
##  [931] -0.7014983424 -0.4695476031  0.0131784857  0.3816783142  0.2876361337
##  [936] -0.6136438118 -0.5953254484 -0.2168638597  0.1771154411  0.0351820539
##  [941] -0.1351407005 -0.1199729290  0.0899450160 -0.1931314680 -0.2981815754
##  [946]  0.4338482734 -0.4512755103 -0.2304300567  0.3025321629 -0.1462888022
##  [951] -0.1779719976 -0.6766877926 -0.1871247365 -0.3209555253  0.3150606952
##  [956]  0.0992019622  0.0442357253 -0.0611450796 -1.3530922106 -1.2846734808
##  [961] -0.8754736314 -0.3375831594 -0.1390483350 -0.2233326261  0.5648173453
##  [966] -0.8357323704 -0.0074943071 -0.2025707417 -0.6259496254 -0.6496417961
##  [971] -0.1512887441 -0.7546811217 -0.6730206350 -1.1432391587 -0.1103952736
##  [976]  0.5654403731 -0.1793627356 -0.3856579128 -0.7007717728 -0.4409847304
##  [981] -0.6625970476 -0.1362343935 -0.3374179995 -0.5354574422  0.2200226039
##  [986]  0.4942801492 -1.0898576897 -0.1503183542  0.6657759917 -0.5957366394
##  [991] -0.4892230439  0.0527059075  0.2075689723 -0.3767121728 -0.1765078543
##  [996] -1.0376638619 -0.2759491273 -1.2090994796 -0.3478376086 -0.3199257480
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

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced
```

```r
fit2
```

```
##       mean          sd    
##   0.50057107   0.16683482 
##  (0.05275780) (0.03729998)
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
## [1] -0.74794359  0.29763658 -1.08339403  0.89367970 -0.06004144  0.55558352
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
## [1] 0.057
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8643831
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
## t1*      4.5 -0.05605606   0.8726419
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 1 2 4 5 7 8 
## 2 1 3 1 1 1 1
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
## [1] 0.0292
```

```r
se.boot
```

```
## [1] 0.9270959
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

