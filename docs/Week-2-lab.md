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
## 0 1 2 4 5 7 
## 3 2 1 2 1 1
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
## [1] 0.0089
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
## [1] 2.700698
```

```r
UL.boot
```

```
## [1] 6.317102
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
##    [1] 4.5 3.9 5.9 4.6 4.3 3.2 3.8 4.5 4.9 4.2 3.6 4.8 5.7 4.7 4.3 5.8 4.2 4.7
##   [19] 5.2 4.6 3.8 5.8 3.9 3.9 4.7 4.9 4.1 2.3 6.5 4.8 4.7 3.7 5.5 4.5 4.2 5.6
##   [37] 5.3 4.8 6.5 4.8 4.6 3.8 4.8 5.4 4.4 4.3 6.8 3.0 5.4 3.9 4.3 3.5 5.6 4.8
##   [55] 4.2 3.3 3.8 5.6 2.2 6.4 7.1 4.4 2.6 6.7 3.8 4.3 5.2 5.2 4.4 5.9 5.1 5.3
##   [73] 7.4 4.0 4.4 5.6 4.0 5.0 4.0 3.2 4.8 4.3 3.9 2.9 3.6 4.0 5.0 4.5 5.3 4.6
##   [91] 5.8 5.3 3.3 4.1 4.7 4.5 5.0 3.9 4.0 5.9 4.0 5.1 5.2 3.0 3.8 4.9 6.5 4.8
##  [109] 4.7 4.0 5.5 4.6 3.2 3.5 4.1 4.0 4.5 3.4 5.6 3.7 5.2 5.4 5.6 6.0 5.2 5.3
##  [127] 3.7 3.8 5.5 4.4 4.5 4.2 5.0 3.0 4.1 6.2 4.0 4.6 5.9 5.2 4.2 5.4 3.0 3.6
##  [145] 4.5 4.5 4.5 3.7 6.2 4.9 5.6 2.9 4.8 4.3 5.5 5.4 3.8 4.3 3.4 5.7 4.2 5.3
##  [163] 6.0 5.9 4.7 3.5 4.5 4.2 4.0 4.1 4.1 5.7 4.0 4.7 5.2 4.5 4.8 4.6 6.3 4.8
##  [181] 4.8 5.8 3.8 5.8 4.1 4.3 3.8 2.4 5.5 4.7 4.0 5.0 3.8 4.8 4.3 5.4 3.8 5.1
##  [199] 4.1 3.6 5.5 4.6 3.2 6.6 4.8 2.3 4.5 4.5 6.5 4.5 4.9 5.1 4.9 3.2 4.9 3.8
##  [217] 5.0 5.7 4.8 5.4 5.7 5.4 5.5 5.3 4.1 4.2 5.1 4.9 4.8 3.8 4.6 4.9 4.7 3.7
##  [235] 5.6 3.6 4.1 4.5 4.7 3.8 4.2 6.3 5.7 4.8 4.5 4.8 3.3 3.8 4.1 5.2 4.2 4.3
##  [253] 4.8 4.8 4.4 4.0 4.5 4.3 3.1 5.2 3.2 3.2 4.4 4.3 5.2 4.0 5.8 3.5 5.2 4.8
##  [271] 4.3 3.9 5.2 4.1 4.5 3.0 6.1 6.6 4.4 5.0 3.8 3.2 3.7 5.7 3.9 5.5 5.0 4.7
##  [289] 4.0 4.4 4.5 4.9 5.1 4.2 3.8 4.8 3.7 5.0 3.4 2.8 4.5 4.8 5.7 3.3 4.3 4.3
##  [307] 4.8 4.2 5.0 5.0 3.6 6.0 3.6 5.2 3.7 3.0 4.8 4.4 4.5 4.8 5.6 4.1 6.1 4.2
##  [325] 3.6 4.1 3.2 5.2 3.5 4.9 3.6 3.8 4.9 4.4 4.5 5.4 5.3 4.8 3.3 3.7 5.7 4.6
##  [343] 4.6 4.4 3.7 3.4 4.7 5.9 5.0 3.2 5.8 6.1 4.2 5.1 4.7 4.2 4.6 5.4 4.0 3.8
##  [361] 3.5 5.0 3.0 4.1 3.5 4.2 3.4 4.2 4.4 4.9 4.3 4.2 5.1 4.4 4.6 4.0 4.4 5.7
##  [379] 3.7 5.9 3.6 5.6 4.3 5.2 3.7 6.2 4.5 6.1 4.9 4.7 5.3 4.5 4.4 4.2 3.5 4.9
##  [397] 5.1 3.9 6.0 5.3 5.3 3.6 5.9 6.0 7.1 4.7 5.5 5.2 4.5 2.9 3.2 5.6 4.4 3.8
##  [415] 4.6 4.3 5.2 4.3 6.1 4.5 5.0 4.3 3.8 4.5 4.7 6.4 3.7 4.8 6.0 4.9 5.1 6.0
##  [433] 5.2 4.5 3.3 5.6 5.2 4.4 3.9 5.0 5.7 5.9 4.8 3.8 4.3 4.2 4.0 3.7 5.0 4.8
##  [451] 4.7 3.4 4.9 3.6 6.3 3.9 3.9 4.0 4.8 3.7 3.7 3.4 5.2 4.7 6.7 3.3 3.3 6.0
##  [469] 2.9 4.8 4.2 4.3 4.4 4.2 4.7 3.9 4.7 4.7 4.9 4.7 2.5 4.5 4.2 5.0 4.5 5.1
##  [487] 4.3 3.2 4.1 4.5 2.8 3.3 4.9 3.4 3.9 5.7 5.0 5.2 3.9 3.9 3.7 4.4 4.7 2.5
##  [505] 4.4 4.9 3.2 3.9 3.8 4.0 4.1 4.1 2.7 4.3 4.7 3.1 3.3 4.6 3.7 5.7 6.3 5.5
##  [523] 3.0 5.2 4.4 6.0 4.8 4.3 4.3 4.5 5.3 3.8 4.1 4.6 5.6 4.6 3.9 5.1 4.3 5.1
##  [541] 5.0 4.2 4.9 4.9 4.7 3.9 4.2 2.9 4.5 4.6 3.4 4.6 5.3 3.1 4.2 4.9 4.8 4.1
##  [559] 5.7 5.1 3.0 3.1 2.4 6.1 4.2 5.5 3.8 3.2 5.2 5.2 3.3 5.6 4.9 3.9 5.8 4.5
##  [577] 5.5 4.1 5.6 2.8 5.2 4.1 4.1 3.9 4.7 4.5 5.1 5.4 4.1 5.0 3.3 5.1 4.8 4.4
##  [595] 5.1 3.9 4.5 3.6 4.8 4.7 4.5 6.0 2.9 5.7 3.9 5.6 3.5 5.4 4.9 2.8 5.3 4.1
##  [613] 5.6 4.0 3.2 5.2 6.8 6.6 3.5 3.6 5.6 4.4 5.0 4.3 5.1 5.0 4.1 4.6 4.5 3.2
##  [631] 3.2 4.9 6.4 5.0 3.0 5.0 3.5 5.3 6.0 5.4 5.3 5.3 3.7 4.9 2.8 4.3 6.5 2.3
##  [649] 5.6 5.2 2.0 4.8 4.9 4.5 4.0 4.5 4.9 5.4 6.1 3.6 4.5 3.4 5.0 4.6 3.3 5.4
##  [667] 4.9 4.4 4.4 5.1 4.1 5.1 4.0 3.8 3.7 4.6 3.7 4.3 3.8 4.5 4.3 5.8 2.8 5.2
##  [685] 4.8 6.2 3.0 5.4 5.7 4.2 4.8 5.0 3.5 3.1 5.9 3.3 6.7 3.0 5.3 5.1 3.2 4.1
##  [703] 4.6 6.1 4.3 5.0 5.2 4.3 4.8 5.2 4.1 2.8 4.6 4.9 4.0 4.0 5.0 6.3 5.1 3.7
##  [721] 3.0 3.8 4.9 5.2 3.6 4.9 5.6 4.1 4.6 4.3 5.0 5.7 4.0 3.5 6.0 5.4 3.8 4.4
##  [739] 4.5 5.0 4.4 5.1 2.3 4.7 4.4 4.7 5.3 3.7 4.4 4.8 4.8 5.8 3.5 4.4 5.2 5.7
##  [757] 5.1 6.0 4.5 3.5 6.8 3.6 3.9 4.6 2.9 5.4 4.9 5.3 4.5 4.7 4.8 4.3 4.7 3.1
##  [775] 4.7 6.1 2.9 4.1 3.5 2.8 6.1 6.2 3.3 5.0 3.8 4.9 4.8 5.2 4.6 4.1 2.6 4.5
##  [793] 4.9 3.4 5.6 5.4 4.6 5.4 4.3 3.3 4.3 5.2 3.2 4.9 4.7 5.6 4.8 4.2 4.3 5.0
##  [811] 3.3 6.3 4.6 3.1 2.5 4.8 4.7 4.1 4.4 3.2 5.3 4.1 3.1 5.9 5.4 3.8 5.4 4.4
##  [829] 3.5 3.4 5.7 3.8 5.5 3.1 3.9 5.7 4.9 5.0 4.3 4.8 3.4 3.8 3.6 5.3 5.8 4.1
##  [847] 3.8 5.2 3.8 1.9 5.0 4.3 5.1 3.8 4.9 2.2 5.5 2.9 4.3 3.7 3.7 3.9 4.4 3.4
##  [865] 6.2 5.3 4.9 3.7 4.4 5.3 5.3 5.4 2.8 5.3 5.0 4.0 3.9 5.1 3.0 4.8 4.5 4.4
##  [883] 4.8 4.7 3.2 2.6 6.1 3.9 3.2 5.0 4.8 3.4 5.8 3.6 4.9 5.4 5.5 4.8 3.6 4.5
##  [901] 4.4 3.6 3.8 2.8 4.3 4.8 3.5 3.4 4.0 3.7 3.9 5.4 4.9 5.3 3.4 3.4 4.7 4.9
##  [919] 5.9 3.4 5.0 4.2 3.7 4.1 4.5 5.5 4.9 3.7 4.8 4.8 4.2 4.8 4.5 5.5 3.5 5.0
##  [937] 3.3 2.8 4.0 3.6 3.8 3.7 6.8 4.5 3.6 5.4 5.0 5.7 5.2 5.2 3.6 4.0 6.4 5.0
##  [955] 5.1 4.0 3.1 5.1 4.9 4.4 4.4 4.4 4.8 4.3 5.2 3.1 6.4 2.6 5.8 3.5 5.8 4.1
##  [973] 4.6 1.9 2.8 4.0 2.8 4.1 5.1 4.1 3.4 5.9 3.2 5.2 4.3 4.4 4.4 4.6 4.3 6.5
##  [991] 4.3 4.3 5.1 4.5 5.1 3.4 5.1 4.7 4.7 4.5
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
##   2.8   6.3
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
##    [1] 3.7 4.5 3.2 2.9 3.2 5.1 3.7 6.1 5.9 5.9 4.2 5.0 4.1 2.6 4.9 4.4 4.1 3.7
##   [19] 4.4 4.4 4.2 5.6 5.5 4.6 5.3 4.9 5.1 3.8 3.3 3.7 4.0 6.0 5.1 3.9 4.1 5.7
##   [37] 4.9 4.9 4.5 4.8 3.6 5.3 4.4 4.7 4.2 4.9 4.0 3.5 4.4 5.0 4.4 5.3 4.6 4.3
##   [55] 3.9 5.3 6.1 4.6 5.1 5.4 3.3 4.0 4.8 4.3 4.3 5.6 4.0 3.4 4.8 5.4 4.3 5.4
##   [73] 3.7 3.7 4.8 5.9 3.8 3.6 4.1 4.6 4.2 4.9 3.8 4.6 5.5 5.4 5.8 4.4 5.3 2.9
##   [91] 5.2 4.4 4.2 3.8 2.5 5.0 4.8 4.5 3.4 7.1 4.4 5.7 3.9 6.1 4.4 3.5 5.7 5.2
##  [109] 4.4 5.0 5.4 4.3 5.0 4.5 3.7 4.3 4.3 5.2 4.4 3.4 4.2 3.8 3.8 4.6 4.4 5.4
##  [127] 4.1 3.8 3.5 5.2 5.4 4.8 5.0 4.4 5.2 3.2 2.6 4.3 4.5 4.9 6.6 4.8 4.0 4.8
##  [145] 4.4 5.0 3.9 2.9 4.2 4.9 4.8 6.3 2.7 4.9 4.8 4.4 3.0 4.8 4.2 4.1 5.9 5.1
##  [163] 6.0 3.9 4.8 5.7 5.4 4.3 4.2 3.8 4.5 3.6 5.8 3.7 2.7 4.5 3.8 4.5 4.5 3.8
##  [181] 4.8 4.8 4.3 3.5 4.1 4.8 4.6 4.9 5.4 5.3 5.4 4.0 6.5 4.4 4.3 4.0 4.0 5.2
##  [199] 5.4 3.9 3.2 5.6 2.9 4.9 5.2 3.8 5.2 3.1 4.9 3.7 4.8 4.8 4.8 4.1 5.7 4.3
##  [217] 3.2 2.7 4.9 5.6 4.3 3.5 5.8 3.9 4.3 3.9 5.1 4.1 6.0 6.1 5.4 5.8 3.4 4.8
##  [235] 4.9 4.4 4.9 5.0 5.7 3.2 5.4 3.9 4.9 4.5 4.0 4.1 5.8 3.8 4.4 3.6 3.9 3.3
##  [253] 4.7 6.1 3.6 3.0 4.1 4.7 3.7 4.8 4.4 5.1 2.6 5.4 4.9 2.9 3.7 4.9 5.3 4.6
##  [271] 4.3 5.5 3.9 4.0 5.3 6.0 3.3 4.1 3.7 6.0 4.0 5.9 4.6 4.0 3.3 5.5 4.0 4.6
##  [289] 5.7 4.7 4.6 3.7 3.8 5.0 6.1 3.9 6.2 3.4 4.2 5.4 4.0 4.2 4.2 4.1 4.2 4.7
##  [307] 4.8 4.1 4.2 5.5 3.1 4.7 5.4 3.8 4.8 3.7 5.7 5.4 3.1 3.5 4.0 4.1 5.6 3.5
##  [325] 4.5 4.9 4.9 3.7 4.6 5.8 5.4 3.9 4.8 5.3 4.7 4.2 3.3 5.5 4.9 6.4 4.8 5.3
##  [343] 4.5 4.3 3.7 4.6 4.6 3.3 3.5 3.6 4.9 4.2 4.8 4.3 3.5 2.8 4.7 3.8 3.4 3.4
##  [361] 5.2 4.5 6.1 3.6 3.5 3.3 3.4 3.8 4.5 6.8 3.9 4.1 4.0 2.5 3.1 3.2 3.6 4.6
##  [379] 6.6 5.9 4.7 4.6 4.4 5.5 3.0 3.9 4.2 4.8 3.8 6.5 4.7 4.6 3.6 2.9 3.9 3.7
##  [397] 4.7 6.1 5.6 6.3 3.8 4.8 3.9 4.1 3.9 5.1 3.7 5.2 3.3 2.2 4.0 5.1 4.0 4.8
##  [415] 4.6 6.5 5.1 3.7 5.1 6.4 5.0 4.5 3.5 4.9 6.0 3.6 4.5 4.3 4.0 4.0 4.6 6.3
##  [433] 4.6 4.8 4.5 5.4 4.7 3.5 4.7 5.5 6.1 5.7 3.0 3.4 5.7 4.5 3.0 3.4 2.2 5.5
##  [451] 4.9 3.9 5.0 4.6 5.7 5.2 4.3 3.9 4.6 3.5 5.9 5.6 5.0 4.0 4.0 3.3 5.5 3.6
##  [469] 4.3 4.9 2.3 5.2 4.8 4.2 3.0 4.9 4.8 4.9 3.6 4.6 5.1 2.5 4.1 5.4 3.5 4.8
##  [487] 5.1 4.2 3.5 4.0 3.8 4.8 3.9 3.7 4.7 4.6 6.3 4.1 3.8 3.0 3.6 3.3 3.7 5.0
##  [505] 4.0 4.7 3.9 3.5 4.2 4.4 3.8 4.5 3.7 3.6 3.1 4.1 4.5 4.8 3.7 4.8 5.4 5.0
##  [523] 4.1 5.3 2.7 6.1 4.6 4.8 3.5 4.4 5.4 4.5 4.4 4.5 5.6 5.7 4.4 3.9 4.4 6.5
##  [541] 4.0 3.9 4.2 5.7 3.4 5.6 4.6 4.8 5.0 3.7 6.3 4.0 5.4 4.4 5.5 4.1 4.6 5.8
##  [559] 4.9 4.8 3.3 5.3 3.4 4.3 5.3 4.4 6.1 5.4 4.4 4.7 3.1 4.2 4.3 4.1 5.1 5.9
##  [577] 3.0 3.9 3.2 5.8 5.4 4.8 6.7 3.6 5.4 3.7 5.4 4.3 5.3 4.1 3.7 3.8 5.0 3.8
##  [595] 4.8 3.1 4.8 4.3 5.7 5.5 3.0 5.1 5.2 3.6 4.6 4.0 5.7 4.7 3.5 3.7 5.0 4.7
##  [613] 3.8 4.6 3.9 6.0 4.8 3.3 4.4 2.9 3.5 3.7 6.0 4.5 4.2 5.3 4.5 3.8 5.5 5.9
##  [631] 3.2 2.7 5.6 6.0 6.2 4.4 2.3 4.2 3.3 3.4 5.5 4.1 6.5 1.8 5.3 4.1 4.5 3.6
##  [649] 5.0 5.9 3.5 4.1 3.1 5.1 4.1 3.9 4.9 3.9 4.6 4.3 4.8 3.7 5.2 4.5 4.4 3.9
##  [667] 5.1 3.6 4.4 4.5 3.5 5.3 4.5 5.5 4.9 5.6 4.6 4.8 3.3 3.0 4.8 5.3 4.9 4.4
##  [685] 4.1 3.8 4.4 3.5 4.4 6.1 4.0 2.9 3.8 4.5 2.4 4.5 4.9 3.7 5.1 4.1 3.8 4.5
##  [703] 4.1 5.0 3.8 4.1 4.9 5.2 5.2 4.3 4.9 4.7 3.5 5.2 5.4 3.4 4.3 4.6 3.8 4.9
##  [721] 4.5 5.1 5.6 5.4 3.4 5.2 5.3 5.6 2.7 2.6 4.9 4.8 4.8 2.5 5.1 4.6 4.0 4.3
##  [739] 3.4 4.9 3.6 3.4 4.4 4.8 3.7 4.0 5.9 3.7 3.0 5.5 3.2 5.1 4.4 4.8 4.7 3.3
##  [757] 3.8 4.2 4.7 5.6 2.0 3.5 3.9 4.9 6.8 4.4 3.0 4.4 3.9 6.1 4.9 4.8 3.9 4.7
##  [775] 4.2 4.0 4.6 5.2 5.1 3.8 3.1 5.0 3.8 4.8 3.5 4.3 3.6 2.8 3.9 5.6 4.7 3.3
##  [793] 3.1 4.9 4.9 5.4 5.5 3.8 5.0 5.3 5.2 6.1 5.9 4.3 5.1 4.5 4.1 4.5 3.7 5.6
##  [811] 4.6 6.3 3.3 6.1 4.4 4.6 2.7 5.6 6.2 6.6 4.0 3.8 3.9 5.5 3.3 5.2 5.1 2.7
##  [829] 4.9 4.8 5.0 4.7 4.6 4.4 4.8 4.5 6.5 5.9 3.8 3.4 3.8 4.3 6.5 3.4 4.9 3.4
##  [847] 4.1 5.2 4.5 4.9 4.4 3.7 6.2 4.6 3.9 4.9 2.0 5.9 3.8 4.2 2.9 5.4 4.1 4.9
##  [865] 5.1 2.9 5.3 3.4 4.7 2.6 2.4 6.5 4.7 5.5 4.3 4.4 4.5 4.9 4.9 4.8 3.5 5.6
##  [883] 5.3 4.8 4.8 5.1 5.3 4.7 3.1 4.4 4.8 4.1 2.8 3.9 5.0 5.5 4.8 5.5 4.3 4.4
##  [901] 5.4 4.3 6.5 4.4 4.9 3.9 6.1 2.5 5.6 3.7 4.1 4.4 3.4 5.0 3.9 6.4 5.2 4.4
##  [919] 2.7 4.3 4.7 5.4 5.0 4.5 5.3 5.3 4.7 5.4 4.0 5.5 3.9 4.1 5.3 3.3 3.8 3.8
##  [937] 5.3 5.5 4.1 5.2 4.9 3.0 4.0 2.4 3.4 4.4 5.2 3.3 4.3 5.2 3.3 3.6 4.4 3.8
##  [955] 5.1 4.4 4.0 4.9 4.0 4.9 2.6 4.7 5.3 4.2 5.3 4.2 4.6 4.3 5.0 2.9 3.9 3.0
##  [973] 5.1 4.5 3.9 4.3 3.8 5.8 4.5 3.8 5.3 5.2 4.5 4.4 4.5 7.1 4.8 3.6 4.7 4.8
##  [991] 4.4 4.3 6.0 4.1 3.7 2.9 4.4 5.3 4.5 4.0
## 
## $func.thetastar
## [1] -0.0254
## 
## $jack.boot.val
##  [1]  0.53848580  0.28813559  0.21684492  0.13352273  0.04753247 -0.06406250
##  [7] -0.24467456 -0.29390244 -0.40708447 -0.53901099
## 
## $jack.boot.se
## [1] 0.9656953
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
##    [1] 5.5 2.9 4.0 3.7 3.8 4.9 4.0 3.9 3.9 5.9 5.2 5.5 4.4 4.5 5.7 4.3 4.5 4.0
##   [19] 4.2 4.7 4.2 3.6 5.0 5.1 5.2 5.7 4.8 5.3 5.4 5.4 4.3 3.6 4.8 4.1 4.3 4.4
##   [37] 3.8 5.5 5.7 5.6 4.9 3.6 4.5 4.7 3.0 3.3 3.4 5.2 4.0 4.8 4.2 4.5 5.1 4.3
##   [55] 4.4 4.2 5.1 4.9 4.1 4.1 4.5 5.5 4.4 3.8 4.0 3.9 4.0 6.3 5.2 3.7 6.4 4.1
##   [73] 3.3 4.1 4.8 4.5 5.0 5.3 4.4 4.9 4.4 5.7 4.4 4.0 5.6 3.8 6.2 4.7 4.5 3.2
##   [91] 5.1 4.5 4.5 4.2 4.4 2.5 3.0 4.0 5.8 4.0 4.3 4.2 5.1 4.6 4.8 3.2 5.4 5.0
##  [109] 4.8 3.6 4.5 3.3 4.3 4.8 5.1 4.8 4.8 3.3 6.6 4.1 3.8 5.6 4.4 3.5 5.5 5.1
##  [127] 5.5 4.7 5.3 4.4 4.6 4.8 3.7 5.3 4.2 4.6 6.1 5.5 3.2 4.1 5.2 4.7 6.3 3.1
##  [145] 4.6 4.4 3.9 2.8 4.6 3.4 4.0 5.0 4.3 3.0 4.5 5.2 3.2 4.6 4.5 4.3 3.5 4.3
##  [163] 4.1 4.8 6.1 5.7 3.7 4.9 4.1 3.1 5.8 5.7 3.4 4.4 4.9 4.6 5.2 4.0 3.8 4.0
##  [181] 4.3 3.5 4.8 4.5 5.1 4.1 5.7 3.6 5.5 4.1 4.3 5.8 4.5 4.2 4.5 4.9 4.5 4.9
##  [199] 5.2 3.9 2.1 3.6 4.7 5.2 5.1 6.4 3.0 5.5 5.8 5.2 5.3 5.4 6.9 4.0 5.5 3.7
##  [217] 4.6 5.2 4.2 5.7 5.1 2.2 4.7 3.2 4.3 5.6 5.4 5.1 5.6 3.7 4.5 3.9 4.4 5.3
##  [235] 4.5 4.9 4.0 4.9 5.3 4.8 4.7 5.1 4.4 5.4 5.7 4.3 5.4 2.5 5.1 3.9 3.8 5.0
##  [253] 4.0 4.0 6.0 4.5 4.2 3.6 3.9 5.3 4.6 6.0 4.9 4.8 4.1 6.7 5.2 3.7 5.5 4.3
##  [271] 4.8 4.2 1.9 4.6 3.7 4.5 6.4 4.2 5.2 3.3 4.2 4.8 3.6 4.5 3.9 4.0 3.4 4.6
##  [289] 4.1 4.4 5.1 6.1 4.5 5.0 2.6 6.0 3.6 4.2 4.4 3.9 3.4 5.3 3.9 5.8 3.0 3.7
##  [307] 5.7 3.6 4.0 3.8 2.9 3.4 5.7 5.0 5.1 4.8 3.4 4.4 5.9 4.3 5.8 4.1 5.2 5.2
##  [325] 5.6 4.9 4.4 3.2 6.3 5.7 3.7 3.5 4.3 3.7 6.1 6.2 5.1 4.1 4.8 3.6 3.0 5.8
##  [343] 5.6 4.2 5.0 4.3 4.7 4.1 5.3 5.4 2.9 4.2 4.4 4.1 4.5 4.2 3.5 4.3 4.3 5.0
##  [361] 4.2 5.5 4.9 4.4 3.8 3.2 4.6 4.7 5.2 4.3 4.9 3.7 3.6 5.0 5.4 3.5 4.5 4.5
##  [379] 4.9 3.6 4.2 3.1 5.4 4.6 5.0 4.9 3.3 4.3 4.4 3.9 4.2 6.0 4.6 2.5 4.3 3.6
##  [397] 5.4 4.9 4.9 4.5 3.1 2.4 5.3 4.1 5.1 5.3 4.6 5.3 3.6 4.7 6.0 4.5 5.3 3.2
##  [415] 4.5 4.1 5.4 5.0 4.9 5.5 5.2 5.4 5.2 5.1 5.8 3.1 5.5 4.3 3.5 2.9 3.9 4.6
##  [433] 3.8 2.7 2.8 4.6 4.0 4.0 5.7 3.9 4.4 4.1 4.2 3.8 5.2 3.4 3.5 4.8 3.8 4.9
##  [451] 5.0 5.7 4.4 3.9 5.4 4.7 3.4 6.1 3.6 4.2 4.3 3.9 5.2 5.0 5.9 4.4 4.2 3.3
##  [469] 3.9 5.4 3.8 4.2 3.6 3.0 2.1 3.9 4.6 3.8 5.7 5.7 3.4 3.1 4.3 3.0 4.6 4.5
##  [487] 5.2 4.8 4.3 4.3 2.7 3.5 3.3 4.8 5.5 4.7 4.6 4.7 1.9 3.9 5.6 3.9 4.3 3.9
##  [505] 4.3 3.2 3.9 4.7 5.0 4.3 4.4 5.2 5.9 4.4 5.4 5.1 6.0 3.3 2.7 4.5 6.3 4.5
##  [523] 3.9 4.7 3.2 5.6 2.9 3.9 4.0 4.2 5.1 4.0 4.5 5.5 4.8 4.3 5.1 3.8 3.4 4.8
##  [541] 5.4 5.1 4.4 6.2 4.5 4.9 4.9 4.9 2.9 2.9 4.4 5.0 5.0 4.6 5.0 3.7 4.5 4.0
##  [559] 3.6 3.8 5.1 4.1 3.7 4.8 5.9 5.0 4.7 6.4 4.8 5.8 4.0 3.9 5.2 3.2 4.7 5.8
##  [577] 4.7 3.9 4.8 5.0 5.1 5.1 4.9 3.5 4.2 3.2 2.7 4.5 3.6 4.4 3.7 6.1 4.7 5.2
##  [595] 5.2 4.8 3.3 2.9 4.7 5.9 4.6 3.4 3.8 3.1 3.8 5.5 4.3 5.1 4.4 3.9 5.1 4.0
##  [613] 3.6 4.1 4.1 4.1 6.1 5.1 2.5 4.6 4.2 4.8 4.4 4.1 6.2 4.7 6.0 3.7 4.9 2.9
##  [631] 5.3 5.1 4.5 5.1 4.4 5.6 5.9 5.3 3.9 5.6 5.9 3.7 5.2 5.0 3.3 5.0 5.6 4.6
##  [649] 6.1 3.2 4.2 4.2 5.3 5.2 4.7 4.9 5.4 5.4 4.9 5.7 6.0 6.4 4.4 4.2 5.7 4.5
##  [667] 4.0 4.1 5.4 4.5 5.4 3.6 4.9 5.0 4.0 5.0 5.2 3.2 5.7 5.0 5.4 3.5 3.1 4.6
##  [685] 4.6 4.0 5.1 4.1 5.4 5.6 5.0 4.8 4.4 4.6 3.6 3.4 3.8 6.0 4.4 3.4 5.3 4.6
##  [703] 4.2 4.1 4.3 3.6 4.3 4.9 4.9 3.9 5.6 5.1 3.9 4.0 3.7 6.0 3.6 5.6 3.9 5.6
##  [721] 3.8 6.8 3.3 5.0 4.6 4.7 5.1 5.0 3.9 4.0 4.8 4.7 2.9 5.3 4.1 3.7 6.1 4.6
##  [739] 5.6 5.2 6.2 4.9 3.2 5.6 5.4 3.5 4.5 4.8 4.1 4.8 4.3 4.6 4.3 6.8 4.4 6.3
##  [757] 5.6 2.7 4.7 2.2 5.0 2.9 5.1 5.0 5.4 5.6 4.6 4.1 4.7 6.3 3.9 6.0 4.9 4.9
##  [775] 5.3 4.7 3.7 3.8 5.0 4.3 3.7 5.3 4.3 4.7 4.7 5.7 3.6 4.4 4.2 4.9 6.7 3.6
##  [793] 5.1 4.3 2.6 5.0 4.4 4.6 2.2 4.2 5.3 4.6 5.3 5.2 4.8 4.5 6.1 4.7 4.5 5.0
##  [811] 5.6 4.3 5.2 5.1 4.0 4.9 3.8 4.6 3.3 5.2 4.9 5.1 5.4 4.6 4.3 6.0 3.2 5.5
##  [829] 4.6 4.7 4.4 4.2 5.2 3.3 3.9 4.9 5.7 4.6 3.7 3.9 5.0 3.2 4.8 4.2 3.1 3.5
##  [847] 3.0 5.6 5.9 5.3 4.7 4.8 4.6 4.3 4.3 4.0 4.9 4.7 3.6 6.0 5.5 4.3 6.0 5.3
##  [865] 3.4 5.1 4.0 4.7 4.1 3.9 2.9 5.8 3.8 6.2 5.3 3.4 4.9 3.7 5.0 4.2 5.1 5.6
##  [883] 5.9 2.9 4.5 3.4 4.8 3.6 4.3 4.1 5.5 4.2 4.6 5.6 3.5 5.1 5.0 4.2 3.3 4.5
##  [901] 4.3 2.7 5.6 4.1 3.6 3.7 4.1 4.3 4.7 5.1 3.0 3.3 4.7 4.5 5.4 4.2 3.5 4.4
##  [919] 3.9 3.5 2.8 4.8 6.0 5.1 3.1 3.9 5.0 4.5 4.5 5.3 3.7 3.7 3.9 3.3 4.7 4.9
##  [937] 4.9 3.6 3.7 4.8 4.0 5.5 3.4 2.8 4.3 2.8 5.1 4.0 4.6 4.3 4.1 3.9 4.1 5.4
##  [955] 5.4 2.7 4.7 4.6 4.7 4.4 4.0 4.3 4.5 4.6 4.7 5.8 4.8 5.9 4.5 4.8 4.7 5.6
##  [973] 5.7 5.1 5.6 3.9 3.6 5.8 4.3 5.1 4.4 4.8 4.3 4.8 2.5 3.8 3.7 5.2 4.3 5.2
##  [991] 4.1 4.6 5.0 3.9 5.3 4.9 3.2 4.7 3.9 5.3
## 
## $func.thetastar
##   72% 
## 5.028 
## 
## $jack.boot.val
##  [1] 5.5 5.4 5.3 5.1 5.1 5.0 4.8 4.8 4.6 4.6
## 
## $jack.boot.se
## [1] 0.9079648
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
## [1] 0.746142
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
##   2.778576   4.013860 
##  (1.175683) (1.861253)
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
## [1]  0.3120231 -0.2896893  2.1359834  0.5680727  0.7030460 -0.3047134
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
##    [1]  0.389596183  1.440697523  1.159667239  0.406957904  1.371970372
##    [6] -0.168443671  1.502516295 -0.166420200  0.228184266  0.554994024
##   [11]  0.703987181 -0.034624135  0.485107258  0.269884081  0.123738856
##   [16]  1.244458141  0.411775523 -0.045277938  0.877624797  0.603646506
##   [21]  0.667328431  0.553839173  0.448366774  0.712669443  0.514987537
##   [26]  0.559947046  1.451221612  1.447092192  0.189436743  0.831615827
##   [31]  0.826222651 -0.074748031  1.109281346 -0.151255196  0.328829844
##   [36] -0.127176994 -0.529424006  0.343668450  0.488486576 -0.347007842
##   [41] -0.123771474  0.359910406  0.287444580  1.292818375  1.245104715
##   [46] -0.789751276  0.837875453  0.825595294  0.025888409  0.870887795
##   [51]  0.592197264  0.033613879  1.157683112  0.467783735  0.315664444
##   [56]  0.859444432  1.091974630  1.425527856  0.359515581 -0.047029956
##   [61]  0.545934969  0.116046157  1.545794304  0.837102245  0.434421557
##   [66]  0.416934023  0.129518986  0.610769770  1.155947128  0.740683119
##   [71]  1.157496950  0.567878727  1.531910174  0.818854537 -0.003406912
##   [76]  0.357857847  0.818038030  0.809967270  0.950818205  0.879806104
##   [81]  0.227300963  0.441777026  0.671410365  0.797532782  0.356915500
##   [86]  0.740683119  0.841517639  0.146414770  0.421141184  0.769645785
##   [91]  0.393676143  0.275231202  0.588697847  1.693577488 -0.604573880
##   [96]  0.149935414  0.592325073  0.219677566  0.196351009  0.161201120
##  [101]  0.771336502  1.418512593  1.799529916  1.508446820  0.815692163
##  [106]  0.376041069  0.590768131  0.711143112  1.363219071 -0.291812414
##  [111]  0.694342207  0.621051330  0.970346867 -0.174284974  0.340308386
##  [116]  0.778927606  0.790052942  1.692256116  0.035530171  0.249758977
##  [121]  0.463322516  0.386295401  0.739795401  0.468465010  1.624121588
##  [126]  0.321910619  1.788575403  1.391357320  0.893857580  0.439376833
##  [131]  0.838538042  0.357725472  0.484210957 -0.615036850  1.431717724
##  [136]  0.119177203  1.067866093  0.677091936  0.692883050  0.102783646
##  [141]  0.349901959  0.731827561  1.818828804  0.324694053  0.042433201
##  [146]  0.762566243  0.927423030  0.999123505  0.672754453  0.771723367
##  [151]  0.074620172  0.684085232  0.482455295 -0.105898080  0.980895214
##  [156]  1.778465748  1.229327855  1.128500941  0.357576918  0.088892262
##  [161]  0.259301247  0.915521705  0.229087728 -0.462831266  0.285708686
##  [166]  0.456286177  1.306971599  0.597639179  1.090108716 -0.416506198
##  [171]  0.252265677  0.698738463  0.493540100  1.052659593  1.070087234
##  [176]  1.090025055  1.046944058  0.584335190  0.387188121  0.053933949
##  [181] -0.278898482  0.414066797  1.034577634  1.266252062  0.522053834
##  [186]  0.732852901  0.750584278  0.339479551  0.315119834  0.542685302
##  [191]  1.474801033  0.429215823  0.744657842  1.968689120  0.030891416
##  [196]  0.421325473  0.709654247  0.402265848  1.268493676  1.181376008
##  [201] -0.073032019  0.587012651  1.414148137  1.260590395 -0.060309385
##  [206]  0.604199790  0.536638676  0.235738377  0.560199511  0.913510052
##  [211]  0.681444359  0.763762295  0.055360084  1.122221138  0.977575547
##  [216]  0.448383471  0.263315338  0.997979293 -0.769262942 -0.055192718
##  [221]  0.511585298  0.203877559  0.727420462  1.849751743  0.244420845
##  [226]  1.005521738  0.745323634 -0.057340302  0.339223001  0.399852810
##  [231]  0.532039701  0.475825325 -0.414740139  0.147153651  0.522292887
##  [236]  0.837004002  1.132678986  0.404232855  1.191440210  1.160890867
##  [241]  1.143499832  0.439142941  1.902510259  1.467550133  0.216746384
##  [246]  0.608091961  0.842307318  1.708973241  2.346909023  0.845376029
##  [251]  1.147388455  0.554933790  0.692236924  0.569287151  0.729396079
##  [256] -0.015363266  1.271097064  0.516274754  0.858535061  0.920485234
##  [261]  0.844755479  1.838825978  2.037843380  0.215787477  0.501398649
##  [266]  0.615756658  1.415389273  0.431655657  1.597865721 -0.104632234
##  [271] -0.243332626  0.539844441  1.067277622 -0.080435500  0.597655628
##  [276]  0.870684652 -0.416300489  0.765134517  1.626020897  0.688440930
##  [281]  0.569777229 -0.281146982 -0.245523909  0.436860388  0.934982348
##  [286]  0.664197796  0.271997622 -0.635552295  0.217297593  0.672450108
##  [291] -0.032776982  0.306182108  0.879122937  0.170734773  0.887465636
##  [296]  1.206941480 -0.609859571  0.752884644  0.464802726  0.760828076
##  [301]  0.080938994  1.480488871  0.475176946  1.844709040  0.092085739
##  [306]  0.735974919  0.870118226  0.711847223  1.026166584  0.396848150
##  [311]  0.700786169  0.680458647  0.461730125  0.505153778  0.678679212
##  [316]  1.629057822  0.446790088  0.208502263 -0.013412206  0.503917529
##  [321] -0.583867851  0.322717058  0.248267018  0.407247712  1.141073372
##  [326]  0.974842309  1.930259931 -0.034896878  0.940941785  0.086651112
##  [331] -0.050908294  0.163679659  0.328407019  0.960468115  1.181515559
##  [336]  0.789387854  0.458529598  0.821339360  0.359053825  0.157799724
##  [341] -0.236008998  0.647647511  0.437019534  1.334761300  0.650800083
##  [346]  1.529018707  1.487090558 -0.195449274  0.726357460  1.326380419
##  [351]  0.947230201  0.508511967  1.378879635  0.097795998  1.244169336
##  [356]  0.966891302  0.038101556  1.095722844  0.736288878  1.311072600
##  [361] -0.197091589  1.609286038  1.265401925  0.810317946  0.715731641
##  [366]  1.044756868  0.021884766  1.157681187  0.885046845  0.646047392
##  [371]  1.520224609  1.169886753  0.185455698  0.659988201  1.256737127
##  [376]  1.237521965  1.045561230  1.094217544  0.593875808  1.124048110
##  [381]  0.419081081 -0.020414039  1.731274092  0.753012432 -0.338339953
##  [386]  0.166447657  0.326198270  1.517780698  1.155849310  0.782010593
##  [391]  0.992049284  1.376021967  0.812948452  0.192728610  0.055259560
##  [396]  1.348699764  0.705515334 -0.079959275  1.222214292  1.087264855
##  [401]  0.137326782  0.807631326  1.104077867  0.142814000  1.458975848
##  [406]  0.038205303  0.882207500  0.778213793  0.924388266  0.840731748
##  [411]  0.830286249  0.550361939  0.652022711  0.353791254  0.799373544
##  [416]  0.457560923 -0.382324255  1.134963252  0.692883050  0.554250703
##  [421]  0.271855691  0.448563613  0.937628103  0.798603172  1.049478270
##  [426]  0.381925587  0.175431602  0.360270503 -0.308227335  1.045669520
##  [431] -0.544849992  0.478867845  0.697853497  0.563653183  1.765160554
##  [436]  0.607530538  0.635629914  0.469897848 -0.114319166  1.093269361
##  [441]  1.840076286  0.947250889 -0.189461955 -0.833471157  1.059606514
##  [446]  0.803670713  1.806238539 -0.170667650  1.782308889  0.729780349
##  [451]  0.885501911  1.049173368  1.079707871  1.593255803  0.535345717
##  [456]  0.845230213  0.519668783  0.678635288  0.265253958  0.344556142
##  [461]  0.476435358  0.904422285  0.332516063  0.292060681  0.362200647
##  [466]  1.977773373  1.905856864  0.310326270  0.931275977  0.471997674
##  [471]  0.368796605  0.677630630  0.757103985  0.301205377  0.947170369
##  [476] -0.432831529  0.711107501  1.114622679  0.537741454  0.084806807
##  [481]  0.735733452  0.664961087  1.133771166  1.027925603  1.109439967
##  [486]  0.087405618  0.769673638  1.247633952  0.979135467 -0.097407236
##  [491]  1.912656098  0.243672489  1.732605757  0.928996933  1.648604135
##  [496]  1.283078189  0.140149419  0.617013326  0.661832114  0.750857951
##  [501]  1.037360460  1.080475194  1.223483825 -0.184562355 -0.006378467
##  [506]  1.029039424  1.369154815  1.071951772  1.539316455  0.979656001
##  [511]  2.050784414  0.960993167  1.415225335  0.501748520  0.910508200
##  [516]  0.948652916  0.473937972  0.223878115  1.683297997  1.708783068
##  [521]  0.220387926  1.191408527  0.755051030  0.657397079  0.634531067
##  [526]  1.223232161 -0.113000653  1.372771844  0.573240936  0.448814009
##  [531]  1.226819232  0.604221869 -0.185346531 -0.157453390  0.645189942
##  [536]  0.606467414  0.880460684  0.489713806  1.265401925  1.039540476
##  [541] -0.147175557  0.431496278 -0.089064423  0.770136272  0.820083933
##  [546]  0.260009551  0.475587460  0.604061935  0.746142010  0.519901579
##  [551]  0.793716897  0.518766757  1.346452556  0.516020743  0.281523371
##  [556]  0.099996704  0.702392046  0.242511986  0.462019499  1.104421340
##  [561]  0.302257228 -0.020507974  0.824540571  0.561684355  0.574898831
##  [566]  1.028195564  0.746643769  0.061802411  1.212227009  0.669194561
##  [571]  0.104421059  1.079422012  0.781402847  0.394241112  1.554958802
##  [576]  1.825232695  0.609610995  0.836209013  1.332028199  1.699037490
##  [581]  0.119177203  1.304070915  0.190644902  0.075023666  0.455777653
##  [586]  0.649680835  1.499726802  0.776173193  0.216212922  0.128974571
##  [591]  0.891475779  0.483910869 -0.160845845  0.131086718  1.622863243
##  [596] -0.242164073  1.774877605 -0.188645734  0.263820162  0.857658745
##  [601]  2.217877232  0.646284775  0.316374547  0.246409318  0.914330941
##  [606]  1.220148575  1.137921953  1.080208962  1.562773354  0.435211135
##  [611]  1.594460646 -0.483853838  1.315681806  0.435549524  0.522292887
##  [616]  1.175983956  0.358527594  0.117379791  0.814669467  0.088484840
##  [621]  0.393629352  0.781402847  1.062899768  0.440096079  0.514722756
##  [626]  0.517641898  0.219026496  0.798499428  0.250765913  0.030170114
##  [631]  0.514982271  0.533055112  0.751036565  0.478863671  0.634326182
##  [636]  0.652783102  0.068875206  0.262845482  0.446414919  1.015707333
##  [641]  0.322432246  1.266252062  0.293385698  0.861168779 -0.513747634
##  [646]  0.603707554  1.089135916  1.090645249  0.229942090  0.695205296
##  [651]  0.589926044  1.446968286  0.966891302  0.092087758  0.778554441
##  [656]  1.895605209  0.490842614  1.326317126  1.209291679  0.709365003
##  [661]  0.387129271 -0.017517369  1.885462718  0.519744847  0.459141637
##  [666]  0.580443322  0.810355912  1.272880608  1.743268472  0.289626147
##  [671]  0.738995791 -0.374550480  0.899701821  1.868746539  0.891753377
##  [676] -0.018762641  0.360057806  1.334353553  0.950242116  1.705980351
##  [681]  0.906748583  1.739466483  0.369641916  0.650284381  0.763089899
##  [686] -0.020230909  0.745957250  0.847934423 -0.053065505  0.793230733
##  [691] -1.088163739 -0.255680086 -0.097568796  0.878697486  0.671709353
##  [696]  0.877706858  0.401297747  0.771873067  1.287281639  0.289578361
##  [701]  1.177334447  0.481010390  0.173960056  1.397292145  0.372002093
##  [706]  0.957891940  0.501210152 -0.426826545  0.304139905  1.216585457
##  [711]  1.465525153  1.604755257  1.369907007  0.859564848  0.853167877
##  [716]  1.262966400  0.701867819 -0.161409127 -0.247535235  0.003878402
##  [721]  1.309576752  0.922604428  0.969625582  1.044305284  0.752590287
##  [726]  0.858022204  1.788018125  0.987746519 -0.163624096  1.888588181
##  [731]  0.324371611  0.578875564  0.549852524  1.186903307  0.377136798
##  [736]  1.003075940  0.420159312  0.853446428  0.511952469  0.798603172
##  [741]  0.393148151  0.230032234  0.635795202  0.249734777  1.082426728
##  [746]  0.412474426  0.685302209  1.137500440  1.261491833 -0.001675974
##  [751]  0.501210152  0.945832609  0.973557936  1.064113240  1.028510801
##  [756]  1.868746539  2.223031304  0.699548719  1.160787494  1.232991060
##  [761]  0.217646192  0.395637874  0.328806991  0.048958481  0.306042095
##  [766] -0.028717350  0.110864750  0.553528684  0.406957290  1.167999450
##  [771]  0.848150758  0.855306482  0.654726216 -0.020967690  1.514242144
##  [776]  1.306726414  0.528916184  0.788734491  0.277948060  1.197818590
##  [781]  1.229436645  0.498158150  1.114256684  1.903672401  0.373283375
##  [786]  0.766244882  0.143425198  0.723719775  0.120612468  0.306533830
##  [791]  0.808940547  1.606334705  0.776173193 -0.189246185  1.770719842
##  [796]  1.528576851  0.351960306  0.943262245  1.100018743 -0.027925540
##  [801]  0.460977239  0.697624886  0.646497558  0.995398070  0.499115577
##  [806]  0.481057536  0.769163113  0.361919748  0.873051649  0.751560267
##  [811]  0.684788470  0.258766166  0.922643954  0.258326425  0.247639313
##  [816]  1.213195316  0.701727529  0.438698704  0.598930907  1.281983474
##  [821]  0.763681472  1.788510218  1.044355164  1.063321994 -0.398404271
##  [826]  1.689400141  1.403600418  0.827715871  0.501880427 -0.438378799
##  [831]  1.338402681  0.837929416 -0.014721909  0.999854871 -0.142897528
##  [836]  0.680502884  0.702346238 -0.038235172  0.855306482  0.712515073
##  [841]  0.110864750  0.930399277  1.155033511  0.347461480  0.194193111
##  [846]  0.692456756  0.329611369  1.307335714 -0.553245865 -0.197428181
##  [851] -0.032338095  1.118624946  1.030687197  0.706142628  0.287003543
##  [856]  0.767580869  1.681515995  0.508117672  0.892093921  0.455906753
##  [861]  1.061264666  0.257016244 -0.345997810  0.689324867 -0.091405473
##  [866]  0.710652121  0.656566151  0.386428268  0.390427082  0.492504924
##  [871] -0.121604020  0.242670695 -0.121097670  0.702592448  1.011770057
##  [876]  0.571186802  0.963099827  0.748417191  0.730954834  1.163033704
##  [881]  0.009852670  0.864357788  0.503507448  1.151645279 -0.084899907
##  [886]  1.799626922  0.632279189  0.772022232  0.268184048  0.709127039
##  [891]  0.878205031 -0.064249229  1.172806255  0.095335533  0.282306979
##  [896]  0.493899694 -0.386073162 -0.811300190  0.471636279  1.268168763
##  [901]  0.638333317 -0.282515406  1.528162437  1.097887322  0.231854020
##  [906]  1.643939478  0.497424845  1.056240138  0.343366395  0.295340525
##  [911] -0.199169527  0.695739055  1.226070290  0.836205896  1.066071071
##  [916]  0.241868162  0.297013504  0.894180991  1.002696063  0.285319622
##  [921]  0.229873994  0.215546090  0.370698718  1.471742768  0.661573252
##  [926]  0.726982873  0.610590465  0.852609217  0.502491315  0.841552150
##  [931]  0.073297468  0.905623102 -0.418789117  0.620494550  0.372138899
##  [936]  0.626973878  1.281666744  0.988361023 -0.389806126  0.774245816
##  [941]  0.583878298  0.335772028  1.201066163  0.732264760  1.868746539
##  [946]  1.776579685  0.589884519  2.636074864  0.072540021  0.774117623
##  [951]  0.664049672  1.091185654  1.381194860  1.076539537  0.612372089
##  [956]  1.969583718  0.722391665  0.704367322  1.083674654  0.602232651
##  [961]  0.558473480  0.811042203  0.459255611  0.731038915  0.658734631
##  [966]  1.430923383  0.541091855  0.393734290  1.037606721 -0.066886312
##  [971]  1.650000390  1.475906825  0.433379954  1.626891315  0.834159147
##  [976]  0.282578289  0.209313349  0.503024582  0.234701397  0.762566243
##  [981]  0.678385758  1.407944881  0.164123321  1.029852871  0.650382420
##  [986]  1.582700442  1.358313910 -0.319661263  0.945652502  0.825779503
##  [991]  1.064678196  0.303624219  1.017614912  0.414509903  1.991730371
##  [996]  0.510071729  0.432792345  0.878439916  0.735539313  0.922304082
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
##   0.69225243   0.42405533 
##  (0.13409807) (0.09481845)
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
## [1] -0.05950368  0.86689857  0.98044933 -0.47956987 -0.38740922  0.32125531
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
## [1] -0.0276
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8756969
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
## t1*      4.5 0.009409409   0.9043712
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 2 4 5 6 7 9 
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
## [1] 0.0258
```

```r
se.boot
```

```
## [1] 0.9161297
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

