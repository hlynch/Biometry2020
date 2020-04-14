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
## 1 3 4 5 8 9 
## 2 1 2 1 1 3
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
## [1] 0.0154
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
## [1] 2.753153
```

```r
UL.boot
```

```
## [1] 6.277647
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##   2.5%  97.5% 
## 2.8000 6.2025
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
##    [1] 4.6 4.4 4.9 5.4 4.0 5.6 3.9 2.5 4.5 5.2 4.8 4.0 4.2 4.1 4.4 5.1 5.1 5.8
##   [19] 3.3 6.1 4.3 3.1 5.7 4.0 4.3 4.0 5.3 4.9 4.3 4.3 4.2 6.1 4.4 4.9 4.7 4.2
##   [37] 6.5 4.4 4.5 4.3 4.8 4.6 4.9 3.5 3.8 3.4 4.2 5.3 5.3 6.3 4.9 4.3 5.5 6.2
##   [55] 5.4 4.8 4.8 5.0 4.1 6.3 5.8 6.0 5.5 6.1 5.8 4.3 4.7 2.4 4.5 5.0 4.3 4.4
##   [73] 5.4 4.4 4.4 4.7 4.7 5.6 3.8 5.0 4.8 4.4 5.7 4.5 4.0 5.5 2.7 5.0 4.8 4.9
##   [91] 4.3 5.0 4.9 5.2 5.5 3.6 4.2 4.1 5.2 3.8 4.8 4.3 5.5 3.8 3.0 4.6 5.1 4.8
##  [109] 3.9 5.0 4.4 3.8 5.0 5.7 4.4 5.2 4.0 5.0 4.5 4.5 4.9 3.3 4.2 5.8 4.4 5.9
##  [127] 5.1 4.2 6.8 4.2 5.0 4.5 4.7 2.8 3.6 3.2 4.2 4.8 4.9 3.4 4.4 4.6 3.5 3.8
##  [145] 4.5 4.7 4.9 5.1 3.3 3.0 5.2 6.0 3.1 5.1 4.4 3.2 3.8 5.7 5.0 2.6 4.5 4.6
##  [163] 4.4 4.2 4.5 5.7 4.3 3.4 4.2 5.0 3.4 3.8 3.8 3.3 4.0 5.6 5.3 6.1 5.2 5.0
##  [181] 4.4 3.8 4.8 3.3 5.6 4.4 4.9 3.7 5.4 4.3 5.2 5.6 2.4 3.4 4.6 3.6 5.2 4.4
##  [199] 6.1 5.2 3.6 5.5 5.6 4.3 4.9 4.1 3.6 4.3 3.5 4.5 4.3 5.3 5.1 4.9 3.3 4.4
##  [217] 5.0 4.5 3.3 3.8 5.4 4.7 5.1 4.5 4.7 5.1 4.4 4.6 4.3 3.2 5.3 4.5 4.5 5.7
##  [235] 4.1 3.8 3.2 3.3 4.2 5.7 4.4 4.5 4.6 5.8 4.7 5.0 3.2 3.8 4.3 3.2 1.5 4.8
##  [253] 5.0 4.2 4.2 4.2 4.6 4.5 6.6 4.6 4.0 5.1 6.1 3.3 4.1 4.9 5.0 2.9 5.3 5.9
##  [271] 5.9 4.7 4.5 3.8 6.0 4.3 5.7 4.8 2.5 4.6 3.4 5.7 4.2 4.7 5.4 3.9 5.8 5.1
##  [289] 4.8 4.0 5.7 4.2 2.8 4.1 6.0 5.0 4.1 3.5 5.3 4.9 3.9 5.6 4.4 3.9 4.7 5.2
##  [307] 3.2 4.2 3.5 4.7 3.0 4.4 5.6 4.4 4.5 3.7 5.2 3.4 4.8 3.3 4.8 5.6 4.0 5.3
##  [325] 5.5 6.0 2.4 3.8 5.2 6.2 5.9 5.5 2.9 4.6 2.2 4.0 3.9 4.3 3.0 4.2 5.2 3.7
##  [343] 6.0 5.4 4.6 4.1 4.4 4.7 5.1 3.6 5.5 3.9 5.0 6.0 3.7 4.3 4.7 5.4 4.5 4.8
##  [361] 5.8 3.9 5.2 4.7 3.0 5.3 4.9 4.7 4.5 4.8 4.4 4.4 2.4 3.3 4.8 5.1 6.5 5.3
##  [379] 5.7 4.9 4.7 4.9 5.6 6.6 5.6 4.9 4.5 4.2 4.4 4.4 3.4 6.1 5.7 5.9 4.5 3.2
##  [397] 5.4 4.6 5.4 4.0 4.6 5.1 5.4 4.4 4.4 4.4 4.1 4.4 3.5 4.5 4.1 5.1 3.8 3.7
##  [415] 4.1 3.7 3.4 5.3 4.9 4.9 2.0 5.4 4.3 4.2 3.9 4.1 4.1 4.9 3.6 5.0 2.9 4.4
##  [433] 3.1 4.2 4.6 4.0 5.2 5.1 4.5 5.2 4.6 3.7 4.7 4.2 5.8 5.3 6.6 4.5 5.3 5.3
##  [451] 5.8 5.7 4.4 4.4 5.5 3.9 3.9 5.3 5.1 3.1 5.3 4.9 4.1 3.2 4.1 5.7 3.8 3.9
##  [469] 2.6 4.0 4.2 3.3 5.0 3.5 3.2 4.2 4.2 5.6 5.6 5.2 4.4 4.8 3.4 3.3 5.9 4.5
##  [487] 3.3 3.7 4.7 5.6 5.4 4.9 4.2 4.0 3.8 4.6 3.9 5.5 5.2 4.4 4.6 5.8 5.2 5.3
##  [505] 3.1 5.3 4.6 5.6 4.2 6.4 6.3 3.2 4.4 3.5 5.1 4.1 4.8 3.7 5.9 4.6 3.6 5.0
##  [523] 4.0 3.6 5.8 4.5 4.7 3.8 4.9 4.4 3.9 3.6 4.5 4.4 4.5 4.9 5.1 4.0 4.9 4.0
##  [541] 3.5 5.0 6.4 4.5 3.2 5.1 5.3 3.9 3.4 4.0 4.2 5.3 4.5 5.9 4.3 3.4 4.6 5.3
##  [559] 3.8 5.2 4.2 6.1 5.3 3.4 4.0 2.7 4.8 5.5 5.4 2.4 4.3 5.5 6.2 4.3 6.2 6.3
##  [577] 3.8 5.6 4.7 6.0 4.7 3.4 5.2 4.6 3.9 4.0 5.6 4.4 4.2 4.9 4.8 5.2 3.9 3.7
##  [595] 4.1 5.4 4.0 3.5 5.4 4.6 4.5 5.5 4.5 4.2 3.0 4.1 5.1 4.5 5.5 3.6 4.1 5.2
##  [613] 3.3 3.8 4.1 4.8 5.5 5.4 3.2 4.3 3.6 4.6 3.3 4.6 4.3 4.6 6.6 4.0 3.5 4.0
##  [631] 3.9 4.5 5.5 4.3 4.1 4.6 3.0 3.6 4.0 5.1 5.3 4.0 3.4 4.0 5.2 3.8 5.3 5.2
##  [649] 3.7 5.4 4.6 5.2 5.7 5.7 4.5 5.3 5.5 4.5 6.5 2.6 4.1 4.2 2.1 5.1 4.3 5.1
##  [667] 5.2 3.9 4.7 4.4 5.5 3.5 4.4 6.2 3.6 5.1 3.3 4.1 3.8 4.5 3.7 5.2 5.5 3.7
##  [685] 4.9 3.1 4.7 5.2 4.1 4.6 2.6 5.2 4.0 4.0 3.5 5.3 4.0 3.9 4.5 4.1 4.6 3.2
##  [703] 3.7 4.8 4.0 6.2 4.4 5.1 3.6 3.8 5.8 3.0 3.5 4.5 4.8 5.5 5.1 4.8 3.9 4.1
##  [721] 5.6 5.0 4.2 3.1 3.7 4.8 5.3 4.7 4.6 2.5 6.0 4.3 5.2 4.4 5.7 6.4 4.0 6.2
##  [739] 2.8 2.8 3.0 5.3 3.1 5.2 3.8 6.4 3.7 3.2 4.7 5.6 4.7 4.4 6.2 4.8 4.3 4.3
##  [757] 4.2 5.9 5.4 3.6 4.8 6.0 4.8 3.1 4.6 4.8 3.6 6.0 4.8 5.7 4.9 4.6 4.5 5.0
##  [775] 4.2 3.6 5.1 5.9 4.6 5.1 5.8 5.3 5.4 1.8 5.0 3.3 5.0 4.4 4.2 4.7 3.4 5.0
##  [793] 5.3 6.0 4.9 3.9 3.9 3.9 5.7 4.5 4.3 6.4 3.7 2.4 4.2 4.1 6.3 4.1 3.6 5.3
##  [811] 4.8 3.8 4.2 6.5 3.5 6.4 3.0 4.1 5.2 3.4 5.1 5.5 5.1 4.9 5.5 3.6 3.7 4.5
##  [829] 3.2 3.7 4.9 5.5 4.4 3.6 4.9 4.3 4.4 4.5 3.9 3.3 4.4 5.6 4.8 3.9 5.5 3.8
##  [847] 3.3 4.5 6.2 4.9 4.9 4.4 4.6 4.7 5.4 5.3 5.4 3.8 5.4 3.3 4.3 2.6 5.5 2.8
##  [865] 4.2 4.9 4.9 6.1 3.7 6.0 4.5 4.0 5.5 4.5 4.8 5.8 6.0 5.1 3.9 4.8 3.2 3.4
##  [883] 5.2 3.5 3.4 4.2 2.8 3.5 4.9 4.4 4.0 5.5 4.7 4.1 4.4 5.0 3.7 4.0 4.0 3.6
##  [901] 5.6 5.2 4.5 6.2 3.4 5.6 3.7 4.6 5.4 4.6 3.6 3.5 3.1 3.9 4.7 5.3 5.1 5.7
##  [919] 2.4 5.7 4.7 3.3 4.6 4.9 4.2 4.4 5.2 5.3 4.1 3.8 4.2 3.6 4.3 4.6 4.5 5.5
##  [937] 4.5 3.9 4.6 2.8 5.6 6.5 4.5 5.1 5.1 2.2 4.6 2.9 4.5 4.1 4.1 4.5 3.3 3.2
##  [955] 4.8 4.8 4.1 3.4 5.9 4.1 4.1 6.6 3.2 5.6 5.4 3.5 5.2 4.4 3.1 4.1 2.3 3.7
##  [973] 4.6 3.7 4.2 4.4 4.6 5.5 4.4 4.4 3.3 4.4 3.2 5.0 3.7 5.2 5.2 4.9 4.5 4.4
##  [991] 4.7 5.1 4.6 4.8 4.9 3.7 4.8 4.0 4.4 4.4
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
##    [1] 5.7 3.0 3.9 4.5 3.3 5.4 5.7 4.0 5.7 5.8 4.4 4.9 4.8 5.1 4.0 4.6 3.0 4.2
##   [19] 4.4 5.1 4.4 5.1 3.5 4.7 4.0 3.6 4.5 4.7 5.8 4.7 4.4 5.5 4.0 2.8 3.0 5.4
##   [37] 5.4 4.6 3.7 4.8 3.9 5.3 2.7 5.0 4.5 4.6 3.9 4.6 5.2 4.1 2.8 3.1 4.6 3.7
##   [55] 5.1 4.4 4.6 2.4 3.9 4.1 5.1 4.4 2.8 4.5 4.7 3.0 3.8 4.7 3.4 3.7 4.9 5.8
##   [73] 6.0 4.5 4.5 4.8 3.4 4.8 4.2 3.9 3.1 6.4 3.0 3.5 4.3 3.1 5.1 4.8 3.8 4.9
##   [91] 5.3 3.9 5.6 3.8 6.7 4.7 4.4 5.1 5.3 5.2 4.0 4.1 4.0 5.6 4.3 6.0 3.7 5.4
##  [109] 5.9 5.8 5.1 4.7 2.9 2.8 4.4 5.3 5.0 3.3 5.5 6.1 6.0 3.6 3.7 5.3 4.2 5.2
##  [127] 3.0 4.4 4.3 3.4 3.1 4.3 3.7 6.1 3.1 3.5 4.1 4.1 4.3 4.5 3.1 6.2 4.2 5.1
##  [145] 5.2 3.9 5.1 3.4 4.3 3.0 6.6 5.9 3.7 6.2 3.5 4.5 4.3 5.0 5.0 4.2 4.3 3.1
##  [163] 5.8 5.2 4.4 6.0 6.2 5.3 5.5 3.9 5.0 5.0 4.4 6.2 4.1 5.0 5.0 5.4 5.6 3.7
##  [181] 5.6 6.8 4.0 2.4 6.5 4.8 5.3 4.7 3.8 4.6 3.3 3.8 3.3 4.3 4.7 3.3 5.5 5.1
##  [199] 5.9 4.6 5.9 3.4 5.0 5.6 5.0 5.3 4.4 5.0 3.5 4.2 5.3 4.8 5.0 3.9 4.4 4.6
##  [217] 3.3 4.3 3.8 4.5 4.8 5.9 4.5 5.9 4.2 5.5 4.6 5.2 4.0 4.7 4.8 3.9 6.3 3.7
##  [235] 4.9 5.1 5.1 5.3 6.6 3.9 4.5 3.4 3.8 3.9 3.3 4.6 4.6 3.1 5.1 4.7 3.1 2.9
##  [253] 3.5 5.1 2.2 3.7 4.5 5.7 4.3 5.4 5.4 3.0 6.0 3.1 4.3 5.7 6.4 4.0 4.9 4.0
##  [271] 4.3 3.7 3.7 4.1 5.2 3.3 3.1 4.3 3.7 5.8 2.9 2.6 4.1 4.3 4.8 6.0 5.2 4.9
##  [289] 5.5 3.5 3.9 3.7 6.3 5.4 2.7 4.1 6.2 3.4 5.7 4.3 3.5 4.3 4.5 4.4 5.3 4.5
##  [307] 4.5 4.5 4.5 3.8 4.7 3.8 6.5 3.4 5.7 5.2 3.7 5.5 6.3 4.7 5.0 4.5 5.4 4.4
##  [325] 4.1 4.8 6.0 3.3 5.3 5.4 5.1 4.3 3.1 4.4 3.5 5.4 5.1 4.8 6.3 3.9 4.6 2.0
##  [343] 4.5 4.1 3.8 3.7 4.9 5.7 4.7 5.5 3.1 4.5 4.6 4.6 4.1 3.4 5.4 5.3 3.5 4.8
##  [361] 5.7 4.2 4.3 5.8 4.1 3.9 4.3 3.2 3.1 5.9 3.9 5.1 5.5 3.5 3.8 3.3 3.5 4.5
##  [379] 4.8 3.9 4.6 6.0 5.2 5.9 5.2 5.0 3.5 4.1 5.6 4.4 4.5 5.4 4.7 2.7 4.3 3.1
##  [397] 4.3 3.8 6.0 4.5 4.9 3.6 4.2 6.4 5.1 5.0 4.8 4.5 4.7 4.4 4.4 6.3 3.7 4.6
##  [415] 5.9 5.2 4.0 4.3 4.3 5.1 4.6 4.3 6.2 3.7 5.1 5.8 3.7 4.9 4.4 4.3 2.8 3.9
##  [433] 5.9 5.2 4.3 6.1 4.8 5.6 4.0 6.4 5.3 4.5 5.3 4.1 4.9 3.8 4.9 5.4 3.7 4.9
##  [451] 2.6 4.4 2.7 4.0 3.9 4.2 5.5 4.2 4.8 3.2 3.6 3.6 2.8 4.1 4.6 4.7 5.3 3.2
##  [469] 3.9 6.4 3.6 4.4 4.5 4.4 5.7 4.4 4.5 4.1 3.2 5.0 4.8 3.0 3.8 3.9 3.8 5.7
##  [487] 3.3 3.5 3.3 4.3 4.7 5.6 4.2 4.6 2.9 5.6 4.8 3.6 4.9 6.4 4.0 2.5 3.7 5.0
##  [505] 4.7 5.2 4.1 4.4 4.5 4.0 4.1 3.5 3.6 5.2 4.0 5.5 4.7 4.3 3.1 5.1 4.5 4.7
##  [523] 4.3 3.9 5.0 5.5 4.7 2.8 5.9 3.7 4.1 4.7 5.2 4.5 5.8 5.3 3.7 3.5 3.8 7.0
##  [541] 5.3 5.0 6.1 6.0 4.4 3.5 3.7 6.8 2.5 4.3 3.1 3.9 5.3 3.5 4.3 6.1 5.5 5.0
##  [559] 6.2 4.4 5.5 5.1 3.4 4.7 3.6 5.2 4.5 5.8 5.9 6.4 5.1 5.0 3.6 5.7 3.1 3.6
##  [577] 3.8 4.8 3.7 5.1 3.0 3.2 4.3 4.3 3.7 4.8 4.8 4.1 4.3 3.1 6.3 4.3 3.7 4.5
##  [595] 4.7 4.6 4.4 4.7 5.9 5.1 4.8 4.1 5.8 3.9 4.0 4.0 3.4 6.5 5.5 3.6 3.9 6.3
##  [613] 4.5 4.3 5.7 3.8 4.4 3.4 4.0 4.4 3.2 4.9 4.2 4.3 4.6 5.3 5.4 4.5 5.0 4.5
##  [631] 5.7 4.5 4.7 4.9 3.1 5.3 4.7 4.7 3.4 3.3 5.0 3.6 5.8 4.0 5.8 4.4 3.7 4.0
##  [649] 5.3 5.3 3.6 5.6 4.7 5.1 4.7 5.6 3.4 4.0 4.0 5.0 4.7 5.0 4.3 5.7 3.3 5.7
##  [667] 4.8 3.9 4.1 4.6 4.4 5.2 4.4 2.8 5.7 5.2 3.0 3.8 4.7 4.6 4.5 5.2 6.2 5.2
##  [685] 4.9 5.6 4.5 4.1 7.2 5.5 4.8 2.4 5.3 5.3 3.6 4.1 3.9 4.7 4.7 4.9 4.7 5.7
##  [703] 3.3 4.1 3.5 4.4 3.8 6.0 4.2 4.9 4.9 5.1 4.8 4.5 4.1 6.5 6.8 4.0 4.9 5.8
##  [721] 4.9 5.8 5.4 4.2 3.8 3.5 6.5 4.7 2.6 3.9 4.5 4.7 5.0 5.8 2.9 3.0 3.0 5.9
##  [739] 5.1 5.4 5.3 3.6 6.3 4.2 3.8 5.9 3.6 2.9 4.9 4.1 5.0 3.9 3.8 3.8 5.5 4.9
##  [757] 5.5 5.0 3.2 4.1 5.3 3.1 3.3 4.9 6.3 4.2 2.9 4.6 4.2 4.6 3.9 6.2 4.0 4.5
##  [775] 5.0 5.3 3.6 2.4 4.9 3.3 5.2 4.6 4.8 3.3 4.7 5.1 5.4 4.2 5.5 4.2 4.2 4.0
##  [793] 3.9 2.3 3.9 3.8 3.6 5.3 4.8 4.9 3.0 5.1 3.8 3.7 4.8 6.2 4.3 4.3 5.8 4.5
##  [811] 4.0 4.9 5.1 3.4 6.4 4.2 5.2 3.1 4.4 3.9 2.7 3.9 3.1 5.3 6.0 6.9 3.8 5.2
##  [829] 4.3 4.7 3.6 4.9 4.5 4.3 4.7 4.3 4.4 4.9 5.0 5.7 4.8 5.2 4.9 4.4 4.9 6.0
##  [847] 4.5 6.3 4.2 4.5 5.2 4.8 3.2 5.1 3.0 5.4 4.7 6.2 4.8 4.6 5.1 3.7 3.3 5.0
##  [865] 5.4 6.0 5.5 5.2 5.2 5.0 2.3 5.2 4.5 3.9 5.6 3.5 5.9 5.1 3.3 4.3 5.3 4.2
##  [883] 4.4 4.6 4.2 4.2 4.8 1.7 4.8 4.9 5.0 5.5 4.4 5.0 4.2 3.9 5.5 4.2 3.5 5.7
##  [901] 3.5 5.2 4.4 2.3 5.0 4.5 3.4 5.3 5.3 5.2 3.5 3.7 4.8 4.7 3.2 4.6 3.7 3.9
##  [919] 4.5 4.6 5.3 5.3 2.9 4.7 6.1 5.5 6.3 4.2 3.6 3.5 4.6 3.3 3.5 5.1 3.5 4.0
##  [937] 5.0 5.9 6.3 4.4 5.5 4.4 5.0 3.2 5.0 4.4 3.9 4.7 4.9 4.3 4.3 6.3 2.9 4.8
##  [955] 3.4 4.9 5.0 3.9 4.1 4.3 2.5 4.9 4.7 4.9 4.2 4.6 4.5 4.9 5.1 4.5 5.5 3.7
##  [973] 3.5 4.4 3.9 5.3 3.9 4.8 5.1 3.8 2.6 5.5 5.8 5.2 4.3 4.6 4.9 5.2 4.9 3.5
##  [991] 2.7 5.3 4.9 4.8 4.6 3.7 5.8 5.6 4.3 4.4
## 
## $func.thetastar
## [1] 0.0304
## 
## $jack.boot.val
##  [1]  0.56897507  0.41436620  0.36299435  0.18641618  0.08966667 -0.09122340
##  [7] -0.12831858 -0.21352785 -0.40430769 -0.46959064
## 
## $jack.boot.se
## [1] 1.001478
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
##    [1] 4.3 4.0 3.8 4.3 3.8 4.5 3.1 5.6 2.6 4.8 4.3 6.3 5.2 3.7 4.6 4.4 3.7 4.5
##   [19] 2.9 4.5 4.8 3.6 4.0 3.6 5.3 4.4 5.7 3.7 3.3 5.7 4.4 4.3 3.4 3.9 4.3 4.2
##   [37] 5.2 4.9 3.0 6.3 4.8 5.1 4.8 3.3 3.3 5.3 6.0 4.4 3.1 6.8 4.4 4.8 3.8 2.6
##   [55] 5.7 3.5 4.8 6.3 3.0 4.1 4.7 3.8 6.4 5.8 2.7 5.8 4.3 4.8 3.6 3.3 5.6 3.7
##   [73] 5.2 6.2 4.3 3.7 5.0 3.8 4.2 6.3 4.9 3.9 5.7 2.5 4.7 4.8 3.5 4.3 3.0 4.2
##   [91] 4.7 3.1 4.2 5.7 4.7 4.7 6.3 5.1 2.9 3.9 4.4 4.4 4.9 4.0 3.6 3.6 4.1 5.4
##  [109] 5.3 4.0 5.2 4.2 4.0 4.1 4.3 5.7 3.9 5.6 4.6 3.2 4.5 5.3 4.1 2.6 6.1 3.3
##  [127] 6.2 4.4 5.1 4.5 6.0 4.7 3.1 4.6 4.9 4.4 4.5 5.1 2.9 5.4 4.5 4.8 4.4 4.5
##  [145] 6.2 4.6 3.8 3.6 4.4 5.3 5.2 4.8 3.1 3.8 5.3 5.6 4.0 5.2 2.5 3.9 4.3 5.1
##  [163] 5.1 3.7 4.4 4.2 4.1 3.9 5.2 4.8 3.8 5.6 3.9 6.5 2.5 3.9 3.8 3.9 4.5 3.8
##  [181] 3.4 3.4 4.4 2.9 5.8 4.4 5.4 4.4 4.8 4.1 3.7 4.8 5.4 4.9 4.9 4.9 5.1 4.1
##  [199] 5.2 4.5 4.0 4.9 3.0 3.9 3.9 5.2 4.0 5.0 4.7 3.1 5.8 5.1 3.2 4.3 3.4 4.0
##  [217] 4.2 4.6 4.2 4.2 3.5 4.5 4.1 4.6 5.9 5.7 5.7 6.3 5.1 3.9 4.6 3.8 4.0 3.9
##  [235] 5.0 3.7 4.5 5.2 5.0 3.7 6.4 3.8 5.8 4.3 5.9 6.6 5.3 5.2 4.0 3.9 5.0 5.1
##  [253] 4.5 4.9 6.0 4.6 5.4 4.5 4.3 4.5 3.9 4.2 5.2 4.7 5.8 3.7 4.5 5.4 5.4 5.6
##  [271] 4.9 4.9 4.0 6.4 6.3 5.0 3.4 4.2 4.0 5.8 5.5 5.1 4.8 5.8 4.4 3.8 5.2 3.9
##  [289] 5.7 3.7 3.9 3.5 4.8 2.7 3.8 3.8 5.6 4.1 5.3 2.9 3.0 3.6 3.0 4.8 3.7 3.7
##  [307] 2.9 5.7 3.6 4.2 3.5 4.1 4.4 5.1 3.5 4.3 3.8 3.5 4.0 4.8 5.1 6.0 5.0 3.8
##  [325] 3.9 4.3 4.7 4.7 3.3 4.2 3.4 4.5 7.0 2.3 4.3 5.5 3.5 5.0 1.5 5.3 3.9 4.3
##  [343] 3.9 3.8 3.8 1.5 5.5 4.8 4.8 4.4 4.2 4.5 6.2 5.2 5.1 3.3 5.5 3.3 5.6 4.1
##  [361] 4.7 3.2 4.2 4.4 4.5 5.1 4.2 4.5 6.1 4.7 3.9 4.5 4.3 3.4 4.3 4.2 4.5 5.2
##  [379] 5.9 5.4 3.2 5.4 5.9 5.4 4.8 4.5 4.8 4.5 4.8 5.9 3.2 4.7 5.8 5.3 4.3 3.0
##  [397] 4.4 4.6 4.5 3.9 4.3 5.7 2.8 3.3 4.8 5.0 4.3 5.3 3.2 5.3 5.0 5.6 5.2 3.1
##  [415] 4.6 3.6 5.1 3.3 4.2 4.7 4.4 3.2 4.6 4.5 5.6 4.9 4.8 3.9 5.2 5.5 3.3 3.1
##  [433] 3.3 4.1 3.6 5.1 4.5 4.5 5.4 4.9 4.3 5.8 4.7 5.9 3.6 4.6 4.1 4.8 4.3 2.7
##  [451] 3.4 5.3 3.1 4.1 5.1 4.7 2.9 3.5 5.8 4.0 6.0 5.4 5.3 3.7 4.1 5.2 3.9 5.6
##  [469] 5.6 3.2 2.9 3.9 3.9 4.9 4.9 4.9 4.5 4.8 5.8 4.2 4.0 3.9 2.9 3.7 5.2 5.4
##  [487] 3.0 5.6 4.6 4.1 3.6 4.9 5.0 4.4 4.6 3.8 4.2 5.1 5.1 4.1 4.9 3.7 4.3 4.0
##  [505] 4.0 3.6 4.9 3.7 3.7 5.2 4.1 5.3 5.3 5.8 4.2 3.2 3.5 4.4 4.6 5.0 2.6 3.7
##  [523] 4.3 3.7 3.5 3.6 4.6 4.7 5.0 4.4 3.8 5.5 5.0 3.5 5.8 4.5 3.6 5.3 3.8 4.3
##  [541] 4.7 2.7 6.6 4.6 3.4 3.9 5.4 5.1 5.0 2.9 4.6 3.1 4.9 4.2 3.1 3.3 3.9 3.2
##  [559] 5.0 4.6 4.1 4.1 4.0 4.9 4.8 6.0 4.5 5.2 4.0 4.7 4.5 3.5 5.4 4.5 4.7 5.3
##  [577] 4.4 4.9 4.9 3.4 4.0 3.1 4.1 4.7 3.1 6.4 6.5 4.6 3.6 4.6 4.1 5.1 3.7 3.8
##  [595] 4.8 4.2 6.3 4.7 3.3 4.8 5.1 5.9 4.7 4.4 5.2 4.8 3.6 5.8 4.2 5.6 4.1 5.0
##  [613] 5.1 3.9 3.8 3.8 4.4 3.9 4.4 4.3 4.3 4.4 5.0 3.3 4.0 5.8 4.4 4.7 4.7 5.2
##  [631] 3.4 4.0 3.4 3.1 4.0 4.1 4.1 4.0 6.8 4.2 5.8 4.5 2.8 3.6 4.9 5.4 5.4 3.4
##  [649] 3.1 3.6 5.3 5.6 3.1 4.4 6.1 4.3 4.2 4.9 5.1 4.9 3.9 3.2 5.2 3.2 4.8 4.2
##  [667] 6.5 3.8 4.6 4.4 5.3 5.1 3.6 4.6 5.9 3.6 3.0 4.9 4.7 4.5 3.2 4.9 5.3 5.0
##  [685] 6.2 4.0 2.6 5.1 4.7 5.0 4.6 3.2 2.9 5.4 6.2 2.5 2.7 6.0 4.9 3.5 3.9 4.9
##  [703] 4.2 5.0 4.3 3.2 4.7 3.3 2.9 4.5 4.0 3.8 4.0 5.8 3.2 3.6 3.8 4.6 5.3 4.8
##  [721] 4.0 4.4 3.9 6.2 3.7 7.3 4.9 4.1 5.2 4.5 4.5 3.6 3.5 5.4 4.9 3.2 5.2 4.6
##  [739] 3.5 4.8 4.8 4.8 3.2 3.9 3.8 4.1 3.2 4.0 4.2 7.0 5.1 3.0 5.6 6.0 4.7 4.3
##  [757] 5.0 5.3 4.6 4.3 3.6 3.7 5.2 4.2 3.6 4.3 4.1 4.6 3.4 3.6 5.1 4.4 3.3 4.6
##  [775] 5.3 3.1 5.9 3.4 5.5 4.4 5.5 4.5 5.2 3.1 4.5 4.1 5.2 3.5 3.6 4.8 4.1 3.5
##  [793] 3.7 3.7 3.9 3.3 4.9 4.7 5.7 5.6 3.6 4.8 5.8 4.4 3.4 4.9 4.5 5.5 5.0 4.5
##  [811] 6.2 4.0 3.9 4.7 5.2 6.3 4.2 5.3 4.5 5.3 4.1 2.5 5.7 4.3 4.5 5.0 4.8 3.9
##  [829] 2.4 4.4 4.9 4.7 5.4 4.9 5.2 4.9 4.4 2.8 6.2 6.4 4.3 5.6 4.6 5.1 4.9 5.4
##  [847] 5.9 3.5 5.3 3.4 4.3 4.0 4.9 5.6 6.0 5.5 4.9 4.1 4.6 3.6 3.9 4.1 4.6 4.6
##  [865] 3.9 5.1 2.5 5.4 4.4 4.0 4.2 3.1 3.4 4.3 3.2 3.1 3.7 5.1 2.8 4.2 4.9 2.9
##  [883] 4.7 4.6 5.4 2.7 7.1 3.5 3.7 4.0 5.1 5.2 3.8 3.8 4.5 5.0 5.1 4.6 6.0 3.0
##  [901] 5.9 4.8 3.6 4.5 5.0 3.8 4.8 4.9 6.0 3.0 4.1 4.8 6.2 3.9 5.0 5.0 5.5 4.8
##  [919] 3.7 4.6 3.2 5.1 3.4 3.8 3.8 4.4 5.1 4.9 5.3 5.4 6.4 5.2 4.7 4.4 3.4 6.7
##  [937] 3.8 4.1 4.3 3.8 2.9 5.0 5.7 5.1 4.4 4.4 4.3 3.9 6.0 4.7 3.8 3.6 3.1 3.8
##  [955] 5.4 3.7 4.9 4.5 6.1 3.5 3.2 4.8 4.7 4.2 4.5 6.1 4.2 5.9 5.2 5.1 4.8 5.2
##  [973] 4.2 3.7 5.3 4.9 2.9 5.0 4.5 4.5 5.4 4.8 4.5 6.8 3.6 4.1 5.4 5.5 4.5 3.3
##  [991] 5.9 4.5 5.0 2.7 3.1 4.7 3.3 5.4 2.9 6.1
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.4 5.4 5.3 5.2 5.1 4.9 4.8 4.8 4.6 4.4
## 
## $jack.boot.se
## [1] 0.980867
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
## [1] 0.5503048
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
##   5.007830   8.944695 
##  (2.169085) (4.075266)
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
## [1] 0.2784519 0.2607583 0.8423341 0.3392877 0.3447115 1.5299831
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
##    [1]  1.1991498655  1.5153053150  0.1738500179  0.3670350153  0.2345953598
##    [6] -0.1319104439  0.3364025865  0.5540267916  0.6747722248  0.6619552334
##   [11]  0.5319863778 -0.3714034195  0.5245830849  0.8205992616  0.0065214613
##   [16]  0.5661355421  0.6479624558 -0.2961551053  0.8242407486  0.1388605843
##   [21]  0.2194822955  0.6040211000  0.2364525597 -0.2568987078  1.5115311283
##   [26]  0.8591128847  0.4730402167 -0.0554714424 -0.1800244684  0.7655480644
##   [31] -0.3790734595  0.1912263421  0.9886493189  0.1516725133  0.3696628442
##   [36] -0.3283922601  0.1481614202 -0.1469427219 -0.2313147661  0.5490691953
##   [41]  0.2639406724  0.0301617162  1.6356474051  0.2745350807  0.1553035426
##   [46]  1.0890453848  0.3442012644  0.4635554898  0.9558365799  0.5102359443
##   [51]  0.7247485520  0.4765456765 -0.1049696885  0.6459915273  0.0949412452
##   [56]  0.5686599628  0.2915420973 -0.0623104401 -0.0978253045  1.6510021775
##   [61] -0.0059512207  0.7859463530  0.3839998832  0.3878277885  0.2050749974
##   [66] -0.2606398498  0.6860206540  0.2410387076  0.4159304340  0.2545878734
##   [71] -0.4818563175  1.5531405103  0.3278334735  0.3761615941  0.0248592997
##   [76]  0.0956291720  0.5757868619 -0.1774623715  0.6635294644 -0.4201288736
##   [81]  0.7161055063 -0.6596957915  0.3812519621  0.0452515060  0.6591701352
##   [86]  0.5867649645 -0.1677885580  0.6161465342  0.7577526286  1.1029465769
##   [91] -0.1695074941  1.5598277475  0.3471025883 -0.1872298026  0.4401886300
##   [96]  0.2428089911  0.8427447695  0.0435166847  0.7613139669  1.4090230348
##  [101]  1.4553178130 -0.3012490772  0.2670742548  0.8729636055  0.1707600095
##  [106]  0.9658726284  0.4146067802  0.6623759939  0.1256995942  0.5785955429
##  [111] -1.0658719316  0.0033958758  0.2160087134  0.6572425048  0.4405636605
##  [116] -0.2010533255 -0.1740675683 -0.5400606485  0.9237053287 -0.1984349627
##  [121] -0.6181333515  0.1891709616  0.7910671237  0.0785298755 -0.5078192882
##  [126]  0.6157931126  0.2282849584  0.9343551214  0.4172924417 -0.1710088089
##  [131] -0.2083948089  0.9231847716 -0.3044090116  1.0380891171  0.7222764939
##  [136]  1.2212147237 -0.3270814482  0.2400758185  0.2220720474 -0.3929569267
##  [141] -0.5018110020  0.3138928188  0.4065612088  0.8253581145  0.3710934844
##  [146]  0.3802144663  0.2014669107  0.4631617658 -0.1495896805  0.0710452923
##  [151]  0.6112568733  0.3666506191  0.4912345072  0.2223567012  1.0364768956
##  [156]  0.6127284036  0.7012393734  0.7411932061  0.6115203949  0.2457144190
##  [161]  1.0547010897  0.2970440645  0.1066311937  0.6750766128  0.1569100712
##  [166] -0.0244197489  0.7716707183  0.4821612711  0.0550313638  0.4469070748
##  [171]  0.6379014235  0.1058586473  0.2350577027  0.1586831991  0.8522449960
##  [176]  0.8031926423  0.3114391514 -0.0942974211  0.8917190609  0.3338260629
##  [181]  0.1149922791  1.2797642139  0.8046183914  0.5368812832  0.1657204638
##  [186]  0.5643677734  0.5781772270 -0.0814695030  0.4960967388  1.1539357283
##  [191]  0.2082443593  0.5526628420  1.2509297811  0.8765918715  1.3495211564
##  [196]  0.3739174917 -0.3788708063  1.1267112484  1.2676234249  0.2565449554
##  [201]  0.0298103301  0.7967327326  0.1652504197  0.0532427730  0.5523537723
##  [206] -0.3146616740  0.0198428977 -0.0013523967  0.0222272476  0.6039562002
##  [211]  0.7876725914  0.2872285082  0.4963494137  0.4855170768  0.1498130153
##  [216]  0.2803484007  0.4932187346  1.1610094822  0.3626599765  0.5768616102
##  [221]  0.1558074970  0.3782099945  1.0178102070  0.5515086938  0.6756802381
##  [226]  0.1589781588  0.2588924040  0.2087866623 -0.0563835068  0.0680334435
##  [231] -0.1651480411  0.2567845393  0.7650215395 -0.0654232917  0.1027215279
##  [236] -0.5421589879 -0.2978983174  1.0604408437  0.2849128068  0.9421996478
##  [241]  0.0715544551 -0.2440300412  0.8176515179  0.7131662164 -0.0013706557
##  [246]  0.7769758621  0.3054992432  0.0899177969  0.8505491990  0.6626513866
##  [251]  0.8184088708  0.2164845710  0.2532101477  0.4541870445 -0.2359938080
##  [256] -0.1503948249  0.7569126196 -0.6207563722  0.4729942235  0.2041235888
##  [261] -0.2423208794 -0.3835177326 -0.9269415742  0.0404355497  0.1372413262
##  [266] -0.6329148689  0.9228171632  1.0551354648  0.5006125500  0.8321233960
##  [271] -0.3866500620  0.6028136094  0.4858948010  1.7307397443  0.5805153321
##  [276] -0.0943918220 -0.6481041902  0.2643747932 -0.3217793573 -0.1268161884
##  [281] -0.3973162695  1.5599605462  0.6724025549  0.9854753935 -0.2197865755
##  [286]  0.0581067703  1.2320190890 -0.0680481358 -0.4172799215 -0.2425394776
##  [291]  0.6728962508 -1.0850183473 -0.9288531110 -0.0786767747 -0.7383675200
##  [296]  1.2333385420 -0.2886739860  0.7440076451  0.4337717436  0.0050437625
##  [301]  0.0264281763  0.1295227844 -1.0790800999  0.7956088370 -0.2074519680
##  [306]  0.3914832978  0.7372066040  0.7319090887  0.5757319485  0.4861772241
##  [311]  1.5376854093  0.5595614525  0.0547522486  0.0719718100 -0.0762810966
##  [316] -0.6140187211  0.1737308264  0.2835141250  0.9104460502  1.0665674751
##  [321] -0.2223650087 -0.3146616740  0.7199805514  0.2178610168  1.3700276129
##  [326]  0.2734000400  0.3727217873  1.0491008442 -0.5858244035  1.6886457721
##  [331]  0.2520419490  0.4600467108  1.0836168497  0.4555901186  0.9107929641
##  [336]  0.5657317541  0.1530346552  0.2857756135  0.0473951490 -0.0371358838
##  [341]  0.2016908576  1.6085062174  0.5635595302  0.0233941809 -0.5737218668
##  [346] -0.0674119000 -0.2342539668 -0.0199346918  0.0697939985  0.3855603048
##  [351]  0.3771842736  0.6901882798  1.2538483586  0.7179134146  0.6459411056
##  [356]  0.5178144174  0.5253699264 -0.1790792567 -0.1337214897 -0.3931393863
##  [361]  0.5035297054  1.8940449986 -0.1457199350  0.8553696757  0.6662228730
##  [366]  0.5719388874  0.6245620583  0.5796118737 -0.2715093510  0.4852373257
##  [371]  0.2927515848  0.5265820749 -0.3117941364  0.5676565542  0.5561307710
##  [376]  0.6154010551  0.2277889734  0.8581769902 -0.2026360919  0.4077624933
##  [381] -0.0820668553  0.0212846681  0.0594760192 -0.5881374309  0.3452311696
##  [386]  0.9265296674  0.8176515179  0.0129959618  0.4217353827  0.5258607478
##  [391] -0.0657969535  0.3996196682 -0.0944648252  0.5896057530  0.1764775772
##  [396]  0.4382729207  0.7353858557  0.1209917697 -0.3036428184  0.7941194090
##  [401] -0.3488710517  0.1619179510  0.9644398708  0.5552044959  0.9462301425
##  [406] -0.0820668553  0.5560316789  0.7215093396  0.4398488361  1.0801328934
##  [411] -0.2451581671  0.4711932922 -0.2431086717  1.0101047445  0.7145448287
##  [416]  0.5552356254  0.3540998969  0.6092583915  0.8263779764 -0.7759903978
##  [421] -0.1275400932 -0.3530049619 -0.1455675909 -0.1626227394  0.5457507643
##  [426] -0.1772784719  0.6842146439  0.9862435499 -0.1902471943  0.4677377891
##  [431]  0.1494344341  0.1978926992  0.9228171632  0.9581748719  0.7370177200
##  [436] -0.8660732237 -0.2568987078  0.4569779741  0.6247599697  0.7271215918
##  [441]  0.2231304567  0.8198251439  1.4136504326 -0.5017335564  0.6356133721
##  [446] -0.7596463319  0.7445236851  0.1290570240  0.6655460283  0.7858870692
##  [451]  0.0057374119  0.6086802080 -0.3800144428  0.1123346369  0.3505921340
##  [456]  0.1544606966  0.7391091762  0.7847422920  0.5305718558  0.3138473934
##  [461]  0.3455523558  0.1086714421 -0.7378746004  0.4833065344  0.3314935336
##  [466]  0.5421139585  0.5980491460  0.2829988493  0.3678441210  0.7308741546
##  [471]  0.7139827783  0.8939629621  0.0479605339  0.7166407376  0.9960461222
##  [476] -0.4450965055 -0.1328779829  0.1796847788  0.2889851559  0.8402110703
##  [481]  0.2338004661 -0.0730674592  0.1561745376  0.9084076981 -0.1649604707
##  [486]  0.1685135038  0.5528055366  0.2242446356  0.4677377891  0.1069207597
##  [491]  0.3736578930  0.4631628429  0.0443986898  0.0444292260  0.7641904057
##  [496] -0.1080383070 -0.0388113894  0.5169987955  0.1125617049 -0.2997945690
##  [501]  0.4035413027  0.6696218957 -0.4702044177  0.1825425369  1.0101047445
##  [506]  0.1079310608 -0.0946400956  0.0552118914  0.5580603443 -0.0761844106
##  [511]  0.3585879732  0.8420573157  0.1797331548  1.0730662192 -0.0059384003
##  [516]  1.3918361761 -0.0404184489  1.0417470110  0.4995027553 -0.0978364194
##  [521]  0.9100575951  0.6509513255 -0.7166732627  0.5275761216  0.4112747799
##  [526]  0.6244135779  0.1647962108  0.5037440491 -0.0233409803  1.0859149157
##  [531] -0.1757399725  0.5858693243  0.6216710680  0.2503968782 -0.2243716639
##  [536]  0.0106054000  0.2420993075  0.7634614593 -0.2442484174  0.4793956204
##  [541] -0.4057862346  0.3704957608  0.4383910184 -0.4259628273  0.2298760001
##  [546]  0.5771331543 -0.3396234945  0.5363631341  1.0255547755 -0.4382384415
##  [551]  0.1518299271  0.5028274849 -0.1811640652 -0.0930127555  1.3595860170
##  [556]  0.5165061082  0.6262929380  0.5604494733 -0.4838841297  0.4884284795
##  [561]  0.0299283888  0.9195205625  0.8913244632  0.9004274621  0.5668387862
##  [566]  0.2718843978 -0.2280595741 -0.0488354483  1.0421446489  0.8609552929
##  [571]  0.1062398509  0.1386536400  0.5264684232  0.0889784183 -0.0594648106
##  [576]  0.4245877740  1.3748677847  0.3764727198  0.5574302851  0.1440909171
##  [581]  0.7122688821 -0.2552507249  0.1396370366 -0.0365621681 -0.0187534540
##  [586] -0.0002842684  0.2453815202 -0.5953651249  0.1264389679  0.8552326764
##  [591] -0.3178148590  0.2404097291  0.0106732000  0.4627050481  0.7585214060
##  [596]  0.2783382992 -0.1781227001  0.3796740343 -0.1539158574  0.2448184136
##  [601]  0.7044121663  0.3946689617  0.3376475999  0.7252639193  0.8255208224
##  [606]  0.3165602975  0.3364025865  0.2109886077  0.1129120290  0.1204051684
##  [611]  0.2290041629  0.7753274380  0.5559375032  0.3487602837 -0.3153548383
##  [616] -0.4842397214 -0.3042473431  0.2277889734  0.5144992164 -0.2644106340
##  [621]  0.2561745164  0.1350715374  0.6055364708  0.8440298016  0.4737636863
##  [626]  0.4481795174  1.2602641496 -0.2190711109 -0.1052601880 -0.1300697366
##  [631]  0.2116688319  0.2945576302  0.3133586212  0.2500798710 -0.2024306620
##  [636]  0.6674386241  0.1412029870 -0.9217041857  0.2731005293  0.3283001846
##  [641]  0.3270191583 -0.9502701064  0.1480267277  0.4821092670 -0.1703117911
##  [646]  0.0550451806  0.5050219041 -0.1734791447  0.4161656847 -0.1184833965
##  [651]  0.7659332348  0.7345139029  0.2682068236  0.7590601049 -0.4586376183
##  [656] -0.0051603582  0.3228194683  0.2676251573 -0.3211878206 -0.3653332522
##  [661]  0.5838223209 -0.1258888322 -0.0438436745  0.5335099724 -0.4270424405
##  [666]  0.1616198008 -0.1140943667  0.2283293049  0.0484941014  0.2166122215
##  [671]  1.4934990482  1.0393726625  1.0458982986  1.2003635163  0.3381565713
##  [676]  0.1827251595  0.3079989185 -0.2183668112  1.7469359110 -0.9673132922
##  [681]  1.0072744997  0.7238773137 -0.5752606951  0.6172349914  0.7751357967
##  [686]  0.3883836805  0.2743380651  0.4792736736  0.8515935598  0.6020895629
##  [691]  0.2407554167  0.4236315136  0.5495511545  1.0040759707 -0.4685843025
##  [696]  0.0292903957 -0.7438372928  0.5676565542  0.4009685669  1.9110277294
##  [701]  0.8273500481  0.9333877984 -0.1053894642 -0.0168693712  1.8229729594
##  [706]  0.4054153458  0.2922595876  0.3327074497  0.6004815376  0.5323812614
##  [711]  0.2547402954  0.0372625875  0.2481752501 -0.0395551300  1.1863747374
##  [716]  0.6967670391 -0.0157035450  0.2477768281  0.4111320318  0.6892773915
##  [721]  0.4205975099  0.4372963316 -0.2579827897  0.5417790208  0.0502652627
##  [726]  0.0829983650  0.2022823713  0.8819424507 -0.2514211147  0.4670515395
##  [731]  0.3278334735  0.0674812790  0.8682888402  1.3459223557  0.9094195098
##  [736] -0.1777312278  0.5452536308  0.6091722261  1.0970285662  1.1138469461
##  [741] -0.5576891696 -0.2814039641  0.4337717436 -0.0958685362 -1.3718849186
##  [746]  0.2055996437  0.1865339739  1.0614023655  0.5099282965 -0.1283037797
##  [751]  0.6725546534 -0.3063481512 -0.1283651357  1.1878806521  0.2292275820
##  [756]  0.8108914409  0.9069937250 -0.0635620002  0.5072296149 -0.0274873529
##  [761]  0.3368285766  0.1375997200  1.0266772574 -1.0354448127 -0.2446376018
##  [766]  0.3123629890  0.3309622024 -0.0923747444 -0.4705856237 -0.6982491054
##  [771] -0.0956080639  0.8295439351  0.0779739979  1.0209406699  1.0266164423
##  [776]  0.9258041295  0.5405097692  0.2109472786  0.7556008851  1.3234191857
##  [781]  0.5949751381  0.7406895221  0.8311342677  0.7293255004 -0.2194205281
##  [786]  0.6013415923  0.7019401437  0.7374418826  0.3756929583  0.8534134919
##  [791]  0.1020212620  0.3919568081  0.9351614824  0.2964454258 -0.0958788756
##  [796]  0.9204276321  0.0020743190  0.8763239473  0.6055408966  0.6595995142
##  [801]  0.5965061400  0.8684517191  0.4789910453  0.0034739981 -0.1480004111
##  [806]  0.5137229077  0.4960174685 -0.0344023794  0.2613677426  1.0479626467
##  [811]  1.0378295354  0.9626319998  0.9989740966  0.2882536463  0.8729607268
##  [816]  0.3744299670  1.1727648188  0.7209135911 -0.3087441754  1.1403840865
##  [821]  0.5642090566  0.8060010649  0.9613950356 -0.4850629843 -0.3781966406
##  [826] -1.0408903467  0.6777771685  0.9639801429 -0.3475468728  0.8386130147
##  [831]  0.6610737102  1.0160782690 -0.2019114593  1.1294004714  0.6777771685
##  [836]  0.2735300861  0.1531969067 -0.0099921971  0.1037603661  0.8850033827
##  [841]  0.4638542987  0.0120972426 -0.2034888461  0.4053997843  0.8362833483
##  [846] -0.3990137527  1.1699120338 -0.5144526804  0.8399258260  0.4063294336
##  [851] -0.5938766014 -0.3968141303  0.9971193455 -0.7643179395 -0.8740825584
##  [856]  0.6248439021  0.1947276705  1.3833841940 -0.4546742435  0.0888303222
##  [861]  0.0605985130  0.0055779286  0.3326229974 -0.7563711984  0.1633127443
##  [866]  0.5082678242  0.3729467824  0.3168778901  0.0917881716 -0.8180568067
##  [871] -1.7797311513  0.7986721947  0.7185005544  0.6864766228 -0.3316594168
##  [876]  0.3562895777  0.4956796818 -0.1405954277  0.2960632280 -0.0946059666
##  [881]  0.4394273730  0.4339972439 -0.3474293416  0.3568110003  0.2431054512
##  [886]  1.6154466679  0.4974770943 -0.3412913117  0.2157083838  0.5975518493
##  [891]  1.5032209902  0.5849303773  0.1609106429 -0.3768661207  0.4923864499
##  [896]  0.7089002968  0.1425679417  0.1550407324  1.0646474809  0.2351275259
##  [901]  0.3202655378 -0.5434569875  0.7506063403  0.0333842514  0.0101821868
##  [906]  0.4617154408  0.4481401962  0.5290953401 -0.5285751871  1.2816897480
##  [911]  0.6698770427  0.3678441210  0.2558171933  0.0404298741  0.4077058245
##  [916] -0.8966423290 -0.0581020722  0.4468786987  0.2109032430  1.1443987468
##  [921]  0.3409813733  0.4307712220  0.9691377290  0.5258247179 -0.4829900511
##  [926] -0.0305969965  1.3619623591  0.3802591748  0.1693193556 -0.5307908529
##  [931]  0.0577217730 -0.1243212393  1.0415124220 -0.1662114403  0.4580723770
##  [936] -0.0846526062  0.2714802838  1.0009261573 -0.0920901718  0.6402422211
##  [941]  0.4059745323 -0.0476589252  0.7850824648 -0.0536905323 -0.0418565397
##  [946]  0.7672108923 -0.4162633748  0.3361803730  0.3902624587  0.9666377790
##  [951]  0.5588754699 -0.3526745036  0.8447521907  0.7954364644  0.9074468089
##  [956]  0.0900401974  0.8510939663  1.0666892918 -0.0572837767  0.2589252754
##  [961]  0.9726536610  0.5943310692  0.5577045407 -0.1123001187 -0.1335416643
##  [966] -0.1143799771 -0.1841034719  0.7065932153 -0.6064823176  0.3136179048
##  [971]  0.5403535638  0.1334801321  0.7774783905 -0.4003712329 -0.0139936121
##  [976]  0.9323154293  0.3397390847  0.2016384066  0.1372413262  0.7118721405
##  [981] -0.3142631406 -0.9465682245  0.0483688084  1.2680437569 -0.1741738624
##  [986] -0.1576409248  0.0665777634 -0.8382370361  0.0709040807  0.2439412795
##  [991]  1.2600016447  0.5987967653  0.7494930145 -0.4957098198  0.3738490693
##  [996]  0.3746323302  0.9401759189  0.5166144847 -0.0671139775  0.5541463740
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
```

```r
fit2
```

```
##       mean          sd    
##   0.55986341   0.24738116 
##  (0.07822879) (0.05531198)
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
## [1]  0.6311468  0.1405383  0.1189826 -1.0197673 -0.6812570  1.0305982
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
## [1] 0.0105
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8608373
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
## t1*      4.5 -0.03603604   0.8939805
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 2 3 4 5 7 8 
## 2 1 1 2 2 1 1
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
## [1] -0.0278
```

```r
se.boot
```

```
## [1] 0.8933001
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

