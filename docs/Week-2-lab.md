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
## 2 4 5 6 7 8 9 
## 1 2 1 1 1 2 2
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
## [1] 0.0445
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
## [1] 2.688258
```

```r
UL.boot
```

```
## [1] 6.400742
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##   2.5%  97.5% 
## 2.7975 6.5000
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
##    [1] 5.7 5.2 3.5 5.7 3.3 5.1 4.7 4.9 3.9 4.1 4.7 6.6 4.2 3.1 4.5 5.3 4.9 4.5
##   [19] 3.8 3.9 5.3 3.1 3.3 4.9 5.0 5.9 4.1 4.6 5.9 3.7 5.0 4.9 5.0 2.5 5.5 5.6
##   [37] 3.9 4.1 5.9 4.8 3.9 4.3 5.0 3.7 4.9 4.5 4.1 4.9 4.7 3.6 3.5 3.5 5.5 3.4
##   [55] 5.1 6.1 6.5 3.7 4.9 2.9 4.9 4.0 4.6 5.5 4.2 4.9 4.1 4.3 4.3 3.3 4.5 4.0
##   [73] 5.0 3.3 3.5 3.8 3.7 3.9 5.1 4.0 4.4 4.8 2.1 4.0 4.5 4.6 5.1 4.4 5.6 5.1
##   [91] 3.8 4.5 3.7 3.5 5.3 4.2 4.5 6.3 4.0 4.1 3.9 4.7 3.0 3.8 3.7 5.1 6.4 2.9
##  [109] 5.9 4.4 3.5 6.3 5.7 4.9 3.3 4.6 4.0 3.6 4.0 4.6 4.5 4.1 4.8 4.4 4.5 2.8
##  [127] 6.2 3.1 5.4 4.2 4.8 4.4 3.8 4.2 4.2 4.8 4.6 4.1 5.4 4.9 4.0 3.8 6.0 5.8
##  [145] 1.9 4.0 4.7 3.4 5.7 3.6 4.1 6.0 2.8 5.4 5.1 5.9 6.7 4.2 4.5 4.4 4.7 4.9
##  [163] 3.4 2.8 3.4 4.1 4.7 3.8 3.7 6.6 5.0 4.8 6.1 4.5 2.3 3.1 4.9 4.7 3.6 3.6
##  [181] 3.4 3.9 3.8 4.2 5.5 4.6 5.5 6.9 4.3 5.1 6.0 5.0 4.3 3.2 3.5 3.8 3.7 2.0
##  [199] 4.4 4.6 4.6 3.1 4.9 4.1 4.2 3.8 3.8 5.1 4.9 4.6 3.7 5.1 3.5 4.9 5.1 5.1
##  [217] 5.3 4.8 3.7 4.0 6.2 5.2 5.7 2.8 4.8 5.6 5.2 3.9 4.0 4.7 2.7 5.2 4.1 4.4
##  [235] 4.8 4.6 5.1 4.8 5.3 3.1 4.1 3.0 4.1 3.2 4.2 4.1 3.2 4.2 5.6 3.1 3.9 3.1
##  [253] 5.1 5.5 4.0 3.6 3.4 3.8 5.4 4.1 4.4 5.7 3.4 3.3 3.2 4.7 4.3 4.3 4.9 5.2
##  [271] 6.3 4.0 3.2 4.4 4.0 3.2 5.0 4.0 4.6 4.3 4.6 5.7 3.7 4.0 4.2 5.9 6.1 4.8
##  [289] 5.1 4.3 4.7 4.9 4.6 4.5 4.4 3.5 4.5 5.3 3.6 4.5 3.8 6.7 5.3 3.4 5.1 3.5
##  [307] 4.5 3.8 4.3 4.4 3.8 4.1 4.9 3.4 4.9 5.6 3.2 5.4 4.8 4.6 5.8 4.3 3.3 6.1
##  [325] 4.1 5.4 3.2 4.9 4.8 4.8 5.6 3.7 5.3 4.1 4.8 4.7 4.6 4.9 4.8 2.4 5.4 5.0
##  [343] 6.3 3.6 4.6 5.1 5.0 4.5 4.7 4.3 4.2 4.5 3.7 3.1 4.7 4.6 2.3 4.6 4.1 4.3
##  [361] 5.4 2.6 4.6 4.9 5.3 2.0 5.9 2.6 4.9 5.5 4.3 4.4 4.0 5.8 4.4 4.6 4.8 4.2
##  [379] 4.8 4.8 4.6 6.1 5.1 4.0 3.7 5.4 5.3 4.7 3.6 5.5 5.0 3.6 4.1 4.8 3.7 5.1
##  [397] 3.7 4.6 3.8 4.2 5.5 4.6 5.5 3.9 4.5 5.6 5.6 3.8 4.6 4.6 5.3 4.4 3.4 2.9
##  [415] 5.0 4.0 4.7 4.7 5.7 3.9 4.6 3.4 5.4 3.3 3.8 3.0 4.2 3.7 5.4 4.6 4.7 5.1
##  [433] 4.9 5.4 4.4 4.6 3.6 4.3 4.3 4.7 4.0 3.5 3.8 4.2 5.9 3.7 4.3 4.8 5.5 4.7
##  [451] 4.5 4.4 5.6 4.5 5.2 4.3 4.0 6.6 3.6 4.1 5.1 4.2 4.3 2.7 3.8 5.3 4.0 3.2
##  [469] 5.2 3.3 3.8 5.1 5.1 4.2 5.0 4.2 4.5 4.5 4.4 5.0 4.3 5.0 4.1 2.5 4.0 3.8
##  [487] 4.8 4.0 5.1 2.9 6.6 4.0 5.5 4.1 3.3 3.4 3.7 4.6 2.7 3.4 3.5 5.1 3.3 3.6
##  [505] 5.0 3.6 3.4 4.6 4.6 4.7 5.2 6.1 5.1 4.4 4.9 3.9 3.1 4.0 6.0 4.6 4.8 4.7
##  [523] 4.5 3.8 4.9 4.5 2.8 2.6 5.0 5.4 3.4 4.8 2.8 4.4 5.1 5.6 5.1 4.5 5.0 6.9
##  [541] 4.2 4.9 5.1 3.9 5.0 5.2 3.4 4.6 5.0 4.5 5.5 5.5 5.2 5.4 3.7 5.7 2.7 3.9
##  [559] 4.0 4.8 2.9 4.6 4.9 4.4 4.6 6.2 4.6 3.4 2.8 4.2 4.7 4.9 4.3 3.8 3.1 5.7
##  [577] 4.1 4.3 4.1 3.9 6.6 4.3 4.9 5.6 4.9 4.1 5.0 4.2 3.6 5.1 4.3 3.1 4.6 5.2
##  [595] 4.8 3.9 5.0 4.8 6.3 4.6 4.0 4.7 4.8 3.3 5.0 5.1 5.0 2.7 4.4 3.4 4.8 2.6
##  [613] 3.4 2.8 3.4 5.0 3.7 4.4 4.1 4.2 5.9 3.1 5.5 4.0 4.2 3.7 4.0 4.6 4.7 3.9
##  [631] 4.8 5.2 3.6 5.4 4.3 4.7 3.0 4.7 4.9 5.5 4.4 5.4 3.9 3.5 3.9 4.0 4.4 4.4
##  [649] 4.7 4.2 3.1 4.2 4.0 5.8 5.3 4.9 4.5 4.1 4.7 4.5 2.9 4.2 3.1 4.7 3.9 3.1
##  [667] 4.8 4.5 4.7 5.0 4.0 4.5 5.0 4.3 6.3 4.4 5.7 5.4 2.9 2.9 5.7 2.3 4.0 5.0
##  [685] 4.7 4.4 4.0 6.6 4.1 4.2 4.1 5.7 5.6 4.0 5.1 4.6 5.0 3.6 4.0 5.8 3.7 2.9
##  [703] 3.8 3.3 4.3 4.7 3.4 5.0 2.7 5.8 5.9 5.2 4.3 4.4 5.1 3.5 4.0 4.3 2.9 5.6
##  [721] 5.5 4.2 4.8 3.1 2.0 3.4 4.8 6.0 4.3 4.9 6.4 3.4 4.5 5.2 5.7 4.6 4.8 4.3
##  [739] 4.9 4.5 3.1 4.3 3.6 4.2 5.9 5.3 3.5 5.0 5.2 4.3 5.4 4.3 4.6 4.0 5.0 3.2
##  [757] 3.5 5.5 7.6 3.4 5.7 5.0 4.4 5.2 5.0 4.6 3.9 5.0 5.2 5.4 3.7 4.6 5.4 6.0
##  [775] 6.4 5.1 4.8 5.0 3.9 3.6 4.8 3.1 4.9 3.7 2.9 5.1 4.1 3.1 3.6 3.8 5.3 4.7
##  [793] 3.6 5.5 5.3 2.9 4.8 4.2 2.7 4.3 3.7 3.8 4.7 3.9 3.9 4.1 5.4 2.6 4.0 4.3
##  [811] 4.2 3.7 5.7 3.3 5.9 4.7 3.8 3.5 4.9 2.8 4.8 3.9 3.2 3.6 4.2 5.7 4.6 2.4
##  [829] 3.9 4.0 6.4 6.0 3.7 4.0 4.1 4.5 4.6 5.8 3.0 4.8 4.7 4.6 4.0 4.6 4.9 3.9
##  [847] 6.3 6.5 5.6 5.5 4.3 3.8 5.7 5.1 5.1 5.9 4.1 3.7 2.5 3.5 4.3 2.6 4.6 4.5
##  [865] 5.7 3.8 5.5 4.9 4.6 3.7 3.1 4.5 4.4 3.6 4.3 5.6 3.9 3.4 5.1 3.9 4.6 4.5
##  [883] 4.3 4.3 3.7 3.1 6.3 4.3 3.7 4.4 6.4 4.4 5.8 4.9 3.1 5.0 5.7 3.8 4.1 4.8
##  [901] 4.7 3.8 4.3 5.4 3.1 4.5 5.4 3.5 4.0 4.4 4.2 4.4 3.9 3.9 4.6 3.9 4.9 3.7
##  [919] 4.4 5.8 3.7 5.2 5.2 4.8 5.4 4.4 4.5 5.4 3.9 4.5 4.2 4.4 5.2 3.9 4.1 5.7
##  [937] 3.1 4.4 4.0 4.2 4.2 3.8 3.4 4.3 5.1 6.7 4.8 4.9 4.6 4.4 5.3 2.7 6.1 3.6
##  [955] 3.0 5.5 5.0 4.0 5.5 4.5 4.9 4.7 3.9 4.8 5.5 5.3 4.6 5.0 4.8 3.3 5.4 4.4
##  [973] 4.6 3.3 4.5 2.9 3.3 4.7 4.9 4.1 5.1 4.8 3.9 5.8 4.3 4.8 4.9 4.2 5.0 5.4
##  [991] 3.6 4.0 4.8 5.1 3.3 5.5 4.3 4.9 4.0 4.2
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
##    [1] 6.4 4.8 5.2 5.3 4.0 5.8 4.5 2.4 3.0 4.9 4.2 4.1 3.3 4.1 5.2 4.9 5.5 3.3
##   [19] 6.1 5.4 4.8 5.4 5.7 4.1 5.4 5.0 5.4 4.4 5.3 3.7 4.2 5.9 5.0 5.0 4.3 3.4
##   [37] 3.8 3.7 4.6 3.6 5.7 5.5 4.2 4.2 4.3 4.3 4.9 5.0 4.2 4.2 5.2 3.7 3.8 6.6
##   [55] 5.3 4.0 4.0 4.5 5.4 4.3 4.8 4.0 3.2 3.3 4.4 3.2 6.8 4.4 5.4 5.2 5.1 4.7
##   [73] 4.2 4.0 4.0 2.7 5.6 4.3 4.7 4.0 4.7 2.5 4.5 3.6 4.3 5.6 5.3 5.4 3.1 4.6
##   [91] 4.4 6.6 3.9 3.5 4.6 5.8 4.2 5.6 5.0 2.5 6.7 5.4 4.2 5.2 5.2 6.1 3.7 4.9
##  [109] 3.9 5.2 4.2 2.8 4.0 4.7 3.7 3.7 3.5 4.2 4.6 3.2 5.1 4.6 4.5 5.5 4.2 4.2
##  [127] 4.8 5.2 4.8 2.8 6.1 4.2 5.0 3.9 5.0 4.9 4.1 5.1 4.8 2.5 4.3 3.1 4.6 3.8
##  [145] 4.3 5.7 3.7 2.8 4.5 5.2 3.2 2.9 4.2 4.9 5.4 4.9 3.0 3.8 3.5 5.5 4.5 3.5
##  [163] 3.9 3.8 4.0 3.9 5.1 3.1 3.9 4.8 4.3 5.3 4.6 4.4 3.7 4.3 4.3 4.7 4.0 5.7
##  [181] 4.3 4.2 6.3 6.0 4.5 4.0 3.4 3.5 3.6 4.3 4.6 4.7 3.4 4.8 5.6 4.9 4.1 4.7
##  [199] 6.1 4.2 6.1 4.8 4.8 5.2 5.9 5.6 5.5 7.0 3.0 4.7 3.4 3.8 3.8 3.8 4.1 4.2
##  [217] 5.0 4.7 4.6 4.3 4.4 3.5 4.9 4.6 3.1 4.0 3.9 3.7 5.8 5.1 6.2 4.7 2.5 2.0
##  [235] 5.3 5.0 5.2 4.7 5.2 4.5 5.4 4.4 5.6 3.3 4.3 4.4 3.5 5.3 4.0 6.1 3.9 5.1
##  [253] 5.6 5.3 1.7 5.1 4.0 5.9 3.4 4.7 3.1 4.2 3.5 5.0 4.7 4.7 6.4 5.1 4.7 5.8
##  [271] 2.4 3.6 5.6 3.5 3.9 5.3 4.5 3.2 5.3 5.8 5.1 4.3 4.3 4.7 4.2 5.4 5.7 5.1
##  [289] 4.8 4.1 4.8 4.5 5.2 6.2 5.9 5.1 4.8 3.2 4.0 5.1 5.4 3.5 4.5 3.9 5.2 6.4
##  [307] 4.4 3.4 3.6 5.0 3.0 4.6 4.9 5.6 4.7 3.5 5.0 4.1 6.4 3.9 4.4 4.8 4.1 5.0
##  [325] 6.1 3.8 5.1 5.1 5.2 5.2 5.4 4.0 4.8 2.8 5.5 4.2 6.0 3.6 4.7 2.4 5.2 4.0
##  [343] 5.0 4.4 5.1 5.1 4.5 3.7 4.5 4.3 4.1 5.7 4.5 3.4 3.0 6.3 5.7 6.1 2.8 3.9
##  [361] 4.4 5.1 5.7 5.5 4.2 3.5 4.2 4.5 4.0 4.6 3.3 3.8 4.4 5.9 4.8 3.4 5.3 3.7
##  [379] 3.9 5.5 2.9 4.6 4.4 2.0 4.7 5.8 3.6 3.7 5.7 5.6 4.1 4.4 3.7 3.9 3.3 4.5
##  [397] 4.8 2.8 4.5 3.7 1.9 4.5 3.8 5.1 5.4 3.1 4.4 4.9 4.9 5.4 5.7 5.3 3.5 5.1
##  [415] 3.4 5.4 5.2 4.1 2.9 4.7 4.0 6.4 2.3 4.3 4.0 3.0 6.1 3.2 5.9 3.5 5.7 4.5
##  [433] 5.0 4.7 5.2 4.3 4.5 4.8 4.3 4.3 5.1 3.9 6.0 5.2 4.6 3.8 5.7 4.9 2.2 3.5
##  [451] 5.3 3.5 3.4 5.3 4.0 5.6 5.3 3.4 5.7 5.4 3.6 4.6 3.6 6.4 6.5 5.3 4.3 2.4
##  [469] 4.8 4.8 5.2 5.2 3.5 4.2 3.5 4.5 4.7 5.9 5.6 3.9 4.9 6.9 5.9 6.9 2.6 5.1
##  [487] 4.0 4.2 5.3 4.6 4.1 4.0 4.2 4.4 4.0 3.5 4.8 3.7 5.4 4.6 4.9 4.8 4.0 4.2
##  [505] 3.7 3.2 4.6 4.4 4.5 5.5 4.0 4.7 3.8 4.5 4.3 5.2 3.4 4.7 3.7 5.4 5.1 3.3
##  [523] 6.0 3.4 4.3 5.1 2.7 5.8 4.7 3.7 3.5 4.9 3.9 4.5 3.8 4.5 4.4 4.7 4.1 4.7
##  [541] 4.8 4.4 4.7 3.8 6.1 4.5 5.3 3.4 5.5 4.9 5.1 3.4 2.7 4.2 5.0 3.2 3.4 5.9
##  [559] 5.4 4.0 3.3 3.7 4.3 4.4 4.2 5.5 3.8 5.8 4.6 5.0 4.0 3.3 5.1 4.8 4.7 3.9
##  [577] 4.2 3.4 5.1 5.2 4.9 4.5 6.1 4.5 6.8 2.9 5.3 5.3 4.9 3.9 5.8 4.3 5.3 4.6
##  [595] 3.4 4.9 4.6 5.1 5.7 3.0 3.1 4.9 4.0 3.8 4.1 5.7 5.7 5.2 4.5 4.4 5.0 4.8
##  [613] 4.2 4.7 4.2 4.3 2.8 3.6 4.5 3.9 5.2 5.2 4.2 4.1 5.0 5.1 3.5 4.3 5.6 3.9
##  [631] 3.7 4.8 4.0 4.4 5.3 5.4 5.6 6.8 4.6 4.2 4.0 3.5 4.4 6.7 4.5 4.6 4.0 6.1
##  [649] 3.4 5.8 2.9 4.2 4.7 2.8 5.1 4.1 4.4 5.1 4.5 4.3 4.5 2.8 4.2 3.1 4.7 3.6
##  [667] 4.0 5.4 5.4 4.2 3.9 5.8 3.5 5.9 4.0 3.9 4.1 4.6 4.5 3.1 6.0 3.2 4.6 5.3
##  [685] 4.9 2.9 4.2 4.0 4.9 4.4 3.6 4.2 5.1 4.6 4.9 4.6 4.6 4.5 3.8 5.7 4.6 5.0
##  [703] 4.8 5.2 3.9 4.0 5.0 4.3 3.1 4.5 4.4 4.0 4.2 4.4 4.3 1.8 2.3 3.2 5.5 5.2
##  [721] 3.3 4.8 5.7 4.5 5.5 4.7 6.0 4.7 4.4 5.1 3.9 5.0 4.2 3.7 4.0 4.4 3.5 4.5
##  [739] 4.7 4.2 4.0 3.7 4.2 6.1 5.6 5.7 5.0 5.0 3.5 3.9 6.0 5.7 5.1 4.9 4.9 3.4
##  [757] 3.3 4.0 5.1 5.9 3.0 4.5 5.0 4.3 5.9 3.6 5.3 5.3 5.7 3.8 3.9 6.5 5.0 4.2
##  [775] 3.8 5.9 4.0 4.7 6.9 4.1 3.2 3.8 3.5 4.6 4.8 3.4 6.0 4.2 3.7 4.5 3.1 4.8
##  [793] 4.0 3.9 5.7 6.3 3.8 3.6 3.3 5.3 5.7 5.2 3.3 4.8 3.6 5.1 4.7 5.2 5.2 4.6
##  [811] 4.7 4.9 4.7 3.4 5.5 5.3 4.8 4.5 5.3 4.0 6.0 3.5 2.9 4.6 4.3 3.9 5.8 4.7
##  [829] 4.3 5.1 3.0 4.3 5.3 4.9 4.3 6.0 4.3 5.2 4.0 2.5 3.4 4.9 4.7 4.9 5.0 5.4
##  [847] 4.4 4.3 6.4 5.3 5.7 4.7 5.0 5.2 4.8 3.6 6.1 5.9 3.0 3.6 4.8 3.7 3.3 6.1
##  [865] 3.4 4.0 3.9 3.7 4.7 4.1 5.2 4.6 3.4 4.0 4.8 2.4 3.8 5.4 5.1 5.7 4.4 4.6
##  [883] 2.7 5.0 4.6 3.9 3.3 2.9 3.3 5.6 4.5 4.3 4.7 5.4 2.2 5.6 4.6 3.7 5.0 4.3
##  [901] 4.9 4.6 5.2 4.0 4.7 4.0 4.1 5.4 3.9 4.1 3.7 5.0 4.1 4.4 4.8 4.9 5.7 4.3
##  [919] 3.2 5.0 4.3 4.5 3.8 4.6 4.5 4.5 5.0 5.0 3.5 4.7 3.5 4.2 4.2 5.3 4.6 3.7
##  [937] 3.9 3.7 3.6 2.8 5.1 4.6 4.5 4.5 4.6 6.0 4.3 2.6 4.3 4.8 3.4 3.5 4.8 3.8
##  [955] 4.1 5.6 4.7 5.6 3.8 4.9 5.9 4.2 5.1 3.9 2.5 3.9 6.1 4.1 4.0 4.9 4.5 4.1
##  [973] 4.1 5.2 3.8 4.2 4.2 4.2 5.9 3.0 6.1 2.5 3.1 5.4 4.6 4.6 4.0 6.4 5.2 2.1
##  [991] 3.4 4.7 3.4 4.1 4.5 3.4 4.1 4.8 4.7 3.6
## 
## $func.thetastar
## [1] -0.0081
## 
## $jack.boot.val
##  [1]  0.49859155  0.33916914  0.31626506  0.13658537  0.09316770 -0.04606742
##  [7] -0.16114458 -0.28005952 -0.39482759 -0.53114286
## 
## $jack.boot.se
## [1] 0.9653413
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
##    [1] 3.8 4.0 3.8 4.0 4.8 3.1 5.4 4.8 4.5 5.4 6.8 3.9 1.9 5.1 3.3 3.8 6.8 4.5
##   [19] 4.3 4.1 4.2 4.4 5.9 4.8 5.1 4.8 4.5 5.4 3.8 3.7 6.2 5.6 4.1 4.7 4.3 5.1
##   [37] 5.5 5.2 2.5 3.6 4.1 3.9 3.4 3.7 3.8 4.7 2.7 4.1 3.3 4.0 3.5 5.3 4.2 5.4
##   [55] 5.5 3.6 3.7 6.6 5.4 5.1 3.5 4.5 3.5 3.5 3.2 6.4 5.9 4.8 2.6 5.5 3.5 5.7
##   [73] 3.4 3.5 5.6 3.3 5.0 4.4 3.9 5.0 4.4 5.5 3.5 5.0 4.5 3.5 6.9 6.6 5.1 4.6
##   [91] 4.5 4.7 5.2 6.1 4.0 4.0 5.4 2.9 5.7 5.0 4.4 4.4 4.5 5.2 2.8 4.7 4.3 4.6
##  [109] 4.8 4.4 5.4 2.9 3.4 3.7 7.1 6.3 4.2 3.4 5.1 3.9 2.5 4.6 2.3 4.3 3.3 3.5
##  [127] 4.3 5.6 4.1 4.2 4.9 3.0 6.1 4.6 4.6 3.9 4.4 4.8 4.1 3.3 4.6 2.7 3.5 4.6
##  [145] 4.8 2.2 3.9 5.1 5.6 4.9 4.8 3.2 3.2 4.7 5.9 3.2 6.6 4.1 3.0 3.8 3.8 4.3
##  [163] 4.5 2.6 3.8 4.7 5.7 4.4 5.9 4.4 6.7 4.8 5.0 4.5 4.2 5.9 3.6 3.4 5.0 3.6
##  [181] 3.9 3.5 6.3 5.6 3.6 4.0 4.6 4.9 4.0 4.0 4.3 4.4 3.5 5.5 5.1 4.9 4.2 3.3
##  [199] 3.4 4.2 5.4 5.3 4.2 4.1 4.8 4.8 4.2 5.9 4.7 5.7 4.7 5.6 5.4 5.8 3.7 4.7
##  [217] 4.0 5.3 4.9 2.6 3.7 6.2 4.0 4.9 3.9 4.8 3.9 5.4 4.9 6.6 3.5 4.5 5.4 5.3
##  [235] 5.8 5.1 4.3 2.6 5.7 4.8 5.0 3.9 5.0 5.5 2.7 3.4 3.0 5.1 5.2 3.7 4.5 5.4
##  [253] 3.6 3.5 5.2 3.6 5.3 3.5 3.4 4.0 6.8 3.8 4.5 5.5 4.3 5.0 4.7 5.3 4.9 4.6
##  [271] 5.0 2.6 4.0 4.6 3.0 5.3 2.5 4.1 4.0 6.6 4.9 4.2 4.4 2.9 4.8 4.2 4.3 2.3
##  [289] 5.3 5.2 4.8 5.4 3.5 3.7 4.1 4.0 4.7 4.6 3.9 3.2 4.6 3.6 4.6 4.6 4.0 4.3
##  [307] 4.6 4.9 3.3 3.6 5.6 5.0 4.9 5.1 4.0 4.9 4.6 3.9 3.5 3.7 4.3 2.7 3.5 3.8
##  [325] 4.1 6.5 3.6 3.9 5.4 3.2 5.3 4.9 4.6 5.7 4.7 2.6 5.3 5.2 4.1 5.8 4.0 5.5
##  [343] 4.4 3.2 4.9 3.9 4.6 5.9 3.7 3.7 5.5 3.6 4.9 4.2 5.8 4.9 3.5 4.7 5.1 3.3
##  [361] 5.8 4.6 4.0 3.6 5.3 3.7 3.5 5.5 4.3 5.7 4.7 4.0 6.2 5.1 3.4 5.3 5.2 5.0
##  [379] 2.5 4.8 5.8 3.7 3.5 4.5 5.9 3.5 4.7 5.7 3.3 4.1 4.2 3.6 4.1 5.3 5.4 4.0
##  [397] 4.3 3.5 3.2 4.2 3.7 3.9 2.8 5.5 3.4 3.8 5.2 3.9 3.5 2.8 5.0 2.8 4.7 3.4
##  [415] 6.1 4.4 4.3 3.9 3.6 3.9 4.3 4.9 3.7 4.7 4.3 4.2 4.0 5.0 4.9 3.2 5.7 5.0
##  [433] 6.5 5.5 4.0 3.4 6.4 5.6 5.6 5.3 3.3 4.9 3.2 4.5 2.6 4.7 3.2 3.1 3.8 6.2
##  [451] 5.4 4.6 4.3 6.6 5.3 6.2 6.7 5.1 5.2 6.3 4.2 3.7 4.0 4.7 4.5 4.7 5.2 5.7
##  [469] 5.4 5.2 3.2 4.4 4.3 3.9 3.5 5.1 5.4 3.4 5.7 4.5 4.8 4.0 2.5 4.4 4.3 4.6
##  [487] 3.1 5.7 4.8 3.7 4.4 5.4 5.0 3.0 4.3 3.6 3.0 6.4 5.0 4.9 5.1 4.9 5.0 3.7
##  [505] 4.6 3.8 5.3 4.9 5.0 5.3 4.3 4.0 4.9 4.3 4.8 5.2 4.0 5.7 5.8 3.8 5.1 6.3
##  [523] 3.2 5.6 6.0 2.8 5.1 3.7 3.4 4.6 5.1 5.0 4.8 4.2 3.2 4.9 3.3 4.3 3.7 3.2
##  [541] 4.5 4.3 6.5 4.9 5.9 3.8 3.8 4.3 4.5 2.6 4.4 5.2 5.0 6.9 3.1 5.7 3.9 3.5
##  [559] 5.4 4.0 2.9 6.0 3.4 5.4 6.1 4.4 4.4 5.2 3.8 5.9 4.2 4.5 4.7 5.6 3.4 3.8
##  [577] 3.8 4.2 5.4 4.4 6.0 6.0 5.9 3.6 3.8 5.2 3.9 2.1 3.4 5.9 4.3 3.8 4.1 3.8
##  [595] 4.9 3.6 4.7 4.1 4.4 5.0 5.3 2.8 4.5 4.0 4.1 5.2 4.0 4.8 3.0 3.7 5.5 5.5
##  [613] 5.2 4.5 5.0 3.9 6.0 4.8 5.4 5.7 4.1 4.5 4.9 3.6 5.3 3.0 4.1 3.1 5.7 6.4
##  [631] 4.4 4.4 5.9 6.3 2.5 5.4 5.1 3.4 3.1 4.3 5.3 4.8 6.1 4.6 5.4 5.4 4.2 4.5
##  [649] 4.3 4.4 3.6 5.1 4.1 3.8 5.3 5.5 4.4 4.7 5.4 3.6 5.2 4.0 4.1 6.1 3.7 3.6
##  [667] 3.5 4.8 4.9 6.2 4.9 4.7 4.6 5.1 4.3 3.4 4.8 4.2 3.8 4.1 5.5 3.2 4.4 6.2
##  [685] 6.0 3.7 4.0 2.9 5.1 4.8 5.0 3.7 3.4 3.4 3.3 4.7 4.5 4.8 5.0 4.4 4.2 4.8
##  [703] 4.4 5.5 4.6 3.9 3.7 4.6 4.7 4.6 4.8 5.0 4.0 3.0 4.8 3.6 4.1 5.2 3.3 5.6
##  [721] 4.8 4.4 3.3 5.2 4.0 5.0 3.8 4.5 5.5 4.6 4.9 5.7 6.0 4.8 5.9 4.4 5.7 3.6
##  [739] 4.9 4.3 3.8 3.6 6.2 5.6 4.0 5.0 3.1 4.9 3.4 5.7 4.7 4.0 6.0 4.0 4.7 5.1
##  [757] 5.3 4.4 4.7 5.7 5.1 3.1 4.3 4.0 4.9 3.8 4.7 3.6 5.0 5.4 4.0 4.9 4.9 3.7
##  [775] 6.1 5.3 4.3 4.0 4.1 5.3 4.9 5.3 6.0 4.6 5.4 5.7 4.6 3.6 3.8 2.7 5.3 4.2
##  [793] 3.0 5.3 3.9 5.9 5.0 3.4 4.5 4.4 5.0 3.7 5.2 5.2 4.5 4.5 4.9 5.5 6.1 5.5
##  [811] 4.1 5.2 3.2 5.3 4.2 5.3 4.5 3.2 3.7 4.4 4.7 5.0 3.5 6.0 5.3 4.6 6.4 3.9
##  [829] 4.1 5.1 5.5 4.7 4.2 5.4 5.6 4.7 3.7 4.2 5.7 5.8 5.1 4.5 4.9 3.9 4.1 5.4
##  [847] 4.0 5.4 3.3 7.5 4.4 3.0 3.9 4.5 4.6 5.6 4.7 5.6 4.3 4.9 3.5 5.2 4.4 4.7
##  [865] 5.3 4.7 4.5 2.9 5.3 3.9 3.7 4.7 6.6 4.2 3.0 4.3 5.1 4.5 4.0 4.8 3.9 6.2
##  [883] 8.1 4.8 5.5 4.2 4.6 3.9 5.3 5.0 3.9 3.8 5.4 3.4 4.2 4.9 4.1 4.1 5.1 4.0
##  [901] 4.7 5.4 4.5 5.1 4.7 3.9 4.6 5.5 4.1 4.5 4.6 3.5 4.5 3.9 3.5 3.7 3.7 3.8
##  [919] 3.2 4.4 5.6 5.6 4.1 6.4 3.4 4.3 4.9 5.5 5.6 4.6 4.0 2.9 4.6 3.7 6.3 5.9
##  [937] 4.4 4.6 3.8 3.9 7.0 4.8 5.4 5.7 4.6 4.0 4.5 3.9 5.1 6.3 3.8 5.4 2.7 4.3
##  [955] 4.1 4.9 5.0 4.5 5.4 4.0 3.5 4.0 5.0 3.5 4.4 5.0 3.9 3.8 3.0 4.1 5.6 3.6
##  [973] 3.8 5.3 3.9 4.1 4.2 4.1 4.7 3.8 4.2 4.7 5.2 4.1 6.0 5.5 4.8 2.5 2.7 4.2
##  [991] 3.3 2.7 5.0 5.7 2.3 3.3 3.5 5.4 5.1 5.5
## 
## $func.thetastar
## 72% 
## 5.1 
## 
## $jack.boot.val
##  [1] 5.5 5.5 5.3 5.3 5.1 5.0 5.0 4.8 4.7 4.4
## 
## $jack.boot.se
## [1] 1.014692
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
## [1] 0.7375143
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
##   5.166786   9.472817 
##  (2.240039) (4.313165)
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
## [1]  0.07428751  0.49983812  1.18290411  0.76774277  1.40480808 -0.11342198
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
##    [1]  0.7930379070  0.7225111729  0.3472528506  0.3500761716  0.8331838586
##    [6]  0.3562611234  0.7517815671  0.2044907194  0.8317940235 -0.1341598413
##   [11]  1.3797841838  0.5994498784  0.7479860366  0.3637786252 -0.4001381865
##   [16] -0.0987110435  1.0379332024  1.3450212554  0.3927933191  0.7417626368
##   [21] -0.0779296283  0.8238857303  0.6319098748  0.7502127292  0.7611266832
##   [26]  1.3090722524  2.1291924997  0.2687727149  0.2611719985  1.3591805631
##   [31] -0.8122379396  0.9130782444  0.3923083747  1.3939780012 -0.4897478555
##   [36]  1.2150218093 -0.0117975373  0.7080048844  0.8437652794  1.0306751869
##   [41]  0.3686449751 -0.0221083281  0.7150065629  1.0347840783  0.0235435123
##   [46]  0.7098165275  0.6839244933  0.3996491931  0.8792090110  1.3429156409
##   [51] -0.8464847470  0.0331198405  0.7277436820  1.8022215716 -0.3257225488
##   [56]  1.3732120943  0.2461435207  0.7283938724  0.2584496662  0.0145170161
##   [61] -0.3540006208  0.7467242838 -0.0172977058  1.9139118131  0.8689502411
##   [66]  1.2792477004 -0.0085338485  0.7756665363  0.0154509507  0.2899146731
##   [71]  1.3494576851  0.8439045691  1.1469326376  0.7422788433  1.0806497552
##   [76]  0.7383097256  2.3671848081  0.6392295709  0.7390398535  0.7972498480
##   [81]  0.8621155883  0.8275457578 -0.0617107144  0.7561894628  0.3050185087
##   [86] -0.3967741054  0.7910784760 -0.0001412753  0.2935826892  0.8464796507
##   [91]  1.3303252239  0.7574368294  1.7767736265 -0.0504290393  1.2806180233
##   [96] -0.4196418737  1.3139107705  1.2996738409  0.3377232018  1.2441929517
##  [101]  0.7756343447  1.0057392986  1.2788921554  2.0284956079  0.7302703527
##  [106]  0.6399777129  1.1848247750  0.7260736759  0.8419126272  0.7328919571
##  [111]  0.1048813824  2.2250103025  1.7729791955  0.6961591389  1.8484769484
##  [116]  0.7682101094  1.3409421627  0.7237179573  1.4076444861  2.0859226577
##  [121]  1.9878700777  1.9129010885  0.3258190332 -0.3622465635  1.2810735483
##  [126]  0.3726144071  0.6375399915  0.3858909267  2.1521182762  0.8836599190
##  [131]  2.0394879739  1.0670524429  0.8389797162  1.3082332252  0.6571128946
##  [136]  1.6242849462  0.3826418042  0.0412555375  0.2808843421  2.3111474984
##  [141]  0.8003330160  0.3828273186  0.7702204596  0.7134904241  0.8516752744
##  [146]  1.3566670623  1.1553495838 -0.8386937375  0.7685177295  0.2423769422
##  [151] -0.0089756161  0.5846179570  0.4064825358  1.4489613853 -0.1056320990
##  [156]  0.7096100329  1.3492631855  1.3061053588  1.3952447668  1.1220724394
##  [161]  0.3337939325  1.3511229166  0.7430657084  2.3091984792  2.0776683319
##  [166]  1.4430650484  0.6896345863  2.3030817080  2.1203734285  0.7617544385
##  [171]  0.7225354450  0.7506079122  0.3200746412  1.2366045062  0.6744345649
##  [176]  2.5156005088  0.5748571680  0.8521574947  1.4378752659 -0.3227275777
##  [181]  0.7254513765  2.4212193397 -0.4067091648  2.1692050410  1.2225383029
##  [186]  1.4296933871  1.0700813269  1.2898854111 -0.0850971965  0.7237068597
##  [191]  0.7276272028  0.4131114517 -0.3930976929  0.3847857848  0.7615179291
##  [196]  1.8294764601  1.2702052103  0.3563587119  0.6312383855  0.5740745676
##  [201]  0.2647014612  0.8013116440  2.1037262269  1.6702589967  1.1022780423
##  [206]  1.3675338359  0.3949573911  2.0386882344  0.6966551581  1.2535908418
##  [211]  0.1918747423  0.7181513522  0.2761210089  0.7966181812  1.8254718877
##  [216] -0.0791595178 -0.7925432967  0.0127890938  1.2740398676  1.2062914201
##  [221]  1.8995708397  0.8031405573  1.0907940169  1.9100520124  0.3532457386
##  [226]  2.0553562394  0.4179717787  1.8953096876  0.7581855840  1.3707655848
##  [231]  1.3932175518  0.3090207443  0.3515500978  0.3996213136  0.6190934023
##  [236]  0.7730876551  0.3373951641  0.8190398437  0.2646721605  0.9672418271
##  [241]  1.1353602452  0.8714580135  0.3816463838 -0.3789694849  1.5603550238
##  [246] -0.0864658703  0.7017396090 -0.1078797950  0.7158022345  1.9034987073
##  [251]  0.0076961441  1.2217278886  0.8274636076 -0.4245623633  1.5149087040
##  [256]  0.5686786004  0.2269144219  0.8238031292 -0.3866842709  1.3565098978
##  [261]  0.8205237611  1.3589019608  1.2016615474  1.3775603294  0.3826049746
##  [266]  1.2397313061  2.2854215384  0.7573484492  1.1842478391  0.8332029584
##  [271]  1.2908805216  0.2063479796  0.7090570258  0.6803107368  0.3075902601
##  [276]  0.3879561829  0.7167958598  1.4016570907  0.7270734304  0.8606310683
##  [281]  0.7687755871  1.4055686343  1.2516606919  0.3202094934  1.2147962692
##  [286]  1.2489110325 -0.4262113531  1.3036044016 -0.0822916609  0.6747680925
##  [291]  0.3962311179 -0.0617110576  1.0616118689  0.4294767820  1.8530067593
##  [296] -0.9302217115  1.9449718201  0.7144122547  2.2873281164  1.2258958844
##  [301]  1.0512252578  0.7567447192  1.1787972619  0.0402697353  0.4001251418
##  [306]  0.8086250327  0.4032522053  0.2770688732  0.8425526486  1.2238468893
##  [311]  1.2657046360  1.6170547489  0.6566318908  1.9548582220  1.3207864910
##  [316]  0.3223532366  1.8445726590  0.3533976549  1.1353602452  0.6365355221
##  [321]  1.5644985032 -0.0574886973  0.8012043068  0.7068464188  1.1113834409
##  [326]  1.8316750414  0.6314170912  1.0907940169  1.9626645500  0.7189979117
##  [331]  0.3212441583  0.7517815671 -0.4168577894  2.3111117889  0.6585799881
##  [336]  1.9067817539  1.3535182331  1.2578376647  0.8639510835  0.7938545179
##  [341]  0.3698557095  1.2356230024  0.7958076114  1.5938565833  0.4406128887
##  [346]  0.2446343169 -0.3914261987  0.0486398424  0.7756051198  0.6979223169
##  [351] -0.8669955597  0.0174714598 -0.0943956697  2.0619856640  0.3345369174
##  [356]  1.8166606838  1.6078186642  2.0884877286  0.3461700668  0.2834867539
##  [361]  0.4374856883  0.4278862380  1.3329453296  2.0470801745  1.2862461656
##  [366] -0.0298692729  2.1447732523  1.2397469388 -0.1518365228 -0.0575977520
##  [371]  0.8055368264  0.7646433153  1.1864651893  1.2788385598  1.3905369693
##  [376]  0.4233335350  1.1634605404  0.7050491347  0.2964037108  0.3347403341
##  [381] -0.1232001773  0.3480744008  0.2683062745  0.6392891166  1.3516376678
##  [386]  0.2751722988  0.3342257870  0.4499477165  0.3341534595  1.0157824326
##  [391]  0.9852113012  0.7858793191  0.8398812181  1.3497783245  0.5816976325
##  [396]  0.7220089035  0.7155401174  1.0644928044  0.1989577227  1.2259813065
##  [401]  0.2108229434  0.3223607309  0.6883013642  2.2233019250 -0.0502298312
##  [406] -0.0002173218 -0.4357165248  1.8747688550  0.6798550539  1.2632068601
##  [411]  0.3448726009 -0.4474868003  2.0389877185  0.2993864718  0.6717220598
##  [416]  0.3916184322  1.8634557997 -0.3969467226  1.1162645160  0.7346518352
##  [421]  1.2789196300  0.0142563524  0.7664368457  0.8246411292  1.2844770765
##  [426]  1.2425459299  2.0459031759  0.7881684872  1.3762311305  0.8801745563
##  [431]  1.1672626402  1.2423471294  0.7966910300 -0.3597590136 -0.5091335871
##  [436]  1.3494638970  2.2974291734  0.8657771212  0.7826781173  0.7891612793
##  [441] -0.0305448771  1.2213960856  1.0199295729  0.0228462417 -0.3499998851
##  [446]  0.3819389770 -0.5149502832  0.0173180989  2.1796030859  1.3008591592
##  [451]  0.7693535599  2.2672721983  0.4592276023  1.1243016673  0.2622190002
##  [456]  1.6925783850  1.3832815072  0.7849314914  0.2311973078  0.7299922676
##  [461]  1.8530067593  1.2632068601  0.6409585195  0.7090570258  0.9760350986
##  [466]  1.2397375841  0.6981756092  0.3840887807  0.0399437797  1.3916543967
##  [471] -0.0408649150  1.2129457317  1.3073534474  1.2207927284  0.0450272780
##  [476] -0.4719207541 -0.5667109265  0.8238857303  1.2872743227  0.2230297599
##  [481]  0.8595400132  1.9551548207  1.8061869969  1.1330697168  0.6773940206
##  [486]  1.1407613128  0.6613381421  0.7620602511  0.7730143016 -0.0507396304
##  [491]  0.3917021879 -0.0824348347  1.1204805973  0.0237833490  0.4457030544
##  [496]  0.3611105720  0.3985223369  0.7356532803  1.0761689901  1.1716523488
##  [501]  0.7765256897  1.3117735538 -0.0837104516  0.3314940284  1.8543878303
##  [506]  1.9478951026  0.8849095714  0.3075198277  1.2212338765  0.7389574953
##  [511]  0.3930350399  0.7207301543  1.1551004351  0.3507291471  0.4338184675
##  [516]  1.2572204571  0.0346718813  2.3023856919  0.8157562019  1.2778577259
##  [521]  0.8760447946 -0.4972194813  0.8418038870  1.2810034534  1.8939587098
##  [526]  0.3901594071  1.1392075248  1.3573944681  0.2802187278  0.2947838908
##  [531]  0.0207803720 -0.0712819495  0.6653630884  0.8040266807  0.7801005604
##  [536]  0.3699254696  1.3017704533  0.7680401932  1.3938628044  0.2981509929
##  [541]  0.7847910420  0.2255992202  0.7719765477  1.0750201356  0.7391966839
##  [546]  0.6435905839  0.8255784023  1.1029523841  0.2607574516  1.8320522807
##  [551]  1.4244620063  0.8497615402 -0.0253987896  1.2745954589  0.6729511891
##  [556]  1.5415359448  1.4132369541  1.0829396555  0.2186037326  0.3902634542
##  [561]  0.0217403410  0.7732276195  0.7553871362  0.9109236688  0.4110910264
##  [566]  0.7462723633  0.6827480233  0.3169359318 -0.2029640988  0.7862181717
##  [571]  0.6719599506  1.2019077317  0.6233962278  0.2914384373  1.2149189777
##  [576]  0.2807822459  0.2696120258  1.8429808297  0.7779698920  0.4062086701
##  [581] -0.9003411062  0.7001917334  2.0611868120  0.1901126123  0.0397502779
##  [586]  0.7583119306  0.3936732003  0.6255978549 -0.7757211731  0.3334303488
##  [591]  0.6760469322  2.0014579397  1.4101536844  1.6245386621  0.3328377933
##  [596]  0.6699209799  0.7069069339  1.2332274950  0.3348796160  0.8285052180
##  [601]  1.7625026090  1.2678205169  0.7787884049  0.7342832693  0.3087262394
##  [606]  0.3790590929  2.2192707374  1.2067341522  0.2911158798  0.8085683442
##  [611]  0.7107905995  0.2985650645  0.1636197239  1.1694616261 -0.0505795865
##  [616]  0.3183569301  0.2444019330  1.0942825650  0.2449078624  1.8813821188
##  [621]  1.7943897204 -0.0207774587  0.2981893976  1.1281675352  0.2703583095
##  [626]  0.7727799824 -0.0300245832  0.7318899362  0.3880725439  0.1374750810
##  [631]  1.2488987190  0.0115780095  0.3038119142  0.6541306799  0.3437720940
##  [636]  0.7661933586  0.3523400504  1.2792477004  0.6392295709  0.8555542480
##  [641]  0.6010341915  0.7816956539  1.4983828105  0.3721136290  1.9911185442
##  [646]  1.3611801067  0.7758199985  0.6908745848  0.3645090622  0.7862356569
##  [651]  0.7934494954  0.4594616407  0.7833567699  0.8133526861  0.0604476077
##  [656] -0.4261013360  0.3649181421  0.5871860274  1.2524845039  0.7732790825
##  [661] -0.3880213519  1.9966169423  1.3332715258  0.8332858304  1.2340339955
##  [666] -0.8573656533  1.3527563849  1.2108560969 -0.0467711130  0.4125159963
##  [671]  1.1229705772  2.2331632731  0.7410450334  0.8086251358  0.6847075679
##  [676]  0.2813871670  1.3518722548  0.4071340629 -0.0291095224  1.3527138223
##  [681]  1.3252691995  1.3611801067  2.0801487181  0.4313631777  0.7401799217
##  [686]  0.3507733720  1.0670494646  0.7811160686  0.7404214047  1.2628007089
##  [691]  2.3859261148  0.3551976025  0.6504008740  0.8207435799 -0.4488841543
##  [696] -0.0984621801 -0.0995037724  2.3469705133  1.2993629772 -0.4263957363
##  [701] -0.1038657077  0.3604377633  1.4263427510  0.6438938167  1.2558739596
##  [706]  1.7605297748  0.3546843906  0.3114373636  1.1940188898  0.3561112107
##  [711]  0.8676294484  0.8036619040  0.3533976549  0.3553607538  1.4138225323
##  [716]  0.8563958520  0.3042407670 -0.0414134435  0.0083940707  0.0450894779
##  [721]  0.0370993057  1.5662801658  0.0167066085  1.8776595023 -0.4582251320
##  [726]  1.7041366193  0.2584207983  0.7082490893  0.2754906266  1.3629921973
##  [731]  1.3984463872  0.7656757092  1.7430331247  2.0031921802  0.3999851887
##  [736]  0.2503978406  0.8372347371  0.7276122801  1.3671166468 -0.0518175112
##  [741]  1.3031465613  1.0225471420  0.4148842184  1.2548574065 -0.0771823848
##  [746] -0.0325869021  1.3613690748  2.0972319663  1.0696275459  0.3752761822
##  [751]  0.7366515416  0.7226169774  1.0884746062  0.7143943071  0.7276990598
##  [756]  0.6729907352  1.2542289609  0.6407504700  0.7297777728  1.9179764036
##  [761]  1.4322977876  0.8040121786  1.2603311565  1.3613630484  0.4042530849
##  [766]  2.0342616268  0.8345199667  1.6454518384  0.7988675371  0.2320576850
##  [771] -0.1035931532  0.7719765477  0.7194638826 -0.0306403286  0.7156690556
##  [776]  0.4063779675  0.8886700702  2.0178819332  1.9420517210  0.8205266312
##  [781]  1.2952120363  0.7645357037  1.1277479740  0.0064468629 -0.0188763222
##  [786]  1.1837536262  0.8399771591  0.8461942085  1.2977807701  1.3624562413
##  [791]  0.7015376957  1.8421963710  2.3255199355  1.3362923107 -0.0206051016
##  [796]  0.0211925057 -0.0436725842  1.5543761915 -0.2712358732  0.6733096300
##  [801]  0.7614362715  0.6357955014  0.3785202886  0.8427224613  0.6445369453
##  [806]  0.3238473297  1.4719153734  2.1574811885  1.2505589961  1.4249832854
##  [811]  0.5698825140  1.3815334125  0.3274127264  0.7924660300 -0.5840296180
##  [816] -0.0970660394 -0.3369862895  0.2049892773  0.7484111638  1.2571002477
##  [821]  1.4103477158  1.5573405853  1.9452942193  1.3876882595  2.2231360519
##  [826]  1.2555954643  0.6452671859  1.2812098458  1.3381619393  0.3849133119
##  [831] -0.0401828296  0.7583713404  1.1569412825  0.8447157666  0.6287321917
##  [836]  1.3679107532  0.4294124850  0.3272573672  0.6243760938  0.6345334588
##  [841]  1.4071224502  0.7086052772 -0.3560554561  0.0323313817  0.6620282459
##  [846]  0.3942444692  0.7921056584  0.6621413258  0.6922501647  1.2582594370
##  [851]  2.3485649023 -0.0022052525  0.5059056891  0.7432596221  0.2821655709
##  [856]  0.4091202767  0.3254116697  0.0298431957  1.0106980090 -0.4210800321
##  [861]  0.4038253630  2.1860088189  0.7945329071  1.2279293232  1.9293416262
##  [866]  0.8047092791 -0.0306403286  0.3707768470  0.3680763915  0.8199810857
##  [871]  1.8393702400 -0.0029764937  0.0101900814  0.3996491931 -0.4125187462
##  [876]  1.2242237877  0.2511391948  0.7227354969  1.4582416134  1.5315895781
##  [881]  1.3565274360  0.2888874394  0.7059443714  2.0311406267 -0.0745478340
##  [886]  1.2732650442  0.3766573811  0.6591928826  1.2444798495  0.3456502688
##  [891]  0.7777464800  1.3472690222  1.1960140642  0.6973385363  0.3221405856
##  [896]  0.5712015812 -0.4495351573  1.1516997911  0.4077012612  0.4191990517
##  [901]  0.9852113012 -0.0171024136  0.0141175712  1.8730468348  1.4158292808
##  [906]  0.5555569682  0.7532812950  0.8306091086  0.6216690245  1.3022670225
##  [911]  0.3786033746  0.8259644562  0.3228374263  1.6701065398  1.3127979620
##  [916]  0.6819956749  1.8203621001  0.6936288718  1.2611376687  0.9699244807
##  [921]  0.7254513765 -0.3815485192  1.2824016394 -0.0630118364  0.7384521228
##  [926]  0.3825971445  0.3381059619  1.8084682969  1.1973183680 -0.7263015584
##  [931]  0.1273569416  0.3261747985  1.4588555175 -0.4999607734  0.3741196160
##  [936]  2.1217235446  0.7234393179  0.8013874679  0.3937082402  0.7271877075
##  [941] -0.0588005952  1.0180378929  1.8933510299  0.9620463087  0.5610872128
##  [946]  1.3139481804  0.8185301261  1.2266861121  0.4116523551 -0.1028290939
##  [951]  0.7568342020  1.5859138482  1.0331580984  1.2471802783 -0.3980288412
##  [956]  0.0045479478  2.4635126666  1.3897901115 -0.4373534562  0.4235492530
##  [961]  1.9730396009  0.8157562019  2.2170627480  0.0271727729  0.4093612994
##  [966]  2.0470524483  0.3439922118  0.8017884555 -0.3954080102  0.6773940206
##  [971]  1.9634543340  1.9216077894  0.3038926076 -0.0902487904  0.8137298679
##  [976]  0.3545617810  1.7739550419 -0.0212809489  1.2289028244 -0.8880284890
##  [981]  1.2003933936  0.7401591891  1.1051342353 -0.0562069189  0.3022237808
##  [986]  1.9079176570  0.3525637700  0.3386351014  0.0475483671  1.1666619340
##  [991]  0.8071457756  0.7907819415  1.2777759396  0.2433642560 -0.0117638214
##  [996]  0.7115373877  0.2913975631  2.1910022288  2.0322783451  1.5757313651
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
```

```r
fit2
```

```
##       mean          sd    
##   0.54542989   0.25394198 
##  (0.08030351) (0.05677983)
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
## [1] -0.2358690  0.4965036  0.3763614 -0.2909760  0.4712937  1.1192304
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
## [1] -0.0194
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8953353
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
## t1*      4.5 -0.00970971   0.8939318
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 2 3 4 7 8 9 
## 1 1 2 1 1 3 1
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
## [1] 0.0047
```

```r
se.boot
```

```
## [1] 0.8906296
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

