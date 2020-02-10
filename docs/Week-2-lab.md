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
## 0 2 3 4 8 9 
## 3 2 1 1 1 2
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
## [1] -0.0029
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
## [1] 2.64039
```

```r
UL.boot
```

```
## [1] 6.35381
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.7   6.4
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
##    [1] 6.0 3.5 3.0 5.0 4.3 4.0 4.5 4.6 5.4 4.2 5.1 2.6 4.9 4.7 4.4 5.1 5.1 5.0
##   [19] 5.3 4.0 4.3 3.1 3.3 3.3 2.8 3.7 5.3 3.4 3.0 4.8 5.5 4.2 3.2 4.7 3.5 4.7
##   [37] 5.6 4.5 2.7 3.6 3.7 5.0 3.4 4.3 4.2 4.8 5.0 5.9 4.9 5.1 4.4 4.8 4.5 4.3
##   [55] 4.8 4.2 4.2 5.0 4.0 4.3 4.5 3.5 2.9 4.1 4.2 4.1 4.2 2.9 1.6 5.9 4.2 4.7
##   [73] 5.0 4.7 4.8 3.5 2.0 4.9 4.4 5.3 6.3 6.7 3.6 3.8 5.3 3.3 3.8 4.9 5.0 4.1
##   [91] 4.7 3.5 4.2 3.8 4.7 2.9 3.4 3.4 4.9 5.7 5.2 4.9 4.3 2.7 4.8 4.6 5.1 4.1
##  [109] 4.9 3.9 5.3 3.7 5.3 5.6 5.2 4.7 5.0 5.6 5.0 4.6 3.6 5.0 4.7 3.9 4.7 4.8
##  [127] 4.1 4.2 5.6 3.3 6.0 4.3 4.2 5.2 5.8 5.3 4.6 4.1 4.0 5.2 3.4 6.1 5.6 2.0
##  [145] 3.2 4.8 4.8 3.6 4.1 4.9 3.9 3.5 5.0 4.3 6.6 5.8 4.8 3.9 3.7 4.6 5.2 5.5
##  [163] 3.4 3.7 3.7 3.2 2.9 4.6 4.4 4.8 3.9 4.7 4.1 5.7 4.9 5.7 2.9 5.9 3.7 3.9
##  [181] 3.5 6.0 4.3 5.6 4.8 5.1 4.9 5.3 4.4 3.9 5.1 4.5 4.6 5.5 4.2 4.2 3.3 3.3
##  [199] 5.7 5.1 5.3 4.0 4.2 5.7 4.8 3.9 2.8 4.1 3.9 4.8 5.4 5.1 5.1 4.8 2.8 3.6
##  [217] 4.8 4.5 4.2 5.6 5.2 5.0 5.3 3.6 5.7 4.8 4.9 3.6 4.3 5.9 3.9 2.5 4.0 4.3
##  [235] 5.2 3.6 5.9 3.8 3.9 4.3 2.9 3.7 3.3 3.9 4.4 5.7 3.6 4.5 5.0 3.7 4.7 3.8
##  [253] 3.8 5.8 3.4 5.4 3.5 5.5 5.8 3.4 4.4 2.1 5.2 3.3 2.3 5.0 4.3 4.6 4.0 2.9
##  [271] 4.3 2.1 5.2 4.8 6.3 5.2 2.5 3.4 3.2 5.0 4.3 4.8 3.7 4.0 4.9 4.7 4.9 4.5
##  [289] 5.9 4.7 3.2 4.2 3.1 5.0 4.1 5.3 4.1 5.6 5.2 5.1 6.2 3.8 3.7 6.1 4.5 4.1
##  [307] 4.9 5.0 4.5 4.0 5.1 4.5 4.4 3.8 4.2 5.1 3.7 3.2 6.4 4.3 2.7 4.0 4.2 5.7
##  [325] 4.8 4.9 5.0 3.9 2.4 5.9 5.3 5.4 3.3 4.4 5.6 5.3 4.3 5.5 5.0 5.0 4.0 4.4
##  [343] 5.6 6.1 4.6 4.8 3.6 4.3 3.8 4.6 4.6 3.8 4.7 4.7 4.1 4.9 4.7 4.4 4.2 3.7
##  [361] 4.6 4.1 5.0 4.5 6.0 4.6 5.3 5.3 5.8 4.6 5.5 5.9 4.7 3.6 4.1 5.1 5.1 3.0
##  [379] 3.0 4.2 4.5 4.7 2.1 4.7 2.7 4.5 4.4 3.2 4.1 5.0 4.2 4.0 4.3 3.5 5.0 2.5
##  [397] 5.7 6.1 5.4 5.2 5.8 4.4 3.7 2.5 5.2 5.2 6.5 5.8 5.2 4.3 2.8 4.4 4.6 5.0
##  [415] 4.6 3.6 4.1 2.8 4.0 5.2 4.1 7.0 4.7 4.9 5.0 5.3 4.4 5.7 4.3 4.0 5.1 3.9
##  [433] 4.3 6.6 6.0 5.4 5.9 5.7 4.8 3.8 5.1 4.6 4.6 5.0 4.4 4.4 3.5 3.4 4.6 3.1
##  [451] 4.5 4.6 4.6 6.3 3.6 4.0 6.1 2.8 3.5 5.2 5.1 3.2 4.9 5.0 4.3 3.9 6.0 2.2
##  [469] 4.2 4.0 4.0 3.6 5.3 4.1 4.2 5.6 3.3 3.0 6.5 3.0 5.6 5.1 2.4 6.0 5.1 5.4
##  [487] 5.7 4.9 5.4 4.8 4.6 4.5 6.6 4.9 4.8 4.6 4.9 5.2 4.8 3.8 4.0 4.3 5.8 5.4
##  [505] 3.7 4.0 4.7 4.5 4.6 4.2 4.7 4.9 3.5 2.7 4.9 4.5 3.7 4.6 5.7 4.9 3.6 6.7
##  [523] 4.5 4.7 3.6 5.6 4.2 4.5 4.0 3.9 4.6 4.6 3.8 3.5 4.0 4.1 5.5 4.2 5.6 4.8
##  [541] 5.1 4.4 4.7 5.2 4.1 2.6 4.9 5.5 4.6 5.2 4.0 5.0 2.7 4.2 5.2 4.0 5.2 3.7
##  [559] 3.8 6.5 4.0 4.2 4.7 4.5 5.7 3.6 6.5 4.8 3.4 3.4 5.3 2.4 5.3 3.5 3.9 4.2
##  [577] 4.1 4.2 3.1 4.8 3.0 4.9 4.7 4.4 5.4 5.8 3.3 6.1 5.3 5.1 5.3 4.8 3.7 3.6
##  [595] 2.9 5.7 3.7 6.1 5.4 5.1 4.2 1.6 4.5 4.7 4.0 4.1 4.5 4.2 5.4 3.6 5.0 4.7
##  [613] 4.2 3.5 5.2 4.4 5.2 6.4 3.6 3.9 5.8 5.4 5.0 4.4 3.5 5.3 3.4 4.0 5.3 3.2
##  [631] 4.3 4.5 5.2 3.5 4.4 4.3 4.8 5.3 5.6 5.7 4.1 4.7 4.2 5.3 4.4 4.5 4.0 3.8
##  [649] 5.1 4.4 3.3 5.2 4.0 3.2 4.7 3.3 3.1 5.0 3.3 3.4 4.6 4.4 5.5 5.4 4.3 4.6
##  [667] 4.9 3.3 4.3 4.5 4.6 5.1 2.2 3.2 3.3 6.0 4.5 5.0 3.9 4.4 5.2 4.7 4.3 4.7
##  [685] 4.8 3.5 4.5 4.9 4.6 5.5 4.0 4.4 2.9 4.8 3.5 4.9 4.5 4.1 4.9 3.8 3.6 3.4
##  [703] 2.3 4.1 3.0 5.5 4.9 5.7 4.6 4.6 4.2 3.4 3.4 4.4 5.3 6.0 4.3 4.4 5.2 3.4
##  [721] 4.4 3.8 6.0 4.1 5.1 3.0 4.0 4.9 3.6 6.0 4.1 5.1 5.7 5.5 4.8 4.2 3.7 3.3
##  [739] 4.2 4.0 5.0 4.9 2.5 5.2 5.9 5.1 4.1 4.6 5.2 4.9 2.9 4.7 5.2 4.9 4.5 3.9
##  [757] 3.8 1.8 5.4 5.5 3.8 3.5 5.8 5.5 4.1 5.8 5.5 4.3 5.0 3.4 2.7 3.7 5.9 3.7
##  [775] 5.1 4.5 5.4 4.7 5.3 4.1 4.0 5.4 4.9 3.4 5.0 5.3 2.5 6.0 3.6 5.9 4.6 3.8
##  [793] 3.5 3.6 4.6 3.6 4.8 3.5 4.7 5.3 4.4 4.0 4.5 5.2 2.9 4.4 5.3 4.5 5.5 4.8
##  [811] 4.2 5.2 4.4 4.6 4.3 4.0 3.5 5.0 4.6 3.6 4.5 3.6 5.4 4.3 3.4 3.4 4.0 3.5
##  [829] 5.3 6.0 4.9 5.2 5.6 4.7 4.2 4.8 4.5 4.8 4.3 4.3 3.2 2.6 4.5 5.0 4.0 5.0
##  [847] 4.5 4.1 5.6 4.7 5.7 4.1 5.0 3.0 5.6 3.3 4.1 4.5 5.5 4.9 3.1 6.2 2.3 4.9
##  [865] 3.0 5.0 5.4 4.0 4.8 3.4 5.0 3.0 4.6 7.0 3.9 6.2 3.2 5.0 5.0 5.4 4.8 4.9
##  [883] 5.1 3.3 5.3 2.3 4.7 4.3 3.3 4.5 5.0 3.6 4.2 5.7 5.0 3.9 5.2 3.7 4.0 4.8
##  [901] 5.3 4.6 5.8 5.5 5.3 5.1 4.2 5.1 4.5 4.5 2.4 5.0 3.5 4.7 5.9 4.6 4.8 4.3
##  [919] 3.2 3.9 5.7 4.0 5.7 4.7 5.3 4.1 4.1 3.3 3.4 4.1 3.9 5.3 4.5 6.3 4.5 4.9
##  [937] 4.6 4.2 4.5 2.8 5.3 5.5 4.2 5.0 4.1 2.5 5.2 4.4 3.0 4.3 6.1 5.9 5.9 4.5
##  [955] 5.2 3.5 5.8 4.0 4.1 4.5 2.7 4.4 4.6 4.6 4.0 4.3 5.5 3.9 4.8 2.5 6.2 4.1
##  [973] 4.5 4.0 3.2 4.4 5.3 4.5 3.1 5.6 5.2 5.2 5.2 6.4 4.9 4.5 6.2 4.9 5.3 3.9
##  [991] 5.6 4.1 5.1 3.2 2.6 4.5 5.2 5.2 5.6 3.4
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
##   2.5   6.1
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
##    [1] 4.8 4.3 3.8 5.5 4.0 3.4 5.2 4.4 2.6 4.5 2.8 5.9 3.9 4.9 5.0 3.9 4.3 6.1
##   [19] 3.9 4.6 5.3 2.5 4.1 3.7 2.3 4.6 5.0 3.8 5.7 4.9 5.9 5.1 5.6 3.6 5.0 4.7
##   [37] 4.7 6.0 4.3 5.8 6.4 4.8 3.5 4.3 4.0 5.4 5.9 4.6 3.9 4.2 5.9 5.6 3.5 3.4
##   [55] 5.0 3.8 5.2 4.1 2.9 3.9 4.2 5.6 4.7 4.2 3.4 4.4 3.0 4.4 5.3 5.1 4.4 4.6
##   [73] 6.2 3.6 4.7 3.8 4.1 4.6 3.9 5.0 4.2 5.1 5.3 4.3 3.3 3.4 5.1 3.2 5.7 4.0
##   [91] 5.9 5.1 4.7 3.4 4.7 5.7 3.9 3.5 5.8 3.9 4.0 5.0 5.3 4.3 3.2 2.7 2.6 4.1
##  [109] 3.3 4.3 3.9 3.8 4.6 5.5 4.8 4.2 3.8 5.0 4.6 4.3 3.5 5.1 3.4 4.0 6.2 4.4
##  [127] 4.7 4.0 2.9 5.4 5.8 6.0 5.7 3.4 4.4 4.7 4.8 3.5 2.7 2.5 4.7 4.6 3.8 4.2
##  [145] 4.6 3.5 4.0 4.9 5.9 5.0 4.5 4.2 3.6 4.3 4.8 3.8 4.5 4.8 4.2 5.4 2.1 3.5
##  [163] 3.8 4.1 2.9 3.2 5.2 5.8 3.4 4.9 4.6 4.4 5.1 4.7 4.5 5.6 4.5 5.4 6.1 4.0
##  [181] 2.8 5.6 4.4 4.1 5.0 5.6 4.0 4.2 4.8 3.2 5.7 8.0 4.5 3.9 5.3 2.6 4.8 4.6
##  [199] 3.6 4.1 5.8 5.8 4.6 4.5 3.7 5.3 5.1 6.8 5.3 4.8 2.4 3.9 4.5 4.6 4.0 4.9
##  [217] 4.6 5.4 4.3 4.8 3.6 5.2 5.3 4.3 4.7 5.2 4.2 3.0 4.4 3.4 4.5 3.7 5.3 5.0
##  [235] 4.5 5.1 4.7 4.1 4.3 4.3 4.2 3.1 3.1 5.4 3.3 5.4 3.8 4.1 3.0 3.5 5.1 3.6
##  [253] 3.4 3.3 3.5 5.9 5.3 3.3 4.3 3.5 3.8 3.8 2.5 3.4 3.7 3.5 4.3 2.5 3.2 5.2
##  [271] 4.0 4.2 4.1 4.6 4.1 3.8 4.4 5.6 4.9 4.7 4.3 5.3 5.5 5.3 3.3 4.0 5.2 5.1
##  [289] 4.9 3.8 5.4 4.3 2.4 4.2 4.5 5.8 5.3 3.7 4.6 3.1 3.8 4.1 6.0 4.4 4.9 5.9
##  [307] 4.9 4.5 6.2 4.5 6.1 5.1 3.9 4.6 6.1 4.5 5.1 3.3 4.1 5.2 4.2 4.0 3.9 4.7
##  [325] 5.6 4.5 5.8 4.5 3.6 3.6 4.4 4.0 5.0 4.8 3.9 5.3 4.8 5.6 6.6 4.7 3.5 5.1
##  [343] 6.1 5.0 5.3 4.8 4.5 5.6 6.1 3.3 5.2 5.8 4.5 4.2 4.6 5.9 4.3 4.5 3.9 3.8
##  [361] 3.2 5.6 4.7 2.9 3.9 4.7 3.4 3.2 4.5 4.4 2.8 5.0 4.3 4.4 4.2 4.6 5.4 5.1
##  [379] 2.0 4.6 6.1 4.4 6.7 4.1 4.2 3.9 4.2 4.4 4.7 5.6 3.7 4.6 4.7 5.0 4.8 6.1
##  [397] 5.1 5.9 3.8 3.9 1.9 3.5 4.0 3.2 5.9 3.9 3.1 4.7 4.6 4.7 4.9 4.2 4.7 4.2
##  [415] 5.3 3.8 4.8 5.0 5.3 5.2 6.3 4.5 4.5 5.6 3.4 4.4 3.9 2.9 3.3 3.3 4.3 3.8
##  [433] 2.9 5.0 4.7 4.4 3.7 5.2 5.2 5.3 4.1 3.0 3.9 5.2 3.6 3.1 4.6 5.1 5.0 4.0
##  [451] 5.2 5.1 4.3 4.8 4.6 4.0 5.7 6.1 3.9 4.9 6.1 5.9 4.4 4.4 6.1 4.1 2.6 6.0
##  [469] 5.0 4.1 4.8 3.5 6.0 5.1 3.8 4.9 5.2 4.3 4.1 4.3 4.7 4.9 5.1 4.7 4.6 4.5
##  [487] 2.8 2.5 4.2 5.4 4.8 4.0 4.5 3.9 5.0 4.0 4.1 5.2 4.9 5.8 3.9 3.6 5.1 3.7
##  [505] 3.7 4.2 5.5 4.9 3.9 4.7 3.8 4.3 3.4 6.2 3.4 3.2 4.5 4.2 3.3 5.0 4.3 6.5
##  [523] 5.0 4.9 4.1 5.1 4.6 3.8 4.7 3.1 4.1 5.4 5.3 4.8 4.9 5.0 4.9 3.4 5.4 5.5
##  [541] 4.8 4.0 2.4 5.6 6.0 3.5 3.0 4.2 5.5 6.1 5.4 6.2 4.9 4.2 3.1 6.3 3.3 6.8
##  [559] 3.6 3.0 3.3 4.3 3.6 6.2 3.4 5.8 5.2 4.3 2.8 4.5 5.2 5.4 5.1 4.6 5.1 5.0
##  [577] 5.4 3.2 3.2 4.2 4.5 5.3 4.5 4.5 2.6 4.1 5.3 4.9 4.8 5.9 4.5 4.0 5.6 6.5
##  [595] 5.2 6.2 4.7 4.2 4.8 2.4 5.7 5.1 4.8 4.1 5.1 5.2 4.6 5.1 3.1 4.1 4.9 6.7
##  [613] 5.1 4.3 4.8 4.7 4.3 4.5 5.1 5.4 3.5 3.7 4.6 4.3 2.2 4.4 5.4 3.0 2.8 5.5
##  [631] 2.8 3.2 3.7 5.3 4.2 3.2 6.5 4.6 3.8 4.8 4.6 4.3 3.5 4.8 3.7 5.4 5.2 4.6
##  [649] 5.4 5.6 4.9 5.0 4.9 3.0 3.3 5.0 4.8 4.4 4.5 4.1 4.7 4.5 3.3 4.4 4.3 5.5
##  [667] 4.6 4.8 4.5 4.0 3.6 3.7 3.4 5.3 3.6 4.7 5.3 4.6 2.5 3.4 5.2 6.2 3.7 4.3
##  [685] 3.6 4.8 4.3 5.9 5.0 6.1 5.8 4.9 3.8 4.6 3.4 4.0 4.1 3.6 4.8 3.1 4.5 4.9
##  [703] 5.5 6.2 3.6 4.0 4.3 3.1 4.3 3.9 5.0 4.6 2.9 3.5 4.1 5.5 3.5 5.4 5.9 4.2
##  [721] 3.7 4.2 4.8 5.2 5.3 4.5 5.4 4.9 3.3 2.8 5.8 3.9 4.8 5.4 6.2 4.5 5.5 4.5
##  [739] 3.3 6.3 4.3 5.0 4.9 4.5 3.7 3.3 2.7 5.4 5.0 4.7 4.8 5.6 2.7 2.8 3.9 4.1
##  [757] 4.3 6.8 3.0 4.3 4.8 5.6 4.1 3.4 4.4 6.3 4.5 3.4 5.5 3.1 3.7 6.0 3.9 5.2
##  [775] 5.5 3.3 4.6 4.2 4.7 4.5 4.9 5.8 2.2 4.7 6.1 7.1 4.7 5.8 4.3 3.4 5.4 2.9
##  [793] 4.4 3.9 4.4 4.9 6.1 4.7 4.8 3.6 3.6 4.8 3.7 5.7 5.3 5.0 5.6 5.0 3.2 3.6
##  [811] 5.2 4.7 4.2 5.4 2.9 3.8 4.9 3.8 5.1 4.6 5.8 2.7 3.6 4.7 3.7 4.6 3.3 3.3
##  [829] 5.9 4.5 5.9 4.7 3.2 4.3 3.8 6.3 6.6 6.1 4.8 4.0 4.9 5.8 4.6 6.0 4.8 5.1
##  [847] 4.4 6.0 4.6 4.9 3.5 5.1 5.0 4.5 6.2 2.8 6.4 5.6 5.3 4.5 3.3 4.5 6.0 3.7
##  [865] 5.1 4.9 2.8 4.9 3.3 2.2 5.2 6.0 5.5 4.9 4.7 4.6 3.8 3.7 5.7 5.1 3.9 4.4
##  [883] 4.4 3.7 5.3 5.3 4.0 3.7 5.6 3.4 3.7 4.1 3.8 4.7 5.3 3.0 4.8 5.6 4.5 4.5
##  [901] 5.7 5.0 4.9 6.0 4.6 4.1 5.8 5.4 3.6 5.9 4.8 3.7 6.0 5.6 5.5 5.4 5.7 3.8
##  [919] 4.0 3.4 3.4 4.4 5.1 3.1 5.5 4.6 4.8 3.1 4.4 4.9 4.8 5.1 6.9 5.3 2.7 3.6
##  [937] 5.4 3.8 4.0 4.0 4.7 5.2 3.2 5.3 5.6 4.7 5.0 6.1 4.8 4.1 4.4 4.3 4.6 4.6
##  [955] 4.8 4.3 6.5 3.8 3.0 5.3 2.8 3.5 4.0 3.6 5.3 3.2 4.4 4.5 3.4 5.3 3.8 3.1
##  [973] 6.0 4.0 5.2 5.0 5.8 3.7 3.4 4.1 5.1 3.7 3.7 5.4 4.2 4.6 4.0 3.1 5.0 5.8
##  [991] 4.4 3.9 3.2 4.4 3.0 4.5 5.4 5.6 3.3 3.7
## 
## $func.thetastar
## [1] 0.0023
## 
## $jack.boot.val
##  [1]  0.50269542  0.35045593  0.26564246  0.20054348  0.04680233 -0.09021739
##  [7] -0.16825397 -0.31732955 -0.35515320 -0.54620462
## 
## $jack.boot.se
## [1] 0.9715941
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
##    [1] 5.2 5.4 4.5 4.7 3.2 5.5 4.5 4.0 4.3 4.7 3.9 4.2 2.9 6.0 3.8 6.5 4.6 3.4
##   [19] 4.9 4.7 5.8 4.9 2.9 6.6 4.1 3.4 4.2 4.9 4.7 4.4 4.0 5.2 5.8 4.7 4.6 4.8
##   [37] 5.4 3.7 6.1 5.6 4.5 5.7 6.9 4.8 4.3 5.2 5.2 4.8 5.7 5.7 4.8 5.7 5.6 4.8
##   [55] 5.6 6.2 7.1 5.4 5.5 3.0 5.3 3.5 4.8 4.3 5.2 5.3 4.2 4.0 4.4 5.1 4.3 3.1
##   [73] 5.9 3.9 4.8 4.6 3.4 3.4 3.2 3.5 4.8 5.2 6.2 4.3 3.5 2.8 6.0 5.1 4.2 4.8
##   [91] 4.3 6.9 4.0 3.7 4.3 3.8 5.1 4.7 3.4 5.8 4.8 3.5 5.6 3.2 3.9 5.5 4.2 4.5
##  [109] 3.5 4.2 5.9 5.4 3.0 3.6 4.5 4.7 4.7 5.2 3.9 4.4 4.7 4.0 5.1 4.4 4.3 3.8
##  [127] 2.8 6.1 2.3 4.9 4.8 4.0 4.5 4.8 3.9 4.5 2.5 3.9 4.6 5.0 4.1 5.2 6.1 4.0
##  [145] 4.9 4.2 3.2 5.1 4.5 3.5 5.3 4.6 5.3 3.8 3.2 6.2 3.9 6.4 2.5 5.1 4.6 4.2
##  [163] 5.5 5.1 4.1 3.0 4.9 5.0 5.1 4.4 5.2 4.6 5.2 4.0 3.7 5.3 4.9 3.1 2.7 3.6
##  [181] 5.6 4.0 3.3 4.2 3.7 5.8 3.6 4.3 5.7 3.8 2.8 3.6 4.8 4.0 5.0 5.5 6.4 2.9
##  [199] 5.2 5.1 5.2 4.4 3.2 4.5 3.9 5.2 3.7 5.8 4.0 3.0 5.5 4.1 4.9 3.4 5.4 4.3
##  [217] 4.5 4.7 4.7 3.4 3.6 5.7 4.4 4.1 4.7 3.5 5.2 5.8 4.1 6.1 4.8 4.9 4.4 3.1
##  [235] 4.8 3.5 4.6 5.2 5.2 4.9 5.0 4.3 3.2 5.4 4.8 5.0 4.8 4.0 6.1 3.2 5.1 3.3
##  [253] 5.2 3.8 4.1 4.9 2.5 4.0 3.5 5.5 6.4 4.4 4.1 3.7 4.7 5.3 6.1 3.9 5.0 3.0
##  [271] 3.6 4.6 4.6 3.9 4.1 4.1 6.5 4.5 4.0 4.0 4.4 4.4 4.4 3.8 5.5 5.3 3.5 4.7
##  [289] 3.7 4.4 3.3 5.8 5.0 4.6 4.9 4.7 3.3 3.1 5.3 5.5 4.4 3.7 4.6 5.3 4.2 3.1
##  [307] 6.5 4.6 4.7 4.6 4.3 3.5 6.1 5.2 1.6 4.6 4.8 3.7 3.2 4.8 5.1 4.0 4.2 5.7
##  [325] 4.5 4.7 3.6 4.6 3.9 4.1 4.0 5.1 4.6 4.0 3.8 5.4 4.2 3.6 4.4 3.7 4.0 4.5
##  [343] 4.8 4.0 5.6 4.7 6.6 2.9 4.1 3.6 4.8 4.1 4.6 3.1 6.6 3.7 5.1 4.8 4.5 4.7
##  [361] 5.1 4.4 6.0 5.0 3.8 5.4 4.8 5.0 5.9 5.9 3.8 4.6 5.0 3.7 4.0 5.0 3.4 5.6
##  [379] 4.6 5.3 4.1 3.3 5.3 4.9 3.6 5.0 4.8 4.6 3.8 3.3 5.2 6.1 6.0 5.2 6.1 5.0
##  [397] 6.1 4.2 4.0 4.4 3.3 4.5 2.3 4.9 4.3 3.9 6.0 3.4 5.3 6.5 4.9 5.6 3.4 3.8
##  [415] 4.2 4.1 4.5 3.6 3.5 4.7 4.5 4.5 5.1 4.5 4.1 4.9 5.1 4.8 4.2 3.3 3.4 3.2
##  [433] 4.1 5.9 4.3 4.7 5.2 4.8 3.4 4.5 5.8 2.1 5.0 5.0 4.1 3.4 4.9 5.5 4.8 5.1
##  [451] 2.9 3.5 4.5 4.9 3.7 5.3 5.7 6.1 2.6 4.8 5.1 2.0 5.0 4.3 4.5 2.6 3.5 4.4
##  [469] 2.5 4.2 5.4 3.5 3.8 4.2 4.5 4.6 4.5 4.1 5.4 5.2 4.2 5.8 3.2 5.0 5.1 4.3
##  [487] 3.9 4.2 2.8 3.4 5.3 3.5 6.2 5.3 3.7 6.3 4.9 6.5 3.0 5.1 3.5 4.6 3.4 3.2
##  [505] 2.1 5.1 5.1 4.6 4.9 4.6 5.0 5.1 4.1 4.7 5.8 5.2 4.9 3.3 5.9 4.2 5.5 5.4
##  [523] 4.1 6.2 3.2 3.6 4.8 4.2 3.1 3.4 3.6 3.8 5.1 3.1 4.9 5.1 4.9 2.8 5.3 3.9
##  [541] 2.6 3.7 4.1 5.5 4.2 4.7 4.4 5.0 4.7 2.8 5.2 5.0 3.4 4.2 2.6 3.9 3.3 5.3
##  [559] 5.0 4.5 5.6 4.7 4.4 4.7 5.5 2.5 5.4 5.5 4.2 3.9 4.7 5.8 4.0 2.9 4.8 5.5
##  [577] 5.1 3.9 5.0 3.9 5.5 4.0 5.2 4.1 4.8 5.6 1.6 6.5 5.6 4.6 3.9 4.9 3.9 5.0
##  [595] 5.9 3.2 5.7 3.4 4.1 4.7 4.2 4.0 4.0 4.4 2.4 5.4 3.3 4.9 6.2 3.9 4.4 5.5
##  [613] 5.2 4.8 5.4 6.0 4.5 4.8 5.7 4.2 2.8 3.5 5.6 4.3 3.9 3.7 3.1 5.2 2.6 3.9
##  [631] 3.4 5.8 3.6 4.4 5.6 2.9 5.4 4.7 2.5 3.7 3.6 4.9 5.5 4.3 5.3 3.8 4.2 5.8
##  [649] 4.7 4.6 3.5 5.5 6.0 5.4 3.3 5.1 3.5 4.4 5.1 5.6 6.2 1.9 3.5 5.2 4.5 5.7
##  [667] 6.3 3.9 4.9 5.1 3.6 5.8 4.2 6.1 4.0 5.4 4.6 3.9 5.0 4.0 5.2 3.2 5.3 5.1
##  [685] 4.0 4.8 5.5 3.8 3.0 5.5 5.8 4.0 2.2 4.6 3.8 3.6 5.0 3.1 4.0 4.5 3.1 3.9
##  [703] 4.8 5.5 4.7 3.3 7.1 4.0 3.5 4.8 4.5 5.6 4.9 3.5 4.5 4.9 4.6 3.6 3.6 5.4
##  [721] 5.4 3.7 5.4 4.4 3.8 5.2 5.6 3.6 5.1 3.8 3.6 2.9 3.3 4.4 4.4 3.5 4.9 4.5
##  [739] 4.6 4.1 3.2 5.3 4.0 4.3 4.3 5.2 4.0 5.3 4.2 2.6 3.6 2.6 4.5 2.7 5.7 4.5
##  [757] 3.9 5.7 5.7 4.8 4.9 4.4 4.8 4.1 4.5 3.6 5.6 6.6 3.8 3.1 5.5 3.0 3.8 4.8
##  [775] 4.5 4.7 5.4 3.2 6.3 4.6 4.5 4.1 4.4 4.9 4.1 4.9 4.6 5.2 4.7 5.3 4.3 4.6
##  [793] 5.4 4.9 3.5 4.9 3.4 2.3 4.4 3.7 3.8 5.1 4.9 4.7 3.7 4.4 4.9 3.5 4.7 3.8
##  [811] 5.1 3.7 5.3 4.1 3.5 4.5 5.2 5.6 4.7 6.5 4.0 4.0 5.0 5.2 4.6 4.6 3.3 3.5
##  [829] 5.6 4.6 4.7 4.1 4.3 4.0 4.5 5.5 4.8 6.5 4.3 4.5 6.1 4.1 4.1 5.0 4.9 4.2
##  [847] 4.3 6.0 5.0 5.4 3.8 3.1 5.6 5.7 2.7 3.3 4.6 6.4 4.5 4.1 4.7 3.5 4.1 5.0
##  [865] 4.2 3.1 4.5 4.5 4.4 4.0 4.3 3.9 5.2 3.8 3.9 3.5 4.6 4.6 5.9 4.0 4.3 4.5
##  [883] 3.9 5.1 3.6 5.0 5.8 5.0 5.1 5.2 4.5 3.7 3.3 4.1 5.2 4.3 4.1 3.6 4.6 3.6
##  [901] 5.0 6.7 4.9 4.9 3.2 4.1 4.7 5.4 4.0 3.8 4.7 4.5 4.5 4.0 4.3 4.8 3.7 5.8
##  [919] 4.2 4.8 4.6 5.6 5.1 4.2 6.2 4.5 2.7 3.9 3.8 6.7 5.2 4.9 3.2 4.7 5.3 4.7
##  [937] 4.5 5.6 5.3 5.7 3.2 4.1 4.0 4.1 3.6 3.3 3.8 4.6 3.7 5.5 3.8 6.4 4.4 5.5
##  [955] 4.9 5.6 4.3 5.8 5.4 5.4 5.1 4.2 3.3 2.7 5.1 4.5 4.4 4.2 4.6 3.5 5.9 4.0
##  [973] 2.4 6.6 4.2 3.3 4.5 3.4 4.6 5.7 4.2 4.8 4.7 5.1 2.6 5.2 4.4 6.0 5.4 2.6
##  [991] 4.5 6.5 5.0 4.9 5.7 4.2 1.3 4.9 3.4 5.1
## 
## $func.thetastar
## 72% 
## 5.1 
## 
## $jack.boot.val
##  [1] 5.500 5.500 5.352 5.200 5.200 5.000 5.000 4.700 4.672 4.500
## 
## $jack.boot.se
## [1] 1.000647
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
## [1] 1.458408
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
##   1.8505809   3.6967563 
##  (0.7647198) (1.7527380)
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
## [1]  1.19994816 -0.05620637  0.39470486  1.29732042  1.23710866  0.73767613
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
##    [1]  1.1439519891 -0.3207833509  0.9988216011  1.4535568256  0.7136744389
##    [6]  2.0848061123  1.1493846236  1.0680963360  2.3307543752  0.6970797194
##   [11]  0.9503898905  0.8160030172  0.2500131596  1.5418307140  0.8056455439
##   [16]  1.8675889067  1.2605988707  0.9200177410  0.9636475673  1.1471584032
##   [21]  0.4050204450  0.3282517611  0.5684088696  0.4686441619  1.2611906059
##   [26]  2.1618918977  0.6735619879  1.2518312845  0.8298533063  1.1042350719
##   [31]  0.5047292018  1.4161682633  0.3113355370  0.9160597712  0.9118832927
##   [36]  1.7281394156  1.5729171597  0.3479277132  1.1958707446  0.2319954404
##   [41]  1.8711798330  1.4582313233  0.1887641340  1.7635634951  1.5503042417
##   [46]  0.4202382313  1.4116652583  1.7750064597  1.8247406360  1.9127616934
##   [51]  0.5018469974  0.7998696609  0.5289169315  1.3910616093  2.1578241337
##   [56]  0.0482347380  0.7493250362  1.0047320983  0.8877278517  0.9494869141
##   [61]  1.1377654480  1.7832061853  1.2379175615  1.3323484611  0.3893017438
##   [66]  1.2062263578  1.1610598399  0.5323759455  0.9760672572  2.1044413380
##   [71]  1.3542798566  1.6278795978  0.4251546926  1.0312672385  1.1332188917
##   [76]  0.8057832916  1.6910756835  0.7291777506  0.7594003604  1.2692954703
##   [81]  0.0353504125  2.3537180254  0.3219377722  1.3359485951  1.8854968855
##   [86]  1.6908109290  0.7870998105  0.3381762420  1.2010591460  1.5566219449
##   [91]  1.6948914905  0.9189710126  1.6157538784  2.1686810436  0.5763528312
##   [96]  0.6505043086  0.8842993061  1.1211137237  2.1470741245  0.9246144002
##  [101]  0.6936945751  1.0871495250  1.1994428433  0.8794158679  1.4279009311
##  [106]  1.1880582189  1.7959661508  1.4656519115  0.5668800607  1.5769357828
##  [111]  1.8985192156  0.2800921634  1.3126974066  0.7445954028  0.7430690992
##  [116]  0.5874094256  0.5635571400  0.3903502868  0.2306771676  1.1980251868
##  [121]  1.7484402244  0.2444007365  0.8946649992  1.5560552422  1.7125986265
##  [126]  0.3185673760  0.6312816391  1.5704644478  0.3881727796  1.8029006221
##  [131]  1.3819059460  1.4232573997  0.9689393573  1.1988199663  0.4656010711
##  [136]  1.2307448680  0.9330374159  0.0316584332  1.5690514559  1.3309214206
##  [141]  1.4049213059  1.6283152528  1.3305896901  0.5965800427  0.7226071883
##  [146]  1.1481139127  0.3067545894  1.1325707166  0.6483813420  1.5070194935
##  [151]  1.8213866445  0.6781360801  0.3817636432  1.5492483882  0.9211535505
##  [156]  0.0385903104  1.3766934642  1.4266982955  0.7240795075  1.2216606406
##  [161]  1.5878219480  1.1309675163  0.9096624769  0.7627370057  0.1406178528
##  [166]  1.0203628975  1.6331377649  1.4653905508  1.2467965793  1.5107967142
##  [171]  1.5508014132  0.4306949775  1.3629748252  1.6734536189  0.4154044070
##  [176]  0.5074057806  0.6653636659  0.6443213687  1.2885020228  0.6670563392
##  [181]  0.2975773880  0.6223398709  1.0445294239  0.1417660889  1.1383446359
##  [186]  1.4296286647  2.5976638051  1.9340197441  1.5517907353  0.8571899345
##  [191]  1.4584075135  1.3045645645  2.0418244936  0.9004603265  1.3163136307
##  [196]  0.7622179356 -0.0856638403  1.4760613920  0.5112895367  0.2236197509
##  [201]  1.2103556058  1.6608414332  1.1285009569  1.0115615760  1.1529550096
##  [206]  0.7967002130  1.0602981995  1.4243082917  1.1237553220  1.5654574830
##  [211] -0.0941588957  1.1789375251  1.4321600548  0.6773316883  1.6115344119
##  [216]  0.4982983971  1.0637908098  1.3186955805  0.0763182809  1.0573521357
##  [221]  0.8651701668  0.8897476158  2.1654453984  0.7316170786  1.0151313484
##  [226]  0.8163197326  0.8275769090  0.6464244533  1.9360667270  2.0327806194
##  [231]  1.3025967613  1.4246959811  1.6299791334  2.1037905787  1.3938454312
##  [236]  0.0450444946  1.4488084329  2.3363509316  1.1529788721  0.7067566997
##  [241]  1.5136844687  2.4004628495  1.0383233492  1.9766729206  0.7461349332
##  [246]  1.2398199824  2.1290744154  1.0734337214  0.5798561780  1.5637562738
##  [251]  0.6795352342  1.9649454379  1.8189716749  1.2656274076  1.1498374423
##  [256]  0.7246713296  0.3728426414  1.3467171491  2.1632811620  1.1365228415
##  [261]  2.2326435017  0.6559790427  1.6204263873  1.0670324504 -0.0385933770
##  [266]  1.4044023044  1.8946051202  0.3583956103  1.3085293099  1.3485894148
##  [271]  0.7364828563  1.2611906059  1.8018419407  1.2187231733  1.6499885591
##  [276]  0.3939260741 -0.0349388663  1.4033760152  0.9914069318  1.0138002350
##  [281]  0.8324652795  1.5567067702  1.2617941428  2.4274976016  1.5165997342
##  [286]  0.6520160834  0.9857148244  0.4563161401  0.8173713677  0.9197664459
##  [291]  1.9816391840  2.0024235849  0.2157952786  1.8872673633  1.9479451124
##  [296]  1.0868730446  1.8636468414  1.0871725724  1.0001043159  0.6129756142
##  [301]  2.1630368465  0.6603186371  1.2182586825  2.1618918977  0.7310636870
##  [306]  0.5617726568  0.4690957334  1.0987189651  1.0620028489  1.6997719081
##  [311]  0.8169418333  0.7184618369  0.9959765927  0.9062276541  0.3747801021
##  [316]  0.4936108414  0.5014795095  0.7929078281  1.1947135519  1.0392328820
##  [321]  2.3933379020  0.7024035246  0.6133501685  0.9835426777  1.0304672503
##  [326]  1.1967480451  0.9943430106  1.3514095402  1.0099886812  1.0023113786
##  [331]  2.1351644961  0.0426505068  1.3927238563  2.1355053357  1.1847670789
##  [336]  1.4981833102  0.2090721283  1.1722701872  0.4657068929  1.7391578905
##  [341]  0.6473451652  2.1956030629  0.4665179526  1.2424574273  2.1746007747
##  [346]  1.2146205825  0.1793819446 -0.0437790374  1.0576845454  1.4099630111
##  [351]  1.2783394345  1.4767715933  0.2284812110  2.4718612446  1.2183238975
##  [356]  0.9284054405  0.8531102054  0.4755496942  0.8100546063  1.2799312125
##  [361]  2.0596442301  0.5395010815  1.0999589882  1.6266505750  1.2280989593
##  [366]  1.9007704305  1.2799677248  0.2483203093  1.9567307207  0.1633329133
##  [371]  0.5937943113  1.0210371661  0.3729393464  1.5582154373  1.0240989365
##  [376]  1.3417134358  1.0991817499  0.8992303122  1.5113713938  1.4678781432
##  [381]  1.0475200445  1.4541642909  0.3571827895  0.5492378011  2.1600448692
##  [386]  1.4948531484  0.4553247081  0.9391976573  0.8695029712  0.6350370188
##  [391]  1.4740009033  1.8289420166  1.4908348906  1.8028016838  1.1325085137
##  [396]  0.7965302874  1.5227778669  0.1683329262  0.7996685523  0.6638829893
##  [401]  0.5823659916  0.4130492681  1.8908688500  0.7762617972  1.7832061853
##  [406]  1.1283667225  1.8908688500  1.4795668149  0.4572053397  1.4928228486
##  [411]  1.0286670335  1.0980336285  1.9852656920  1.5352644549  1.1658544973
##  [416]  1.3701864807  0.1500529703  1.5753161688  0.7953546924  0.5292976705
##  [421]  1.1913740924  0.6851196803  2.4471456286  1.0590130871  1.2940552358
##  [426]  0.4437703591 -0.2060669838  0.0915942428  0.7705024782  1.5504166627
##  [431]  0.7585050687  1.0445294239  0.6739039437  0.5086672121  1.2049344239
##  [436]  1.9289976247  1.4883569204  1.6203031021  0.8942404831  1.3077997001
##  [441]  1.7157690841  1.5327624998  1.2697993143  0.2872258352  1.7438428763
##  [446]  0.5433389432  0.8364998349  1.1510629175  1.1340838451  2.4655572902
##  [451]  1.1572305325  1.6166655883  0.5491748789  1.7154380278  1.4084666271
##  [456]  1.8192162880  1.0914940821  1.1403625704  1.8979348561  1.4732359521
##  [461]  1.2738158798  2.2259099279  0.7718861779  1.0573207808  1.4881628407
##  [466]  1.1899929748  0.5526940657  1.5467870249  0.6203766160  0.7236324010
##  [471]  1.4289780686  0.4022055628  0.2874737795  0.7334947485  0.9308769464
##  [476]  1.1184076371  1.3406757423  1.2492838460  1.8599203179  0.8669517196
##  [481]  1.3241492624  0.1652428752  0.3053502094  1.6689524466  1.0059154611
##  [486]  1.3812883075  1.4543475420  1.6836656974  1.3367357773  1.2326481242
##  [491] -0.1206807332  2.2017744887  1.2304643388  0.8378830451  0.7797079616
##  [496]  0.9142368226  0.0962646975  0.9263017010  0.6047220645  0.9288261856
##  [501]  0.8823215645  0.9298268300  0.5820774939  1.4850211212  1.2981115249
##  [506]  1.7141518818  0.9495549725  0.1056828173  0.3744598265  1.4412410259
##  [511]  0.7008111228  1.1967735024  0.0149045999  1.9982299554  1.3797909337
##  [516]  0.8734333205  0.6807334197  1.5892636589  0.9388985390  0.7886655638
##  [521]  1.4271628818  0.8209005350  1.5792766490  0.6039221462  1.5749021832
##  [526]  1.4369840236  0.6225362622  1.9382797105  1.1492399173  1.0542984703
##  [531]  0.2348972876  1.5970125198 -0.0942995751  1.6637530001  0.6333513095
##  [536]  0.9589013980  1.2039296709  0.4187033857  0.2066239199  1.2587600673
##  [541]  0.6245803531  1.3337430415 -0.0731455273  1.4527221364  1.0530664987
##  [546]  1.5083150131  1.2174268314  1.2107241624  1.7184114825  1.7029831411
##  [551]  0.6783075036  0.7008518568  1.0489784893  0.5426487867  0.9234300767
##  [556]  1.1477105877  2.0389340559  1.0250479357  1.2844991556  1.6241026587
##  [561]  1.1180480682  0.4392854929  1.6763076334  1.3292704597  1.6560554310
##  [566]  1.5070846879  1.3706251699  0.5015245308  1.2374484361  0.1278840586
##  [571]  0.8824013980  2.0924812570  1.0264123312  0.9865667653  1.7141518818
##  [576]  1.0577512839  0.1713149180 -0.0336108181  1.5300811721  0.4126338882
##  [581]  1.0698554512  0.7126889171  1.0846367677  1.0389997040  0.5378827307
##  [586]  1.1749145097  1.2017872014  2.3312141374  0.6119735909  0.5344458437
##  [591]  0.5696152091  1.2629328676  0.9514766903  1.1160018125  0.9620634869
##  [596]  0.4262694694  1.1099554727  0.7357910673  0.7411327131  1.2413223302
##  [601]  1.8365713754  1.4748611425  0.2376807494  1.1006384890  1.8985192156
##  [606] -0.0516150060  0.5626394400  0.5557776346  2.0362109242  0.6740591840
##  [611]  0.8039960269  1.8294134250  0.2111462040  1.6348923192  1.0681815986
##  [616]  1.7686813675  0.1158508696  0.7579866026  2.1423232323  0.5333387430
##  [621]  2.1141225286  0.8291180856  0.8058096564  1.6578369125  0.6182025746
##  [626]  1.0853051402  0.3460967596  1.1852401422  0.6096308733  1.8646267289
##  [631]  1.2223199232  2.1911667787  0.9807485542  0.8741904031  1.3423021570
##  [636]  1.0309047200  1.2899003042  0.7994717138  0.5631161510  1.2350278012
##  [641]  0.2439542899  0.9739951317  1.8006282512  1.0303122179  0.9148398388
##  [646]  1.0470301386  0.0377559749  1.0210982387  1.6266879447  1.3600604747
##  [651]  1.8701676317  0.8848867487  0.6438770982  1.0099886812  1.4055867058
##  [656]  1.1007179424  1.7762376113  1.9427444472  2.1482566072  1.8989935544
##  [661]  1.8240759081  1.8081110991  0.7252577485  0.5398353747  1.5272696055
##  [666]  1.8939071896  0.8796636307  0.9066503954  1.4153870707  2.2581211850
##  [671]  2.1500080428  1.0865039166 -0.4947535310  1.2352150200  1.5160305550
##  [676]  0.6392371136  1.5176624918  1.5320978115  1.1891360108  1.6188888087
##  [681]  0.6875128643  0.8583904090  0.4502315792  0.9326148382  0.8226149154
##  [686]  2.1411345490  1.1960665514  0.7248088611  0.7672609845  0.7527559496
##  [691]  1.8361802517  0.9463167036  0.8731436786  0.4732181621  0.9510606677
##  [696]  1.7517686216  0.1576554091  0.5718255208  1.0977419022 -0.0219564398
##  [701]  1.6956999849  0.9333088179  1.3969569330  1.5684938456  1.5774792984
##  [706]  0.5147256873  1.1193206000  0.5613396590  0.9019723664  0.2695571280
##  [711]  1.0636096090  1.1145150750  1.4287532555  2.0329670751  1.5937370643
##  [716]  2.1508971033  0.9161124290  0.7465266760  0.7081476443  1.2374298144
##  [721]  1.2210952475  1.1982763265  0.9661196896  1.6267928588  1.3656922643
##  [726]  1.7606556565  0.4865980115  1.0508307679  1.4325693991  1.6372628506
##  [731]  0.5818555038  1.3311087089  0.6996915924  1.2246391048  0.8581392958
##  [736]  1.1524848474 -0.0395893807  2.2532076254  0.5069407083  2.1481690210
##  [741]  1.1216760652  1.4036173080  2.2479857259  0.1240001029  1.0383583671
##  [746]  0.4885076802  1.7220295115  1.0517875513  0.3613806779  0.5097556835
##  [751]  1.0384747824  1.2642488263  2.1750224656  2.0264669921  0.8031774431
##  [756]  1.3539357033  2.2437435517  1.3711429240  0.9903893922  0.8260907355
##  [761]  0.2483834711 -0.0226741352  0.5334935857 -0.0336572066  1.5576878307
##  [766]  0.3608305446  0.4980321939  1.7145715656  0.1065469811  1.2110159162
##  [771]  0.5180991627  1.2996626111  1.0780629044  1.5945770931  0.5053028415
##  [776]  0.6593409275  1.3346802241  2.2919304453  1.1877580794  1.6180182931
##  [781]  1.0096274917  1.1766223273  0.8118210853 -0.1391636253  1.0370951045
##  [786]  2.0182780814  1.8972514714  1.4617016051  1.6324056444  1.4018424220
##  [791]  2.3075028616  1.5815146988  0.2371937793  0.8088379169 -0.1091239912
##  [796]  1.0178889463  1.4082407181  1.7853837256  1.1743815904  0.8193495266
##  [801]  1.4035967284  1.1555282763  0.7484260358  0.3176944240  1.2478026567
##  [806]  1.4697224940  0.9222457660  0.0831319216  1.4558502044  1.2822593212
##  [811]  0.5927862492  1.3736927923  1.4242681171  1.8425324755  1.3316353599
##  [816]  0.4586224491  0.9688872262  0.1649861238  1.2813570240  2.1484514121
##  [821]  0.7390477279  1.5227778669  1.2036441599  1.8660177815  0.0774303441
##  [826]  1.6327761793  0.8508725904  0.6187812691 -0.0882705994  0.3468334450
##  [831]  1.0374565841  0.4999073834  0.2752718111  0.6300115095  1.2985742335
##  [836]  1.0564539838  0.3007949272  1.7177901905  1.7546375570  0.8916903761
##  [841]  0.8559930124  1.1912310851  1.4058576888  1.3886767897  0.6893286111
##  [846]  0.3212598914  0.2044746043  1.1206573598  0.8195339149  0.4427796733
##  [851]  1.5231352744  1.2365872592  0.7793054373  1.1534836146  1.4995907534
##  [856]  1.4357291713  0.7319009500  1.0688649309  1.4236137044  0.8535658059
##  [861]  1.8660254479  1.5349663525  1.8845087671  1.7508271763  1.0905665243
##  [866]  1.0942966805  1.3673717379  1.3482458314 -0.0368533082  0.9814060166
##  [871] -0.1267782901  1.1019328415  0.4665179526  0.7067566997 -0.6886125548
##  [876]  0.7590224723  1.2662959417  1.5723204425  1.7543894444  1.7684993425
##  [881]  1.4463104936  1.5819569382  1.8186411595  1.6417278755  1.0459136583
##  [886]  1.2076924402  2.0900064995  1.2295711469  1.8752507940  0.3573127924
##  [891]  1.4193186954  1.2873142533  1.5489397317  0.8328868634 -0.1838697952
##  [896]  0.9539378507  0.6263355272  1.0319315723  1.3329872353  1.2279909716
##  [901]  1.3819059460  1.0115200758  0.8851892708  0.6581730656  1.7448470190
##  [906]  0.5485019208  1.2582719240  0.9208843288  0.7716415520  0.0626485044
##  [911]  0.5994505124  1.1026972641  0.6693267344  0.8493239270  1.1085837919
##  [916]  1.3488452220  0.4980321939  2.1010745854  1.0391614700  0.8345770810
##  [921]  1.0311882972  1.5226864608  0.3003827104  0.8718300550  1.6344994378
##  [926]  1.0692444769  0.9398904647  0.9065302254  0.1627591143  1.9958311689
##  [931]  1.9432347665  1.8702018971  0.8408521360  0.5824128528  0.9551701695
##  [936]  0.8686977079  1.2058841227  0.8090227192  0.5629108714  0.1893169340
##  [941]  0.4678088421  1.4766024342  0.1570838074  2.5056644574  0.8130394669
##  [946]  1.5654789264  0.1420295105  0.5625098979  0.8624975411  2.2701223452
##  [951]  0.8620674224  1.6925193698  0.2256976311  1.6778925564  1.5179200503
##  [956]  1.9943945044  1.7021384409  1.1679305188  1.7360520561  0.7627370057
##  [961]  1.4839927033  1.1139719964  0.2611728152  0.5885515074  1.0222272161
##  [966]  1.6879645740  0.4250448449  1.7231717001  1.4756185625  0.8018443761
##  [971]  2.3107427150  0.9388265402  1.2724765657  2.4798434194  1.5400692840
##  [976]  0.3153688634  0.8506724943  1.7864843646  0.8061387226  0.7562691040
##  [981]  0.8120901716  2.0071419158  1.5922150977  1.9489067568  0.5307408241
##  [986]  1.9989093389  1.6828813249  0.3512720776  0.8080182709  0.0007622836
##  [991]  1.0407853397  1.4868907927  0.5340842195  1.1091961987  0.5392737380
##  [996]  0.3775286960  0.9379575296  0.8517512635  0.7997966838  1.3631068220
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
```

```r
fit2
```

```
##      mean         sd    
##   0.5005957   0.4118744 
##  (0.1302461) (0.0920958)
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
## [1]  0.2499614 -0.1008305 -0.2665211  0.5634764 -0.1459139  1.2104035
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
## [1] 0.0138
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9096797
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
##     original        bias    std. error
## t1*      4.5 -0.0005005005   0.9191951
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 2 3 5 6 7 8 9 
## 1 1 1 2 1 1 3
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
## [1] -0.0025
```

```r
se.boot
```

```
## [1] 0.8923228
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

