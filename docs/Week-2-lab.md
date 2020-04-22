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
## 0 1 2 3 4 5 6 9 
## 1 2 1 1 1 1 1 2
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
## [1] 0.0334
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
## [1] 2.743373
```

```r
UL.boot
```

```
## [1] 6.323427
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
##    [1] 4.2 3.4 5.3 5.9 4.6 5.9 4.1 4.8 3.3 4.3 5.8 3.9 3.4 4.9 3.5 5.1 3.8 3.6
##   [19] 4.0 5.6 3.6 4.7 3.4 5.5 3.4 3.7 4.8 3.1 6.3 3.5 3.5 5.4 3.4 6.0 5.3 3.2
##   [37] 5.2 3.6 2.9 4.7 5.2 4.9 3.8 3.6 3.9 4.8 3.7 3.7 3.7 3.8 4.2 5.1 3.0 4.1
##   [55] 4.4 5.3 3.2 3.5 6.4 3.6 4.0 3.3 3.5 4.0 4.1 5.3 4.9 4.7 4.0 4.3 4.0 4.4
##   [73] 4.9 4.3 5.1 4.3 4.5 4.6 4.8 5.2 6.0 4.7 3.3 5.3 5.4 4.1 3.0 5.2 2.9 6.9
##   [91] 4.8 5.1 4.9 4.5 5.0 4.5 6.2 4.3 4.4 4.8 4.0 4.2 3.2 4.9 4.6 3.5 3.1 5.4
##  [109] 4.1 5.1 4.9 3.5 3.2 5.6 4.4 5.0 6.3 3.8 4.5 5.4 4.6 6.3 6.3 3.7 5.4 4.3
##  [127] 4.9 5.5 4.7 4.9 4.9 4.4 4.4 3.8 4.1 4.7 4.6 4.9 5.4 4.1 3.1 5.0 3.8 4.6
##  [145] 2.2 5.3 5.3 6.2 3.8 4.4 4.5 4.7 3.5 5.8 5.9 3.9 5.7 5.5 3.2 3.6 2.4 6.3
##  [163] 5.5 5.3 4.7 4.5 4.0 5.2 3.7 5.2 3.6 4.8 4.5 5.3 4.8 2.9 2.5 5.3 3.7 3.6
##  [181] 4.8 4.3 3.6 4.7 3.3 5.5 5.1 4.8 4.0 3.9 4.9 4.8 4.1 4.2 2.1 3.5 3.4 3.4
##  [199] 5.3 3.3 4.6 4.9 4.6 3.5 4.4 5.5 4.8 5.3 4.5 3.7 4.9 4.7 3.2 4.0 3.7 6.2
##  [217] 4.9 5.9 5.8 2.4 5.1 4.8 3.9 5.6 6.6 4.4 4.3 4.3 5.2 5.1 4.4 3.9 4.6 3.5
##  [235] 3.6 6.3 3.8 4.4 5.0 4.2 4.3 3.5 4.7 3.8 5.3 4.8 5.0 3.2 4.2 5.2 4.4 4.4
##  [253] 3.9 2.8 2.7 4.5 4.5 4.8 6.1 4.3 5.2 5.4 5.6 3.7 6.4 5.6 6.0 3.6 5.0 3.1
##  [271] 3.5 5.0 5.5 4.7 2.4 5.3 5.2 4.8 5.7 4.5 3.7 4.8 3.7 3.6 4.9 5.4 4.3 3.4
##  [289] 4.9 3.1 4.8 4.1 4.6 4.5 3.9 3.2 4.2 3.7 5.1 4.0 3.5 5.0 4.4 5.5 3.5 5.9
##  [307] 5.6 3.6 5.1 4.9 4.3 4.2 5.6 3.1 4.9 4.9 2.6 2.9 4.7 3.8 4.9 5.3 4.3 2.8
##  [325] 4.4 4.0 4.4 4.0 4.3 4.6 3.6 3.3 4.9 4.0 4.1 3.2 5.1 5.1 3.0 4.1 4.5 6.7
##  [343] 2.8 4.7 2.7 4.4 2.3 4.7 4.5 4.4 5.4 4.5 3.8 5.4 1.6 3.6 4.5 3.0 4.6 4.4
##  [361] 4.7 5.5 5.8 4.4 3.5 3.6 4.5 6.0 4.8 4.7 6.3 4.6 4.7 4.7 3.4 5.2 3.9 4.8
##  [379] 3.8 4.1 3.2 4.5 3.2 4.6 5.0 4.3 5.3 3.3 5.1 4.5 5.1 6.0 2.4 5.2 5.2 5.0
##  [397] 2.9 4.6 3.8 4.6 3.6 5.5 4.2 3.1 4.5 5.1 3.5 4.6 3.2 4.3 6.4 5.0 4.0 6.1
##  [415] 6.3 5.3 5.8 5.6 4.6 5.9 3.4 4.1 3.7 4.7 4.2 2.8 4.7 3.9 3.8 4.6 4.3 5.0
##  [433] 4.9 5.8 3.3 5.1 4.3 4.9 2.9 4.3 5.3 5.4 4.2 4.4 4.7 4.3 5.3 4.3 4.8 4.3
##  [451] 4.0 5.5 2.5 5.9 4.1 6.0 5.0 4.2 3.3 5.6 4.3 4.7 4.2 4.1 3.6 4.7 5.6 4.6
##  [469] 4.7 5.4 5.0 3.8 3.9 5.7 3.1 4.5 5.1 5.2 2.8 5.3 3.9 5.2 5.9 3.8 3.7 5.2
##  [487] 4.9 4.6 4.6 4.6 5.0 3.1 4.4 3.6 5.4 5.1 5.7 3.4 5.5 4.3 3.1 3.5 4.1 4.6
##  [505] 3.5 4.0 3.8 4.1 4.7 4.9 2.6 5.1 4.4 5.0 3.8 4.7 5.3 2.2 5.3 3.9 4.7 3.8
##  [523] 4.2 4.8 4.4 4.4 3.9 5.2 4.9 4.5 4.7 3.9 5.0 3.9 4.6 3.5 4.5 4.6 4.3 4.3
##  [541] 2.0 5.6 3.6 4.2 2.1 3.7 6.1 4.6 5.2 3.9 4.4 3.5 3.9 4.6 5.8 5.7 4.6 5.2
##  [559] 4.9 6.0 2.9 3.4 5.4 5.7 4.5 5.9 4.6 4.2 3.7 4.1 5.9 6.1 4.6 4.9 4.5 5.2
##  [577] 3.2 4.6 4.4 4.1 3.2 4.2 4.0 4.8 5.4 4.2 4.4 5.7 4.4 4.9 4.4 5.1 5.2 3.9
##  [595] 3.6 4.1 4.2 3.6 3.7 4.8 3.9 4.1 3.8 4.0 2.8 4.0 3.2 5.5 5.8 5.7 4.5 5.1
##  [613] 4.4 4.9 4.1 4.7 5.0 4.4 5.5 5.1 5.2 4.8 5.7 4.5 5.4 3.1 4.4 6.7 5.4 4.0
##  [631] 4.2 4.5 5.7 4.0 5.4 3.2 4.4 5.2 5.7 5.8 4.1 5.1 4.3 2.7 4.4 4.4 4.0 3.9
##  [649] 4.8 5.1 2.5 6.3 3.4 3.9 4.2 3.6 4.6 6.6 4.1 3.5 4.9 4.9 4.2 4.0 5.2 3.9
##  [667] 3.4 3.2 4.1 3.1 3.7 4.5 3.6 4.8 3.8 5.2 5.5 5.5 3.8 3.3 6.6 4.0 4.5 4.0
##  [685] 4.4 5.8 4.8 4.4 5.1 3.8 6.0 5.9 5.1 5.2 4.9 4.9 5.4 5.3 4.8 4.4 4.9 3.8
##  [703] 4.1 5.6 4.8 5.9 4.1 5.7 4.0 4.2 5.5 5.5 4.1 6.0 3.6 4.6 3.4 4.8 2.8 5.8
##  [721] 5.3 5.6 7.0 5.1 5.5 2.5 5.6 5.1 4.0 5.3 3.7 3.9 4.5 4.4 3.8 6.1 4.0 4.2
##  [739] 3.9 4.7 4.2 5.0 4.3 4.6 4.7 4.5 5.4 3.8 3.8 4.1 4.5 4.8 5.0 4.9 4.0 4.9
##  [757] 5.5 3.5 4.0 4.4 5.3 5.7 5.2 3.5 3.9 5.0 4.0 4.8 5.5 4.2 5.0 4.0 5.7 3.5
##  [775] 4.5 3.6 4.8 5.2 5.4 5.3 4.2 4.6 4.3 4.9 5.1 5.9 5.1 4.2 5.1 4.0 3.8 3.6
##  [793] 6.0 3.9 5.7 4.9 4.7 4.4 4.3 3.8 3.9 4.4 4.3 5.2 3.5 4.2 3.3 6.1 4.0 4.6
##  [811] 4.6 4.6 5.9 4.7 4.4 5.1 5.5 4.0 2.9 4.8 3.6 3.9 5.7 5.8 3.2 3.8 3.8 3.1
##  [829] 4.3 3.4 5.1 5.0 4.3 5.2 4.5 5.3 4.5 6.2 5.2 5.5 3.3 5.1 6.5 4.0 4.4 4.0
##  [847] 3.3 6.3 4.9 5.4 6.1 3.4 5.6 5.4 4.7 3.9 4.0 4.3 4.9 4.0 5.5 3.8 3.8 4.5
##  [865] 4.1 4.0 5.3 4.9 5.3 5.8 5.4 3.7 3.7 4.2 5.5 4.0 4.4 4.6 5.1 4.7 3.5 2.8
##  [883] 4.6 6.2 5.5 3.9 3.5 4.1 5.9 3.7 3.0 4.5 5.3 3.7 2.7 5.0 5.6 5.0 4.2 4.3
##  [901] 4.6 5.3 5.4 4.6 4.9 4.6 4.1 4.8 5.2 5.5 4.5 3.3 5.4 3.9 4.4 4.4 5.3 5.8
##  [919] 4.8 5.1 5.5 5.4 4.9 5.0 3.8 4.9 3.8 3.8 3.6 6.2 4.3 4.2 4.6 2.6 4.6 5.4
##  [937] 5.9 6.1 3.9 5.0 4.2 3.2 4.1 5.6 5.1 4.7 4.1 4.5 4.9 4.2 4.2 3.2 5.7 5.3
##  [955] 4.2 4.8 5.0 4.6 6.4 5.3 3.5 3.5 2.9 3.3 3.6 4.6 5.0 3.6 5.6 5.0 6.4 5.3
##  [973] 4.9 5.5 3.1 4.6 3.8 4.0 4.6 3.7 3.7 4.2 4.9 4.8 4.6 4.5 4.3 6.2 5.9 4.8
##  [991] 4.0 4.9 3.4 3.5 4.8 2.3 5.7 4.1 4.8 3.8
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
##    [1] 4.2 3.0 5.6 4.7 3.8 4.5 3.4 5.6 4.8 3.4 3.3 5.2 3.8 5.7 3.9 5.5 4.1 4.9
##   [19] 4.3 4.7 4.4 2.8 4.7 3.7 4.2 3.8 3.4 2.5 4.1 3.8 4.0 3.2 2.6 2.4 4.0 4.7
##   [37] 5.2 4.5 5.7 4.8 5.9 3.3 3.7 6.0 2.7 4.3 5.2 4.0 4.3 3.4 4.2 3.9 4.8 5.0
##   [55] 4.2 3.9 4.3 3.8 3.7 2.9 4.1 4.5 4.5 5.1 3.8 4.8 5.2 4.9 4.5 2.4 3.5 6.2
##   [73] 4.0 3.0 4.2 4.2 4.4 3.9 4.4 4.1 4.0 4.4 4.1 3.1 3.9 3.6 5.3 5.2 4.9 3.0
##   [91] 4.2 3.7 4.8 3.9 4.8 6.1 4.5 3.6 5.4 4.4 3.1 4.1 5.5 2.8 5.1 2.7 4.8 4.8
##  [109] 5.3 5.1 2.4 3.8 4.9 5.4 4.3 5.2 4.1 5.4 3.1 3.2 4.0 4.9 4.1 5.3 5.2 3.5
##  [127] 3.8 4.7 5.6 5.6 6.3 3.5 5.0 4.7 3.9 2.8 5.6 5.6 3.6 3.5 5.5 5.5 3.4 5.2
##  [145] 4.1 5.4 4.1 5.9 4.9 4.8 5.5 3.9 3.3 4.7 4.0 3.2 5.7 4.3 4.1 3.8 4.2 3.3
##  [163] 5.4 4.4 3.4 4.3 5.4 3.3 3.9 4.9 4.9 5.1 4.2 3.6 3.6 3.6 5.8 4.0 3.5 4.7
##  [181] 5.2 4.4 4.0 4.3 6.2 2.9 5.2 5.9 5.3 4.3 4.6 5.8 4.9 4.7 5.1 6.0 4.7 5.2
##  [199] 5.8 4.4 5.2 4.6 5.5 5.5 3.9 7.0 4.3 4.3 4.0 3.5 3.0 3.9 6.0 4.5 4.2 3.3
##  [217] 6.0 4.6 3.7 5.1 4.5 4.2 4.5 4.0 5.1 4.3 3.7 4.6 4.2 4.7 4.5 3.8 4.8 5.0
##  [235] 4.5 2.8 6.1 4.2 4.1 3.3 5.7 3.3 3.0 4.1 5.2 5.6 5.4 3.9 5.4 3.8 4.6 3.9
##  [253] 4.1 5.1 5.2 5.4 3.0 5.5 6.6 3.1 4.5 3.5 4.3 2.7 5.3 4.9 5.0 4.0 5.4 4.8
##  [271] 3.7 3.9 4.7 4.7 4.9 3.8 4.9 4.1 3.2 5.0 4.7 3.5 3.7 5.3 5.2 4.2 4.6 4.0
##  [289] 4.7 5.2 3.9 4.4 4.5 5.3 5.1 5.1 3.8 4.9 6.3 5.1 4.5 3.6 3.5 5.5 4.5 4.3
##  [307] 5.4 3.5 5.7 3.8 5.4 4.0 3.3 4.8 5.1 2.7 3.6 5.3 4.7 4.2 5.5 3.6 3.6 3.5
##  [325] 4.1 5.0 3.7 5.3 5.4 3.7 6.8 5.8 5.5 4.2 3.2 2.6 5.0 4.8 5.2 4.0 4.0 2.5
##  [343] 4.3 5.4 3.3 4.6 4.9 4.1 4.1 3.2 5.1 4.8 5.6 4.5 3.9 5.0 4.5 4.7 4.7 4.9
##  [361] 5.3 4.2 4.9 5.3 3.8 4.2 5.1 4.8 4.5 4.5 4.7 5.4 3.6 5.8 5.4 4.9 2.7 3.5
##  [379] 4.0 3.7 4.3 4.0 3.9 4.3 3.0 5.8 4.0 3.4 4.1 4.3 3.9 4.3 3.2 3.4 5.3 4.6
##  [397] 4.9 3.2 3.5 3.1 3.7 5.1 4.1 5.3 5.8 3.3 5.3 4.2 2.9 3.2 3.5 5.2 4.2 6.3
##  [415] 4.5 3.6 4.2 3.5 2.8 4.2 3.1 3.4 4.8 4.3 2.4 4.1 4.6 3.2 4.9 5.3 5.1 4.4
##  [433] 5.2 6.2 4.1 3.4 4.9 5.4 5.3 4.5 4.2 4.8 5.3 4.2 5.2 3.3 4.4 3.4 5.0 3.4
##  [451] 4.9 4.0 4.6 4.0 4.7 4.6 5.5 3.9 2.5 2.8 4.4 4.8 4.4 3.1 3.8 3.7 3.5 4.2
##  [469] 3.9 3.7 5.0 2.9 5.8 4.3 5.1 5.2 4.0 5.5 3.9 4.3 4.8 5.2 3.8 5.0 5.0 3.1
##  [487] 4.6 3.1 4.7 2.1 3.8 4.2 4.9 3.2 4.7 3.8 3.1 5.0 4.1 4.3 3.7 4.7 3.7 4.4
##  [505] 6.1 4.6 3.3 5.2 4.8 3.8 4.1 4.5 4.5 4.8 4.1 3.0 6.3 2.8 3.8 4.5 5.5 5.4
##  [523] 3.7 4.1 3.8 6.2 3.9 3.9 3.4 3.2 5.0 4.7 5.9 6.7 4.2 4.7 3.2 3.7 3.7 5.1
##  [541] 5.2 6.0 5.6 4.8 5.6 3.3 2.4 4.5 5.3 3.8 4.2 3.7 4.6 3.4 5.7 6.0 6.0 7.1
##  [559] 4.0 6.4 5.3 4.3 4.7 5.3 6.0 3.5 2.1 5.4 3.3 4.6 4.0 3.9 3.5 3.5 3.7 4.2
##  [577] 3.9 5.6 5.2 4.4 5.9 5.2 4.6 4.3 5.1 5.0 4.6 4.1 4.1 4.6 5.7 4.7 3.6 4.8
##  [595] 3.4 4.3 4.1 4.4 3.8 4.5 3.2 6.5 5.0 3.5 4.7 4.6 5.7 4.2 5.5 5.2 4.7 6.3
##  [613] 6.1 5.4 4.5 5.3 4.1 3.8 4.2 3.0 4.2 5.4 4.7 4.4 4.9 4.2 4.1 4.1 3.9 4.4
##  [631] 6.0 5.3 3.9 5.0 5.3 4.5 5.9 5.5 4.3 6.4 4.8 2.3 5.1 3.8 4.5 4.6 3.5 3.6
##  [649] 4.0 4.0 5.1 3.1 4.2 5.3 3.1 4.1 4.2 4.7 3.7 4.2 6.7 4.1 5.5 5.2 3.8 3.2
##  [667] 5.2 5.8 3.6 5.1 3.8 5.4 4.6 3.9 4.4 4.0 4.3 4.3 4.7 3.9 3.6 5.2 4.4 3.7
##  [685] 3.8 4.6 3.5 4.7 5.6 5.6 6.4 4.8 4.9 5.4 5.6 4.0 4.1 3.6 5.4 3.3 4.7 4.2
##  [703] 4.0 5.3 5.4 3.6 4.5 5.3 4.3 3.8 5.8 4.9 3.7 4.3 2.9 5.0 5.4 5.2 4.0 2.6
##  [721] 5.1 6.1 4.3 7.0 4.9 5.2 4.0 4.5 3.5 3.6 4.7 3.6 5.3 4.9 3.7 5.7 4.6 4.3
##  [739] 5.5 5.7 5.5 4.9 6.0 3.5 5.5 4.5 4.5 5.4 5.8 4.5 5.5 3.9 4.1 5.4 4.2 4.4
##  [757] 4.6 3.1 6.4 3.0 3.9 4.3 3.5 3.3 5.5 5.7 5.6 5.7 5.0 3.1 5.2 4.5 5.8 4.3
##  [775] 5.9 3.3 4.5 4.2 6.2 5.8 4.3 5.2 3.8 3.8 4.1 4.0 5.0 4.4 4.6 3.9 4.2 5.3
##  [793] 3.6 4.7 4.8 4.3 5.5 5.6 5.6 4.8 2.8 4.2 6.1 3.6 4.2 5.3 3.5 4.4 3.5 4.9
##  [811] 3.6 4.1 3.5 4.2 5.8 5.1 3.9 4.4 4.5 3.6 5.5 5.3 5.6 4.6 4.1 4.8 5.7 5.6
##  [829] 5.2 4.4 5.7 4.9 6.7 4.9 4.9 4.2 4.5 5.9 4.7 5.3 6.0 4.2 4.2 4.2 4.1 4.7
##  [847] 5.6 5.8 5.1 3.9 3.3 5.7 6.1 3.9 4.6 4.8 4.9 5.0 5.5 5.9 5.1 5.6 3.3 3.4
##  [865] 2.3 3.4 5.7 5.5 5.0 3.8 3.3 4.9 6.4 4.9 5.9 4.6 5.1 6.4 3.8 4.2 4.1 5.1
##  [883] 4.9 6.4 5.0 3.3 5.2 3.6 3.5 4.4 4.8 3.6 4.3 5.9 4.0 4.4 2.5 5.8 4.7 4.6
##  [901] 3.9 3.8 4.3 6.0 3.5 4.1 3.3 4.7 5.4 3.7 4.0 4.7 3.9 6.9 4.0 4.7 5.2 2.9
##  [919] 5.2 4.3 5.3 4.8 4.7 4.7 5.3 6.1 4.4 4.2 4.4 5.1 4.7 4.6 3.4 4.6 5.3 4.9
##  [937] 4.9 4.4 4.1 5.4 5.1 4.4 3.7 3.1 3.6 5.0 3.6 5.4 3.3 4.6 4.4 4.0 4.7 4.7
##  [955] 4.3 3.6 5.4 4.4 4.6 3.8 4.0 3.5 5.5 5.2 4.9 6.0 5.0 4.2 4.1 4.9 3.8 3.3
##  [973] 4.3 4.4 4.9 3.4 4.9 4.1 6.0 3.7 4.9 4.1 4.4 5.3 3.0 3.0 3.5 6.5 4.0 4.5
##  [991] 3.0 4.3 4.6 4.3 3.4 5.0 3.8 4.0 4.5 5.7
## 
## $func.thetastar
## [1] -0.0214
## 
## $jack.boot.val
##  [1]  0.43387978  0.32522255  0.30028409  0.11477987  0.04852071 -0.06000000
##  [7] -0.16381766 -0.30769231 -0.38173913 -0.49858357
## 
## $jack.boot.se
## [1] 0.9077487
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
##    [1] 6.5 2.6 4.0 1.7 7.4 4.3 2.5 4.9 5.9 6.3 4.1 2.7 4.8 3.8 4.2 4.8 4.7 4.4
##   [19] 2.9 4.6 5.8 3.9 5.7 3.8 4.1 3.8 4.1 2.8 3.0 5.3 3.5 3.5 5.3 3.9 4.6 5.1
##   [37] 4.7 6.1 3.1 4.5 4.7 4.2 3.9 4.3 4.4 4.4 3.1 4.8 4.5 3.6 3.7 4.2 3.1 3.5
##   [55] 5.8 4.4 3.7 3.6 5.5 3.8 4.1 6.0 5.0 5.0 3.3 3.6 4.9 5.1 3.8 5.7 4.0 5.4
##   [73] 5.0 6.6 4.2 3.8 5.5 4.6 5.0 5.9 5.4 2.9 5.9 2.6 4.4 5.8 5.3 4.7 4.7 5.3
##   [91] 2.8 4.0 4.7 3.7 4.6 5.4 3.5 5.3 4.4 5.9 4.0 4.7 4.1 4.3 5.0 5.4 5.5 3.2
##  [109] 4.3 3.2 5.0 5.7 4.0 5.6 3.3 5.4 3.6 5.8 5.2 2.9 5.2 3.6 6.4 4.0 4.3 3.9
##  [127] 6.2 6.1 4.6 4.0 4.1 4.4 4.4 3.1 4.4 6.0 5.1 5.1 6.0 3.2 4.4 4.8 5.4 4.1
##  [145] 4.7 5.8 5.6 3.7 4.4 2.9 5.4 2.4 4.8 5.5 4.3 3.9 5.8 3.1 6.4 4.6 3.8 3.1
##  [163] 4.9 5.8 3.8 3.9 4.8 3.7 3.0 3.3 3.6 4.3 4.2 3.5 4.4 3.7 4.7 4.2 5.3 4.2
##  [181] 5.2 4.9 5.1 3.8 4.3 3.1 4.0 4.4 4.2 4.3 4.7 5.3 4.1 5.0 4.3 3.1 4.8 4.7
##  [199] 5.3 4.3 4.3 4.4 4.7 4.6 3.4 4.3 4.3 4.9 6.2 4.2 5.2 3.8 2.2 2.5 4.7 5.1
##  [217] 4.2 4.7 5.2 4.4 5.1 6.8 4.8 4.1 6.0 5.1 5.8 2.4 5.1 3.8 5.8 3.5 4.9 5.3
##  [235] 5.5 5.2 4.5 3.7 4.1 6.0 5.5 3.7 4.7 4.2 6.0 4.6 4.5 3.2 3.7 3.7 2.9 3.2
##  [253] 4.6 4.7 5.3 5.0 5.1 2.9 3.3 5.6 4.5 3.8 4.1 6.2 5.6 3.9 3.5 3.9 6.0 4.2
##  [271] 5.0 4.5 5.1 4.3 6.0 4.0 4.7 5.9 4.3 3.2 4.3 4.8 4.8 3.4 4.2 5.3 6.2 3.7
##  [289] 5.4 6.2 5.4 4.5 5.3 2.3 4.7 3.8 5.7 3.7 5.1 3.1 5.4 4.1 5.2 4.1 4.6 5.0
##  [307] 4.7 6.1 4.2 4.9 3.5 4.3 6.0 1.7 5.3 5.0 5.6 3.9 4.2 5.4 4.4 4.5 5.6 4.8
##  [325] 3.7 5.4 4.3 3.6 4.1 4.0 4.7 3.5 4.1 5.0 4.7 4.6 4.6 5.9 4.8 3.2 4.0 4.9
##  [343] 5.3 5.7 3.3 4.6 5.4 5.3 5.4 3.6 5.1 5.4 2.9 4.5 4.5 3.9 6.0 4.5 4.6 5.2
##  [361] 4.2 3.7 5.0 4.0 5.2 3.0 3.1 5.0 5.6 4.4 3.5 5.3 3.4 5.7 6.5 6.1 4.4 6.2
##  [379] 4.2 4.4 4.4 3.6 4.0 4.6 4.6 4.6 3.0 4.2 2.8 6.2 5.6 3.5 4.9 3.3 4.2 4.1
##  [397] 5.6 6.1 5.5 4.3 3.8 3.3 5.1 3.5 4.8 2.3 4.9 5.0 7.9 4.1 5.0 3.5 5.1 3.5
##  [415] 4.4 4.7 3.2 4.5 5.2 4.6 4.7 3.9 4.5 4.7 3.4 4.9 4.7 4.2 4.7 5.3 4.8 5.3
##  [433] 3.8 4.9 6.3 3.4 5.6 4.0 5.0 4.8 2.7 5.6 3.2 4.9 6.7 5.3 6.7 4.8 6.1 4.7
##  [451] 5.7 4.6 4.4 4.7 4.3 4.7 4.6 4.9 5.3 4.4 3.3 3.8 4.1 4.1 4.6 4.1 3.0 5.0
##  [469] 5.0 5.7 4.7 4.3 3.6 4.6 4.1 5.8 4.3 4.3 4.8 4.0 3.7 5.5 4.5 3.8 4.5 5.0
##  [487] 4.4 5.5 6.0 3.0 4.2 5.3 3.1 5.1 5.2 5.1 5.8 2.2 4.0 5.7 4.6 3.4 6.3 4.3
##  [505] 6.1 4.0 4.6 4.9 5.0 4.0 4.6 4.1 3.4 4.6 2.8 4.7 6.1 5.6 4.5 4.5 5.4 4.5
##  [523] 5.7 7.4 3.9 4.5 5.1 5.0 4.5 4.4 3.4 3.6 5.7 5.4 5.7 5.0 5.3 2.8 2.4 4.3
##  [541] 3.9 5.8 5.1 5.1 3.4 4.4 4.8 4.1 4.6 3.9 4.9 4.6 2.9 5.8 3.3 6.1 6.1 3.5
##  [559] 4.8 5.9 5.1 4.2 2.9 5.1 4.5 3.1 3.8 4.2 5.9 2.7 4.9 4.6 4.9 4.0 4.9 5.1
##  [577] 5.5 5.1 4.6 3.3 3.9 5.0 4.7 4.6 4.6 6.7 4.7 4.3 5.3 5.1 4.0 4.1 3.4 4.5
##  [595] 4.2 4.2 4.3 3.0 4.2 3.9 5.6 5.2 4.6 5.6 5.4 5.3 4.2 4.7 5.3 3.6 5.4 4.0
##  [613] 6.2 4.9 4.9 5.8 4.6 5.2 5.0 4.0 5.5 5.0 5.7 3.8 4.2 3.9 4.3 4.8 4.8 4.4
##  [631] 4.2 5.1 6.2 5.8 5.3 4.9 3.9 5.0 4.4 4.7 4.2 4.4 3.2 4.4 5.5 5.3 5.0 5.1
##  [649] 3.6 4.1 5.9 5.2 5.4 5.5 3.6 3.8 5.4 4.3 4.3 5.3 5.7 5.1 4.0 5.4 4.5 3.8
##  [667] 3.5 5.2 5.4 2.9 3.5 4.6 4.6 4.7 4.1 6.7 5.2 5.7 3.5 4.5 4.4 4.9 2.4 2.7
##  [685] 4.5 3.9 4.7 4.3 5.8 4.5 3.9 3.9 3.7 5.1 4.6 5.8 3.9 3.3 4.4 2.9 3.3 4.3
##  [703] 5.6 4.0 2.8 4.0 5.4 3.4 4.8 2.3 3.5 4.6 3.4 4.9 5.9 3.6 4.3 6.0 4.1 3.9
##  [721] 5.8 4.4 4.5 5.3 4.7 3.8 4.3 4.8 5.7 3.5 4.6 3.6 4.5 4.3 2.3 4.2 4.9 3.5
##  [739] 6.5 4.4 5.0 4.0 3.1 4.8 5.5 4.4 5.6 4.8 5.0 3.5 4.5 7.1 4.4 5.0 4.8 2.5
##  [757] 3.7 5.2 4.8 5.5 2.3 3.7 3.5 3.7 4.3 3.3 3.3 4.7 4.1 3.8 4.5 4.9 6.6 5.1
##  [775] 4.3 3.8 4.7 3.8 3.9 5.0 3.0 3.6 5.7 3.6 4.9 4.4 6.2 5.0 3.5 5.0 5.2 6.1
##  [793] 5.5 4.6 3.6 3.9 4.4 2.3 4.6 4.8 4.9 4.8 5.7 4.7 4.4 2.1 5.6 4.9 4.4 4.0
##  [811] 5.2 3.4 5.2 4.5 4.3 4.0 4.4 4.4 6.0 3.7 4.8 4.6 4.9 4.5 4.3 4.5 5.1 4.7
##  [829] 5.3 3.7 5.5 5.5 3.6 3.8 5.8 5.4 4.6 4.5 4.5 3.8 5.9 4.9 4.8 3.8 3.0 3.6
##  [847] 3.6 5.1 5.0 4.5 5.3 4.1 4.7 4.8 4.5 4.6 4.4 5.3 5.0 5.1 4.8 5.0 4.7 4.8
##  [865] 3.9 5.6 4.0 4.2 5.5 4.9 3.8 6.8 4.1 2.6 5.6 6.1 5.4 5.3 4.7 3.7 3.2 3.0
##  [883] 4.7 5.0 3.0 4.4 5.5 4.8 4.2 5.5 4.2 4.8 3.5 4.6 4.1 3.5 5.4 3.7 4.4 6.0
##  [901] 4.3 5.2 3.6 4.0 4.3 5.4 5.7 4.1 4.5 6.1 3.8 5.4 5.7 5.6 3.6 4.8 3.6 4.3
##  [919] 2.0 5.4 5.1 5.6 5.3 5.0 3.7 4.8 5.7 3.6 5.2 4.0 5.4 4.2 5.1 3.3 3.3 6.0
##  [937] 4.4 5.1 5.3 3.5 5.0 4.7 4.3 4.9 4.4 5.0 5.6 3.8 3.1 3.8 5.1 5.5 5.4 3.7
##  [955] 2.6 3.1 3.3 5.5 4.7 2.3 4.9 3.0 3.9 4.9 6.9 5.5 3.2 3.6 4.4 3.9 5.6 6.2
##  [973] 5.0 5.3 5.4 6.4 3.5 4.4 5.0 4.7 4.4 5.3 3.5 3.8 4.6 3.6 3.6 3.5 5.5 5.2
##  [991] 4.0 3.1 4.9 7.1 4.9 4.1 4.2 4.4 5.6 5.9
## 
## $func.thetastar
## 72% 
## 5.1 
## 
## $jack.boot.val
##  [1] 5.5 5.5 5.4 5.3 5.2 5.0 5.0 4.8 4.7 4.4
## 
## $jack.boot.se
## [1] 1.046136
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
## [1] 0.8657736
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
##   3.275939   4.754453 
##  (1.396885) (2.190915)
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
## [1]  0.61300866  1.24835844 -0.21520693  0.06749493  1.14364131  1.11661225
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
##    [1]  0.933294817  0.605001932  0.931564924  0.787489220  1.140187212
##    [6] -0.204941673  0.558960035  0.367925306  2.101690204  1.899189187
##   [11] -0.846451994  1.029519566  0.043591209 -0.390979102  0.707125906
##   [16]  0.060290937  0.894887507  2.097348276  0.847005220  1.054851117
##   [21]  0.674158770 -1.315064250  0.721919586  0.659796921  1.189790762
##   [26]  0.993059406 -1.009924895 -1.125530373  1.114586763  1.419866849
##   [31] -0.378170722 -0.462557035  0.879439662 -0.681403449 -0.518200107
##   [36] -0.255392151  1.563804482  0.296806755 -1.087032779  2.015278709
##   [41]  2.140758228  1.635964865  1.798803225  1.326033975  0.965862711
##   [46]  1.309199207  1.204450012  1.183090622  1.536539983  0.551542471
##   [51]  0.391720760  1.223337607  0.761166865  0.870503430  0.306662140
##   [56]  1.157764330  0.496486123  0.492184283 -0.632606928  0.685285126
##   [61] -0.292699966  0.477553761  1.334308370 -0.646373408  0.383044776
##   [66]  0.269785078  0.189710701  0.837875096  1.446158637  1.102132954
##   [71] -0.782968036 -0.900712104 -0.481776197  0.175589131  1.881213434
##   [76]  0.546593629  1.372872161  0.429608821  1.185045250  0.700093419
##   [81] -0.187451692 -1.634400376 -0.433893190  1.354238747  0.781832567
##   [86]  0.365949846  0.866517829  0.233973745 -1.868306118  0.898605943
##   [91]  1.411680946  1.447447567  2.084065843  0.586051693  1.610231721
##   [96]  0.024342648 -0.319830692  1.000275667  1.502836051 -1.536693112
##  [101]  1.229751397  1.602094755  0.235990971 -0.184685477  0.473065184
##  [106] -0.405690025  0.415608029  1.243536152  0.550477283  1.357625878
##  [111]  0.132833884 -0.376952518  0.923908016 -0.390632313  0.324066583
##  [116]  0.381992104  1.670360815  0.208931572  0.866364704  0.683624339
##  [121]  0.797976856  1.161773224  1.231350631 -1.302503188  0.939329186
##  [126] -0.648343330  0.579068695  0.187451041  2.436063344  0.740503481
##  [131] -0.856725713  2.041312740  0.889721348 -0.350626195  1.792137648
##  [136]  1.444804762  0.591665994 -0.418832874  1.805714180  0.963725532
##  [141] -0.170551227  0.697579684  1.619075953  1.076620379  0.573720573
##  [146] -0.929762358  1.092482189  0.853428968  2.152417500  2.261450428
##  [151] -0.052842475  0.738885292  0.663014957  1.023520842  0.873026589
##  [156]  0.742285110  0.913876520  0.007600748 -0.832836396  0.093186490
##  [161] -1.040765830 -1.127216388 -1.034816749 -0.546391193  0.922987224
##  [166]  0.507099498  1.305575782  1.553922373  0.648560905  1.272532395
##  [171] -0.481776197  1.321501517 -0.703090178  0.745875300  1.039273828
##  [176]  0.764241101  1.332613443  0.382421102  0.413680790 -0.176565612
##  [181] -0.066349361 -0.270825476  0.655321559  0.756454653 -0.010138941
##  [186]  0.488089499  0.908532450  0.816009983 -0.190075625  0.953533938
##  [191] -2.058601447 -0.685938367  0.121924633 -0.633434151  0.867284328
##  [196]  1.065997949  0.844166863  0.685509552 -0.458589691 -0.023117565
##  [201]  0.587399458 -0.550105718 -0.815443247 -2.292123009  0.956206346
##  [206]  0.767398177  0.201170079  0.709461253  1.897067156  0.773447973
##  [211]  0.519307185  0.476900674  0.670173668 -0.637691537 -0.287374066
##  [216] -1.085296966  0.987084810  1.439507021  1.816907239  1.098121237
##  [221]  0.900553746  1.481107376  0.817584003  0.272113955 -1.288487295
##  [226]  0.279649177  1.201639074 -1.117235825  0.592198168 -0.676006632
##  [231] -0.889997711  0.239170198 -0.471215750  0.418202945  0.753406781
##  [236] -0.624320385 -0.746676620 -0.912178680 -0.464381237  0.795674232
##  [241] -1.205445664  0.355043937  0.777469799  1.367637194  1.014567693
##  [246] -0.400290734  0.347662908  0.736587745  0.467464007 -0.868867131
##  [251] -0.479725716  1.052925724  1.522739986 -0.617839455 -0.106437916
##  [256]  1.340306060 -0.457786603  1.100596485  0.782454660  0.876952204
##  [261] -0.130356625  1.210754296  1.868947827  1.136054546 -0.372392207
##  [266] -0.469506447 -0.579035676  0.679360519  1.014722556  0.983385153
##  [271]  1.069139093  1.203595205  1.181212919 -0.797583995  0.649669610
##  [276]  1.673325590  1.098365576  0.592687066 -0.437262328  0.989287709
##  [281]  0.735484497  1.272441637  0.773226303  0.225809136  0.754341569
##  [286] -1.348875622 -1.200349164  0.683330675  0.852003219 -0.827165590
##  [291]  0.347493212 -0.940160711 -0.192480382  1.055577035 -0.731826758
##  [296]  1.367052358  0.865773553  0.355668540  0.816400952  1.035668384
##  [301]  0.169217961 -1.029269203  0.484950263  1.568508064  0.609892247
##  [306]  1.179559543 -0.393265189  0.568101912  1.183561670  1.834722523
##  [311]  1.121888584  0.866661186 -0.005302090 -0.101950539 -0.756801959
##  [316]  1.363892300 -0.738546819  0.088303668  0.973028543 -0.673191355
##  [321]  0.644588550  0.819074356  0.615363669  0.909238227  0.427110519
##  [326]  1.720959182  0.543938574  0.818942789 -0.752988487  0.061102420
##  [331]  1.040487534  0.918347996 -0.170276493 -0.605985960 -1.138750526
##  [336]  1.128590265  2.025063872 -0.518593876  1.101534201  1.246984637
##  [341] -0.961003224 -0.625779614  0.606838980  1.011617733  1.343655577
##  [346]  0.955211843 -0.265783102 -0.707403813  0.244078850  1.512731128
##  [351]  0.354515953  1.282461730  1.148709126  0.851171625  0.685848702
##  [356]  0.801673894  1.222451333 -0.011228368  0.085042597 -0.671110968
##  [361] -1.030000766 -0.848512816  0.453573360  1.844457105 -1.113724445
##  [366]  1.400226360  0.502646608  1.172746710  0.520812728  1.007463459
##  [371] -0.042775460 -0.655892204 -0.811893410  1.157213524 -0.887677796
##  [376]  0.983602347 -1.270469779  0.369714198  0.763518434  1.147345396
##  [381]  1.083095161  0.769891807  1.692622960 -1.019508164 -0.581194522
##  [386]  1.007463459 -0.139846807 -1.176961301 -0.477626670  0.777188430
##  [391]  1.296758881  0.103579127 -0.814895785  0.690294899 -0.605678879
##  [396] -0.543361593 -0.912259150 -0.272572850  1.270075424  0.873152944
##  [401] -1.197580376  1.368946008  1.234608810 -0.251928243  0.455901258
##  [406]  0.813172064  0.356584089  1.114840629  0.630329788  0.857889940
##  [411]  1.198768453  1.080812409  1.862136802 -1.002359275 -0.951844905
##  [416]  0.240042271 -0.037907298  1.130499812  1.284575200  0.783538709
##  [421]  0.534555905  0.748548928  0.768745347  0.554141430 -1.301009140
##  [426] -0.890381775 -0.576685130  1.207522485  0.576766636 -0.381936305
##  [431]  1.078694496  0.861828434  0.391500118  0.582904968 -0.502312660
##  [436]  1.176612980 -0.879134490  0.081374057  1.123789317  0.183368475
##  [441]  0.615948193 -0.473808523  0.913756462  0.691871355  0.529671542
##  [446]  0.609978374  1.622762552  1.694586232  0.625428576  0.655542454
##  [451] -1.009924895  0.990723093  0.418381787 -1.320104626  0.616451411
##  [456] -0.790157534  1.494134926 -0.821352605 -1.492541398  1.150562386
##  [461]  0.958898495  1.071743141 -0.048886583 -0.522645304 -2.028836947
##  [466]  0.926154231  0.758377522 -0.408300673  0.320663239 -0.379166662
##  [471] -0.476426546 -0.726539734 -0.391570303 -0.487189318  0.015990327
##  [476]  1.048298155 -0.712799368  0.952094550  0.373382912 -0.936075903
##  [481]  0.463168204  0.141013234  0.549012879 -0.045295688  1.013737240
##  [486]  0.999637701  0.975209104  0.517158134 -1.590908493 -0.730900688
##  [491]  1.341333005 -0.500542933  0.630278828  1.388050790  0.204077487
##  [496]  0.327302648  0.616998946  0.449671668 -1.099883995  0.715529270
##  [501]  1.362594386  0.185209338 -0.536616999  0.843215517 -1.300893026
##  [506]  0.747663135  0.332195507 -0.239989477  0.721058171  0.803865701
##  [511] -1.178496187  2.055412463  1.969873972  1.106974102  1.689926596
##  [516] -0.447872351 -1.001485117  1.234400301  1.072481379 -0.765895402
##  [521]  1.272532395  0.821901463  0.594751381 -0.088502486  1.941371072
##  [526] -0.888013931  1.615257448  1.087236115  0.961133626  1.168122000
##  [531]  1.680721498  1.186133224  1.165481872  0.632295093  0.244926131
##  [536]  0.568567128 -0.883506635  0.910431045 -1.028584142  0.611864921
##  [541] -0.004666852  0.563342029  0.178226471 -0.529451486  0.430724627
##  [546]  0.669092150  0.160080655 -1.030749808  1.048913588  0.943203812
##  [551]  1.128615601  0.914129033 -0.924742768  0.247610851  0.376164911
##  [556]  1.016069285 -0.619895626  0.901717512 -0.169859648  1.636976712
##  [561]  1.384086982  0.881683114  0.456862188  0.568019113  1.097244469
##  [566] -0.417801551  0.318787643  0.561399467 -0.693729606  0.807476175
##  [571]  0.335044172  1.170362157  1.167838226  0.720066470  0.772171337
##  [576]  0.804228527  1.697037724  1.936555018  1.384820504  1.862136802
##  [581]  0.414538437  0.886614587 -0.002294285  0.866305068  0.188781245
##  [586] -0.606136605  0.833085117 -0.767865103 -0.134320034  0.622342113
##  [591]  0.897417121  0.999417375 -1.291363205  1.962056868  1.140346059
##  [596]  1.366468244  0.724871379  0.406283352  1.171126732 -0.401428580
##  [601]  0.242694414  1.762091485  0.404987932 -0.888529082 -1.021040315
##  [606]  1.236839520  0.747585708 -0.703692935 -0.403252001  1.459454504
##  [611]  0.877457740 -1.143036667  0.832426987  1.249508959 -0.125500536
##  [616] -0.271074984  0.368413220 -0.918573763  1.621117900  0.955899846
##  [621]  1.183465576  0.550686251  0.874629805 -0.128958203 -0.657050632
##  [626]  0.734839863 -0.158214458  0.765527066 -0.328136459  0.036870611
##  [631]  1.230406079  0.898537669  1.137627691  0.852392417 -0.937463289
##  [636]  0.593261986  0.413515965 -0.311346673  0.417822387  1.470461569
##  [641] -0.213449151  1.344130063  1.669374785  1.017162916  0.788754461
##  [646]  1.887664643  0.688169723  0.514555841  0.026719644  1.173592246
##  [651]  0.265942731  0.878524875 -0.995317437  1.203959960  0.569573251
##  [656]  0.929285995  1.128549970  1.014869325  1.079272190  0.910737331
##  [661]  0.860384325 -0.620437126 -0.329911791  1.585563055 -1.652677082
##  [666] -0.704035132  0.796067388 -0.497464182  0.642869237  0.616034943
##  [671]  0.701908350  0.934423944  0.940131108  1.648736778 -0.105745514
##  [676] -0.898459076 -1.307288601  0.802605147 -0.527707059  0.514905549
##  [681]  0.215345855  0.615948193  0.536571248  0.853715605 -0.064580693
##  [686] -0.620184479 -0.190250138  0.587585283 -0.454494988  0.384706936
##  [691]  1.169422535 -0.543998451  1.216279970 -1.594037525  0.695217192
##  [696]  1.302769757 -0.049757680  1.877329081  0.739077466  0.522706178
##  [701]  1.301020520 -1.916286814 -0.376701455  1.172893600  1.510437081
##  [706]  1.179135879  1.293723597  0.887052144  0.769719440  0.399992895
##  [711]  1.235321787  0.650641990  0.985350389 -1.077411960 -0.863581344
##  [716]  0.237526308 -0.176605021  0.505923133  0.841187757  0.178507214
##  [721]  0.105222642  0.496976639 -0.356409466  0.773191720  0.930907685
##  [726]  0.699627236  1.512731128  1.116430529  0.978793816  1.444463386
##  [731]  0.149415714 -0.300285487  1.756081552  0.450856183  0.958157485
##  [736]  0.892766673 -0.395780512  0.770861436 -0.576835087  1.198076961
##  [741]  0.689244072 -1.076341608  1.982255958  0.923908016 -0.158054166
##  [746] -0.143649104 -0.015764086 -0.542286934 -0.356971811 -0.157333471
##  [751] -0.529619504  0.716068883  1.270607988  0.163572344  1.267971334
##  [756]  0.544579792  1.125420118 -1.160162625  1.098365576 -1.858881274
##  [761] -0.040193556  1.170365040 -0.507253270  2.057424270  0.336487214
##  [766] -0.218025795 -1.379234628  1.303431847 -1.628826949 -0.007097882
##  [771]  0.912108934  1.432850715  0.116471024  0.871346673  0.866231809
##  [776] -0.142624336  0.648158012  1.268075835  0.726141620 -0.499221026
##  [781] -0.799432854  0.120441855  0.809851104  0.484794141  2.211299551
##  [786]  1.102888486  0.805169347  0.849595303 -1.647865224  1.650300125
##  [791]  1.515164935  0.867962137  1.124203254 -1.612716789 -0.199565565
##  [796]  1.312050176  0.876930783 -0.276868221  1.126432594  1.856050487
##  [801] -0.488568042  1.953259332  0.866310396 -0.762885350 -0.411094487
##  [806]  0.906749365  0.816400952  0.326625499  1.266717289  0.977815743
##  [811] -0.755660401 -0.488891764  1.052213156  0.890640090 -0.248267626
##  [816]  1.071982202  0.656224464 -0.876776325  1.135328127 -0.721454095
##  [821]  1.336918996  0.178297948  1.230721804 -0.721454095 -0.794282347
##  [826]  0.716167352  0.616160018  1.401241939  1.321018251  1.464546317
##  [831]  0.875186354  0.650055331  0.825476725 -0.212112766 -1.190659203
##  [836]  0.643580742  1.340306060  1.434424561 -0.693958080  2.061677477
##  [841]  0.626420188  0.643794436  0.394795292 -0.716259359  0.471953588
##  [846]  0.695634522 -1.131285012  0.106493883  1.188533492  0.708706852
##  [851]  0.548524427  1.660784823  1.216882535  0.127939821  0.518918990
##  [856] -0.429950889 -0.001927454  0.561973132  0.034539602  1.419042983
##  [861] -1.165689584  0.190364981 -1.334928135 -0.557312948 -1.009808382
##  [866]  0.787862691  0.106479453 -0.461663643  1.013044821  0.674232358
##  [871] -0.361619865  0.525561763  0.851303612 -1.424763044  1.292373252
##  [876]  0.725503668  0.762693410  0.408081560 -0.957905050  0.545823431
##  [881] -0.574106548  0.702066416  0.645621777  0.618711352 -0.782161529
##  [886] -1.039559674  0.230125024  1.555643997 -0.505922823  0.200402784
##  [891]  1.002694542  1.099251270  0.880598107  1.140287868  1.936758581
##  [896]  1.180439466 -0.942250059  0.557072318  0.533501991  0.967542292
##  [901]  0.135692534  0.668221066 -0.928120602  1.200002577  0.292885482
##  [906]  0.859841991 -0.004353478  0.910737331 -0.605985960 -0.721732302
##  [911]  1.510801987  1.078381757  0.287367940  0.853505692  1.848037167
##  [916]  0.866429188 -0.549982560  0.739023785 -1.742500257 -0.602566378
##  [921] -0.326846322  0.777350189  0.369898950 -0.353809486  1.326532644
##  [926]  1.256455722 -0.229591206  0.999614192 -0.208911008  0.586254598
##  [931]  0.279778974 -1.719710765  0.207881313 -0.754134187  1.957467175
##  [936]  0.648721365  0.478450322  0.906741512  0.906929034  0.255944904
##  [941] -0.443166622 -1.175507137  0.502896044  1.444396291  0.951536106
##  [946]  1.321197293  1.138259679  1.075431744  0.405985629 -0.512089445
##  [951]  1.958478225  0.159159093 -0.479725716 -0.785458643  0.385869976
##  [956]  1.766904024 -0.820421778 -0.568874408 -0.336324743  1.204450012
##  [961]  1.574338108 -0.613388004  1.356001831 -0.664219683  0.631175355
##  [966]  0.828944548 -0.569794168  0.021226945 -0.068268630 -0.211173946
##  [971]  2.401856261  0.660004831  0.747877476 -0.830999128  1.219193092
##  [976] -1.462544010  1.253223107  1.034714612  1.544387001  0.820218568
##  [981]  1.862136802  0.733777832  0.650687381 -0.156880076 -0.967808532
##  [986]  2.113203867  1.122514876  0.343684703 -0.872932362  1.284072428
##  [991] -0.473922408  0.903713853  0.051011547 -0.691827592  0.933032225
##  [996] -0.599323263 -0.183326504  0.977094415  0.625951446  0.885356946
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
##   0.68902708   0.35480154 
##  (0.11219810) (0.07933422)
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
## [1] -0.2775080 -0.2603954  0.2107337 -0.3176093 -0.3513032 -0.4545804
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
## [1] -0.0057
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9148412
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
## t1*      4.5 -0.02662663   0.9122259
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 1 2 3 4 5 6 8 
## 2 1 1 1 1 1 2 1
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
## [1] -0.0061
```

```r
se.boot
```

```
## [1] 0.9274875
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

