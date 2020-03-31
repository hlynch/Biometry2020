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
## 0 1 2 3 4 5 8 
## 1 1 2 1 1 3 1
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
## [1] -0.0245
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
## [1] 2.757777
```

```r
UL.boot
```

```
## [1] 6.193223
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.7   6.1
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
##    [1] 3.9 5.1 3.2 2.0 3.4 5.9 3.8 4.6 4.5 4.6 4.9 4.6 5.4 4.8 4.9 5.0 5.6 4.6
##   [19] 3.3 3.1 5.1 3.8 4.3 5.2 5.1 3.5 4.3 4.3 4.3 3.9 5.9 3.3 4.1 4.9 4.6 6.3
##   [37] 3.8 5.0 5.8 4.9 4.5 4.5 5.4 3.4 4.1 5.7 6.1 2.8 4.1 5.2 3.3 5.4 3.4 5.4
##   [55] 5.2 4.8 5.1 5.9 4.4 5.7 3.2 5.3 3.8 3.3 3.6 5.5 3.8 3.0 3.1 5.1 3.2 5.7
##   [73] 4.5 4.5 4.2 4.4 3.4 4.2 6.4 5.3 4.5 4.7 3.1 4.0 4.4 4.7 5.3 5.9 4.7 5.9
##   [91] 5.2 5.7 4.9 4.7 4.3 5.4 3.5 4.5 4.8 5.7 4.1 4.4 4.8 4.2 3.9 3.9 3.3 4.9
##  [109] 4.6 5.2 4.7 4.1 3.6 4.5 4.2 6.1 3.9 3.3 5.1 6.2 4.2 4.1 3.9 6.0 5.3 3.8
##  [127] 6.3 4.5 4.6 6.0 3.0 4.3 3.1 5.5 5.2 3.1 3.7 4.2 5.4 4.3 4.9 4.3 4.7 3.5
##  [145] 4.1 3.6 5.6 5.1 3.4 4.2 5.0 4.9 3.7 5.0 4.5 3.4 6.3 3.9 4.3 3.7 3.8 3.7
##  [163] 3.2 4.1 3.6 4.6 4.2 5.6 5.6 6.1 4.7 4.5 5.4 4.3 5.1 4.8 3.9 3.8 6.1 3.0
##  [181] 5.8 5.0 5.0 4.4 4.1 4.6 4.0 4.3 4.0 5.3 5.4 3.4 4.3 4.5 5.3 4.1 6.0 4.4
##  [199] 5.1 5.8 4.8 6.0 4.5 3.5 5.3 3.7 4.8 4.1 3.7 5.0 2.5 3.7 4.9 4.2 3.8 5.0
##  [217] 3.0 3.9 4.6 4.6 4.6 5.3 5.9 4.4 3.9 3.6 4.6 6.5 3.9 4.5 5.1 5.4 3.3 4.0
##  [235] 5.5 4.7 5.0 4.2 4.6 4.4 5.5 3.6 6.1 5.9 3.2 4.0 4.5 4.9 4.7 3.3 3.3 4.4
##  [253] 2.5 4.2 6.5 4.4 4.8 3.0 4.0 5.2 4.2 6.3 3.8 3.5 5.2 5.4 4.6 2.8 3.8 6.2
##  [271] 3.1 5.5 4.3 5.1 5.5 5.3 5.5 3.7 5.6 4.5 4.0 3.7 3.5 5.1 4.8 4.5 4.6 4.2
##  [289] 4.8 3.5 4.0 4.3 3.4 3.2 2.3 2.4 5.0 3.6 3.5 3.9 5.3 4.8 5.0 4.8 3.7 6.5
##  [307] 5.1 4.4 4.9 5.4 4.7 4.9 5.1 5.0 4.7 4.3 2.7 5.4 4.7 4.8 4.0 3.9 6.1 5.7
##  [325] 5.0 4.4 6.1 4.6 3.5 4.4 5.3 5.9 3.6 4.4 3.8 5.4 3.3 4.8 4.2 5.8 3.5 4.3
##  [343] 6.2 4.9 4.6 3.9 5.7 4.3 4.1 4.3 3.9 4.5 4.4 3.5 4.3 5.9 2.5 4.8 5.4 5.1
##  [361] 3.7 5.4 5.4 3.8 5.2 4.6 3.1 5.1 4.5 5.1 4.3 4.7 3.0 2.4 3.2 4.5 3.5 5.0
##  [379] 5.2 4.0 4.6 5.3 4.4 3.5 5.0 4.5 4.5 4.4 4.6 3.2 2.9 5.0 5.1 4.8 4.5 5.0
##  [397] 4.7 5.8 3.1 4.7 5.8 4.7 5.3 3.4 5.4 4.6 3.9 2.3 4.2 3.1 4.0 4.7 3.2 4.8
##  [415] 3.4 4.4 4.5 5.3 6.1 6.2 4.5 5.0 5.1 4.6 3.4 5.5 5.9 6.6 4.6 4.9 4.6 4.3
##  [433] 5.5 4.2 3.8 5.0 4.4 4.8 4.9 3.9 5.2 4.7 5.0 4.6 3.7 4.3 4.0 5.5 3.3 4.6
##  [451] 5.3 4.6 5.1 4.2 5.7 5.4 3.2 3.3 3.2 5.2 5.1 5.0 3.3 5.6 5.8 3.7 4.7 3.1
##  [469] 5.4 5.3 4.5 5.1 4.5 4.6 5.4 4.7 4.3 4.4 4.0 5.1 3.9 4.6 5.8 4.3 4.3 2.7
##  [487] 2.3 4.4 4.3 4.9 6.1 3.3 5.6 5.7 6.1 3.8 4.1 4.4 5.1 4.7 3.3 4.0 4.4 3.5
##  [505] 4.6 4.6 2.6 4.5 3.4 4.0 4.3 3.7 4.9 4.1 4.3 5.3 5.0 3.8 4.0 5.9 2.6 5.1
##  [523] 5.4 4.3 3.8 4.8 4.6 3.6 2.9 5.4 4.7 4.8 2.7 5.0 2.5 5.0 4.0 6.3 4.2 4.2
##  [541] 3.2 4.6 4.5 3.7 4.3 6.2 5.1 4.8 2.7 4.9 5.7 3.7 5.7 4.5 5.0 3.9 4.0 6.4
##  [559] 5.7 3.4 4.9 4.2 4.5 4.7 4.9 4.6 4.3 4.9 5.5 4.7 4.7 5.4 2.5 4.4 2.6 3.0
##  [577] 3.7 4.4 4.1 5.7 4.7 5.2 3.2 4.0 3.9 5.7 4.5 4.0 5.2 3.4 4.4 4.9 2.9 5.0
##  [595] 4.5 4.0 4.6 4.4 4.3 4.2 3.4 6.2 4.5 5.3 3.4 4.3 4.5 4.1 5.6 3.8 4.4 3.9
##  [613] 2.3 4.1 5.3 4.9 3.4 3.4 3.6 3.7 6.5 2.7 5.2 3.9 6.8 2.7 4.0 3.0 5.3 4.7
##  [631] 4.3 2.7 4.8 3.9 4.5 4.1 5.0 3.5 5.5 5.2 5.3 4.3 2.5 5.3 3.9 5.2 4.4 6.6
##  [649] 5.1 3.9 3.8 4.8 5.3 5.8 5.7 4.5 5.6 4.6 4.7 4.7 3.8 5.0 5.8 5.2 4.5 3.4
##  [667] 5.6 3.4 4.9 3.5 4.6 3.6 4.6 4.1 5.0 3.8 4.7 3.2 3.7 5.9 4.1 5.8 4.6 5.0
##  [685] 4.4 4.8 2.8 5.9 3.3 4.4 4.2 3.6 3.4 4.2 5.2 3.5 4.6 4.0 3.6 3.0 5.3 3.3
##  [703] 3.7 5.1 2.9 3.7 5.3 4.5 6.4 6.2 5.6 4.4 5.9 5.4 5.2 5.1 5.3 2.8 4.8 3.8
##  [721] 5.1 5.0 5.3 5.5 4.5 3.4 4.9 5.0 5.5 2.7 3.9 5.2 4.5 5.3 6.2 4.9 5.1 3.3
##  [739] 4.9 5.3 5.7 5.1 3.9 3.3 5.3 4.3 3.3 5.7 5.1 4.4 4.5 3.9 5.7 4.6 3.8 3.9
##  [757] 4.8 4.3 5.6 4.2 2.9 5.9 6.1 4.8 3.2 3.2 5.1 4.8 4.9 4.8 4.5 5.1 3.4 4.1
##  [775] 4.3 5.8 4.1 4.6 3.0 3.8 4.7 5.8 4.4 4.6 4.1 4.2 5.1 5.2 3.5 4.2 4.4 3.4
##  [793] 4.3 4.5 5.4 3.7 4.8 6.1 3.9 4.2 4.8 3.1 5.3 4.1 4.5 5.1 4.8 5.4 4.9 2.7
##  [811] 3.6 3.8 5.1 2.8 4.8 4.5 5.8 4.4 4.6 6.6 4.7 4.2 4.6 6.6 5.9 4.5 5.6 5.3
##  [829] 4.3 4.1 3.8 2.8 2.1 5.3 2.9 4.7 3.1 4.9 4.3 4.6 5.4 4.7 3.9 2.1 4.6 4.8
##  [847] 4.4 6.5 3.4 4.4 4.1 3.4 6.1 5.2 4.1 5.5 3.8 4.6 3.9 4.8 4.5 5.4 4.9 6.8
##  [865] 5.8 2.6 4.2 3.5 3.0 5.9 5.5 4.0 4.4 5.6 5.2 5.0 3.9 5.3 3.2 5.9 4.0 4.8
##  [883] 3.9 5.2 5.6 2.8 5.0 3.3 4.3 4.9 6.1 4.0 5.8 4.7 3.9 4.9 5.0 2.8 3.3 4.2
##  [901] 5.0 5.3 4.7 4.3 6.0 2.3 3.5 4.0 5.4 6.1 2.5 4.9 4.2 4.1 5.5 4.3 4.6 6.6
##  [919] 4.2 4.1 2.5 4.7 3.7 3.6 4.3 5.1 3.3 5.0 3.9 4.0 3.7 6.4 2.7 5.8 4.1 5.4
##  [937] 5.7 3.0 4.4 4.6 4.3 4.5 5.0 4.4 4.3 3.9 5.5 3.5 3.8 5.5 5.0 5.2 5.2 6.3
##  [955] 5.9 2.9 4.8 4.6 5.0 3.3 4.0 4.1 5.5 4.0 4.2 4.8 4.5 5.0 4.6 3.6 5.2 4.6
##  [973] 3.1 4.3 6.0 4.6 4.8 5.4 4.2 3.8 6.0 5.8 3.7 4.2 5.2 4.8 4.3 3.4 3.3 4.5
##  [991] 3.2 3.8 4.6 3.4 3.4 4.4 4.8 4.8 4.4 4.4
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
##   2.7   6.2
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
##    [1] 5.6 3.6 4.2 4.2 4.6 6.0 4.5 5.3 5.6 6.4 4.2 4.5 4.4 4.1 6.2 5.1 5.2 5.0
##   [19] 4.0 6.0 4.8 5.0 5.1 3.3 6.1 4.7 4.7 4.0 5.8 4.4 5.6 4.5 4.2 4.2 3.8 4.4
##   [37] 4.6 4.6 5.4 4.3 4.0 5.0 3.8 3.8 5.1 6.2 5.2 6.0 5.2 5.0 3.4 3.3 3.7 2.7
##   [55] 5.7 6.1 5.6 5.2 5.6 4.3 4.1 4.1 3.8 4.8 4.7 3.3 4.4 3.7 6.0 4.4 3.1 3.5
##   [73] 3.5 4.4 3.7 4.5 4.0 3.9 5.0 3.8 6.0 4.9 6.1 6.0 4.6 3.7 2.9 3.9 4.5 4.2
##   [91] 6.9 4.8 3.7 6.6 4.4 4.2 5.2 4.5 5.2 4.2 4.6 5.8 2.5 4.9 3.9 4.3 3.6 5.2
##  [109] 5.3 4.4 5.5 5.5 5.4 4.6 3.6 3.6 3.1 4.0 5.9 3.7 4.6 5.0 5.0 1.4 4.9 6.2
##  [127] 5.4 5.1 4.4 5.2 4.5 6.0 4.4 5.3 4.9 3.6 3.9 3.6 3.9 4.2 4.6 4.6 2.8 3.4
##  [145] 6.3 2.9 2.6 4.9 3.6 5.0 3.5 3.9 4.1 3.5 3.2 4.6 3.9 4.2 4.9 4.9 4.7 3.5
##  [163] 3.9 2.2 4.4 4.0 4.4 3.6 5.7 4.1 4.5 4.2 5.2 4.0 4.4 3.7 5.6 4.6 3.5 4.8
##  [181] 3.9 5.4 3.9 3.4 4.6 4.8 2.8 4.1 3.5 4.6 4.7 5.0 4.4 4.6 4.2 5.2 4.7 4.5
##  [199] 4.5 3.2 4.4 5.1 5.5 5.1 4.3 4.4 3.3 3.8 5.6 3.4 5.0 4.2 6.0 5.3 5.2 5.2
##  [217] 5.4 6.8 3.7 4.0 2.5 4.3 6.0 5.6 4.3 3.0 3.5 4.3 4.2 5.4 5.1 3.4 3.8 4.4
##  [235] 4.4 5.4 6.3 5.5 4.5 2.7 5.0 4.8 5.2 4.9 4.3 5.0 4.0 5.1 4.3 4.3 5.4 3.5
##  [253] 3.9 5.2 4.2 3.9 5.1 4.4 3.3 5.2 6.3 5.4 4.0 5.1 4.5 3.6 3.6 4.9 3.6 5.5
##  [271] 3.2 5.8 2.1 5.8 6.2 4.1 4.5 4.6 3.8 5.4 3.5 5.5 4.8 5.0 2.6 5.3 4.4 5.2
##  [289] 6.1 5.2 4.6 4.1 4.2 4.7 3.6 5.3 5.2 4.3 3.6 4.9 4.6 6.7 6.6 4.2 6.4 4.7
##  [307] 3.5 4.9 4.2 4.8 3.5 4.0 4.1 3.8 5.0 4.9 3.6 4.8 5.1 4.8 4.8 4.9 3.3 4.4
##  [325] 5.5 4.2 4.9 4.1 5.0 3.4 4.0 4.0 5.1 4.7 3.6 3.9 5.5 3.9 4.6 4.1 3.7 6.5
##  [343] 3.2 1.6 6.2 4.8 4.7 4.8 4.9 4.2 3.3 5.4 3.9 3.0 4.3 4.3 3.3 5.2 5.8 3.3
##  [361] 2.9 5.8 3.6 5.0 3.4 3.6 2.9 4.6 4.8 4.2 5.1 4.5 3.4 3.3 4.5 3.0 4.6 3.8
##  [379] 3.5 3.7 3.0 5.4 4.7 3.3 4.1 3.3 2.9 4.7 4.2 5.6 4.9 6.1 5.2 4.0 4.8 4.8
##  [397] 4.6 5.6 4.4 4.2 3.1 4.5 4.0 5.3 3.1 6.0 6.4 4.4 5.3 5.1 4.1 4.8 4.7 4.5
##  [415] 4.8 3.5 4.5 4.1 3.7 4.4 4.8 4.9 4.9 4.3 5.3 4.3 4.2 3.1 4.8 5.3 3.4 5.1
##  [433] 4.7 4.4 4.0 4.6 3.8 4.7 4.8 4.6 5.1 4.7 4.5 4.1 4.4 3.7 4.0 5.1 3.7 4.7
##  [451] 3.9 4.9 3.8 5.0 5.2 4.7 6.1 4.2 4.4 4.0 5.3 2.6 4.2 4.7 4.1 4.9 5.1 4.4
##  [469] 4.7 4.1 3.3 3.5 4.1 6.7 5.2 4.4 3.0 3.4 3.7 4.6 6.4 5.4 3.4 3.5 5.0 3.3
##  [487] 5.1 4.9 4.8 5.2 3.5 3.3 4.2 4.2 4.5 4.4 4.1 4.4 6.2 5.0 6.2 5.0 3.9 5.2
##  [505] 6.1 4.9 6.2 5.4 4.8 5.1 5.2 3.9 5.0 4.5 4.3 4.4 4.9 4.1 6.1 3.7 6.3 5.6
##  [523] 4.3 3.8 3.6 5.0 4.3 3.2 5.5 4.0 5.6 5.2 4.3 3.4 5.8 3.6 4.0 3.9 5.1 3.8
##  [541] 3.6 5.5 2.8 4.8 4.1 6.6 4.6 4.8 6.4 3.6 5.4 4.0 4.7 3.7 5.4 5.8 4.8 3.7
##  [559] 4.9 4.5 4.6 2.9 2.4 6.1 3.1 3.8 4.2 3.0 4.2 3.9 4.8 4.3 4.1 3.3 5.1 5.0
##  [577] 4.3 6.2 5.8 3.6 6.7 3.9 4.3 5.0 4.1 5.6 4.4 3.1 4.0 5.0 3.1 4.7 4.1 4.8
##  [595] 4.8 5.6 6.6 4.0 3.8 4.3 4.3 5.6 4.1 4.3 3.4 4.8 4.7 6.1 3.5 6.0 4.0 4.5
##  [613] 6.2 4.2 2.3 4.0 5.7 3.5 6.7 5.4 6.3 5.0 3.6 5.4 3.8 4.8 5.9 4.0 5.1 3.5
##  [631] 5.8 7.5 3.5 3.8 4.9 4.6 3.2 2.9 3.7 4.6 3.6 4.9 5.5 6.3 4.7 4.2 2.8 4.7
##  [649] 3.0 3.2 5.3 5.1 4.0 3.9 4.0 3.7 4.2 4.6 5.1 4.9 5.1 2.3 4.6 6.0 5.0 4.7
##  [667] 4.7 4.7 4.6 5.5 3.6 4.2 4.1 2.9 4.1 4.9 4.0 4.0 3.2 3.7 4.4 4.5 4.9 4.3
##  [685] 5.1 5.5 6.6 4.1 4.9 5.5 5.1 4.5 5.9 5.1 4.0 5.0 2.6 4.4 5.5 5.3 4.1 5.4
##  [703] 5.0 4.0 3.8 3.2 5.1 3.7 5.7 6.5 5.6 3.1 4.5 4.5 4.9 3.7 4.8 3.9 4.6 5.9
##  [721] 3.3 4.0 4.6 3.6 5.8 4.8 4.4 5.1 4.4 4.4 5.6 2.7 4.9 4.4 3.8 4.3 3.3 3.6
##  [739] 5.5 5.6 5.7 4.5 4.8 2.5 4.5 3.8 3.9 4.7 5.1 4.4 4.7 3.6 4.9 4.0 4.2 3.5
##  [757] 4.6 5.2 5.0 4.3 6.6 4.4 4.6 3.9 5.7 3.2 4.6 3.6 4.0 4.0 3.5 2.8 5.3 4.6
##  [775] 5.6 4.7 4.4 3.5 4.3 4.7 4.0 3.6 4.1 4.0 4.2 3.6 3.7 5.5 3.4 3.5 4.7 4.8
##  [793] 5.4 4.4 4.0 5.2 5.4 6.0 4.4 5.5 3.7 3.3 2.5 5.0 5.2 3.2 6.4 4.0 5.4 3.5
##  [811] 4.7 6.0 5.3 3.3 2.9 3.6 6.6 4.8 5.0 4.5 5.2 5.0 4.3 3.9 3.5 3.7 6.1 3.6
##  [829] 2.6 3.2 5.0 5.5 6.1 3.9 4.6 3.8 3.3 4.6 3.8 3.5 5.1 4.9 5.8 3.9 5.0 4.6
##  [847] 4.5 5.6 3.9 4.8 5.0 6.3 5.6 4.0 4.2 7.2 3.5 3.4 4.7 4.0 3.7 2.2 2.5 6.0
##  [865] 2.3 4.8 4.4 4.8 4.0 4.9 5.3 3.7 3.9 3.7 5.3 5.3 5.7 5.2 4.5 4.8 5.1 3.9
##  [883] 6.4 3.8 4.4 4.9 4.3 2.2 5.7 4.5 3.7 4.6 5.9 4.5 3.0 5.2 3.6 3.7 3.3 3.9
##  [901] 4.7 5.7 3.8 4.2 5.5 4.6 3.3 4.6 3.7 4.6 4.3 5.2 4.5 3.9 4.6 6.0 4.5 4.1
##  [919] 5.3 4.6 2.4 4.9 4.3 5.0 4.8 5.7 4.5 4.8 4.6 5.3 6.1 5.1 5.7 4.4 5.0 3.5
##  [937] 3.1 5.1 4.1 4.6 5.2 4.4 5.1 4.2 4.1 3.2 4.0 6.0 4.3 4.3 4.6 3.3 3.3 3.7
##  [955] 5.2 4.3 3.6 4.2 4.7 5.7 5.5 5.4 4.6 4.7 4.3 4.4 5.5 3.0 4.1 4.5 3.9 3.9
##  [973] 5.7 5.4 5.3 1.8 3.4 3.3 4.1 3.0 4.9 5.3 4.7 3.8 6.9 5.4 4.1 5.3 3.2 4.3
##  [991] 4.5 5.8 6.2 4.8 2.5 4.8 5.0 3.8 3.0 3.9
## 
## $func.thetastar
## [1] -5e-04
## 
## $jack.boot.val
##  [1]  0.58238994  0.33241758  0.32015707  0.15486726  0.03578275 -0.09243697
##  [7] -0.15337079 -0.25153203 -0.42212389 -0.51636905
## 
## $jack.boot.se
## [1] 1.002914
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
##    [1] 4.6 4.7 3.8 5.2 4.4 3.0 4.2 4.7 6.1 3.7 5.4 6.0 4.6 3.0 3.4 3.9 6.3 5.8
##   [19] 5.1 4.3 3.9 5.6 5.6 3.1 4.2 3.8 5.5 3.7 3.4 6.5 6.1 4.2 5.4 5.7 3.7 5.0
##   [37] 4.4 5.7 4.9 4.0 2.5 4.0 5.0 4.8 5.4 4.5 2.8 3.7 4.9 4.1 4.7 5.6 3.9 4.7
##   [55] 4.2 5.5 6.2 4.6 5.2 2.7 6.0 2.8 4.4 4.7 4.2 5.2 5.7 2.2 4.0 4.5 3.2 5.9
##   [73] 5.1 4.9 2.8 4.8 3.4 4.3 4.7 4.3 5.0 3.6 3.7 4.0 5.0 4.9 4.3 4.1 5.7 5.9
##   [91] 4.5 3.8 5.6 3.9 3.5 5.4 5.0 4.3 4.3 5.3 4.0 4.3 5.3 3.4 6.4 2.9 5.0 4.5
##  [109] 4.3 5.7 4.8 3.1 4.3 4.8 5.2 5.1 4.4 4.7 4.2 6.7 4.9 5.2 5.4 3.6 5.2 4.6
##  [127] 3.8 2.6 3.1 3.3 5.4 3.4 5.0 3.7 6.1 3.4 5.0 4.4 6.0 5.5 3.8 4.0 4.0 4.1
##  [145] 3.7 2.8 3.8 3.1 4.9 5.1 3.7 3.2 5.6 6.0 5.0 4.1 4.7 4.8 4.6 6.3 2.7 6.5
##  [163] 2.8 4.2 4.1 5.7 3.5 2.9 3.1 3.5 4.3 3.8 4.3 5.6 5.6 6.0 4.8 4.5 3.5 4.0
##  [181] 4.1 5.1 5.4 3.2 5.6 4.5 4.5 4.9 5.0 5.1 3.8 4.2 3.9 5.1 4.7 3.4 4.9 3.6
##  [199] 5.0 6.7 4.0 4.8 4.2 2.9 5.0 4.5 3.3 5.4 3.8 4.7 4.6 6.0 5.2 5.7 5.6 4.4
##  [217] 5.4 5.1 4.4 4.5 5.1 5.4 4.5 2.5 5.7 3.0 5.0 2.7 4.9 4.4 5.9 2.2 5.6 3.0
##  [235] 4.0 4.7 5.5 5.2 4.4 3.6 4.2 5.0 4.4 4.4 4.3 4.2 6.0 4.6 6.2 3.2 4.2 5.2
##  [253] 3.5 5.6 3.4 6.4 5.9 3.9 3.0 3.7 3.4 4.9 5.1 5.5 2.7 5.0 6.3 2.9 5.1 5.9
##  [271] 3.3 3.8 3.3 3.9 4.9 4.1 4.5 3.2 6.4 5.5 5.0 4.6 5.8 2.7 5.5 4.4 3.7 4.1
##  [289] 4.5 3.3 4.8 4.3 5.2 4.4 3.9 4.8 5.2 4.5 4.8 4.3 3.7 5.3 4.1 4.8 3.7 5.2
##  [307] 4.8 5.7 4.2 3.9 5.5 4.6 4.5 3.6 3.8 4.2 3.8 4.5 4.4 4.4 7.2 4.4 5.7 4.8
##  [325] 3.0 3.1 5.8 3.7 3.9 4.8 4.7 4.1 3.8 3.7 5.4 3.4 4.6 5.2 4.8 3.9 3.2 5.6
##  [343] 6.8 4.9 4.1 4.6 5.1 4.0 5.4 4.6 4.0 5.8 3.7 5.2 3.7 4.8 4.1 5.7 3.1 5.0
##  [361] 5.9 4.9 2.9 4.0 3.6 2.3 5.2 4.8 3.8 5.9 5.5 4.1 2.7 2.2 3.8 6.1 4.3 4.7
##  [379] 3.5 5.9 4.6 4.8 5.3 3.2 5.7 4.9 5.9 4.2 4.5 4.7 3.9 3.0 3.8 5.3 5.3 4.2
##  [397] 5.6 4.0 4.5 4.7 5.2 4.7 4.5 6.7 3.3 3.2 3.0 3.1 3.8 5.4 4.8 4.0 4.9 4.5
##  [415] 5.4 3.6 3.7 4.4 5.3 5.5 5.0 5.3 4.2 5.0 3.6 5.7 4.7 3.0 5.0 3.5 3.3 6.0
##  [433] 4.3 5.9 5.6 6.8 5.4 5.1 5.6 4.6 5.0 5.8 4.9 7.2 4.1 4.1 4.7 3.8 4.5 3.9
##  [451] 5.4 4.0 5.1 3.6 4.6 3.9 3.9 6.1 2.8 4.2 3.2 3.9 4.9 5.4 3.4 5.7 3.9 3.6
##  [469] 5.0 4.6 5.1 3.8 4.3 4.6 5.0 3.9 4.7 6.1 3.6 5.2 3.9 3.2 4.0 3.6 4.9 4.3
##  [487] 4.8 4.0 6.0 3.8 3.1 4.1 6.7 3.8 3.3 4.6 5.3 3.0 3.5 5.3 2.7 3.9 5.0 5.1
##  [505] 6.0 2.5 4.1 5.0 5.9 5.3 3.1 5.3 4.1 5.8 4.1 5.0 5.2 4.7 4.8 4.8 4.4 3.9
##  [523] 3.9 4.6 3.8 5.4 3.5 2.7 3.1 5.4 5.5 5.5 4.8 4.9 3.9 4.4 4.5 5.6 4.2 3.1
##  [541] 6.7 4.1 4.9 5.0 3.7 5.0 3.3 5.6 2.1 4.5 2.5 4.1 4.0 4.1 4.7 3.6 3.1 5.5
##  [559] 3.6 3.1 4.4 4.4 4.7 3.9 4.6 4.5 3.4 4.6 6.0 3.9 3.6 6.0 3.7 3.1 4.3 3.3
##  [577] 5.0 4.8 4.4 5.5 4.8 5.2 4.4 6.4 4.6 3.7 4.2 5.9 3.3 6.5 3.8 3.8 3.7 5.8
##  [595] 4.2 3.3 4.5 3.8 3.6 4.3 3.7 5.0 6.8 6.3 3.2 4.2 3.3 3.7 4.7 3.7 3.8 5.6
##  [613] 3.3 4.9 4.4 3.9 4.6 5.0 3.9 3.9 4.1 4.3 4.2 3.5 4.5 4.6 4.7 3.7 4.6 5.5
##  [631] 3.3 4.1 5.0 4.6 3.5 4.4 6.3 4.7 5.8 3.7 3.8 4.0 3.4 4.5 4.8 5.4 5.4 4.7
##  [649] 5.1 5.6 4.5 3.4 6.1 4.1 4.6 4.0 4.7 3.2 6.2 4.2 4.4 3.3 5.5 6.3 3.1 3.8
##  [667] 4.7 2.7 4.4 4.6 4.8 5.9 2.7 4.4 4.2 3.7 2.5 4.9 4.6 3.2 5.8 4.2 4.7 3.1
##  [685] 5.3 3.3 3.4 5.5 3.8 3.2 3.8 6.3 4.6 4.6 4.4 4.5 4.4 5.2 4.2 5.3 4.1 4.4
##  [703] 3.0 3.9 3.6 4.9 4.8 4.8 3.6 4.6 4.8 4.6 4.7 5.7 4.6 5.8 6.4 3.1 5.0 3.8
##  [721] 3.0 4.4 5.0 4.8 3.0 2.8 3.3 4.5 3.5 4.1 5.1 5.9 3.5 4.2 5.9 5.8 4.1 4.8
##  [739] 4.4 5.0 4.6 5.1 2.5 4.1 4.8 5.0 5.1 4.3 4.8 3.6 2.6 4.7 5.0 3.8 4.5 4.1
##  [757] 5.2 5.2 3.5 4.1 2.9 4.9 5.9 4.1 4.3 4.5 3.5 2.6 3.7 4.4 5.2 5.1 4.7 3.2
##  [775] 4.5 5.1 5.1 4.0 4.1 4.6 3.0 3.9 5.3 4.8 4.6 4.6 5.6 5.0 3.7 3.3 3.8 5.0
##  [793] 5.1 6.1 5.0 4.7 4.2 4.1 4.7 4.0 5.0 6.2 5.0 3.7 4.5 4.7 6.2 5.7 4.5 5.4
##  [811] 3.7 3.0 4.2 5.2 4.3 3.8 5.3 4.9 2.5 5.0 5.4 5.5 3.9 4.4 5.0 5.5 4.2 4.9
##  [829] 4.5 4.0 4.1 5.8 4.2 5.9 5.3 5.0 4.4 3.9 4.6 5.2 4.4 5.6 4.2 5.6 3.8 3.7
##  [847] 5.3 6.2 2.6 3.7 5.1 4.9 3.5 5.3 5.3 3.3 4.5 4.0 5.0 4.6 5.0 5.1 3.4 4.9
##  [865] 2.8 4.6 3.8 3.6 4.1 6.4 3.9 4.8 4.4 4.0 5.8 3.9 4.5 5.2 4.6 4.8 4.2 5.2
##  [883] 5.5 3.6 2.4 5.4 3.9 3.8 4.6 5.6 4.8 4.8 5.5 6.7 4.9 3.8 3.4 4.4 5.0 4.2
##  [901] 4.9 5.0 5.7 4.4 3.5 4.1 4.5 4.2 3.9 5.0 4.6 5.9 3.7 4.1 4.2 4.8 3.7 5.5
##  [919] 5.4 4.6 5.6 4.7 4.6 4.1 3.9 3.7 6.0 4.7 5.2 4.0 4.8 4.4 5.5 3.8 4.0 3.5
##  [937] 6.4 5.3 3.2 3.6 4.6 4.6 5.3 5.3 1.6 3.6 4.6 4.4 5.2 4.8 5.2 2.7 5.3 3.8
##  [955] 5.7 5.3 4.5 5.8 4.3 5.8 5.0 4.0 5.5 4.9 3.9 5.6 3.8 3.6 4.0 4.4 6.7 3.6
##  [973] 4.1 3.1 4.5 6.1 5.2 4.8 5.0 3.0 5.6 3.2 4.6 4.3 4.8 3.9 4.4 5.9 4.3 4.9
##  [991] 4.4 6.2 5.1 5.3 3.9 3.0 3.0 5.6 3.6 3.1
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.600 5.452 5.400 5.200 5.100 5.000 4.904 4.700 4.600 4.500
## 
## $jack.boot.se
## [1] 1.066046
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
## [1] 0.1494392
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
##   3.014816   6.355150 
##  (1.280699) (2.937331)
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
## [1]  0.9382089 -0.1261748  0.6532820  0.3160619  0.4296186  0.6000152
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
##    [1] -0.1638575667 -0.2109406966 -0.2843366256  0.0448905749  0.3356060146
##    [6]  0.1266038140  0.4844384821 -0.0390488391  0.4977272741  0.1181945757
##   [11] -0.8330658624 -0.4417204454  0.3696219782  0.2079004034 -0.1617167540
##   [16]  0.2811462046 -0.2916568125  0.0052759509  0.4525442941 -0.3514584096
##   [21]  0.1750317328  0.5021475746  0.6424148147  0.5840772193 -0.3576439855
##   [26]  0.5786317349  0.2189712455  0.4419185917 -0.1995521644 -0.2061686902
##   [31]  0.3696685000  0.3077097935 -0.7515597024  0.7492038583  0.4821548569
##   [36]  0.1222889464  0.4551319397  0.2289445905  0.2667338314  0.3967289882
##   [41]  0.4060147375 -0.7772237265 -0.6017473272  0.1225911816  0.3621842182
##   [46]  0.0695699053 -0.0062650620 -0.6579006912 -0.4377583046 -0.0605795387
##   [51]  0.1166745803 -0.6479234617 -0.6245300973  0.0984321036  0.1793241595
##   [56] -0.1815961617  0.4720825602  0.7087128861 -0.4286588965  0.1520303046
##   [61] -0.6475544774  0.5259207984  0.2461192409 -0.0768438681 -0.2796542229
##   [66] -0.1792316942  0.3380284434  0.2064084201 -0.5225897453 -0.1938463770
##   [71] -0.2511715722  0.2604103673  0.3238415556  1.2352423336  0.4637696007
##   [76]  1.2703689897 -0.0421160371  0.0629193432  0.6764602817 -0.0141263924
##   [81]  0.3751723500 -0.0132554248  0.1896151448  0.1470198509  0.1479411767
##   [86] -0.4856572460  0.2485462583  0.1326675289  0.5227157611  0.0167229138
##   [91]  0.0748582217 -0.3179908038 -0.3619332633 -0.0665284181  0.2337935023
##   [96] -0.4688593586  0.1501284585 -0.1829305509  0.5988506066 -1.3893771527
##  [101] -0.3839935155 -0.2331553282 -0.4772817675 -0.6555369235 -0.1375380821
##  [106]  0.3236476319  0.4050954648  0.3307248243 -0.6227008831  0.2518428557
##  [111]  0.1760961106 -0.3226338691  0.5029267229 -0.6097348227 -0.1382519659
##  [116]  0.2175699181  0.4207803625  0.0738771480  0.2695345703 -0.2015728892
##  [121] -0.3149303874  0.5028162859 -0.0511443498  0.1039405977  1.2836212617
##  [126] -0.3913405110 -0.0018800545  0.0445031095 -0.8997507363  0.1428405212
##  [131] -0.1662160232 -0.0984636746  0.2781013051 -0.2505036840 -0.7854727202
##  [136]  0.4451102410 -0.5711825774  0.2464813643 -0.9625284210  0.2957548644
##  [141]  0.3139781329 -0.0706158904 -0.1254650265 -0.0886238567  0.3796760064
##  [146]  0.2989155350 -0.2313770569  0.4274944702 -0.0730076656  1.1085883439
##  [151]  0.6157553680  0.5264447091 -0.2988160190  0.5390212011  0.1140913915
##  [156] -0.1137947465 -0.7028977163  0.4993747786 -0.1721716799 -0.5217737018
##  [161]  0.2352249886 -0.3774914470 -0.7052620386  0.2029919344  0.7622529000
##  [166] -0.0401340325  0.2610473992  0.4447373965 -0.5197505631  0.1947445000
##  [171] -0.0300676746  0.2711325894 -0.1605769622  0.1239437791  0.0806561838
##  [176] -0.4660020411  0.5146560107 -0.7701613279  0.9278400968  0.6793893495
##  [181]  0.1675697156 -0.1665440755 -0.4169066376 -0.3240250700  0.0804262860
##  [186] -0.4278502009  0.4065376559  0.2643280334 -0.7242249286 -0.0919663955
##  [191] -0.5701632360 -0.0085979300  0.4892615113 -0.0471938140 -0.0948273510
##  [196]  0.0668741400 -0.1073043431 -0.4918629388 -0.3364366104 -0.1528692951
##  [201] -0.5203250518 -0.7047936528 -0.0686652239 -0.3121580271  0.4715510649
##  [206]  0.4171365296  0.3291708663  0.0143302600  0.2112630962 -0.2182045622
##  [211] -0.5971340220 -0.1266629189 -0.2626101465  0.1913804156 -0.6960913247
##  [216] -0.1988150303  0.4016432314 -0.5691851991  0.1039809192  0.6443746049
##  [221]  0.0627655876 -0.0098512989 -0.1301534317  0.1357837354  0.8783466983
##  [226]  0.7546255051  0.1903163755 -0.4849719801  0.1142191963  0.0668741400
##  [231] -0.0961690686 -0.3179074464  0.2988928309 -0.4341475039 -0.7108181336
##  [236]  0.0941505725  0.0635923267 -0.4900550603  0.2040911989  0.0364906677
##  [241]  0.3527853652  0.2525795601  0.6595835387  0.2503663983  0.4837948326
##  [246] -0.0953943728 -0.9010270268  0.3983332644  0.2056723243 -0.1865028428
##  [251] -0.1641837721  0.2868764596  0.3469391790  0.9044843718 -0.1534774790
##  [256]  0.0274649735 -0.6720921611 -0.3839613813  0.3590902088 -1.4760872502
##  [261] -0.2106072094  0.2198038771  0.4079657659 -0.0017988742  0.4246388189
##  [266]  0.3661836521 -0.0118332510 -0.2926099194  0.1266302083  0.0953139670
##  [271] -0.0208838976 -0.2182955602  0.7210879318  0.1432333273  0.4853655118
##  [276]  0.3495326045  0.3076929786  0.1375807116  0.7012022292  0.4275174584
##  [281] -0.5655047377  0.1308675891  0.2245117533  0.5665807874  0.0786317625
##  [286]  1.3545023645  0.0607355439 -0.2184737429  0.4676119721  0.3367429261
##  [291] -0.5603038334  0.1277785074  0.3062128232  0.0066073353 -0.1258267240
##  [296] -0.0803462235 -0.1071590323 -0.1713276797 -0.1062296132  0.2174690391
##  [301] -0.1260147931  0.1343016253 -0.3527639191 -0.3580544324 -0.7712014661
##  [306]  0.6307140398 -0.1843638728 -0.3549787544  0.5858441519 -0.5621186982
##  [311]  0.5621119759  0.0672646330  0.0036449793 -0.6209876988 -0.1556923952
##  [316] -0.3862376344 -0.5913098731 -0.2556050834 -0.1690608012  0.2270028567
##  [321] -0.3098340157  0.1011221369 -0.4027123343 -0.2368986027 -0.0453945492
##  [326] -0.2788844738 -1.4057359009 -0.1563194639  0.1823682534 -0.1655774042
##  [331]  0.3077097935 -0.2707470135 -0.5752095993 -0.0197932182  0.2524776374
##  [336]  0.1197824376  0.6425375367 -0.0454954811 -0.3368843745  0.2066477821
##  [341]  0.3983332644 -0.3687441950 -0.2112274930  0.5153853052 -0.0656510301
##  [346] -0.0398637878 -0.0390488391  0.4740294558 -0.3246277543 -0.7960149008
##  [351]  0.1948323517  0.5903051583  0.4123253054 -0.1126473436  0.0086974568
##  [356] -0.6640238444 -0.5184677888 -0.3810458679  1.0473307930  0.0561258630
##  [361] -0.0887756439  0.9264690256  0.5581522808 -0.0374199134 -0.2156156780
##  [366]  0.1195450053  0.0772856574 -0.5687860330 -0.5978936103  0.1556558456
##  [371]  0.4183105771  0.5118933050  0.2037343784 -0.2910653638  0.8004572996
##  [376]  0.5417298381  0.3588788422  0.3251712869 -0.4198601983  0.6761668722
##  [381]  0.4698905374  0.0724339711 -0.1145503010 -0.2309053530 -0.2640358610
##  [386]  0.3060542015  0.6761152052 -0.2785386927  0.4136054094 -0.2759053643
##  [391] -0.1668240168  0.2329632266  0.4253549121 -0.2945017729 -0.1239022658
##  [396] -0.1329020714  0.1001600629  0.2369994692  0.8508176327  0.2108663358
##  [401] -0.0524192085  0.3443444027  0.5214278883 -0.0929125209  0.6485423984
##  [406]  0.4844426825  0.2549006189  0.3416902457  0.6818235905 -0.2579663821
##  [411] -0.1989171988  0.8744180451 -0.0537191859  0.9551792514  0.2736670793
##  [416]  1.0911861736 -0.3318504867  0.2350514435 -1.0218239617  0.1060142078
##  [421] -0.4703768828 -0.3770422454  0.0378603385  0.4434222288 -0.0631702312
##  [426] -0.4300732891 -0.5398544429 -0.1396363523  0.5216055959 -0.2914446393
##  [431]  0.2454967285  0.4512475946  0.3609531663  0.1297831977  0.6490889267
##  [436]  0.1845305097  0.3589922033  0.4017675094  0.1110796805 -0.4392796825
##  [441] -0.0989074453  0.1302201920  0.5940318045 -0.4937483380  0.9454957411
##  [446]  0.0095028229  0.3539036630  0.8808809203  0.5665807874  0.2156318921
##  [451] -0.3375856367  0.7928186121 -0.2186281706 -0.0676049600 -0.4629423627
##  [456] -0.5752164305  0.5694617441  0.0242918328 -0.0127995349 -0.0106821802
##  [461]  0.6226633652 -0.1106570221 -0.4928184266 -0.0037370515  0.0066216386
##  [466]  0.4751837748  0.5313123496 -0.4430280872  0.9814163496  0.1530389130
##  [471] -0.0786402001  0.4700327741 -0.4126188390 -0.1979333445 -0.5160974619
##  [476] -0.0432490957 -0.5970547816  0.0502425982  0.2718852084 -0.6164615376
##  [481]  0.3207159163  0.2759751878 -0.1652262201  0.4034223896  0.7622064062
##  [486]  0.9429335645  0.0453888248 -0.1691357858  0.6175955028 -0.1271475073
##  [491] -0.0008911605 -0.6254435879  0.2386327856 -0.8253127477  0.4961521298
##  [496]  0.3022967076  0.9116627868  1.3089863003  0.3761835971  0.2786671065
##  [501] -0.0270575996  0.4917190898  0.0477218285 -0.8499722320 -0.0675301823
##  [506] -0.4674938079 -0.3481686836  0.4162819786  0.7024504924 -0.7596277083
##  [511]  0.3860187566  1.1420708653 -0.9682823892  0.1763366969  0.0904495271
##  [516] -0.2932735099  0.4654672030  0.5919308459 -0.2718787062  0.0664602164
##  [521] -0.0451351478  0.5112902414 -0.2414578274  0.2260031575 -0.0247997514
##  [526]  0.1314569745  0.0664602164  0.2796046514 -0.6218594605  0.4761751118
##  [531]  0.1168978985 -0.5988101904  0.1156077903 -0.6579554379 -0.1320100544
##  [536]  0.2125413263 -0.2289685373 -0.1028825286 -0.2038019205  0.2249878854
##  [541] -0.1588315524 -0.4216189455 -0.2998313789  0.9010179440  0.5119051434
##  [546]  0.8731427747 -0.1527436830  1.5868303196 -0.6074129269  0.7050365611
##  [551]  0.2146020896  0.4778398242 -0.0925022735  0.1083460435  0.6587346683
##  [556] -0.1880329034 -0.4480658056 -0.4024094663 -0.3782278830  0.5406852259
##  [561] -0.5215104927  0.5153831896  0.4070832504  0.1946333877  0.9709480590
##  [566]  0.9429335645  0.5354925028  1.0892860425  0.6332343785  0.4016385340
##  [571] -1.8071893064 -0.4549122337  1.0743692582 -0.1819958209  0.2495136982
##  [576] -0.1913318453  0.7715411434 -0.3190375048 -0.4290728266  0.1333845140
##  [581]  0.4419387059  0.1805554906  0.3323994374 -0.4642037224 -0.3749677391
##  [586] -0.0444157310 -0.1830268436  0.7516752429 -0.1145503010  0.1832180279
##  [591]  0.0672646330  0.2066075623 -0.0730188065 -0.1169295665  0.7764926499
##  [596]  0.2480085534 -0.0632277297 -0.1431160016  0.1023258260 -0.3780554995
##  [601]  1.0998100986  0.1561743508  0.9020963103  0.4695293212  0.0670856254
##  [606]  0.0403962832  0.4573740721  0.7663445528 -0.3225026547 -0.1305114057
##  [611]  0.8794517010  0.0122977862  0.1933346293  0.5769073653  0.3599172200
##  [616] -0.4596304100  0.7870363140 -0.2854873694  0.0270671031 -0.3081569935
##  [621] -0.1070301664 -0.8827081879 -1.1364201113  0.0851040897  0.7282048826
##  [626]  0.3501515525 -0.1319007737  0.5030328221  0.6121991726 -0.6034880948
##  [631]  0.4421070583 -0.9257308947  1.3529000891  1.0548071426  0.3928202079
##  [636]  0.0843673421 -0.1337968111  0.1062626951 -0.2058192782  0.1244610785
##  [641]  0.5164579469  0.1754581591 -0.0503076382  0.3259634708  0.0604627881
##  [646]  1.2632956458  0.1633220879 -0.3898350207 -0.5092460656  0.3004881960
##  [651]  0.0535080897 -0.4080303888  0.3183433025 -0.0119253218  0.4744455463
##  [656] -0.3862006947 -0.3201651888  0.2344906597 -0.4629423627  0.4512475946
##  [661]  0.3473162231  0.3608596811 -0.0524368918 -0.2122000596  0.5195693152
##  [666]  0.4813610280  0.3928691033  1.0469753785 -0.5788965647  0.4844384821
##  [671]  0.2785458457  0.3708962796 -0.5872450805 -0.2353729016  0.8070303810
##  [676] -0.2075588535 -0.0498423943 -0.5595196719 -0.0101666823 -0.0115849954
##  [681] -0.0553591112  0.8272317665 -0.8757089863 -0.1541066881  0.2781045694
##  [686]  0.0440753533  0.7876612187 -1.1647764283  0.6991845229  0.1426254939
##  [691]  0.2147147733 -0.1306955701 -1.0319552453  0.1528880474  0.5124226484
##  [696] -0.0474876758  0.0945817008 -0.1491549112  0.6793419828  0.7871984688
##  [701] -0.1564900444  0.0712692649  0.1148185736 -0.0467059361 -0.4717597161
##  [706]  0.0497864506 -0.7736923280 -0.1311597368 -0.1968446719  0.0300837496
##  [711]  0.3199799459  0.4070470131  0.4083024663  0.1535624188  0.2811645341
##  [716]  0.1552138556  0.2756422210 -0.4740789856 -0.1482976322 -0.5996738243
##  [721]  0.8334210366 -0.2372479852  0.4259060208 -0.2298099728  0.4385517377
##  [726]  0.0789589573  0.6488583831 -0.3247594207 -0.3365886197  0.0664602164
##  [731]  0.1564753908 -0.3861328782  0.0248191267 -1.0514262840  0.0754361239
##  [736] -0.0040999280  0.2471870773 -0.4625416839 -0.0402406178 -0.1942056814
##  [741]  0.1466327403  0.3814291088  0.1805818158 -0.2896082297 -0.0593155198
##  [746]  0.0825157523  0.4894444633 -0.3003471975 -0.1416050957  0.8535077458
##  [751] -0.2454105784 -1.0351769240  0.2436625011 -0.6584647758  0.4707202289
##  [756]  0.4338463130 -0.1000718543  0.0277304029  0.5284994281  0.1159434683
##  [761] -0.2099937196  0.5253290914  0.0682207784 -0.4300635195  1.1687969176
##  [766] -0.0593155198 -0.3019312389  1.1148575567 -0.4266239055  0.2001570313
##  [771]  0.2995569897 -0.2505673101  0.0411486038 -0.6414927592  0.5052947797
##  [776] -0.0870391545 -0.5620103569  1.0283582988  0.4027943656 -0.4141479442
##  [781] -0.6167865372 -0.0115687242  0.1509062730 -0.6291769453 -0.1564567839
##  [786]  0.2961896995 -0.5396831776  0.0799440846 -0.6567649479  0.6461550493
##  [791] -0.4864641986 -0.4754741659 -1.0686263372 -0.2244904161 -0.1410570589
##  [796]  0.5801233647  0.5520510750  0.0358379977 -0.1891728656 -0.3747416215
##  [801]  0.0806561838 -0.4125978560 -0.4107867482  0.1101359657  0.2847964667
##  [806] -0.1756318518 -0.0602110985 -0.6865718707  0.3189981044  0.2646526447
##  [811]  0.1275160183  0.1429291657  0.6549393190  0.7744402637 -0.1961780823
##  [816] -0.7479009759 -0.4300635195  0.9122891310  0.3247049728  1.3520389185
##  [821] -0.1208459993  0.1616787594  0.0072364367 -0.3619332633  0.0058034037
##  [826] -0.3458838221 -0.0135054157  0.2825645012 -1.0194278386  1.0770275042
##  [831]  1.2415668339 -0.1640530199 -0.9562128353  0.6494033957  0.4747488856
##  [836]  1.7327116142  0.4415614543 -0.3174067924  0.4016369405  0.1326563148
##  [841]  0.2478691716 -0.1711070538  0.6795644185  0.1343612721  0.2736679537
##  [846]  0.3054431813 -0.0193750781  0.6347096548  0.4583463440 -0.0079750337
##  [851] -1.1525299167  0.0858344081 -0.2564475803 -0.1079823548  0.8411609701
##  [856] -0.1781123294  0.0190311588 -0.5267862538 -0.6352420699  0.0502033818
##  [861]  0.1095380080  0.5501614501 -0.0318085935 -0.6636515061 -0.0031118638
##  [866]  0.5057863587 -0.1210151625  0.3625721469  1.4626550339  0.4669310013
##  [871] -0.3438195284 -0.0518441819 -0.1823174837  0.2910831211  1.0889382726
##  [876] -0.6782072086  0.2338881552  0.3437848867 -0.2458466256  0.3109816656
##  [881]  0.0881532622 -0.7719610630  0.1832854632 -0.1110412911  0.8543367341
##  [886]  0.5621119759  0.3274794384  0.1998999399 -0.2216343691  0.0789162604
##  [891]  0.6850554464 -0.3175328632  0.4423926980  0.6258344579 -0.1199421954
##  [896]  0.1972451195 -0.1703146232 -0.6385403731 -0.1459489087  0.1078750861
##  [901]  0.3822523438 -0.8470223102 -0.1318744976 -0.5307046227  0.0388176120
##  [906] -0.0901006214  0.2998764033 -0.4312494909  0.5315027228 -0.1500692873
##  [911] -1.0285300531  0.5040663821 -0.4043007566 -0.3419706701  0.6889132855
##  [916] -0.1031834573  1.5683285474  0.4166803223 -0.6705099293  0.2614656916
##  [921] -0.2194919746  0.5941238316 -0.0385046142 -0.4330981066 -0.9326607403
##  [926]  0.4029165558 -0.7821763801  0.6164080706  0.3827608214  1.4816061853
##  [931]  0.2693554054 -0.4424900919  0.3768348041  0.4066192856  0.5029451615
##  [936]  0.3712717567  0.4154185269  0.2905858204  0.7574273931  0.1455326819
##  [941] -0.3645126007  0.9586309999 -0.7594870077 -0.0497290938  0.1151123133
##  [946]  0.0416805402  1.0197296526 -0.8609497882 -0.0272966314 -0.1601853282
##  [951] -0.0477696222 -0.5509084745 -0.0218525495 -0.3207526834  0.2989155350
##  [956]  0.3741521368 -0.0827359259  0.4157305356  0.2128293220 -0.7938716679
##  [961] -0.3920413895 -0.0334258484 -1.2385076471  0.3372441934  0.1931582068
##  [966] -0.4336027581  0.4701648147  0.5550196257  0.0751623096 -0.4262805139
##  [971]  0.6865420238 -0.0830891916  0.3887872169  0.5561894489  0.1384585120
##  [976]  1.0436675691 -0.6007624922  0.7036625419  0.4561734087  0.4520265138
##  [981] -0.2288875074  0.2859153948 -0.4289466873  0.3681296530 -0.2799140559
##  [986]  0.3324644395  0.3063598607  0.3513730978 -0.0773173575 -0.1100955977
##  [991]  0.4033569465  0.7424794738  0.1195370055 -0.0374199134  0.6331012675
##  [996]  0.0268183282  0.6921290053 -0.4102615511  0.4985819991  0.7835635919
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
##   0.47438748   0.23930433 
##  (0.07567468) (0.05350631)
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
## [1]  0.3419889 -0.1459430  0.2720352  0.3328103  0.4932827  0.7696521
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
## [1] -0.0018
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9059346
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
##     original     bias    std. error
## t1*      4.5 0.01611612   0.9130816
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 2 3 5 7 8 9 
## 1 2 1 2 1 1 2
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
## [1] 0.0045
```

```r
se.boot
```

```
## [1] 0.9077961
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

