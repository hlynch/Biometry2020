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
## 0 1 2 3 4 5 8 9 
## 1 1 2 1 1 1 2 1
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
## [1] -0.0264
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
## [1] 2.759762
```

```r
UL.boot
```

```
## [1] 6.187438
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.8   6.1
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
##    [1] 3.8 4.7 4.2 5.1 6.0 4.5 5.1 4.0 6.2 4.5 6.4 3.2 4.6 4.4 4.7 5.3 4.7 5.4
##   [19] 4.3 2.9 5.9 4.1 5.3 3.0 4.1 4.9 3.8 2.9 4.5 5.9 3.9 4.1 5.3 4.7 3.5 4.0
##   [37] 2.9 5.5 5.9 3.8 4.1 4.7 5.1 5.2 3.3 5.7 4.0 6.9 5.1 3.0 5.9 4.1 6.5 3.6
##   [55] 4.6 4.4 4.5 3.1 3.1 5.0 3.8 5.9 4.9 4.3 4.3 3.8 4.4 3.9 3.1 5.0 5.0 5.3
##   [73] 4.8 4.4 5.0 3.4 5.0 4.8 4.9 5.1 3.4 4.4 5.8 5.3 3.3 3.5 4.1 4.3 3.8 2.2
##   [91] 5.3 2.6 5.8 5.0 3.0 5.8 4.8 4.0 4.9 4.4 5.8 4.2 4.5 3.7 4.5 5.1 5.4 4.1
##  [109] 5.0 5.2 4.5 3.3 4.5 5.9 4.5 4.3 3.4 3.7 5.6 3.0 4.5 3.7 3.7 2.5 4.8 4.5
##  [127] 4.4 4.1 4.9 2.6 4.8 4.5 3.2 4.1 6.7 5.0 5.0 5.5 3.5 6.5 5.7 4.1 3.7 4.7
##  [145] 5.6 4.5 5.1 3.8 3.8 4.6 3.7 4.4 3.7 4.1 3.2 6.2 4.8 5.3 6.1 4.5 4.5 5.4
##  [163] 4.9 4.7 4.6 4.0 3.8 5.0 3.5 5.0 4.6 5.0 3.1 2.2 5.6 4.5 3.2 6.5 4.0 4.6
##  [181] 6.0 3.5 4.5 3.9 5.0 5.0 5.0 3.0 4.0 5.0 3.1 4.3 3.9 4.1 3.6 3.7 4.9 4.6
##  [199] 6.0 3.5 5.4 4.0 5.5 3.8 4.4 5.0 5.5 3.1 4.0 5.1 3.9 5.1 4.2 4.5 3.0 3.6
##  [217] 6.6 4.8 3.5 4.7 6.6 4.1 4.5 4.1 4.8 5.4 6.5 3.1 3.8 6.0 5.1 4.7 3.7 3.8
##  [235] 4.0 4.5 6.2 2.9 4.3 4.4 5.5 2.2 3.7 5.4 4.8 3.4 3.8 4.0 4.5 3.0 4.4 4.7
##  [253] 4.3 4.2 4.6 3.7 5.8 5.1 4.5 6.0 4.3 5.7 4.0 5.4 5.3 3.4 4.8 4.0 3.9 4.8
##  [271] 2.9 4.8 4.2 4.6 5.7 3.3 3.7 3.7 4.2 6.2 3.8 4.2 4.9 3.6 5.1 5.3 3.9 4.4
##  [289] 4.2 3.3 2.6 3.5 4.0 5.7 5.2 6.1 5.5 3.9 6.4 3.7 2.5 3.4 5.1 4.1 4.3 4.7
##  [307] 3.8 5.3 2.5 3.8 5.0 4.3 4.1 5.2 3.5 3.1 5.0 2.7 4.4 5.2 3.9 3.5 5.0 4.0
##  [325] 4.1 4.2 3.0 2.5 5.4 4.9 3.5 4.9 4.3 3.8 4.3 4.4 4.0 4.1 4.5 4.1 5.9 6.7
##  [343] 4.5 4.8 4.7 1.2 5.4 4.1 5.5 3.7 3.2 4.2 3.7 4.8 4.8 3.2 7.0 4.7 4.6 5.3
##  [361] 3.7 4.8 4.8 3.9 4.8 4.6 4.8 4.1 3.1 5.1 3.5 5.9 5.1 6.8 5.2 4.7 3.6 5.3
##  [379] 4.9 5.2 3.5 4.7 5.4 3.9 4.9 3.7 4.9 3.9 3.7 3.6 4.4 3.9 2.9 3.3 4.6 4.2
##  [397] 5.3 5.2 3.1 4.1 3.3 6.3 6.3 4.3 4.3 7.0 4.1 4.8 4.5 5.4 3.2 5.2 5.3 4.4
##  [415] 5.1 5.1 5.5 3.9 5.2 5.1 6.0 4.1 4.9 4.7 4.5 4.4 5.3 4.8 3.7 5.8 5.5 6.3
##  [433] 3.8 4.4 3.6 4.3 3.4 5.3 4.3 4.0 3.6 6.1 4.6 4.7 5.4 2.8 5.1 5.1 3.6 4.1
##  [451] 4.9 4.6 4.7 3.5 4.4 3.7 3.3 5.0 3.5 4.9 5.1 6.4 3.4 3.7 4.0 3.6 3.8 3.4
##  [469] 5.5 4.3 3.9 4.7 4.4 4.5 4.2 4.4 3.0 4.6 5.0 3.7 4.9 3.4 4.2 5.9 2.7 3.6
##  [487] 4.7 4.3 2.6 3.4 3.7 4.0 4.8 5.9 6.0 3.6 4.6 4.0 4.3 4.3 5.0 5.6 3.1 5.0
##  [505] 3.9 4.0 5.8 5.2 5.6 5.9 5.2 4.2 5.1 2.6 3.6 3.9 3.7 3.1 4.1 4.7 5.9 2.0
##  [523] 4.4 5.2 4.0 5.3 4.9 4.9 3.3 5.7 3.9 3.3 5.1 4.7 5.4 4.0 2.7 6.5 4.8 4.6
##  [541] 4.8 4.1 3.9 4.1 3.6 3.2 3.9 4.4 4.7 6.3 4.6 3.9 3.4 5.0 4.7 5.1 4.6 2.4
##  [559] 2.8 3.7 4.8 4.3 4.7 4.2 4.2 3.3 4.2 4.7 5.1 2.9 5.2 4.9 4.1 4.3 4.6 5.1
##  [577] 3.4 4.1 4.0 4.3 6.3 4.3 5.3 4.4 4.7 5.0 3.4 3.3 5.4 4.0 3.5 5.4 4.6 4.9
##  [595] 3.2 4.7 6.0 5.4 5.4 3.2 4.7 3.5 6.6 3.9 2.9 3.8 3.8 6.5 5.1 3.9 4.3 4.4
##  [613] 3.7 2.4 5.1 4.4 5.4 3.6 4.0 3.5 4.5 5.5 5.2 4.0 4.0 3.6 3.1 5.6 3.7 4.9
##  [631] 4.5 3.8 3.1 4.5 4.6 5.3 4.0 5.5 5.1 3.8 4.1 4.5 3.4 4.8 3.8 3.4 4.1 5.8
##  [649] 3.6 3.8 5.5 3.7 5.2 6.2 3.3 3.9 5.9 5.8 3.7 3.9 4.3 3.2 4.7 4.8 4.5 3.5
##  [667] 3.2 4.5 4.9 5.4 4.2 5.2 4.3 3.4 4.2 4.4 5.6 4.6 4.1 3.5 5.7 5.4 4.8 4.9
##  [685] 5.3 5.0 3.5 5.8 4.8 5.8 5.0 4.8 3.9 3.8 5.1 5.9 5.0 4.3 6.7 5.7 5.4 4.1
##  [703] 4.1 6.5 5.3 5.3 5.8 5.0 4.8 3.3 5.1 4.1 5.0 4.4 6.3 3.7 4.1 5.4 5.2 2.6
##  [721] 4.4 3.4 4.7 4.7 3.7 4.9 4.5 3.8 4.7 5.5 4.4 4.5 6.0 4.5 5.8 5.6 5.2 5.0
##  [739] 5.7 5.3 5.2 5.5 5.7 5.7 3.6 4.5 3.9 3.8 4.8 3.8 6.4 4.4 4.2 2.5 4.6 5.3
##  [757] 5.4 4.0 4.3 4.1 4.3 4.3 4.8 4.6 3.7 4.0 7.3 5.1 4.3 3.9 4.2 4.6 5.5 3.8
##  [775] 4.3 4.7 2.8 3.7 4.0 4.2 4.7 4.0 4.1 3.6 4.9 6.1 3.7 5.5 2.4 4.5 3.3 4.0
##  [793] 2.7 5.0 2.8 6.0 4.4 3.8 5.8 3.3 5.9 3.6 3.0 4.8 5.4 4.8 4.0 5.2 4.2 4.9
##  [811] 4.6 4.0 3.9 4.9 6.8 3.5 4.9 5.2 4.8 4.4 5.2 3.0 5.7 5.3 3.0 4.6 5.0 2.9
##  [829] 5.5 5.6 4.4 4.3 6.0 2.6 5.4 4.8 5.0 5.3 4.1 3.3 4.8 3.4 5.3 4.3 5.5 2.0
##  [847] 2.9 5.4 4.7 4.3 4.5 2.5 5.4 5.3 5.8 5.7 7.3 5.3 3.8 3.9 4.7 5.3 6.1 3.7
##  [865] 4.6 3.1 4.8 5.5 3.3 3.9 6.0 4.0 4.8 7.0 3.1 4.5 3.7 5.7 3.3 3.9 4.1 4.0
##  [883] 3.6 3.7 3.7 3.8 6.0 4.1 4.8 5.0 4.3 3.5 4.4 3.8 6.0 4.3 3.5 4.3 5.5 4.3
##  [901] 4.9 3.1 4.6 6.1 5.1 3.3 4.3 4.2 5.3 5.2 4.4 4.5 4.1 3.5 4.4 4.4 3.5 5.7
##  [919] 3.3 3.8 5.9 3.0 5.2 3.3 4.3 4.2 4.5 3.7 4.8 5.2 4.0 5.6 3.8 4.4 4.4 2.4
##  [937] 5.1 5.1 3.8 5.4 5.9 3.6 5.5 4.8 3.3 3.6 4.0 3.3 5.0 3.8 3.8 5.3 3.6 5.0
##  [955] 4.7 4.8 4.8 4.9 5.0 3.4 4.3 3.9 4.2 4.4 4.0 5.6 3.6 2.6 3.2 3.9 3.4 4.2
##  [973] 4.3 5.6 5.0 4.8 3.3 3.1 3.0 3.8 6.4 4.5 2.9 3.5 3.9 4.2 4.4 6.0 4.5 5.2
##  [991] 4.9 4.6 6.3 4.9 3.9 5.2 5.1 4.5 5.6 4.5
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
##   2.7   6.4
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
##    [1] 5.0 5.2 5.1 5.6 4.1 6.3 4.5 4.9 2.3 4.6 4.8 5.4 4.7 5.7 4.4 4.1 4.1 3.4
##   [19] 4.7 4.6 3.5 5.7 4.0 4.6 5.6 3.1 5.8 3.6 4.3 3.8 4.0 4.4 4.2 4.1 5.6 4.1
##   [37] 4.8 3.4 4.3 5.8 4.2 5.8 4.2 6.3 4.1 4.7 4.5 4.5 5.1 3.9 5.1 5.0 4.9 3.7
##   [55] 3.7 5.5 5.6 4.6 3.0 3.0 3.8 5.5 3.4 5.4 6.3 4.9 4.2 4.8 4.8 4.2 4.4 4.3
##   [73] 5.8 4.8 5.7 5.8 3.5 5.2 4.8 3.6 4.8 4.8 4.2 3.3 4.0 3.1 3.5 5.8 4.9 4.0
##   [91] 5.0 4.1 4.1 3.8 4.4 5.7 4.6 5.4 3.3 5.2 3.7 3.8 4.7 5.9 4.1 4.2 4.3 4.2
##  [109] 4.1 5.2 2.0 5.8 4.5 6.0 4.0 5.0 4.7 4.3 4.0 5.3 4.4 4.7 4.3 3.6 5.3 4.8
##  [127] 3.9 5.0 5.2 5.1 4.9 2.8 5.7 3.9 3.7 6.4 4.8 4.3 6.6 4.7 2.9 5.8 6.5 4.8
##  [145] 4.9 4.4 4.6 4.1 2.4 4.8 4.5 3.9 4.3 5.3 4.9 3.3 4.9 4.7 7.2 6.8 6.3 3.2
##  [163] 4.2 3.6 5.0 3.5 7.1 2.8 4.9 5.4 4.2 3.0 3.8 4.4 4.9 4.0 3.7 5.5 5.2 5.9
##  [181] 3.5 3.9 4.1 4.6 3.9 3.7 5.6 5.6 3.0 4.8 7.0 4.6 3.3 5.4 4.7 4.3 3.4 4.2
##  [199] 3.3 4.1 5.3 4.9 4.3 6.2 4.3 4.7 4.9 3.7 3.9 4.2 4.3 5.3 4.5 5.2 5.6 4.9
##  [217] 5.5 2.8 5.1 3.2 4.2 4.0 4.5 4.3 2.9 4.3 4.4 4.2 4.5 5.1 3.7 5.3 3.3 4.0
##  [235] 4.9 4.6 3.8 2.0 5.6 4.0 5.9 5.0 5.2 4.2 3.7 3.9 5.7 4.3 5.0 5.0 3.0 5.2
##  [253] 5.9 4.0 2.4 4.1 4.0 6.2 3.3 3.7 3.9 4.5 6.3 4.7 5.1 6.6 4.1 4.3 4.9 3.5
##  [271] 4.6 4.1 5.0 4.2 4.6 3.3 2.2 4.8 5.0 5.0 3.3 3.0 5.0 4.6 3.4 3.8 5.8 5.7
##  [289] 3.5 4.8 5.3 5.8 3.8 3.4 5.8 4.1 3.8 5.2 5.1 5.0 5.1 4.6 4.2 3.6 5.0 4.6
##  [307] 3.6 4.0 4.9 3.8 4.7 5.2 6.0 3.8 4.8 4.5 4.2 4.6 3.1 3.8 5.0 3.9 4.2 5.6
##  [325] 4.0 2.9 3.3 4.4 3.4 4.1 2.2 4.2 3.9 5.6 5.4 5.5 4.6 6.6 4.6 4.2 4.5 4.5
##  [343] 3.6 4.7 4.8 5.0 4.0 4.4 3.2 4.5 2.4 3.5 6.1 5.7 5.5 3.7 5.9 5.1 4.2 4.9
##  [361] 4.6 4.0 3.3 4.0 4.4 4.3 3.4 4.6 4.8 4.6 4.1 4.8 4.6 4.1 3.4 4.5 5.5 4.8
##  [379] 3.5 2.8 4.6 3.9 4.4 5.3 5.1 4.0 5.3 4.3 4.2 3.2 5.5 4.0 3.6 4.2 2.8 5.3
##  [397] 4.8 4.5 4.9 3.7 4.2 4.5 4.8 6.0 4.1 4.4 3.4 6.4 5.3 5.5 3.9 6.2 5.1 3.3
##  [415] 4.8 4.1 4.0 6.4 4.3 6.1 4.9 3.4 5.2 4.5 4.3 4.9 5.2 3.7 5.7 3.4 4.7 4.1
##  [433] 4.2 5.4 3.9 3.0 4.7 3.7 4.1 5.4 3.4 5.7 4.9 4.4 4.7 3.6 4.5 3.8 3.7 3.1
##  [451] 4.4 5.0 3.8 5.2 5.8 4.9 4.9 3.5 5.1 4.9 6.1 3.9 6.0 3.8 4.3 5.3 3.8 5.7
##  [469] 4.9 5.0 2.0 4.1 4.8 5.4 6.5 4.3 2.9 5.2 2.7 3.9 3.7 4.1 4.4 5.9 4.7 5.6
##  [487] 4.9 4.6 5.1 3.8 4.2 4.2 4.4 4.0 4.9 4.0 5.4 5.2 4.5 2.3 3.4 5.8 4.1 4.1
##  [505] 3.9 4.6 4.1 5.5 4.2 4.1 3.3 6.8 3.9 4.2 3.7 3.7 6.4 3.2 5.8 5.7 5.1 5.4
##  [523] 3.6 4.0 6.0 4.4 4.7 3.0 6.0 4.5 2.8 4.6 4.1 4.6 5.1 3.9 4.6 5.0 4.6 3.7
##  [541] 4.7 2.3 3.6 5.6 4.8 5.4 3.5 3.2 3.9 4.3 3.8 4.0 6.6 4.4 4.9 4.9 4.2 3.4
##  [559] 3.9 4.0 6.4 4.3 4.6 5.6 4.4 4.5 4.7 3.9 3.1 4.2 3.8 3.0 4.5 5.7 4.5 4.7
##  [577] 4.6 4.2 5.3 4.8 5.8 3.6 4.3 5.5 4.1 3.9 4.8 6.3 4.0 4.5 4.7 1.9 3.8 6.4
##  [595] 5.1 4.5 5.2 4.7 5.7 2.9 4.9 6.4 3.6 6.0 6.3 4.3 4.4 4.8 5.5 5.0 3.9 3.1
##  [613] 4.5 5.4 4.5 6.1 4.7 3.6 4.1 4.1 5.4 3.6 4.8 4.8 5.1 5.1 4.0 3.6 4.9 3.9
##  [631] 3.7 4.6 5.1 3.7 4.8 5.3 3.1 5.9 3.8 3.6 3.6 4.9 4.3 4.2 4.0 3.6 3.8 4.7
##  [649] 4.0 4.0 5.0 4.5 4.8 5.4 3.6 5.4 5.8 4.9 2.7 5.7 3.5 4.6 2.8 4.5 3.0 3.2
##  [667] 5.3 4.9 4.9 5.9 4.0 6.0 4.4 3.9 3.1 3.0 4.2 5.5 4.6 4.9 5.2 4.9 4.6 5.0
##  [685] 4.4 5.3 4.1 4.7 5.3 4.4 3.8 5.6 4.4 3.5 4.6 4.3 3.3 4.9 4.1 5.1 6.0 3.1
##  [703] 4.4 4.2 4.3 3.3 4.0 4.9 4.1 4.9 4.8 4.6 6.3 5.4 3.9 4.6 3.9 4.2 4.5 3.5
##  [721] 4.7 4.2 5.5 4.7 5.7 4.4 5.5 5.1 4.1 4.6 5.6 5.2 4.1 3.6 4.9 4.0 3.2 4.0
##  [739] 4.2 3.5 4.5 3.9 3.3 5.0 2.4 4.9 3.3 5.7 4.4 6.7 3.8 4.2 3.9 4.1 5.0 3.8
##  [757] 4.6 5.2 4.2 3.6 5.0 4.8 4.0 6.7 4.3 4.2 3.9 3.8 3.7 4.3 5.0 4.0 5.1 4.4
##  [775] 4.3 6.7 2.9 3.2 4.5 5.9 4.2 5.3 4.1 4.9 4.3 4.3 4.9 4.5 3.4 4.5 5.9 4.4
##  [793] 4.4 5.1 6.4 3.6 2.5 3.5 3.5 3.4 3.9 6.0 4.7 2.7 3.5 4.7 3.2 6.3 5.9 4.7
##  [811] 6.0 5.9 5.7 4.4 4.0 4.2 3.6 4.7 4.6 5.2 4.1 4.5 5.8 5.1 6.7 3.4 4.2 3.6
##  [829] 3.4 6.3 4.3 3.5 2.4 4.1 3.7 5.2 4.9 6.8 4.4 5.2 4.9 5.2 4.6 5.4 6.5 3.6
##  [847] 4.5 4.6 5.6 4.5 5.6 5.2 5.2 3.0 5.0 3.8 4.5 3.4 4.3 5.9 4.7 5.8 3.8 3.8
##  [865] 5.6 5.5 5.7 4.8 2.2 4.2 4.4 4.8 5.4 5.2 4.7 4.0 3.9 4.5 3.9 6.6 3.7 6.2
##  [883] 5.0 2.9 4.4 5.1 5.7 5.8 4.3 4.2 5.6 3.3 4.4 4.3 4.8 5.3 3.9 4.2 4.6 2.9
##  [901] 5.1 6.3 3.8 4.0 4.6 4.3 3.4 4.0 5.0 5.7 5.3 5.4 3.9 2.9 5.8 3.9 4.1 2.9
##  [919] 3.2 6.0 5.8 5.2 5.8 3.5 5.0 3.2 5.2 4.0 5.6 7.0 3.0 4.8 4.0 5.0 5.1 4.3
##  [937] 3.0 4.5 6.1 5.4 5.6 4.4 3.8 2.8 4.6 5.6 5.2 5.3 3.9 5.6 4.2 4.5 4.9 4.9
##  [955] 3.5 3.5 6.5 5.1 3.6 4.2 5.1 5.8 4.4 4.5 4.9 4.0 4.7 4.5 4.6 4.2 5.2 4.7
##  [973] 5.9 4.0 6.2 4.3 3.9 5.2 2.7 4.1 4.6 3.9 3.8 5.3 4.2 5.0 3.8 2.8 4.3 4.8
##  [991] 5.3 3.9 4.3 4.3 4.9 3.5 4.5 5.3 4.4 5.4
## 
## $func.thetastar
## [1] 0.0185
## 
## $jack.boot.val
##  [1]  0.54346591  0.38397790  0.30835913  0.17832898  0.09646739 -0.05406162
##  [7] -0.15967302 -0.25093168 -0.31615854 -0.49337176
## 
## $jack.boot.se
## [1] 0.9525436
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
##    [1] 3.0 3.5 3.9 4.7 4.7 5.4 5.8 4.4 3.7 3.8 3.4 5.3 3.9 3.2 3.3 3.2 3.1 4.5
##   [19] 5.0 6.4 5.8 5.2 3.9 3.7 4.5 2.9 6.2 4.0 4.4 4.9 3.9 3.7 5.5 3.7 4.0 4.8
##   [37] 4.4 2.2 4.7 4.1 6.1 3.9 4.1 4.0 6.3 3.9 4.8 3.7 3.8 4.3 5.1 4.2 6.1 4.7
##   [55] 4.9 6.0 3.2 3.8 4.0 5.1 3.6 3.6 4.4 5.1 4.8 2.7 4.2 3.2 4.4 5.1 4.9 3.8
##   [73] 4.7 5.2 4.0 5.0 4.6 5.1 3.8 5.1 5.3 4.0 5.7 4.1 3.6 3.1 3.4 4.8 4.3 4.7
##   [91] 4.0 5.5 4.7 4.7 3.6 4.2 4.7 5.3 5.7 4.0 4.1 5.0 5.4 4.0 2.7 5.3 3.5 4.3
##  [109] 3.8 4.9 5.3 3.1 6.3 6.0 4.6 4.2 2.5 6.0 4.7 3.9 5.4 3.9 3.3 3.8 3.5 3.3
##  [127] 4.0 5.1 5.2 4.4 4.7 4.0 3.0 4.9 3.9 4.8 4.3 4.8 4.3 3.0 4.8 3.2 4.4 5.5
##  [145] 4.1 3.3 5.2 4.2 6.3 6.0 4.5 4.8 6.4 3.6 3.2 4.0 5.0 3.5 4.8 3.3 5.6 5.4
##  [163] 3.7 3.3 3.5 4.2 4.3 5.4 3.3 4.8 4.5 4.7 5.5 4.1 4.4 4.8 4.2 4.8 4.2 5.5
##  [181] 2.4 3.9 5.3 5.1 6.1 2.9 4.0 2.9 4.8 5.5 4.5 2.2 3.2 4.3 5.4 4.3 5.6 3.5
##  [199] 4.2 3.3 3.9 4.6 4.0 4.7 4.4 4.8 4.1 4.4 4.4 4.0 4.0 3.2 4.6 5.3 3.6 4.8
##  [217] 5.7 3.1 6.0 5.7 4.6 3.6 2.6 5.5 5.3 5.6 4.5 4.8 4.6 3.1 4.6 3.7 3.2 5.2
##  [235] 3.7 4.2 3.9 4.5 5.0 4.3 4.6 3.1 5.2 5.5 5.1 4.3 3.8 3.0 5.6 3.7 4.7 4.9
##  [253] 4.4 2.3 4.1 2.6 5.8 3.4 3.5 2.9 4.6 5.4 3.9 4.5 5.5 3.8 4.8 5.2 3.8 3.4
##  [271] 4.8 3.9 6.1 5.0 4.8 4.8 4.3 4.8 5.0 4.5 6.0 4.6 2.6 6.3 4.7 4.4 3.2 2.1
##  [289] 4.8 5.0 4.7 4.2 4.9 4.1 4.2 6.4 6.1 3.4 3.5 4.4 2.7 3.8 3.0 5.7 5.8 5.6
##  [307] 5.3 3.2 5.1 3.7 3.9 5.0 4.4 4.9 4.7 2.7 5.3 4.4 4.1 3.9 5.9 4.5 4.1 6.3
##  [325] 5.5 3.5 5.4 4.1 3.5 5.0 4.1 3.8 5.5 3.5 4.9 4.1 2.6 4.1 5.5 6.3 4.2 4.7
##  [343] 4.4 5.9 5.7 4.7 5.4 6.1 4.5 3.9 5.0 6.8 4.4 4.9 4.7 5.6 4.8 6.5 3.5 4.4
##  [361] 3.2 3.6 5.7 4.6 5.1 2.8 5.6 4.2 5.0 5.3 3.8 3.8 4.6 4.5 3.7 4.8 5.2 6.2
##  [379] 3.9 3.2 4.5 2.9 5.2 5.3 5.4 3.2 5.9 3.4 4.6 2.7 4.2 4.7 4.1 4.4 5.3 4.1
##  [397] 4.3 3.8 4.1 4.3 4.9 5.0 4.5 5.1 3.8 6.5 6.5 4.8 2.6 4.3 3.0 4.4 3.3 5.2
##  [415] 4.7 5.0 5.3 5.6 5.0 3.7 6.3 4.1 4.5 5.4 4.5 4.5 4.2 3.3 4.4 5.6 5.2 4.8
##  [433] 3.0 4.1 5.3 4.1 3.8 5.8 5.0 5.1 3.4 4.7 5.8 5.4 5.4 6.1 5.7 6.3 4.0 3.5
##  [451] 4.9 5.4 3.5 4.4 4.6 5.7 4.3 5.5 4.5 4.6 5.0 3.7 3.5 3.0 4.9 4.2 5.2 4.8
##  [469] 5.3 4.2 6.0 4.1 4.0 4.2 4.6 5.7 5.1 5.7 4.0 4.6 6.0 4.8 5.5 6.3 3.2 5.2
##  [487] 4.7 4.6 4.2 3.5 4.2 3.4 3.8 4.9 3.6 3.7 3.0 5.6 4.3 4.1 4.6 4.6 6.8 4.9
##  [505] 3.8 6.8 4.0 5.6 4.0 5.4 5.7 4.7 5.1 5.3 6.4 3.7 4.1 4.3 4.4 5.9 4.6 3.7
##  [523] 3.4 5.4 4.6 5.2 3.9 4.0 4.3 4.9 2.5 4.1 4.9 5.5 4.6 5.2 4.4 2.2 4.8 3.6
##  [541] 4.2 3.3 4.0 3.6 5.5 5.1 4.3 4.8 3.9 4.2 3.4 3.2 4.5 4.9 4.0 4.3 4.5 5.0
##  [559] 3.6 5.2 4.2 4.4 4.6 3.2 4.7 4.2 5.0 4.6 3.9 4.1 5.0 3.4 3.5 4.4 4.5 2.6
##  [577] 4.5 3.9 5.1 5.2 3.2 3.6 4.0 4.8 4.5 3.7 5.0 4.1 3.9 1.9 3.6 4.2 5.2 4.6
##  [595] 3.8 3.7 5.0 3.9 5.2 5.6 3.5 3.7 5.5 4.4 5.0 5.8 5.1 5.3 2.4 3.7 5.1 5.2
##  [613] 3.2 5.6 5.5 6.3 3.4 4.5 4.1 5.9 5.1 4.8 4.3 5.1 5.9 4.6 3.4 5.8 2.6 6.3
##  [631] 3.7 4.7 5.2 3.7 3.6 4.6 4.6 5.4 4.6 2.8 4.3 3.5 3.4 5.3 5.4 3.9 3.4 3.5
##  [649] 4.1 4.9 3.3 4.2 4.9 2.0 2.4 5.5 5.9 5.2 4.5 5.5 3.8 4.8 4.6 5.2 4.8 4.6
##  [667] 4.8 5.4 3.2 5.0 3.9 4.6 4.8 4.5 3.2 4.1 3.8 5.9 4.9 6.2 4.0 4.6 4.1 3.9
##  [685] 4.0 3.5 4.3 5.9 5.1 5.8 3.1 4.9 4.4 6.1 4.5 3.2 4.3 4.3 3.8 5.0 5.2 4.5
##  [703] 5.5 5.7 3.4 4.5 5.3 4.3 5.2 5.2 4.7 5.7 3.5 5.3 3.2 5.2 4.1 3.8 6.7 4.8
##  [721] 2.8 5.2 3.6 4.4 3.7 5.2 5.3 4.6 2.8 5.1 4.8 5.3 3.5 4.6 3.3 5.4 4.1 3.5
##  [739] 3.4 4.3 4.7 4.9 4.8 4.6 4.1 5.2 3.3 5.0 5.9 4.5 5.7 3.8 2.9 6.0 5.6 4.0
##  [757] 4.5 4.9 6.1 7.2 5.2 4.7 3.9 4.6 4.8 4.1 5.1 4.7 5.0 4.7 5.0 4.5 3.9 4.3
##  [775] 4.3 6.0 6.9 3.9 4.0 3.8 3.0 3.9 4.5 4.1 3.3 4.1 5.7 4.2 4.1 4.5 3.6 3.8
##  [793] 5.2 4.8 5.9 6.0 6.0 3.9 3.9 5.6 4.0 4.1 4.3 4.3 4.0 5.7 3.8 4.2 4.2 3.9
##  [811] 4.4 4.0 3.5 5.9 4.5 4.5 3.7 4.0 2.9 4.6 4.3 5.3 4.9 4.6 3.9 4.8 3.7 4.4
##  [829] 4.8 5.5 5.0 4.7 4.1 4.2 5.4 4.9 5.9 5.1 3.5 3.9 4.0 4.7 5.4 4.8 3.8 4.6
##  [847] 6.1 3.8 3.0 4.2 5.3 3.9 4.5 5.4 2.8 5.3 5.4 4.4 5.1 5.0 5.0 4.4 3.4 3.7
##  [865] 5.3 5.1 5.4 3.1 5.0 3.8 4.3 4.5 4.7 4.5 5.0 4.6 4.1 3.3 5.4 4.5 4.1 2.5
##  [883] 4.1 3.6 4.2 5.0 4.0 4.2 5.0 4.8 4.0 5.2 4.6 3.2 3.7 5.0 4.4 4.0 5.0 2.9
##  [901] 4.8 5.6 4.4 5.7 3.3 4.3 4.4 6.1 3.5 3.0 5.4 3.6 3.9 5.2 4.3 4.2 4.4 4.6
##  [919] 3.9 2.9 4.6 4.7 4.5 3.8 4.1 3.6 4.6 5.6 5.4 5.5 3.6 4.3 6.0 4.4 4.1 5.0
##  [937] 6.5 4.3 5.9 4.5 4.0 3.0 5.8 4.5 5.3 4.9 3.2 4.8 5.6 4.8 5.2 2.3 5.1 4.9
##  [955] 4.5 4.0 4.6 4.5 6.1 3.0 4.4 6.2 5.4 5.2 4.8 4.0 4.5 3.8 4.2 4.4 3.5 4.5
##  [973] 3.4 4.4 4.0 4.0 5.6 4.3 4.0 4.4 4.8 4.2 3.4 5.8 5.4 4.0 4.1 4.3 5.1 4.5
##  [991] 4.6 2.6 4.1 5.0 3.9 5.8 4.6 6.3 4.8 3.6
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.4 5.4 5.3 5.1 5.1 5.0 4.8 4.7 4.5 4.5
## 
## $jack.boot.se
## [1] 0.9748846
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
## [1] 1.218362
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
##   1.3478014   2.4967154 
##  (0.5442186) (1.2160698)
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
## [1] 0.2660319 0.5202124 0.1872670 1.5599780 0.6891409 0.7980852
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
##    [1]  9.655713e-01  1.613660e+00  1.232400e+00  9.478591e-01  6.563075e-01
##    [6]  5.065514e-01  5.025901e-01 -1.809303e-02  2.374754e+00  1.207194e+00
##   [11]  1.438253e+00  1.521751e+00  8.259270e-01  1.008945e+00  4.808107e-01
##   [16]  1.013523e+00  8.887065e-01  8.427529e-01  5.937180e-01  8.673893e-01
##   [21]  1.476253e+00  3.494945e-01  1.120672e+00  7.529311e-01  1.113549e+00
##   [26]  8.490284e-01  1.118133e+00  1.451553e+00  7.983034e-01  5.055390e-02
##   [31]  8.908355e-01  9.789039e-01  5.630028e-01  1.138873e+00  1.812122e+00
##   [36]  1.033212e+00  4.907294e-01  1.708757e+00  5.051352e-01  1.147219e+00
##   [41]  1.603919e+00  1.852874e+00  1.146063e+00  5.360756e-01  6.987125e-01
##   [46]  4.308924e-01  1.972338e+00  1.459631e+00  8.791373e-01  2.674243e-01
##   [51]  1.629356e+00  2.712681e-01  1.331935e+00  5.811451e-01  1.640710e+00
##   [56]  7.657309e-01  5.164523e-01  6.223503e-01  1.508664e+00  1.477534e+00
##   [61]  5.982301e-01  1.418036e+00  1.285584e+00  6.137042e-01  1.667447e+00
##   [66]  2.528432e-01  9.358939e-01  7.872170e-01  1.285141e+00  1.853308e+00
##   [71]  8.012571e-01  1.439343e+00  9.256145e-01 -6.387186e-02  1.244885e+00
##   [76]  1.186462e+00  8.155790e-01  1.443150e+00  6.126101e-01  8.768639e-01
##   [81]  1.159634e+00  2.073842e+00  1.328417e+00  3.106329e-01  1.634444e+00
##   [86]  8.143059e-01 -1.012540e-01  1.080445e+00  6.850186e-01  7.118733e-01
##   [91]  8.320867e-01  1.046742e+00 -4.566800e-02  6.599575e-01  1.326256e+00
##   [96]  1.456416e+00  7.797307e-01  7.816517e-01  7.365783e-02  9.367171e-01
##  [101]  8.925460e-01  3.629180e-01  1.144357e-02  1.157212e+00  1.425668e+00
##  [106]  9.712222e-01  1.010265e+00  1.521018e+00  8.706898e-01  9.299689e-01
##  [111]  6.274310e-01  8.907324e-01  1.176084e+00  1.697967e+00  3.869188e-02
##  [116]  4.481136e-01  3.892897e-01  2.118869e+00 -1.431509e-01  5.466282e-01
##  [121]  1.038907e+00  4.884167e-01  7.162577e-01  7.776818e-01  3.656269e-01
##  [126]  9.725436e-01  1.073402e+00  6.678025e-01  1.229508e+00 -2.249406e-01
##  [131]  6.345998e-01 -3.194526e-02  1.475028e+00  1.409485e+00  4.340312e-01
##  [136]  4.777436e-01  1.062585e+00  3.576351e-01  9.691841e-01  6.596685e-01
##  [141]  7.207291e-01  6.458741e-01  3.813600e-01  1.534577e+00  1.135353e+00
##  [146]  1.410871e+00  8.917986e-01  1.005705e+00  8.421859e-01  1.683174e+00
##  [151]  1.366447e+00  1.105653e+00  4.730019e-01  1.076400e+00  1.240766e+00
##  [156]  1.127923e+00  7.138789e-01  1.416137e-02  1.120279e+00  5.052703e-01
##  [161]  4.577610e-02  5.300532e-01  3.483445e-01  9.379690e-01  1.240896e+00
##  [166]  2.576717e-02  1.079877e+00  2.680638e-01  1.121535e+00  1.133110e+00
##  [171]  1.263891e+00  1.084650e+00  1.238484e+00  1.266480e+00  7.849801e-01
##  [176] -1.647031e-01  7.730618e-01  6.960721e-01  1.667240e+00  1.115471e-01
##  [181]  1.034071e+00  1.167663e+00  1.968055e-01  1.383594e+00  8.983513e-01
##  [186] -1.447116e-01  2.025725e+00  5.477038e-01  1.924199e+00  3.195554e-01
##  [191]  1.312828e+00  1.535486e+00  7.158187e-02  8.399312e-01  8.109535e-01
##  [196]  1.514892e+00  4.067154e-01  6.983135e-01 -2.555858e-01  9.089166e-01
##  [201] -4.475382e-02  1.462535e+00  8.615801e-01  9.449621e-01  1.730882e+00
##  [206]  2.380807e-01  1.407506e+00 -2.587231e-02  7.686173e-01  1.140256e+00
##  [211]  4.145608e-01  9.576705e-01  2.703825e-01  5.969649e-01  1.011624e+00
##  [216]  1.084359e+00  8.735699e-01  1.055591e+00  1.555450e+00 -9.474952e-03
##  [221] -7.221351e-01  6.554367e-01  1.482670e+00  1.209850e+00  1.309987e+00
##  [226]  1.021334e+00  1.212482e+00  6.493595e-01  1.071356e+00  9.454680e-01
##  [231]  7.917974e-02  1.387236e+00  2.366263e-01  2.447368e+00  1.103793e+00
##  [236]  1.016528e+00  1.731990e+00  3.498733e-01  6.451635e-01  1.886050e-01
##  [241]  8.818328e-01  1.947828e+00  1.248831e+00  6.206422e-01  3.315839e-01
##  [246]  9.320870e-01  1.111243e+00  7.439119e-01  9.700979e-01  8.255923e-01
##  [251]  1.003720e+00  7.802682e-01  1.105614e+00  9.503037e-01 -3.150633e-01
##  [256]  6.008351e-01  7.745431e-01  6.083047e-01  8.269373e-01  1.210901e+00
##  [261]  5.993767e-02  1.215733e+00  6.692785e-01  1.490754e+00  1.228715e+00
##  [266]  1.183742e+00  7.246170e-01  1.978759e+00  5.352816e-01  1.406458e+00
##  [271]  3.392343e-01  7.281201e-01  8.651322e-01  1.346687e+00  8.622400e-01
##  [276]  1.315383e+00  1.510185e+00  2.018742e+00  7.727037e-01  5.289771e-01
##  [281]  1.072578e+00  1.935272e+00  3.746739e-01  6.585287e-02  8.122240e-01
##  [286]  1.246566e+00  1.203307e+00 -3.450319e-02  3.494805e-01  3.516492e-01
##  [291]  1.170851e+00 -1.287896e-01  1.528330e+00 -6.356738e-01  1.529185e+00
##  [296]  1.133836e-01  1.060365e+00  1.149286e+00  8.304564e-02  4.478156e-01
##  [301]  2.077511e+00  3.474356e-01  1.657446e+00  8.765720e-01  1.333447e+00
##  [306]  8.653329e-01  6.939100e-01  1.246962e+00  1.098114e+00  1.668201e+00
##  [311]  6.652277e-01  5.005105e-02  1.117015e-01  1.140865e+00  6.023109e-01
##  [316]  2.387949e+00  9.655713e-01  1.642550e-01  8.010916e-01 -4.321735e-02
##  [321]  1.418813e-01  8.252820e-01  1.393977e+00  1.428112e+00  7.872170e-01
##  [326]  3.129921e-01  3.134391e-01  8.850220e-01  5.299097e-01  1.199357e+00
##  [331]  1.290569e+00  7.508653e-01  5.345632e-01  8.462912e-01  1.266586e+00
##  [336]  2.021224e-01  1.202761e+00  6.701539e-01  1.206483e+00  4.671965e-01
##  [341]  1.314255e+00  8.426502e-01  2.001976e+00  2.147511e-01 -4.461841e-01
##  [346]  7.415751e-01  3.757365e-01  1.085824e+00  1.181533e+00  1.089663e+00
##  [351]  9.508036e-01  9.778220e-01  6.612173e-01  1.162781e+00  1.079163e+00
##  [356]  3.024129e-01  5.462483e-01  1.180014e-02  9.305274e-01 -2.368514e-01
##  [361]  1.851432e+00  1.006574e+00  2.739914e-01  1.703419e+00  1.219880e+00
##  [366]  1.267959e-01  1.036202e+00  1.894835e+00  9.937547e-01  6.159171e-01
##  [371]  4.610405e-01  6.940393e-01  7.405513e-02  1.172383e+00  6.659660e-02
##  [376]  4.748378e-01  8.509953e-01  1.275497e+00  6.372602e-01  1.178737e+00
##  [381]  5.187094e-01  4.495217e-01  1.440912e-01  9.532739e-01  1.321110e+00
##  [386]  5.318331e-01  1.205900e+00  1.080452e+00  1.065199e+00  1.507467e+00
##  [391]  1.367901e+00  1.270327e+00  1.735720e+00  1.136125e+00  6.290911e-01
##  [396]  1.406458e+00  6.326144e-01  1.139949e+00  4.697505e-01  1.446807e+00
##  [401]  5.605806e-01  5.094779e-01 -3.561606e-03  6.354823e-01  1.276520e+00
##  [406]  1.216898e+00  1.124361e+00  1.326177e+00  2.077461e-01  1.032384e+00
##  [411]  1.043701e+00  1.720780e+00  1.559483e+00  1.935272e+00  6.367261e-01
##  [416]  8.410032e-01  6.606455e-01  1.332776e+00  1.075539e+00  1.166174e+00
##  [421]  1.249433e+00  8.916024e-01  1.226325e+00  4.566763e-01  2.053428e-01
##  [426]  5.606716e-01  5.831233e-01  2.528006e-01  1.283257e+00  3.395643e-01
##  [431] -1.507756e-01  1.425588e+00  1.729354e+00  1.274804e+00  1.342513e+00
##  [436]  8.215768e-01  9.100759e-01  1.755499e+00  1.717809e+00  9.058144e-01
##  [441]  3.408506e-01  8.581449e-01  1.122831e+00 -4.634543e-01  1.526535e+00
##  [446]  1.475516e+00  2.459136e-01  1.493684e+00  6.607730e-01  8.134635e-01
##  [451]  3.070045e-01  1.547290e+00  3.875926e-01  8.887679e-01  8.106837e-01
##  [456]  6.201497e-01  8.441471e-01  9.270854e-01  1.084853e+00  6.160971e-01
##  [461]  9.043058e-01  8.576246e-01  1.382029e+00  3.616281e-01  1.187017e+00
##  [466] -1.966899e-01  1.040393e+00  1.295196e+00  1.880120e-01  3.395252e-01
##  [471]  7.656470e-01  1.637488e+00  1.384572e+00  8.045381e-01  2.401988e-01
##  [476]  1.243769e+00  6.007949e-01  4.273401e-01  1.740610e+00 -5.986583e-01
##  [481]  1.025483e+00  2.104939e+00  7.162577e-01  2.658800e-01  2.261753e-01
##  [486]  3.204561e-01  5.735073e-01  9.968686e-01  5.597440e-01  1.360511e+00
##  [491]  5.804333e-02  3.195629e-01  7.997836e-01  6.441587e-03  8.766636e-01
##  [496]  7.333006e-01  1.357085e+00  7.619309e-01  1.275539e+00  1.172596e+00
##  [501]  1.700381e+00  6.476523e-01  3.672814e-01 -8.369032e-01  1.353369e+00
##  [506]  1.015470e+00  7.185985e-01  4.829944e-01  5.588759e-01  9.630921e-01
##  [511]  3.376392e-01  2.382026e-01  1.444946e+00  1.174508e+00  1.352157e+00
##  [516]  5.404183e-01  8.796335e-01  1.154160e+00  1.052584e+00  1.702422e+00
##  [521] -3.655239e-01  9.966112e-01 -1.396228e-02  7.050235e-01  4.041291e-01
##  [526]  2.096253e-01  1.337420e+00  1.511345e+00  8.295225e-01  1.409767e+00
##  [531]  1.207037e+00  1.230116e+00  5.157991e-01 -1.606620e-01  1.292045e-01
##  [536]  1.465591e+00  1.080570e+00  7.179811e-01  3.498626e-01  1.650945e+00
##  [541]  1.394286e+00  1.346199e+00  7.120343e-01  4.130187e-01  5.063532e-01
##  [546]  5.763910e-01 -9.063263e-02  1.939406e+00  9.454680e-01  6.745589e-01
##  [551]  1.153400e+00  1.174508e+00 -5.900037e-02  7.755901e-01  1.037945e+00
##  [556]  1.210554e+00  1.136128e+00  8.250558e-01  7.757854e-01  8.521916e-01
##  [561]  1.497452e+00  4.501990e-01  6.894344e-01  1.140821e-01  1.073189e+00
##  [566]  1.310996e+00  1.102848e+00  1.911244e+00  9.928132e-01  4.615295e-01
##  [571] -1.071093e-02  1.035854e+00  6.047453e-01 -2.355508e-02  7.775472e-01
##  [576]  1.000925e+00  1.067071e+00  5.763257e-01  1.426645e+00  9.261839e-01
##  [581]  1.083377e+00  3.868406e-01  1.266985e+00  5.485300e-01  9.537957e-01
##  [586]  1.118984e+00  9.435620e-01  5.534014e-01  1.600988e+00  4.074928e-01
##  [591]  1.111027e+00  1.664978e+00  1.247197e+00  1.419844e+00  1.691357e+00
##  [596]  8.532271e-01  8.829620e-01  1.743176e+00  1.257932e+00  1.338112e+00
##  [601]  1.366447e+00  6.553496e-01 -9.317243e-03  1.873321e-01  6.863797e-02
##  [606]  1.913300e+00  1.436504e+00  8.875095e-01  4.314602e-01  7.124286e-01
##  [611]  6.460094e-01  8.348525e-01  1.482837e+00  1.091167e+00  1.079967e+00
##  [616]  1.260775e+00  7.132736e-01  1.116595e+00  4.461565e-01  1.025130e+00
##  [621]  5.978723e-01  1.115726e+00  1.516664e+00  1.821778e+00  8.632069e-01
##  [626] -5.828197e-01  1.520698e+00  9.345320e-01  1.421205e+00  1.535454e+00
##  [631]  6.262360e-01  1.061397e+00  2.004489e+00  6.352051e-01  1.440245e+00
##  [636]  6.026646e-01  7.679542e-01  1.248194e+00  3.458934e-01  1.455340e+00
##  [641] -2.124085e-01  9.349762e-01  1.575920e+00  1.223262e+00  2.008784e+00
##  [646]  1.378156e+00  8.999982e-01  8.747341e-01  2.065243e+00  4.405605e-01
##  [651]  3.720100e-01  1.705334e+00  1.927290e+00  1.971701e+00  1.433253e+00
##  [656]  1.281645e+00  1.967873e+00  1.740868e+00  8.071908e-02  3.938205e-01
##  [661]  1.207163e+00  1.484899e+00  2.789462e-01  1.698607e+00  8.035819e-01
##  [666]  1.030517e-01  1.075958e+00  1.639139e+00  1.179759e+00  1.611012e+00
##  [671]  1.264492e+00  8.953743e-01  7.746556e-01  7.337445e-01  9.059827e-01
##  [676]  9.940564e-01  3.665538e-01  6.078399e-01  4.120544e-01  9.447172e-01
##  [681]  8.837664e-01  6.610480e-01  2.810628e-01  1.477844e+00  2.289338e-01
##  [686]  1.749664e+00  8.354655e-01  1.010232e+00  1.812287e-01  1.041674e+00
##  [691] -2.886420e-01  2.675035e-01  5.640958e-01  9.904990e-01  9.387365e-01
##  [696]  5.661839e-01  4.432725e-01  1.375520e+00  3.421890e-01  1.000258e+00
##  [701]  1.724718e+00  2.971965e-01  1.310622e+00  8.116034e-01  4.520159e-01
##  [706] -6.396408e-01  9.618854e-02  1.445381e+00  2.098986e+00  1.394639e+00
##  [711]  5.328900e-01  7.673274e-01  8.728758e-01 -3.342471e-02  5.860045e-01
##  [716]  4.537430e-01  7.082981e-01  1.537568e+00  1.716861e-02  1.140256e+00
##  [721]  1.266929e+00  8.200293e-01 -2.206036e-01  4.180101e-01  5.864482e-01
##  [726]  2.330158e+00  9.610023e-01  7.377065e-01 -4.606202e-01  4.381390e-01
##  [731]  3.609012e-01 -3.994145e-02  9.937547e-01  1.128583e-02  3.543883e-01
##  [736]  1.136337e+00  5.799279e-01  8.597063e-01  1.271462e+00  1.295909e+00
##  [741]  2.249626e+00  5.880301e-01  4.814241e-01 -3.709709e-01  6.806497e-01
##  [746]  2.709820e-01  1.013804e+00  9.639950e-01  1.421977e+00  6.290911e-01
##  [751]  1.104142e+00  1.447675e+00  8.582595e-02  1.084666e+00  1.667147e+00
##  [756]  1.971269e+00  7.854464e-01  4.814878e-01 -9.871652e-02  1.283257e+00
##  [761]  4.022030e-01  5.231835e-01  1.877064e-02  5.580352e-01  1.158675e+00
##  [766]  7.929167e-01  8.885520e-01  1.779088e+00  1.450728e+00  3.056582e-02
##  [771]  4.079647e-02  1.172920e+00  6.806831e-01  2.033048e+00  1.270006e+00
##  [776]  3.058859e-01  4.067154e-01  9.454680e-01  1.511716e-03  5.160652e-02
##  [781]  1.207718e+00  1.215098e+00  1.725738e+00  3.337670e-01  1.363019e+00
##  [786]  8.436720e-01  1.137268e+00  4.603762e-01  6.608034e-01  5.885250e-02
##  [791]  1.306297e-01  4.491177e-01  8.141315e-01 -3.118639e-01  1.848675e+00
##  [796]  5.773613e-01  1.185306e+00  6.861534e-01  4.339629e-01  5.314167e-01
##  [801]  1.115594e+00  1.452003e+00  1.342244e+00  6.027739e-01  2.789082e-01
##  [806]  1.078820e+00  1.213678e+00  1.302215e+00  6.541560e-01  1.087655e+00
##  [811] -6.869219e-01  6.894704e-01  4.578868e-01  5.203011e-01  4.807997e-01
##  [816]  1.461710e+00  1.484569e+00  8.186383e-01  3.081047e-01  1.818190e+00
##  [821]  1.441012e+00  1.181533e+00  8.622365e-01  9.320442e-01  4.072889e-01
##  [826]  1.572538e+00  1.962941e+00  1.071272e+00  1.436511e+00  1.028188e+00
##  [831]  1.050247e+00  1.119376e-01  1.064579e+00  5.682112e-05  1.070194e+00
##  [836]  9.105354e-01  1.208570e+00  5.976154e-01  1.243777e+00  1.044359e+00
##  [841]  1.135435e+00  7.697132e-01  2.310314e-01  1.322988e+00  1.825890e+00
##  [846]  6.894344e-01  3.203513e-01  1.116844e+00  1.825890e+00  1.227657e+00
##  [851]  1.121947e+00  2.016966e+00  1.062032e+00  9.320891e-01  4.469733e-01
##  [856]  1.574573e-01  8.774077e-01  7.069994e-01  1.035584e+00  1.309987e+00
##  [861]  6.526392e-01  9.063501e-01  6.137333e-01  1.074579e+00  6.138911e-01
##  [866]  7.700102e-01  1.253292e+00  2.045756e+00  9.977783e-01  2.674243e-01
##  [871]  1.208114e+00  1.378244e+00  6.085133e-01  3.456508e-01  6.950532e-01
##  [876]  9.985791e-01 -1.783769e-01  7.447517e-01  7.593095e-01  4.779721e-01
##  [881]  8.672226e-01  6.766021e-01  5.338457e-01  7.092088e-01  6.250175e-01
##  [886]  1.238554e+00  1.461952e+00  1.241602e+00  3.986952e-01  1.664791e-01
##  [891]  1.447282e+00  6.516148e-01  3.785488e-01  3.396872e-01  1.035689e+00
##  [896]  7.615247e-01  1.237346e+00  7.826392e-01  1.081672e+00  1.627494e+00
##  [901]  8.215545e-02  4.391697e-01  8.188930e-01  3.270705e-01  7.117942e-01
##  [906]  6.379439e-01  1.354226e+00  7.788699e-01  1.194312e-01  1.180598e+00
##  [911]  2.064118e+00  8.714129e-01  1.116735e+00  8.939249e-01  1.015360e+00
##  [916]  6.266814e-01  4.170129e-01  9.613255e-01  1.088740e+00  1.628294e+00
##  [921] -6.461407e-02  9.023850e-01  9.967585e-01  1.369986e+00  7.476078e-01
##  [926]  5.894766e-01  5.924987e-01  7.104541e-01  2.386743e+00  6.636787e-01
##  [931]  1.358230e+00  1.383273e+00  5.407374e-01  1.099253e+00  1.560025e-01
##  [936]  7.196150e-01  1.680026e+00  9.514282e-01  1.141095e+00  1.116248e+00
##  [941]  1.248169e+00  1.559502e+00  1.001752e+00  3.781037e-01  2.569982e-01
##  [946]  7.223169e-01  1.385066e+00  3.584647e-01  1.641371e+00  2.746818e-01
##  [951] -7.287189e-01  1.034677e+00  7.663392e-01  9.730381e-01  1.183945e+00
##  [956]  3.221108e-01  1.081988e+00  4.083712e-01  6.739306e-01  5.819790e-01
##  [961]  1.037556e+00  1.735186e+00  1.370980e+00  1.195649e+00  6.848113e-01
##  [966]  1.784145e-01  5.188170e-01  1.944512e+00  1.832698e+00  1.118587e+00
##  [971]  3.981828e-01  1.134179e+00  1.244507e+00  2.100475e+00  3.781398e-01
##  [976]  9.645438e-01 -2.898742e-01  2.148218e+00  1.220731e+00  1.463919e-01
##  [981]  1.353745e+00  5.141678e-01  1.555754e+00  8.091090e-01  7.826312e-01
##  [986]  5.452383e-01  2.064958e+00  8.374409e-01  1.050442e+00  1.513875e+00
##  [991]  9.355892e-01  4.782605e-01  1.093618e+00  3.428159e-01 -9.513572e-02
##  [996]  1.153552e-01  5.121430e-01  4.243526e-01  1.509363e+00  2.004489e+00
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
##      mean         sd    
##   0.5398299   0.4674424 
##  (0.1478183) (0.1045206)
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
## [1] -0.32450823  2.32837806  0.22809965 -0.07488543 -0.69883542 -0.29150659
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
## [1] -0.0085
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9011492
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
## t1*      4.5 -0.01711712   0.9208002
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 2 3 5 6 8 
## 1 1 2 3 1 2
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
## [1] -0.04
```

```r
se.boot
```

```
## [1] 0.9188494
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

