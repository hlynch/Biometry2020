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
## 0 1 3 5 6 
## 2 2 2 1 3
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
## [1] 0.0138
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
## [1] 2.707685
```

```r
UL.boot
```

```
## [1] 6.319915
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
##    [1] 4.1 5.5 4.4 3.7 3.8 4.8 3.7 3.9 3.7 3.2 4.8 4.1 4.2 4.6 3.0 6.2 4.7 5.2
##   [19] 5.3 4.6 6.3 4.2 3.6 6.1 4.3 4.1 4.9 6.1 4.9 6.2 4.9 3.8 6.0 5.2 3.1 5.3
##   [37] 4.0 4.1 4.0 3.5 3.0 3.0 6.0 5.0 6.0 3.8 6.1 4.2 5.4 5.3 5.2 4.3 5.8 5.0
##   [55] 3.6 3.7 5.2 6.4 4.5 3.3 5.2 3.6 4.2 4.7 4.7 3.3 5.4 4.4 4.7 3.3 4.0 4.1
##   [73] 6.1 4.6 4.1 3.5 5.7 4.3 5.4 4.0 4.1 3.3 3.6 4.1 3.8 4.8 4.1 5.4 4.8 3.4
##   [91] 4.3 4.9 3.4 5.3 3.7 5.6 3.8 2.0 4.9 3.5 3.3 4.6 3.8 2.9 5.2 4.8 5.1 5.1
##  [109] 4.0 3.9 3.6 4.7 3.8 4.9 5.1 4.7 3.6 4.4 4.4 5.4 2.9 6.2 2.9 4.9 2.3 5.1
##  [127] 4.3 4.8 4.9 5.7 4.3 4.4 4.7 3.4 3.8 4.2 5.1 3.9 4.8 5.7 5.8 4.8 5.2 4.7
##  [145] 4.1 5.7 5.3 5.3 5.3 4.8 4.6 4.2 4.8 5.7 3.7 6.2 4.7 5.4 5.2 3.0 4.5 3.9
##  [163] 5.5 5.5 3.8 3.6 5.3 4.1 4.4 4.8 4.0 3.2 2.7 3.7 5.2 3.9 3.9 3.4 5.8 2.2
##  [181] 4.1 3.3 5.9 3.9 4.5 4.9 3.2 4.7 5.0 4.7 4.8 5.3 3.3 4.2 3.5 3.9 3.4 5.3
##  [199] 5.2 4.1 5.7 4.5 2.7 4.0 5.9 4.2 5.3 4.4 5.0 5.0 5.0 4.1 3.9 4.9 5.5 4.4
##  [217] 3.1 3.4 5.9 4.0 2.5 3.5 5.5 3.9 3.3 4.7 4.4 4.3 5.6 4.3 4.9 5.5 4.9 4.9
##  [235] 6.0 3.6 2.9 6.0 4.4 3.7 4.9 4.8 3.5 4.3 2.9 5.4 3.7 2.8 5.6 4.7 4.6 4.0
##  [253] 4.8 4.4 3.9 3.5 4.9 5.5 4.7 5.8 5.3 4.1 4.4 3.8 3.8 3.4 5.4 3.7 2.4 4.5
##  [271] 5.5 4.2 4.1 4.7 3.7 5.8 5.3 4.7 4.5 4.3 4.2 4.2 3.6 4.1 4.0 4.3 5.9 3.5
##  [289] 5.1 5.2 5.0 4.9 3.5 3.5 5.0 4.9 4.3 4.1 3.2 4.8 6.4 5.7 6.6 4.6 5.7 3.5
##  [307] 5.7 4.9 3.8 4.4 5.7 4.0 4.2 2.5 4.3 3.1 4.8 4.7 4.0 4.5 5.1 3.6 5.0 3.9
##  [325] 4.5 6.4 4.5 4.2 5.4 4.9 3.4 5.8 4.6 4.3 6.1 4.7 4.9 4.2 4.3 2.9 3.9 5.1
##  [343] 4.6 5.7 3.6 4.5 5.0 5.3 4.8 5.7 4.8 3.9 5.2 5.4 4.4 4.1 5.0 7.2 3.7 3.1
##  [361] 4.6 4.6 5.1 4.4 4.8 5.9 5.0 6.4 5.0 4.5 5.7 5.1 5.2 4.5 5.9 5.3 5.6 4.2
##  [379] 5.5 5.7 7.2 4.7 4.7 4.6 2.6 5.6 2.8 3.5 4.8 4.0 4.4 4.3 5.1 4.7 5.2 3.6
##  [397] 6.2 4.2 5.4 4.8 4.4 4.9 4.2 4.3 3.8 3.6 4.5 4.9 4.0 4.6 5.1 5.7 5.2 3.4
##  [415] 4.2 4.2 5.3 4.3 3.5 4.4 5.8 3.7 4.9 4.5 4.4 3.3 5.5 4.2 5.0 5.9 3.6 3.3
##  [433] 4.5 5.1 3.5 3.3 5.2 6.1 3.3 2.9 4.0 4.9 5.5 5.6 4.0 6.0 3.2 3.6 5.0 4.1
##  [451] 4.6 4.0 4.0 3.3 4.9 4.7 4.5 5.1 6.2 4.6 4.7 5.9 3.7 4.1 4.4 3.4 4.5 6.0
##  [469] 3.8 5.7 4.4 4.0 5.4 5.1 4.8 7.0 3.3 4.6 3.2 4.0 6.9 3.7 4.9 4.2 4.0 4.5
##  [487] 3.6 4.5 5.9 5.6 3.3 5.5 4.7 4.7 3.6 2.8 4.8 5.0 5.9 5.4 4.3 3.5 4.2 3.7
##  [505] 1.9 4.3 5.7 5.1 5.0 5.9 5.1 4.5 5.0 4.9 4.4 4.5 4.0 5.3 4.5 5.7 3.7 4.4
##  [523] 4.1 5.9 4.5 4.8 4.4 4.0 4.5 5.0 3.3 5.4 4.5 4.1 2.8 5.1 3.3 4.1 5.4 4.8
##  [541] 3.1 4.3 4.7 3.9 4.3 5.6 4.7 6.1 5.1 5.0 2.8 3.9 4.4 4.9 6.3 4.5 5.0 3.2
##  [559] 5.1 4.0 4.7 5.9 4.7 3.4 4.8 4.1 4.5 3.9 3.6 3.9 3.3 5.1 5.3 4.3 5.2 4.3
##  [577] 4.4 4.0 3.9 4.1 3.4 5.1 3.7 4.6 5.1 4.5 4.1 4.3 4.9 4.0 4.6 4.4 5.6 3.9
##  [595] 6.1 3.6 4.9 4.5 4.5 5.5 5.2 6.1 4.1 5.1 4.2 5.5 5.5 5.1 5.4 5.6 4.0 3.6
##  [613] 3.9 4.9 5.1 3.3 3.2 4.3 4.9 4.2 4.8 2.2 3.5 3.6 5.6 3.4 5.5 4.4 4.3 4.9
##  [631] 5.6 5.9 5.6 4.6 3.7 5.7 4.8 4.6 3.6 3.5 4.9 4.1 4.5 4.6 4.2 5.1 3.6 5.7
##  [649] 4.1 4.6 5.8 4.2 3.7 4.7 3.9 4.0 4.5 3.9 3.2 5.4 5.0 4.2 4.7 4.3 5.4 3.9
##  [667] 4.3 4.2 4.0 4.2 2.5 2.6 5.4 5.2 6.5 6.1 4.4 4.4 4.8 5.5 3.7 2.8 5.3 4.2
##  [685] 5.8 4.9 3.5 3.8 5.6 4.9 3.9 5.6 2.9 4.8 5.4 4.7 3.9 5.0 4.2 4.9 5.6 6.6
##  [703] 3.7 5.6 3.6 5.2 5.4 5.9 4.2 5.2 4.8 4.5 4.2 3.9 3.8 4.0 3.7 4.9 3.5 5.4
##  [721] 3.8 3.9 5.2 3.3 5.1 3.9 3.4 4.2 5.4 5.0 4.0 4.0 4.9 6.0 5.2 4.4 3.0 3.6
##  [739] 4.2 5.3 4.2 2.6 3.1 5.3 4.3 3.0 5.0 4.8 4.8 5.1 3.8 4.4 4.9 4.7 4.8 5.2
##  [757] 3.8 5.4 3.9 2.9 5.5 6.5 4.9 4.9 3.7 4.5 2.6 2.7 2.9 3.4 4.8 3.7 6.3 3.6
##  [775] 3.2 4.3 6.5 4.1 5.7 2.8 4.2 4.1 4.8 4.8 2.6 5.0 4.9 3.5 4.4 5.9 5.9 6.0
##  [793] 5.2 2.6 5.6 2.2 6.0 5.3 4.7 3.7 4.7 4.4 3.5 5.3 4.5 4.9 5.7 4.0 4.1 3.6
##  [811] 6.3 5.0 4.3 4.1 5.1 4.4 2.5 5.1 5.4 4.2 6.2 3.5 4.6 3.6 5.2 5.3 4.4 3.9
##  [829] 2.6 4.9 4.2 4.3 4.8 3.8 4.0 3.0 4.7 4.9 3.6 4.9 5.7 2.7 6.3 4.1 5.4 4.5
##  [847] 3.2 4.7 4.0 3.6 3.6 4.4 3.6 3.0 3.9 3.9 3.6 4.8 4.8 3.9 4.1 4.4 6.0 5.9
##  [865] 4.8 6.3 5.2 4.1 5.0 4.0 4.4 3.9 6.4 4.5 4.6 3.0 4.2 5.6 4.8 4.5 4.2 2.7
##  [883] 5.3 5.5 4.5 5.1 4.2 4.0 4.3 3.8 3.4 3.5 3.6 4.4 4.0 4.4 3.0 3.8 4.9 4.3
##  [901] 2.9 2.8 3.8 4.2 4.2 4.0 4.4 4.8 4.3 3.8 3.5 5.3 3.0 3.5 4.4 4.3 4.6 4.4
##  [919] 3.5 3.4 4.9 4.9 4.2 4.1 5.4 4.8 4.7 4.4 5.7 5.9 5.5 4.9 5.0 4.2 6.6 4.2
##  [937] 3.9 4.1 3.7 4.6 4.9 5.4 4.2 4.1 4.9 4.6 5.5 5.7 4.8 5.1 5.5 3.4 3.3 5.1
##  [955] 4.8 4.3 4.0 4.5 4.6 5.0 3.1 4.1 4.4 3.8 4.1 4.5 5.1 4.9 3.3 5.6 6.3 4.9
##  [973] 6.3 4.8 5.1 5.8 5.8 5.9 4.3 4.9 4.1 4.4 4.6 4.1 3.8 3.5 3.9 6.1 4.4 4.6
##  [991] 4.1 2.4 2.3 4.0 3.7 4.1 6.3 6.3 3.7 3.5
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
##    [1] 4.4 5.9 4.7 2.8 4.4 4.5 4.5 5.4 4.6 5.4 3.8 4.2 3.9 5.9 5.2 3.3 4.8 2.6
##   [19] 3.6 5.6 3.8 5.7 3.7 3.7 5.3 5.2 4.0 5.0 3.5 5.0 4.0 4.2 4.0 4.1 6.6 5.8
##   [37] 4.5 3.3 4.6 4.4 3.9 5.1 4.0 4.7 5.0 4.7 4.4 3.0 4.4 4.7 4.7 4.7 3.7 5.7
##   [55] 4.7 6.4 6.3 4.7 4.1 5.2 4.4 6.2 4.4 3.6 5.4 4.8 5.3 4.6 4.7 4.1 3.7 4.3
##   [73] 4.5 4.5 4.4 4.7 4.0 4.1 2.2 5.6 5.0 4.1 6.3 4.3 5.5 5.4 4.7 5.6 4.3 4.4
##   [91] 5.1 5.2 5.1 3.9 4.2 4.4 4.0 4.6 5.2 3.3 5.0 4.9 2.8 6.4 3.9 4.2 5.1 4.7
##  [109] 5.0 3.7 6.6 4.4 4.4 4.4 4.6 5.0 4.5 4.4 4.8 5.8 4.0 3.8 2.8 4.5 5.1 4.6
##  [127] 3.9 5.2 4.5 4.2 3.8 4.7 4.2 5.1 4.9 5.1 4.8 3.9 4.8 2.8 4.6 3.7 2.4 4.1
##  [145] 6.2 5.4 5.3 2.4 3.7 5.0 3.6 4.5 5.7 6.6 4.9 3.6 5.5 5.7 5.5 5.7 4.4 4.0
##  [163] 4.8 4.6 4.4 4.6 4.6 2.6 5.1 3.7 4.6 3.7 4.7 4.6 4.9 3.2 5.0 3.5 3.8 4.2
##  [181] 5.8 4.7 4.0 3.5 4.1 5.4 4.6 3.8 5.3 5.0 3.4 5.0 4.8 5.5 4.9 2.9 3.5 4.6
##  [199] 3.8 2.9 2.1 3.4 3.7 5.2 3.5 4.5 4.7 4.3 6.2 3.7 3.2 4.7 5.9 4.9 4.5 6.2
##  [217] 4.8 3.5 5.7 4.1 4.5 4.6 4.9 5.5 4.8 3.5 4.7 5.1 3.0 3.8 3.7 3.1 3.7 5.0
##  [235] 4.4 5.0 4.2 4.7 3.3 3.7 6.0 4.8 3.3 3.6 4.1 4.0 3.7 2.4 4.6 4.3 3.2 4.1
##  [253] 4.6 5.2 3.6 4.8 4.9 6.0 4.5 4.1 4.1 5.3 5.2 5.1 4.8 3.9 5.1 4.7 4.1 4.5
##  [271] 4.5 4.3 4.5 3.8 4.8 4.7 4.6 4.9 4.6 5.5 5.5 4.6 5.3 3.1 4.1 4.7 4.4 4.9
##  [289] 4.5 5.0 5.3 4.8 4.2 4.3 3.7 5.0 4.8 4.2 4.2 4.3 3.9 4.8 5.3 5.4 5.5 3.4
##  [307] 4.4 3.3 5.9 3.8 3.2 4.2 4.2 3.4 5.4 3.3 4.2 5.0 4.2 3.5 3.5 6.0 4.5 4.9
##  [325] 5.5 4.5 5.9 4.1 3.0 3.8 4.6 4.0 5.6 3.9 5.0 4.8 4.3 3.4 5.4 4.6 3.7 4.4
##  [343] 5.3 4.8 5.2 4.6 4.6 3.4 4.2 4.5 3.5 3.0 5.2 5.0 6.3 4.2 5.1 5.1 6.1 3.8
##  [361] 4.6 4.7 6.2 5.7 3.7 5.1 4.0 2.8 3.2 4.6 3.9 5.0 4.4 4.2 4.2 4.6 3.0 5.0
##  [379] 4.0 3.5 4.3 3.3 3.9 5.3 5.3 3.2 3.7 4.3 5.1 4.7 4.9 4.7 4.5 3.8 2.8 5.3
##  [397] 6.5 5.6 5.7 4.8 3.8 5.4 5.2 5.6 5.6 4.8 4.9 4.7 3.3 3.4 4.7 4.5 5.6 5.6
##  [415] 5.0 5.7 4.9 4.0 3.4 6.6 4.6 5.5 5.2 3.6 4.3 5.8 4.6 4.1 4.6 4.5 4.9 3.6
##  [433] 5.4 4.3 5.2 4.5 6.9 3.5 4.3 3.8 4.8 4.2 2.7 4.7 4.0 5.3 4.8 4.7 4.1 6.2
##  [451] 4.9 4.5 4.3 3.5 4.6 4.6 5.2 4.8 4.3 3.4 4.7 4.7 3.0 4.2 6.0 4.1 4.2 4.8
##  [469] 4.4 4.1 5.9 6.5 5.1 3.4 4.8 4.6 3.8 4.6 5.3 5.9 4.3 4.6 3.2 5.1 6.4 6.5
##  [487] 6.0 5.6 5.3 2.2 5.3 5.8 3.5 4.7 3.7 3.9 3.8 5.8 5.3 6.5 4.2 4.2 5.7 4.1
##  [505] 5.8 5.2 5.4 4.2 3.7 5.1 4.0 3.6 4.6 5.3 5.1 4.9 4.8 4.1 3.3 5.7 3.7 4.1
##  [523] 4.7 3.7 5.1 5.7 3.6 3.9 4.0 4.3 3.7 4.2 6.4 6.9 5.8 4.7 4.3 4.2 4.6 5.2
##  [541] 3.9 6.6 5.0 4.2 4.7 3.7 2.5 4.2 6.9 3.7 4.5 4.5 4.6 4.5 4.0 5.6 3.4 4.5
##  [559] 5.1 5.4 4.1 5.6 4.0 4.7 6.1 5.8 5.5 5.9 2.7 4.8 3.3 3.3 3.5 5.5 3.2 5.0
##  [577] 4.3 6.2 4.5 4.2 5.3 2.5 4.9 4.3 6.1 5.8 5.7 4.5 5.7 4.8 5.9 4.6 3.9 3.3
##  [595] 3.5 5.4 2.9 4.9 4.9 2.7 5.1 4.7 4.7 5.2 4.5 3.3 4.1 4.6 4.2 2.7 4.2 4.2
##  [613] 2.7 5.1 4.9 5.3 4.2 3.7 4.6 4.7 2.2 6.1 4.9 4.5 4.8 5.5 6.1 3.4 5.9 4.3
##  [631] 5.5 5.5 4.3 3.0 4.1 4.5 3.5 3.0 2.5 4.3 6.0 5.2 4.8 6.2 4.5 4.5 4.8 3.5
##  [649] 6.8 5.5 5.2 4.3 3.8 5.3 3.6 5.5 4.3 4.9 4.9 4.4 6.4 6.4 5.3 5.2 5.5 4.9
##  [667] 3.5 4.3 3.8 4.9 5.9 5.7 4.2 4.5 3.0 4.8 4.6 4.8 4.6 4.5 4.3 5.7 4.7 3.2
##  [685] 5.0 3.1 2.8 4.6 3.8 5.1 3.9 4.2 5.1 2.8 3.1 4.6 5.0 3.2 3.0 5.6 3.8 4.4
##  [703] 3.9 4.2 3.1 3.8 4.7 7.3 5.0 4.1 3.8 5.8 6.2 5.4 4.7 5.7 5.1 4.9 5.2 5.2
##  [721] 3.4 4.6 4.0 5.1 4.6 5.4 5.3 5.2 3.3 4.0 4.2 4.6 5.8 3.7 4.0 5.7 5.6 2.9
##  [739] 4.1 5.8 5.3 4.9 5.4 4.3 5.1 2.4 4.6 5.1 4.7 3.6 4.5 3.8 5.6 4.8 5.0 4.4
##  [757] 4.1 3.7 3.4 4.5 3.8 4.0 5.2 4.3 4.2 4.5 4.8 4.4 4.1 4.3 4.9 3.5 4.2 4.5
##  [775] 4.4 4.2 4.0 3.0 5.9 6.2 4.4 4.9 3.7 4.4 4.0 7.0 6.1 5.6 5.1 4.3 4.4 5.6
##  [793] 4.6 4.7 5.6 6.0 4.3 3.5 5.3 3.0 6.0 4.4 4.0 4.4 6.7 4.1 5.3 6.3 5.0 4.5
##  [811] 4.9 4.6 3.0 3.9 3.1 4.7 3.8 5.3 6.1 3.5 4.2 5.1 3.7 4.7 4.8 5.8 4.2 2.5
##  [829] 5.6 3.6 4.7 3.9 4.8 4.8 4.9 4.0 4.3 3.8 4.6 4.3 4.4 4.4 3.7 5.8 3.6 4.7
##  [847] 4.3 5.6 5.3 3.1 3.7 4.9 3.8 4.8 4.2 5.2 3.8 3.8 5.8 3.5 4.3 4.7 4.0 3.2
##  [865] 6.2 4.9 3.6 4.3 3.9 4.8 4.0 4.6 3.8 5.3 4.3 3.3 4.5 5.3 5.2 5.1 4.4 4.7
##  [883] 4.1 3.0 3.0 3.3 5.6 5.5 3.9 4.2 5.6 5.0 4.5 5.8 5.4 4.8 4.1 4.7 5.1 3.0
##  [901] 4.8 5.6 2.7 4.2 3.3 3.7 4.7 5.3 5.0 4.6 5.8 4.4 5.1 4.7 3.1 5.4 5.1 3.6
##  [919] 4.3 4.7 3.2 5.8 4.6 4.8 3.7 3.8 3.4 3.1 3.8 5.1 4.5 5.2 5.7 5.1 2.7 3.1
##  [937] 3.9 5.0 4.0 4.7 4.7 4.4 3.9 3.1 3.6 4.2 6.6 4.6 2.8 3.9 3.6 4.7 5.0 4.3
##  [955] 3.6 4.3 5.6 4.5 3.8 3.7 3.8 4.9 5.0 2.7 4.6 5.5 3.6 4.9 5.7 6.9 3.7 3.9
##  [973] 4.3 2.5 4.8 3.6 4.2 5.2 5.0 3.5 6.8 4.8 3.1 4.9 4.5 3.9 5.0 6.0 5.8 6.5
##  [991] 4.8 3.8 4.0 4.2 5.8 4.8 3.8 6.3 5.1 4.6
## 
## $func.thetastar
## [1] 0.0398
## 
## $jack.boot.val
##  [1]  0.50664820  0.44858757  0.25470588  0.24537815  0.08630137 -0.01820809
##  [7] -0.08181818 -0.27191011 -0.24227642 -0.48389058
## 
## $jack.boot.se
## [1] 0.9220278
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
##    [1] 4.6 6.0 5.5 4.3 4.7 4.9 4.7 4.7 4.8 3.0 3.9 6.4 4.6 5.2 4.3 4.5 3.3 6.1
##   [19] 5.3 4.6 5.9 3.7 3.2 5.7 4.9 3.7 5.4 4.0 5.2 5.1 4.7 4.7 3.7 3.7 5.6 5.0
##   [37] 4.9 5.1 5.5 4.5 4.5 4.0 3.8 5.1 4.1 3.4 4.5 5.9 3.8 4.8 5.7 4.2 4.3 4.6
##   [55] 3.4 3.0 5.4 3.7 4.0 2.6 2.6 4.9 3.8 4.7 2.8 3.0 3.7 6.0 2.6 5.6 4.4 3.5
##   [73] 4.4 6.3 4.1 4.8 2.4 5.7 4.9 5.0 3.3 4.1 4.6 2.5 5.7 5.5 4.0 3.5 3.7 4.6
##   [91] 4.7 5.3 3.9 4.5 5.8 4.8 3.3 5.2 4.3 3.8 3.6 5.6 4.9 6.1 3.3 5.2 4.6 6.5
##  [109] 4.5 4.7 4.0 5.3 5.5 4.1 3.4 4.8 4.7 3.5 3.0 4.3 3.2 5.1 3.7 4.2 5.1 6.6
##  [127] 3.9 6.3 4.7 4.1 5.2 5.2 4.1 4.1 3.9 4.7 5.4 4.3 3.8 4.5 4.7 4.1 4.8 5.1
##  [145] 3.9 5.3 6.2 4.1 3.6 4.8 5.3 4.7 4.8 5.4 4.4 5.2 3.7 4.6 4.1 3.9 4.7 4.8
##  [163] 5.4 4.9 4.2 5.7 3.1 4.7 4.0 4.1 5.1 4.7 5.5 4.1 4.4 4.8 3.8 3.6 6.0 3.7
##  [181] 5.8 6.4 4.0 3.7 6.1 4.0 4.7 5.4 3.9 5.0 3.9 4.9 4.6 4.4 3.5 5.1 4.5 5.1
##  [199] 5.5 6.2 5.6 4.8 5.3 3.9 4.5 3.6 4.6 2.6 3.6 4.7 4.2 5.7 5.0 5.3 4.5 3.8
##  [217] 3.7 5.3 2.4 4.4 5.1 5.2 6.1 3.6 5.1 4.0 5.6 3.5 2.5 4.1 4.9 3.3 2.6 4.8
##  [235] 5.5 4.7 3.9 4.7 4.0 3.9 5.0 3.8 5.5 6.1 4.5 3.8 4.3 2.8 5.7 4.1 5.2 2.9
##  [253] 5.7 4.3 4.1 3.6 4.1 4.2 4.3 4.4 4.1 5.0 5.9 4.5 5.3 4.9 5.3 4.1 5.0 5.5
##  [271] 3.0 4.5 5.4 6.0 4.4 4.2 5.4 3.7 4.2 4.3 4.2 3.7 4.0 3.6 4.0 4.9 4.4 5.0
##  [289] 3.1 5.5 4.8 2.9 6.3 5.2 3.6 4.5 4.8 5.8 5.2 4.3 4.7 6.2 1.6 5.4 4.1 4.4
##  [307] 4.6 4.3 6.0 3.6 3.4 3.5 6.6 3.5 4.2 2.8 4.9 5.4 4.2 4.9 3.5 5.0 6.3 3.9
##  [325] 4.1 3.4 4.2 3.7 6.1 4.1 3.1 3.6 6.5 5.8 5.6 4.3 4.9 4.3 4.0 5.0 3.3 5.0
##  [343] 4.8 5.4 5.7 3.0 5.1 3.3 5.2 4.8 4.2 6.0 4.1 5.2 4.4 5.1 5.6 3.6 5.8 3.4
##  [361] 3.1 6.0 4.7 4.3 6.2 4.8 4.6 3.0 5.9 5.0 4.9 3.9 4.6 5.4 5.3 5.6 2.6 3.1
##  [379] 4.4 3.0 5.5 4.2 6.0 4.0 4.5 6.9 5.1 5.4 5.2 4.2 5.5 4.2 4.3 4.4 4.7 4.8
##  [397] 4.0 4.1 4.2 4.7 5.6 4.9 4.7 5.9 7.0 3.3 4.5 4.4 4.9 5.7 6.1 4.4 4.4 3.5
##  [415] 3.2 3.0 5.4 2.8 3.1 4.6 5.6 5.1 3.9 4.3 5.4 3.7 5.0 5.0 4.5 4.5 6.1 6.0
##  [433] 5.9 3.4 6.7 3.4 4.3 5.0 4.3 5.0 4.4 5.4 4.5 4.7 5.2 5.8 5.7 4.3 4.1 4.9
##  [451] 3.3 5.4 5.8 5.5 6.0 4.9 4.2 3.5 5.6 5.6 4.0 3.8 3.8 4.1 3.9 4.2 4.8 4.4
##  [469] 5.0 5.1 3.7 4.0 4.0 4.2 5.3 3.6 4.2 4.3 4.7 4.9 4.0 5.2 2.7 3.0 6.4 6.1
##  [487] 4.1 4.6 5.2 3.8 3.6 4.9 3.2 3.0 6.4 5.2 3.2 5.7 4.6 2.6 4.0 4.3 5.4 4.2
##  [505] 4.6 4.3 4.1 5.7 3.8 4.4 3.8 4.2 3.7 3.8 4.2 3.2 4.7 4.4 4.3 3.7 5.1 5.6
##  [523] 3.5 5.4 5.2 3.9 5.0 4.5 5.2 4.9 3.9 4.1 4.2 3.8 5.8 3.8 6.2 5.1 5.6 3.6
##  [541] 3.8 5.2 5.3 3.0 3.9 5.0 4.5 3.5 5.5 3.3 5.7 4.0 4.5 3.1 5.0 2.1 4.9 3.8
##  [559] 4.4 2.7 3.2 5.2 4.0 5.4 3.4 4.8 3.3 6.4 5.3 5.5 3.4 3.5 3.3 4.6 4.6 5.1
##  [577] 3.7 5.4 4.2 5.4 4.8 2.9 2.7 3.2 3.8 2.8 4.6 3.3 6.5 2.1 4.9 4.8 3.0 4.9
##  [595] 4.6 3.5 6.4 4.4 5.4 5.8 5.3 4.2 6.2 5.7 6.0 4.0 5.2 5.5 3.6 4.4 5.0 4.0
##  [613] 4.3 2.2 6.5 3.6 6.7 5.2 5.1 5.9 5.2 5.9 4.5 3.5 6.1 6.9 5.3 4.2 4.5 2.9
##  [631] 4.1 3.6 5.8 4.4 5.7 3.6 4.6 4.8 4.4 4.4 4.8 3.7 4.3 5.3 4.5 4.4 4.1 4.5
##  [649] 6.3 3.8 4.1 3.9 5.7 4.6 4.6 6.1 6.2 4.7 5.2 3.7 3.2 5.4 5.3 5.1 3.7 4.3
##  [667] 3.3 4.6 3.5 5.0 4.0 5.6 4.0 3.9 3.9 5.0 5.0 5.1 5.1 5.2 4.3 3.9 4.2 6.4
##  [685] 3.2 4.7 4.4 5.6 4.0 5.3 6.3 4.7 3.8 4.9 3.2 5.0 4.6 3.8 4.0 3.4 5.1 5.3
##  [703] 6.3 4.6 5.4 2.4 5.0 6.4 4.9 5.6 6.0 3.4 5.8 4.8 4.1 4.5 4.5 4.8 4.7 5.1
##  [721] 3.6 4.5 4.3 5.1 4.3 3.7 3.0 3.4 3.6 4.0 3.9 3.0 3.9 3.7 5.0 4.2 5.3 4.7
##  [739] 5.5 4.6 4.0 5.6 3.3 5.0 3.5 4.8 3.9 4.1 2.9 4.4 2.8 4.8 4.8 6.2 3.9 4.2
##  [757] 4.7 3.5 5.1 5.0 5.2 4.6 5.1 6.7 5.9 4.6 2.8 4.7 4.3 5.2 4.2 4.4 4.6 4.6
##  [775] 5.7 4.5 5.2 5.5 5.4 4.6 6.8 4.6 5.2 6.0 5.0 3.9 4.2 4.1 4.0 3.8 3.8 4.8
##  [793] 4.4 2.8 2.9 5.1 4.6 3.9 6.5 7.3 4.4 5.6 4.7 6.0 2.9 3.4 4.8 5.0 3.5 4.6
##  [811] 3.9 4.4 4.0 2.1 4.0 4.1 5.7 4.4 4.3 4.7 2.6 5.4 2.5 4.4 6.3 4.6 5.6 4.2
##  [829] 4.9 3.7 3.8 4.6 4.7 3.8 4.1 4.3 3.9 3.7 4.7 5.1 5.8 2.0 5.2 3.8 4.4 4.8
##  [847] 5.4 3.7 5.2 5.8 4.2 4.8 5.7 3.4 4.9 4.7 4.5 4.2 4.2 4.5 3.1 4.5 3.6 3.2
##  [865] 4.6 4.4 3.7 5.9 4.6 5.2 5.4 4.1 5.6 5.2 3.6 4.3 3.5 5.0 3.5 4.8 5.2 4.9
##  [883] 2.8 3.6 5.1 4.9 6.2 4.1 4.4 4.2 5.4 4.8 4.4 4.5 4.8 4.8 5.1 3.6 4.2 4.5
##  [901] 3.4 4.9 5.7 4.3 3.7 3.4 5.5 5.5 5.1 5.9 5.0 4.1 4.8 5.5 4.4 5.0 2.1 4.9
##  [919] 4.0 5.0 4.9 3.8 4.6 3.7 4.1 4.6 3.5 3.4 3.9 3.1 4.2 4.2 2.5 4.8 6.5 4.4
##  [937] 6.3 3.8 4.7 3.8 5.0 5.8 4.3 4.0 4.3 4.9 4.0 4.5 4.8 4.9 5.6 5.2 4.1 3.3
##  [955] 4.2 5.2 4.4 3.9 6.2 4.8 4.5 5.2 4.8 5.1 3.3 3.3 5.1 5.4 3.4 3.0 5.7 4.9
##  [973] 6.7 4.9 4.0 4.4 6.3 4.5 4.3 3.8 3.4 3.7 4.3 5.0 6.3 4.7 4.7 4.4 4.9 6.4
##  [991] 4.6 4.6 5.1 4.0 4.7 4.7 5.0 4.7 4.5 4.2
## 
## $func.thetastar
## 72% 
## 5.1 
## 
## $jack.boot.val
##  [1] 5.500 5.400 5.392 5.300 5.200 5.100 5.000 4.800 4.600 4.500
## 
## $jack.boot.se
## [1] 0.9908803
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
## [1] 0.2710317
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
##   2.741522   4.424167 
##  (1.159221) (2.052682)
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
## [1]  0.72325600  0.92216282  1.68209278  0.05131418 -0.28623243  1.27703529
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
##    [1]  4.322402e-01  4.554689e-01 -8.468548e-01 -3.219503e-01  2.424820e-01
##    [6]  8.724896e-02  1.031860e-01  3.612718e-01 -5.540836e-02  5.346473e-01
##   [11]  3.369559e-01  3.484142e-01  7.831166e-01  5.257814e-01  3.420148e-01
##   [16]  8.360263e-01  5.010267e-01  6.388938e-01  8.897863e-01  1.299262e+00
##   [21] -1.694509e-01  4.802298e-01 -2.233967e-01  5.054100e-01 -1.180384e-01
##   [26]  5.110042e-01  7.790291e-01 -5.275531e-01  5.139775e-01  2.709924e-01
##   [31] -1.737462e-01 -4.016809e-01  9.774982e-01  3.507102e-01 -7.327271e-01
##   [36] -2.254102e-01  9.472541e-01  4.975028e-01 -2.458151e-01  4.465275e-01
##   [41]  9.132828e-01 -9.602120e-01  3.187009e-01 -2.258309e-01  2.719218e-01
##   [46]  4.401399e-01  1.203611e-01  1.679715e-01  7.282247e-01  1.922390e-01
##   [51]  4.991349e-01  1.518800e-01  1.517657e+00  3.426182e-03 -3.896973e-01
##   [56]  9.103411e-01  8.791292e-01  5.174089e-01  1.106051e+00 -7.450325e-01
##   [61] -3.398780e-01  1.169658e+00  8.541230e-03 -2.459413e-01  1.678781e-01
##   [66]  9.094257e-01  7.133070e-01  5.236164e-02  6.714631e-01  5.442151e-01
##   [71]  5.308781e-02  6.835020e-01  3.156841e-01  5.703708e-01  4.507183e-01
##   [76]  1.723111e+00  1.398301e+00 -3.417086e-01  1.028418e+00  4.065741e-03
##   [81] -2.858964e-02  8.419814e-01  5.503626e-01  6.705272e-01  3.282642e-01
##   [86]  6.939077e-01  6.497454e-01  3.910004e-01  1.313986e-01  1.256367e+00
##   [91]  4.464744e-01  1.018404e+00  8.922451e-01  8.344740e-01  6.122309e-01
##   [96]  6.023198e-01  3.322604e-01 -5.341724e-02  7.852154e-01  6.039403e-01
##  [101]  4.280559e-01 -8.717728e-01  3.495653e-01 -9.191965e-01  1.155773e+00
##  [106]  1.648250e-01 -4.938337e-01  8.604192e-01  1.712914e-01  3.270721e-01
##  [111]  5.169550e-01 -1.063605e-01 -2.348677e-01  4.341801e-02  6.879173e-01
##  [116]  9.610208e-01  5.054100e-01  1.324900e-01  2.627745e-02  2.209350e-01
##  [121]  1.941670e-01 -7.906402e-02  3.612718e-01  6.717333e-01  2.834802e-01
##  [126] -6.586324e-01 -5.344802e-02  9.214554e-02  2.813319e-01 -2.956362e-01
##  [131]  2.128064e-01 -1.250913e-01  1.095437e+00  4.272319e-01  4.234091e-01
##  [136]  2.619513e-01 -1.600203e-02  1.436174e+00  7.034666e-01  1.096630e+00
##  [141] -1.239666e-02 -1.811811e-02  2.037831e-01  5.410175e-01  6.105406e-01
##  [146] -2.406579e-01  5.006908e-01  4.439676e-02  7.443112e-01  6.026800e-01
##  [151]  4.462266e-01 -1.741409e-01  2.215961e-01  2.608910e-01  2.153528e-01
##  [156]  6.231705e-03  7.701136e-01  7.945716e-01 -4.518511e-01  3.879854e-02
##  [161] -5.525245e-01  7.255437e-02  4.746218e-01  3.252159e-01  4.432553e-01
##  [166]  5.098673e-01  4.444302e-02  3.776089e-01 -2.456256e-01  5.737996e-01
##  [171]  3.118391e-01  8.003087e-01  1.872918e-01 -7.081308e-02  4.250254e-01
##  [176]  9.643794e-01  1.508678e-01 -3.328417e-01 -4.976254e-02  4.499440e-01
##  [181]  9.472273e-02  8.715228e-01 -2.873805e-01  1.026483e-01 -6.070547e-01
##  [186]  2.983129e-02  7.939815e-01  4.903228e-01 -4.208128e-02 -4.331543e-01
##  [191] -1.473138e-01  4.287126e-01  5.742647e-01 -2.094017e-01  5.884992e-01
##  [196]  4.497145e-01  3.590624e-01  9.902398e-02  8.959705e-01  6.820725e-01
##  [201]  4.044112e-01  5.973500e-01  1.395246e-01  4.179888e-01  7.506452e-01
##  [206]  3.941111e-01  1.882975e-01 -2.170738e-02  3.739244e-01  5.020297e-01
##  [211]  2.030499e-01  2.335726e-01  1.010724e+00 -1.607023e-01  1.591607e-02
##  [216]  9.969381e-02  9.060550e-01  3.372452e-01  4.620438e-01  7.012133e-01
##  [221] -1.184961e-01  7.862728e-01  5.283527e-01  5.089383e-01  4.749754e-01
##  [226] -1.331026e+00  2.306610e-01 -1.192346e-01  2.028265e-01 -3.450975e-01
##  [231]  6.251421e-01  4.278785e-01 -8.135522e-01  5.077315e-01  7.569950e-02
##  [236]  1.933718e-01  4.352837e-01  7.828020e-01  2.878054e-01  1.129502e+00
##  [241]  8.866713e-02  4.131786e-01  4.736977e-01  3.882502e-01  5.692558e-01
##  [246] -3.480849e-01 -3.157366e-01 -2.334969e-01  2.601065e-01 -1.819694e-01
##  [251]  4.032259e-01 -9.468359e-03  6.898416e-01  4.136058e-01  3.761561e-01
##  [256]  6.834306e-02  6.060989e-01 -2.284694e-01  4.159623e-01  2.817592e-02
##  [261] -3.727146e-01  2.212875e-01  7.444598e-02  6.872434e-01 -5.573170e-02
##  [266]  5.114010e-01 -4.159612e-01 -6.221381e-01  1.708114e-01  1.255433e-01
##  [271]  1.137532e+00 -2.110947e-01  7.848128e-01  4.569847e-01  9.275777e-02
##  [276] -3.669412e-02  2.066804e-01  5.026446e-01  4.040856e-01 -1.891041e-01
##  [281]  2.435359e-01  4.406469e-01  1.775003e-01  4.914489e-01  2.701683e-01
##  [286]  6.616593e-02  8.944265e-01  1.144016e-01  5.105964e-01  7.844161e-03
##  [291]  3.860370e-01 -6.649988e-01  6.402055e-01  5.261887e-01  7.501294e-01
##  [296]  1.159708e-01  5.508299e-01  4.203538e-02  4.774286e-01 -6.739643e-01
##  [301]  1.694643e-01 -2.163583e-01  2.915625e-01  5.801946e-01  4.534582e-01
##  [306] -1.638456e-01  4.702140e-01 -3.322872e-01  1.002404e+00 -8.024659e-01
##  [311]  1.096959e+00  1.240617e-01 -5.948043e-01  4.542882e-01 -5.574310e-01
##  [316]  2.350448e-01  5.231833e-01 -1.571424e+00  5.537464e-01  1.916387e-01
##  [321]  1.616290e+00  2.266410e-01  1.142089e+00  9.242007e-05  9.115044e-01
##  [326]  6.141474e-01  6.189622e-01  6.863336e-01  5.598159e-01  6.503220e-01
##  [331] -2.279836e-01  2.000816e-01  2.391743e-01 -2.061352e-01  3.379934e-01
##  [336] -2.195697e-01  4.481568e-01 -2.482391e-01  4.230318e-01 -4.633950e-02
##  [341]  5.285702e-01 -1.332321e-01 -2.661044e-01  2.000816e-01 -6.520770e-01
##  [346]  4.644054e-01 -3.395739e-01  1.191980e+00 -2.131151e-01  7.795259e-01
##  [351]  1.316941e+00  3.625161e-01  9.351848e-02  4.901763e-01  5.848705e-01
##  [356]  8.480331e-01  4.213946e-01  4.657211e-01  3.536526e-01 -4.611794e-01
##  [361]  6.439814e-03  3.875994e-01  4.964632e-01  4.553276e-01  2.325874e-01
##  [366]  4.866913e-01  3.858483e-01  3.441732e-01  6.860229e-01 -1.099734e+00
##  [371]  1.015190e-01  7.507555e-01  8.069373e-03  5.676182e-02  1.755595e-01
##  [376]  1.689262e-02  8.465031e-01  1.491081e-01  2.898279e-01  5.015830e-01
##  [381]  3.671170e-01  8.017646e-01  4.799811e-01 -4.103019e-01 -6.321396e-01
##  [386] -4.712871e-03  6.662624e-01  7.394289e-01  6.829043e-01  1.717331e-02
##  [391]  5.648465e-01  3.778490e-01  2.707832e-01  1.533734e-01  1.504153e-01
##  [396]  4.413187e-01  4.551624e-01  4.098936e-01  5.114282e-01 -1.113366e-01
##  [401]  5.190568e-02  2.144466e-01 -9.973858e-02  4.338239e-01 -2.814143e-02
##  [406]  1.105966e-01  6.901883e-01 -3.559748e-01  1.978352e-01  3.252159e-01
##  [411]  9.715059e-01 -2.343459e-01  6.721074e-03  4.490626e-01  1.805282e+00
##  [416]  1.439252e-01  6.962611e-02 -8.106929e-01  1.779845e-01 -1.277108e+00
##  [421]  4.185729e-01  4.409251e-01  6.057954e-01 -7.190637e-01  4.234149e-01
##  [426]  3.813839e-01 -1.529172e-01 -4.689541e-03  7.952052e-01  5.106920e-01
##  [431]  5.116578e-01  1.254691e-01  2.196562e-01 -8.824216e-01  4.333444e-01
##  [436]  5.077315e-01  3.211122e-02  8.861304e-01  2.213282e-01 -1.058373e+00
##  [441]  7.875115e-01  7.382169e-01  5.348551e-01  2.012446e-01 -5.524109e-01
##  [446] -2.051435e-01  4.897391e-01  5.877595e-01  8.041859e-01  2.156132e-01
##  [451] -2.247458e-01  6.352864e-01  5.163768e-02  4.684106e-01  5.662999e-01
##  [456] -2.737194e-01  1.129614e+00  5.901966e-01 -4.380389e-01  2.810820e-01
##  [461]  1.780173e-01 -2.865158e-01  2.539985e-01 -6.350764e-01  1.371376e-01
##  [466]  1.147796e-01  1.163313e+00  4.477634e-01  1.765567e-01 -1.569190e-01
##  [471]  2.646625e-01  7.202338e-02 -1.844095e-01  1.001265e+00  4.764526e-01
##  [476] -2.057189e-01  6.977254e-02 -1.448894e-01 -7.336507e-01  1.732606e-01
##  [481]  7.567368e-01  8.201321e-01 -2.717895e-01  5.218697e-01  7.071223e-01
##  [486]  2.834296e-01 -5.882455e-01  7.820447e-01  8.076787e-01  5.624382e-01
##  [491]  6.298576e-01  5.265070e-03  1.092390e+00 -7.796289e-01  2.441199e-01
##  [496]  3.413600e-01  3.804820e-01 -5.010560e-02  5.817822e-01  2.115642e-01
##  [501] -2.765488e-01  8.210013e-02  5.548279e-02  2.841553e-01  1.887517e-01
##  [506]  4.044938e-01 -1.094772e+00 -2.355974e-01  3.794307e-01  1.391891e-01
##  [511]  1.679715e-01  4.195010e-01  4.898372e-01  4.657705e-02 -1.033740e-01
##  [516]  8.827240e-02  1.964691e-01  3.444141e-01  1.977424e-01 -6.497679e-02
##  [521]  8.451090e-01  6.330213e-01  5.085463e-01  7.675737e-01 -1.433285e-01
##  [526]  9.725565e-01 -8.204040e-02  4.921100e-01  3.785680e-01  2.905673e-01
##  [531]  7.448708e-01 -2.144816e-01  1.476196e+00 -7.280925e-02  1.382120e+00
##  [536]  8.949351e-01  5.665103e-01  1.016527e+00  4.383142e-01  1.077951e-01
##  [541]  1.298430e-01  6.124744e-01  6.171243e-01  3.561508e-01  1.000189e-01
##  [546]  1.818535e-01  3.324668e-01  4.429459e-01  2.697090e-01  1.069477e-01
##  [551]  1.198357e-02  9.831179e-01  1.040466e-01 -4.982625e-02  6.828287e-01
##  [556]  1.375732e-01  2.716603e-01  7.991872e-01  9.054974e-02 -1.791476e-01
##  [561] -1.793760e-01 -1.825022e-01  3.838018e-01  1.035156e+00  5.506520e-01
##  [566]  7.103866e-01  5.094011e-01  2.950808e-02  6.050951e-01  1.459764e-01
##  [571]  4.448983e-01  4.574303e-01  7.071223e-01  1.124636e-01  3.833288e-01
##  [576] -6.252377e-01  4.359492e-01  5.818298e-01 -6.829304e-01  1.930343e+00
##  [581]  3.647192e-01 -4.168816e-01  1.390237e-01  8.822184e-01 -1.146857e+00
##  [586]  3.352310e-01  5.275647e-01 -2.253599e+00  5.293397e-01  1.334488e+00
##  [591]  1.452692e+00  3.714099e-01 -4.368781e-01  3.020893e-01  1.506359e-01
##  [596]  4.197353e-01 -2.884756e-01  7.658277e-01  1.765953e-02  1.971202e-01
##  [601]  3.906449e-01  2.477542e-01  8.012850e-01  2.129157e-01  5.084773e-01
##  [606]  8.831069e-01  5.015657e-01  6.096610e-01  9.482321e-01 -5.648360e-01
##  [611]  8.069373e-03  3.556786e-01 -5.145637e-01  8.681320e-01  7.899653e-01
##  [616]  3.393494e-01  1.071793e-01  5.096488e-01  3.374457e-01  1.290612e-01
##  [621]  5.548040e-01  3.829661e-01  4.049377e-01  1.941980e-01 -2.499483e-01
##  [626]  5.460433e-01 -9.191965e-01  2.992449e-01  9.043818e-01 -4.293986e-01
##  [631]  9.124626e-02  4.558761e-02 -1.957019e-01  4.461443e-01  4.659732e-01
##  [636] -2.268889e-01 -4.897867e-01  9.024929e-01 -1.978283e-01  5.898231e-01
##  [641]  1.337565e-01  2.105209e-01  4.683406e-01 -9.153225e-01  4.746218e-01
##  [646]  7.406261e-01  1.104587e+00  4.881118e-01  3.730941e-01  3.935409e-02
##  [651]  7.940905e-01  4.664229e-01  3.817749e-01 -2.404918e-01  7.176523e-01
##  [656]  2.477542e-01  3.990045e-01  2.790502e-01  4.535732e-01 -1.022050e+00
##  [661] -1.307928e-01 -1.637216e-02  4.551624e-01  1.382507e-01  9.068467e-01
##  [666]  1.014005e+00  6.175186e-01 -1.981672e-01 -4.498091e-01  1.461477e-01
##  [671] -5.904050e-01  4.282155e-02  3.354926e-01 -1.995314e-01  4.771026e-01
##  [676]  1.030126e+00 -1.458916e+00  8.790536e-02  7.115735e-01  4.539824e-01
##  [681]  3.687245e-01  9.326583e-02  6.131593e-01  8.745126e-01 -2.558167e-01
##  [686]  4.408843e-01  2.916799e-01  4.284118e-01  3.679410e-01  4.753498e-01
##  [691]  2.296460e-01  2.322553e-04  4.927902e-01  2.676496e-01  4.784728e-01
##  [696]  2.427026e-01  6.534806e-01 -1.397853e+00 -1.945853e-01  5.663150e-01
##  [701]  7.048965e-01  7.032578e-01  1.209597e+00  2.181939e-01  4.608998e-01
##  [706]  5.035982e-01  9.377440e-02  8.402116e-02  9.448241e-01  5.816979e-01
##  [711]  3.435212e-01 -5.447943e-01  4.545318e-01  1.005029e+00  1.204807e-01
##  [716]  6.781605e-01  6.512341e-01  5.584277e-01  7.485342e-01 -1.979408e-01
##  [721]  1.001323e+00  5.265704e-02 -2.456256e-01  3.451101e-01  5.309369e-01
##  [726]  5.750397e-03 -2.587005e-01 -1.221090e-02  9.168863e-02 -5.736537e-01
##  [731] -3.053068e-01  4.404860e-01  3.368504e-01  1.059789e-02  2.005248e+00
##  [736]  4.762928e-01  5.438978e-01  3.997293e-01  9.788704e-01 -4.923710e-01
##  [741]  3.543586e-01  2.163004e-01  6.789873e-01 -7.688199e-01 -7.521053e-01
##  [746]  3.870977e-02  1.292732e-01  1.977424e-01  6.795511e-01  8.969946e-01
##  [751] -1.314414e-01  3.627953e-01  7.614920e-02  2.620087e-01 -3.419795e-01
##  [756] -3.406649e-01 -3.523938e-02  3.147115e-01 -4.613789e-01  7.862728e-01
##  [761]  2.966618e-01  1.083874e-01  3.899004e-01  4.044871e-02 -2.258309e-01
##  [766]  4.587536e-01  1.335201e-01 -3.965181e-01  2.085281e-01  8.757750e-02
##  [771]  1.160617e-01 -4.136668e-01  6.056274e-01  6.994211e-01 -1.685148e-01
##  [776]  2.105667e-01  5.231325e-01  2.004307e-01 -7.941969e-02 -3.339366e-01
##  [781]  4.129232e-01  2.342900e-01  1.178240e+00 -7.906402e-02  9.742157e-01
##  [786]  1.487195e+00  5.932134e-01  6.326263e-01  6.926918e-01  2.767523e-01
##  [791] -2.457698e-01  2.367461e-01 -2.930160e-01 -7.628654e-02  1.450299e+00
##  [796] -1.529843e-01 -6.047808e-01  1.908278e-01  7.214035e-01  6.867394e-02
##  [801]  9.202736e-01  1.754835e-01  3.709685e-01  4.799811e-01  3.317917e-01
##  [806]  4.254838e-01 -2.964429e-01  3.546368e-01  1.921670e-01  3.613876e-01
##  [811]  3.696988e-01  6.482899e-01 -3.253207e-01 -3.665169e-02  9.210861e-01
##  [816]  3.533522e-01  7.311497e-01 -2.887920e-01  1.911318e-01  1.473008e-01
##  [821]  1.003970e+00  1.931759e-01  7.097387e-01  6.994389e-02  5.077315e-01
##  [826]  6.831838e-01  2.661805e-01  6.566937e-01  3.747007e-01  1.597028e-01
##  [831]  3.916253e-01  2.188471e-01  5.333418e-01  1.140522e+00  8.914128e-01
##  [836] -1.266442e-02  6.447079e-01 -4.920230e-02 -3.031870e-01 -9.801956e-02
##  [841]  1.219252e+00  1.594501e+00  3.175620e-01  3.342490e-01  1.079632e-01
##  [846]  2.246702e-01 -1.282866e-01  1.692566e-03 -2.723287e-02 -3.921583e-01
##  [851]  9.354037e-01 -2.686276e-01  7.030288e-01  2.512901e-01  3.597126e-01
##  [856]  7.604551e-01  4.463468e-01  2.104550e-01 -9.622713e-02  2.453607e-01
##  [861]  2.441333e-01 -5.146517e-01 -1.580859e-01  2.258908e-01 -9.451357e-02
##  [866]  8.280094e-01  8.067801e-01  5.363174e-01 -1.855320e-01  3.229048e-02
##  [871] -3.129047e-01  2.568015e-01  9.572118e-01 -1.943772e-01 -1.251741e+00
##  [876]  1.012604e+00 -7.213226e-01  1.359843e-02 -1.973260e-01 -3.979662e-02
##  [881]  4.532194e-01  9.997312e-01  2.136450e-01  4.477234e-01  1.728028e-01
##  [886] -2.174873e-01  5.435312e-02 -1.190920e-01  1.138871e-01  7.621504e-01
##  [891]  8.480331e-01 -1.009730e-01  3.419022e-01  1.161216e+00 -1.640941e-01
##  [896]  2.472567e-01 -5.462853e-01 -3.880209e-01  6.998994e-01  4.453589e-02
##  [901]  7.220811e-01  4.695952e-01  9.619984e-01 -1.253287e-01  6.481210e-01
##  [906]  2.907884e-01  2.594197e-01 -2.860585e-01  1.014085e+00  4.634976e-01
##  [911] -1.003207e+00 -1.855123e-01 -3.379564e-01 -1.555113e-01  9.395254e-02
##  [916]  4.792325e-01 -5.514581e-02  8.694157e-02  2.113985e-01  8.274690e-01
##  [921] -1.466297e-01  1.129012e-01  1.405868e+00  1.780296e-01  1.806387e-01
##  [926]  3.462729e-01  1.014094e+00  1.560138e-01 -1.124944e-01 -1.348251e-01
##  [931]  4.477234e-01  7.931102e-01 -3.240213e-01 -6.145822e-01  2.509584e-01
##  [936]  2.180768e-01  1.093630e+00  8.699278e-01  7.113516e-01  7.288603e-01
##  [941]  4.485530e-01 -5.922838e-01  6.083502e-01  8.861304e-01  9.547846e-01
##  [946]  5.347000e-01  3.508928e-01  1.379521e+00  5.432686e-01  7.523734e-01
##  [951]  6.853721e-01  2.269645e-01  4.400819e-01  1.072865e-01  5.109136e-01
##  [956] -5.368763e-01 -5.087746e-01  4.204847e-01  5.660781e-01  4.948807e-01
##  [961] -1.434921e-01  1.330605e-01  5.438978e-01  5.550683e-01  6.727561e-01
##  [966]  4.333984e-01  1.001501e+00  1.040689e+00  3.366929e-01 -2.772358e-01
##  [971]  6.182349e-01  1.744232e-01 -1.552795e-01  9.944391e-02  7.817614e-01
##  [976]  8.448356e-01  5.448177e-01  1.750758e-01  2.019449e-01  4.369111e-01
##  [981]  7.654563e-01  6.710820e-01  5.672158e-01  8.738299e-01 -2.935285e-01
##  [986] -3.223994e-02  3.450104e-01  3.088737e-01  5.814491e-01 -8.491293e-01
##  [991]  2.227006e-01  1.129502e+00  9.801561e-01  5.370204e-01  8.683268e-01
##  [996] -2.483617e-01  2.297078e-01  4.287126e-01 -2.030151e-01 -7.041065e-01
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
##   0.61967390   0.33921122 
##  (0.10726801) (0.07584673)
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
## [1]  0.3803586 -0.1097397 -1.4537419 -0.3029218 -0.2807832  0.6209265
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
## [1] 0.8994664
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
## t1*      4.5 0.05105105   0.8829585
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 1 3 5 6 7 8 9 
## 1 1 1 1 2 1 2 1
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
## [1] -0.0452
```

```r
se.boot
```

```
## [1] 0.8773065
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

