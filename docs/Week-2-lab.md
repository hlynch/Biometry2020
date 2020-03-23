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
## 0 3 4 5 7 9 
## 2 2 3 1 1 1
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
## [1] -0.0509
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
## [1] 2.650998
```

```r
UL.boot
```

```
## [1] 6.247202
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
##    [1] 4.6 5.8 4.6 3.3 4.1 5.7 5.5 3.2 4.8 4.4 4.2 4.8 4.2 3.6 4.0 2.6 5.5 3.2
##   [19] 4.9 6.7 4.0 3.8 6.4 6.0 5.1 4.7 6.2 3.4 5.0 5.0 4.6 5.6 4.3 4.7 2.5 4.6
##   [37] 3.0 5.3 4.4 4.8 5.3 4.0 5.7 3.4 5.4 5.3 3.4 4.3 4.2 6.0 4.7 5.4 4.4 4.8
##   [55] 4.3 5.4 4.0 4.1 4.9 6.0 5.9 4.7 4.7 4.5 4.3 5.8 4.9 4.0 4.5 5.1 5.0 4.4
##   [73] 5.0 5.0 4.3 5.0 3.8 4.6 2.9 3.8 3.9 5.2 3.3 3.4 6.1 4.2 3.9 5.6 4.7 5.4
##   [91] 4.2 4.7 6.0 4.0 3.9 4.2 5.5 3.8 5.1 5.4 3.9 4.5 3.7 4.2 4.5 4.6 5.5 5.5
##  [109] 4.9 5.0 3.9 4.8 3.2 5.6 3.9 6.5 5.2 4.0 5.6 4.6 4.5 4.8 3.5 3.5 5.3 4.1
##  [127] 4.7 5.8 3.4 4.9 4.5 4.7 4.4 6.7 4.9 3.1 5.2 4.6 5.0 2.4 3.8 5.5 4.9 5.7
##  [145] 5.0 4.9 5.2 3.9 6.4 4.7 2.9 3.4 5.4 3.7 4.5 4.6 5.5 4.2 4.2 5.8 5.1 3.7
##  [163] 2.6 4.7 5.6 5.1 4.4 3.1 6.4 4.3 5.7 3.5 3.7 3.5 4.5 5.1 5.6 3.7 4.7 5.7
##  [181] 4.2 2.9 5.6 4.0 4.8 3.9 4.3 4.6 4.8 5.3 5.3 3.5 3.7 6.0 5.5 5.7 4.8 4.1
##  [199] 3.6 4.2 5.9 4.1 4.0 4.7 4.1 5.0 5.7 5.2 5.6 5.8 3.6 3.6 4.5 3.5 5.3 3.3
##  [217] 5.8 4.3 3.6 4.1 5.4 3.8 4.6 5.1 5.4 5.1 3.8 3.2 5.3 4.6 4.6 4.3 4.4 5.4
##  [235] 5.9 4.1 6.1 5.7 3.6 5.2 5.0 7.4 3.9 4.2 3.8 5.7 4.7 5.1 3.6 5.1 4.2 3.5
##  [253] 5.4 4.0 5.7 3.4 2.8 5.6 4.8 3.6 5.2 2.5 3.9 3.9 5.8 3.8 4.1 5.5 5.3 3.9
##  [271] 4.7 4.2 5.0 3.9 4.9 5.0 5.3 4.6 5.2 5.0 3.2 4.5 4.6 5.1 5.1 4.9 3.7 4.3
##  [289] 5.1 5.7 6.4 3.2 4.0 4.0 6.1 4.0 5.1 5.3 6.0 4.7 3.8 4.1 5.9 3.7 3.6 4.6
##  [307] 4.3 5.6 3.6 2.5 4.2 4.7 3.7 3.0 2.4 5.6 4.3 4.9 4.9 3.6 5.0 4.8 6.0 5.7
##  [325] 4.9 4.1 4.5 4.4 4.0 4.4 4.2 5.2 4.1 3.6 4.6 5.1 6.3 3.5 3.8 4.3 3.3 4.2
##  [343] 4.2 2.8 3.8 4.9 6.5 4.7 5.0 3.8 3.8 3.8 4.3 4.6 3.8 4.1 4.7 4.2 5.0 4.3
##  [361] 3.2 5.5 4.4 5.6 5.0 5.1 5.0 3.4 5.5 3.2 5.2 5.3 3.0 3.8 3.8 5.0 5.0 3.7
##  [379] 5.1 4.9 4.0 4.3 4.8 3.5 4.3 5.5 4.3 4.1 5.9 4.9 4.4 4.5 5.2 3.8 4.7 5.9
##  [397] 5.3 5.9 5.3 5.7 5.2 3.0 4.1 4.4 5.3 3.5 5.3 4.0 3.1 4.2 4.4 4.9 5.1 5.0
##  [415] 3.5 4.2 4.0 5.5 5.0 6.1 3.8 3.2 4.5 4.1 5.8 4.3 3.4 5.6 6.0 5.0 6.2 3.1
##  [433] 6.4 4.8 4.2 3.1 3.1 4.9 5.2 3.5 4.2 4.5 5.5 3.8 4.5 4.1 4.5 4.2 4.6 4.1
##  [451] 3.8 4.8 5.1 5.0 5.2 6.1 4.5 5.6 4.9 2.9 3.6 4.6 5.0 3.3 4.7 5.5 5.9 4.5
##  [469] 3.7 4.3 5.5 3.7 4.5 3.0 6.1 5.5 4.5 4.9 5.0 4.1 3.4 4.6 2.4 4.0 6.2 5.4
##  [487] 5.2 5.9 6.2 3.8 5.0 4.0 5.3 3.5 3.6 5.6 4.8 3.2 3.9 5.3 5.0 4.6 4.2 3.1
##  [505] 4.3 4.7 4.7 5.0 3.1 4.0 4.4 3.8 5.5 3.6 3.3 3.4 3.6 5.1 5.8 3.6 5.1 3.1
##  [523] 5.1 4.2 4.3 4.0 3.3 4.6 4.5 6.6 3.7 4.9 3.2 5.2 4.1 3.8 3.3 4.0 6.1 5.3
##  [541] 5.2 5.0 3.4 3.2 3.8 5.1 5.5 4.6 3.9 5.0 3.3 4.4 3.4 5.9 2.5 5.0 3.3 5.2
##  [559] 4.6 3.8 4.5 3.9 4.9 4.3 4.0 4.2 3.9 6.2 4.8 4.7 4.8 5.6 4.9 4.2 4.1 4.6
##  [577] 5.2 3.0 6.1 5.2 4.4 7.0 4.1 4.3 6.3 5.0 5.2 4.1 7.4 6.9 4.8 3.4 3.7 4.9
##  [595] 5.0 5.9 3.3 4.6 4.9 3.8 5.0 4.6 5.4 3.7 3.8 5.0 3.7 4.7 4.3 3.2 4.7 4.5
##  [613] 4.5 4.1 3.6 3.4 4.4 4.7 4.3 4.2 4.9 4.9 5.7 4.4 5.2 6.0 4.8 4.6 2.3 5.3
##  [631] 4.3 6.8 3.7 4.7 5.0 3.8 5.4 4.5 5.0 4.1 4.9 4.1 5.4 4.4 5.4 4.4 4.4 3.3
##  [649] 4.0 5.3 4.1 4.5 6.1 4.0 3.5 3.9 5.0 4.7 4.4 3.3 3.9 4.4 4.0 5.4 4.3 3.5
##  [667] 4.3 3.9 4.9 2.4 4.7 4.1 4.0 2.7 4.9 4.6 2.8 5.0 6.5 4.0 4.2 4.8 4.0 3.0
##  [685] 2.9 3.8 3.7 4.6 4.6 4.6 3.0 5.3 4.8 5.0 4.0 4.8 4.3 6.0 5.6 3.7 5.6 4.4
##  [703] 4.0 3.1 3.2 3.6 4.8 4.6 4.4 4.5 4.7 3.9 5.1 5.3 5.0 6.2 4.2 3.8 3.7 3.2
##  [721] 5.9 6.0 3.9 4.9 3.7 4.2 5.0 4.7 4.4 4.8 5.6 3.6 5.3 5.2 2.9 3.6 5.5 3.9
##  [739] 5.6 3.4 4.7 3.7 6.7 4.2 4.5 3.6 3.4 5.4 4.5 5.0 4.5 4.8 3.6 4.5 5.5 5.5
##  [757] 2.8 5.0 5.0 3.6 3.2 3.6 3.3 5.0 3.8 3.4 4.6 4.1 3.7 6.1 5.1 5.8 3.5 5.1
##  [775] 4.5 3.6 4.5 4.6 3.3 4.5 3.3 4.7 4.8 4.2 4.5 5.9 4.8 4.5 5.4 3.0 3.7 5.1
##  [793] 3.9 3.9 3.7 4.5 3.7 4.0 4.9 5.6 3.2 3.3 6.2 4.2 5.5 3.7 4.3 4.3 2.6 4.7
##  [811] 4.3 3.9 4.7 4.8 4.1 3.7 4.1 2.3 6.5 5.7 5.6 4.6 4.3 4.6 5.0 4.9 4.6 3.6
##  [829] 4.6 5.5 3.4 4.7 3.0 3.6 4.3 5.7 3.8 4.4 3.5 3.7 3.2 5.2 5.2 4.2 5.9 4.1
##  [847] 2.8 4.8 2.8 3.2 5.7 7.1 4.9 3.8 3.8 4.7 3.4 5.9 3.6 4.9 4.3 5.0 2.9 5.2
##  [865] 4.0 5.0 3.7 5.2 4.1 5.4 4.6 6.1 5.2 5.6 4.6 5.5 4.9 3.5 6.0 5.3 6.5 3.9
##  [883] 6.8 4.6 4.6 3.8 5.7 4.4 3.7 4.6 4.5 4.0 5.2 4.7 4.2 2.6 4.3 3.7 3.3 2.5
##  [901] 3.3 5.5 5.9 4.2 4.6 3.7 4.2 4.6 6.2 5.1 3.6 5.3 4.8 5.4 4.8 3.7 4.4 5.0
##  [919] 3.4 3.8 4.8 4.1 3.9 4.1 4.1 2.7 5.3 5.5 3.7 3.8 7.2 4.5 4.2 5.4 3.9 3.7
##  [937] 4.7 4.5 3.6 3.9 6.3 4.0 4.3 5.0 3.3 4.1 3.4 3.6 6.1 4.8 3.2 5.5 5.0 4.9
##  [955] 4.1 4.7 4.3 4.8 2.6 3.5 6.1 4.8 5.2 5.7 6.6 4.7 4.4 4.3 4.0 4.7 5.2 4.2
##  [973] 3.8 3.8 5.9 4.9 4.5 4.6 4.5 3.6 4.9 2.0 5.3 4.6 5.4 3.9 4.7 5.1 4.9 4.7
##  [991] 4.7 5.8 4.9 4.0 4.2 6.6 4.9 4.0 4.8 4.5
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
## 2.8975 6.3000
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
##    [1] 4.4 5.4 4.0 3.9 3.6 4.4 5.1 5.3 4.7 3.8 4.7 3.1 6.1 2.1 4.1 5.5 3.8 4.2
##   [19] 6.0 4.9 3.1 3.9 3.9 3.5 3.0 3.0 4.2 4.5 3.7 4.5 3.8 3.6 5.3 4.4 5.4 5.6
##   [37] 4.9 2.2 6.1 3.4 3.4 4.8 4.8 5.0 3.2 3.3 5.4 6.2 6.5 4.4 4.8 4.5 5.3 5.0
##   [55] 3.5 4.9 3.0 5.4 4.8 6.1 3.2 4.3 2.9 4.0 3.9 4.3 5.0 5.0 5.3 3.8 3.2 3.7
##   [73] 5.2 4.0 5.5 5.1 3.8 5.6 4.3 6.2 4.0 4.9 4.3 4.4 5.6 4.2 3.3 3.8 4.2 4.3
##   [91] 3.4 3.7 5.4 4.5 5.8 4.9 4.2 4.7 2.6 4.7 4.1 5.4 4.5 5.5 4.8 4.4 5.2 3.6
##  [109] 4.1 6.6 4.5 3.9 2.6 4.1 4.5 4.1 5.1 5.1 4.1 4.9 5.9 4.0 4.0 4.4 5.8 3.4
##  [127] 5.3 5.2 4.2 5.3 5.7 4.1 6.3 2.4 4.4 5.7 4.3 4.4 6.5 4.7 5.3 4.4 5.3 6.0
##  [145] 4.4 4.9 5.3 3.7 4.5 4.9 4.6 5.2 4.2 3.7 3.0 2.8 4.6 3.4 3.3 3.8 4.1 5.2
##  [163] 4.1 2.7 3.2 3.9 2.9 5.2 4.3 4.8 4.5 6.5 6.1 3.9 5.8 5.8 3.1 5.6 5.7 4.6
##  [181] 4.0 3.9 4.1 3.3 4.6 4.5 4.5 3.8 4.2 5.7 4.8 2.5 3.5 5.8 4.1 3.6 2.9 3.7
##  [199] 5.0 2.9 4.1 4.7 3.6 3.7 5.1 5.5 3.2 3.0 4.2 4.2 4.9 4.9 4.0 4.6 4.5 4.4
##  [217] 3.7 4.7 3.1 3.9 5.4 5.4 5.9 4.0 5.3 3.9 4.6 2.4 5.7 3.4 3.5 4.8 3.7 4.0
##  [235] 2.9 5.2 4.1 4.7 4.8 3.6 4.7 5.6 2.8 4.4 3.2 4.7 3.9 5.2 4.3 4.2 3.4 3.9
##  [253] 2.9 5.2 6.4 5.3 4.1 5.3 4.9 2.3 5.7 4.6 5.4 4.0 4.2 4.5 5.7 5.6 3.5 5.2
##  [271] 4.4 4.7 5.5 5.4 5.2 4.9 4.8 4.3 5.0 4.3 4.4 3.8 3.9 4.7 4.8 5.3 3.0 5.3
##  [289] 4.2 4.3 5.3 4.9 4.9 4.9 3.8 4.7 4.6 4.4 4.4 4.2 4.7 4.7 3.8 2.6 4.7 4.4
##  [307] 3.9 4.3 4.9 3.8 7.4 4.0 6.1 5.1 4.5 5.1 4.7 4.8 6.3 4.4 4.4 5.9 3.6 4.6
##  [325] 4.6 5.9 5.5 5.7 4.0 4.4 4.2 2.9 4.2 5.9 3.5 4.9 5.6 4.8 3.7 3.8 3.8 4.7
##  [343] 5.5 5.1 5.4 3.7 5.1 3.7 3.9 5.4 4.9 3.7 4.4 5.0 3.0 4.5 5.5 4.5 4.8 3.8
##  [361] 4.5 5.3 4.5 4.7 4.7 5.8 5.2 4.4 4.7 5.0 2.6 5.3 4.3 4.4 4.9 4.2 5.3 3.4
##  [379] 3.6 3.1 5.1 4.6 5.6 5.1 4.3 4.5 4.5 4.8 4.2 4.2 4.9 4.3 4.5 3.8 3.2 4.2
##  [397] 4.3 3.4 4.3 4.1 4.4 4.1 3.0 4.7 4.1 5.7 5.2 5.1 4.4 4.6 3.5 4.8 4.3 5.1
##  [415] 4.2 4.5 4.3 4.0 4.1 5.5 5.0 4.4 3.8 5.4 3.5 4.9 4.8 5.6 4.2 4.0 4.8 4.5
##  [433] 3.2 3.6 3.2 6.1 4.7 4.3 3.5 2.7 4.7 4.7 3.0 3.9 3.7 3.6 4.6 5.0 6.0 5.1
##  [451] 5.2 4.7 3.6 3.9 4.1 3.8 4.5 4.3 3.3 6.2 2.3 3.5 5.6 5.4 5.8 4.1 5.4 3.4
##  [469] 4.3 4.6 4.4 3.7 3.5 6.5 5.3 6.1 4.1 3.9 5.1 7.5 4.8 4.8 5.4 4.1 3.8 4.2
##  [487] 4.7 5.3 4.4 4.5 3.0 3.0 5.0 4.1 4.6 3.8 4.0 3.8 3.4 4.7 4.3 5.2 5.9 3.5
##  [505] 3.4 3.5 3.5 5.2 5.1 5.0 5.8 4.7 4.1 4.8 3.5 3.2 5.1 5.2 6.5 5.4 2.8 5.2
##  [523] 4.4 4.4 4.6 4.6 5.7 5.4 4.5 4.1 3.5 5.6 4.3 4.4 5.1 3.0 4.0 5.4 4.3 4.8
##  [541] 5.0 3.6 3.4 5.0 4.4 5.0 4.2 5.0 2.7 3.9 3.9 4.5 4.9 4.0 4.6 4.6 3.9 4.4
##  [559] 5.2 6.1 4.6 3.7 3.9 5.5 4.9 6.6 3.0 4.6 5.0 5.8 4.0 4.1 3.9 5.2 4.0 4.5
##  [577] 5.0 5.0 4.6 3.8 6.1 5.3 5.2 4.6 4.6 5.1 3.1 5.9 4.3 4.0 4.1 2.8 6.5 3.7
##  [595] 3.0 5.0 3.8 3.4 5.1 3.7 4.4 4.9 3.6 3.6 6.1 3.9 2.7 2.1 6.4 2.1 3.4 3.8
##  [613] 5.8 3.0 4.6 5.9 4.2 3.8 3.0 5.1 6.3 4.1 4.7 4.4 5.7 5.0 4.3 5.4 4.0 5.8
##  [631] 3.7 3.0 5.5 4.1 5.0 5.5 5.0 3.6 2.3 4.8 4.6 5.4 5.8 4.1 4.8 4.5 4.4 5.9
##  [649] 4.9 3.6 4.8 3.9 5.7 5.7 4.5 4.1 5.0 6.3 3.8 5.8 5.9 2.9 2.6 5.6 3.8 5.3
##  [667] 4.0 3.1 4.8 4.2 3.3 4.1 4.6 3.9 5.1 4.0 3.9 5.3 4.5 5.2 4.5 6.4 5.7 4.1
##  [685] 4.8 3.1 4.0 4.4 4.8 6.3 4.1 4.3 4.9 3.8 4.3 5.4 4.2 3.2 5.0 4.9 3.8 4.4
##  [703] 3.1 5.0 4.2 3.8 4.5 3.0 3.3 4.0 4.5 5.8 3.2 4.6 3.7 4.8 5.1 5.1 6.0 5.5
##  [721] 4.8 4.8 3.8 4.8 5.1 4.5 5.2 4.3 5.6 3.9 3.9 2.6 6.1 5.4 4.7 4.2 5.4 3.9
##  [739] 6.1 4.6 4.0 5.2 5.8 5.5 4.3 5.1 3.7 6.0 4.5 3.8 3.4 4.2 4.1 4.0 2.9 5.5
##  [757] 3.9 5.2 4.9 5.5 3.6 3.1 5.0 3.4 3.6 5.4 3.0 5.9 6.6 4.9 3.5 4.3 2.9 4.6
##  [775] 4.5 5.0 3.5 4.4 6.4 5.5 3.8 4.2 5.7 4.0 4.4 3.5 4.6 4.4 4.8 4.4 4.4 3.9
##  [793] 5.6 4.2 4.8 3.8 4.9 4.8 5.4 5.2 4.3 4.8 5.2 4.0 3.0 3.6 4.6 4.3 3.5 5.1
##  [811] 5.6 6.8 5.2 5.7 4.0 4.7 6.0 2.5 4.4 2.7 4.2 3.0 4.2 5.1 3.2 5.9 5.8 4.4
##  [829] 5.2 6.1 3.6 4.8 5.1 4.9 4.5 5.1 5.7 4.3 4.6 5.5 3.8 4.0 2.2 3.7 4.1 4.2
##  [847] 3.8 5.2 4.0 4.0 3.9 5.0 4.3 4.3 5.8 2.6 4.5 4.4 4.8 4.1 3.6 4.2 2.2 4.9
##  [865] 6.3 2.4 4.6 5.7 4.8 4.1 4.5 4.1 5.5 4.3 5.3 5.3 3.9 5.3 5.5 2.7 4.2 5.2
##  [883] 3.6 5.6 4.6 5.2 4.5 5.3 2.6 4.7 4.3 3.6 4.8 5.1 2.7 4.6 4.4 5.1 2.7 4.4
##  [901] 3.3 4.6 5.6 3.9 3.6 4.2 5.3 4.8 5.7 4.7 4.3 3.5 2.2 5.8 5.0 4.2 4.8 4.0
##  [919] 5.4 4.1 4.2 5.2 5.3 3.5 3.9 5.6 3.8 5.6 5.0 5.3 5.5 3.3 5.1 4.1 4.6 2.9
##  [937] 2.2 4.4 2.0 4.1 4.8 5.3 2.5 4.6 5.0 4.5 4.4 6.9 4.1 6.0 5.5 5.1 3.8 5.4
##  [955] 3.9 3.6 3.4 3.1 4.3 4.3 3.5 3.4 5.0 3.7 4.2 5.6 4.0 2.8 4.0 3.9 4.6 6.0
##  [973] 3.8 5.7 4.4 5.6 5.1 3.4 5.7 2.6 4.2 3.8 5.7 3.0 5.5 4.0 5.9 5.3 3.4 4.9
##  [991] 6.0 4.9 5.7 4.4 2.9 4.2 4.5 4.6 4.6 5.6
## 
## $func.thetastar
## [1] -0.024
## 
## $jack.boot.val
##  [1]  0.493530997  0.414029851  0.162039660  0.157267442  0.004011461
##  [6] -0.147297297 -0.230894309 -0.289602446 -0.376323120 -0.540462428
## 
## $jack.boot.se
## [1] 0.9702977
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
##    [1] 3.9 4.7 3.1 3.4 4.7 4.3 3.8 4.1 5.3 4.2 5.5 3.2 4.9 4.2 4.0 3.5 3.9 3.5
##   [19] 3.0 5.6 5.1 4.6 3.7 4.1 4.4 4.3 5.3 4.9 3.8 5.2 5.4 3.8 4.3 4.2 4.8 4.8
##   [37] 4.5 5.1 4.2 3.5 5.6 4.9 6.5 4.6 4.4 5.2 3.1 4.7 4.8 4.3 4.1 4.1 6.1 4.3
##   [55] 5.4 4.7 2.9 3.8 5.1 4.0 2.5 4.7 4.6 4.9 4.2 2.7 2.8 4.3 5.4 4.7 3.1 4.9
##   [73] 4.8 5.1 5.3 3.5 5.8 4.9 4.4 4.7 4.6 5.5 4.5 3.6 4.3 5.5 6.7 5.0 5.0 4.7
##   [91] 4.8 5.1 5.6 4.3 3.4 5.1 5.2 5.1 5.1 4.9 3.2 4.7 2.3 4.3 5.5 3.2 4.1 6.3
##  [109] 3.1 5.6 4.2 3.8 4.6 4.5 4.1 5.0 3.6 4.4 3.8 4.1 4.2 5.5 2.9 5.1 6.6 4.9
##  [127] 3.7 5.2 4.7 4.6 5.1 3.5 5.4 2.9 4.6 2.7 4.0 2.5 3.9 5.4 4.5 5.2 4.1 5.4
##  [145] 4.7 4.9 5.7 2.7 3.6 5.0 5.0 6.1 2.4 4.4 3.4 4.8 4.2 4.0 3.7 3.4 5.8 5.3
##  [163] 6.2 4.6 3.4 6.2 4.3 4.8 4.6 5.4 5.2 4.5 4.9 5.0 4.4 5.2 3.4 5.5 4.5 5.4
##  [181] 4.4 2.8 3.2 4.8 1.7 4.2 5.8 5.2 4.5 5.1 3.9 4.3 5.8 5.3 4.9 3.4 7.0 3.4
##  [199] 5.1 4.6 4.9 4.8 6.2 6.3 3.2 3.5 4.7 3.6 4.9 5.1 3.4 5.5 5.9 4.0 5.0 3.4
##  [217] 4.4 5.2 3.8 6.5 3.7 3.4 3.2 4.4 2.8 4.9 4.4 3.7 5.4 5.0 3.1 4.0 4.7 4.1
##  [235] 4.7 5.0 4.8 3.8 5.7 3.9 4.1 3.6 3.5 3.6 4.7 3.8 5.1 4.9 5.7 3.6 5.0 3.3
##  [253] 5.3 4.9 5.1 4.0 4.3 2.9 4.1 3.7 3.7 5.6 6.6 4.9 6.3 3.9 4.7 5.3 5.5 4.3
##  [271] 5.9 4.5 3.6 2.7 4.4 5.5 5.3 4.5 3.7 4.3 3.3 4.5 5.0 4.3 4.4 3.9 5.4 3.1
##  [289] 3.6 2.9 3.6 5.8 4.0 6.2 5.0 4.7 5.7 3.9 4.4 6.6 4.9 5.6 5.8 3.5 5.4 4.8
##  [307] 3.2 3.6 4.4 3.3 4.4 5.3 3.1 4.1 5.9 3.7 5.2 4.7 4.1 5.4 6.3 4.5 5.2 3.4
##  [325] 5.1 5.8 5.3 4.7 4.6 3.6 3.3 4.8 5.0 3.9 6.4 5.4 3.4 3.7 3.3 3.8 4.5 4.9
##  [343] 3.5 3.6 3.2 4.5 5.0 4.2 4.6 5.6 5.0 5.6 4.0 4.9 3.1 3.6 3.2 3.3 5.4 4.0
##  [361] 4.3 3.2 4.5 5.4 3.5 4.0 4.4 5.6 6.1 4.0 6.3 3.5 6.0 4.2 5.3 3.4 3.6 6.3
##  [379] 5.6 4.6 4.2 2.8 4.9 4.8 5.0 4.6 5.6 4.9 4.8 4.3 4.1 3.6 3.9 4.1 4.3 4.5
##  [397] 4.8 4.0 3.6 4.0 4.8 2.8 5.3 4.8 2.7 3.3 3.5 3.6 3.2 5.6 6.7 4.8 4.1 3.5
##  [415] 5.9 3.3 4.6 4.5 4.6 5.5 3.6 4.1 3.9 4.4 3.6 3.0 4.3 5.6 5.4 3.4 4.6 3.5
##  [433] 4.5 2.9 4.8 4.0 5.2 3.6 2.0 4.2 5.2 4.9 5.6 5.3 4.5 4.5 5.3 5.5 2.8 6.1
##  [451] 4.3 4.3 3.1 5.1 5.1 3.5 4.6 4.2 4.2 5.0 4.3 3.5 4.3 3.6 5.1 5.0 3.9 5.6
##  [469] 3.5 5.0 3.5 4.7 6.1 3.5 4.4 4.3 5.0 5.4 4.0 4.2 4.7 4.8 5.0 4.2 5.3 5.2
##  [487] 5.5 5.6 3.2 5.3 2.5 4.4 3.7 6.7 4.1 5.1 4.9 3.1 6.5 4.1 3.8 3.8 4.6 4.6
##  [505] 5.1 5.0 2.8 5.6 4.1 4.3 4.4 4.8 5.0 6.0 3.9 4.1 5.1 4.0 3.1 3.7 4.7 6.2
##  [523] 5.2 3.8 4.4 4.7 4.3 3.0 4.8 5.1 3.4 4.2 5.2 5.6 6.2 6.4 4.1 3.7 4.6 3.5
##  [541] 3.2 2.8 4.3 5.6 4.0 4.0 4.3 4.8 2.6 4.0 3.0 3.9 3.5 3.9 4.4 5.0 4.8 5.2
##  [559] 4.2 5.7 5.9 4.6 4.5 6.0 4.4 5.3 4.7 3.9 4.1 5.2 4.0 3.5 4.0 5.2 4.1 5.0
##  [577] 4.3 4.2 6.4 3.3 2.6 2.2 5.0 6.2 6.3 7.5 4.6 3.2 5.3 5.3 4.6 3.4 5.6 5.9
##  [595] 4.2 5.7 4.4 5.8 3.1 5.5 5.7 4.1 3.5 4.8 4.3 3.2 4.6 5.8 4.6 4.0 5.2 5.9
##  [613] 5.2 4.6 4.2 6.4 4.9 4.8 3.6 5.5 4.6 4.8 4.7 3.4 5.3 3.8 5.6 2.5 4.0 3.9
##  [631] 5.5 5.0 4.0 4.3 4.8 5.6 4.8 5.8 4.1 4.4 2.5 5.6 4.5 4.6 4.0 4.6 2.7 4.2
##  [649] 3.4 4.3 4.9 5.0 4.1 4.7 4.4 5.0 3.7 5.0 4.7 5.1 5.6 3.4 6.6 4.9 3.5 5.8
##  [667] 4.9 2.7 4.6 4.4 5.0 6.0 5.9 4.6 5.1 3.8 4.4 5.1 3.5 3.6 4.6 3.0 4.4 3.2
##  [685] 3.5 3.9 4.7 3.3 4.6 5.4 3.5 5.7 5.8 5.9 5.8 2.9 3.9 4.3 5.4 4.4 4.6 5.1
##  [703] 5.6 3.0 5.7 4.0 3.6 4.4 3.7 4.5 4.7 4.7 5.9 4.5 3.4 3.6 4.3 4.4 3.8 4.5
##  [721] 4.2 5.5 3.3 5.0 4.5 5.3 5.7 4.7 5.4 4.3 3.8 4.2 4.6 1.9 4.7 3.7 4.7 5.6
##  [739] 4.2 3.1 4.6 2.8 5.4 3.2 4.0 5.1 4.9 3.7 4.1 4.9 5.3 5.8 4.5 3.6 5.6 5.0
##  [757] 4.0 4.6 3.5 5.7 4.3 5.3 4.1 5.5 3.1 5.4 5.1 3.6 3.5 3.6 4.7 5.9 5.1 4.1
##  [775] 3.8 3.4 3.9 4.4 4.5 4.6 5.5 2.7 5.2 3.8 2.4 5.5 4.5 3.6 4.3 4.9 4.6 4.0
##  [793] 4.3 4.2 2.8 5.1 5.3 3.6 3.9 5.1 2.8 4.3 5.3 4.3 4.6 4.9 3.7 7.1 3.1 4.2
##  [811] 4.6 3.2 3.4 4.3 3.2 3.9 4.1 2.9 4.3 4.6 4.7 3.0 3.1 4.7 5.1 4.4 4.5 4.2
##  [829] 2.4 3.5 4.8 3.1 4.9 5.8 3.8 4.9 5.5 4.1 4.5 4.7 4.8 5.9 3.0 3.5 4.0 5.4
##  [847] 4.7 4.6 4.1 5.1 4.7 4.3 3.6 4.7 4.7 5.0 4.7 4.9 3.7 3.2 4.2 5.6 5.1 5.2
##  [865] 5.9 3.3 3.2 4.8 2.9 3.7 5.1 5.4 4.5 3.8 4.1 4.6 4.7 4.4 3.8 5.3 5.8 5.3
##  [883] 3.1 2.3 5.9 4.8 5.5 5.2 5.1 5.0 5.6 4.9 3.1 5.5 3.5 3.1 3.7 3.7 3.9 1.8
##  [901] 5.2 4.7 3.0 3.9 4.5 4.9 5.1 6.2 3.9 4.0 4.0 4.3 4.5 4.0 3.8 4.6 3.7 5.0
##  [919] 3.8 4.0 4.4 4.6 4.4 5.2 5.7 2.6 3.1 2.2 4.1 3.5 4.4 2.9 5.7 4.0 4.4 5.3
##  [937] 6.5 4.8 3.5 4.8 4.7 3.7 4.2 5.6 4.8 5.2 4.8 6.6 5.6 3.8 4.3 4.4 4.8 5.7
##  [955] 3.8 3.2 3.1 4.3 3.6 5.8 3.7 5.3 3.5 4.6 5.8 5.0 4.7 4.7 4.9 6.1 4.9 5.5
##  [973] 3.9 4.6 2.9 4.3 5.8 4.4 5.7 5.4 3.8 3.8 4.7 4.1 5.1 5.2 4.4 2.9 4.4 5.8
##  [991] 5.8 4.7 4.2 3.7 5.2 5.2 4.7 5.1 4.9 4.7
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.500 5.300 5.300 5.200 5.100 5.000 4.900 4.700 4.700 4.448
## 
## $jack.boot.se
## [1] 0.9355545
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
## [1] 0.5846178
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
##   3.833008   6.729235 
##  (1.645027) (3.085774)
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
## [1] 0.2867095 0.2949235 0.7158902 0.2723924 1.8546401 0.4384653
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
##    [1]  0.683646432  1.035930597  0.753069558  1.989888024  0.949653660
##    [6]  0.493768088  0.415393295  0.580525750  0.640278351  0.554419079
##   [11]  0.172137678  0.702636803  0.527054492 -0.381821428  1.278187254
##   [16]  0.377237879 -0.295046137  0.617727759  0.022716793  0.696713375
##   [21]  0.360849539 -0.238120916  0.541590931  0.056021314  0.489369730
##   [26]  0.692654832  0.523581949  0.146361244  0.812755050  0.507986750
##   [31]  1.123586682  1.143623926  1.633393893  0.199877897  0.282741999
##   [36]  0.757622390  0.813834607  1.283883336  0.662876417 -0.474996283
##   [41]  0.724746384  0.947207279  1.140181447  1.176210542  0.784819356
##   [46]  1.076115170  1.036138759  0.792260464 -1.147510138  0.157926153
##   [51]  0.364381357  0.924271143  0.641524551  0.315212383  0.703934494
##   [56]  0.464040777 -0.254306545  0.902052454  0.636505596  0.204935873
##   [61] -0.050151555  0.396090552  0.636361012  0.252863624  0.795977155
##   [66]  0.021441071  0.274539173  1.220863541 -0.158946368  0.994593189
##   [71]  0.324832585  0.600320530  0.989978597  0.683561619  0.519538794
##   [76]  0.245418756  1.151824633  0.333544494  0.534465791  0.304600731
##   [81]  0.609505970  0.542164081  0.209703255  0.318592392  0.212869685
##   [86]  0.601002500  1.697286576  1.051816490  0.349782198  0.066673545
##   [91]  1.199089680  1.451396574  0.496648367  0.763598651  0.689026875
##   [96]  1.175760599  0.411506757  0.556504612  2.325406594  0.251434853
##  [101]  0.198341547  0.756583057  0.356918099  1.707716131  0.635732650
##  [106] -0.112698947  0.650498823  0.407273790  0.677197012  0.634369876
##  [111]  0.489463193  0.756859868  0.510316999  1.154308844  0.484976555
##  [116]  0.796103925  0.122555730  0.835669543  0.615912066  0.293692552
##  [121]  0.438425823  0.893896419  0.801671867 -1.142142109  0.613349078
##  [126]  0.584617771  1.059112578 -0.478906098  1.818124544  0.357773633
##  [131]  0.892928472  0.461786620  0.446986091  0.464780241  0.672283750
##  [136]  0.298100879 -0.063573915 -0.770849374  0.550472103  0.758183152
##  [141]  0.329437484  0.443565216  0.829943174 -0.111720492 -1.137959420
##  [146] -0.709590063  1.517111983  0.135544211  0.820072645  0.551708071
##  [151]  1.227539903  0.704167462  0.941065861  0.201506787  0.074648052
##  [156]  0.509096763  0.704534272  0.033635592  0.409321963 -0.353593145
##  [161]  1.897390478  0.244461351  1.608158549  2.098371834  0.122076373
##  [166] -1.346577884 -1.104532681  0.627008005  0.894456975 -0.055897272
##  [171]  0.870759416  0.400741644  0.543459597  1.152111652  1.925333671
##  [176]  0.654929816  0.772097062 -0.053069747  0.639832646  0.608256936
##  [181] -0.556087815  1.386035862  0.609975140  0.286582050 -0.079219107
##  [186]  1.298891704  0.340293942  0.368497157  0.094415122  0.642183116
##  [191] -0.689637566  0.973736829  0.360404753  1.188966658 -0.410498177
##  [196]  1.100822795 -0.440341106  0.121538780  0.383125854  0.667387396
##  [201]  0.364582666  0.338742496  0.223361824  0.799815595  0.734838803
##  [206]  0.262854533 -0.529996306 -0.525828605  1.264845888  0.580785407
##  [211]  0.744962906  0.127917489  0.541714410  1.242672156  0.450108833
##  [216]  0.566855699  0.605631739  0.315115819  0.836636321  0.957952977
##  [221] -1.055457824  0.396027409  0.435579969  0.142399520  0.057357669
##  [226] -0.543183564  0.633119875  0.533544171 -0.329119149  0.384414893
##  [231]  0.361962537  0.183862346  0.632116876  1.294472507  0.362533484
##  [236]  0.428197846  0.375740095  1.038805029  0.095531525  1.432841198
##  [241]  0.881141630  0.834419570  0.606246650 -0.326264114  0.738812096
##  [246]  1.333054505  0.467606018  0.712342504  0.565620311  0.702450671
##  [251] -0.031708547  0.346504292  0.196213166 -1.074657330  0.549999300
##  [256]  0.518711328  0.629477897  0.431338204  0.573747128  0.445978506
##  [261]  0.226885608  0.137884018  1.214205440 -0.003389333  0.663979030
##  [266]  1.065221677  0.182052615  0.764095414  0.296963348  0.821381490
##  [271]  1.673431479  0.829943174  1.817625087  0.509237121  0.514674994
##  [276] -1.057945484 -0.210476730 -1.347857449  0.124585443  0.176654859
##  [281]  0.792437117  0.524672821  0.861326594  0.465107261  0.541140076
##  [286]  0.881253378  0.604649016  0.707587345  1.171708231 -0.295272661
##  [291]  0.209183916 -0.500605349  0.519641240  0.437459004  0.504247627
##  [296]  0.067684315  1.220099553  1.087905022  0.848676878  0.566498578
##  [301]  0.385234822  0.453938755  0.500684907  1.980686794  1.029044753
##  [306] -0.471766932  0.588589883  0.775022050  0.221931360  2.229456223
##  [311]  0.591039687  0.874657214  1.193509629  0.827159457  0.259035954
##  [316]  0.340140372  0.295439457  0.139500111 -0.412968271 -0.145408171
##  [321]  0.116439964  0.851400317  1.355590395  0.863960264  0.919035444
##  [326]  1.401943432  0.438439129  0.960965459  0.480799853  0.850875700
##  [331]  0.536207840  0.311595122  0.780937204 -0.097903207 -0.511329518
##  [336]  1.154568918 -0.823834080  0.393447508  0.298427969  0.393761594
##  [341]  0.848228256 -0.776835824  0.362033126  0.759866822  0.335652610
##  [346]  0.054147469  0.137655866 -0.917943053 -0.885893564  0.322639512
##  [351]  0.480933661  1.228549462 -0.230562159  0.904464129  0.453706352
##  [356]  0.340022828  0.960103559  1.279606032  0.328372851  0.122205097
##  [361] -0.051208398  0.558157661  0.545678831  0.399137014  1.227351728
##  [366]  0.595943521  1.068720600 -0.008307194  0.677222750 -0.208734505
##  [371]  0.271943899  0.040748676 -0.990866916 -0.523817340  0.150299118
##  [376]  0.231501451  0.173583766  0.396874758 -0.383464207  0.883615259
##  [381]  0.874657214  0.879381646  0.845037533  0.572075592  0.813834607
##  [386]  0.353382412 -0.354776725  0.015648748  0.988723894 -0.026439899
##  [391]  0.650028884  0.421823463 -0.129776270  0.603672906  0.820904257
##  [396] -1.538111736  0.829835215  0.626196631  2.210794447  0.703305959
##  [401]  0.312907601 -0.863456982  0.968800423 -0.078787277  1.064746724
##  [406]  0.793208763  0.332221265  0.309657398  0.589189506  0.649374321
##  [411]  0.724810417  0.294308435  1.131402485 -0.133945273  0.524019506
##  [416]  0.801377315  0.501357089  0.816881294 -0.251686164  1.030284320
##  [421] -0.028682344  0.794805023  0.178520150 -0.353593145  0.698801510
##  [426]  0.351484783  0.696576480  0.662251150  0.322741858  0.415908551
##  [431]  0.481653235 -0.739920782 -0.017558501  0.172071017  0.368762902
##  [436]  0.053259411  0.513916445 -1.732803872  0.346326498  1.560160767
##  [441]  1.324182412  0.392649517  0.286690378  1.961547005  0.122145852
##  [446]  0.504361678  0.886288243  0.883488570  1.107109477 -0.050490640
##  [451]  0.746506223  0.509991194  0.403687621  0.624976323  0.725125385
##  [456]  0.429171634  0.590050166  0.529065727 -0.347317980  0.969710923
##  [461]  0.738262038  0.344760616  0.840003763  2.416527974  0.506756480
##  [466]  0.237066337  0.696378214  0.606483822  0.526546493  0.731384483
##  [471]  1.003850019  0.362384434  0.469835529  1.060504180  0.502267512
##  [476]  1.576029030  0.353355735  0.080504046  0.269998987 -0.210544766
##  [481]  0.903982466  0.391986908  0.353516072 -0.810424459  0.583212831
##  [486]  1.190900504  0.485617855  1.101926414  0.699066941  0.823068169
##  [491]  0.398680330  2.112651131  0.105334366  0.575250269  0.233846797
##  [496] -0.283011634  0.444568609  0.167057102  0.264824946  0.444531298
##  [501]  0.449525707 -0.370446637  0.169097019  2.133611101 -0.525777676
##  [506]  0.359178471  0.753509676  1.306154343 -1.222185252  0.317939517
##  [511]  0.276172698 -0.063573915  0.387783808  0.467242240  0.586376364
##  [516]  0.313141880 -0.158177861  0.986499787  0.382565334  0.550648843
##  [521]  1.643341533 -1.841879573  0.628705631 -0.008398195  0.489120168
##  [526] -1.162021079  0.077735252  0.611426344 -0.218780061 -0.905848997
##  [531] -0.322473940  0.592819300  0.943923467  0.954196132 -0.614282423
##  [536]  0.339555746 -0.201282038 -1.665892116  0.275815760  0.017181252
##  [541]  0.668940593  0.687656594  1.372127265  0.442037404  0.819045083
##  [546]  0.473723467  0.513743095  1.934687144  0.939033313 -0.432957438
##  [551]  1.141945692 -0.907414222  0.463659474  0.498077038 -0.921842170
##  [556]  0.318962630  0.293873491  0.741186956  0.527129868  1.084987516
##  [561]  0.674913346  1.112085213  0.703573563  0.504503374 -0.802249914
##  [566]  0.888298259  0.540575244  0.222011214  0.019067215  0.876892679
##  [571]  0.758321587  2.168070646  0.724781985  1.381322181  1.585893499
##  [576]  0.362138392 -0.008894297  0.680159123  1.975796943  0.743728930
##  [581]  0.561480961 -0.156690395  0.170691805  0.845124099  0.822979456
##  [586]  1.422978598  0.695833757 -1.226326819  0.457688664  0.431443888
##  [591]  0.292864780  0.475233727  0.743737214  0.265606661  0.180086450
##  [596]  0.696637211  2.270197931  0.678584356  1.183170828  0.850753803
##  [601]  1.165749335  0.643717710  0.871279310  0.937073600  0.383145782
##  [606]  0.702816546 -0.655495165  0.531745717  1.108621775  1.048130935
##  [611]  0.390911675  0.688762105 -0.343031015  0.671884644  0.548453171
##  [616]  0.504920091  0.369079409  0.491789663  0.700942191  0.134222442
##  [621]  0.279532460 -0.157746235  0.215365563  0.492483469  0.957416810
##  [626]  0.203936416 -0.277878168 -0.331453418  0.476053524  0.593302731
##  [631]  1.211342907  0.863384039  0.435346165  1.248369428  0.705897589
##  [636] -0.112227034  0.455329233  0.784058811  0.465107261  0.224533661
##  [641]  0.466418938  0.427543815  0.347641142 -0.700612200  0.609947982
##  [646]  0.830006198  0.683464395  0.348552250  0.767659301  0.661552301
##  [651]  0.580164896  0.524672821 -0.211349140  0.154661048  0.081709452
##  [656]  0.526168993  0.414373261  0.563367863 -0.046618831  0.830006198
##  [661]  0.401424950  0.729840829 -0.083120573  0.590109194  1.124742070
##  [666]  0.614072809  0.469164081  0.374836948 -0.107427679 -0.111088999
##  [671]  1.026088915  0.043203829 -0.090343016  0.417257550  0.068683026
##  [676]  0.581180569 -0.342987769 -0.245859966 -0.824034574  1.505528671
##  [681]  0.661300477 -0.169204315  1.010064071  0.472649371  0.087122880
##  [686]  0.838240484  0.339104598 -0.051422785  0.476465232 -0.165503626
##  [691]  0.152068753  0.125140497  0.622988136  0.894786358  0.500399709
##  [696]  0.677222750  0.844320364  0.280888699  0.079991697  0.361905662
##  [701]  0.164883820  0.601448690  0.758417252 -0.784266243  0.659624270
##  [706]  0.427982133  1.110898991  0.727066246  0.482110052  0.553496343
##  [711]  0.636220821 -0.392344010  1.638192023 -0.225713434  0.655116452
##  [716]  0.685495159  0.066842522  1.026716036  1.655507654  0.636424606
##  [721]  0.428233552 -0.512120956 -0.789248925  0.391425590  0.231501451
##  [726]  0.761212906  0.717898656  1.617271539 -0.257091486  0.341887590
##  [731]  0.642730964 -0.692541080  0.916463211  0.201506787  0.685766024
##  [736]  0.247126205  0.598474982  0.726975729  1.212678325  0.447645641
##  [741]  1.902967313  1.154308844  0.957944264  0.942901374  0.577100997
##  [746]  0.706821559  1.216714803  0.194951243  0.760879377  0.380866368
##  [751] -0.139757410 -0.131877812  0.657457016  0.213675217  0.452457935
##  [756]  0.422453288  1.123886386  0.493135726  0.208040781  0.803982991
##  [761]  0.564288711  1.028375756 -0.042519431  1.820715685  0.339555746
##  [766]  1.167161714  0.342326692 -0.241135081  1.033172509  0.173262742
##  [771]  0.639982309  1.223738399  0.149095962  0.931992389  0.202874510
##  [776]  0.542593750  0.132557575  0.605212315  0.273185089  0.764095414
##  [781]  0.737056287  0.633550411  0.749136890  0.767850739  0.419196036
##  [786]  0.365738604  0.798427258 -0.042611353  0.415914015  0.767173567
##  [791]  0.993417466 -0.013971183  0.879787517  0.593724269 -0.032971067
##  [796]  0.886039006  0.716734087 -0.019448768  1.168128905  0.438197061
##  [801] -0.312401040  0.735886875  0.232401178 -0.446247022  1.039934117
##  [806]  1.146225168  0.600381938  0.521953774  0.551043731  0.452424434
##  [811]  0.610313646  1.314764178  1.332515032 -0.104955173  1.163491896
##  [816]  0.775255677  0.769122393  0.454829649  1.194355301  0.704534272
##  [821]  0.871372867  1.242052072 -0.074620493  1.130845757  0.932656565
##  [826]  1.006445990 -0.153841386  0.521763020  0.440820034  0.211963756
##  [831]  0.355239986  0.425185667  0.574059924  0.418511217  0.807906721
##  [836]  0.010162875  1.323713582  0.915923546  0.320410328  0.897596184
##  [841]  0.845064091  0.721213021  0.667243418  0.130082484  0.487114160
##  [846]  0.686347428  0.562366234  0.499855235  0.325441840  0.357444404
##  [851]  0.247197911 -0.133624485  0.676744223  0.865698525  0.546662988
##  [856]  1.294472507  0.696883101  1.749386232  0.696449991  0.521594518
##  [861]  1.409312757  0.512101624  0.433666725  0.621741496  0.453938755
##  [866]  0.492856843 -1.427093010  0.147196028  0.737883682  0.672825467
##  [871]  0.405458707 -0.443596379  0.261647187 -0.519531995  1.166588069
##  [876]  0.684049956  0.399675939  1.624266960  0.858186173 -0.089321463
##  [881]  0.495070003  1.372009617  0.464106622  0.211027059  0.911752417
##  [886]  1.386611309  2.453253647  0.105334366  0.564921377  0.316940833
##  [891]  0.570632983  0.265905737 -1.104745978  0.418668140  0.220538276
##  [896] -0.666961348  1.410505414  0.421383100 -1.113120334  0.734175579
##  [901]  0.978081302  0.178740347  0.221119768  0.514856059 -0.290321644
##  [906]  0.381955628 -0.577028881  0.516799838  0.283262590  0.817781007
##  [911]  0.027570575  1.179329291  0.652445216  0.902652258 -1.006620513
##  [916]  0.674644735  1.315875769 -0.698415490  0.973324420  0.122287496
##  [921]  0.361128101  0.741774505 -0.220363343  0.351533915 -0.789791364
##  [926]  0.842471410  0.380882155  0.618156624  0.372745460  0.812218675
##  [931]  0.420405475  0.167860879  0.256309821 -0.234300566  0.766464886
##  [936]  0.647786017  0.147160175  0.730282962  0.518913828  0.981919954
##  [941]  0.891751506  0.488185398  0.927181926  0.692407309  0.406248813
##  [946] -0.845878506  0.425530669  0.571162162  0.888950212  0.084479163
##  [951]  1.128550238  0.353057050  1.672609060 -0.807611462 -0.200241437
##  [956]  0.388867389  0.731803834  0.610480224  1.061829540  0.826827345
##  [961]  0.814192364 -1.009616070 -0.667773687  0.252047474  1.982456476
##  [966] -0.344902841  0.185156704  0.582960329  0.339076568  0.671269636
##  [971]  0.066842522  0.939752386  0.432122385  0.529938090 -1.165102740
##  [976]  0.644430700 -0.086482127  0.297484464 -0.167010858  0.257494820
##  [981]  0.489963385  0.034924458 -0.734267945  0.418185571  0.987818454
##  [986] -0.870002458  0.792138576  0.918912637 -0.506576088  1.202027440
##  [991] -0.437684229  0.892505060  0.695447280  0.380517858  0.421157216
##  [996]  0.879093859  0.595490496  1.074529243 -0.959761932 -0.247343063
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
```

```r
fit2
```

```
##       mean          sd    
##   0.56960276   0.28107255 
##  (0.08888295) (0.06284551)
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
## [1] -0.8720188  0.1180424  0.3392822 -0.4398478 -0.7681597 -0.1693603
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
## [1] -0.0415
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8810925
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
## t1*      4.5 -0.05065065   0.8937419
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 2 3 4 6 7 
## 1 3 3 1 1 1
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
## [1] -0.016
```

```r
se.boot
```

```
## [1] 0.8902227
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

