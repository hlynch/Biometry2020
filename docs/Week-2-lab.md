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
## 0 1 2 4 7 8 9 
## 1 2 2 1 2 1 1
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
## [1] -0.0151
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
## [1] 2.739493
```

```r
UL.boot
```

```
## [1] 6.230307
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##   2.5%  97.5% 
## 2.7975 6.2000
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
##    [1] 6.1 3.8 4.7 5.4 4.8 4.5 6.2 4.6 4.8 5.5 5.6 4.0 3.1 5.4 3.9 4.8 4.6 4.8
##   [19] 4.0 3.8 3.8 5.8 3.1 4.4 3.7 2.3 5.1 5.5 4.0 3.6 5.5 6.5 5.0 5.5 4.9 2.7
##   [37] 5.8 4.6 3.1 4.9 4.2 4.7 3.4 4.3 4.8 5.5 4.7 5.8 5.4 5.1 4.7 3.7 4.3 4.7
##   [55] 3.7 5.1 5.0 2.0 3.3 4.0 4.2 3.0 3.8 4.7 5.1 2.7 5.5 5.1 4.1 3.9 5.6 3.8
##   [73] 3.6 4.2 4.7 5.2 5.2 4.7 5.8 5.5 4.8 4.3 5.4 4.7 4.9 4.2 3.6 2.3 3.7 5.0
##   [91] 4.3 3.9 4.0 2.5 3.6 5.0 4.7 3.9 4.3 4.7 4.9 3.1 7.1 3.9 3.4 5.5 5.3 5.0
##  [109] 4.8 3.6 3.0 6.0 5.3 5.0 4.3 5.6 6.6 4.2 4.2 5.0 4.7 4.9 6.3 4.6 5.8 5.0
##  [127] 3.6 4.3 3.8 4.1 5.8 4.7 6.0 3.0 4.1 4.4 4.8 2.9 2.9 4.2 4.8 4.5 4.7 3.5
##  [145] 4.4 2.6 4.8 2.4 5.1 3.9 4.6 5.8 5.1 3.6 4.2 4.4 4.2 5.7 4.1 4.0 3.8 2.5
##  [163] 5.1 4.4 5.1 4.3 3.5 6.3 3.7 3.5 4.0 4.2 3.7 4.5 5.1 3.3 4.0 3.6 4.4 5.2
##  [181] 4.2 3.4 4.0 6.2 5.3 5.0 4.0 4.7 3.5 3.3 3.8 4.9 5.6 3.3 4.8 6.1 4.6 5.0
##  [199] 3.7 4.8 2.9 3.8 4.4 4.0 5.3 3.5 4.8 4.5 4.7 4.1 4.8 4.1 4.8 5.9 4.4 4.3
##  [217] 4.6 4.2 4.5 3.2 4.7 4.0 5.8 4.0 4.5 5.5 5.1 4.6 4.2 4.9 4.0 5.5 4.4 4.4
##  [235] 4.9 4.4 3.7 5.4 3.7 4.8 4.4 4.0 2.7 4.0 6.1 4.2 4.1 5.6 4.7 3.2 4.6 2.0
##  [253] 5.6 4.8 4.5 4.6 5.4 5.1 5.1 4.4 3.3 4.5 5.2 4.1 5.4 4.3 4.6 4.5 6.8 5.1
##  [271] 2.2 6.1 4.7 6.1 2.5 5.3 3.8 6.2 3.4 4.3 6.3 3.3 5.1 5.3 3.8 4.3 4.8 4.3
##  [289] 4.1 3.2 6.2 4.5 5.1 5.5 5.5 3.8 5.1 4.2 5.4 4.8 4.4 3.9 4.4 5.0 5.0 3.9
##  [307] 5.3 4.3 4.5 3.8 5.8 6.0 5.1 4.9 5.7 3.6 5.5 3.2 4.0 5.3 3.5 4.6 3.9 4.6
##  [325] 4.2 3.9 6.0 3.4 4.5 3.8 3.6 3.5 4.7 4.7 3.6 4.7 3.9 4.3 5.7 4.3 5.4 5.6
##  [343] 4.2 5.6 4.4 2.4 5.0 3.9 6.1 3.3 5.6 4.4 5.8 1.7 5.6 5.3 5.4 5.3 6.8 6.1
##  [361] 3.6 4.3 4.3 5.5 4.9 6.1 5.0 3.7 4.9 3.9 5.2 2.8 4.0 4.2 5.3 4.3 5.0 4.0
##  [379] 5.9 7.5 5.7 7.1 4.1 4.8 4.8 4.2 3.7 6.1 4.3 4.6 4.6 4.4 4.6 6.0 4.3 4.1
##  [397] 4.3 6.0 4.2 5.5 5.3 3.1 4.7 3.5 4.8 3.7 2.8 3.3 3.4 5.6 4.1 5.6 5.4 5.4
##  [415] 5.0 5.3 4.9 3.7 3.3 4.4 3.7 3.6 4.0 4.6 5.4 6.0 5.8 4.2 3.9 4.2 5.4 6.3
##  [433] 5.6 4.9 5.8 6.0 5.0 3.7 5.2 5.6 5.3 3.8 3.8 3.6 6.2 3.1 5.0 4.8 4.1 4.0
##  [451] 3.5 4.3 4.6 3.8 5.0 3.8 5.2 3.4 4.0 4.0 4.4 6.3 3.0 4.1 4.1 5.0 6.1 4.7
##  [469] 3.7 5.1 3.3 4.2 1.3 4.5 4.9 4.1 5.2 4.3 5.8 4.3 4.4 4.0 5.4 4.8 5.0 2.1
##  [487] 4.3 5.1 4.8 4.1 3.5 3.6 4.3 4.5 5.1 5.1 4.5 4.7 4.2 4.8 5.5 4.7 4.6 4.6
##  [505] 3.8 4.0 4.0 5.8 4.4 4.4 6.5 6.1 5.7 6.6 4.8 5.3 4.5 3.8 5.2 3.4 3.7 4.1
##  [523] 4.6 4.0 5.6 5.6 5.0 4.2 3.7 4.2 4.9 3.3 3.6 5.4 4.2 2.4 5.3 3.7 5.8 3.5
##  [541] 4.0 4.2 4.6 2.6 4.4 3.2 4.1 6.0 3.8 5.1 2.3 5.1 3.8 4.4 3.5 3.2 5.6 4.3
##  [559] 5.1 4.7 5.1 4.7 4.0 3.7 3.8 5.8 4.5 4.4 4.3 4.6 2.9 2.2 3.3 5.1 4.2 3.8
##  [577] 4.5 4.9 4.4 4.1 4.8 4.0 2.7 5.5 4.1 4.7 5.6 5.1 5.7 5.1 4.4 3.4 4.4 6.0
##  [595] 4.7 6.7 5.6 4.1 5.3 4.3 4.4 4.5 2.5 6.6 4.1 4.3 5.1 4.4 5.1 4.2 4.6 4.4
##  [613] 4.2 6.0 6.1 3.3 6.1 2.7 4.9 4.6 4.7 4.9 3.3 4.4 4.2 6.7 3.1 5.2 6.5 4.5
##  [631] 4.1 3.9 4.5 4.0 3.1 3.7 3.5 4.8 4.5 2.6 4.2 6.6 3.4 3.7 3.1 5.7 2.8 5.3
##  [649] 4.4 4.4 3.6 4.5 4.0 4.6 5.5 5.1 5.3 2.5 4.0 4.9 3.3 5.8 3.9 4.7 3.2 4.1
##  [667] 5.8 4.5 4.2 4.0 4.3 4.2 3.6 5.0 2.8 5.2 3.0 3.5 4.4 4.4 3.0 3.1 5.6 5.7
##  [685] 4.9 3.8 4.3 5.3 4.6 5.3 4.4 4.6 3.6 6.4 5.3 3.7 5.1 4.1 4.3 4.6 4.6 3.4
##  [703] 3.7 4.1 4.5 3.4 4.9 5.8 4.5 4.3 5.2 2.8 4.4 6.0 4.9 5.1 3.9 5.1 4.9 4.2
##  [721] 4.9 3.3 2.8 5.1 4.7 4.2 3.4 4.9 3.6 4.3 4.0 5.1 4.1 3.2 4.4 4.4 4.2 4.9
##  [739] 5.4 4.9 3.3 4.9 3.6 4.6 4.2 3.3 4.6 3.9 3.9 3.7 3.2 5.1 5.3 5.4 4.7 5.6
##  [757] 4.3 4.5 3.6 4.5 4.1 3.5 4.0 5.5 4.9 5.8 4.6 4.4 5.0 5.7 4.0 4.1 4.9 4.2
##  [775] 3.8 5.5 5.9 4.6 4.3 5.0 4.4 3.9 4.9 5.8 4.5 4.2 4.9 3.8 4.2 2.9 3.7 5.7
##  [793] 4.3 2.9 3.6 4.3 4.0 4.1 5.1 4.0 3.5 4.8 4.6 5.7 4.8 3.3 4.5 5.5 4.5 4.8
##  [811] 4.8 3.8 4.1 4.2 5.3 4.4 4.1 6.0 3.7 3.4 4.6 4.1 3.9 4.8 4.7 5.5 4.7 4.2
##  [829] 4.6 4.1 4.6 5.5 5.0 3.5 6.1 4.2 2.2 6.0 4.0 4.6 4.3 3.3 5.3 3.8 3.9 5.3
##  [847] 3.4 3.6 4.3 5.8 5.6 3.4 4.9 5.8 5.0 5.4 4.6 4.4 3.5 2.6 4.0 5.3 5.8 5.4
##  [865] 5.3 5.0 4.0 4.8 4.1 4.3 4.4 3.9 4.4 5.7 3.9 4.1 4.3 4.0 3.8 3.1 4.5 5.2
##  [883] 5.2 5.1 2.8 3.9 3.6 4.6 3.6 4.3 5.7 4.3 3.3 3.9 2.9 5.7 4.6 3.9 3.3 2.5
##  [901] 5.6 6.1 5.0 4.8 5.0 4.3 5.0 4.6 5.8 2.9 5.8 4.1 5.4 3.9 5.3 5.3 3.6 4.6
##  [919] 3.9 5.4 4.1 4.9 4.5 3.7 5.1 3.3 5.1 3.8 5.2 2.9 4.0 5.3 5.2 5.5 4.1 4.3
##  [937] 5.5 4.4 5.1 4.7 4.1 4.2 5.3 4.0 5.3 5.2 4.5 4.6 4.3 3.9 3.8 4.5 3.0 2.8
##  [955] 4.5 3.2 2.7 3.4 6.1 4.2 4.8 4.5 4.0 4.0 5.4 4.3 4.3 5.1 2.0 4.1 3.6 3.4
##  [973] 3.9 5.0 4.0 3.6 4.4 4.9 4.4 4.7 4.5 4.2 3.5 4.9 3.2 4.0 4.2 4.1 4.5 4.9
##  [991] 4.1 4.2 3.4 5.1 4.0 3.1 4.6 6.4 2.6 3.4
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
##   2.6   6.2
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
##    [1] 4.0 2.4 4.1 3.3 5.7 4.3 4.2 5.0 3.7 3.8 5.0 5.2 5.6 4.5 5.2 4.4 5.6 3.6
##   [19] 5.1 4.8 4.5 5.6 4.0 6.4 5.0 4.7 2.9 4.3 5.6 5.6 4.3 4.7 5.4 3.6 3.8 5.5
##   [37] 4.3 3.8 4.3 4.4 3.8 4.4 5.0 4.4 2.6 4.8 5.1 4.1 4.2 1.7 3.4 2.9 4.9 5.6
##   [55] 3.2 3.6 4.9 3.6 4.2 6.0 3.9 4.3 2.9 3.7 5.0 5.2 5.0 4.6 2.5 4.0 5.0 3.5
##   [73] 5.3 4.0 3.9 4.9 3.7 4.1 3.9 5.5 2.8 4.5 3.7 4.3 3.3 2.7 4.8 5.0 4.7 4.8
##   [91] 3.6 5.7 3.5 3.8 5.6 3.1 7.3 5.1 6.6 4.1 5.2 4.4 4.9 4.9 3.3 4.4 4.0 4.4
##  [109] 3.7 5.6 3.5 4.5 3.8 4.4 5.3 3.7 4.2 4.8 5.3 3.9 4.7 5.2 5.3 3.5 4.2 5.5
##  [127] 4.2 3.8 5.6 5.1 3.4 3.6 2.7 5.5 5.3 6.4 5.1 3.9 4.4 4.2 4.7 4.7 5.1 4.3
##  [145] 5.3 4.3 4.5 4.7 6.1 4.1 2.6 4.8 4.8 4.2 5.2 5.4 4.1 5.9 3.9 5.4 4.4 3.2
##  [163] 2.7 5.0 5.0 4.9 5.0 4.3 1.7 4.8 5.0 5.1 4.2 4.2 3.9 4.2 4.0 2.3 4.7 5.3
##  [181] 6.0 5.1 4.1 2.8 4.3 3.7 4.2 3.5 6.8 4.8 6.3 3.8 5.3 4.3 5.4 4.7 4.6 5.1
##  [199] 5.9 5.0 4.1 3.3 2.8 4.6 4.8 3.2 5.7 5.0 3.9 3.6 4.9 2.6 4.8 3.7 4.8 4.8
##  [217] 5.1 5.4 5.1 4.1 4.6 4.1 5.7 4.7 5.7 4.5 4.8 5.0 4.7 5.7 5.2 3.6 3.5 3.4
##  [235] 3.9 5.7 5.9 4.5 4.4 4.2 4.6 5.1 4.8 4.3 3.8 3.9 5.2 3.4 4.0 3.4 2.6 5.9
##  [253] 4.8 4.1 6.5 4.1 4.4 3.9 4.4 6.6 5.0 5.5 4.6 3.8 5.2 3.6 5.1 6.0 4.5 4.3
##  [271] 2.7 3.5 3.2 4.2 5.2 6.0 3.6 4.7 5.4 4.0 4.7 6.2 5.4 4.5 5.9 3.9 3.6 4.6
##  [289] 3.7 4.0 4.0 3.8 3.3 3.3 5.5 5.8 5.8 3.6 6.6 5.4 4.0 4.3 5.2 5.5 3.8 3.2
##  [307] 5.5 4.0 4.7 2.8 5.6 4.6 4.0 4.9 5.1 3.9 5.5 4.8 4.4 3.7 4.3 4.0 4.8 6.0
##  [325] 3.3 4.2 3.4 4.6 4.2 3.7 4.5 4.3 4.8 4.4 4.4 4.4 5.4 3.7 3.7 3.9 2.5 4.4
##  [343] 4.3 4.8 4.3 5.1 4.4 3.9 5.3 3.4 4.8 4.0 4.9 3.1 3.2 4.6 5.8 4.8 5.3 3.7
##  [361] 5.8 5.5 4.0 4.2 5.7 4.1 4.5 5.9 6.0 5.0 3.7 3.9 5.1 3.5 3.9 4.5 3.4 4.8
##  [379] 5.5 2.7 6.1 5.8 4.6 4.9 5.1 3.6 4.8 5.3 5.3 4.8 4.7 5.3 5.2 4.9 3.5 5.1
##  [397] 4.7 3.8 4.2 2.8 5.2 4.5 5.9 4.8 5.4 4.3 3.7 6.1 5.1 4.5 4.0 4.5 4.0 4.4
##  [415] 5.0 3.6 4.3 3.2 5.2 6.3 4.5 2.8 3.5 6.2 4.5 3.8 4.6 3.9 2.9 6.0 4.5 4.1
##  [433] 5.6 3.5 4.5 3.9 5.0 5.5 4.9 5.4 4.8 6.1 4.6 5.0 4.1 4.6 6.2 3.9 3.3 4.1
##  [451] 2.3 3.9 3.2 6.2 5.1 5.5 4.3 4.1 3.8 5.2 6.1 3.1 4.0 2.7 5.2 5.3 3.7 5.5
##  [469] 3.8 3.3 5.0 5.2 5.7 4.7 4.7 5.1 4.3 5.3 5.5 5.3 3.8 4.5 4.5 5.1 4.3 5.2
##  [487] 5.3 6.0 5.4 6.5 4.7 4.5 5.9 3.4 3.3 4.5 3.5 2.9 4.6 3.3 4.8 5.2 4.9 4.9
##  [505] 4.0 4.2 5.0 5.4 5.9 5.4 3.4 5.8 3.4 3.1 4.8 5.1 3.7 4.4 4.3 5.4 5.7 4.2
##  [523] 5.1 3.3 4.0 4.6 3.9 5.1 5.7 6.1 4.4 4.8 6.7 5.0 4.1 5.1 5.3 4.2 3.6 3.5
##  [541] 4.2 3.3 5.1 4.3 4.6 4.7 5.3 4.0 5.0 4.4 4.3 4.0 3.9 5.3 3.8 4.2 4.4 5.3
##  [559] 4.8 4.9 4.1 4.7 4.9 5.3 4.9 4.6 4.5 4.3 5.8 3.6 4.9 5.0 5.6 4.4 5.1 4.5
##  [577] 3.0 4.0 7.2 4.0 5.1 4.8 3.7 4.4 4.6 5.1 7.3 3.6 4.5 5.0 5.3 3.0 3.8 3.6
##  [595] 5.3 4.3 4.3 4.6 4.6 4.2 4.7 5.5 4.9 6.9 4.8 3.9 4.1 2.9 4.6 4.7 4.1 4.4
##  [613] 4.5 4.5 5.9 4.6 5.5 4.9 4.5 4.5 4.9 4.5 5.2 5.5 4.3 5.2 3.3 4.9 4.2 3.3
##  [631] 5.8 5.0 3.3 4.2 3.8 4.6 4.7 4.5 3.2 5.0 4.4 4.8 2.2 3.8 4.2 3.6 4.3 4.2
##  [649] 6.3 4.9 3.5 6.3 3.7 3.7 4.7 5.3 5.4 5.4 3.5 5.4 4.6 4.9 4.7 5.7 5.7 4.9
##  [667] 4.7 7.1 5.2 7.0 3.3 4.5 3.9 3.6 4.3 4.7 2.7 5.3 4.9 3.8 5.6 2.6 3.5 5.5
##  [685] 5.0 4.2 4.4 4.3 7.2 4.1 4.4 5.5 4.7 4.9 6.1 4.1 5.4 3.6 4.8 3.8 2.6 3.5
##  [703] 5.0 5.3 5.1 2.5 5.2 4.3 4.1 4.6 5.3 4.2 5.5 5.8 3.8 4.9 4.8 5.1 4.7 4.7
##  [721] 5.1 3.6 6.1 5.1 5.1 3.2 5.2 3.2 4.7 6.3 5.2 5.4 5.4 3.5 3.8 3.5 2.5 4.3
##  [739] 4.1 4.6 5.7 5.9 3.8 3.0 4.1 5.9 4.4 3.6 2.5 4.1 4.4 6.2 4.9 4.6 4.5 6.2
##  [757] 3.7 3.1 4.8 5.5 4.9 4.0 4.5 4.9 4.4 4.6 6.4 4.3 4.0 3.9 6.2 3.9 3.4 4.5
##  [775] 3.6 3.8 5.0 3.1 4.3 5.0 3.9 5.1 4.3 5.4 5.2 3.5 5.3 4.4 2.0 4.3 6.2 3.7
##  [793] 4.7 4.6 3.0 3.3 6.2 5.7 5.0 5.0 4.9 4.1 6.3 5.1 5.0 5.5 3.9 2.9 4.0 4.0
##  [811] 5.8 3.9 4.7 4.0 5.9 5.8 5.1 4.7 3.9 3.1 3.8 4.2 5.2 5.2 3.8 4.9 6.9 5.7
##  [829] 3.7 4.3 5.2 3.3 3.1 4.7 2.9 6.4 4.9 6.0 4.4 3.9 3.9 4.6 5.4 6.3 4.1 3.2
##  [847] 2.8 5.3 2.9 3.3 4.7 4.4 3.0 5.7 3.5 5.7 3.0 3.2 4.2 4.2 3.6 4.6 4.7 5.1
##  [865] 3.6 3.9 3.5 6.2 4.5 3.9 4.2 5.0 4.0 4.5 3.7 4.5 3.1 3.6 3.6 3.9 5.4 3.9
##  [883] 6.5 5.0 4.9 3.9 3.7 4.6 3.2 5.6 7.1 4.9 4.9 2.1 5.8 4.7 3.4 5.8 5.7 5.2
##  [901] 3.5 4.2 5.9 2.3 3.4 3.9 4.1 3.8 4.3 5.3 3.9 2.5 3.9 5.4 2.7 4.6 4.2 5.5
##  [919] 3.7 6.2 3.9 4.5 3.3 4.1 4.0 5.6 5.0 3.4 4.5 4.7 4.4 4.7 3.7 5.5 5.2 3.4
##  [937] 3.2 3.0 5.2 5.3 4.2 4.8 3.4 4.2 3.6 3.3 5.5 3.9 4.4 4.0 3.7 6.1 4.6 5.5
##  [955] 3.7 3.6 5.5 5.2 4.4 5.4 4.6 4.5 4.4 3.3 3.7 3.6 5.2 4.4 4.5 5.1 5.0 5.0
##  [973] 4.0 3.6 5.8 5.1 4.8 3.7 5.4 5.1 3.7 3.6 4.4 5.0 4.4 6.3 4.3 6.0 2.9 4.7
##  [991] 6.1 3.1 3.2 5.3 3.7 3.9 5.9 3.7 4.7 3.1
## 
## $func.thetastar
## [1] 0.0154
## 
## $jack.boot.val
##  [1]  0.537869822  0.482608696  0.320281690  0.170454545  0.067055394
##  [6]  0.007278481 -0.133701657 -0.230153846 -0.468154762 -0.506764706
## 
## $jack.boot.se
## [1] 1.038856
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
##    [1] 4.0 5.5 6.1 4.6 3.6 4.8 5.0 4.6 4.8 5.5 4.9 5.8 4.9 3.3 3.9 3.6 4.4 3.0
##   [19] 4.8 2.8 4.0 5.2 4.4 4.3 4.3 4.1 4.2 4.3 4.0 1.6 5.8 4.0 4.2 5.1 5.7 5.5
##   [37] 4.7 6.3 5.2 3.6 3.4 4.8 4.7 4.0 7.0 4.4 4.4 5.2 6.3 4.3 4.8 3.8 3.8 6.4
##   [55] 4.2 5.5 2.7 5.4 2.7 2.2 3.9 5.9 5.4 6.6 4.1 5.6 2.4 3.9 6.3 4.7 4.3 4.8
##   [73] 3.6 5.5 4.5 3.9 4.6 4.0 4.5 5.2 4.7 2.5 4.8 6.2 4.0 4.5 5.4 4.0 4.1 5.1
##   [91] 5.3 4.0 4.2 4.0 4.6 4.0 4.0 5.3 5.1 5.1 4.2 3.7 4.3 4.4 4.7 5.8 3.3 3.7
##  [109] 4.1 4.8 4.3 3.9 4.8 5.2 3.1 3.5 5.4 3.7 4.2 2.6 3.5 4.0 4.0 4.4 3.1 3.0
##  [127] 3.8 4.1 3.7 5.3 3.9 3.8 3.8 3.2 4.7 5.5 5.0 4.6 4.6 5.7 4.5 5.4 5.9 2.6
##  [145] 5.2 4.4 4.8 3.3 5.3 4.8 5.6 4.2 2.9 5.1 3.3 3.4 5.7 5.2 4.0 4.4 3.2 5.7
##  [163] 4.0 3.8 5.4 3.6 5.8 3.7 4.2 4.9 3.2 4.6 4.7 5.6 4.7 3.7 5.1 5.7 5.6 5.7
##  [181] 3.8 5.3 4.2 4.3 5.3 5.1 4.9 4.7 4.2 3.0 3.2 4.7 5.3 4.7 6.3 5.1 6.6 4.3
##  [199] 4.2 5.1 5.3 4.6 5.5 4.4 3.8 5.4 4.8 4.6 4.5 4.1 5.2 6.3 5.4 4.7 4.4 4.9
##  [217] 4.5 5.5 4.6 4.7 4.4 5.1 5.3 4.2 6.0 3.7 4.8 4.4 3.7 5.7 6.0 6.6 4.1 4.7
##  [235] 4.5 5.3 4.4 6.1 3.1 5.3 4.1 4.5 3.7 5.8 3.6 4.3 4.8 6.0 5.2 2.7 4.8 4.0
##  [253] 2.4 4.6 3.6 3.4 4.6 3.6 4.7 3.6 4.7 4.5 3.5 4.9 4.3 4.4 4.5 3.0 5.0 3.6
##  [271] 5.1 5.3 4.0 5.1 5.0 5.8 5.6 4.1 4.9 4.7 4.1 3.1 4.7 3.4 3.8 4.8 4.1 4.8
##  [289] 7.3 4.4 4.8 4.2 6.1 4.3 2.9 4.6 4.1 3.6 4.2 4.5 6.3 4.1 4.4 4.1 2.8 2.0
##  [307] 3.9 2.3 4.8 4.2 4.5 3.7 5.1 4.2 3.2 3.0 4.0 4.3 4.6 4.3 5.3 3.4 5.1 5.4
##  [325] 3.1 4.6 3.7 4.8 4.1 5.0 4.9 5.8 4.6 5.2 3.4 4.9 4.0 5.1 5.1 4.8 5.0 3.5
##  [343] 4.2 4.1 4.1 4.8 5.1 4.7 4.5 5.4 4.3 3.7 4.9 5.1 4.1 3.6 3.4 4.5 4.1 3.9
##  [361] 3.8 6.6 5.7 4.2 4.1 6.0 5.3 4.5 5.1 4.1 5.6 4.8 2.6 4.5 4.7 3.2 5.0 3.9
##  [379] 2.8 5.1 3.5 3.7 4.4 2.7 3.4 5.2 3.5 5.1 5.7 5.0 4.4 5.1 3.3 3.7 5.6 3.3
##  [397] 5.1 3.6 5.2 4.2 3.7 5.6 2.9 5.3 6.2 6.0 3.7 5.6 6.2 3.7 4.5 5.1 2.8 2.5
##  [415] 3.7 3.2 4.5 3.9 5.7 4.4 3.8 4.2 3.9 3.8 3.6 3.7 4.3 4.5 5.7 4.5 3.3 3.2
##  [433] 4.4 4.0 2.7 5.3 5.3 3.9 3.5 5.3 4.2 3.9 4.1 4.1 4.0 4.9 5.3 5.7 3.5 5.6
##  [451] 4.2 4.4 3.9 5.3 3.7 4.6 4.4 5.2 4.0 4.3 3.9 5.4 4.8 4.0 4.9 3.6 4.9 4.8
##  [469] 3.9 4.1 3.7 5.1 3.7 5.0 5.8 4.4 5.3 5.8 4.1 6.6 5.0 5.6 4.2 4.2 3.8 5.6
##  [487] 4.8 4.9 6.1 3.3 5.5 2.8 3.8 4.0 5.2 3.8 5.0 4.6 3.4 4.7 5.5 4.7 5.0 3.7
##  [505] 5.8 4.5 3.8 5.6 4.9 3.9 6.6 4.1 3.9 4.1 5.1 4.5 3.8 3.4 4.1 5.0 3.7 5.3
##  [523] 5.8 3.7 4.5 3.1 4.5 4.0 4.6 4.8 4.0 3.1 5.6 5.5 1.8 4.3 2.1 5.2 3.7 3.1
##  [541] 5.6 3.3 4.5 3.3 5.5 4.4 3.8 5.7 5.2 3.7 4.5 5.1 4.3 4.9 5.0 3.6 5.0 3.0
##  [559] 3.3 4.8 4.7 4.8 3.2 3.3 4.3 3.7 3.7 4.4 3.6 4.6 4.4 4.9 5.5 5.3 4.0 5.3
##  [577] 5.6 5.0 4.6 3.9 5.0 4.1 5.1 5.0 3.2 4.2 4.6 5.6 4.5 4.7 4.3 3.6 3.6 4.6
##  [595] 4.9 5.2 5.3 5.1 5.2 4.4 5.2 3.9 5.4 3.5 6.0 3.3 5.4 3.7 4.6 3.3 3.8 4.9
##  [613] 4.2 4.3 5.1 5.4 4.8 4.2 4.7 3.7 4.0 6.0 4.5 5.6 4.5 5.5 5.8 3.7 3.5 4.9
##  [631] 3.9 6.7 4.8 2.7 3.9 5.2 5.8 3.9 4.2 3.0 4.8 4.5 5.5 2.7 5.3 4.8 4.3 5.0
##  [649] 4.2 5.1 6.0 5.2 2.7 5.3 4.1 5.2 4.1 3.7 4.6 5.3 4.9 4.4 4.7 5.2 3.2 6.0
##  [667] 4.3 4.1 5.8 4.2 4.0 4.5 3.9 5.0 4.5 4.5 4.9 3.9 3.7 5.4 4.3 4.9 5.2 5.3
##  [685] 3.9 4.5 5.0 5.1 5.6 4.4 5.2 5.3 4.0 2.9 4.1 4.9 4.5 6.0 4.1 3.5 3.2 4.4
##  [703] 4.9 3.1 4.6 5.4 4.0 5.5 5.0 4.3 4.1 4.4 2.4 4.6 4.9 3.7 4.7 3.9 4.5 3.5
##  [721] 5.8 4.3 4.8 4.4 5.1 6.3 3.8 3.7 3.6 5.2 2.8 4.2 5.6 4.1 4.6 5.4 4.7 5.9
##  [739] 4.7 4.2 2.9 5.4 4.4 4.1 4.7 6.0 4.3 4.5 5.3 4.3 4.9 4.0 4.8 5.2 4.1 4.2
##  [757] 4.2 3.7 5.9 5.1 4.6 2.6 4.1 3.6 4.1 6.5 3.4 4.0 2.9 3.6 4.2 4.8 3.3 3.1
##  [775] 4.2 4.0 3.2 5.5 3.2 5.1 4.3 4.4 4.7 4.1 4.7 3.2 4.9 4.0 5.2 3.9 5.1 6.2
##  [793] 6.0 4.9 4.9 2.9 4.8 4.9 4.3 5.0 5.6 4.5 3.8 3.4 3.4 4.6 5.9 6.2 5.6 4.8
##  [811] 4.1 4.1 3.4 2.9 4.9 6.3 4.4 4.1 3.4 4.5 3.4 6.7 4.4 3.8 3.7 5.0 3.1 6.8
##  [829] 4.3 3.1 4.8 5.2 5.3 4.2 3.2 4.1 3.9 4.5 3.5 4.6 4.8 5.2 3.4 4.3 3.4 5.8
##  [847] 5.3 4.8 4.3 5.7 3.6 5.2 5.2 4.1 4.2 4.1 1.3 5.7 4.2 4.1 4.3 4.3 5.2 3.9
##  [865] 4.0 3.7 7.0 5.9 4.2 4.8 5.5 3.8 5.7 5.8 5.1 5.1 6.0 3.8 5.1 4.9 5.4 4.1
##  [883] 4.6 3.5 5.0 4.1 5.0 4.5 3.6 4.6 5.2 5.6 4.4 3.8 2.9 4.2 6.3 5.0 5.0 4.0
##  [901] 3.8 4.1 4.8 5.7 3.2 4.5 3.3 4.8 3.5 3.6 5.0 3.1 5.3 4.3 4.4 4.1 5.9 5.3
##  [919] 3.9 5.5 2.7 3.4 3.1 3.1 4.2 4.8 3.3 3.8 4.8 4.8 4.0 4.4 5.0 6.0 2.8 3.9
##  [937] 3.7 5.6 2.9 3.8 2.3 3.8 5.5 3.5 4.4 6.0 5.9 3.7 2.8 4.3 5.5 4.6 4.5 4.8
##  [955] 5.1 5.1 4.3 3.3 4.4 4.9 3.1 3.4 4.1 4.4 4.9 4.0 4.4 5.0 5.6 4.6 4.6 4.8
##  [973] 4.3 5.3 2.2 2.9 2.6 4.4 3.5 3.5 5.2 4.2 4.2 4.5 3.6 3.1 5.2 4.2 3.3 3.4
##  [991] 4.5 6.2 4.8 5.5 5.0 4.5 4.6 5.7 4.9 3.3
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.30 5.40 5.30 5.10 5.00 4.90 4.80 4.78 4.50 4.40
## 
## $jack.boot.se
## [1] 0.9634438
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
## [1] 0.9946638
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
##   3.972902   5.684245 
##  (1.707387) (2.603996)
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
## [1]  1.6112344  1.0459881  0.4788180  1.1673952  0.9495311 -0.4499067
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
##    [1]  1.0090838833  0.0162666201 -0.6424930478 -0.0101971597  0.7026611632
##    [6] -2.1278283923  1.0885860754  0.9666438043  0.8945033882  0.5571411527
##   [11]  0.9029575447  0.9864588548 -0.7057620958 -0.3293677863 -0.7821464456
##   [16]  0.6982748695  0.4207432090  0.4284274958  1.2969402186  0.6897159893
##   [21]  0.0643267760  1.9173818949  0.0180400351 -0.7715869718  1.1897831688
##   [26] -1.3289241924  0.9575826621  0.0220093253  1.3922468974 -1.3402926400
##   [31]  0.0223705058  0.8557685874  0.6196466727  1.5512471981  0.7813280149
##   [36]  0.7661221577  0.0004864566  0.9499553244  0.0329045576  1.3596638925
##   [41]  1.3115386515 -0.3307550597  0.1056208160  0.9688942704  1.3919003126
##   [46]  1.1431417787  1.2297282453  0.6841961436  0.6984725028  0.8570768693
##   [51]  1.5083589008  0.8788335689  1.6290422462  0.9671209547  0.6500151318
##   [56] -0.3424068903  1.5865696231  0.8160337970  0.4009302433  0.9424994007
##   [61]  0.7949327329  1.1824105853  0.9026475446 -0.2828858005  1.1738584575
##   [66] -0.7691357278 -0.2389607749  0.0317550099  1.2620591240  1.0878325806
##   [71]  1.3612879573  1.3523971600  0.7646379541  0.0732408445  0.3757102219
##   [76] -1.3037869629  0.4003095658  1.0198422519 -0.7601660968 -1.1223506510
##   [81]  0.8281919918  0.5080975234  0.3900705122  1.2676879696  0.6560930613
##   [86]  1.0817407493 -0.0005736267  0.0742984250  0.9150788299  0.8450128811
##   [91] -0.7814382750  1.0633073597  1.3240292544  1.6765651075 -0.7309208777
##   [96]  0.5002845749  0.4450240008 -0.8499947397  0.9354997460  1.1288323011
##  [101] -0.3200741876  0.9850473578  0.7766681981  0.5120284453  0.4688484065
##  [106]  0.7389413612  1.1875969527  0.8041150224  0.8202714834  1.8603028000
##  [111]  0.8678953849  0.8023519096 -0.2659484698  0.8326267400  0.9432935164
##  [116]  0.0929109028  0.9039677178 -0.8202706994  0.8706103752  1.0708844583
##  [121]  0.4049460443  1.0286869866  0.0546242999  0.8311455353  0.6002967949
##  [126]  0.6803786343  0.7912027849  0.7115981566 -1.4263801184 -0.8489955605
##  [131]  1.2316021492 -0.4007280181  0.1148588220  0.0223705058  0.9304451092
##  [136]  0.8076897015  1.0279742917 -0.0017056698 -0.4015916321 -0.3879360444
##  [141]  1.0701636997  0.7433854389 -1.1522554064  0.8490920637  1.1354021892
##  [146] -0.2875298335  1.3990468681  1.1076508023  0.7734660492  0.8085035209
##  [151]  0.9697356916  1.0758574762  0.6297720851  1.0232093589  0.9957188210
##  [156]  0.9210051237  1.0672369192  0.8366402159  0.6865860143  1.0844629696
##  [161]  1.2237756546  0.8361987774  0.3870004084  0.5259745233  1.4672784568
##  [166]  0.9652412895  0.6009533454  0.7201742736  0.8956494955  0.6316984120
##  [171] -0.2387172344  1.7986656390  0.4344857967  0.7095279928  0.9443458270
##  [176]  1.4586824805 -1.4711439516  0.8103433249  0.2956484853 -0.3448433613
##  [181] -0.0003464674  0.4003254778 -0.0006224711 -0.3969826943  1.7763770729
##  [186]  1.0184569720  0.0917051616  1.0053588185  0.7845487113  0.1159669929
##  [191]  1.4088723186  0.9988106989  1.1996584367  2.3485508140 -1.4009515608
##  [196]  1.1296937429  0.8579868165 -0.2075533780  1.2822956309 -0.4015002322
##  [201]  0.5654554811  1.3589885167  0.8919196139  0.9397476379 -0.8286217439
##  [206]  0.1808355938  0.7646379541  1.2173938666  1.3537966359 -0.2894726474
##  [211] -0.8521955496 -0.6993664208  0.8610069887  1.1511215379  1.6197606389
##  [216]  1.0456575042  0.1654302778  1.0810599697 -0.7303695009  0.4131083380
##  [221]  0.7003845269 -1.2810428630  0.3820164413  0.7210185621  1.5967199625
##  [226]  0.1837424822 -1.4844622667  1.4251091710  0.3953782377  1.4442066279
##  [231]  0.6813856998  0.5572392250  1.1234943831  1.4119989581  1.2795592270
##  [236]  0.4031670789  1.0337332376  0.4856240140  0.4039095566  0.9877667423
##  [241]  0.5190347556  0.5689398593  0.4075117610  1.1793193555  1.5637568522
##  [246]  0.4761841641  0.9131754830 -0.2737935932 -2.1505730176 -0.7499451363
##  [251]  0.4611043597  0.0909474385  1.3267635249  0.0353785858 -2.1002247636
##  [256]  0.7245383036  1.2808919097  0.7485728239  0.8879200341  0.8370503355
##  [261]  0.6299768712  0.6804130087 -0.4076551120  0.8317220914 -0.7314380806
##  [266]  0.5008287549  0.4523988140  0.1360241619  1.1609747453  1.4560659035
##  [271]  0.7277050769  1.6531407933 -0.8162015613 -0.8533632935  0.4463398028
##  [276]  0.9343161842 -0.6793516633  0.9227478062  1.1917548909  1.8702556454
##  [281]  1.3406988639  0.7313867343  1.0856390346  0.8080999316  0.5272546759
##  [286]  1.2633712313  0.2173289577  0.4897885607 -0.3302016774  0.4406331325
##  [291]  0.3863574774  0.3168950769  1.0559624230 -0.0316281294  0.8194132358
##  [296]  1.0298078768  0.7039878052  0.6448137210 -0.7750153073  0.6994133675
##  [301]  0.9488146695 -0.7119685769  0.9526185243  0.9744402115  0.8680962108
##  [306]  1.0105033327 -0.2947946324  0.9999568389 -0.6805935185  0.8104133217
##  [311]  0.8489236157  1.2780384158  1.3333141735  1.8126062980  1.0578238778
##  [316]  0.4693536771  0.9059984288 -0.3713016373  0.8566805864 -0.6503493123
##  [321]  1.3691340824  1.5606057705  0.0103343972  0.0626694889  0.7228187365
##  [326] -0.0199130977  0.8951966173  1.4723253853  0.8401709363  0.0744233545
##  [331]  0.9439022621  0.6780625762  0.0929204383 -0.2278480115  0.0223675427
##  [336]  0.7525920917  0.8751674428 -1.1375697559 -0.2971810050  0.5307789843
##  [341] -0.3336304152  0.4531735642  1.2781213445  0.4946902549 -0.3140582695
##  [346]  0.9382587571  0.5020731358  0.5194677297  0.0484924218 -1.1347655859
##  [351]  1.2527373959  0.6580284711  0.4841871487  1.2486496043  1.2267086606
##  [356] -0.8416512351 -1.2922916456  0.6199461643  0.8268321922  0.9168377613
##  [361]  0.9859599238 -0.3998365924  0.4683192419  0.6549923942  0.9645215764
##  [366]  0.7074713168  0.5799121453  1.1522812703  0.3713273679  0.7196086643
##  [371]  0.8555796561 -0.3947520579  1.3459887334 -0.7112014389 -0.8516514502
##  [376]  0.7351484722  1.2188583136  0.0636188309 -1.2753661729 -0.7787691509
##  [381]  0.7949327329  0.1371149808  0.7423341131 -1.3893939320  1.3937385281
##  [386]  1.0986412387  1.0132651427  1.3584099919  1.2777055829  0.8566283731
##  [391]  1.2327410508  0.7474090252  0.6874625654 -0.3397174637  0.4620487306
##  [396] -1.3068979361  1.7367390226  0.5693044953  0.7866828797 -0.3684727815
##  [401]  0.9283196145  1.1742397277 -0.3852975722  0.4578138698  0.0102273611
##  [406]  0.7399828240  0.6612815324 -0.4146927225  1.0504978020  0.9289898933
##  [411]  0.8260450231  0.9019277083  0.4698883050  0.7479771039  0.2742558855
##  [416]  0.8812005310 -1.1554912140  0.6196680428  0.4551380756  0.9160634767
##  [421]  0.5946522633  1.3408042465  0.7595557909  0.9142733754 -0.3701483827
##  [426]  0.2507997048  0.8738794957  0.8399834550  0.8011903044  0.8477540751
##  [431]  1.1850571906 -1.3303376705  0.0355037916 -0.7363289479  0.9502335008
##  [436]  0.4725927553  1.1117431697  0.7270114798  0.8845665743  0.0022668879
##  [441] -0.2934862010 -0.2945504696  1.0243858687  1.4757724685  0.9835942389
##  [446]  0.3888008134  0.5172555554 -1.1646748272 -0.0077677456  1.1338254517
##  [451]  0.7271812872  0.8801489871 -0.2744584824  1.0013657128  1.1103918851
##  [456]  0.9501799798  0.9961920653  1.3740172162 -0.3262970658  0.3294258365
##  [461]  0.5838271087 -1.3217781496  1.1462375940  0.8921322622  0.6940333285
##  [466]  0.9444639739  1.2040522519  0.7569232271  0.9426181943  0.0400273274
##  [471]  1.2108286541  0.9235177666  1.0861416014  0.8968371681  1.0044664368
##  [476]  0.9391811952  0.8023300932 -0.8162015613 -0.2113416234 -0.2738025432
##  [481] -0.6568841361  0.6940574871  1.4987019293  1.0243948556 -0.2101653818
##  [486]  0.8622952440  1.0040704012 -0.6507364246  0.8446277233  0.6337773284
##  [491] -0.6416142334  1.1994616368  0.7009180192 -1.3202599283 -0.2647297268
##  [496]  1.2163765118 -0.6383913872  0.7403681344  1.5459363560  0.1009146364
##  [501] -1.1324861030  0.9202268747  0.0689326113  0.4918062158 -0.7228226689
##  [506] -0.0441403010 -0.3330699247  1.0831297860  0.5885331123  0.6845286934
##  [511] -0.7974437901  1.1644689844  0.8732704232  0.7504292558  1.0665934196
##  [516] -0.3273658395  0.8604948252 -0.3844512628 -0.8343192067  1.0566244875
##  [521] -0.7677176222  0.3803582607  0.7201123606  0.9219254778 -0.3435533254
##  [526]  0.0961258699  0.8125717082  0.8451982546  0.6737201775  0.6979282027
##  [531]  2.2928401041  0.4397193245  0.4777774217  1.1085030148  1.0939556919
##  [536]  0.9128079122  0.7380719350  0.4875925530  0.9685628966  0.8686509903
##  [541]  0.3764443785  0.9782463561  0.5754524879  0.6249565129 -0.7264012917
##  [546]  0.7315490130  0.9675298388  0.0512882389 -0.3852975722  0.9066116210
##  [551]  0.5783662770  0.5787247450 -0.3276495747  1.2445937689  1.0741884966
##  [556]  0.7166899837 -1.1375697559  1.6487326744  0.6457798664 -0.6752995241
##  [561]  0.4114672210  0.1479420657  0.3307119425  0.5662321329  0.8982970042
##  [566] -0.2079166333  0.4952702131  0.8702782645  1.5859719883  0.3527806940
##  [571] -0.7460604535 -0.3201098059  1.0173699141  0.4841554525 -0.6562997877
##  [576]  0.5440075339  1.5016115701  1.1345921941  0.4804551005  0.3780042299
##  [581]  0.8454787691  0.6023783448  0.5085865095  0.8591412883  1.3418060809
##  [586] -0.0028018142 -0.2697240396  1.1448172775  0.8999971670  0.9356700927
##  [591]  0.9553078330  0.7189887933  0.7818968946 -0.3945169411  0.8869512687
##  [596]  0.5302671712  1.4716344880 -0.2901409145 -0.9588195380  0.4315141055
##  [601] -1.8494046420  0.0180608325  0.8560651888 -0.8191683043  0.7361581305
##  [606]  0.8603825609  0.9837845861  1.2722511219  0.9552756916 -0.8410101468
##  [611] -0.8510011306  0.6239480823  0.9565243393 -0.0063389864  1.0279500687
##  [616]  1.0368693960  1.7443263473  0.1530864436 -0.6630530245  1.3537044155
##  [621]  0.4649738086 -0.8050827296 -0.0040680297  0.8397433903  0.8921979429
##  [626]  1.0180415014  0.9741213848  0.0438036157  0.0922720351  0.3918647655
##  [631] -0.8198711291  0.8618769281  1.3911287837 -0.3544978846  1.4205608640
##  [636] -1.3227474965  0.6162302469  0.6692611013  0.6301028495  1.5441974919
##  [641]  0.8730281199  0.2972031820  0.7808584204 -0.8528009978  0.0501072796
##  [646] -1.3091756938  0.5834192152  0.0012895016  0.9337808242  0.0094676379
##  [651]  0.8055934231  1.1936218294 -0.7839125423  0.3934470240  1.0578973421
##  [656]  1.4323097417  0.2555677102  1.6059167919  1.1483654567  0.4500681737
##  [661]  0.8309078894  0.9000411488  0.0543170389  0.6058355377  0.6825623623
##  [666]  0.9790983494  1.3223908231  0.2884543438  1.4425100035  1.6009420553
##  [671]  1.3230938886 -0.3278500154  0.6009895173 -0.6028943190  0.6236759182
##  [676]  0.6326091140  0.3913062764  0.4735477686  0.8160337970  0.9864588548
##  [681]  0.4705724206  0.7040065378  0.0006098027  0.9380031837  0.8875469699
##  [686] -0.7344822662  0.0589532321  0.8318791176  0.8173178788  0.8923732616
##  [691]  0.9230623850  1.1562708598 -0.8059444173  0.8304760805  0.9576623777
##  [696]  0.0568044473  0.9788220961  0.1414367277  0.6709028060  1.1408383786
##  [701] -0.3938830495  1.1876637634  0.8856598445  0.8157716810  0.8841073480
##  [706] -0.0052779984  0.6334842534 -0.0025688596 -0.0106375176  1.0390687981
##  [711]  0.8871130970  0.8201644675  0.5053301746 -0.3935526827 -1.2068996951
##  [716]  0.0628208225  0.9362743227 -0.3508355594  0.5180281154  1.3037168555
##  [721]  1.3111696281  1.3980611855 -2.2199763747 -0.3937209454  0.6749671649
##  [726]  0.7701867558  0.9557289086  0.9319486660  0.1760326429  0.9830208971
##  [731]  0.4496043124  0.7639693844 -0.2658243818 -1.3111683014 -0.5589247275
##  [736]  0.6072316193  1.0026614231 -0.2212791655  0.8993152795  1.7750634838
##  [741]  0.8386122880  1.2110383625  1.3995975324  1.1260311081  0.7844029965
##  [746]  0.8357904818 -1.9683311787  0.7992240147  1.2576702318  1.1655952355
##  [751]  1.1636141066 -0.2737935932  2.4391942370  0.0168826401  0.4323467814
##  [756]  0.7474674582  0.9729389961 -1.2957133900  1.6699645437  0.6386220752
##  [761]  0.4019306508  0.7770241194  0.9877276839  1.0181960282  0.0629777564
##  [766]  1.0329068434  1.3086447206  1.0125589689  0.4600829219  0.7627482680
##  [771] -0.8413892684 -0.2695975634  2.4380769306 -0.6860645899  0.7989915638
##  [776]  0.3951274664  1.1262630856  0.7061348015  1.0424084134 -0.3689086755
##  [781]  0.3832243920  0.5915032708  0.6254510230 -0.7189293038  0.7042318793
##  [786]  0.9898206540 -0.7025341434  1.9076581800  0.8652993353  0.4193274269
##  [791] -0.2919147689  0.8581294894 -0.4087830124  0.8470228104  0.7691759546
##  [796]  0.0841899632  0.5626716220  0.6954253300  0.0086194321 -0.3636646150
##  [801]  1.1952488049  1.1786124857  1.3868993927  1.0270502597  0.1300170581
##  [806] -0.6518943112  1.0546849817  1.3629399274  0.9956724316  0.1022892765
##  [811]  0.9071298926  0.0695011693 -0.7655905648  0.6086990341  1.3016733837
##  [816]  1.1116852121  1.0646564353  0.0643850382  0.5916542696  1.5195492142
##  [821]  1.2565863492  1.4608209492 -0.6946219286  1.6741389030  0.8042570299
##  [826]  0.8469876111  0.7462315126  0.3376046003 -1.7968677141  0.8845415153
##  [831] -0.7177263551  1.1283642743  0.3521114179  1.3510063366  0.6434702489
##  [836]  0.9523535299 -0.3957385248  1.6156847055  1.2772885779  1.0356043995
##  [841]  1.1876323992  0.7338926919  1.4162224960  0.8904741645 -0.6860645899
##  [846] -0.2311709146  2.2752889841  0.4929362675  0.4053518516  1.1795462838
##  [851] -1.6217245216 -1.3726981973 -0.7789814521  1.0164573530 -0.6929027600
##  [856]  0.9705065725  1.4399034706  0.6040788008 -0.3337080434  0.6533416666
##  [861]  1.0511909233 -0.2504054372  0.8448769737  0.9353491847  1.0656089029
##  [866]  0.8545258586  0.6868930748  0.6982748695  1.0574209280  0.7973002560
##  [871]  0.4396115777  0.1208252995  1.1753702732  0.8011903044  0.2745234347
##  [876]  0.8842278093  0.9819292566 -0.6893961334  0.0333328016 -0.0114180514
##  [881] -2.1106447065  0.0120092857  0.9053699666  0.5263990881 -0.8487360386
##  [886]  0.8689216521  0.0141230781  0.5892134931 -0.0056707400  0.9740597858
##  [891]  1.2777173039  0.0620447877  0.3941512439 -0.4256166640  1.2428357677
##  [896]  1.0595393484  0.6905429685 -0.0501997889  0.4128841115  1.0790908816
##  [901]  0.7777610202  0.7646397653 -0.7384790460  0.6417864781  1.1973059885
##  [906]  0.7978031644  1.0913776387  0.8187338571  1.1494270547  2.0231508004
##  [911]  1.0359893178  0.9039677178  0.3811629330  0.6867089520  0.5479946821
##  [916] -1.4232416201  0.4445127328 -0.3275641603  0.5880540112  1.2905737029
##  [921] -2.3903367159  0.8849941815 -0.6191201499  1.0209075173  0.7669810494
##  [926] -0.4115558637  0.6242471780  0.6176577455  1.1747382624 -0.3497509763
##  [931]  0.4578989381  0.8671980140  0.9475908120 -0.3119361121  0.8262422302
##  [936]  0.0485599976  0.8787547966  1.3760552113  0.7539460978  1.1230011972
##  [941] -0.7654941147  0.8591888829  0.7994791567  0.0665476338  0.0984056285
##  [946]  0.9204068533  1.8899873724  0.1008547669 -1.2613012349  0.8289736162
##  [951]  0.7509259489  0.0011941665  0.7643209119  0.8159022711 -0.7109795281
##  [956]  0.6173965479  0.6293587231  1.4793525543  0.5912709611 -1.4490665017
##  [961]  0.1813712227 -1.2264807601  1.2336837301  1.0900705592  0.0090329833
##  [966]  0.7474971404  0.3823269648  1.0711756551  0.9957710404  1.4850089806
##  [971]  0.7804307701 -0.3206860921  0.7289382526  0.7938709556  0.9992874666
##  [976]  0.0630312234  0.3659722017  0.9091598546  2.1513358186  0.9131486451
##  [981]  1.3455159448 -0.8440682557  0.8463972722 -2.0939604838  0.6459926274
##  [986]  0.5766187001  1.5950605157  0.4155160818  0.9682283477 -0.0027880269
##  [991] -1.1493692006  1.1056255186  1.7390973332 -0.2601275167  0.9165830631
##  [996] -0.8509140187 -0.0090033676 -0.2110396493  0.8576231268  0.7958902395
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
##   0.69893970   0.36447409 
##  (0.11525683) (0.08149757)
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
## [1] -0.8385366  0.2978760 -0.2203714  0.7182985 -0.2837672 -0.9118182
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
## [1] -0.0095
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.935358
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
## t1*      4.5 -0.03513514    0.912624
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 1 2 4 5 6 7 8 
## 2 2 1 1 1 1 1 1
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
## [1] -0.0139
```

```r
se.boot
```

```
## [1] 0.9021811
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

