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
## 0 1 2 5 6 7 8 9 
## 1 1 1 1 1 3 1 1
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
## [1] -0.0201
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
## [1] 2.71778
```

```r
UL.boot
```

```
## [1] 6.24202
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.8   6.2
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
##    [1] 3.1 4.7 6.0 4.6 4.1 6.5 3.4 3.9 5.3 5.1 3.7 4.2 5.0 4.5 5.2 3.5 3.1 3.9
##   [19] 6.4 4.1 4.9 3.4 4.7 5.2 3.4 3.7 5.6 3.9 5.5 5.0 5.7 5.1 3.3 3.4 4.9 5.4
##   [37] 4.4 4.5 3.0 5.2 4.7 4.9 4.6 4.5 5.0 4.0 5.2 4.0 4.6 6.1 5.4 4.2 4.2 4.1
##   [55] 4.8 6.0 3.8 3.5 4.6 4.1 4.7 4.9 5.0 6.4 4.8 4.2 2.4 5.5 3.2 5.0 5.0 4.5
##   [73] 3.8 2.8 3.9 5.3 4.0 5.3 4.4 5.8 3.7 5.2 3.0 4.3 4.8 5.5 4.6 5.4 4.7 3.9
##   [91] 4.3 4.4 4.3 4.5 5.9 5.3 6.0 5.3 5.2 4.0 4.5 3.9 3.2 5.3 4.9 4.1 5.3 4.6
##  [109] 5.3 4.4 4.7 5.0 4.1 3.5 3.1 4.5 5.1 3.3 4.5 5.4 6.5 4.3 6.5 5.3 4.9 3.1
##  [127] 3.3 5.4 6.4 3.8 5.8 3.8 4.5 6.0 3.6 6.5 4.1 3.4 5.0 4.4 5.0 3.4 5.3 5.2
##  [145] 5.5 2.3 5.6 3.4 5.4 5.3 4.2 5.3 6.2 5.1 4.1 5.5 5.7 4.2 3.1 3.6 3.5 3.3
##  [163] 3.8 4.1 5.2 4.5 4.0 3.4 3.7 6.4 3.3 4.8 5.0 4.8 3.7 4.0 5.4 5.4 5.0 2.6
##  [181] 5.1 5.6 4.2 3.7 6.7 3.9 3.2 5.8 4.5 5.2 4.5 2.9 3.8 5.6 5.0 5.6 4.6 4.2
##  [199] 5.8 3.3 4.1 4.5 5.6 4.8 5.1 3.9 3.0 4.7 4.4 2.6 5.7 4.9 4.1 5.5 4.8 5.1
##  [217] 4.7 4.8 5.9 5.8 3.5 5.6 3.4 5.4 5.1 4.7 4.1 4.8 1.7 3.1 5.6 3.5 4.7 3.5
##  [235] 3.8 3.6 5.2 5.6 4.5 4.8 4.7 4.8 5.1 5.1 4.6 4.7 5.7 4.4 5.0 3.8 3.4 3.5
##  [253] 3.2 4.1 5.1 4.8 4.3 3.8 5.4 3.7 4.4 4.4 3.0 4.8 4.7 4.0 2.7 3.9 4.9 4.9
##  [271] 4.2 4.2 4.9 4.5 4.0 4.6 3.5 4.7 4.7 5.8 3.7 4.5 4.2 5.9 5.1 5.2 4.6 3.3
##  [289] 3.9 3.9 5.2 5.8 4.7 4.9 4.0 4.0 2.6 3.7 4.1 3.6 3.8 3.4 6.7 5.0 4.4 4.3
##  [307] 4.2 4.2 4.8 4.7 5.2 6.8 3.9 5.4 4.2 5.6 2.7 4.4 5.5 5.4 5.7 5.1 2.9 2.9
##  [325] 4.6 3.4 4.0 5.0 6.5 4.4 2.7 4.1 5.1 5.2 5.3 6.3 3.6 3.8 2.4 4.1 4.4 5.3
##  [343] 5.1 5.1 5.0 5.8 2.6 4.3 4.3 3.0 5.3 5.7 4.1 3.3 5.5 4.2 4.2 3.8 3.8 3.6
##  [361] 5.9 2.8 3.8 3.5 5.0 3.3 4.3 3.4 5.4 5.5 5.2 4.8 4.7 2.7 4.6 2.8 5.3 4.9
##  [379] 5.0 3.1 4.8 5.8 5.5 4.5 4.1 4.1 3.8 3.6 6.0 4.4 4.5 4.1 4.7 4.2 4.3 3.3
##  [397] 4.4 4.6 3.7 3.6 5.3 4.3 4.9 6.2 4.7 5.7 3.5 3.9 5.0 4.1 4.3 4.4 5.3 4.8
##  [415] 4.5 3.8 5.2 3.5 5.5 6.4 3.2 3.2 4.0 3.2 4.4 4.3 3.1 3.7 5.2 4.9 4.4 3.1
##  [433] 4.7 4.7 4.6 5.2 5.5 6.5 4.0 3.7 2.7 4.4 3.3 4.8 5.7 4.7 4.2 5.4 3.8 3.7
##  [451] 4.4 5.9 4.0 5.3 4.4 6.1 6.1 4.9 4.7 5.2 4.7 4.3 5.3 4.4 4.5 3.8 5.5 5.1
##  [469] 5.2 4.5 2.1 4.0 3.0 5.3 4.9 3.2 4.6 5.3 5.0 5.6 3.8 4.8 4.1 5.2 4.3 5.7
##  [487] 5.4 4.6 4.6 5.0 4.9 4.9 5.4 5.2 2.4 3.3 5.4 3.0 4.8 4.5 3.7 4.9 5.1 3.0
##  [505] 4.3 4.6 7.4 5.3 5.2 4.7 2.4 4.2 4.7 4.1 4.5 3.8 5.5 4.9 6.4 4.5 4.4 5.3
##  [523] 3.6 3.3 4.8 3.2 5.3 3.9 5.3 4.1 5.1 4.3 3.0 4.4 4.9 5.8 4.7 3.3 4.9 4.9
##  [541] 5.4 4.0 5.5 4.4 4.0 2.9 4.1 4.6 4.0 3.7 2.7 5.6 3.4 3.4 4.2 2.0 4.3 5.2
##  [559] 3.8 4.4 3.2 5.0 4.1 3.5 4.4 3.7 3.2 3.8 3.9 4.7 6.2 4.7 6.2 6.5 5.5 3.9
##  [577] 6.0 3.8 5.0 5.4 4.3 6.1 4.2 4.6 5.1 3.4 4.5 3.8 3.3 4.0 6.9 4.3 3.9 5.1
##  [595] 5.9 3.7 5.7 3.4 5.3 5.5 4.2 3.8 2.0 4.5 3.0 5.0 3.0 4.6 4.1 4.4 5.6 5.7
##  [613] 3.4 4.2 5.6 4.8 4.6 6.7 4.1 4.9 4.8 5.5 4.2 4.5 4.7 3.7 3.6 4.7 5.7 5.3
##  [631] 4.2 6.6 5.5 5.1 4.6 4.6 4.9 2.8 5.2 4.1 4.7 3.5 4.9 5.6 5.4 5.8 4.8 4.5
##  [649] 4.2 3.6 4.8 3.8 3.3 3.5 5.1 4.6 3.6 5.1 4.1 3.4 4.8 3.7 5.9 6.3 3.9 4.0
##  [667] 4.6 3.9 4.5 6.1 3.9 3.8 3.5 4.8 4.8 4.3 3.2 4.1 4.6 4.1 5.9 3.6 5.7 4.8
##  [685] 2.4 4.2 5.5 5.3 5.3 5.4 5.1 4.0 4.2 7.1 4.5 5.7 4.8 3.6 4.9 5.7 5.5 4.3
##  [703] 5.4 3.5 2.8 4.2 3.7 4.6 4.1 3.7 4.3 2.9 4.0 4.5 4.8 4.3 4.7 4.6 6.5 4.7
##  [721] 3.7 3.6 2.3 6.2 5.4 5.5 5.7 5.6 3.3 4.1 4.7 3.8 3.8 6.0 5.1 5.1 7.0 3.7
##  [739] 2.8 3.7 4.6 4.8 4.2 5.5 2.1 3.3 5.1 3.9 6.4 5.8 3.8 4.5 3.2 4.3 1.9 5.3
##  [757] 4.4 5.2 4.9 5.5 4.9 3.9 4.4 3.6 3.2 4.1 5.0 6.8 5.8 4.1 4.1 5.2 4.5 4.5
##  [775] 5.1 3.7 6.1 4.0 4.2 3.7 4.6 4.2 4.9 4.3 5.1 4.6 6.5 3.4 3.9 5.5 4.2 4.1
##  [793] 4.5 3.1 5.4 4.7 4.9 4.0 4.5 5.4 4.3 3.6 5.4 2.7 5.2 5.7 6.1 3.1 3.7 3.4
##  [811] 3.2 5.2 6.0 5.1 4.3 4.5 5.6 4.8 2.6 5.0 3.2 5.1 5.4 5.8 3.7 3.9 4.9 2.9
##  [829] 4.8 4.1 4.7 5.6 4.4 5.8 4.1 3.3 4.1 6.0 3.3 4.2 5.0 3.6 5.6 3.9 5.7 4.9
##  [847] 4.0 5.3 3.3 4.0 3.5 4.9 4.8 3.9 6.0 4.7 3.4 4.9 4.5 5.9 2.9 4.0 6.2 3.8
##  [865] 4.8 3.1 3.4 4.0 2.0 4.1 3.2 4.9 5.7 4.6 4.0 6.0 4.5 4.6 6.2 4.7 5.7 5.4
##  [883] 3.9 5.8 3.6 5.0 3.2 4.3 2.5 4.3 4.7 4.8 3.9 4.3 6.0 5.5 3.5 4.2 3.7 3.9
##  [901] 4.1 5.6 4.0 5.4 4.0 4.8 4.6 5.5 4.1 5.3 4.4 4.7 4.8 4.1 4.1 3.6 4.3 5.3
##  [919] 4.2 4.3 4.6 4.6 4.8 4.3 5.2 5.9 3.7 5.9 4.8 2.3 3.8 4.7 4.7 3.5 4.7 3.3
##  [937] 2.6 4.7 3.9 4.5 5.5 3.4 2.7 3.1 2.7 4.6 4.1 4.6 5.6 1.7 3.7 6.3 4.7 5.3
##  [955] 4.5 4.3 3.2 3.7 4.9 7.2 3.5 4.1 4.4 4.7 4.3 3.8 4.4 6.1 4.1 3.9 5.7 4.3
##  [973] 4.8 4.9 2.2 2.7 3.0 3.0 2.9 5.3 4.0 5.5 4.1 4.9 4.9 4.5 5.6 3.4 3.4 4.9
##  [991] 5.8 5.7 6.0 3.7 4.7 3.2 5.0 5.7 5.0 5.5
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
##    [1] 4.7 3.1 4.2 3.8 5.2 4.4 4.4 3.8 2.9 4.6 4.2 6.7 4.1 3.5 6.6 4.5 5.0 4.8
##   [19] 3.8 5.3 4.6 3.3 5.1 6.5 4.2 3.8 4.8 5.6 3.8 3.7 4.8 3.2 3.0 4.5 4.2 5.0
##   [37] 4.0 5.1 6.7 4.1 4.7 4.4 3.8 4.7 5.4 5.8 6.0 3.3 2.9 6.4 4.4 4.3 4.1 4.2
##   [55] 4.5 4.0 4.6 5.5 5.1 6.8 3.0 3.2 3.3 6.2 5.6 6.1 4.1 4.7 6.5 5.4 4.4 4.6
##   [73] 4.9 5.0 5.1 4.6 3.9 3.7 5.1 4.6 3.7 5.6 4.0 4.2 5.2 4.2 6.9 5.0 5.1 5.5
##   [91] 5.1 5.4 4.1 4.6 2.3 5.8 5.0 5.4 4.7 5.0 3.8 5.5 3.4 3.4 4.2 4.3 4.3 5.5
##  [109] 4.2 4.8 4.7 3.8 4.7 4.5 3.9 5.6 4.2 4.5 4.1 5.9 4.4 3.8 4.3 3.2 3.0 5.5
##  [127] 4.8 4.9 4.8 3.2 4.3 4.9 5.3 6.6 4.1 6.0 2.3 5.1 4.5 4.6 4.9 5.3 5.2 4.7
##  [145] 4.5 4.0 4.7 3.4 6.3 3.8 6.1 3.8 5.3 5.3 4.5 4.8 4.9 5.1 2.3 5.0 2.9 4.0
##  [163] 5.7 4.3 3.9 4.5 4.1 3.0 6.9 4.4 6.1 5.1 5.9 3.9 4.7 3.5 4.7 5.9 4.7 3.7
##  [181] 3.4 4.7 4.6 4.6 5.5 5.2 4.2 3.4 5.9 3.8 3.6 4.7 4.8 5.1 2.9 5.0 4.4 5.0
##  [199] 5.2 5.0 5.3 4.0 4.8 6.8 4.7 4.8 4.9 5.2 5.4 4.2 4.4 4.6 4.5 5.3 4.4 4.0
##  [217] 3.7 5.5 4.5 6.4 5.6 3.2 4.6 3.9 3.9 3.5 4.7 6.2 3.8 3.8 5.1 4.0 3.0 4.5
##  [235] 3.7 4.2 6.3 6.0 5.7 3.9 5.1 4.8 5.0 4.1 3.3 5.3 4.4 5.2 2.9 4.8 4.5 4.9
##  [253] 3.7 7.0 4.1 5.7 5.2 4.2 5.5 4.2 4.6 2.5 5.0 5.3 3.3 4.8 4.0 6.6 5.0 4.2
##  [271] 4.7 3.0 5.7 4.0 4.5 5.8 5.6 6.9 4.9 3.7 4.8 4.0 5.4 3.8 5.3 5.6 4.5 4.9
##  [289] 5.4 4.5 3.9 3.7 4.6 2.7 4.2 5.2 5.0 4.5 6.2 4.3 2.6 4.0 2.9 5.1 5.1 3.6
##  [307] 5.2 5.3 3.9 4.8 3.7 4.4 5.2 3.6 5.0 3.0 3.1 3.4 4.5 4.7 3.6 5.7 5.5 4.6
##  [325] 4.5 3.6 4.5 4.7 5.9 4.4 3.7 3.2 4.5 4.2 4.7 5.0 5.2 4.5 4.3 3.9 5.0 5.5
##  [343] 3.7 4.2 4.4 4.7 5.2 5.9 5.0 3.8 4.3 5.8 3.7 3.8 4.3 3.9 3.8 4.5 5.0 5.4
##  [361] 4.2 3.2 4.6 3.1 4.7 3.7 5.6 5.6 4.3 3.4 3.7 3.9 4.8 6.4 3.4 5.8 4.2 4.9
##  [379] 3.2 5.2 3.5 3.9 3.9 4.4 4.7 6.2 4.3 3.3 4.1 4.3 4.4 4.6 4.8 4.0 5.2 3.7
##  [397] 4.9 3.6 5.6 3.3 4.3 6.8 4.8 5.1 4.6 4.5 6.1 3.5 5.0 1.8 5.4 5.9 6.1 4.6
##  [415] 5.3 4.1 3.5 3.6 5.7 4.7 4.5 3.5 3.2 4.2 4.2 5.7 4.7 5.1 3.4 4.4 3.6 4.7
##  [433] 3.9 3.2 3.1 3.7 4.3 4.2 4.2 4.6 4.8 3.4 5.3 4.9 3.8 4.9 2.9 4.3 4.8 4.4
##  [451] 4.7 4.5 4.6 4.7 4.9 5.2 4.3 4.6 4.8 4.6 5.5 3.7 3.8 5.2 4.3 5.7 5.8 3.7
##  [469] 5.6 4.0 3.6 5.1 4.4 2.1 6.1 4.3 4.7 4.6 4.6 3.3 5.4 4.7 5.7 4.7 4.6 6.2
##  [487] 4.7 4.7 5.7 4.9 4.4 4.6 3.3 4.5 4.7 5.9 6.7 3.6 5.5 4.3 5.1 2.2 4.0 3.3
##  [505] 4.2 4.7 4.5 3.4 4.2 4.9 2.8 3.4 4.9 4.8 5.9 4.6 3.6 5.0 4.0 5.5 3.1 5.2
##  [523] 3.5 4.8 4.1 6.3 4.0 3.5 3.5 4.5 4.8 4.6 3.0 3.9 4.5 5.0 5.3 5.0 5.1 4.2
##  [541] 4.4 5.3 3.4 5.4 4.8 4.3 4.8 4.2 5.8 3.6 4.6 3.1 4.0 4.4 4.5 3.0 5.8 5.3
##  [559] 4.7 4.7 5.0 4.7 4.4 5.4 4.8 6.6 6.1 6.7 4.3 2.2 6.6 3.4 4.2 5.6 4.1 4.5
##  [577] 5.2 3.9 4.0 5.0 4.0 4.8 5.8 3.9 4.6 3.0 4.6 5.2 4.0 4.7 5.1 4.4 4.0 3.5
##  [595] 4.9 3.9 2.0 4.2 6.2 4.9 5.1 4.5 6.9 3.0 4.2 4.0 6.7 4.7 4.0 2.8 2.4 2.5
##  [613] 3.8 4.0 3.1 4.7 4.3 3.2 4.1 4.3 5.5 4.3 4.4 3.4 4.6 4.8 5.7 5.6 4.9 5.3
##  [631] 5.5 5.0 6.2 4.5 5.5 3.8 4.3 2.7 4.2 4.5 4.8 4.1 3.9 4.3 4.9 4.7 4.1 5.1
##  [649] 6.0 5.3 3.2 5.0 3.7 2.4 4.8 5.3 4.9 4.4 3.9 4.6 5.4 5.4 5.0 4.4 3.9 4.2
##  [667] 5.7 3.6 4.5 4.7 4.4 4.8 5.2 5.1 5.6 4.3 5.6 4.8 2.4 5.5 4.9 4.3 5.8 4.8
##  [685] 5.0 4.1 4.4 4.4 4.3 4.9 5.2 4.9 6.1 1.9 5.0 3.3 4.2 4.7 3.4 4.7 6.1 4.7
##  [703] 3.9 4.7 3.2 4.8 5.1 4.6 3.8 3.7 4.0 4.2 3.3 4.8 4.3 5.3 5.3 4.4 2.3 4.7
##  [721] 4.0 4.3 5.6 4.5 5.5 4.7 4.1 5.4 4.1 4.2 4.2 4.4 5.7 4.6 5.0 5.2 2.7 4.0
##  [739] 3.8 6.3 5.4 4.1 5.6 4.2 4.9 4.5 3.8 5.6 6.6 3.9 4.5 3.5 2.3 5.2 4.6 3.8
##  [757] 4.5 4.6 4.2 5.4 6.2 4.6 4.3 4.4 4.1 5.8 4.7 3.9 3.8 4.6 5.4 5.8 4.1 4.5
##  [775] 4.5 2.7 5.1 5.9 5.7 4.3 3.1 4.9 3.7 2.0 4.3 3.8 6.2 3.5 5.1 3.1 5.1 4.6
##  [793] 2.0 3.9 5.2 4.3 4.5 4.7 3.4 4.7 4.6 4.3 4.6 4.0 5.3 5.9 3.3 2.2 4.3 5.2
##  [811] 3.6 4.2 5.5 3.7 3.7 4.3 4.3 4.4 6.0 3.2 4.1 4.0 3.1 4.7 3.4 4.1 4.5 4.8
##  [829] 5.2 4.5 3.9 4.5 4.1 3.9 4.7 5.2 4.6 6.4 3.4 4.3 2.0 4.4 4.1 3.7 4.2 6.0
##  [847] 4.2 4.3 2.8 5.7 4.0 3.1 5.6 3.2 4.6 3.5 6.8 4.5 4.9 5.1 6.7 4.6 4.8 5.4
##  [865] 5.4 5.0 4.7 4.4 4.3 5.3 6.1 3.2 4.5 3.5 4.5 3.7 4.8 2.4 3.5 4.8 5.5 5.0
##  [883] 3.8 4.5 5.7 4.8 4.1 4.4 3.2 3.6 4.1 4.9 3.8 2.6 2.9 5.8 5.1 3.7 4.7 5.0
##  [901] 5.1 5.6 5.6 3.7 3.2 5.5 5.2 5.8 5.2 5.9 4.1 4.9 4.5 4.3 4.8 5.8 3.8 3.8
##  [919] 4.6 3.3 4.3 5.2 4.2 4.7 5.1 5.8 4.2 5.0 3.8 3.4 4.1 4.4 5.7 2.3 3.9 5.5
##  [937] 5.0 4.9 3.5 3.9 5.9 5.6 4.1 5.8 4.4 6.4 5.6 4.9 3.6 4.1 5.1 4.8 4.1 6.2
##  [955] 4.9 4.6 4.4 4.5 5.4 4.9 6.3 3.0 5.6 4.1 5.3 5.9 5.2 3.3 4.6 3.9 6.4 4.2
##  [973] 4.3 5.5 5.2 5.6 4.3 5.1 3.6 4.4 3.5 2.8 5.0 6.0 5.4 3.7 4.6 3.9 6.1 4.4
##  [991] 3.6 5.5 4.0 4.1 3.7 4.5 5.6 4.6 3.9 3.9
## 
## $func.thetastar
## [1] 0.0382
## 
## $jack.boot.val
##  [1]  0.4923706  0.4031250  0.2633238  0.1979290  0.1695015 -0.0462141
##  [7] -0.1759104 -0.2005731 -0.4151335 -0.4707317
## 
## $jack.boot.se
## [1] 0.9509005
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
##    [1] 4.8 4.3 4.2 4.4 3.3 5.0 3.7 3.9 4.2 3.5 4.1 5.1 3.4 5.3 5.2 4.0 4.8 3.1
##   [19] 3.7 5.6 5.1 4.0 3.7 4.3 5.1 3.9 4.4 5.5 2.6 4.1 4.6 5.5 3.6 4.4 3.9 3.0
##   [37] 4.8 2.9 5.0 6.2 3.3 3.8 4.8 4.9 5.2 3.4 3.2 1.6 4.2 3.6 2.4 5.0 5.3 4.6
##   [55] 5.6 4.5 5.1 5.7 4.4 3.4 5.0 4.8 3.7 4.3 3.2 3.4 3.7 4.0 3.6 3.7 3.8 4.6
##   [73] 4.0 5.5 3.5 3.3 6.0 5.0 4.1 6.2 4.1 4.0 2.5 5.4 4.1 4.2 5.4 4.2 3.4 4.4
##   [91] 6.4 4.8 5.6 5.3 4.1 4.7 5.0 3.6 5.3 4.4 3.7 5.0 2.8 3.9 4.6 5.0 6.1 4.6
##  [109] 5.4 4.1 6.0 4.5 4.7 5.4 4.5 5.3 4.2 5.7 7.1 3.9 3.9 7.2 3.8 4.2 4.9 2.4
##  [127] 4.1 5.0 3.7 5.1 2.0 4.1 4.9 5.6 6.1 5.1 4.1 4.0 3.7 3.9 5.0 5.8 5.1 5.0
##  [145] 5.3 6.2 2.9 5.3 4.5 5.8 3.9 3.7 4.6 3.5 3.7 6.6 3.1 5.9 5.5 4.0 5.0 6.2
##  [163] 3.5 4.4 5.3 4.1 3.7 3.6 3.7 5.3 5.3 2.2 5.8 2.7 4.4 4.9 2.6 4.6 4.8 3.5
##  [181] 3.6 4.0 4.1 6.2 5.1 4.6 5.9 2.2 4.5 3.2 5.3 5.3 5.4 5.1 3.8 4.9 4.4 3.3
##  [199] 5.1 3.4 4.8 3.9 5.0 5.3 4.5 4.5 5.8 3.7 5.2 6.4 4.7 5.0 4.7 2.2 3.2 4.0
##  [217] 4.7 4.4 4.4 3.8 2.7 4.7 4.5 5.6 4.0 4.1 4.6 3.9 3.2 4.5 4.6 4.3 4.6 5.1
##  [235] 3.8 4.9 5.4 4.6 6.0 6.2 3.4 3.9 3.6 3.7 2.8 5.6 4.3 3.7 5.1 5.6 5.8 4.7
##  [253] 5.3 4.1 4.9 5.6 3.5 3.3 5.6 6.2 4.9 3.4 3.8 5.6 3.7 5.0 5.4 5.4 3.8 5.7
##  [271] 4.2 7.3 5.0 4.1 5.1 5.5 4.8 2.8 5.4 5.1 2.9 2.7 4.3 6.0 4.0 4.2 3.2 5.0
##  [289] 3.7 4.6 4.0 3.0 4.2 4.8 4.5 5.0 6.1 4.4 4.8 5.7 5.2 3.4 3.5 6.4 4.3 4.0
##  [307] 3.9 4.0 5.3 4.0 4.9 5.3 4.4 3.3 4.6 2.2 5.5 5.1 5.2 4.2 6.5 4.8 3.6 3.5
##  [325] 5.4 6.1 5.2 3.2 3.1 4.4 3.1 5.2 4.8 4.0 4.2 6.2 4.8 4.6 3.8 4.9 3.3 4.5
##  [343] 5.5 4.8 3.8 5.5 3.3 4.1 4.3 4.3 4.1 4.7 5.4 3.7 3.3 3.9 5.1 5.7 4.7 4.9
##  [361] 4.4 4.0 6.1 2.7 4.8 5.8 3.5 5.2 5.7 4.8 3.5 1.8 4.0 5.0 4.1 4.8 5.5 4.8
##  [379] 4.4 3.1 5.2 3.6 4.9 4.5 4.2 4.7 5.5 3.6 2.7 5.2 5.8 4.5 4.2 6.9 3.1 3.6
##  [397] 4.8 3.3 5.5 4.5 4.6 4.6 3.1 5.5 3.6 5.0 3.9 3.8 5.6 3.1 3.6 4.2 6.7 3.3
##  [415] 3.0 3.5 4.0 5.6 3.8 4.0 5.0 5.5 4.6 2.9 4.8 5.0 5.2 6.1 5.2 5.4 4.0 4.8
##  [433] 3.9 5.2 3.3 6.5 5.4 3.2 3.8 4.7 4.9 5.5 5.1 3.9 4.7 4.4 6.5 3.4 4.9 3.8
##  [451] 5.2 3.6 3.9 4.4 4.4 3.6 5.9 4.1 5.5 5.4 2.9 4.7 2.0 4.5 5.7 3.8 4.0 5.2
##  [469] 5.4 4.2 4.5 4.4 3.6 4.5 3.8 5.8 3.9 5.4 3.7 4.4 6.4 4.1 4.7 5.9 4.6 2.9
##  [487] 4.9 4.9 4.8 4.4 4.4 4.3 6.3 3.8 6.0 5.6 5.1 5.0 4.2 4.7 3.7 5.4 3.6 4.0
##  [505] 3.9 4.6 5.7 5.3 4.5 3.8 4.5 4.8 5.7 4.3 3.4 3.7 4.3 3.8 4.1 3.5 4.9 4.6
##  [523] 4.5 4.0 4.5 5.3 3.4 4.9 5.7 5.4 6.0 5.2 4.9 5.0 3.2 4.2 6.3 5.1 3.0 4.1
##  [541] 4.6 3.1 3.9 4.9 4.0 5.3 4.6 4.7 5.1 4.2 4.6 4.0 4.2 3.8 4.7 4.2 4.9 5.0
##  [559] 5.6 3.2 4.5 4.1 4.2 3.9 3.0 5.4 4.6 3.8 4.2 4.4 5.6 5.1 6.3 4.5 3.3 3.3
##  [577] 6.0 3.2 3.1 4.5 4.5 2.4 3.7 5.6 4.1 5.1 7.2 5.5 5.4 2.8 4.9 4.6 5.4 4.3
##  [595] 6.0 4.5 5.9 4.8 5.1 5.3 4.3 4.5 4.4 6.0 2.7 4.4 3.8 5.7 6.3 4.3 4.9 3.2
##  [613] 5.7 5.1 4.0 3.7 4.4 5.1 4.0 4.1 4.8 4.5 3.9 5.2 5.9 3.4 4.1 4.0 3.6 5.3
##  [631] 3.7 4.5 4.4 4.3 2.7 4.6 3.6 4.6 4.3 3.9 4.2 2.6 5.6 4.1 5.6 4.1 5.9 4.9
##  [649] 4.9 4.6 1.9 2.9 5.2 5.2 5.5 5.0 5.3 4.4 5.8 5.4 5.1 3.6 4.0 5.2 5.8 4.3
##  [667] 5.6 5.4 5.1 5.2 4.6 4.6 4.2 5.8 3.1 4.4 5.3 3.4 2.3 5.6 4.5 5.0 4.7 4.9
##  [685] 4.9 4.6 4.9 5.0 4.3 3.2 4.1 3.9 5.8 5.6 5.4 4.1 6.3 5.8 3.0 5.2 2.8 4.8
##  [703] 3.7 5.1 4.8 4.9 4.0 4.2 3.8 4.5 6.0 4.1 5.1 3.6 5.0 3.2 3.3 4.1 3.1 4.1
##  [721] 3.2 4.5 4.9 3.1 4.2 4.1 4.6 5.6 4.8 4.2 3.6 4.1 4.6 4.0 5.4 4.0 5.3 3.2
##  [739] 5.2 4.4 3.6 4.2 2.7 3.9 6.7 5.5 4.0 5.7 7.1 5.1 4.5 3.9 3.2 3.6 5.6 5.1
##  [757] 5.5 5.1 3.7 3.1 2.5 4.4 5.6 3.6 3.7 4.6 5.1 3.6 5.8 3.8 5.2 4.4 3.7 5.6
##  [775] 4.8 5.1 4.6 4.7 6.5 4.0 4.8 4.0 2.8 5.1 4.0 2.7 3.9 5.1 4.2 4.9 5.2 4.8
##  [793] 4.6 3.5 5.1 4.7 2.9 5.9 4.9 5.5 1.9 3.9 5.3 4.9 4.9 2.1 5.5 4.4 2.7 4.7
##  [811] 4.5 6.8 5.5 5.1 4.7 6.1 4.3 4.6 3.6 4.4 5.1 4.7 5.6 4.1 4.9 4.9 3.2 5.3
##  [829] 4.0 4.7 3.6 4.1 4.9 5.6 4.3 3.5 3.1 4.2 3.8 4.4 4.8 4.7 5.7 4.0 3.6 3.1
##  [847] 6.1 3.9 4.4 3.8 3.2 5.2 6.0 4.2 4.1 3.2 4.6 4.3 4.9 4.3 4.4 4.9 4.9 3.5
##  [865] 3.6 5.7 5.5 4.7 3.6 5.7 5.0 3.3 4.7 4.2 5.0 4.9 4.5 3.7 2.4 5.1 3.3 4.2
##  [883] 4.7 5.1 5.3 4.4 4.3 4.5 4.6 4.3 4.0 5.5 3.9 5.1 2.9 3.6 4.0 4.8 6.4 3.3
##  [901] 3.7 3.6 5.0 4.0 3.5 4.2 6.1 2.9 3.1 3.2 5.3 3.6 5.6 4.6 4.9 4.0 2.0 4.2
##  [919] 5.1 3.8 4.8 6.3 5.2 5.4 3.8 4.0 4.0 2.0 4.5 4.1 4.0 3.0 3.1 5.7 4.2 5.5
##  [937] 4.8 4.4 2.9 4.7 5.2 4.4 5.4 5.7 4.7 4.3 5.2 3.3 4.8 4.9 3.9 3.5 5.3 4.2
##  [955] 4.1 6.0 5.1 5.2 3.9 3.6 5.5 4.1 5.4 2.8 5.6 4.4 3.3 4.9 6.6 4.6 5.7 4.0
##  [973] 4.5 5.6 3.4 4.8 5.3 4.3 4.1 5.7 3.6 3.7 3.8 5.1 3.2 4.3 3.4 5.5 3.7 4.4
##  [991] 4.8 5.8 3.7 3.5 4.4 4.1 4.0 3.3 5.1 3.9
## 
## $func.thetastar
## 72% 
## 5.1 
## 
## $jack.boot.val
##  [1] 5.50 5.40 5.40 5.30 5.10 5.08 4.90 4.70 4.70 4.40
## 
## $jack.boot.se
## [1] 1.040684
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
## [1] 0.09718536
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
##   3.871240   6.583455 
##  (1.662068) (3.018082)
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
## [1]  0.7078051 -0.1579458  0.5104537  1.6263846  0.7387317  1.7482810
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
##    [1] -0.3133952487  0.0839485441 -0.0233344217 -0.0294810102  0.3326931095
##    [6] -0.6390277752  0.9842185131 -0.2711625147 -0.0360951267  0.5895629405
##   [11] -0.1014038335  0.1372913319  0.0079045935  0.3651146479  0.0281907653
##   [16]  0.1338880865 -0.1880595639  0.4126359740  0.3447250584 -0.2467303477
##   [21]  0.4241268210  0.2076300373 -0.7795322166  0.0313566500  0.2156195932
##   [26]  0.7772345503 -0.2712648448  0.0063467317  0.3646301172  0.5521108517
##   [31] -0.4258218430  0.3130775901  0.5318447897 -0.4868400459 -0.0076724616
##   [36] -0.1512618399  0.0237915088 -0.2199282595  0.5740705931 -0.1155171344
##   [41]  0.7490003241 -0.4690690955 -0.5032924090 -0.1464673505  0.4436391578
##   [46]  0.2195829803 -0.0222218487 -0.5305339178 -0.3744532285 -0.1186548907
##   [51]  0.1054227976  0.4527264384 -0.1667812735 -0.9591241895 -0.3085131085
##   [56] -0.2952796833  0.6916367454 -0.6847253024 -0.0047292055  0.0889546440
##   [61]  0.8397435268 -0.4356247600 -0.2301013350 -0.1848114463  0.5629888373
##   [66] -0.6417131963 -0.1784789744 -0.0928465924  0.1379861774  0.0016441689
##   [71]  0.2620425909 -0.4269267491  0.3282027140 -0.0232060280 -0.0795206184
##   [76] -0.4481620817 -0.3210206693  0.4560927945  0.4436391213 -0.1907686214
##   [81] -0.1086439177  0.1751931558 -0.1695636523 -0.4630416349 -0.0595139513
##   [86]  0.3663177729 -0.6364230920 -0.1855195135 -0.3215856610 -0.0167093958
##   [91]  0.2621281875 -0.0108635433 -0.3616534896 -0.3262165547  0.6864489689
##   [96] -0.2187854144  0.1130117423  0.1504408793  0.1379861774  0.3652498608
##  [101] -0.4382164224  0.0468281296 -0.2080213589 -0.2983036675  0.2399986496
##  [106]  0.2850602952  0.2446976007  0.3583413514  0.2080899579  0.1182546837
##  [111]  0.1042394667  0.0358819860  0.6697504261 -0.2589711991  0.5572775504
##  [116] -0.2698843372  0.2658755253 -0.2508512090  0.6389388191  0.2806279213
##  [121]  0.1048753077  0.6965366643  0.3010994605  0.2456837023  0.3268697680
##  [126]  0.4465613378  0.4155062360 -0.4481404605 -0.3644210773 -0.0347677174
##  [131] -0.0173608667 -0.0284186174 -0.2740459311  0.0590519599 -0.0504859886
##  [136] -0.5264871422 -0.1082808267  0.0645213988  1.0487787419  0.4787176581
##  [141]  0.5284484750  0.3275652174  0.0055479533 -0.2322413398 -0.0938265988
##  [146]  0.8445609662  0.0246821177  0.3482376873 -0.1051951862 -0.2728480166
##  [151]  2.2983394832 -0.0721695537 -0.0011598478  0.0021763903  0.3833116921
##  [156] -1.0996522890 -0.1282991041  0.4527264384 -0.2035100304 -0.1051148598
##  [161]  0.0784416042 -0.0236657984 -0.3064843381 -0.2241988732  0.3454093529
##  [166]  0.6588290853 -0.1468501076  0.7640416649  0.0213778178 -0.4883763753
##  [171]  0.0586301953  0.8056682377  0.2090247977  1.1059353630 -0.7721822690
##  [176] -0.0397229414  0.1115682023 -0.1081220628  0.1523992698 -0.1253923086
##  [181]  0.1904327785 -0.4265167547 -0.1757560858  0.5000138410  0.6715738943
##  [186]  0.3156374563 -0.1827437407 -0.1099052329 -0.0721695537 -0.0237293960
##  [191] -0.6821830692  0.2370216147 -0.2277947060  0.6721861315 -0.0772264743
##  [196]  0.0016620504 -0.0371260478  0.1514837368 -0.3469161180 -0.6784495717
##  [201]  0.5479222929  0.1342552989  0.0088185167 -0.1353656550 -0.9580769977
##  [206] -0.6535708324  0.0591919154  0.0381750086  0.2814979849 -0.2720573135
##  [211]  0.1796028241 -0.2183743897  0.0136809346 -0.4573400917  0.6338365259
##  [216] -0.2162401764 -0.6765147606  0.5425835088  0.2260198651  0.1514593063
##  [221] -0.2905354334  0.7327437044 -0.1373836382 -0.0050235137 -0.1436951745
##  [226] -0.5765335300 -0.3104177529  0.1262459692 -0.0913521795  0.1539855406
##  [231] -0.1299550744  0.3644833360 -0.4446019801 -0.4349257585 -0.2430070100
##  [236]  0.0758890905  0.5656866111  0.0478627590 -0.0020979712 -0.3259303061
##  [241]  0.2866708582  0.4284996018  0.0789355751  0.2014241276  0.9112112381
##  [246] -0.5165458187  0.1419081399  0.1913616137 -0.4461459233 -0.3293149602
##  [251]  0.1059212739 -0.4856165434  0.3399821808  0.1376516669 -0.3517214199
##  [256] -0.0666683550 -0.8623097111 -0.0324217613  0.1305437365 -0.2102321368
##  [261]  0.5652518094  1.2378127929  0.7860133967 -0.8572213115  0.2818297471
##  [266]  0.5384372368 -0.3781611914  0.2831713357  0.2780852797  0.1177098589
##  [271]  0.7054547782  0.3647303028 -0.0166939686 -0.4559229866 -0.4386740394
##  [276]  0.2206088586  0.5551092268  0.3161811749  0.2295592515 -0.5914080197
##  [281] -0.4262711985  0.5878673114  0.7695309600 -0.1366137239  0.2929261550
##  [286]  0.5493716046  0.1412881263 -0.2293830849 -0.0173712621  0.1590692428
##  [291] -0.2320698131  0.1362667273  0.0707852997  0.7013889091  0.6188708983
##  [296]  0.3670471085 -0.0681323321  0.1136844514 -0.2569321288  0.2513974146
##  [301]  0.0737147104  0.0753761646 -0.1086265851 -0.1392072219  0.0313566500
##  [306] -0.0620306072  0.1325006121  0.0408081256  0.0617247387  0.0897486021
##  [311]  0.6276831338  0.7617448834 -0.3842639852  0.3267460003 -0.1074910413
##  [316] -0.7957953068  0.3781676400  0.5031531965 -0.3184724829 -0.0465437453
##  [321] -0.6410910925 -0.0689937019 -0.0055186289 -0.7913300528 -0.0658071000
##  [326] -0.0596008297 -0.9075899797 -0.2630210486  0.5430061253  0.6721861315
##  [331]  0.5177217161  0.1900484309  0.1352335076 -0.0887222789  0.2920203232
##  [336] -0.5686549461 -0.0767490585 -0.1808217231 -0.0160136730  1.5489183147
##  [341] -1.2876899449 -0.4137949814  0.2976920186  0.4648674710  0.1305934744
##  [346] -0.1364801738  0.4189137548  0.5697473827 -0.4499970940 -0.0447338572
##  [351] -0.1775056726 -0.4450029878  1.0615102420 -0.0315185404  0.4576381155
##  [356] -0.6017530535 -0.3850209565  0.7131722939 -0.1870125031  0.3326889719
##  [361] -0.3162385836 -0.0946229287  0.2289084816  0.6523905213  1.0072410888
##  [366] -0.7164415334 -0.6633213207  0.8248782123 -0.1511194647  0.0405299412
##  [371] -0.9522539925 -0.0860968329  0.7485896474  1.0924399879  0.3338658762
##  [376] -0.0191417078  0.3999676354 -0.4419378332  0.2782873149 -0.0315185404
##  [381]  0.1762617650  0.0644848222  1.0414958434  0.7724797527 -0.4172212974
##  [386]  0.6561859591 -0.0473298983 -0.1890373408 -0.2517054061 -0.5855422521
##  [391]  0.7689224114  0.1395387920  0.6026061100  0.0975270096  0.1379861774
##  [396]  0.5543248284  0.5628616227  0.1856614692  0.5214806269 -0.4246480701
##  [401]  0.9468656838  0.2560229347  1.4961004385 -0.1958584275  0.3767709199
##  [406]  0.4413770272  0.9555780220 -0.5324534609  0.2180426383  0.2905460709
##  [411]  0.3164531158  0.4911323494  0.4454023060 -0.1184048698  0.1419456974
##  [416]  1.2013951627  0.2870401808 -0.8964850471  0.3543162253 -0.7541495249
##  [421] -0.0109246163  0.1502778697  0.0405217630  0.6729682908  0.1317208203
##  [426]  1.1353409989  0.2019469280  0.7541526095 -0.6460456306  0.3572909724
##  [431]  0.1476169333 -0.9334482337  0.3862262531 -0.2469580500 -0.1306073878
##  [436] -0.1255448608  0.0565892395  0.1597269685  0.1032925928  0.0712959389
##  [441]  0.8282090421  0.8964526657  0.2763494369  0.4775589298  0.0070683868
##  [446] -0.4637866595  0.4274744031  0.2263107763  0.0057376571 -0.5639522907
##  [451] -0.4862876369  0.5509695748  0.3826837601  0.1376640959  0.3855091808
##  [456] -0.0549201913  0.1713596343 -0.0364926733  0.7500903006 -0.3883026920
##  [461]  0.4831217607 -0.4558568122 -0.2326263929  0.0183759783  0.9648888651
##  [466] -0.7541495249  0.1485835513  0.0005556318  1.0022938013 -0.6307914790
##  [471]  0.1564621213  0.3005484495  0.9065997515  0.5137186225 -0.5163055965
##  [476] -0.7955748841  0.4903341675  0.0949880361 -0.4163170456 -0.1479326414
##  [481]  0.6734628711  0.1773085932  0.1055255287  0.9392456174  0.0234198149
##  [486] -0.0504859886 -0.4949944835  0.5140395798 -0.0946484828  0.0669022936
##  [491]  0.2358658181  0.1696477361 -0.4908694958 -0.5062031124  0.1168765400
##  [496]  0.1150691301  0.1273739390  0.4232118444 -0.3309118535 -0.3357852959
##  [501] -0.4212367043  0.1428185400  0.1273602659  0.1466088062  0.9649350078
##  [506] -0.6677751247  0.8599717443  1.0793755355  0.4455307634  0.2367150749
##  [511] -0.3028779550 -0.2141858312  0.5283253878  0.2712422931  0.2927641129
##  [516]  0.0397911169 -0.3392515905 -0.2043044993 -0.3633917966 -0.7361463747
##  [521] -0.2403069830 -0.1796566150  0.1630512494 -0.4195000880  1.3673924044
##  [526]  0.0848747705 -0.3140455646  0.0256726947 -0.0851149155  0.1598545429
##  [531] -0.4204909008 -0.8717465331 -0.3162159239  0.4753103924 -0.4092852110
##  [536]  0.7708806470  0.1579573617 -1.1075464028  0.0281483098 -0.3816913511
##  [541]  0.3900209425 -0.4786617199  0.2108443086  0.7672155593  0.6321981175
##  [546]  0.3773741121  1.3579852188  1.0498732277 -0.1021305194 -0.1499513803
##  [551]  0.1605198801  0.0341487191 -0.5117280755 -0.4150154940  0.5645146383
##  [556] -0.3798114669  0.1685191458 -0.1275082785  0.1818044353  1.0831570870
##  [561]  0.1842946823  0.1100975869  0.3677566501 -0.1804267369  0.4700735574
##  [566] -0.1054090001 -0.6759494786 -0.0681192394 -0.2046955079 -0.2483112509
##  [571]  0.3454795563  0.1566486585  0.1205220460 -0.6360063170  0.3441789078
##  [576]  0.3085573210 -0.5811540387  0.6614646190 -0.6224163073 -0.0414209553
##  [581]  0.2887404934 -0.2983545557 -0.3229666881  0.6317160311 -0.1549742078
##  [586] -0.8322896723  0.0275623864 -0.2453886904 -0.6949609441  0.4806996047
##  [591]  0.3622249302 -0.6006638272  0.7641647330  0.3842969984 -1.4142573296
##  [596]  0.6969704997 -0.5349706106 -0.2469407158  0.0566629722  0.0958202554
##  [601] -1.1060475730  0.2958742351 -0.1487671026 -0.0552290339 -0.9014986859
##  [606] -0.2710981328 -0.3692206040  0.2353000981  0.2927869852  0.0079453948
##  [611] -0.0202139639 -0.8825983953 -1.2414694172  0.3971632177  0.9316852762
##  [616]  0.1183680633  0.2723568474  0.5033086182 -0.0996939351  0.5389538899
##  [621]  0.1012409889 -0.3577991898 -0.4773814716 -0.6602126778  0.5472825254
##  [626] -0.3819488692  0.0558767926  0.1075606220  0.4870084622 -0.1172686955
##  [631] -0.3768597241 -0.1115753704  0.6542179084  0.2436403696  0.4573835132
##  [636] -0.1744922967  0.1427004290  0.0489339660  0.9783235424 -0.0804604318
##  [641] -1.6836965259  0.3508113518 -0.0047292055  0.1533827718 -0.0668804011
##  [646] -0.9136748457  0.3517434272  0.6184579341  0.1551037325  0.3337282385
##  [651] -0.2061356726  0.4729305686  0.1830625053  0.9819657334  0.0534473158
##  [656]  1.0968022257 -0.4184831038 -0.2307627198 -0.0326612117  0.5499987299
##  [661]  0.3919760516 -0.2914167623 -0.0592016873 -0.9617851980  0.0477451738
##  [666]  0.1436678755  0.1130719955  0.4849715441  0.7037688771  0.0099942302
##  [671]  0.4365135846  0.7014722969  0.5275299820 -0.1555503762  0.1209815258
##  [676]  0.1818641860  0.3833495747  0.7677796124 -0.3078485826  0.0770658757
##  [681]  0.3155024118  0.2040165266 -0.0973034606 -0.7729800760 -0.8118331260
##  [686] -1.6738464037 -0.1122585090 -0.0817754132 -0.6853865935  0.0949219154
##  [691]  0.1040317783 -0.3067936500  0.3119519135 -1.0902690800  0.6658503896
##  [696] -0.2941929721 -0.7661119503  0.3481930098 -0.8514074128  0.3999693382
##  [701]  0.2126978004  0.4581636933  0.8035823439  0.3220169424  0.3110491475
##  [706]  1.3244976198  0.3839630856 -0.7476216896  0.3319988357 -0.5387762724
##  [711] -0.1625666114  0.6550736530 -0.0514749417 -0.4933467598 -0.1293641017
##  [716] -0.1953278402 -0.2891512889  0.3414612573  0.3132725921  0.7134871287
##  [721] -0.0668673294  0.5572351776  0.2330029164 -0.9563300138  0.0862545021
##  [726] -0.2450870059  1.2550125555  0.0501053287  1.1463020822  0.1100459393
##  [731]  0.1058483940  0.9372067896 -0.2253007212  1.8737527016  0.3154241216
##  [736]  0.0905923341 -0.5824659406 -0.9522539925  1.3087014181  0.0730406749
##  [741]  0.8248782123  0.7099535044 -0.6270119367  0.1269308144 -0.0169995947
##  [746]  0.8403133842  0.1899102911 -0.2871934384  0.3318943887  0.0874331932
##  [751]  0.6063202326  0.4570625396  0.4085615140 -0.1682895362  0.2945981202
##  [756] -0.7819939606  0.0071532738 -0.2384185792 -0.3648893669 -0.0498454991
##  [761] -0.2142492774  0.1572268099  0.8163492797 -0.0513690243 -0.3897083473
##  [766]  0.3684747002 -0.4704793797 -0.1612048019 -0.4380467393 -0.3966149813
##  [771]  0.0562920399 -0.6722889314 -1.0495452207 -0.0583158300  0.5026000763
##  [776]  0.5271130674 -0.1348922356  1.0840572305 -1.1557593519  1.2394994432
##  [781]  0.1167002363  0.1262245825  0.2603972989 -0.3284226584  0.0821343439
##  [786] -0.3174517693 -0.4279409247 -0.4751407124 -0.0922665489 -0.5233635135
##  [791] -0.1613635425  0.5546617477  0.9977681139 -0.0962218865  0.4171521129
##  [796] -0.7139668838  0.7922355638  1.1875750053 -0.1033945249  0.0839667178
##  [801]  0.3175170421  0.3387519046 -0.5268221546  1.5008197143  0.6466680944
##  [806]  0.1248061109  0.3939080400  0.1458973741  0.0005170502  0.1210850965
##  [811]  1.2031344762 -1.1296769256  0.2160739236  0.5766158031  0.0720477622
##  [816] -0.5438661484 -0.1450483776  0.3119003000  0.4784040660  0.3190008751
##  [821] -0.3529304081  0.0163980401  0.3670471085 -0.8601589176  0.3415964585
##  [826] -0.1149468711  0.5342450271 -0.5379945199  0.5231050763 -0.1492960736
##  [831]  0.3129603491 -0.8929608360  0.5588670329 -0.1265574035  0.3259251126
##  [836] -0.4862892971 -0.3309118535 -0.1942939266  0.4378778510 -0.1453042705
##  [841] -1.4073733163  0.3543977777 -0.4266011988 -0.4730232076 -0.0900791822
##  [846] -0.2077420792  0.0745256373  0.3085573210  0.5150501711  0.0170843791
##  [851] -0.3033125940  0.5398298180  0.8030588462  0.4234020545 -0.2215951896
##  [856] -0.0767608525 -0.4908694958  0.0286084942  0.4259608019  0.2007882683
##  [861]  0.3345788067 -0.1383194716 -0.4554836630  0.0471174656 -0.3093307806
##  [866] -1.6443055839  0.5886678732 -0.1103788003 -0.0032576360 -0.4577890779
##  [871] -0.1831417362 -0.1327665240  0.2918626299 -0.2047959712 -0.6900010522
##  [876] -0.3224128576  0.0071067748  0.3855867772  0.6295661582  0.1027880489
##  [881]  0.3889386542  0.4362646709 -0.6722080926 -0.1758714011  0.3013364477
##  [886] -0.6847792531 -0.0663379758  0.3840752693  0.1068349503  1.0705441859
##  [891]  0.5630153981 -0.0295032495  0.3236177405  0.0903913023 -0.0504859886
##  [896] -0.5460065048  0.1981333640  0.3433538653  0.1965816940 -2.1214991292
##  [901]  0.6443121211  0.5200225937  0.5471289226 -0.8127836352  0.4063305079
##  [906]  0.3488663605  0.2704555002 -0.2447065587 -0.0477307543  0.3262244106
##  [911] -0.0253867353  0.0904699133 -0.2218173600  0.5368440028 -0.7065507283
##  [916]  0.3339218367 -0.6286619415 -0.2803150545  0.3724352496 -0.0497049890
##  [921] -0.0256023933  0.7204596919  0.1059418070  0.2190487011  0.3786743426
##  [926] -0.1162092937  0.6648629596  0.3959501524 -0.4284223266  1.1606564285
##  [931] -0.4584179386  0.1086658056 -0.4961742456  0.1504741803 -0.7555601140
##  [936]  0.4926934359  0.0705844731  0.5327394451  0.7816895581 -0.0987285093
##  [941]  0.9882874942 -0.1067024730  0.0163742869  0.1260301403 -0.4859996710
##  [946] -0.2201891633  0.0293496565 -0.3195782394 -0.0826881205  0.6639765548
##  [951] -0.0422087616 -0.5073877063  0.1767041611 -0.0439681548 -0.0435282434
##  [956]  0.4137505982 -0.2483454642 -0.2801164673 -0.6670424517  0.4511953777
##  [961]  0.8043887417  0.3428589599 -0.4071218044 -0.2505730936  0.3827675197
##  [966]  0.1112283490  0.3587858607 -0.1538929525  0.4780473350  0.1504847159
##  [971]  0.1068914209 -0.8378077207  0.5044783080  0.6380573179  0.0581505998
##  [976]  1.0457429475  0.6843453452 -0.1684576211  0.1106978020  0.5174822690
##  [981] -0.3221179142  0.6952629273  0.1493466589  0.0896054793  0.3732991934
##  [986] -0.5566929609  0.1527909113 -0.2345053456  0.2774883429  0.0108206400
##  [991] -0.4099259556 -0.4640321721  0.5793579737 -0.6794757988 -0.9805617630
##  [996] -0.1084787286 -0.5015732848 -0.2430648192  0.5673809660 -0.0203864128
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
##   0.58802197   0.27873563 
##  (0.08814395) (0.06232359)
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
## [1] -0.1927104 -0.1233473 -0.9248511 -0.3958816  0.1818407 -1.1126661
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
## [1] 0.0342
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9145528
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
## t1*      4.5 0.01211211   0.9022814
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 2 3 4 5 6 8 
## 1 2 1 1 1 1 3
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
## [1] -0.0056
```

```r
se.boot
```

```
## [1] 0.9315774
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

