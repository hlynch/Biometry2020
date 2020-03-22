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
## 0 1 2 3 4 6 8 
## 1 1 1 3 1 2 1
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
## [1] 0.0012
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
## [1] 2.740617
```

```r
UL.boot
```

```
## [1] 6.261783
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
##    [1] 4.7 3.6 4.3 3.1 3.8 4.8 3.6 6.0 3.8 4.3 3.9 4.7 5.0 5.8 5.2 5.7 5.9 4.6
##   [19] 4.2 4.2 5.5 4.2 5.8 4.0 5.5 4.4 4.3 3.8 3.7 4.6 3.7 4.9 5.3 4.7 4.5 4.3
##   [37] 4.7 4.9 3.3 5.0 6.0 5.3 5.7 5.2 5.3 6.2 3.2 2.3 6.0 3.5 5.7 4.7 2.4 6.1
##   [55] 3.8 3.7 3.5 4.3 4.6 4.4 3.6 3.8 6.0 3.3 4.8 4.5 4.1 3.2 5.0 3.4 3.7 4.9
##   [73] 3.0 4.9 6.1 5.3 5.6 3.2 4.7 4.3 3.5 3.9 5.2 5.6 4.5 5.2 5.9 3.6 3.0 6.1
##   [91] 3.7 4.8 3.8 5.5 4.8 2.2 2.8 3.0 4.2 5.3 6.6 3.9 2.7 3.6 5.6 4.5 3.2 4.8
##  [109] 3.8 4.5 3.2 4.0 5.6 3.6 3.0 4.7 4.6 3.8 3.6 3.5 3.7 5.2 4.9 5.8 4.2 4.7
##  [127] 3.9 3.5 3.3 4.7 4.5 3.6 4.8 3.0 4.8 4.3 3.4 5.3 4.9 4.5 5.9 4.2 3.6 4.6
##  [145] 3.6 3.1 5.0 4.4 6.4 5.4 4.2 4.3 5.4 3.7 5.0 1.7 2.2 4.1 5.8 4.8 4.6 4.2
##  [163] 4.5 3.5 5.3 5.2 5.3 4.9 4.0 6.0 6.1 3.6 2.7 3.5 3.2 4.7 4.0 3.5 3.9 2.9
##  [181] 3.2 4.2 3.3 2.8 4.4 5.3 4.7 3.4 5.0 5.8 5.5 5.0 4.7 3.5 3.4 3.5 5.3 2.8
##  [199] 5.7 3.7 4.4 4.6 4.7 4.7 3.3 3.8 3.9 3.3 5.7 2.3 4.8 4.4 2.9 4.7 4.9 5.8
##  [217] 3.6 3.4 6.7 5.2 5.8 2.9 3.8 5.2 4.4 5.1 4.1 5.5 3.8 2.2 4.6 2.9 6.0 7.5
##  [235] 5.1 3.9 5.0 4.2 3.8 6.9 4.6 4.1 3.4 4.9 4.3 1.9 4.6 3.9 5.8 5.8 3.7 5.7
##  [253] 4.9 3.2 3.2 3.5 4.7 3.3 4.4 6.2 3.0 4.5 6.3 5.3 5.6 5.6 5.9 3.8 5.4 5.2
##  [271] 3.5 4.8 3.9 3.5 4.8 4.7 4.3 4.6 4.3 5.7 3.8 4.3 5.3 4.0 4.1 4.3 3.5 4.6
##  [289] 4.4 5.0 3.9 2.0 5.0 4.5 4.6 3.9 4.6 4.9 5.5 4.4 4.1 3.8 4.5 3.9 4.8 2.0
##  [307] 5.3 4.5 5.9 5.2 5.2 3.4 4.4 5.0 3.3 3.4 4.1 5.8 4.2 5.4 2.8 2.8 4.4 4.1
##  [325] 5.5 3.7 5.5 4.2 3.8 4.1 4.6 5.6 4.4 4.6 3.9 3.6 4.4 4.7 5.2 6.0 4.5 5.2
##  [343] 4.5 5.0 5.5 3.1 4.1 4.7 4.2 3.5 4.8 4.5 3.1 4.8 4.0 5.3 4.0 4.7 4.5 4.0
##  [361] 4.1 4.1 5.6 3.0 3.6 3.7 5.0 3.5 4.6 3.9 4.2 4.7 3.1 4.4 4.2 6.1 4.7 3.3
##  [379] 5.6 2.4 4.6 4.0 4.3 5.1 5.4 5.7 5.4 5.7 3.5 3.8 4.3 5.0 5.5 3.6 5.2 4.2
##  [397] 5.2 4.7 4.6 5.2 3.6 3.8 5.2 2.7 6.0 5.2 3.9 4.4 5.8 5.4 6.4 3.6 5.6 4.2
##  [415] 3.6 5.2 3.7 4.7 4.5 5.6 4.2 5.7 4.0 4.1 3.9 4.4 5.8 5.2 2.8 5.7 3.6 3.3
##  [433] 6.8 3.5 4.5 5.4 3.5 3.3 4.9 2.9 3.7 4.1 4.3 5.0 4.5 4.4 5.0 6.2 3.9 4.5
##  [451] 3.3 3.9 3.9 5.9 4.9 3.9 4.3 3.9 4.3 5.4 5.3 4.1 4.4 5.0 4.0 5.6 4.1 3.5
##  [469] 4.4 5.3 3.2 5.0 5.0 5.0 4.6 4.0 5.6 3.1 3.7 2.7 5.3 5.0 4.9 3.5 4.5 3.9
##  [487] 3.7 3.0 5.9 5.5 3.6 3.2 3.6 4.4 4.7 3.0 4.2 4.6 4.7 5.3 5.9 4.5 3.4 4.0
##  [505] 4.7 7.1 3.1 4.9 3.9 5.2 4.4 4.1 5.2 3.6 4.2 2.8 4.5 5.9 4.5 3.7 4.3 5.5
##  [523] 6.6 2.9 4.9 3.6 3.1 4.9 5.1 5.4 4.2 5.7 5.2 4.5 3.0 5.3 3.6 5.2 4.0 4.1
##  [541] 4.3 4.8 4.0 5.0 5.1 5.7 4.4 4.7 4.5 6.9 4.1 3.4 5.6 4.2 4.1 4.4 5.5 4.5
##  [559] 6.7 4.9 5.5 3.5 5.0 3.8 6.6 4.4 3.5 5.4 3.4 3.3 4.8 2.6 4.6 4.0 3.6 5.6
##  [577] 4.8 3.0 3.7 5.8 4.7 4.7 4.9 5.7 3.6 4.6 4.3 5.2 4.3 4.5 3.5 4.9 3.6 3.7
##  [595] 4.8 4.1 5.4 5.3 4.0 4.5 5.5 4.4 5.8 4.5 5.6 5.6 3.1 3.2 3.9 6.0 4.3 3.2
##  [613] 5.1 6.1 3.5 4.2 3.4 3.5 5.3 4.0 4.8 4.0 3.5 4.5 4.5 6.2 5.9 3.9 3.5 5.1
##  [631] 3.8 4.9 4.5 6.3 2.3 3.6 4.4 4.8 5.3 4.5 3.9 4.1 5.1 3.7 3.6 5.7 4.5 4.8
##  [649] 5.4 6.2 3.1 3.8 3.4 4.4 4.2 3.4 4.6 2.6 4.7 4.2 3.6 4.7 5.3 3.9 6.2 5.2
##  [667] 3.9 4.7 4.9 5.6 4.3 4.2 4.1 4.3 3.1 2.5 4.3 4.9 5.9 3.5 5.3 6.0 4.3 5.1
##  [685] 5.2 5.0 6.4 4.6 5.3 4.2 4.4 5.5 4.9 4.4 4.2 5.5 3.6 4.8 4.0 6.3 6.2 5.9
##  [703] 4.7 3.4 4.6 4.8 7.0 3.6 4.9 4.3 4.4 4.2 4.6 5.0 3.1 3.1 5.8 4.9 4.3 5.5
##  [721] 3.2 4.6 3.7 4.2 4.7 3.0 5.6 5.5 5.6 4.4 3.3 5.0 2.9 3.0 5.8 5.0 6.4 5.3
##  [739] 3.2 5.0 5.0 4.1 5.4 2.9 4.8 4.1 5.4 4.8 4.9 5.5 6.1 4.1 5.4 5.1 5.7 3.9
##  [757] 5.2 4.1 4.0 5.1 5.1 4.0 4.1 5.2 4.2 6.6 4.8 6.2 4.7 4.2 5.5 4.9 4.4 4.4
##  [775] 5.1 4.9 1.9 4.4 4.1 3.3 3.8 6.0 5.0 4.8 5.1 4.5 5.2 4.9 3.7 5.0 4.3 4.3
##  [793] 5.2 3.9 4.2 4.4 4.0 4.6 4.2 4.7 5.4 4.5 3.8 5.4 3.2 4.1 5.6 5.9 3.8 4.0
##  [811] 5.6 4.5 2.8 4.3 1.9 4.6 4.8 2.8 4.7 5.6 4.3 4.7 5.2 5.1 4.0 5.9 4.2 4.6
##  [829] 7.2 4.7 4.4 2.6 4.0 5.8 3.9 5.3 4.4 4.8 4.5 3.3 3.8 4.5 5.5 4.7 3.4 4.2
##  [847] 5.2 3.6 3.7 3.9 5.9 5.3 5.0 4.7 4.3 3.8 4.2 4.9 4.2 4.6 3.7 3.5 4.3 5.3
##  [865] 4.7 4.6 5.5 4.1 3.8 4.5 5.3 4.4 4.3 4.7 4.9 4.7 5.7 3.3 3.7 5.2 2.7 5.2
##  [883] 5.8 4.3 3.7 5.6 5.1 3.4 3.9 5.9 4.8 5.6 3.6 4.4 4.4 4.4 3.5 3.9 5.0 5.1
##  [901] 4.2 4.7 4.2 3.3 4.8 5.7 6.2 5.5 4.1 4.4 3.1 3.7 3.2 5.7 4.0 4.6 4.5 2.5
##  [919] 5.0 4.0 4.3 3.8 5.2 5.0 4.7 4.8 4.4 6.0 2.7 4.1 5.5 4.0 4.1 6.0 2.3 3.9
##  [937] 5.3 5.6 5.0 5.3 3.2 4.0 3.2 5.2 4.0 3.6 3.1 3.5 4.8 4.4 4.6 5.5 4.6 3.7
##  [955] 4.3 5.8 4.1 4.4 5.3 4.2 3.3 3.8 4.7 6.5 3.5 3.7 3.7 5.5 4.8 3.7 4.9 5.1
##  [973] 5.7 5.3 4.0 4.1 4.8 4.2 4.6 4.9 4.2 5.0 5.5 5.3 3.4 3.6 5.1 4.9 3.8 3.6
##  [991] 2.6 4.0 4.9 5.0 4.3 4.7 3.9 5.0 3.6 5.3
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
##    [1] 5.0 4.9 3.1 5.7 5.0 3.3 3.8 4.3 4.1 3.3 5.2 4.2 4.2 4.6 3.3 5.9 5.5 3.2
##   [19] 3.2 4.9 3.5 6.5 4.1 6.0 3.8 5.0 5.1 4.3 4.4 4.5 5.4 4.2 5.5 4.3 6.0 3.9
##   [37] 4.4 3.4 2.7 5.6 2.8 4.8 4.1 5.1 3.5 4.1 5.5 4.1 4.0 4.4 4.0 4.8 3.5 4.9
##   [55] 3.8 5.1 6.1 4.3 3.1 3.2 5.4 4.8 4.6 5.4 4.1 3.5 5.0 4.4 4.6 5.0 4.7 4.3
##   [73] 3.5 4.5 6.1 4.8 5.0 3.8 4.6 3.5 4.1 4.5 4.6 4.1 4.4 4.2 3.7 4.1 4.8 3.8
##   [91] 5.3 5.8 6.0 5.0 3.8 4.4 5.7 4.6 5.6 4.9 5.4 4.8 3.8 3.5 5.1 4.6 4.1 4.9
##  [109] 3.9 5.5 5.2 3.7 5.6 3.5 3.2 4.7 5.4 4.9 5.4 4.5 3.4 5.0 6.2 3.0 5.2 4.5
##  [127] 4.9 5.7 4.7 5.4 5.0 4.7 3.2 3.1 4.6 4.6 5.9 3.4 5.8 5.3 4.7 4.4 5.1 4.7
##  [145] 4.2 4.1 4.4 4.5 5.1 4.8 5.4 5.4 3.9 5.4 2.4 3.4 4.3 4.4 4.0 3.8 4.6 5.7
##  [163] 3.7 4.7 4.4 4.2 4.4 4.8 4.7 5.8 5.2 4.9 4.6 5.0 5.2 4.8 3.3 3.6 4.5 3.6
##  [181] 3.7 2.0 5.0 4.1 3.2 4.2 3.0 5.1 4.4 5.0 3.5 5.9 6.0 5.0 2.7 3.8 4.8 3.8
##  [199] 4.6 2.6 4.3 4.8 4.7 3.4 4.5 2.8 2.7 4.5 3.8 3.5 5.5 4.9 5.0 4.7 3.5 3.7
##  [217] 5.1 4.7 3.4 4.4 4.3 4.6 4.9 3.5 4.2 5.8 3.8 5.9 4.4 3.5 4.1 3.7 3.5 4.9
##  [235] 3.7 4.9 4.6 4.9 6.2 4.5 4.3 5.7 3.6 4.6 3.4 4.5 4.6 3.8 3.1 3.5 5.0 3.7
##  [253] 5.6 5.3 6.2 3.8 3.1 4.4 2.9 4.5 4.0 5.3 4.4 4.1 5.6 3.4 3.5 2.3 4.8 3.8
##  [271] 4.7 5.4 6.0 6.0 5.2 3.8 2.9 3.6 4.0 4.5 4.9 3.3 5.4 4.1 4.3 4.0 5.1 5.2
##  [289] 3.8 5.2 5.5 4.2 5.3 3.8 3.5 2.9 3.3 4.2 3.3 5.0 3.7 3.9 5.3 4.3 2.9 3.0
##  [307] 4.6 5.1 4.1 4.8 2.8 5.2 3.0 4.9 4.5 4.1 4.0 5.3 5.0 4.5 6.6 4.9 3.8 6.3
##  [325] 4.1 5.7 5.0 5.4 4.6 4.5 3.5 5.1 4.3 5.0 3.9 5.8 4.5 3.5 4.8 3.7 3.0 5.4
##  [343] 4.4 6.7 4.0 6.1 4.2 4.9 4.4 4.0 4.6 4.9 5.2 4.2 4.2 4.6 5.5 3.4 4.5 5.7
##  [361] 3.9 6.2 3.1 5.9 5.7 5.2 3.7 5.1 5.1 5.7 4.1 5.0 5.6 6.1 4.5 4.7 4.0 4.0
##  [379] 4.2 5.2 6.3 3.8 4.9 4.7 5.5 3.5 3.7 5.1 3.4 3.0 5.0 4.5 3.9 4.3 3.7 5.8
##  [397] 4.0 5.0 4.2 3.8 3.8 3.3 4.7 4.9 6.0 5.0 6.9 3.8 4.1 4.7 4.2 3.7 4.0 3.5
##  [415] 4.1 3.9 6.0 4.9 4.9 2.7 5.9 5.8 4.1 5.4 4.1 4.2 4.3 5.2 5.9 3.1 4.5 5.5
##  [433] 3.1 4.1 3.4 3.5 2.7 3.0 3.7 5.4 6.1 5.1 3.8 3.7 3.5 5.7 5.4 4.7 5.3 5.1
##  [451] 2.7 5.1 4.6 4.1 5.8 5.7 3.4 4.9 6.4 3.9 4.5 4.7 4.0 2.5 4.6 4.4 6.2 5.4
##  [469] 5.1 4.1 3.2 4.5 4.6 5.5 4.9 3.8 4.9 4.3 4.6 5.0 4.4 5.3 5.6 3.7 5.9 4.2
##  [487] 3.9 5.2 4.1 4.7 2.6 4.1 4.9 3.8 3.7 2.8 6.5 4.8 3.8 5.2 4.9 5.3 4.3 4.4
##  [505] 6.7 4.6 4.9 4.4 3.3 4.6 5.8 4.2 5.0 4.5 3.9 5.7 4.3 4.3 5.0 5.0 3.8 4.1
##  [523] 3.4 5.6 4.2 2.6 4.3 2.2 4.5 3.8 5.8 3.6 3.8 5.1 4.2 3.2 4.4 3.9 3.0 4.8
##  [541] 6.4 4.2 4.9 4.5 5.9 4.9 5.8 4.3 4.9 4.6 5.9 5.5 5.3 4.1 4.5 4.1 3.8 4.4
##  [559] 3.8 4.1 5.5 4.5 4.3 3.3 4.1 4.9 3.9 4.7 4.1 5.6 3.8 3.8 5.0 2.9 6.0 5.8
##  [577] 5.1 5.4 5.9 3.5 3.2 3.5 3.9 4.0 3.8 6.2 3.9 3.4 5.6 4.0 5.0 3.7 4.6 2.8
##  [595] 6.2 5.5 5.1 5.9 4.6 4.6 5.9 3.4 3.8 3.4 5.0 4.8 3.4 3.5 5.6 4.7 3.7 2.6
##  [613] 5.8 5.1 5.4 4.6 3.3 4.2 4.0 3.7 5.8 5.0 4.1 3.8 5.6 3.9 3.9 3.3 5.7 5.6
##  [631] 3.1 3.7 6.7 4.8 5.1 5.0 3.1 3.2 3.9 4.8 4.8 6.5 5.0 2.0 4.0 3.9 3.6 3.9
##  [649] 5.2 3.1 3.0 3.3 4.4 3.7 4.6 3.4 3.3 3.7 4.3 4.7 4.3 5.5 4.2 5.8 5.6 5.6
##  [667] 5.8 4.7 6.0 4.7 4.1 4.3 4.5 4.7 5.2 5.0 5.1 4.0 6.2 5.0 4.1 5.5 5.2 4.6
##  [685] 5.4 5.7 3.5 4.3 4.3 2.8 4.0 5.9 6.0 4.7 3.3 5.0 5.9 4.6 6.0 3.8 2.9 5.0
##  [703] 4.1 3.9 4.9 4.7 4.4 3.9 3.4 3.4 4.0 3.0 5.6 2.5 4.5 4.8 7.7 2.7 3.6 5.1
##  [721] 6.7 5.6 4.6 4.7 5.6 3.5 2.5 6.0 3.8 4.4 3.9 5.0 5.9 4.7 3.1 3.9 3.5 6.6
##  [739] 4.0 4.5 6.5 4.8 3.7 3.0 3.4 4.5 5.0 3.9 6.6 4.1 4.3 4.1 4.5 6.3 4.7 5.7
##  [757] 5.9 3.0 4.9 3.4 3.8 4.7 3.6 5.6 4.3 3.8 4.6 4.2 3.4 6.3 5.7 3.0 4.4 5.4
##  [775] 4.6 3.7 6.0 5.7 3.3 6.0 4.2 4.3 4.9 4.1 3.7 4.0 2.6 5.8 5.1 4.9 5.8 4.5
##  [793] 3.8 3.8 5.5 2.9 5.5 4.0 4.9 5.8 4.7 4.4 5.3 2.7 3.9 3.5 4.4 3.9 3.8 4.7
##  [811] 6.1 5.9 5.3 4.3 4.7 5.5 5.2 3.4 4.4 3.3 4.6 3.0 5.3 4.8 3.2 4.7 3.5 4.4
##  [829] 5.6 4.1 4.0 3.8 4.2 6.1 3.7 4.9 5.0 5.9 4.6 3.9 4.8 4.6 4.9 4.3 5.2 4.3
##  [847] 4.3 5.9 5.1 4.2 5.6 5.4 3.9 5.5 3.7 5.2 2.1 5.3 4.3 4.2 3.8 5.0 4.8 4.2
##  [865] 5.6 3.4 4.1 2.8 3.8 5.6 4.6 3.4 4.7 4.4 4.3 3.7 5.3 5.5 5.6 4.7 4.5 4.1
##  [883] 3.3 5.7 4.5 4.0 6.0 4.1 3.8 3.7 3.4 4.5 4.8 3.7 5.3 3.6 4.2 4.7 3.7 2.9
##  [901] 4.6 3.9 4.9 4.1 3.6 5.0 4.9 5.0 6.2 5.3 4.7 5.0 3.7 4.9 4.6 4.4 3.4 6.4
##  [919] 4.6 5.3 4.5 4.3 4.3 5.0 6.4 3.6 4.5 3.5 5.4 4.7 5.2 4.8 3.5 4.1 4.6 3.2
##  [937] 3.7 4.2 4.8 4.3 4.7 5.3 4.4 4.8 3.7 5.7 4.9 4.9 4.8 4.3 3.5 3.6 5.4 3.8
##  [955] 3.7 3.7 5.9 4.3 3.3 5.3 4.7 2.0 4.6 5.2 3.5 4.4 5.7 3.6 3.8 4.8 3.0 3.5
##  [973] 3.5 2.8 6.1 5.0 3.9 5.5 3.5 4.5 4.4 3.3 4.4 3.9 3.9 4.6 5.5 3.9 4.9 3.7
##  [991] 4.5 4.2 4.5 5.0 4.2 5.6 4.9 4.8 4.5 5.1
## 
## $func.thetastar
## [1] -0.0126
## 
## $jack.boot.val
##  [1]  0.49433428  0.38997361  0.31460674  0.17341390  0.09501466 -0.15993976
##  [7] -0.27560241 -0.24858757 -0.41424731 -0.49048913
## 
## $jack.boot.se
## [1] 0.9988514
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
##    [1] 5.3 5.6 4.5 4.7 5.2 4.5 4.4 4.2 4.4 5.4 5.5 6.6 3.6 4.8 4.6 3.5 5.3 3.8
##   [19] 5.4 4.8 4.8 4.7 4.7 5.5 4.6 3.5 4.1 4.4 5.4 5.3 3.8 5.6 3.6 3.0 5.0 4.2
##   [37] 4.2 2.5 4.0 3.8 3.7 5.1 4.2 3.5 3.8 5.5 5.0 3.7 4.1 4.9 3.9 6.2 4.5 4.6
##   [55] 4.8 5.4 3.6 4.9 5.4 3.8 4.7 4.5 5.1 4.6 5.0 4.1 4.3 4.6 2.7 5.1 3.4 3.6
##   [73] 4.6 3.9 5.7 4.6 4.9 5.2 3.3 3.2 2.9 5.5 3.4 4.3 4.7 2.6 3.7 4.1 4.6 4.2
##   [91] 6.1 4.2 5.7 4.2 2.6 5.1 2.9 4.7 2.9 5.3 4.9 4.1 5.5 2.9 4.2 5.1 5.7 4.5
##  [109] 5.6 4.7 3.6 5.4 4.8 3.4 4.8 4.9 3.7 4.7 5.1 4.6 6.6 3.7 2.5 3.2 6.5 5.8
##  [127] 4.1 4.5 4.1 4.7 4.5 4.9 4.5 4.8 5.6 3.7 5.4 5.7 3.5 3.7 3.9 3.3 4.5 4.0
##  [145] 3.8 4.0 4.3 4.0 5.0 5.2 4.7 3.9 5.0 3.3 5.1 3.5 5.3 4.7 5.3 4.3 2.7 5.5
##  [163] 4.2 3.0 3.5 5.1 5.1 3.8 3.8 5.6 4.1 5.4 3.8 4.2 4.8 5.9 3.2 5.6 4.8 4.2
##  [181] 5.6 4.7 4.5 3.8 3.0 5.0 4.8 5.0 5.5 4.9 2.8 2.1 4.7 3.5 3.7 5.2 4.4 5.7
##  [199] 3.3 4.8 3.9 4.5 4.8 3.3 5.1 4.0 5.1 4.4 5.3 5.6 4.1 3.3 5.2 3.5 4.4 5.6
##  [217] 4.1 5.1 5.4 4.0 5.2 4.4 4.3 4.7 5.3 4.8 5.2 6.4 3.5 3.0 4.2 5.9 4.6 5.0
##  [235] 4.3 3.8 4.6 6.0 5.1 4.1 3.7 5.3 5.4 3.7 5.2 6.0 2.8 4.6 5.9 3.7 5.5 6.1
##  [253] 5.0 4.6 4.5 4.6 4.8 4.4 4.5 4.2 2.6 4.2 4.8 4.7 4.1 5.1 5.4 4.0 4.9 3.7
##  [271] 4.3 3.6 4.9 4.8 5.5 6.1 4.4 5.9 5.2 4.2 4.7 5.5 4.9 4.4 5.2 4.3 2.8 4.3
##  [289] 6.8 4.8 5.9 4.9 5.8 5.1 3.1 4.8 3.1 4.4 4.8 4.8 3.9 3.4 4.8 6.5 4.7 4.3
##  [307] 5.7 3.9 5.9 5.1 5.4 4.9 5.1 6.7 4.5 2.4 5.7 3.4 5.4 4.8 3.5 4.6 5.1 5.4
##  [325] 5.2 3.7 5.4 5.5 4.9 5.3 4.8 4.0 5.1 4.8 4.4 4.1 5.2 3.7 4.7 2.9 4.5 3.3
##  [343] 4.2 5.7 4.2 5.9 4.9 5.2 4.6 3.4 6.2 4.6 4.6 4.5 4.2 4.3 3.9 4.0 3.5 4.9
##  [361] 3.5 3.8 4.6 5.0 5.3 2.7 4.6 5.1 3.1 5.0 5.3 4.1 5.4 4.6 4.5 4.1 3.0 5.2
##  [379] 3.2 4.1 5.2 3.1 4.2 4.0 3.5 5.0 5.4 3.7 5.4 4.1 2.8 4.6 5.3 4.4 4.1 5.1
##  [397] 4.0 4.4 4.8 4.3 4.6 4.5 3.5 5.2 5.7 4.3 5.8 3.3 6.5 3.9 4.8 4.4 5.6 4.2
##  [415] 4.9 2.9 5.4 4.7 4.7 4.8 3.8 3.4 5.7 4.1 4.0 5.7 6.0 4.2 6.0 4.5 2.5 3.7
##  [433] 5.2 3.2 4.8 4.9 3.3 5.8 3.9 3.6 5.0 4.2 5.0 3.3 5.6 4.0 4.5 6.1 5.4 4.7
##  [451] 5.1 3.7 6.6 3.4 3.5 5.0 4.9 4.8 5.0 4.6 5.0 4.8 4.0 2.9 4.4 3.6 5.0 4.1
##  [469] 5.2 4.9 4.2 4.5 5.2 5.9 6.0 2.7 4.1 3.5 5.7 4.6 3.8 4.8 4.8 5.3 4.6 6.0
##  [487] 5.2 4.4 4.4 5.6 4.6 2.3 4.8 4.0 5.2 4.9 4.3 5.6 5.6 4.5 4.6 4.0 4.4 3.3
##  [505] 4.7 4.7 2.1 6.0 5.3 5.1 4.5 4.0 3.4 2.9 4.5 3.7 4.6 4.0 4.4 4.6 6.1 4.5
##  [523] 3.5 4.4 3.2 4.4 3.7 7.2 3.4 4.2 5.2 4.0 4.3 3.6 2.7 5.3 5.5 6.1 5.1 3.4
##  [541] 3.9 5.0 4.1 4.2 6.3 4.7 6.4 5.7 5.0 2.7 3.6 2.9 5.8 4.6 7.2 4.5 4.0 4.7
##  [559] 4.6 3.7 4.0 3.3 2.9 5.0 3.5 3.8 4.3 5.3 3.1 3.7 4.5 4.6 4.7 5.2 3.6 4.5
##  [577] 3.9 3.8 5.3 4.4 3.2 4.4 4.4 6.1 4.3 2.8 4.4 5.9 3.8 4.7 5.1 3.6 3.5 5.1
##  [595] 3.9 5.4 4.5 2.9 6.0 3.0 3.5 6.4 3.2 4.1 3.7 4.9 3.8 3.6 6.5 4.8 4.6 3.6
##  [613] 2.6 3.9 6.5 5.7 3.2 5.1 4.2 4.0 3.6 5.3 5.4 2.6 4.1 4.4 4.0 4.0 4.4 5.2
##  [631] 4.2 4.4 3.9 4.8 3.1 3.7 5.3 4.2 2.3 4.3 3.5 4.5 5.9 3.9 4.9 4.7 4.9 4.7
##  [649] 5.4 4.5 5.0 4.7 5.3 3.8 6.0 5.4 3.7 6.3 3.3 4.4 3.3 4.6 5.5 4.6 6.4 3.1
##  [667] 3.9 5.1 5.0 4.5 4.8 4.0 3.3 5.3 4.8 5.1 4.6 3.8 4.9 6.8 4.8 5.2 4.3 4.8
##  [685] 4.9 5.0 3.5 5.2 5.4 3.9 4.2 4.4 5.2 4.2 3.9 4.0 3.1 6.3 3.7 4.0 4.7 3.5
##  [703] 4.5 4.0 3.4 4.0 5.0 2.7 6.3 4.2 5.1 4.7 3.6 5.1 5.4 2.1 6.0 6.0 4.7 3.8
##  [721] 5.7 3.4 4.0 5.1 4.5 4.4 5.7 2.8 3.5 6.3 2.7 4.0 4.1 4.0 4.9 5.6 5.0 4.3
##  [739] 4.5 5.0 3.1 4.6 4.6 4.3 5.4 3.8 4.3 4.0 5.5 3.3 4.9 3.8 4.8 3.8 5.2 5.2
##  [757] 3.8 5.3 4.1 4.2 3.8 5.7 3.7 4.0 5.3 4.8 4.1 3.7 3.7 5.2 5.5 4.8 5.2 5.5
##  [775] 5.1 3.6 3.4 2.1 3.2 4.3 4.6 5.5 5.9 3.7 5.0 4.1 5.6 5.3 4.4 4.7 4.0 4.2
##  [793] 3.0 4.3 4.5 4.1 4.9 3.3 5.4 4.1 6.3 5.7 4.6 5.2 5.8 4.4 3.6 4.6 6.8 6.0
##  [811] 3.4 3.4 4.3 3.3 5.4 2.9 4.7 4.7 4.3 5.8 5.4 3.5 6.6 5.1 4.6 4.2 3.9 4.9
##  [829] 3.5 5.2 5.7 4.5 6.1 5.4 3.6 5.4 3.6 5.2 5.1 5.1 3.7 3.4 2.6 6.2 4.7 2.6
##  [847] 4.0 4.3 4.3 4.6 4.8 5.7 4.8 3.9 5.9 3.6 2.6 3.2 4.8 5.1 4.6 4.2 3.4 3.2
##  [865] 3.7 3.0 4.6 7.2 3.4 5.5 5.1 4.1 3.4 5.4 4.4 3.9 4.7 3.8 4.8 2.4 5.9 3.6
##  [883] 4.5 3.9 3.9 3.6 4.2 3.2 5.4 5.7 4.6 4.0 5.0 3.9 3.7 5.7 5.0 5.3 3.7 3.8
##  [901] 4.4 5.2 4.7 2.6 4.3 5.2 5.1 3.8 3.4 4.4 4.4 5.1 4.6 4.1 5.2 4.2 3.2 4.6
##  [919] 4.3 4.0 2.6 5.3 4.5 5.1 3.0 4.7 4.8 5.6 5.0 5.8 3.0 4.7 2.3 4.3 5.0 5.4
##  [937] 5.0 4.4 5.0 5.0 5.1 3.5 4.6 3.7 4.6 4.0 6.0 5.4 3.8 5.9 3.8 6.0 4.6 4.4
##  [955] 4.2 5.3 5.5 4.5 2.9 3.3 5.9 4.7 4.7 5.1 3.9 4.9 4.0 3.5 6.4 4.6 3.4 4.5
##  [973] 5.4 3.6 4.9 3.7 5.0 5.1 4.3 4.7 2.6 4.6 3.9 6.3 6.5 4.3 4.8 3.1 4.7 4.2
##  [991] 4.6 3.5 4.2 4.6 3.2 4.7 3.0 4.9 3.4 2.9
## 
## $func.thetastar
## 72% 
## 5.1 
## 
## $jack.boot.val
##  [1] 5.500 5.400 5.300 5.100 5.100 5.000 4.996 4.700 4.600 4.400
## 
## $jack.boot.se
## [1] 1.00806
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
## [1] 0.1316489
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
##      shape       rate   
##    9.439873   14.787833 
##  ( 4.149179) ( 6.675632)
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
## [1]  0.2777550  0.4063802 -0.3019270  0.7076923  0.7483140  0.3183859
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
##    [1] -0.4629874755 -0.4766280249  0.7457271318  0.0552261671  0.0003637105
##    [6]  0.1684171115  0.9888172695 -0.7947394489  0.7841257942  0.3882416972
##   [11] -0.3956992058 -1.5104634031 -0.8146788938 -0.4481179670 -0.4423752556
##   [16] -0.4703114588 -0.2137928562  0.5820285445 -0.3174055776  0.3149303616
##   [21]  0.2582248581  0.1083435114  0.2044176945 -0.0744401694  0.1828466903
##   [26]  0.2002862616  0.2802869416 -0.7375160836 -0.1036221179  1.0222628611
##   [31]  0.1514994227 -0.4483515944  0.5253460531 -0.0888237640 -0.0180049906
##   [36]  0.0290349831 -0.3915623425 -1.3096191064  0.6646264932  0.3676843114
##   [41]  0.0074282204 -0.0659164845 -0.3235862540 -0.3431053353 -0.6511486805
##   [46]  0.0333629301  0.4163001564 -0.1889787081 -0.0987917900  1.0250768949
##   [51]  0.1234395165 -1.6868453477  1.0097689193 -1.6683042518  0.7401792338
##   [56] -0.0101196074  0.0009177295  0.1056896386 -1.0408979356 -0.1410777898
##   [61] -0.3721228576 -1.7415797193  0.3494686059 -0.5036330335  0.2566453206
##   [66] -0.1787867229  0.0023463885  0.1170131672 -1.2114278710  0.2294068899
##   [71] -0.6499061676 -0.7672466245 -0.9872711662  0.6901516092  0.2591528657
##   [76]  0.1059132401 -0.5847136257  0.0935873010  0.2971881090 -1.3944986536
##   [81] -2.1557319797  0.5095780335  0.2436440328 -0.4008056145 -0.0223465634
##   [86] -1.3031683749 -0.2850890714 -1.4598074285  0.0458396478  0.6543291697
##   [91]  0.0844600098  0.2087158917  0.2037774103 -0.5561610011 -1.3475755509
##   [96] -0.0863591004 -0.7072403034 -0.8886009204  0.6082469532  0.1976961630
##  [101]  0.2592934335  0.1222744016  0.9257471899 -0.7999127151  0.2355013971
##  [106] -0.4275275662  0.2293277860  0.8386758620  0.3057005116  0.0802677539
##  [111]  0.4656609197 -0.1046777671 -0.0009231949  0.3032088827 -0.9439266702
##  [116]  0.1144080483 -1.2488543573  0.7135862174  0.5216933171  1.3184210175
##  [121] -2.0819879644  0.3431890481 -1.4249036974  0.2733154322  0.7179430997
##  [126] -0.4631223986  0.1859878997 -0.8392442693 -0.1426658424  0.3976552838
##  [131]  0.0643456886 -0.7569111640 -0.1124033044  2.2198463236  0.0101682721
##  [136] -2.1534352237  1.2778635511 -0.1808828273  0.4070152983 -0.1375213711
##  [141] -0.3819258603 -0.1114445553  0.0772650100 -0.0289438800  0.2344580733
##  [146]  0.0906031164 -0.0525620726  0.2648177566  0.5240298929  0.3170355785
##  [151]  0.4329136086  0.5917278314  0.0952367500 -0.0931162079  1.4470937582
##  [156]  0.2893075026 -0.0709871015 -0.0215775810 -1.0001604439 -0.1625128358
##  [161]  0.5717285453 -0.4655041224  0.2557080777 -1.5270383055  1.2082889259
##  [166]  0.0695697655  0.7078614880 -0.4767048495  0.0836679817  0.1493599699
##  [171]  0.0164341504  0.2126608514 -1.2374205297 -2.0698139314  0.8103836427
##  [176] -0.2104735285 -0.1099482163  0.6549102508  0.3631744680 -0.4267039235
##  [181]  0.3045849773 -0.9439266702 -0.8301170488  0.3254377189 -0.8460613749
##  [186]  0.1743310846  0.7192050597  0.9879476746  0.2799976483 -2.0786956399
##  [191]  2.0170109028  0.5985986683 -0.1867260844  0.5606989197 -0.7417935958
##  [196]  0.7021601977 -0.1736922012 -0.1108576526 -1.5411644640  0.0415062037
##  [201]  0.0629640745 -1.3911022043  0.3940626217  0.0600542651  0.3930032243
##  [206]  0.3401222707 -0.1515833559 -2.1007727838  0.1532636176  0.1457325117
##  [211]  0.2271312047  0.1981115764 -0.8130311094  0.4130376268 -1.5524071669
##  [216] -0.0162981328 -0.8420465370 -0.1049593110 -1.3828086200 -0.8287488372
##  [221]  0.0740294919 -1.3983270104  0.0274705541 -0.8285527278  0.6352956534
##  [226] -0.0693182719  0.5888361055 -0.1012888596  0.2048374215 -0.2269177769
##  [231]  0.5894502969  0.8680569609 -0.1500419074 -2.0981072994  0.5995900507
##  [236]  0.4551124345  1.8298079881 -0.1254078821 -0.0112888486  0.4877002570
##  [241]  0.2799138145 -0.1004311356  0.8971307134  0.1229855912  0.2002919556
##  [246]  0.1916374806 -0.8657158251 -0.1210847229  0.5543400148 -0.4654669508
##  [251]  0.2808868696  0.6925562496 -0.9880510033  0.7070493534  0.0819033326
##  [256] -1.1767224332  0.2194874552  1.0643680856  0.2950596312 -0.1384415236
##  [261]  0.3044496857  0.7451651531  0.2631424390  1.1125118959  0.5203208317
##  [266]  0.2439708762 -1.3497767078  1.0622021500 -1.4201932236  0.1281386960
##  [271] -1.0947954626  0.2017872357  0.0355454948 -0.4041827839 -1.4904945333
##  [276] -0.7759514156  0.3990973837  1.0821146933 -0.0642870551  0.2795603532
##  [281]  0.6126368977  0.5459954334 -0.4514005622  0.0572739222 -0.3473236062
##  [286] -0.3524165251 -0.8175696669  0.0911561564  1.2764159739 -2.1624323242
##  [291] -0.3463393122 -0.1973524709  0.0292208142  0.2437599516  0.2700557748
##  [296]  0.7525606161  0.8460637767 -0.2791801190  0.0741073191  0.4308724966
##  [301]  0.2583005163 -0.7739004789 -2.1491865628 -0.4559796863 -0.3417951942
##  [306]  0.0771720861  0.1415575277 -0.0013667085  0.4709603124 -1.2475941985
##  [311]  0.6179496460 -1.5038677356 -0.3961462763  0.1755053870 -1.5937191549
##  [316] -0.9013300607  0.0766818156 -1.9194519788  0.5857442011 -0.8313458340
##  [321]  0.4033873454 -0.9417483819  0.2290295402 -1.3442444368 -0.0828571564
##  [326]  0.0565749287  0.4473238243 -0.1086879200 -0.1941873259 -1.2025820121
##  [331] -1.3019223673  0.4985631068 -0.1464422309  0.3001273594 -2.0819879644
##  [336]  0.1948989955  0.0218289468 -0.1643842098  0.5771124949 -0.1105978222
##  [341] -0.2201926786  0.1985414700 -2.0520768590  0.7060012456 -0.0710051076
##  [346] -0.1499004733  0.5110101624  0.7738211359 -0.1321147484 -0.0846332815
##  [351] -0.7506976959 -0.3806765388 -0.8742046589  0.1049032765 -0.7446877918
##  [356] -0.6403600272  0.1941491182  0.1450582000 -0.8849470963 -1.3132056702
##  [361]  0.4197114517 -0.1532104548  1.2132805650 -0.1522750651 -0.3383668154
##  [366]  0.4676902099  0.5771124949  0.4511825921 -1.0385987186  0.1223580199
##  [371]  0.3952047048  0.2003861782  0.2146262301 -0.3512624527 -0.0515290157
##  [376]  0.4014785488  0.1931695136  0.6578178575 -0.7723346222  0.4226066678
##  [381] -0.2545069553 -0.8841264313 -1.4504341724 -0.0679957964  0.3442930028
##  [386] -0.2023872722 -0.0205676017 -1.4642216156  0.9372721775  0.0180018075
##  [391]  0.1841550569  0.3610147040 -0.5668818959  0.4044835345 -0.7033415934
##  [396]  1.1245069007  0.6643505322  0.0861806792  0.1971110468  0.5711518817
##  [401] -1.1980056663  0.1392384301  0.0882936489 -0.4359918339 -0.8008662452
##  [406]  0.2572411058 -1.5407333561 -0.0282501533 -0.7738438101 -0.8491888976
##  [411]  0.6209385067  0.5457685705 -0.1227357298  0.3916316608 -0.7640304965
##  [416]  1.0768919947  0.4034532170  1.3254499889  0.4982092379 -0.0647388977
##  [421]  0.0175523881 -0.2912519282  0.1192423446 -0.1241481425  1.2457028615
##  [426]  0.3265133854 -0.3812679029 -0.4974772056  0.5049614315  0.0806767022
##  [431] -0.4428283637  0.0444212101  1.2265001320  0.3302733155 -0.7981739589
##  [436] -1.5168265265 -0.1699827574  2.2063042989  0.2716368232  0.2246878590
##  [441]  0.3939174747  0.4021585340 -0.6646424448 -0.1596133387 -0.1186815982
##  [446]  0.0238722301 -0.1097368165 -0.8805204887  0.0389048439  0.0614038523
##  [451] -1.5405260099 -0.5821576551 -0.8530101119  0.1752136285 -0.0285961769
##  [456]  0.4295188257 -0.4681563180  0.1351608793  0.0248455691  0.4485674008
##  [461]  0.3339498728  0.3095976152  0.2237195717  0.4561708714  0.6144183041
##  [466] -1.0038051822 -0.7503976954 -0.0770232069  0.3297510120 -1.6147260859
##  [471] -2.2142655939 -0.4531774472  0.2126608514 -1.6625439141 -0.1066374037
##  [476]  0.6669063000  0.4670343968 -0.0456118240  0.3432189637 -0.1140250911
##  [481] -0.9361015435 -0.0615343219 -0.2023942772  0.1450926007  0.2470673438
##  [486] -1.2831981208 -0.0850636513 -0.1231994516  0.1881586715  0.3853986685
##  [491] -0.0517253728  0.0689862698  0.6557035630  0.3720721381  0.5708940423
##  [496]  0.0914911745  0.2873085962 -0.5078527313 -0.8394960607  0.8879552963
##  [501]  0.2834456296 -1.9496340911  0.1028159987 -0.8715841138 -0.8155260063
##  [506]  0.4888736762  0.3882664663  0.3575632815 -0.0773290537  0.1793643472
##  [511]  0.4125972152 -0.0136734873 -0.8113443976 -0.4250394469  1.4511544988
##  [516]  0.0710603634  0.9202375050  0.0836897187  0.5483415915  0.1183474352
##  [521]  0.3060548557  0.4003496688  0.6439696675  0.0443587735 -0.3613186518
##  [526] -0.1941873259 -0.8797966271 -0.5623670592  0.2049722000  0.2448285347
##  [531] -0.1143896871  0.2704573660 -0.1902468356  0.4464352554  0.6104718869
##  [536]  0.4562123448 -1.2359846358  1.1764591472  0.3921805558  0.4908013529
##  [541]  0.5782142258 -0.1960827083  0.0487744089  0.0903484146 -0.0682532291
##  [546] -1.3199217286  0.7108335450 -0.0326039363  0.2802869416  0.0636713478
##  [551]  0.3070509107 -1.3922313823  0.5657158275 -1.5879061091  0.6365368043
##  [556]  0.2035079499 -1.0206826076 -0.4217375642  0.3499963215 -0.9195931942
##  [561]  0.5335740699  0.5251660505 -0.0289417038 -0.9177314827  0.2485747767
##  [566] -0.1457144725 -0.1842468694  0.3996694328  0.0302596251  0.0678890312
##  [571]  0.0901463217 -0.7692769654  0.0820601205 -0.0240022536 -0.6922781318
##  [576]  1.0737176221  0.2367268770 -0.3167924379  0.1550900035 -1.0532713370
##  [581]  0.1858893030 -0.7843480077 -0.0118198730  0.2318216740  0.1940999269
##  [586] -0.0165271988 -0.0811010394  0.3317056169 -0.2469358820  0.2020111977
##  [591] -0.0877675039  0.4496031927  0.4565790068  0.2010866281  0.1174437975
##  [596]  0.0193263699 -0.0766130557 -0.3525578540 -0.6694763801 -0.0385421720
##  [601] -0.4813515361  0.2016468121  0.1394312334 -1.3199217286 -0.1921906986
##  [606]  0.0359537456  0.3856634525  0.0167388496  0.3684332352 -0.0583897334
##  [611]  0.2433141770  1.1097583440 -0.4535510914  0.3409137686 -2.0521398361
##  [616]  0.0652265131  0.4523567778 -0.3706436786 -0.2350338348  0.6524468307
##  [621] -0.0126781530 -0.8183567232 -0.0327658644 -1.2068615828 -1.9124208379
##  [626] -2.4878423889  0.2334337605 -0.5817734055 -0.0256669555  0.3141171404
##  [631]  0.2795603532 -0.0600296637 -0.1766709294 -0.9432724644  0.1350607063
##  [636]  0.0536980775  0.7111754752  0.6581399967  1.2327945725 -1.3722962409
##  [641] -0.1887104345  2.1476420471  0.2649287815 -0.5401877046 -1.2904663913
##  [646] -0.9425609904  0.1432589438  0.3468298048  0.3795705932  0.7039346008
##  [651] -0.7270917239  0.0565749287  0.7731004142 -1.4901639790 -0.8825523607
##  [656] -0.1384415236 -0.0834225585  0.2306457264  0.5143447352  0.0793074721
##  [661]  0.1280940507  0.1930166561  0.5983895988  0.4709306163 -1.1690191710
##  [666]  0.3231320998 -0.9459392388  0.0705408198  0.1645966959 -1.1392902370
##  [671] -0.5362918806  0.2490528304  0.2923094266  0.3383338097  0.1245164730
##  [676]  0.4036984449 -0.7033415934 -0.4130271877  0.1956269099 -0.7761498781
##  [681] -0.0028024291  0.0514498987  0.3038566227  0.5395720500 -0.0625168563
##  [686]  0.0114494714  0.2887993750  0.2441235538 -1.3489540432  0.2865941841
##  [691]  0.0524417365 -1.4700181508 -0.1269149883  0.5799459805 -0.3474568381
##  [696] -0.5151772176  0.0673717594 -0.3013277125  0.9125120862 -1.2401648353
##  [701] -1.4140684949  0.1143858445 -0.4995742464 -0.0327255968  0.8892441336
##  [706] -0.4637609580  0.1442813319 -0.2461362630 -0.0297820705  0.0300470185
##  [711] -0.7755778095 -0.0843500251  0.1549875324  0.2142323945  0.3222078588
##  [716]  0.2908601390 -0.0750673197 -1.6072961699  0.3691404740 -0.1615991981
##  [721]  0.0985484211 -1.4712045436 -0.0972177826 -0.7871788803  0.5380358844
##  [726]  0.2935883793 -0.7836116366 -0.4667965500 -0.2324373545  0.5811939594
##  [731] -0.9063177975 -0.1090748309  1.0819464724  0.2172521549 -0.5761787801
##  [736]  0.7532130774 -0.6962750490 -0.1745557995  0.1982570113  0.3860891278
##  [741] -1.4017786656  0.2569836522 -0.0782933055  0.8582829537  0.0927990413
##  [746] -0.3251069735 -0.4891076433  0.5011257722  0.5820504271 -0.1077800364
##  [751] -2.2179106792 -0.0002931578 -0.1688788020 -0.1348359322 -0.9011445146
##  [756]  0.3876316508  0.5791327860  0.5132242947  0.1636242909  0.3655032468
##  [761]  0.1155107179  0.4750795359  0.1237806616 -0.0846242413 -0.1810138291
##  [766]  0.9093225689  0.1535178607  0.2098548174  0.3857550824  0.6346257352
##  [771] -0.4951195342  0.3873967387  0.0315329298  0.9452604651 -0.0690665577
##  [776]  0.0798800747 -0.8346274896  0.6970826907 -0.1182888675  0.5516210315
##  [781] -0.2540092406  0.3714602274  0.3122874495  0.2962551368 -1.3522637061
##  [786]  0.0463717024 -0.3387895010 -0.1156877213 -2.0513915622 -0.4048977706
##  [791] -0.5799674871  0.4323059502  0.0463238136  0.8268280691  0.3104020465
##  [796] -0.8572213865  0.0328530209  0.0057487635 -0.0274096890 -2.0596344527
##  [801]  0.5466287123  0.8553588464 -0.1575757510 -0.9198761567 -0.0326035421
##  [806]  0.0088486868  0.1081262575  0.0593648066 -0.3137143762  0.1316488766
##  [811] -0.1682549214  0.0587498368 -0.0673731737 -0.1427207170 -2.0534609608
##  [816]  0.1578840005  0.5394344188 -0.2372400329  2.1789376158  0.9336043736
##  [821]  0.1468059916 -0.2843021120 -0.1370318111  0.5381089286  0.2334736395
##  [826]  0.2696906076 -1.2729773663 -0.8203805468 -0.2570515555 -0.4765001193
##  [831]  0.2619405353  0.8479063204 -0.9194159012 -1.4422257341  2.1212850078
##  [836] -0.2855319900 -0.0377890544 -0.2342718593 -0.3154791661  0.3649628601
##  [841]  0.6732728361 -0.9673726398  0.0702700561 -0.0813968793  0.0563142753
##  [846] -0.0752620773 -0.3938282420 -0.7486557009 -0.8552453318  0.1492860464
##  [851] -0.3142928386 -0.0056639479  0.3074311896 -0.7076965803  0.6984848366
##  [856] -0.3733067192 -0.1871333532 -0.4225139695  0.9166088442  0.8682289193
##  [861] -1.9809322160 -0.8510685299 -0.1014375319  0.8465196219  0.0991382887
##  [866] -0.7363345566  0.3540303995 -0.3013277125  0.2388746846  0.0238927222
##  [871]  0.0909543363 -1.3487263437 -0.5798931865  0.2175005957  0.3992806189
##  [876]  0.4949073406 -0.5142462515  0.3782408193  0.0878613865 -0.0943392046
##  [881]  0.7532971793  0.2770059593 -0.0290151799 -0.1601851158  0.2738756663
##  [886]  0.2168344340 -0.0962102090  0.2599261533  0.3013166790  0.4621385957
##  [891]  0.1826681019 -1.2639106984 -0.6351385142  0.5912162919  0.6667056112
##  [896] -0.0518825580  0.7244201858 -0.8801475693  0.2917568421  0.1522325412
##  [901]  0.0163557547 -0.7291320362  0.1416188682 -1.5943548498  0.6678555749
##  [906] -1.2612088475 -2.1638201742 -0.9935135052 -0.0032188496  0.2957185612
##  [911] -0.6061481042  0.1728244132 -0.9220227835  0.5716104337  0.5525805379
##  [916]  1.3733230891 -0.9338471439 -0.3863616723  0.2481756922 -0.8889627769
##  [921] -0.1480657258 -0.3372194570 -0.1474509402  0.1646969341  0.0470316498
##  [926] -0.1554324380  2.1770727397  0.1040693620 -0.1493240090 -0.8676099369
##  [931] -2.3622354633  0.3397077165  0.7618844850  0.1433524149 -1.3866080968
##  [936] -0.5163511933 -0.9604672716  0.0741829895 -2.1261607831  0.2807876288
##  [941] -1.2361941291  0.3086476031  0.2185487243 -0.2806475661 -0.2157129303
##  [946]  0.7179430997  0.0189838584 -0.0372072252  0.6501630437  0.1812133861
##  [951] -1.4327136074  0.1743332591  0.9881739872 -0.4496111397 -0.0331441267
##  [956]  0.1579737509 -1.2374199264  0.5771124949  0.7045947541  0.6101930639
##  [961] -0.5139219108 -0.0354323341  0.3946000533  0.3002097989 -0.0074150191
##  [966]  0.3106597105 -0.4665938001 -0.4568013352 -0.6989564329  0.5990766012
##  [971] -0.7207670158  0.3785888681 -0.0758235404 -0.0675118813 -0.8938062786
##  [976]  0.2585793612  0.4386019686 -0.9515624133 -1.0639953174  0.3051099606
##  [981]  0.3077312506  0.6687677239  0.2464448032 -0.9179287876  0.0678767110
##  [986] -0.0167034454 -0.2578009156 -1.3600183143  0.0003637105  0.0678182865
##  [991] -0.0940920899  0.1487805601  0.1991803557  0.1328808725  0.8632666459
##  [996] -0.8938609021  0.5535262527  0.5474164802 -0.1128724663 -0.4801382161
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
##   0.63834901   0.19898779 
##  (0.06292546) (0.04449137)
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
## [1]  0.62358170 -0.48956643 -0.15081062 -0.89153424 -1.04018621 -0.03817116
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
## [1] -0.0299
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9247005
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
## t1*      4.5 0.01011011   0.9317533
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 2 6 8 9 
## 2 4 1 1 2
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
## [1] -2e-04
```

```r
se.boot
```

```
## [1] 0.9203407
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

