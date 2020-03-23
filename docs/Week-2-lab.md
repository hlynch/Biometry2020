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
## 0 1 3 4 5 6 8 
## 1 1 2 1 1 2 2
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
## [1] 9e-04
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
## [1] 2.70174
```

```r
UL.boot
```

```
## [1] 6.30006
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
##    [1] 4.3 3.6 5.3 4.9 4.5 4.1 4.7 6.1 4.4 6.2 4.0 4.6 4.3 3.4 4.9 2.5 5.5 4.0
##   [19] 4.9 5.7 4.1 3.8 4.0 4.2 3.5 4.6 2.1 4.5 4.8 4.8 6.1 4.8 3.7 2.8 5.8 4.5
##   [37] 5.0 2.4 5.0 3.2 5.1 3.9 4.4 5.6 4.2 4.2 5.2 3.8 3.8 5.0 4.6 3.9 5.3 5.4
##   [55] 4.9 4.2 4.2 3.6 3.4 4.6 4.0 3.9 5.7 3.6 3.9 6.5 3.5 5.0 3.9 3.7 4.6 6.4
##   [73] 2.7 4.2 5.2 5.1 5.2 4.8 4.2 4.8 4.8 2.7 3.9 5.5 4.0 5.2 4.7 5.5 4.9 3.1
##   [91] 3.6 3.8 5.2 6.5 5.6 4.5 4.6 3.5 4.0 6.5 5.5 5.3 4.7 3.0 3.5 4.9 5.3 6.0
##  [109] 5.4 3.3 6.4 3.9 4.2 4.5 3.9 5.4 3.9 5.2 4.3 4.5 4.1 5.5 3.5 3.8 6.0 4.4
##  [127] 5.1 4.8 5.5 5.7 2.4 3.9 5.7 6.2 6.2 3.8 5.8 3.8 4.7 5.1 4.3 4.9 3.5 2.7
##  [145] 3.3 5.8 4.3 4.6 2.7 4.7 6.0 4.3 1.8 5.3 4.9 5.0 5.4 3.6 5.2 5.6 2.0 3.9
##  [163] 3.3 4.1 3.2 3.2 5.3 5.2 3.7 3.0 4.8 5.0 3.1 3.4 4.4 5.8 1.9 4.2 4.5 4.4
##  [181] 4.5 3.0 4.4 4.8 3.1 4.2 3.3 5.1 5.3 4.0 4.1 3.1 3.3 4.6 5.0 4.4 4.2 4.0
##  [199] 6.0 3.8 5.2 3.4 3.4 5.5 3.8 4.8 4.0 4.3 3.7 5.8 6.2 3.9 4.9 5.5 3.9 3.2
##  [217] 4.1 3.1 4.2 4.7 4.3 3.9 4.5 6.0 4.8 3.4 4.9 5.0 5.0 4.8 3.6 3.9 6.2 2.9
##  [235] 4.6 4.5 5.2 6.2 5.4 4.3 4.1 4.5 4.2 6.2 4.5 5.3 4.5 5.4 3.5 3.0 4.7 2.2
##  [253] 5.0 4.8 5.1 5.5 3.7 5.3 5.4 5.4 3.4 4.4 4.5 4.5 4.0 5.1 3.1 4.8 3.9 5.7
##  [271] 4.3 3.1 3.6 5.4 3.6 4.4 3.8 3.1 5.1 3.4 4.7 4.4 6.4 4.4 4.7 5.9 3.7 4.4
##  [289] 3.3 5.4 5.2 3.9 2.5 2.8 5.5 4.1 4.6 3.9 3.8 3.9 3.3 3.1 3.8 5.0 4.1 3.9
##  [307] 2.3 4.7 5.6 4.8 2.2 3.8 4.1 4.2 4.2 3.9 4.1 5.0 3.8 5.1 5.5 4.8 3.4 5.0
##  [325] 3.6 4.9 3.6 4.3 4.1 6.2 5.0 5.3 2.7 4.5 4.5 6.6 4.5 5.3 4.5 7.3 6.5 3.0
##  [343] 2.5 4.6 4.7 5.6 5.9 6.5 3.4 5.4 4.3 2.9 4.9 4.2 3.4 3.8 4.9 5.4 4.6 4.3
##  [361] 4.5 5.7 3.4 4.6 3.3 5.8 3.4 4.5 5.6 7.2 5.2 4.2 4.6 5.8 1.8 5.2 4.4 6.8
##  [379] 3.2 3.3 6.0 5.7 5.2 5.5 5.4 3.8 4.0 2.9 2.6 5.9 5.1 3.8 4.6 4.7 4.7 4.9
##  [397] 4.3 4.6 2.9 4.2 4.7 3.9 5.8 3.9 4.3 3.9 3.9 4.5 2.9 3.6 2.6 5.0 4.9 4.5
##  [415] 5.1 2.8 5.6 4.8 5.6 3.8 5.2 3.9 3.1 4.2 3.8 5.3 3.6 3.4 3.8 3.4 5.1 4.6
##  [433] 3.9 4.5 3.9 3.4 3.6 4.2 4.9 5.8 5.1 5.5 2.5 5.3 4.7 2.7 5.0 4.9 5.5 3.8
##  [451] 4.8 5.2 5.6 4.3 5.8 4.7 3.7 4.8 5.7 5.3 4.9 4.4 6.0 4.7 6.7 3.8 2.7 4.6
##  [469] 3.2 3.6 2.3 4.7 5.1 3.9 5.0 3.9 5.8 4.3 2.6 4.4 3.7 4.7 4.1 3.9 5.1 4.9
##  [487] 4.2 5.5 6.2 5.1 5.0 5.1 3.2 5.5 6.3 4.0 4.5 5.7 5.4 5.0 4.4 5.0 3.9 3.9
##  [505] 3.8 4.2 6.0 4.9 5.1 5.0 4.6 4.5 4.4 4.0 5.2 5.5 3.4 4.8 4.1 5.3 4.5 4.7
##  [523] 3.5 6.1 4.8 3.0 4.0 4.5 5.0 3.6 3.5 3.6 3.8 4.7 4.6 3.7 5.3 3.8 3.8 4.9
##  [541] 5.0 4.5 3.9 4.3 5.1 3.2 3.3 4.7 3.4 5.1 4.6 6.7 4.4 3.8 5.1 3.8 5.2 3.6
##  [559] 4.9 5.8 3.2 4.0 5.0 5.5 5.1 3.4 4.8 4.8 4.9 4.5 4.0 4.6 2.9 4.5 4.0 5.6
##  [577] 5.7 4.3 3.7 4.4 3.9 3.2 4.2 3.5 5.0 6.1 3.1 4.3 3.8 5.8 3.9 4.8 4.7 5.3
##  [595] 4.9 4.8 5.6 1.8 3.8 3.7 4.0 4.8 5.2 4.1 5.3 4.9 4.5 3.7 5.4 3.5 5.7 5.0
##  [613] 4.4 5.1 6.3 4.6 3.4 4.2 3.2 3.1 5.1 5.3 3.9 6.0 4.3 3.0 6.1 2.3 4.9 4.5
##  [631] 3.8 4.5 4.5 4.1 4.2 4.8 6.8 6.3 3.4 3.4 3.2 3.4 5.4 2.5 3.8 4.3 2.3 4.4
##  [649] 3.9 3.4 5.5 5.3 4.3 4.0 4.4 5.1 5.8 5.3 4.0 3.8 5.3 5.2 2.9 4.1 4.3 4.6
##  [667] 4.3 4.6 5.7 4.2 4.6 3.9 4.1 3.9 4.1 3.9 3.2 4.3 4.8 3.0 5.2 5.8 3.1 4.6
##  [685] 4.4 3.6 5.2 4.1 4.9 3.9 5.3 5.1 4.4 5.5 3.9 4.2 5.7 3.3 6.3 6.2 3.9 4.0
##  [703] 4.4 5.9 5.0 4.1 3.5 5.4 4.7 2.7 2.7 3.9 4.7 3.8 4.6 4.4 4.1 3.2 4.7 5.3
##  [721] 4.4 4.7 5.7 5.5 5.5 4.7 6.5 4.3 4.3 5.0 6.3 3.7 5.6 5.1 5.7 4.8 5.5 4.8
##  [739] 4.6 4.9 4.0 4.3 3.4 4.5 5.7 3.0 4.3 4.8 4.8 3.5 5.3 3.6 4.3 3.5 5.0 5.1
##  [757] 4.1 4.7 4.6 4.7 3.3 3.6 3.8 5.6 6.3 5.3 4.7 3.6 4.8 3.4 4.5 3.6 3.5 5.3
##  [775] 4.4 4.9 4.6 4.1 4.0 5.6 6.1 5.4 5.9 4.9 4.5 5.8 3.5 4.5 3.8 4.4 5.5 4.7
##  [793] 2.4 3.7 5.7 4.5 3.4 3.5 4.5 4.5 6.6 4.8 3.0 3.6 5.1 4.1 4.0 4.3 4.9 4.8
##  [811] 5.6 5.3 3.7 4.4 3.6 4.9 4.0 3.6 4.0 4.6 6.2 5.0 4.8 3.8 5.5 4.7 5.5 4.1
##  [829] 3.9 4.3 6.1 5.2 5.3 6.4 4.3 4.0 2.8 4.3 4.8 4.2 3.8 5.7 2.4 5.3 5.2 3.1
##  [847] 3.9 4.3 4.7 4.5 6.5 5.6 4.7 3.7 5.1 5.3 2.5 3.5 3.0 3.8 3.3 2.2 4.3 3.3
##  [865] 6.9 4.3 3.9 5.8 5.4 4.7 4.2 4.5 4.7 3.3 5.1 4.8 6.8 3.0 4.7 4.3 5.6 4.3
##  [883] 4.1 4.7 6.0 3.6 4.8 4.2 6.0 6.7 4.0 4.6 2.9 3.9 5.5 6.7 4.4 2.8 4.2 3.5
##  [901] 6.0 3.8 5.9 6.3 4.8 7.0 4.1 4.2 4.3 3.8 4.9 4.2 5.7 4.3 4.7 6.0 4.3 4.4
##  [919] 4.1 4.4 3.6 4.3 5.1 5.6 4.2 3.7 5.1 5.9 3.8 3.9 4.7 4.8 4.9 5.6 5.1 3.7
##  [937] 4.5 5.6 5.7 5.4 5.0 4.2 5.2 3.3 3.8 4.9 5.3 4.5 2.9 4.8 3.6 4.3 4.7 5.7
##  [955] 3.9 4.5 4.7 4.1 5.3 4.5 4.4 3.8 6.2 4.7 5.3 6.3 4.9 2.9 5.0 5.3 4.3 3.0
##  [973] 4.3 4.4 4.6 4.4 4.2 3.3 4.8 4.1 3.5 2.8 4.3 4.9 4.5 4.5 5.6 4.4 4.6 2.9
##  [991] 5.1 3.9 4.4 4.6 5.1 4.3 4.9 3.7 4.1 5.5
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
##   2.6   6.3
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
##    [1] 4.2 5.1 5.0 3.9 4.0 4.5 4.4 3.0 3.3 2.7 4.8 3.3 5.6 6.8 4.1 4.3 3.2 3.9
##   [19] 6.6 2.3 5.0 3.7 3.2 5.0 4.1 5.6 4.6 3.9 3.7 4.2 4.6 2.8 4.0 4.5 3.8 4.2
##   [37] 3.6 4.7 5.1 3.9 6.2 3.3 4.5 4.9 4.9 5.5 5.4 5.5 5.2 3.5 4.4 5.0 4.0 4.3
##   [55] 4.7 4.7 4.6 4.6 5.5 4.7 4.7 4.9 5.3 4.3 4.8 3.8 3.1 3.4 5.2 6.0 4.9 4.4
##   [73] 3.8 5.5 3.7 2.9 4.4 5.2 5.9 2.6 6.0 3.4 4.4 3.9 3.7 4.5 4.9 4.7 4.1 5.1
##   [91] 3.0 5.4 5.2 4.9 4.7 4.9 4.5 5.3 4.3 4.1 5.1 3.5 3.9 5.5 5.3 6.0 5.9 5.0
##  [109] 5.1 5.9 3.4 5.1 4.3 3.6 5.6 3.2 5.3 5.1 4.7 2.6 4.6 5.7 3.3 5.0 5.3 4.2
##  [127] 4.3 4.8 4.8 3.8 4.2 4.8 4.6 4.7 4.8 4.9 4.3 2.7 3.8 5.3 3.8 4.7 5.2 3.7
##  [145] 4.3 4.4 4.2 3.8 5.5 4.4 6.0 6.1 4.7 5.4 4.3 5.3 4.5 4.0 5.1 4.9 5.1 5.3
##  [163] 4.8 6.5 4.5 2.6 4.0 6.0 5.0 3.5 4.1 5.8 6.0 4.1 5.1 2.6 5.7 4.1 4.6 4.3
##  [181] 2.6 4.5 3.9 3.9 5.1 4.8 4.4 6.0 5.6 5.0 5.9 4.5 4.0 4.8 5.1 4.4 5.9 5.6
##  [199] 5.3 4.9 5.0 5.1 5.6 5.0 3.9 3.5 3.2 3.9 5.2 5.0 3.5 5.9 4.2 6.5 4.6 5.6
##  [217] 4.8 4.2 6.8 4.7 4.0 5.4 3.0 3.4 4.2 6.7 5.4 5.6 4.3 5.1 5.2 5.9 6.0 5.2
##  [235] 4.6 5.1 5.1 5.1 2.2 2.9 4.4 3.7 4.8 5.0 3.8 4.4 5.6 3.1 5.4 4.5 6.4 5.5
##  [253] 3.4 6.1 3.7 5.2 3.8 3.4 3.8 5.9 3.6 3.1 3.8 4.8 3.5 5.1 4.0 3.0 6.1 4.6
##  [271] 4.8 5.0 6.1 4.2 3.0 6.0 6.6 3.2 5.2 5.3 5.6 5.3 5.7 5.0 3.4 3.8 4.9 4.6
##  [289] 4.5 6.2 4.6 2.9 3.8 4.9 4.2 4.1 4.3 3.7 5.1 4.0 5.6 4.7 4.7 4.0 4.7 5.0
##  [307] 4.5 4.1 4.3 2.6 3.4 4.8 4.7 3.6 5.2 5.7 2.5 5.6 5.3 4.3 3.1 4.6 5.7 5.4
##  [325] 2.0 4.9 5.5 4.3 4.0 5.4 4.3 4.0 5.0 3.2 4.6 4.7 5.1 4.4 4.4 5.6 4.2 3.7
##  [343] 5.4 5.4 3.3 5.6 6.6 4.3 4.8 4.2 2.8 4.8 5.3 3.7 6.6 3.7 4.4 4.4 5.4 4.1
##  [361] 4.1 3.7 4.0 3.4 4.6 4.5 4.4 4.2 5.5 4.8 3.5 3.4 5.2 4.5 5.1 3.3 2.6 5.7
##  [379] 5.2 3.9 4.1 4.7 5.9 4.9 4.2 3.2 4.5 5.3 2.3 5.1 6.2 6.1 5.8 3.3 4.3 3.3
##  [397] 5.2 3.7 4.4 4.8 2.8 3.4 5.1 5.3 4.6 4.9 4.8 4.8 3.2 3.1 4.2 4.7 5.6 4.3
##  [415] 4.5 3.6 4.7 4.9 4.4 5.4 4.1 6.0 4.5 3.8 3.6 5.3 2.4 4.6 5.2 4.3 5.4 3.9
##  [433] 4.7 4.9 4.8 5.5 5.4 5.5 5.1 3.4 4.9 4.1 4.4 5.3 4.7 2.5 5.0 3.3 5.4 5.5
##  [451] 5.5 4.9 4.4 4.4 5.0 5.3 4.9 3.4 4.8 4.4 3.5 5.6 4.7 4.1 3.8 4.7 5.2 5.2
##  [469] 4.6 4.6 4.7 5.5 3.2 4.6 3.7 3.4 4.4 7.0 5.5 4.0 4.9 5.2 3.6 4.1 4.6 5.5
##  [487] 3.5 3.8 4.4 6.1 4.8 3.9 4.9 6.6 5.0 3.9 4.7 5.1 4.6 4.0 3.7 1.6 2.8 4.0
##  [505] 3.2 5.1 2.1 5.5 3.9 7.2 3.8 4.2 4.2 4.7 4.5 5.3 4.0 4.3 5.7 6.3 4.6 3.7
##  [523] 4.6 4.8 4.1 4.4 3.4 5.5 4.5 4.0 4.9 4.3 5.8 3.8 4.9 4.0 4.4 5.7 4.9 2.4
##  [541] 3.1 3.2 3.2 3.4 4.9 4.2 5.8 4.1 4.8 3.9 3.3 4.9 3.1 4.3 5.0 5.7 3.5 3.9
##  [559] 4.4 4.3 3.5 4.3 4.7 2.8 6.5 3.5 5.1 2.7 3.8 4.6 3.9 4.2 3.6 3.2 4.9 4.5
##  [577] 4.0 4.7 4.8 4.8 5.1 5.2 5.2 4.6 4.6 4.6 4.8 6.9 3.3 5.8 3.8 4.9 6.3 5.2
##  [595] 5.8 3.6 3.9 4.5 5.2 3.6 5.9 5.0 5.2 4.7 5.7 4.1 6.5 2.6 4.7 3.7 6.3 4.8
##  [613] 3.9 3.7 5.3 4.6 4.6 4.4 5.4 5.3 6.3 4.9 4.4 4.4 3.7 4.2 3.5 5.1 4.3 5.1
##  [631] 6.6 4.6 5.6 3.8 5.2 4.6 4.9 4.9 4.7 3.3 3.6 4.0 5.0 3.2 3.9 3.6 3.5 5.7
##  [649] 4.8 3.8 6.3 4.8 3.7 4.9 4.5 3.8 5.2 5.9 3.2 3.7 4.0 3.0 2.7 3.9 2.7 3.8
##  [667] 3.8 3.7 4.1 4.6 4.7 3.5 5.4 3.9 4.5 5.0 4.8 4.5 5.3 5.4 3.4 4.9 6.1 4.9
##  [685] 4.1 5.4 7.0 3.2 4.5 5.1 4.4 6.1 4.5 6.1 6.0 5.1 5.2 5.1 5.7 5.4 2.7 3.8
##  [703] 4.9 4.6 5.3 5.0 5.2 4.7 5.3 4.1 2.7 4.9 3.8 4.5 4.7 5.1 5.5 4.0 6.0 5.9
##  [721] 5.9 2.9 4.2 1.7 4.3 3.6 3.6 3.5 3.8 4.9 3.5 4.4 4.5 5.2 4.8 5.7 4.8 3.8
##  [739] 4.8 4.8 5.0 3.6 5.5 5.9 5.7 4.9 4.3 5.9 4.5 3.1 3.8 4.5 5.2 4.7 4.4 4.4
##  [757] 3.8 4.3 5.6 5.1 4.0 5.2 4.7 4.5 4.2 5.2 3.8 5.1 4.0 5.1 3.8 4.9 3.1 5.5
##  [775] 5.3 4.4 4.7 4.8 5.3 4.4 4.6 3.4 4.9 2.8 5.4 4.9 4.7 4.3 4.3 4.9 4.8 4.6
##  [793] 5.3 3.7 3.3 4.8 3.7 5.0 4.1 4.6 4.8 5.1 4.4 4.5 3.8 3.2 4.2 4.5 4.2 4.0
##  [811] 3.0 4.8 4.0 3.9 5.2 3.9 4.0 6.6 4.0 2.7 4.6 4.0 4.6 5.5 4.9 5.6 5.4 4.4
##  [829] 4.0 3.7 4.3 4.6 6.0 3.0 5.3 4.6 3.2 4.7 6.3 3.3 4.2 4.3 5.5 3.7 4.1 3.6
##  [847] 4.1 4.5 3.0 5.3 2.8 3.8 5.2 4.0 4.5 5.5 3.8 4.8 4.6 4.1 3.9 6.4 4.1 4.1
##  [865] 5.1 5.2 5.0 4.5 4.4 5.4 3.5 5.7 5.7 5.0 5.0 3.2 5.7 5.3 4.1 4.8 2.9 5.5
##  [883] 3.3 6.1 3.9 4.8 6.3 4.4 3.8 4.4 3.9 4.0 4.1 5.0 2.1 5.3 3.0 3.5 2.3 4.7
##  [901] 5.0 4.4 3.3 2.7 4.3 4.4 5.2 4.6 3.5 6.5 5.0 5.1 3.8 4.2 5.7 5.5 4.7 3.3
##  [919] 5.4 5.4 3.7 2.6 4.4 4.7 4.0 3.6 3.1 5.8 4.6 4.5 3.5 3.5 4.8 4.5 4.7 3.2
##  [937] 4.9 4.7 4.5 5.7 5.6 3.1 3.7 3.3 4.3 2.9 3.3 5.9 4.8 4.5 4.8 4.7 3.5 5.4
##  [955] 6.4 4.6 4.5 4.0 4.1 4.7 4.9 5.9 4.3 5.4 5.0 5.7 4.8 6.0 4.5 3.1 4.8 3.8
##  [973] 5.8 3.8 5.9 4.4 4.1 3.3 6.4 3.8 3.6 4.8 4.3 4.9 4.2 4.7 4.3 5.1 5.5 3.9
##  [991] 5.3 3.1 5.0 4.1 4.8 4.4 5.1 2.5 3.3 3.4
## 
## $func.thetastar
## [1] 0.0261
## 
## $jack.boot.val
##  [1]  0.50835735  0.40989011  0.32893258  0.18221574  0.07657658  0.01478261
##  [7] -0.16270270 -0.26590909 -0.36264706 -0.51569767
## 
## $jack.boot.se
## [1] 0.9785633
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
##    [1] 4.1 4.3 5.9 4.1 5.9 4.6 4.2 4.8 4.0 3.9 4.0 3.5 5.4 4.8 5.0 5.1 4.5 4.6
##   [19] 5.5 5.4 3.6 4.3 5.3 4.7 4.9 4.5 4.6 4.4 6.0 4.4 4.7 6.0 5.6 5.4 3.3 3.8
##   [37] 3.6 4.9 4.5 4.7 5.2 4.3 3.2 3.3 6.1 3.1 4.4 3.4 4.1 2.5 4.1 3.1 5.7 3.1
##   [55] 4.6 4.8 4.0 3.9 5.2 6.0 5.2 5.3 4.2 3.1 5.0 5.7 3.7 2.9 4.2 6.3 5.0 2.8
##   [73] 4.8 6.0 6.7 3.6 4.9 5.1 5.5 4.1 3.7 4.1 4.2 4.3 5.0 2.9 3.4 4.1 3.4 5.0
##   [91] 4.3 3.4 4.8 5.4 4.4 3.7 4.2 5.6 2.9 5.3 4.1 3.4 3.2 4.9 4.8 6.2 5.6 3.4
##  [109] 4.5 3.6 4.8 5.7 5.8 4.1 4.9 2.8 4.6 4.0 5.2 3.9 4.6 4.3 4.6 5.3 5.3 4.8
##  [127] 4.7 5.3 4.0 5.7 6.2 3.8 2.8 4.5 4.9 5.9 4.8 5.6 4.3 2.2 5.2 3.8 3.1 4.6
##  [145] 5.4 3.1 4.0 5.2 5.2 3.8 4.8 4.3 3.6 3.9 3.9 6.5 4.4 3.5 2.8 3.5 3.6 3.8
##  [163] 5.7 3.4 4.1 4.4 5.5 4.6 4.9 2.4 4.8 4.3 3.7 3.5 5.0 3.9 3.6 5.3 4.4 3.4
##  [181] 4.4 4.4 6.3 5.0 4.3 4.8 6.6 5.2 5.0 4.7 3.8 3.9 5.1 5.2 4.8 5.3 5.0 3.6
##  [199] 5.3 4.7 3.3 5.7 5.4 4.5 5.6 5.6 5.3 4.4 4.8 4.9 4.7 4.8 4.7 5.4 4.3 4.7
##  [217] 5.6 4.8 4.9 4.0 5.9 5.2 4.1 2.6 7.0 3.5 5.0 4.7 6.1 2.7 4.0 5.5 5.0 4.0
##  [235] 3.6 4.8 5.8 2.8 4.7 5.4 3.8 3.6 5.2 3.8 3.8 5.7 3.8 4.0 4.7 5.8 4.0 4.6
##  [253] 5.1 4.4 3.2 4.9 5.0 3.6 4.6 4.8 4.7 4.5 3.7 4.7 3.5 3.2 3.7 4.5 4.3 7.0
##  [271] 4.2 5.2 3.3 4.3 6.1 4.1 4.8 4.8 5.3 4.0 4.4 5.7 4.2 4.8 2.9 4.1 3.7 5.3
##  [289] 5.4 5.2 4.8 2.4 4.6 5.0 5.6 4.5 4.1 4.3 5.0 3.9 5.3 3.3 4.4 6.3 5.8 4.8
##  [307] 6.3 3.3 4.3 4.0 4.7 6.5 5.3 3.5 4.9 5.7 3.6 3.8 5.7 4.3 3.3 4.7 3.8 4.8
##  [325] 4.7 4.5 4.0 4.6 4.6 4.2 4.4 5.3 5.4 4.2 3.7 5.4 5.4 4.8 3.7 4.7 7.5 6.2
##  [343] 5.7 5.8 3.8 4.0 4.0 4.8 4.3 5.3 4.5 5.4 3.9 4.9 4.0 3.6 3.5 5.0 5.3 4.4
##  [361] 5.9 3.2 3.9 4.7 3.9 5.5 5.3 5.1 4.0 4.9 3.7 2.0 5.5 4.9 2.8 5.2 5.6 3.2
##  [379] 3.5 3.9 4.2 3.6 4.6 4.0 4.9 5.8 5.5 4.5 4.2 3.6 4.6 4.5 4.0 2.6 7.2 4.3
##  [397] 4.6 5.9 5.0 6.7 3.9 2.6 6.8 3.3 4.4 4.7 2.5 4.7 4.8 5.0 5.4 2.8 4.1 2.9
##  [415] 3.9 2.9 4.6 4.0 5.3 6.0 4.5 5.1 3.7 5.7 5.2 4.8 3.6 5.1 3.8 4.5 3.1 5.8
##  [433] 4.3 6.1 5.2 3.7 4.1 4.5 4.3 4.1 3.2 5.3 6.6 4.6 5.2 4.1 4.5 3.4 3.8 4.5
##  [451] 4.2 5.5 3.1 4.6 5.4 5.8 2.9 5.5 5.6 5.3 3.7 3.6 3.0 3.9 5.5 5.0 4.8 4.5
##  [469] 4.7 5.0 4.3 5.5 4.3 4.4 4.7 6.7 5.4 4.9 4.5 4.4 5.1 4.2 4.3 4.7 3.7 4.8
##  [487] 3.6 4.3 4.3 5.2 4.3 5.0 6.0 5.2 5.8 4.1 5.7 4.7 4.4 4.7 5.5 3.6 3.9 3.8
##  [505] 4.9 3.8 6.1 4.2 4.5 5.3 4.3 5.2 5.4 4.6 4.4 2.8 3.4 4.5 4.7 3.3 4.0 4.7
##  [523] 4.2 5.0 4.5 3.7 3.8 5.6 3.8 4.4 5.9 4.0 5.1 6.2 3.7 4.1 5.5 3.6 5.4 3.8
##  [541] 4.6 3.9 2.6 4.6 5.8 6.3 4.7 2.6 2.5 3.8 6.3 4.1 4.6 3.5 5.4 4.9 4.2 6.2
##  [559] 3.8 5.4 5.9 4.4 5.3 4.7 4.7 4.5 5.1 3.6 3.9 3.5 4.8 5.7 4.9 4.2 3.3 6.7
##  [577] 5.0 5.1 2.5 4.3 3.1 4.4 5.1 5.0 3.4 4.8 5.2 5.0 5.4 3.4 4.8 4.8 4.4 5.9
##  [595] 5.5 3.0 3.6 3.5 4.7 5.2 4.1 3.7 5.6 3.9 4.6 3.7 4.6 3.2 5.2 3.8 5.6 5.2
##  [613] 4.6 5.2 2.8 4.8 3.3 5.8 5.2 4.0 2.8 3.5 3.3 5.2 3.9 4.6 3.8 4.5 3.9 4.4
##  [631] 3.7 4.0 4.1 4.7 4.2 4.1 4.2 4.1 3.4 5.1 5.6 4.1 4.0 4.5 3.9 3.6 5.3 4.2
##  [649] 1.9 3.9 4.7 4.4 6.6 4.7 4.4 4.0 5.3 4.9 3.7 3.0 3.8 4.1 3.8 4.9 3.2 4.2
##  [667] 5.4 2.5 5.1 4.1 5.1 4.4 4.7 3.4 4.7 7.0 3.4 4.7 4.8 3.4 5.5 5.3 4.3 4.5
##  [685] 3.7 3.8 4.1 3.0 3.4 4.3 5.1 4.4 3.9 3.1 4.4 3.2 3.8 4.5 5.3 3.8 4.2 5.6
##  [703] 4.8 4.3 6.2 4.1 4.3 4.4 4.1 2.9 3.3 4.6 5.2 4.5 3.2 4.7 3.8 6.1 4.8 4.7
##  [721] 3.8 4.3 3.7 4.7 5.2 5.5 4.6 3.7 3.7 3.6 4.1 5.1 4.1 5.1 4.1 5.0 4.7 5.3
##  [739] 3.2 5.5 6.2 5.1 5.4 3.2 3.6 5.4 4.4 3.6 4.6 5.3 4.6 4.4 4.8 4.7 4.4 4.7
##  [757] 4.5 4.0 5.3 3.2 4.5 4.1 5.9 4.2 6.0 4.7 4.1 3.8 5.2 4.9 5.2 3.6 3.7 3.4
##  [775] 4.8 4.5 4.5 4.5 4.4 4.1 4.0 3.9 4.8 5.3 5.2 6.2 4.1 5.5 4.3 5.0 5.3 4.2
##  [793] 5.3 3.3 4.5 5.0 4.3 5.2 4.7 5.7 5.1 5.1 2.0 3.5 4.9 4.7 3.6 6.2 4.7 3.1
##  [811] 5.2 2.3 4.0 4.3 4.1 5.6 4.3 6.0 5.5 3.4 5.5 3.4 4.8 6.2 3.8 4.8 4.3 3.9
##  [829] 3.7 5.7 5.6 5.2 5.4 5.5 6.1 4.1 3.5 3.3 4.1 5.5 4.2 4.4 5.8 4.5 3.0 4.3
##  [847] 4.4 3.4 5.8 5.1 4.9 3.1 3.3 4.0 4.5 3.7 4.6 4.8 5.0 3.2 4.2 6.1 4.4 4.9
##  [865] 5.3 5.2 4.2 3.3 6.2 3.8 5.5 3.6 4.4 4.5 3.5 7.0 4.1 3.6 3.5 4.7 5.7 3.9
##  [883] 4.0 4.0 4.5 4.3 5.5 4.6 5.5 4.9 4.7 5.0 3.0 4.7 5.6 4.7 5.3 2.5 3.4 2.6
##  [901] 5.6 5.6 3.5 4.9 5.3 4.4 4.0 3.6 4.5 4.7 4.2 3.8 3.7 5.0 5.4 4.3 4.9 5.7
##  [919] 3.7 4.1 5.1 5.1 4.2 3.6 5.3 3.6 4.0 4.2 3.0 4.2 3.2 3.4 6.2 3.9 4.7 3.1
##  [937] 4.2 5.9 3.5 5.5 3.6 3.1 4.0 4.0 4.2 5.9 3.9 3.7 2.7 4.0 4.4 4.1 5.4 6.7
##  [955] 4.4 5.7 4.0 4.7 4.0 4.5 3.8 4.4 5.4 2.4 4.7 3.7 2.8 4.6 6.1 4.5 2.9 4.2
##  [973] 2.8 5.4 5.1 4.9 2.9 3.2 3.9 5.0 4.8 3.7 4.4 5.4 5.3 5.0 4.3 3.6 3.7 6.5
##  [991] 4.2 5.0 6.1 4.7 3.6 5.2 4.4 3.7 3.9 4.0
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.4 5.4 5.3 5.3 5.1 4.8 4.8 4.7 4.5 4.4
## 
## $jack.boot.se
## [1] 1.073732
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
## [1] 0.320635
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
##    9.669401   19.239900 
##  ( 4.251775) ( 8.683371)
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
## [1]  0.668637473 -0.001216743 -0.092584375 -0.130665512  0.910914930
## [6]  0.071471573
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
##    [1] -1.332036901  0.044929159  1.072410475  0.909573856  0.560870685
##    [6]  0.591141680 -0.777252443 -1.136503319  0.856967192 -0.835130666
##   [11]  0.662172210  0.637539328  0.028308765  0.988321889 -1.414743200
##   [16]  0.955081861 -1.336090020  0.972753974  0.970649089  1.394431919
##   [21]  0.141835625 -0.284521444  0.262551715 -0.511255039  0.415101682
##   [26] -0.675388549  0.171078325  0.884409922 -0.082505509  0.891751123
##   [31]  0.008908579  0.637283800  0.322008190  0.644858510 -0.228125706
##   [36]  0.842009292  0.105675222  0.825936298  0.655698404  0.977335930
##   [41] -1.191085366  0.341098218 -0.799968118  0.285532992 -1.384499229
##   [46] -0.601056981 -0.435582593  0.474334760 -0.164503003  1.392278017
##   [51]  0.130282425 -0.541818234 -0.495546144 -0.070266143  0.050170260
##   [56] -0.149189413  1.692182267  0.019913537 -0.743917135 -0.229718449
##   [61]  0.329473471 -0.592023542  0.241203001  0.823626204 -0.483049277
##   [66]  0.872180062  0.872851823  1.145102336  1.175123204 -0.398347326
##   [71] -0.224357532 -0.034169046 -0.014926430  0.398262739  0.072069982
##   [76] -0.666034106 -1.293010787 -0.265654269  0.087660229 -0.002272248
##   [81] -0.084573877  0.496226219  0.100268351  0.612662216  0.548511984
##   [86]  0.857040911  0.269801425 -0.891398253  1.050782983 -0.158449341
##   [91]  0.867880391 -0.157205320  0.153552594  0.231860421 -0.443286238
##   [96] -0.845885926 -1.273564332  0.535515896  0.364790425  0.648174836
##  [101] -0.746149506 -1.386259567 -0.407629457  0.098358716 -0.106914277
##  [106]  0.587721401 -0.125752241 -0.290158749  0.269183320  0.665613023
##  [111]  0.512114022  0.109839804 -0.578323919  0.557397647  0.219400345
##  [116]  0.720659424 -0.423183489  0.996733654 -0.570509321  0.476432723
##  [121]  1.371224427 -0.633850151 -0.472794088  0.041679676 -0.417736003
##  [126]  0.325202630  0.261470411 -0.023060966 -1.553388237 -0.785336615
##  [131]  0.135698136 -0.054988074 -0.469686747 -0.239024656  0.505402332
##  [136]  0.953345082 -0.110283142  1.070691404  0.463613831  0.586415984
##  [141]  0.271495825  0.217164124  0.244825072  0.470959277 -0.099322761
##  [146]  0.945347397 -0.642500937 -0.154599764  0.722029781 -1.034979140
##  [151] -0.513355262  0.186977976  0.582170019  0.228707814  1.017504349
##  [156]  0.125390441  0.067443890  0.692632524 -1.241317717  1.066910385
##  [161]  0.147455696  0.013003692  0.306986465  0.074721374  1.091608005
##  [166]  0.492145687  1.160691572 -1.908205319 -0.285432555 -0.901363731
##  [171]  0.055167328 -0.326533705 -0.661376091  0.315526270  1.066625713
##  [176] -1.076454695 -1.775676971 -0.709646392  0.967605285  0.732526454
##  [181]  0.985875643  1.249214932  0.816859198  0.064613558 -0.504649453
##  [186]  0.092774745  0.018257039  0.307473346 -0.213232027 -0.264708500
##  [191]  0.233223380 -0.024247849 -0.399643980 -0.483336466  0.608852646
##  [196] -0.772800103  0.077614168 -0.411924644 -0.639918423  0.489345865
##  [201]  1.231128537  0.821463100  0.039590714 -0.733679790  0.360702484
##  [206] -0.453686385 -0.363855711  0.295721850 -0.029924554  0.612783075
##  [211]  0.480180677 -0.193511032  0.381049971  0.041746142  0.319737176
##  [216] -0.325692875  1.348271354  0.711497940 -0.916466580  0.403584244
##  [221] -1.029816044  0.244825072 -0.234963104  0.198248214  0.033383863
##  [226]  0.173961383 -0.086041362  0.460632021  0.120124931  0.177990149
##  [231] -0.436979135  0.580388790 -0.618940772  1.176474007 -0.938696903
##  [236] -1.262537225  1.173956229  0.455614345 -0.789575037 -0.664944753
##  [241]  0.928634350 -0.095610024  0.170026260  0.228676444  0.640058680
##  [246] -0.428316337  1.259067419 -0.777688994  0.094390474 -0.672069954
##  [251] -0.125762008  0.113109072  0.270072968  1.138376444 -0.723595788
##  [256]  0.440260244  1.552921555  0.504093634  1.616280421  0.699600171
##  [261]  0.442320951 -0.909226694  0.592962258 -1.007310470  1.343400891
##  [266] -0.198516908  1.256598497  0.859011206 -0.679305868 -0.518559729
##  [271]  1.743849104 -0.634674244  1.227855441  0.663450330  0.828832038
##  [276] -0.224558735  1.043831557  0.624761369 -1.215266067  0.411462810
##  [281] -0.049850358  0.954782960  0.184473411 -0.711213736  0.220378099
##  [286]  0.532164586  0.528649248 -0.451529277  1.276668773  0.372374098
##  [291] -0.682285211 -0.802586730 -0.250752742  0.029466466  0.276903129
##  [296]  0.180037030  1.153183513  1.336873320  0.393105107  0.167510006
##  [301]  0.161027851  0.891888685  0.459678199 -0.161497756 -0.169904704
##  [306]  0.743055699 -0.484851501 -0.733219543  0.347416341  0.011722559
##  [311]  0.943729230 -0.274261805  0.510263895  0.970445212 -0.168090050
##  [316] -0.212054880  0.905552174  0.132825548 -0.252157127  0.841912676
##  [321] -0.330315249 -0.902588960  0.564803287 -0.372873347  0.500902908
##  [326] -0.455269721  0.671095781 -0.280062421  1.207339899  0.204769573
##  [331]  1.781305803  0.757758896 -0.472794088 -0.574179244  0.724002596
##  [336]  1.544067982  0.429574921 -0.408461910  0.360179106 -0.623991342
##  [341]  1.170198720  0.473278608  0.164562163 -0.063266348  0.237345139
##  [346]  0.442238274  0.709367135  0.037298210  0.706405798 -0.791576344
##  [351]  0.290559606  0.018667083 -1.568243737  0.362231739 -1.356442514
##  [356]  0.403676950 -0.264346194  0.426519663  0.175414256 -0.589815622
##  [361] -0.045211458  0.559568704  0.295628856  0.319132997  0.824485139
##  [366]  1.320974413  1.131529774  0.740079085 -0.319251188 -0.030195376
##  [371]  0.642662571  1.031202798 -0.391686831  1.053854960  0.091781720
##  [376]  1.249956580 -0.011056445  1.091396179  0.144985925  0.755185377
##  [381] -0.470255892 -0.513188336  1.124242055  0.415597176  0.114696590
##  [386]  0.333756356  0.341672085  0.389670560  0.936914019  1.449262338
##  [391]  1.278520744  0.279070410 -0.622324031  0.189840889 -0.002740240
##  [396]  0.531514999 -0.333962046  0.097561684  0.302731582 -1.598637868
##  [401] -0.041165728 -0.435168565 -0.570133982  0.807720098 -0.031376166
##  [406]  0.117977652  0.889678355  0.398026328  0.315810337 -0.373621700
##  [411]  0.432754250  0.249892262 -0.291194956 -0.021031649 -0.025782975
##  [416]  0.220900006  0.500214132  0.117998348  0.593069884  0.045114428
##  [421]  0.494541151  1.541009053 -0.750474271 -0.541155118  1.564505936
##  [426]  0.294179144 -0.030985427 -0.548493120  1.052490743 -0.204248694
##  [431]  0.455337150 -1.000128644  1.019945805  0.053308017  0.067062339
##  [436]  1.099641307 -0.264213490 -0.467651396  0.124291017  0.127800285
##  [441] -1.248936940 -0.651414198 -1.222069292  0.029462136  1.246097388
##  [446]  0.900947202  1.491321430 -1.112752113 -0.572290793 -0.826249400
##  [451]  0.666408519 -0.376680188  0.569865770  1.715982142 -0.054676616
##  [456] -0.636347531  0.283291144 -0.006667344 -0.394170308  0.004042862
##  [461] -0.297844893 -0.349772446  0.352163281  0.375026655  0.767344696
##  [466]  0.655698404  0.996880452  0.257226567 -0.010549421  1.172242432
##  [471]  0.620778110 -0.962435299 -0.125762008  0.188718314 -1.232823025
##  [476]  0.757928016 -0.751021155  0.456851047  0.218258070  0.137017271
##  [481] -0.243333524  0.741993277  0.388665238 -0.436515127 -0.058296828
##  [486]  0.259578477  0.204767617 -1.570349955 -1.170746787 -1.088541151
##  [491]  0.597374262  0.036994024 -0.134546412  0.682851344  1.578255235
##  [496] -0.077621269  0.251658460  0.378653292 -1.301794639  0.005957600
##  [501]  0.191518067  0.119440401  0.287058437 -0.348649034  0.697999309
##  [506] -0.153188076  0.604166360 -0.459961810 -0.283126830 -1.383338474
##  [511] -0.408069134  0.383137231  0.726061008 -0.982304142 -0.867009783
##  [516]  0.118854903  0.655773828  0.603013831  1.420802660 -0.601630116
##  [521]  1.229166386 -0.284691195 -1.160961907 -0.087798462 -0.962962876
##  [526] -0.107909300  0.381049971  0.837444013  1.053647748 -0.237883129
##  [531]  0.446449514 -0.289085033  0.657035646  0.461359149  0.717469357
##  [536]  1.316793157  0.090683094  0.252618508  0.391606949 -0.460517402
##  [541]  0.352672172  0.254298115  0.709989994  1.283886230 -0.338530277
##  [546]  0.180555107  0.623046923 -1.240238801 -1.336144999  0.100178778
##  [551]  0.092250230 -0.025213853 -0.091208761 -0.107410440  0.619036760
##  [556] -0.070627209 -0.708974324  0.774312415  0.837928913 -0.151003480
##  [561]  1.835677292  0.924948343  0.873178785  0.575829126 -0.647632118
##  [566] -0.677079524  0.333317460 -0.249181934  0.779020285  0.145319502
##  [571] -0.872703627 -0.767563977 -0.210680479  0.573460849 -1.643467841
##  [576] -0.622704270 -0.008378258  0.512989572 -1.430786682 -1.171566002
##  [581]  0.529453049 -0.860834017  0.580758822  0.705262531  0.949899508
##  [586]  0.799936737 -1.207839659  0.030010484  1.391561383  0.168053493
##  [591]  0.595116381  0.448235511  0.570380657  1.195721626 -0.918889128
##  [596]  0.101222804 -0.817733597 -0.004711855 -0.059289571  0.231582823
##  [601]  0.091008396  0.410272204 -0.500115114  1.709876969  0.359265991
##  [606]  1.545312764  1.057555413 -0.021189222  0.150311317  0.613798815
##  [611]  0.655193063  0.611173224  0.250240322  0.638932957 -0.447582006
##  [616]  0.539216121 -0.056172008  0.029082314  0.873144212 -0.094853776
##  [621]  0.315258553 -0.593254661 -0.631989078  0.377578130  0.409825556
##  [626]  0.822242715  0.056077319 -1.144607388  1.275844741  0.805110357
##  [631]  0.457522699  0.317375871  0.520029755  1.216602808  0.285334748
##  [636]  0.100178778 -1.788191921  0.284754336  0.178798878  0.698568344
##  [641] -1.249843310 -0.526747679  0.450201951  1.105511584 -0.110283142
##  [646] -0.511759094 -0.538542192 -1.253435470  0.216546310  0.302953497
##  [651] -0.932419736  0.933977232  0.313613175  0.445272307  0.178227284
##  [656] -1.368260906  0.699801882  1.259067419  0.158618783  0.349755774
##  [661] -0.093039454  1.164138902 -0.176108079  0.614263280  0.890699262
##  [666] -1.397255033  0.757431822 -0.011037607 -0.604353790 -1.111179401
##  [671] -1.216829756  0.410299174 -1.326362764 -0.538910447 -0.279977703
##  [676]  0.702329107 -0.150276191  0.885443771 -0.786604199 -1.409361715
##  [681]  1.117277899 -0.050848041 -0.093835031 -0.587176275  0.827055714
##  [686]  0.155734533  0.302731582  0.624312702 -0.539705433  1.285416150
##  [691] -0.327494969 -1.162005148 -0.905436186  0.951119665  0.371690826
##  [696]  0.447556640  0.310688237 -1.367143082 -0.578477031 -0.881285533
##  [701]  0.194182429  0.534410838  0.249768301  0.731541873 -0.438871127
##  [706]  0.849949222 -0.135202501 -0.154641188 -1.237160419 -0.463629179
##  [711] -0.124377847 -0.244250798  0.889678355  0.246465096  0.720960574
##  [716] -0.811880190 -0.651381562 -0.748339697 -0.971820607  0.102678686
##  [721]  0.435709140 -0.236152253  0.864705810 -0.432716502 -0.207086173
##  [726]  0.052531663  0.657947483  0.259524549  1.007620586 -0.059631037
##  [731]  0.676871203  1.284599932 -0.530092505  0.920891362  0.594337629
##  [736]  1.550054370  0.579682667  0.125598321  0.922069682 -0.267555798
##  [741]  1.312325889  0.825494204  1.617172942  0.157360654  0.563187646
##  [746] -0.258159392  0.268948001 -1.285699818  1.379967910  0.721527360
##  [751] -0.667502706  1.303812269  0.318566847 -0.456015698  0.339488034
##  [756] -0.424383743 -0.728831391 -1.332036901  0.189792454  0.295265231
##  [761] -0.686369411  0.712938815 -0.134191422 -1.438473543  0.021339481
##  [766] -1.205511352  0.677712889 -0.723500555  0.487731976  0.379719643
##  [771]  0.075297615 -0.614828164  0.271476562  0.295157337  0.200966261
##  [776] -0.602093073  0.194372109  1.294255463  0.986979611  0.214507041
##  [781] -0.880867034  0.926812002  0.525935710 -0.067287189  0.015581919
##  [786]  1.240480612 -0.557906642  1.365075221 -0.420579004  0.519649771
##  [791] -0.197410518 -0.158639722 -0.471993929 -0.376359316  0.890694143
##  [796]  0.574105891  0.772899292  0.865256720  1.071675963 -0.495546144
##  [801]  0.120487266  0.382891474 -0.416679567 -0.055388612  0.653408769
##  [806] -0.660206415 -0.157922769 -1.476556518 -0.396563682  0.198768269
##  [811]  0.452584317  0.557322071  1.390702576  1.179853353  0.443918327
##  [816]  0.136981068  0.736200623  0.308422052  0.476690337  0.925995964
##  [821]  0.648521013  0.803472555  0.441284715 -0.067169898 -0.240479729
##  [826]  1.247027319  0.853426589 -1.496057281  0.640058680 -0.924470175
##  [831]  0.763508753 -0.828197602 -0.603295314  0.100286864 -1.142542870
##  [836] -0.528373106 -0.202710493  0.249714575  0.896292567 -0.803909188
##  [841] -0.490067961  0.327119948  0.544886820  1.791335430  0.096614555
##  [846] -0.182771374  0.367064630 -0.491909704  0.988270225  0.063270109
##  [851]  0.121616376  0.564803112  0.756991485  0.785231598 -0.271095586
##  [856]  0.569405685  0.187179845  0.934837682  0.320635037  0.745037876
##  [861]  0.432849552  0.190160647  1.056646188  0.595646283 -0.266017302
##  [866]  0.167442118  0.365143839 -0.424988450 -1.643870967  1.094033812
##  [871]  0.883485291  0.158912031 -0.896078673 -0.896078673  0.135860090
##  [876] -0.818365756  0.008295828 -0.692666235  0.936782537 -0.233068588
##  [881] -0.187521515  0.669293292 -0.393323954  1.652060853  0.575742878
##  [886]  0.500721520 -0.444786846 -0.638002585 -0.237625426 -0.502971714
##  [891] -0.033196095 -0.476583218  0.440826197 -0.777329942  0.302242406
##  [896] -0.566207524 -0.188693953  1.022810310  0.084101113 -1.209592170
##  [901]  1.804074240  1.174288630  0.941861034  0.224950342 -0.138882107
##  [906]  0.697945524  0.607130591 -0.149506284  0.377408136  0.219918598
##  [911] -0.387019554  1.535750047  0.491130739  0.049195597 -0.441593785
##  [916] -0.536235896 -1.316430689 -0.085845032 -0.228864352  0.192259987
##  [921] -0.249849484 -0.341282582  0.814688215 -0.037773374  0.632710797
##  [926] -1.001103811 -0.998275325 -0.496953244  0.109052066  0.537040989
##  [931]  1.221750101  0.396857612 -0.265204344  0.188718314  0.387682706
##  [936] -0.036816118 -0.124754177 -0.212847573  1.507200816 -0.523023643
##  [941]  0.209878900  0.220722588  1.410208677 -0.587430152  0.038070248
##  [946] -1.404958519 -1.220036098  0.294611325  0.204674384  0.641273949
##  [951] -1.008886300  0.988753151  0.203124230  0.881963844 -0.088268209
##  [956]  0.141805093  0.312568072  0.014213006 -1.349994406  1.021741499
##  [961]  0.229288552 -0.813390817 -0.266627473  0.491937983 -0.810347566
##  [966]  0.145090323  0.822842383  0.511784850 -0.618009131 -1.002946274
##  [971]  1.207982691  0.158259721  0.482949700  0.565434606  1.154496575
##  [976]  1.576722281  0.028617717  0.996958206  0.479353202  1.095941244
##  [981] -0.967779689 -0.409880352  0.564164275  1.622049908  0.665001456
##  [986]  0.838516913  0.602816804 -0.121071441  0.659948495  0.573460849
##  [991]  0.401412597  1.093581314 -0.500457616  0.491433690  0.355761133
##  [996]  1.192548954 -0.138622066  0.224114046 -1.011357220 -0.001381419
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

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced
```

```r
fit2
```

```
##       mean          sd    
##   0.50254603   0.15528479 
##  (0.04910536) (0.03471603)
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
## [1] -0.3118808 -0.5803309  0.1917947 -0.3971647  0.1679505 -0.7433723
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
## [1] -0.0323
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.938044
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
## t1*      4.5 0.01761762   0.8765767
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 2 3 4 5 6 
## 2 1 1 1 2 3
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
## [1] -0.0475
```

```r
se.boot
```

```
## [1] 0.9173632
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

