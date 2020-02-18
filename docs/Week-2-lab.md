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
## 0 1 4 6 7 9 
## 1 1 1 2 3 2
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
## [1] 0.0173
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
## [1] 2.782202
```

```r
UL.boot
```

```
## [1] 6.252398
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
##    [1] 3.2 4.4 4.1 3.7 3.8 4.5 4.7 4.1 4.6 3.6 4.1 5.5 4.9 4.6 4.8 4.7 4.6 4.7
##   [19] 5.5 2.6 4.4 4.2 5.8 5.2 3.6 5.6 4.3 5.9 5.5 4.3 4.0 6.7 5.1 3.3 4.1 2.8
##   [37] 3.4 5.6 4.1 4.6 3.2 5.5 5.3 5.3 5.7 4.0 4.1 4.3 4.0 4.5 3.6 5.1 2.8 4.1
##   [55] 4.0 3.8 2.5 3.1 5.6 3.4 3.6 5.5 5.4 2.6 3.3 4.9 5.7 3.7 5.3 4.5 3.0 5.4
##   [73] 4.1 4.5 3.1 5.0 3.8 4.7 2.1 3.9 4.8 5.0 3.1 5.2 6.1 4.4 2.1 3.7 3.8 4.3
##   [91] 4.8 3.7 5.6 6.5 4.9 5.2 5.2 2.9 3.9 4.5 3.5 3.3 4.6 4.2 3.3 5.1 4.6 4.5
##  [109] 2.8 6.1 4.0 4.8 3.6 2.7 3.8 5.3 4.5 3.7 3.8 4.4 4.1 4.6 3.3 4.4 6.1 6.9
##  [127] 4.6 4.8 5.1 4.5 4.8 5.1 3.9 2.5 3.8 2.7 3.9 5.0 3.4 3.9 4.3 5.5 5.7 4.3
##  [145] 4.9 5.0 4.5 4.2 3.7 3.8 4.4 4.9 3.0 5.6 5.1 6.0 2.8 4.8 2.7 3.6 6.2 4.6
##  [163] 4.6 5.7 4.7 4.1 6.0 4.4 3.4 4.6 5.0 5.0 5.3 4.2 4.2 3.4 5.6 5.2 5.3 3.9
##  [181] 4.4 3.6 5.0 4.6 3.8 3.7 4.2 3.6 5.2 2.8 4.4 4.6 2.7 3.9 4.7 2.3 5.1 5.1
##  [199] 5.1 6.0 4.4 4.7 3.9 5.0 5.5 5.2 4.9 5.6 4.9 6.2 4.4 3.9 6.3 4.1 5.6 5.3
##  [217] 3.8 4.4 4.1 5.5 3.4 5.2 4.7 5.5 4.6 5.1 5.0 6.8 5.5 5.1 4.0 4.2 3.9 4.1
##  [235] 5.7 3.6 4.6 3.2 4.6 4.6 3.6 3.8 4.5 4.9 5.0 3.7 4.1 3.1 4.4 5.7 4.6 4.4
##  [253] 3.5 4.9 5.7 4.1 4.9 3.6 3.6 5.0 4.1 4.9 4.8 5.6 5.8 4.1 5.1 4.1 4.1 3.6
##  [271] 3.4 5.6 4.6 4.3 2.6 4.5 5.9 4.3 2.3 5.3 5.5 4.8 6.4 3.9 3.6 5.6 3.6 3.5
##  [289] 3.8 5.8 5.0 4.4 4.5 5.4 3.9 4.5 6.3 4.7 4.6 4.8 3.6 4.1 3.6 3.0 5.4 3.4
##  [307] 4.0 4.3 4.8 6.3 3.8 4.3 4.4 4.3 5.5 4.6 3.6 5.0 4.5 4.7 5.7 2.6 5.4 5.0
##  [325] 4.5 3.5 3.9 4.1 4.6 5.0 4.2 5.7 4.9 5.6 2.4 5.1 3.7 4.1 6.7 2.9 5.7 3.2
##  [343] 3.2 2.9 5.8 3.6 4.9 3.6 3.5 4.7 4.9 3.6 2.4 5.3 4.3 3.7 4.2 5.0 2.1 3.7
##  [361] 5.3 2.9 3.7 3.4 5.3 3.9 3.0 4.2 3.7 5.5 4.5 4.3 5.9 5.0 6.2 6.3 4.7 4.6
##  [379] 3.5 4.9 3.3 5.2 4.8 3.4 6.1 5.1 4.8 5.2 4.5 3.5 6.1 3.4 4.0 5.9 4.9 4.6
##  [397] 6.1 4.5 4.1 5.2 5.5 5.8 4.5 5.8 4.8 6.9 5.1 3.4 4.7 4.0 4.4 4.9 3.2 5.0
##  [415] 3.9 5.4 3.4 4.7 5.8 5.4 5.4 3.9 5.9 3.2 5.6 5.4 4.0 5.5 5.1 4.1 3.5 3.7
##  [433] 4.6 3.5 4.5 3.7 5.5 5.1 3.9 3.5 4.2 4.4 4.3 5.5 2.2 6.4 5.0 4.8 3.8 5.4
##  [451] 3.8 5.0 5.3 4.1 2.7 6.6 5.5 4.9 4.5 4.4 5.9 5.2 3.7 5.7 4.8 4.5 3.9 2.3
##  [469] 3.1 4.1 4.1 4.5 4.2 2.8 4.1 4.7 5.4 5.6 4.6 3.8 4.1 5.8 2.4 4.4 5.4 3.6
##  [487] 4.5 4.7 3.4 4.4 5.5 5.3 3.6 3.1 5.4 5.8 3.4 5.1 5.6 4.0 2.8 3.8 5.8 4.2
##  [505] 5.0 3.7 3.2 4.8 5.2 4.1 4.5 4.9 5.1 4.2 5.6 5.5 3.1 3.3 5.7 4.4 5.0 3.9
##  [523] 3.6 4.1 4.5 5.7 6.1 2.9 3.9 5.9 5.2 3.3 3.7 4.4 5.5 3.6 4.2 4.7 3.1 3.5
##  [541] 5.1 4.8 5.4 6.1 4.4 5.0 4.9 3.7 4.3 6.3 4.1 4.9 4.0 3.3 3.7 4.2 5.9 4.0
##  [559] 3.0 5.6 4.4 4.7 3.7 5.6 4.9 5.6 5.0 3.8 3.2 4.7 3.3 2.8 3.9 4.9 3.9 3.2
##  [577] 5.0 4.7 3.3 6.6 5.6 3.7 5.4 4.4 4.5 5.7 7.0 4.9 3.9 4.2 7.2 5.8 4.7 4.6
##  [595] 4.2 5.5 5.6 3.7 3.8 3.9 3.4 5.0 3.4 3.2 4.8 3.8 3.2 5.5 4.7 5.0 3.5 4.8
##  [613] 4.7 3.6 5.3 4.5 5.7 3.6 5.3 5.8 3.4 3.8 6.1 3.8 5.4 4.6 5.2 4.9 4.8 4.8
##  [631] 4.8 3.1 4.0 5.4 5.2 5.1 4.2 4.6 4.4 4.4 2.8 4.8 4.2 5.2 3.1 5.3 3.5 3.5
##  [649] 3.2 3.8 3.6 6.2 3.6 4.3 5.5 5.1 5.2 4.5 5.0 6.1 3.5 5.2 2.2 3.7 2.2 4.2
##  [667] 6.0 3.5 5.7 5.5 4.1 4.3 4.9 3.0 5.8 4.4 4.0 4.4 5.7 5.3 2.1 4.5 4.2 5.0
##  [685] 5.1 5.1 5.1 6.1 5.0 3.3 2.3 4.7 3.9 6.5 4.2 4.3 6.7 4.4 4.8 5.0 5.3 3.8
##  [703] 3.6 5.0 3.7 4.8 4.8 4.2 4.7 4.4 5.0 5.3 3.4 3.5 3.6 5.4 4.2 4.9 5.8 5.1
##  [721] 5.4 5.1 4.7 4.9 4.7 4.2 5.7 5.7 4.3 3.3 6.2 4.1 3.4 5.5 3.9 4.4 5.0 4.0
##  [739] 2.5 3.9 4.5 4.1 5.6 4.0 5.0 5.8 5.2 6.3 4.0 2.4 5.6 3.6 4.2 5.4 4.5 4.3
##  [757] 4.4 4.4 4.4 5.3 4.8 5.0 4.2 4.9 4.3 4.0 4.0 4.5 5.2 2.3 3.4 6.0 3.4 4.1
##  [775] 2.4 3.5 4.6 7.2 5.5 5.8 5.0 3.8 3.5 4.8 5.1 4.7 3.4 5.7 5.3 6.5 6.0 5.5
##  [793] 4.6 4.8 6.1 5.0 4.4 5.0 3.8 3.7 3.9 4.0 3.3 4.5 5.3 5.7 2.3 4.7 3.6 4.2
##  [811] 4.1 4.5 5.2 3.8 4.4 4.4 3.6 4.8 4.1 5.3 4.3 4.1 5.2 4.7 5.7 5.3 2.7 5.1
##  [829] 3.3 5.5 4.6 4.2 4.6 3.7 3.5 5.0 3.6 3.9 4.3 3.8 3.5 4.7 5.9 5.5 4.6 3.4
##  [847] 4.5 4.6 4.9 5.4 3.6 5.6 3.8 3.4 5.7 3.1 4.6 5.9 3.5 6.2 4.3 4.3 4.9 5.0
##  [865] 5.3 4.3 6.0 5.1 4.8 5.5 3.8 4.4 5.4 3.6 4.5 4.2 3.9 5.1 4.1 6.6 3.2 4.8
##  [883] 3.5 4.8 4.8 5.4 5.0 3.7 4.3 5.0 4.2 5.2 2.7 3.9 4.7 5.7 5.6 4.0 3.0 2.2
##  [901] 6.3 5.4 5.6 4.2 5.3 4.8 4.3 4.0 4.1 4.4 4.2 5.1 5.2 5.7 5.5 5.1 4.7 6.1
##  [919] 4.7 4.8 5.7 6.5 6.0 4.3 5.2 4.9 6.1 5.3 5.9 2.7 4.0 6.2 4.7 4.4 4.6 4.8
##  [937] 4.1 5.1 4.3 3.7 3.2 3.5 4.0 4.5 5.3 3.8 3.4 3.5 5.2 4.4 3.5 4.2 4.7 3.1
##  [955] 3.2 3.3 3.7 4.2 3.8 4.7 5.3 5.8 3.4 3.5 6.4 3.5 5.8 4.9 5.1 3.2 4.7 4.5
##  [973] 4.8 3.5 6.1 2.8 5.2 5.3 4.6 4.7 4.0 5.4 4.4 4.9 4.8 3.6 5.2 6.6 3.4 4.6
##  [991] 3.9 4.9 4.0 4.1 5.2 3.6 5.2 5.0 3.2 4.8
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
##    [1] 3.7 3.6 4.8 5.6 4.8 3.7 3.9 4.0 3.4 4.5 3.3 5.1 3.8 4.8 4.9 4.7 3.6 4.2
##   [19] 6.0 4.2 4.1 4.4 4.5 4.3 5.2 5.0 3.1 3.3 3.3 4.2 3.2 5.3 5.4 4.7 2.7 5.1
##   [37] 3.4 2.8 5.0 3.9 4.6 4.3 5.0 5.6 4.9 4.1 3.9 4.1 4.8 4.1 7.0 3.9 4.4 4.2
##   [55] 4.2 2.9 5.2 4.9 5.3 4.0 6.6 3.0 6.1 5.2 4.7 5.2 5.7 3.1 4.6 3.0 3.5 4.7
##   [73] 5.1 4.1 5.5 4.0 5.1 5.8 4.0 4.3 5.9 4.4 4.4 3.5 5.0 5.8 4.4 3.3 3.6 4.5
##   [91] 2.3 5.6 3.4 3.5 3.7 4.2 6.6 4.5 4.5 2.8 4.7 4.7 5.2 4.1 3.7 3.5 5.0 5.3
##  [109] 5.8 4.9 5.1 5.0 4.0 4.6 4.3 4.7 5.4 4.9 3.8 3.4 3.7 3.5 6.3 4.6 3.1 3.4
##  [127] 4.5 5.6 5.8 4.6 2.6 4.5 6.3 6.0 4.7 4.1 5.5 5.2 5.7 4.8 4.1 3.0 5.8 4.4
##  [145] 3.4 4.7 4.0 5.5 2.9 3.6 5.2 6.1 3.2 6.5 3.7 3.9 5.1 5.2 6.3 4.4 3.4 3.8
##  [163] 3.3 3.7 4.9 3.4 4.9 5.0 4.6 4.1 4.0 5.9 4.0 3.9 5.6 5.5 3.4 5.0 4.4 4.9
##  [181] 4.3 5.4 2.8 4.3 2.9 5.6 4.6 4.1 3.9 5.4 3.4 4.9 5.7 3.7 6.5 3.9 3.3 3.8
##  [199] 4.6 3.3 5.1 3.7 4.5 4.8 5.6 4.7 3.4 5.1 4.6 4.4 4.0 7.2 6.1 3.5 5.0 3.2
##  [217] 4.4 5.4 4.6 4.4 5.4 3.7 4.2 2.8 3.2 4.0 3.7 3.5 3.5 4.2 5.6 5.0 2.7 3.8
##  [235] 5.0 4.5 3.9 3.5 4.4 2.2 2.8 6.3 6.1 4.6 5.0 4.7 5.5 3.9 4.3 4.8 5.4 4.8
##  [253] 2.1 3.8 5.4 4.3 4.2 4.5 5.8 4.1 3.7 4.2 4.3 2.3 4.5 4.5 4.2 3.0 4.9 5.1
##  [271] 5.7 3.8 4.2 3.5 4.5 3.9 4.5 5.6 2.6 5.6 2.6 4.8 4.7 4.7 3.8 5.2 5.1 4.7
##  [289] 4.2 6.6 4.7 4.9 4.3 3.8 4.5 4.3 5.3 5.3 3.9 6.0 5.6 4.9 6.2 2.1 4.6 5.4
##  [307] 4.5 4.7 4.9 2.9 6.3 5.0 5.2 4.7 3.3 5.6 5.1 4.0 4.2 4.2 4.6 4.4 4.9 4.0
##  [325] 3.1 6.0 3.7 3.2 5.9 3.0 3.5 4.0 4.9 6.5 3.9 5.0 3.2 5.7 3.6 4.9 3.8 3.7
##  [343] 5.5 1.4 4.5 6.4 3.8 5.0 4.5 6.1 4.0 4.6 5.9 4.8 3.8 5.4 4.1 6.3 3.5 4.6
##  [361] 3.5 3.6 3.4 4.9 4.6 5.0 3.3 5.9 4.5 6.4 3.9 3.2 4.2 5.8 5.1 3.7 5.8 6.9
##  [379] 4.6 5.5 4.7 5.1 3.0 4.1 5.2 4.4 4.8 4.8 3.7 4.2 5.7 5.0 4.3 4.2 2.6 2.7
##  [397] 5.7 3.2 3.4 5.5 4.7 4.8 4.9 4.2 4.5 4.9 5.3 4.1 3.7 4.7 3.7 4.0 5.7 5.0
##  [415] 3.3 6.1 3.9 3.6 5.9 3.4 3.7 4.4 5.2 4.0 4.2 4.8 5.7 2.3 4.1 4.5 3.3 5.0
##  [433] 4.4 5.0 3.7 5.7 5.5 5.1 4.7 5.1 3.1 3.5 5.4 5.7 4.0 4.4 5.3 5.1 4.0 4.5
##  [451] 5.5 4.6 5.5 5.1 4.7 5.2 3.3 2.3 2.8 4.5 6.3 4.8 6.0 5.5 4.6 6.2 5.3 3.8
##  [469] 4.7 3.0 4.7 5.3 4.8 4.7 4.0 6.0 2.5 4.0 4.0 4.6 3.3 3.7 4.2 4.8 4.9 4.4
##  [487] 5.5 4.7 4.4 5.1 5.6 1.6 6.6 3.0 2.6 5.4 4.2 3.1 3.4 3.6 3.7 2.3 2.9 4.4
##  [505] 3.3 5.3 4.0 4.9 4.6 5.3 3.7 5.4 3.1 4.7 4.2 5.9 3.5 4.5 4.6 5.1 4.9 2.5
##  [523] 5.7 3.1 3.8 5.0 4.0 6.7 2.6 3.1 5.2 4.0 4.3 5.8 5.0 4.7 5.3 4.7 5.9 4.2
##  [541] 3.8 4.2 3.7 4.8 3.4 6.7 5.3 5.0 5.1 4.3 3.6 3.8 3.8 5.5 5.8 4.4 5.4 2.8
##  [559] 4.6 4.8 5.1 4.1 3.9 5.9 5.5 3.4 4.2 5.2 5.4 4.7 2.2 6.4 5.4 5.4 3.7 5.2
##  [577] 4.8 4.3 5.4 4.5 5.3 4.2 5.8 3.8 4.1 4.0 5.7 3.4 6.3 3.1 4.3 4.3 4.0 5.7
##  [595] 5.1 5.2 4.0 5.1 4.0 4.7 5.4 4.6 5.6 3.2 2.1 5.9 5.2 4.6 5.1 5.6 4.8 2.6
##  [613] 3.8 4.6 5.5 3.7 3.9 3.2 4.7 4.1 4.0 3.9 4.5 5.7 5.6 4.7 5.4 3.3 4.3 5.7
##  [631] 4.7 4.8 5.7 4.6 4.4 4.9 5.8 5.5 3.8 4.1 3.8 5.9 4.9 5.2 4.6 4.7 4.5 4.7
##  [649] 5.3 3.0 3.8 4.3 5.3 4.6 5.2 4.9 5.2 5.1 4.1 3.3 3.4 3.7 3.9 3.9 3.2 4.4
##  [667] 6.5 5.0 3.5 6.4 4.4 4.8 2.8 4.1 2.7 5.4 3.9 6.4 4.1 4.6 3.5 5.1 4.0 4.1
##  [685] 4.4 4.7 4.9 4.5 5.6 4.5 4.8 2.8 3.5 4.3 4.6 5.6 5.1 4.9 4.8 3.1 3.1 4.5
##  [703] 5.1 4.8 4.9 3.8 4.0 5.9 3.0 3.7 3.3 6.4 4.2 4.6 5.5 5.0 4.9 3.5 5.3 4.5
##  [721] 3.7 4.8 4.7 6.4 5.6 4.9 3.2 4.5 4.5 4.8 5.1 4.8 4.5 4.2 3.4 4.3 4.4 4.1
##  [739] 3.4 5.2 4.0 4.7 4.5 4.9 4.9 4.9 3.8 5.0 4.2 5.0 4.5 3.9 3.0 4.5 4.8 4.5
##  [757] 5.8 2.5 4.9 4.3 2.8 2.8 3.0 3.6 4.3 6.2 4.7 5.6 5.3 6.8 3.6 4.2 5.8 6.3
##  [775] 3.7 5.0 5.2 4.6 4.5 4.0 5.4 4.7 2.5 4.4 3.8 3.8 4.3 4.8 3.6 5.4 4.4 3.8
##  [793] 5.2 4.1 5.2 3.7 4.1 4.8 5.5 4.4 4.1 6.1 6.3 4.6 4.3 3.8 6.5 4.2 5.8 4.8
##  [811] 4.7 4.1 5.0 3.7 5.7 5.2 3.3 3.9 4.7 5.6 2.7 4.2 5.3 4.7 3.6 3.9 4.9 4.3
##  [829] 5.1 6.4 4.0 5.6 4.1 4.3 5.8 4.4 3.6 4.9 5.0 3.9 5.5 3.7 5.1 4.9 4.3 3.4
##  [847] 6.0 4.1 4.1 3.6 4.5 5.5 5.6 3.6 2.2 5.5 5.0 4.3 3.5 5.5 4.7 4.0 4.5 2.9
##  [865] 3.1 5.5 4.1 4.2 5.1 5.3 3.2 4.0 4.4 5.0 4.8 3.2 4.8 4.6 3.9 5.3 4.8 3.3
##  [883] 5.1 4.4 3.6 5.2 4.2 5.6 4.7 5.3 3.2 4.0 5.5 2.9 4.6 3.6 5.2 4.4 5.5 6.3
##  [901] 3.5 3.9 4.2 4.2 4.4 3.8 4.1 6.0 6.8 3.6 6.2 4.8 4.7 4.4 3.8 4.4 5.9 2.9
##  [919] 5.1 3.8 3.9 5.4 3.7 4.8 5.1 4.0 3.7 5.3 2.8 4.3 4.8 5.0 5.7 5.8 4.8 3.6
##  [937] 5.2 5.8 3.7 5.7 4.4 2.7 4.5 4.5 4.4 4.2 4.5 6.4 4.1 2.8 4.3 6.1 4.9 4.3
##  [955] 4.0 4.0 5.0 3.5 3.5 4.5 4.2 4.3 4.1 3.8 4.2 4.4 3.9 4.3 4.8 2.6 4.3 3.3
##  [973] 5.9 4.7 3.7 3.8 4.2 5.2 4.8 4.6 5.7 4.0 3.8 3.7 4.3 4.2 5.8 3.8 4.4 5.1
##  [991] 3.9 4.1 5.6 5.4 4.6 2.8 4.5 5.9 5.8 3.6
## 
## $func.thetastar
## [1] -0.0082
## 
## $jack.boot.val
##  [1]  0.53903134  0.44285714  0.27604790  0.14141689  0.05536723 -0.10263930
##  [7] -0.12972222 -0.28337802 -0.41084011 -0.51910448
## 
## $jack.boot.se
## [1] 1.009294
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
##    [1] 3.0 4.9 4.5 6.5 3.1 7.1 4.3 2.8 5.1 5.3 4.2 3.8 3.6 4.0 4.0 4.3 4.7 2.3
##   [19] 2.8 3.4 5.1 4.3 4.2 4.6 3.9 3.9 4.7 4.6 4.7 4.1 4.2 3.4 6.7 4.0 3.4 5.4
##   [37] 5.9 5.2 2.9 4.2 4.5 4.2 4.0 4.6 4.8 4.6 5.5 3.7 6.1 4.8 2.9 5.2 3.2 3.5
##   [55] 4.6 3.7 1.7 3.6 3.8 3.7 5.6 5.2 3.1 5.4 4.7 5.3 4.5 3.6 4.1 4.1 4.5 4.6
##   [73] 2.6 3.9 4.0 4.9 4.2 4.2 5.2 3.4 3.3 3.5 3.3 3.0 3.9 4.6 3.5 4.3 3.8 5.2
##   [91] 3.6 5.0 5.3 2.4 5.3 5.1 4.6 5.8 4.9 4.7 5.0 5.0 5.5 4.6 4.2 2.7 4.6 3.1
##  [109] 4.0 5.2 3.8 5.0 4.6 5.3 4.6 7.0 4.9 5.5 5.6 4.1 4.5 4.1 4.6 3.8 4.6 5.0
##  [127] 4.0 4.1 3.8 3.9 3.4 3.6 5.1 4.2 6.2 3.0 6.0 4.6 5.2 3.5 4.4 4.4 3.0 3.6
##  [145] 4.6 4.4 4.6 3.6 3.7 5.0 6.0 4.6 2.2 5.0 5.5 4.8 5.5 6.2 4.8 3.8 4.2 4.7
##  [163] 3.8 5.0 4.2 5.0 4.2 5.1 4.5 4.1 3.1 4.5 3.9 3.7 5.5 5.3 2.9 4.8 4.5 6.4
##  [181] 5.1 3.6 3.3 5.0 5.3 3.9 5.7 6.3 3.8 5.0 4.4 4.4 4.8 4.2 4.9 5.1 4.3 3.2
##  [199] 4.2 3.7 3.8 5.8 3.8 5.8 4.0 3.5 3.6 3.8 4.5 4.3 6.0 2.4 5.5 5.8 5.5 4.7
##  [217] 5.4 3.0 6.2 6.4 6.3 4.3 4.8 4.4 3.7 4.3 6.3 5.0 5.0 3.4 4.0 2.4 3.1 5.3
##  [235] 4.8 4.6 5.2 4.7 3.3 4.0 3.7 5.7 4.8 5.3 4.1 4.2 4.3 4.5 5.5 3.5 5.0 4.4
##  [253] 4.7 4.5 4.3 4.3 5.6 4.3 4.1 5.1 4.9 4.3 4.3 3.8 4.9 4.5 5.3 4.0 4.3 5.1
##  [271] 5.9 7.2 5.1 3.2 5.3 4.5 4.2 3.7 4.6 3.4 3.6 3.9 7.0 3.5 5.8 4.6 5.1 5.0
##  [289] 4.2 3.4 6.2 3.3 6.0 5.1 4.5 3.8 5.3 5.8 2.2 4.5 4.6 4.1 4.6 4.8 5.1 5.1
##  [307] 4.9 4.5 3.9 4.7 4.3 3.9 5.2 6.2 5.5 3.1 4.3 4.5 4.1 3.6 4.7 5.1 3.7 4.2
##  [325] 4.5 3.9 3.5 3.3 4.9 6.0 3.0 5.4 4.8 2.9 5.6 5.4 3.9 4.1 4.4 2.5 6.8 5.2
##  [343] 4.0 6.0 3.9 2.0 4.7 4.0 3.2 3.0 4.7 5.6 4.6 4.6 4.5 5.4 4.7 3.7 3.7 5.2
##  [361] 5.2 4.1 3.8 3.7 4.7 3.3 4.7 3.3 5.7 3.2 3.4 3.2 5.5 4.2 3.6 3.9 4.6 4.1
##  [379] 4.6 5.1 4.7 4.0 4.4 5.0 4.2 4.4 4.9 4.0 3.6 4.1 4.1 3.3 4.8 3.5 3.7 3.6
##  [397] 3.6 4.4 3.6 4.5 3.3 3.8 3.5 4.8 3.7 3.8 6.5 5.0 5.2 5.4 5.3 5.3 4.2 4.2
##  [415] 5.1 4.8 4.4 5.1 3.8 4.0 4.5 3.7 4.8 5.4 4.8 5.1 3.3 4.8 5.2 5.1 5.2 3.9
##  [433] 5.6 4.7 3.9 3.6 4.7 5.4 5.5 3.7 5.9 4.6 2.8 4.7 4.6 3.7 3.7 6.3 3.5 3.2
##  [451] 3.8 4.9 4.6 3.2 3.8 5.3 4.4 6.6 5.0 5.6 3.3 6.0 5.5 3.9 4.2 3.3 4.1 5.0
##  [469] 5.4 4.9 4.2 3.6 6.5 2.6 4.9 4.1 4.6 5.3 5.9 5.3 4.3 4.1 3.9 2.8 5.6 6.1
##  [487] 6.1 5.1 4.4 6.9 4.8 3.9 4.6 4.4 3.4 4.4 2.7 5.5 5.0 3.3 3.5 5.8 5.7 3.3
##  [505] 3.3 5.1 4.0 5.7 4.2 4.7 4.5 3.9 5.1 5.7 4.3 4.9 5.6 3.7 4.0 4.0 4.2 4.9
##  [523] 6.2 3.8 5.0 6.5 5.8 4.4 3.3 4.0 4.9 5.0 2.8 5.0 5.0 3.9 4.7 4.6 3.7 5.5
##  [541] 2.6 5.0 3.5 4.9 4.3 3.6 4.7 4.6 4.1 3.7 5.7 3.6 5.1 4.3 4.7 5.2 4.5 5.6
##  [559] 5.6 3.4 3.6 4.1 3.7 3.7 3.9 4.2 4.9 4.6 4.1 6.7 5.1 3.6 3.9 4.9 4.2 4.2
##  [577] 4.9 4.2 4.2 5.7 5.4 5.8 7.1 5.9 3.5 4.6 4.0 4.5 4.0 6.4 6.0 4.8 5.9 3.9
##  [595] 4.1 4.5 4.4 4.3 6.8 3.1 6.0 4.8 4.9 5.0 5.5 4.6 3.2 3.5 5.4 4.7 5.3 4.7
##  [613] 3.7 3.3 4.5 5.1 4.8 5.2 3.4 3.9 2.6 4.8 4.7 6.1 4.3 4.8 5.8 4.4 2.3 4.6
##  [631] 4.2 4.1 4.0 4.0 3.5 3.7 3.7 4.6 5.1 5.0 3.4 4.2 3.6 5.1 5.6 3.0 3.8 4.2
##  [649] 4.1 5.5 4.0 3.6 4.0 4.7 5.2 5.1 4.5 3.4 4.7 2.5 4.8 3.1 4.4 5.9 4.3 4.3
##  [667] 5.0 4.0 4.1 4.6 3.2 2.1 4.3 2.6 4.3 4.2 2.8 4.8 4.7 4.3 4.1 4.9 4.1 4.8
##  [685] 4.1 5.2 2.5 3.8 4.6 3.7 4.0 4.9 3.7 3.0 3.9 5.2 4.1 4.5 4.2 2.8 5.6 4.0
##  [703] 3.2 2.3 5.1 4.7 5.0 3.6 4.1 4.3 4.2 2.4 6.0 3.1 4.9 4.2 5.0 6.4 3.4 3.8
##  [721] 3.8 4.6 3.4 4.7 5.9 4.5 4.1 5.3 3.7 3.2 4.8 4.4 3.2 5.0 4.0 5.0 4.7 3.9
##  [739] 4.7 4.9 6.2 4.8 3.1 4.6 4.6 6.7 5.3 4.9 4.3 4.3 4.5 4.3 5.0 4.1 4.0 3.8
##  [757] 4.9 4.7 6.8 4.5 4.7 5.9 4.7 4.6 2.8 5.3 2.6 4.1 4.2 4.3 5.1 4.2 5.2 5.9
##  [775] 5.2 5.5 3.5 5.3 3.4 4.8 5.0 5.4 4.2 4.2 3.2 4.6 5.9 3.6 4.8 5.9 5.5 4.5
##  [793] 4.8 5.4 3.6 4.9 3.5 3.9 4.6 5.0 6.8 3.7 3.4 4.4 4.1 5.5 3.6 5.4 5.3 5.0
##  [811] 6.7 3.5 5.5 4.7 3.6 3.2 6.3 3.8 4.9 3.6 5.5 4.8 5.1 6.0 5.6 4.9 3.4 4.5
##  [829] 2.1 3.4 6.4 4.8 3.8 3.3 4.0 4.5 4.3 5.6 4.1 3.3 4.4 3.1 4.1 4.1 3.4 3.9
##  [847] 4.3 4.1 3.8 4.8 3.5 3.7 5.3 5.3 2.8 4.3 4.0 4.6 4.0 4.8 3.6 3.1 4.5 4.2
##  [865] 5.7 4.1 4.3 5.3 4.4 3.3 4.5 3.6 5.4 4.9 3.2 5.0 5.1 3.2 4.3 4.2 6.4 5.4
##  [883] 3.8 5.4 4.7 3.7 3.9 1.8 6.3 5.4 4.0 4.8 3.0 4.0 4.7 4.6 5.3 4.6 5.1 5.0
##  [901] 3.3 5.7 5.2 4.9 5.2 5.4 4.7 5.4 5.3 4.9 3.6 4.8 5.7 5.0 2.3 5.2 2.7 5.1
##  [919] 4.2 4.9 6.4 3.6 3.5 5.7 4.8 5.0 4.3 3.8 5.0 6.2 5.2 5.1 4.6 4.1 3.0 6.4
##  [937] 5.0 3.8 4.4 3.7 6.6 5.4 4.6 3.2 4.0 5.7 3.8 4.2 6.2 5.1 5.4 4.9 5.9 2.0
##  [955] 5.0 5.0 4.0 4.9 5.2 4.2 4.3 5.3 4.8 3.5 5.3 4.7 4.1 4.9 3.4 6.0 2.6 5.5
##  [973] 4.2 3.9 3.7 3.0 5.2 4.6 3.4 2.8 6.3 4.2 3.8 4.0 4.3 5.6 4.5 5.9 4.8 4.5
##  [991] 4.3 5.2 3.7 5.3 5.4 4.2 5.0 5.6 3.3 4.4
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.4 5.4 5.3 5.2 5.1 4.9 4.8 4.7 4.5 4.5
## 
## $jack.boot.se
## [1] 0.9931767
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
## [1] 2.024552
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
##   3.732747   8.208215 
##  (1.600344) (3.766848)
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
## [1]  0.1974963  0.4901254  0.3871659 -0.1322566  0.6022018  0.1418273
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
##    [1]  1.7387409671  2.3859135404  0.4494748379  2.1361568147  2.0006049694
##    [6]  1.3612187377  2.2458421045 -0.8961232260  1.6327022613 -0.1422491924
##   [11]  1.3262795267  1.1102527519  0.9316345346  2.4192831076  0.7548684300
##   [16]  1.8236002938  0.8098757313  0.7285412891 -2.5614271040  1.9769877823
##   [21]  1.6916699827  2.2812995216  1.7639477221  0.5497217920  1.3726187996
##   [26]  1.2712721352 -0.1428164763  1.7986548932  1.1786302355  1.9373073903
##   [31]  1.3442601311  1.3171986800  1.6742146386  1.6935125125  2.0640255451
##   [36]  0.1846171986  0.7100346325  2.0141014869  2.3630805252  2.4792533891
##   [41]  0.5038550904 -0.0261514125  0.1072338932  0.0567481008  0.8715661159
##   [46]  1.8124439105  0.6704068774  2.1844845691  1.7292081663  0.4636523220
##   [51]  1.9123817005  2.0228457219  1.8139969987  2.2924740752  0.6707575968
##   [56]  1.4195014102  0.2436514354  1.1494617233  1.0659603837  0.5644942141
##   [61]  2.2693996788  1.9805549578  0.6964136024  2.3552073071  1.7218066240
##   [66]  2.4731892495  2.1361457806  1.3186855235 -0.1676583173 -1.9682992946
##   [71]  2.4268732291  2.0386138218  0.9540763389  1.3703318989 -0.1329702101
##   [76]  0.6859122588  1.3257587817  1.7936638877 -0.4774790350  1.9096664196
##   [81]  0.1602704353  0.4227240644  1.3330794134 -0.1668917193  1.1467764751
##   [86]  1.1024971792  0.2398714011 -0.2044117701  2.0172389078  2.4227109476
##   [91]  0.7349007573  0.6974928521  0.4284771762  1.3645857836  2.4277454984
##   [96]  2.5111184205 -0.8520304451  2.0982252496  0.5323428199  0.8287964645
##  [101]  1.3377761791  0.8250134823  1.3494092145  1.6166579949  0.6616106217
##  [106]  2.5800054878  2.4080748457  0.5011214993  1.0489903615  2.1558424007
##  [111]  1.6469627054  0.8286754509 -0.0924100801  2.2890676260  2.4397396432
##  [116]  2.2977705037  0.4424326856  0.7018413752  1.8443678918  0.8266529115
##  [121]  0.6461092007  1.6899223949  1.2186544580  2.0245523393  2.3163916842
##  [126]  1.4202869254  1.5565409677  2.1764631951  0.3881870641  1.9742619901
##  [131]  1.6046859016  0.1808018427  1.9956634113  0.6286342427  2.0245167448
##  [136]  2.1133799870  2.0152124731  1.3183696478  0.7247045093  1.7930639636
##  [141]  1.3693658884  0.5971410573  2.2833166443  1.4468518803  1.2435111769
##  [146]  0.8547204530  1.3385017472  1.8224500458  0.7724416538  2.4803067005
##  [151]  1.4442466037  1.1052799202  2.3383468742  1.3419777305  0.5393340098
##  [156]  1.9886776763  0.4147209468  1.7078901936  1.1119827959  1.7427805919
##  [161]  1.1651656840  0.6620116096  0.8708065850 -0.6680372048  0.7513140276
##  [166]  1.9315078211  1.3416186287  1.3168863981 -1.3122395511  2.4603173786
##  [171] -0.3308340841  0.8523836329  1.2044976271  0.7593143531  0.4036149855
##  [176]  0.9946844796  2.3448595042  1.4169379441 -0.1000665871  2.5557245467
##  [181] -0.0924100801 -0.8291694544  2.0178702770  0.8128539004  0.5495763302
##  [186]  2.3605580745  0.1292305085  1.6666276269  2.2704318596  1.6102354140
##  [191]  0.7089977120  0.5699063874 -0.9050686928 -0.4491948870  1.4194532892
##  [196]  1.0848304131  0.0468540848  1.3631634782  2.5473658901 -0.3470325649
##  [201]  0.7015161709  0.4208640636  1.2433636927  0.5437780935  1.2852826707
##  [206]  1.1873498106  1.2399672729  0.4100557953 -0.8879061735  1.0287695920
##  [211]  0.8105198989  2.1000718249  1.8992172821  0.7789104413  1.8438326415
##  [216]  2.1554972216  2.0238885281  0.4904473826  1.6154815476  1.8696830006
##  [221]  1.0052238000  0.6913161898 -0.8902039465  1.2326630679  0.6608531277
##  [226]  1.0255584377  1.4745863395  1.0162131936  1.2868925866  1.9993799342
##  [231]  1.6988181681  1.9313789953  0.7631762900  0.7593072322  2.4274208824
##  [236]  2.0702349454 -0.4633393140  0.7261920608  2.2748205138  0.1401347179
##  [241]  0.6866608794  2.1026098397  0.4944432106  0.0354011152  0.6941055634
##  [246]  2.2050039622  2.2310068601 -0.0541004065  1.5016150944  1.2145977167
##  [251]  0.0981827323  1.3522559726  0.9159703839 -0.9009812446  2.5587441521
##  [256]  2.1579546066 -0.2371676416  0.6052684344  0.2935534526  1.8143611993
##  [261]  2.0403284568 -0.4212308661  1.0268562231  0.5759354736 -0.4017020010
##  [266]  1.5873550769  1.3586048044  0.7095283067 -0.4019877012  0.3529592841
##  [271]  2.4346028193 -0.4056515545  1.9285951652  2.2138687323  0.5206566717
##  [276]  0.4741788162  2.3830557798  0.3491999717  2.3655528861  2.3583465585
##  [281]  1.9125994106  2.3269855291  1.2331421748 -0.2520816724  2.2779940540
##  [286]  1.0381532828  0.7986732483  0.8246573119  1.1136595013  0.7104919179
##  [291]  1.8879547671  0.4641117614  2.3590037321  2.2138504126  0.8585515609
##  [296] -0.1917393989  0.0631914233  2.3003297863  0.5758180264  1.3739938978
##  [301]  0.9584158060  0.8670381556  0.3021402360  1.7954586884  1.9529684908
##  [306]  1.9568528576  2.2451073348  0.7777476898  1.1665364399  0.6131851604
##  [311]  2.1388585298  0.7225748832  0.1876527805  0.8094523438  1.6917883539
##  [316]  1.2485381134  1.8227011248  1.0488577462  0.4391714513  2.1544072967
##  [321]  0.3945416544  0.2241206015  0.7070088542  1.9965336440  0.7148598909
##  [326]  2.2734026879  1.2333255377  2.6183738523  0.8036914014  2.0960936284
##  [331]  0.7103869652  0.2808247744  2.4436945754  1.4065308295  1.9141414616
##  [336] -1.1446157404 -0.6341481893  0.7245053820  2.4390402959 -0.6383401451
##  [341]  2.4349204158  0.5628896724  1.2265080187  1.8533680214 -0.8154756721
##  [346]  1.9676400208  2.0893330947  1.9576561474  2.1743082915  1.1739776678
##  [351]  1.2361813887  1.4018538776  0.9092340990  1.5398992149  2.0659429304
##  [356] -0.5282092893  0.8086881720  1.3231207512  2.1911872136  0.2070306748
##  [361]  1.9481122138  2.3693813788  0.3062523231  1.8944413919  1.3976073433
##  [366]  2.1652653340  1.1946433197  0.2654608083  0.8076208413  1.9545311289
##  [371]  1.2010198164  0.0937705624  1.3903124159  0.3437017730  1.0635042002
##  [376]  0.6194994559  1.2880011614  2.4570586836  0.4219278855  1.3396346926
##  [381]  1.9721591998  2.2085284408  1.1851735499  1.2038641714  2.1692420477
##  [386]  1.8839982389  1.0844658638  1.3323082144  0.9135148181  0.6080771067
##  [391]  1.4221131580 -0.3475149148  1.2707856230  2.5543959702  2.3616576293
##  [396]  1.0693214360  1.0754284827  2.1233152714  2.1371167583  2.5782619627
##  [401] -0.8979077662  1.0677618632  2.0555046940  0.5443539132  2.3112437739
##  [406]  1.2000375401  1.2765459963  1.9404998662  2.2127235383  1.2486762882
##  [411]  0.4014702643  0.3763239700  0.9923357904  0.3399339967  1.1140955093
##  [416]  1.2911691330  0.3885207467  2.4046203657 -0.0685746030  2.1643794804
##  [421]  0.0889975461  2.3646753130  1.2243362688  1.3959237123  1.7591874853
##  [426]  1.2794151472  2.0334973601  2.2635286313  1.2499287026  1.3900511574
##  [431]  1.4097607091  1.7407476850  1.3184001185  1.1025533032  1.8055890133
##  [436]  1.0290076161 -0.7071093811 -1.0588330371  0.6170997521 -0.9298250573
##  [441]  1.1552362238  2.2454219967  2.0004238410  1.2685807289 -0.2873425728
##  [446]  0.4264100744  0.6084199552  2.3426437027 -1.5964451569 -1.0561709137
##  [451]  1.3381459490  1.0167214004  0.3812075952  1.1119556322  0.6455818088
##  [456]  2.5463947993  1.4231132694  1.7026373607  0.7881730060  1.4172218360
##  [461]  1.0583759555  0.3760096359  1.3513109189  1.6927955475 -0.0198423420
##  [466]  1.2062338623  2.1923669483  2.2989691537  1.7707633070  0.4318871285
##  [471]  1.9027662264  1.2834691052  0.4562669792  0.9297044744  0.7888385796
##  [476]  1.2474314278 -0.4943058017 -0.1285518042  1.3977276666  0.8001005487
##  [481]  2.0742609865  1.0485711025 -0.9392231466  0.4582040989  1.4932410287
##  [486]  2.0806215252  2.0681579967  2.5526965081 -0.7191284381  0.2824002019
##  [491]  1.0471950146  2.3409757092  1.0461663594  1.9229160468  0.6855855664
##  [496]  0.1087764976  0.8133485168  2.2331027537  2.0064145319  2.1635355979
##  [501] -0.4943058017  1.1781579559  0.3826991480  1.1019169139  2.3093209115
##  [506]  1.2515728933  1.3836550601  0.4015024590  2.1946448164  2.2647495487
##  [511]  1.9589215379  1.2303698626  0.1359411709  1.9214914050 -0.7708830210
##  [516]  1.2998735946  1.8234329054  0.9321843305  0.8171911885  2.4971231784
##  [521]  0.8576357075  2.2958273820  2.4634838087 -0.5046291585  1.0518529410
##  [526]  1.8212364773  0.9383263131  2.0812535797  2.2744010751  1.1539668250
##  [531]  0.0005504124  2.4138213864  1.4143599790  1.2875562723  2.1374011211
##  [536]  0.7876387269  2.3403534548  2.0433591578  0.3182477187  2.3363155776
##  [541]  0.2483138756  2.2903804371  1.9990835124  0.3581619630  2.4230627434
##  [546] -0.3849507442  0.3551807392  0.7939976142  0.7946412083  2.2828207434
##  [551]  2.1903591615  0.5520186437  1.5534260390  1.6888439995  0.5114273478
##  [556] -1.1613786670  0.4219278855  1.6973996223  1.3920793934  0.3852421615
##  [561]  0.7353237220  1.9789490410  1.2317032751  2.0937161851  2.2124521988
##  [566]  1.3647149011  2.0316861427  0.1687671979 -0.2915581340  0.6946230505
##  [571]  1.2203259625  1.3095726933  1.3467692250  0.6962029227  0.4772260518
##  [576]  2.3179608039  1.7707633070 -0.7736028139  1.2669858947  2.3091292751
##  [581]  0.0031283476  2.3667202526  1.1313959060  1.2054336687  0.1491698518
##  [586] -0.5344569991  1.1108428955  0.5790544117  0.6329259901  2.0545610925
##  [591]  1.2559579516  1.8154356142  0.2596488294  1.3257779231  1.7316804373
##  [596]  2.0718580553  1.7037324143  2.3066910037 -0.1547212262 -0.3070911278
##  [601]  2.2734026879  1.1709320163  0.6897428031  1.6848729936  0.6878016412
##  [606]  2.0144180661  2.0438481098  1.4495920939  1.2314736148  1.8047841572
##  [611]  2.3282970335  1.2826016313  1.4606629713  2.2411114892  2.3898155470
##  [616]  2.4512558830  1.1845540928  1.7163374532  1.4458751512  1.7625319103
##  [621]  0.6043726484 -0.1234779735  1.5975691796  0.9918321528  1.2202359157
##  [626]  0.2074497560  1.6312523881 -0.3364424844  2.2783892835  0.3297507005
##  [631]  1.1910428196  1.7420116086  1.2638408452  2.1805864288  1.8514347226
##  [636]  2.3459578660  0.4927059663  2.5293263713  0.9315211408 -0.7527706049
##  [641]  0.5584172022  0.6603851394  0.9941430532  1.7656915152  1.3690072645
##  [646]  2.3788548643  1.8180811218  1.6375767264  2.4436667252  1.1787135038
##  [651]  1.2365668814  0.1892673465  0.9706154959 -0.3532865917  0.3373563532
##  [656]  0.6814558056  2.1228237481  1.6809033059  0.6870909888  1.9362946407
##  [661]  2.0980483598 -0.6520156910  0.6833318697  1.1956019218  1.2444464506
##  [666]  2.3359502884  2.1257247403  1.6711559905  0.7658482055  0.0149967508
##  [671]  1.3822664320  1.4206748478  1.4089080041  2.3297983908  2.2125700297
##  [676]  1.7839428006  0.8093683040  1.4553458454  2.1089219193  0.0711200635
##  [681]  0.9303898843  0.3048575368  1.3383636078  1.3637498647  0.7492745401
##  [686]  2.3677514025  2.0023199737  1.3589577476  0.8560070386  2.4579795349
##  [691]  0.5494829003  1.2096250419  0.4100017783  2.4900375243  1.1802378946
##  [696]  1.4141984272  1.3182072098 -0.0463623379 -0.1729495017  2.3195358995
##  [701]  0.7983013412  0.7688493019  2.0093957909  1.0967596985  0.5656082064
##  [706]  2.4525000666  1.2518599733  0.4414364671  0.9137781650  1.8140240406
##  [711]  1.3856518924 -0.5832940315  0.3377861115  1.9729761828  1.0820043500
##  [716]  1.2755783525  1.3378429472 -1.1525917135  1.1026381619  2.0393440499
##  [721]  1.5906097481  0.7346297414  0.1889856774 -0.0057881732  1.9064605153
##  [726]  1.8984061084  2.3058901829  1.4105310418  1.3554393869  2.3806463526
##  [731]  0.0238463088  1.7236794980  1.4406703742  1.4212463177  0.1962421465
##  [736]  2.2774708380  0.7823754642  1.2101182296  0.0021684961 -0.1873386480
##  [741]  1.1608022877  0.2754810529  2.4365928428  2.3168504725 -0.7107317045
##  [746]  1.1539668250  1.1224025363  1.3494092145  0.4494748379  0.0759038950
##  [751]  2.1565284037  2.1376712246  0.5404611876  1.2197360125  2.0420789580
##  [756]  0.7798004625  0.7622217822  2.2508580278  2.0910974433  0.1855864794
##  [761]  1.0854420068  2.2799289101  1.2343945424  0.7551975199  0.8383120602
##  [766]  0.3513508277  1.2103981020  0.6736065607 -0.3310458785  0.8282232536
##  [771]  0.8060111422  0.5541591997  2.0082368008  2.4130392873  1.9048429997
##  [776]  1.2054187480  0.3919382142  0.7084872429  2.3443836587  1.0658759478
##  [781]  1.2215932993  1.7669849811  2.3087126527  2.2847435752  1.7777304693
##  [786] -1.4097206802  2.0277318744  2.1722585436  2.3075520468  1.4176605921
##  [791]  1.6928402117  1.0703525025 -2.0254083867  0.8781411256  1.3314549680
##  [796]  0.8451771088  2.4070343925 -0.9566609308  1.3630706133  1.9898126465
##  [801]  1.3913408487  1.2188525004  0.7033441461  0.5493533270  1.7943749780
##  [806] -0.1899716652  2.0396249411  1.6406885308  2.1482918525  1.3385677796
##  [811]  0.0333526660  2.4218594321  1.2585147308  0.2860571749  2.5851618558
##  [816]  2.0100255841  2.3581873008 -0.0975168946  1.2094934488 -0.6701497920
##  [821]  1.5841680423  1.6866874208  2.3688482502  0.8960253101  2.0422464172
##  [826]  0.9937040891  2.4501253035  2.3282970335  2.3308174461  1.3964184738
##  [831] -0.4927271677 -0.2397199154  2.3073587748  1.6846166299  1.0890219565
##  [836]  1.2768307808  2.4146663165  2.0228144129  1.3467280507 -0.2334812733
##  [841]  1.9421261291  0.4264100744  1.9478682214  0.3103061995  1.3436952620
##  [846]  1.3077498731 -0.2529168394 -0.4070984701  0.8551753518  0.4313872986
##  [851]  1.1459502243  0.2782361008  0.1309542626 -1.2066734095  1.3582631113
##  [856]  2.0759614588  1.2213969583  1.8337642109  2.0828287960  1.8065821350
##  [861]  1.4413905608  1.9630164112  2.4196067183  2.0615177795  0.1265298354
##  [866]  1.2964508483  0.7808204398  2.3465112330  1.6879493493  1.3236253324
##  [871] -0.4448704962  2.3985997531  0.1756921949  1.9790997954  0.4664514749
##  [876]  2.1588761374 -0.4629628951  0.9020889412  1.4905047801 -1.0013588223
##  [881]  1.4254189064  2.1110241805  2.1220394619 -0.0460652178  1.4985403037
##  [886]  1.9952482314  2.2986369402  1.6393142251  1.7622226943  0.0803693678
##  [891]  0.5761451572 -0.4751214800  1.4077437761  1.2797194346  1.1317940381
##  [896]  2.0294144805  1.1623125034  1.2361190142  2.2536503256  1.3382265585
##  [901]  2.0382278469  1.2781241392  0.0206593509  1.0597472610  0.5736710043
##  [906] -1.0115412151  1.2133257607  2.4805001348  0.2368586388  2.0117575635
##  [911]  1.2535528877  1.3689410997  0.7431371869  2.0039697386  1.2170481597
##  [916]  0.1431484536  1.5794459716  2.1550640905  0.8643918534  2.0331842841
##  [921]  1.2535108017  2.1711526356  1.9630164112  1.7726625377  0.4666153528
##  [926]  1.9344462058  1.0665746391  0.8677434564  2.3010432063  0.9510754287
##  [931]  2.4748603157  1.7625319103  0.9816116167  1.9776014673  1.7737605377
##  [936]  0.3788080354  2.3260947846  0.3757798646  0.4787405443  2.0852978188
##  [941]  2.4287206335  1.9981849302  0.7272892371  0.8701920133  0.4642766933
##  [946] -0.7291114986  0.2096951819 -1.1804434396  0.3639072265  2.1116558157
##  [951]  0.1223011769  1.8257698414  0.5744204449  1.9578063845  2.1523071550
##  [956] -0.4019190427  1.4056628242 -0.0265623616  2.2174906081  1.5205053446
##  [961]  1.9564806762 -0.1706452140  1.2486071425  2.1393950510 -0.4529580078
##  [966]  2.1869493486 -0.5528419587 -0.5637426725  2.1743752733  1.1708511306
##  [971] -0.3930958000  2.2303922542  2.0893741048  0.7059055756  1.7381383196
##  [976]  1.5179306037  0.7707425017  0.9853535732  2.2824997080  0.9475235527
##  [981] -0.2895788904  0.5770496257  2.1027851151  1.6063980906  1.2870514895
##  [986]  0.2563640126 -1.2144152820  2.3086549054  1.5247414944  0.4532836121
##  [991]  0.9743586126  1.2389685714  0.7120194511  1.9046700787 -0.5314574815
##  [996]  0.5706374129  1.4572699757  0.5757687793  1.9278597173  2.4865781070
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
##   0.45475471   0.28961195 
##  (0.09158334) (0.06475626)
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
## [1] -0.4683735 -0.5640362  0.4563706 -1.0371835 -0.5962961  0.2432481
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
## [1] 0.0472
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8904072
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
## t1*      4.5 0.004604605   0.8833796
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 1 2 4 6 9 
## 2 2 2 2 1 1
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
## [1] -0.0123
```

```r
se.boot
```

```
## [1] 0.9144807
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

