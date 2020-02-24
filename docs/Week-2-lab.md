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
## 2 4 5 8 9 
## 4 2 1 1 2
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
## [1] 0.0176
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
## [1] 2.711953
```

```r
UL.boot
```

```
## [1] 6.323247
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##   2.5%  97.5% 
## 2.7000 6.3025
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
##    [1] 3.9 3.9 4.7 4.2 3.6 3.8 4.8 2.7 5.7 5.5 4.8 4.4 5.0 3.2 3.6 3.7 4.4 4.1
##   [19] 3.5 3.2 3.0 6.6 4.7 4.2 5.8 3.1 5.1 5.2 5.2 4.7 3.0 3.0 4.8 5.4 4.5 4.5
##   [37] 5.6 4.7 4.0 2.8 3.4 2.9 4.5 5.3 5.2 5.3 3.1 3.9 4.8 5.9 3.5 4.1 4.8 5.1
##   [55] 5.9 4.0 3.7 4.8 5.5 4.7 4.8 5.1 5.4 5.5 5.4 3.7 3.3 5.8 4.0 5.1 5.3 4.6
##   [73] 3.8 3.4 4.4 3.5 3.5 5.4 3.7 4.6 5.4 3.1 5.3 4.8 4.5 1.6 4.2 4.1 1.6 4.1
##   [91] 4.7 4.7 5.2 6.8 4.4 2.9 5.6 5.2 5.3 4.1 4.8 2.7 4.2 3.6 3.7 4.7 4.3 5.4
##  [109] 5.5 6.7 5.5 5.6 4.2 5.7 5.2 4.6 4.0 7.1 5.5 4.3 5.2 5.1 2.6 5.9 5.2 3.3
##  [127] 2.7 5.3 5.0 3.3 3.7 3.2 4.7 5.1 3.0 5.4 5.5 3.1 4.1 6.0 5.5 3.8 5.4 4.5
##  [145] 5.0 5.7 5.1 4.0 3.3 6.6 4.5 6.9 3.8 4.3 4.1 4.0 3.8 3.2 3.5 4.0 5.2 2.6
##  [163] 4.5 4.0 2.6 4.9 4.8 3.7 5.6 3.8 6.4 4.5 2.7 4.7 4.0 5.3 5.1 5.5 3.3 4.8
##  [181] 5.2 4.1 4.7 4.5 5.7 4.5 5.3 6.6 3.5 4.4 3.8 5.5 6.7 4.6 4.2 5.5 2.6 4.2
##  [199] 5.9 3.8 4.3 5.2 3.4 6.4 4.1 3.4 4.2 5.5 4.9 7.7 4.9 5.2 3.4 4.3 2.6 4.4
##  [217] 5.7 3.6 5.5 4.0 3.9 4.0 4.8 4.2 3.1 4.7 2.4 5.3 2.9 5.3 4.9 3.5 5.7 5.2
##  [235] 3.5 3.9 4.2 3.6 4.6 3.7 2.7 5.7 5.5 5.3 3.5 4.6 3.3 3.3 5.3 5.7 3.7 3.5
##  [253] 4.2 5.4 5.9 3.2 5.2 4.3 2.8 5.2 2.9 5.2 4.1 4.5 4.2 4.6 5.8 4.6 4.0 3.8
##  [271] 3.2 5.1 5.2 4.7 4.8 4.4 4.3 5.2 6.2 5.3 3.5 4.6 4.4 4.6 4.0 5.0 6.4 5.5
##  [289] 3.3 4.5 4.0 5.6 4.9 4.8 5.3 4.6 3.7 3.8 3.9 4.7 5.8 5.2 4.6 3.5 4.2 4.3
##  [307] 5.0 4.7 6.6 5.4 4.6 4.5 3.3 5.7 4.1 4.8 6.0 5.7 3.7 3.6 5.2 6.2 5.0 5.6
##  [325] 3.2 3.6 4.5 3.7 4.4 3.0 4.4 4.6 5.5 4.9 3.8 3.0 4.9 5.5 4.1 4.1 4.4 4.9
##  [343] 3.2 5.1 4.9 5.4 3.2 5.5 2.6 5.2 2.9 4.1 5.1 6.2 5.7 4.2 5.1 5.0 5.2 4.6
##  [361] 3.1 5.7 5.1 3.3 4.7 4.0 5.1 5.5 4.8 4.8 3.8 4.7 3.2 4.3 4.8 5.3 4.5 4.4
##  [379] 5.1 3.6 5.1 3.9 4.6 5.4 5.4 5.6 4.8 4.0 5.0 4.7 4.4 6.1 3.7 4.8 4.1 3.7
##  [397] 2.6 4.1 5.0 3.0 5.5 4.5 5.7 6.5 4.9 6.6 3.4 4.1 3.9 4.9 5.2 5.0 5.8 4.0
##  [415] 3.2 4.2 4.9 3.7 4.5 4.5 4.4 5.7 4.8 4.5 3.5 6.1 5.4 5.9 6.0 3.2 4.8 3.9
##  [433] 4.6 5.3 5.5 5.3 4.9 4.7 4.6 2.7 3.6 4.5 4.1 6.5 5.6 5.3 4.4 4.7 5.2 5.6
##  [451] 4.5 4.1 5.3 6.2 3.2 4.9 4.8 4.9 4.1 4.6 4.1 3.9 4.5 4.2 4.5 2.5 5.7 4.9
##  [469] 3.5 3.8 4.7 5.2 5.4 5.2 4.4 4.5 5.1 4.9 5.5 4.0 2.0 4.0 4.1 3.9 3.9 3.3
##  [487] 5.7 4.3 3.8 4.9 5.3 4.3 5.7 5.2 3.2 5.3 3.7 3.2 4.4 3.4 4.7 3.8 4.2 4.2
##  [505] 4.4 4.6 5.5 5.3 4.3 4.9 5.4 3.2 4.7 4.4 6.6 4.4 4.4 4.1 3.3 5.5 5.5 3.9
##  [523] 5.5 4.9 4.9 4.5 5.7 4.8 5.6 4.9 3.5 5.8 5.2 4.9 5.2 4.6 4.2 5.6 4.5 5.4
##  [541] 3.9 5.6 5.2 4.5 3.5 3.0 3.5 5.3 3.3 5.4 4.9 5.9 7.0 5.2 4.5 4.7 5.6 2.6
##  [559] 3.7 6.9 4.8 5.9 4.8 3.5 4.0 5.3 3.8 6.0 5.6 5.0 4.5 4.3 5.5 4.7 4.7 3.1
##  [577] 5.2 5.7 5.5 4.5 2.8 3.4 3.9 4.8 5.1 4.2 3.3 4.3 3.8 4.2 4.8 4.7 3.8 4.5
##  [595] 3.9 3.6 3.4 4.8 5.2 4.5 4.0 5.5 5.2 5.2 3.8 3.2 3.2 5.1 3.2 4.2 4.4 4.3
##  [613] 4.8 4.3 3.6 5.7 5.1 4.9 5.2 4.0 4.4 4.4 2.8 5.1 6.5 4.7 5.4 5.3 4.0 4.1
##  [631] 5.1 4.8 3.9 4.8 5.7 4.3 3.7 5.6 5.3 5.4 3.9 5.9 4.2 4.7 4.9 4.0 2.4 5.2
##  [649] 4.1 5.3 3.4 6.1 4.8 5.5 3.6 5.8 5.2 4.0 2.8 4.7 4.9 4.6 3.1 6.7 4.9 4.6
##  [667] 6.1 5.4 5.2 5.4 4.4 3.6 5.0 4.7 3.0 4.5 4.2 4.6 4.6 2.7 5.1 5.3 5.1 4.7
##  [685] 2.5 5.0 4.7 3.4 5.2 4.5 4.3 6.2 3.8 5.5 4.4 3.1 4.3 4.8 3.1 5.5 5.4 3.6
##  [703] 7.0 3.9 4.2 4.2 3.9 4.3 5.1 2.3 5.2 4.0 5.0 3.8 5.2 3.8 5.0 3.3 3.1 3.9
##  [721] 2.4 4.2 4.8 4.5 4.5 5.4 5.5 4.8 4.7 4.4 5.6 3.6 3.5 5.9 5.1 5.3 3.7 4.9
##  [739] 3.1 4.3 4.4 4.6 4.2 4.8 4.3 3.8 3.4 3.9 3.8 5.9 5.1 4.9 3.2 3.8 4.8 2.8
##  [757] 4.9 3.5 4.3 5.2 5.0 4.3 4.7 4.2 4.4 6.1 4.5 4.2 4.3 5.2 3.5 5.1 4.3 3.7
##  [775] 2.7 4.0 4.8 3.2 4.2 4.0 3.4 3.4 4.9 5.0 2.5 5.3 4.1 3.3 3.9 3.3 2.7 4.0
##  [793] 4.7 4.9 3.9 4.9 3.2 4.0 4.7 3.9 4.5 4.6 4.5 3.9 3.5 5.2 3.9 5.8 4.4 5.3
##  [811] 4.2 6.6 4.3 4.0 3.5 2.2 3.9 5.5 4.2 5.8 4.3 6.2 4.0 5.3 5.4 4.0 3.6 2.8
##  [829] 3.8 5.8 4.2 4.2 3.6 3.4 5.5 3.9 4.4 5.2 3.0 4.9 3.9 3.9 4.7 4.7 4.8 5.4
##  [847] 4.6 3.2 4.9 7.7 5.1 5.2 3.6 4.0 5.4 3.7 4.7 5.7 4.3 5.8 4.7 3.8 4.8 4.6
##  [865] 6.3 3.7 5.5 5.2 4.9 5.4 4.2 4.7 5.3 4.0 4.2 2.8 4.1 5.7 5.9 3.5 5.0 4.2
##  [883] 3.6 5.3 2.4 3.7 3.8 4.1 4.5 4.5 4.6 4.9 4.9 6.6 3.2 6.5 4.2 3.5 4.2 3.8
##  [901] 3.6 4.1 4.6 4.1 4.1 5.9 3.8 5.9 4.4 4.4 4.3 4.0 4.3 5.5 5.9 5.0 4.7 3.7
##  [919] 3.9 4.1 3.1 3.2 4.4 4.1 5.8 5.0 4.5 2.2 4.4 3.7 4.3 5.1 4.2 3.9 3.9 2.9
##  [937] 4.5 4.6 4.7 4.4 4.6 3.6 5.0 3.7 3.6 4.4 2.7 5.8 4.6 5.2 5.0 5.8 3.9 3.8
##  [955] 5.6 3.9 3.9 3.2 4.0 4.6 4.1 4.4 4.6 3.9 5.4 4.9 4.4 4.8 3.0 4.4 4.5 6.0
##  [973] 3.3 4.8 3.9 5.4 6.0 5.4 3.7 3.8 4.0 3.3 5.1 4.4 3.4 4.8 2.2 5.4 3.3 4.9
##  [991] 4.4 4.6 4.7 5.8 3.9 6.2 4.2 3.7 5.6 5.1
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
##    [1] 4.2 6.2 5.2 5.0 5.1 5.3 4.2 3.7 3.1 4.9 4.0 2.9 4.3 4.6 4.8 2.8 5.5 5.0
##   [19] 5.1 3.3 5.1 2.7 5.8 3.3 5.9 3.7 3.0 5.8 4.8 5.6 4.6 5.7 2.5 4.7 5.3 3.5
##   [37] 5.5 5.2 3.1 4.8 3.4 5.5 3.3 3.7 3.8 6.7 3.1 2.9 5.8 4.8 6.6 2.7 4.4 4.4
##   [55] 6.2 4.1 3.4 5.5 4.3 4.9 3.4 5.6 3.6 5.2 4.8 5.0 5.3 5.5 3.1 5.5 5.4 5.2
##   [73] 4.2 4.9 3.8 4.1 6.1 5.3 6.0 3.8 3.6 5.9 4.0 4.6 4.7 4.7 4.5 3.1 6.6 5.9
##   [91] 3.0 3.8 3.3 4.5 4.3 3.7 3.9 4.3 3.8 5.2 4.7 2.3 4.0 5.4 5.9 5.8 4.2 4.1
##  [109] 5.5 5.6 4.6 3.1 4.9 3.6 4.9 3.6 3.0 4.1 3.8 2.7 4.2 5.3 4.3 4.3 3.8 4.9
##  [127] 4.8 5.3 2.2 4.8 4.9 3.7 5.6 4.3 5.3 4.5 2.7 4.1 3.5 4.9 5.4 5.3 4.2 5.0
##  [145] 5.2 5.0 6.3 4.9 5.3 4.3 4.0 5.4 4.5 3.8 4.4 4.3 3.8 4.8 4.9 3.6 3.6 5.2
##  [163] 5.0 4.7 5.0 4.7 5.7 4.4 4.2 3.0 3.9 4.3 4.2 6.4 6.3 5.4 5.7 4.7 3.1 5.2
##  [181] 5.5 3.0 4.6 2.5 4.2 4.6 4.2 3.7 3.3 5.5 5.0 4.8 5.4 4.5 5.9 3.9 5.4 3.9
##  [199] 3.3 3.4 2.5 3.0 3.2 3.8 4.5 3.3 5.5 3.2 4.6 5.6 6.3 3.1 3.8 5.5 4.8 3.9
##  [217] 5.5 3.1 5.9 4.5 3.8 5.4 5.1 3.7 3.7 5.9 4.5 4.3 4.0 2.8 4.4 4.2 5.4 5.3
##  [235] 5.0 4.6 4.3 4.6 3.9 2.7 3.9 5.4 5.4 4.0 4.8 3.2 4.4 4.7 3.7 5.2 2.5 3.7
##  [253] 3.7 5.5 4.4 2.3 4.2 3.9 4.8 4.4 4.1 6.6 3.6 4.1 4.3 5.3 4.8 2.8 1.7 3.4
##  [271] 3.0 4.2 6.2 4.5 4.6 3.7 3.8 4.8 5.0 4.9 3.3 3.9 3.5 5.7 3.5 3.3 4.3 5.9
##  [289] 3.6 4.4 5.8 5.3 5.4 4.5 4.7 3.2 4.9 3.5 4.4 6.8 2.6 4.9 3.8 3.0 5.5 4.8
##  [307] 4.9 3.4 3.1 5.1 3.6 3.6 3.1 3.6 4.3 5.1 4.8 5.0 4.4 5.4 5.8 3.2 5.9 2.4
##  [325] 4.9 5.0 4.9 4.4 3.2 4.6 5.5 6.2 3.7 5.3 5.6 3.5 4.3 5.0 3.8 5.3 2.2 5.4
##  [343] 4.3 5.5 5.8 4.7 4.7 4.7 4.3 4.4 5.7 4.1 5.7 3.7 4.6 4.7 4.6 5.0 3.5 5.5
##  [361] 3.9 5.1 4.1 4.4 5.1 5.0 5.3 6.1 4.1 4.3 3.3 5.2 5.0 4.3 5.2 3.2 4.2 3.2
##  [379] 6.1 4.6 5.6 5.4 3.6 3.5 5.1 4.0 5.0 3.5 4.2 5.5 4.6 3.0 4.5 3.8 4.6 4.2
##  [397] 6.4 4.0 4.1 5.8 4.4 4.8 5.3 3.8 4.0 4.9 3.5 4.0 4.5 5.5 4.0 5.7 4.5 5.5
##  [415] 5.8 4.0 6.4 2.8 5.1 4.6 4.6 3.7 5.3 3.1 4.0 5.3 4.0 5.3 5.2 4.6 4.8 3.3
##  [433] 5.2 4.3 4.6 3.3 3.8 4.4 4.9 3.2 5.1 4.5 4.2 3.6 5.5 4.4 5.7 4.2 4.0 4.5
##  [451] 3.5 5.3 4.7 3.9 3.2 5.4 4.5 5.9 4.7 3.1 4.9 4.8 3.7 4.9 5.7 2.9 3.4 5.6
##  [469] 4.6 4.5 3.9 4.1 3.4 5.2 4.1 4.2 3.4 3.0 3.3 5.1 4.3 4.4 5.8 6.0 4.2 4.2
##  [487] 4.3 5.9 5.8 6.6 5.4 6.1 2.8 4.3 4.8 6.5 5.6 5.0 5.2 5.3 4.6 5.6 5.0 4.2
##  [505] 5.5 4.6 4.7 5.1 5.6 4.8 5.4 2.6 5.8 4.2 4.1 2.6 4.6 4.7 4.1 4.8 4.8 5.9
##  [523] 3.5 5.6 2.9 3.3 3.8 4.4 5.3 3.5 5.1 4.7 5.3 5.2 3.8 4.0 3.5 3.9 3.1 5.7
##  [541] 4.7 4.7 4.5 3.9 4.4 4.1 6.0 3.1 3.8 2.6 3.1 4.7 3.4 4.3 5.1 4.5 4.6 3.4
##  [559] 5.1 5.2 4.3 3.5 3.7 4.6 4.4 6.3 4.8 5.6 4.6 5.2 4.6 4.9 5.9 5.5 3.5 3.3
##  [577] 4.4 3.3 5.1 3.9 5.4 3.7 4.6 4.5 4.6 4.4 4.6 3.6 3.5 5.5 5.2 4.5 4.3 5.1
##  [595] 4.7 5.6 4.1 6.5 4.1 5.0 5.9 3.6 4.2 4.7 3.0 4.8 5.3 4.0 4.5 4.6 2.4 6.0
##  [613] 4.0 4.5 3.8 5.3 3.5 3.5 3.9 4.0 4.6 4.4 5.5 5.1 4.0 5.4 3.8 5.8 4.2 4.2
##  [631] 4.8 5.2 6.6 4.2 3.9 2.7 5.4 3.4 4.7 4.5 5.2 4.8 3.4 3.9 4.0 4.5 4.2 4.3
##  [649] 5.3 4.3 4.4 4.1 3.6 5.4 5.7 5.0 3.0 3.3 3.3 4.6 6.1 4.0 5.2 3.7 4.4 3.3
##  [667] 3.9 4.8 6.0 4.3 4.1 4.2 4.4 4.8 3.4 5.7 4.4 4.2 5.3 5.2 5.0 5.0 3.3 3.8
##  [685] 6.1 4.4 3.9 6.2 2.5 4.0 5.1 5.4 4.8 5.3 6.1 5.4 5.0 4.6 3.3 5.3 6.1 3.4
##  [703] 4.2 3.1 3.3 5.4 4.1 4.5 4.7 5.2 4.8 4.4 4.4 4.8 3.4 4.0 6.7 6.0 4.1 3.6
##  [721] 3.9 4.8 3.8 3.4 4.9 4.9 3.2 4.6 5.8 4.6 3.4 5.0 4.0 5.5 3.4 5.0 4.9 5.2
##  [739] 5.4 3.5 5.1 3.1 5.5 3.8 4.9 5.6 5.9 4.8 4.6 4.8 2.2 3.9 4.4 3.2 5.5 5.6
##  [757] 5.8 4.4 6.1 6.1 5.5 3.7 4.7 3.2 3.0 5.9 5.1 5.6 5.2 7.1 4.2 5.2 4.9 3.9
##  [775] 5.1 2.8 3.1 4.2 4.2 3.7 4.9 3.1 6.5 3.8 6.8 3.9 4.1 4.7 4.6 3.9 3.2 3.8
##  [793] 3.1 4.7 4.2 2.3 4.7 4.5 2.8 6.1 3.8 4.8 4.6 5.3 2.6 5.5 4.5 5.2 4.2 4.6
##  [811] 5.4 4.5 4.2 4.7 5.1 3.1 6.3 4.7 5.8 5.5 4.1 3.9 3.3 4.6 3.5 4.0 5.0 5.2
##  [829] 4.5 4.6 2.7 5.5 5.8 3.6 4.6 3.8 3.0 3.0 4.2 3.5 5.2 5.6 5.5 5.5 5.6 4.6
##  [847] 3.7 5.5 6.3 4.0 4.8 5.2 4.2 5.5 4.9 5.9 3.9 5.9 5.0 4.8 4.8 6.9 5.0 3.9
##  [865] 4.7 3.4 3.8 4.4 5.0 4.2 3.7 3.3 5.0 5.8 4.9 6.1 4.4 6.1 4.8 4.4 4.4 5.3
##  [883] 2.8 5.3 4.1 4.1 6.0 2.6 4.9 5.9 3.6 5.0 6.8 4.8 3.8 4.4 3.3 4.3 5.8 4.8
##  [901] 5.1 4.6 3.9 4.2 5.1 4.0 4.9 3.7 4.4 5.2 4.5 4.9 3.5 3.2 3.5 3.8 3.7 5.3
##  [919] 3.4 6.1 3.8 4.7 5.0 3.9 6.1 2.5 5.2 6.0 4.4 3.7 5.0 3.2 4.4 4.1 3.1 3.6
##  [937] 3.2 3.6 5.8 4.5 4.7 4.6 4.7 3.2 4.3 4.9 4.7 3.1 5.4 3.9 3.5 5.4 4.4 4.4
##  [955] 4.6 4.7 4.9 4.7 5.2 4.5 4.0 5.0 4.1 5.2 4.4 5.0 6.6 4.1 5.5 5.1 3.7 3.7
##  [973] 4.2 3.4 4.4 6.0 6.7 3.9 4.5 6.1 3.4 3.0 3.9 5.0 6.5 5.2 5.1 4.3 3.4 4.6
##  [991] 4.6 3.1 3.6 3.5 4.2 5.1 5.3 2.8 4.0 4.4
## 
## $func.thetastar
## [1] 0.0018
## 
## $jack.boot.val
##  [1]  0.54101124  0.42195122  0.25513196  0.18618785  0.08888889 -0.07302452
##  [7] -0.18795518 -0.27631579 -0.41383285 -0.49350282
## 
## $jack.boot.se
## [1] 0.9994876
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
##    [1] 5.2 6.1 4.9 3.1 5.0 5.9 5.9 3.7 4.8 3.9 5.1 5.2 3.7 5.6 4.2 5.6 3.4 3.7
##   [19] 5.3 5.6 6.6 6.1 4.9 2.3 5.3 4.8 3.4 3.6 4.9 5.3 4.4 5.4 5.2 5.0 3.7 2.6
##   [37] 6.3 3.1 3.7 4.7 4.3 4.0 6.0 3.3 4.4 4.6 4.1 4.8 4.8 3.9 5.0 4.8 3.5 4.6
##   [55] 4.5 4.2 4.6 4.0 5.5 3.2 4.6 3.6 5.8 3.9 3.8 4.7 4.0 4.6 5.1 3.6 5.5 4.8
##   [73] 3.6 5.3 5.4 3.9 5.1 4.0 5.2 5.4 4.3 3.9 2.5 4.9 3.3 4.8 3.9 4.4 6.0 3.9
##   [91] 4.7 4.0 4.4 5.1 6.2 3.0 5.2 3.6 3.5 4.0 5.9 5.9 4.2 5.6 4.8 5.3 4.3 3.5
##  [109] 3.7 4.6 4.1 6.0 4.5 5.8 2.2 3.8 5.7 3.7 4.1 4.0 5.2 3.8 3.4 3.7 4.4 6.4
##  [127] 4.9 4.5 3.7 5.5 4.1 3.5 3.9 4.4 4.1 4.3 5.4 5.5 4.7 4.1 3.3 3.2 5.5 3.4
##  [145] 4.6 6.0 4.9 3.9 4.2 5.7 4.4 3.8 2.7 4.4 6.4 4.4 4.1 4.9 4.8 3.3 3.3 6.8
##  [163] 4.8 4.5 5.7 6.4 4.2 5.2 4.1 4.7 4.4 3.6 4.9 3.9 2.4 5.1 5.4 5.9 4.3 4.5
##  [181] 5.1 3.0 3.9 3.9 2.4 5.0 4.5 3.1 5.1 4.9 5.5 4.7 4.1 4.9 3.6 3.8 4.4 3.8
##  [199] 4.4 4.8 5.5 3.1 3.2 4.3 4.3 4.0 2.4 4.8 3.5 3.5 4.5 5.5 5.1 4.9 4.1 3.3
##  [217] 3.7 4.5 3.8 4.4 4.2 5.8 5.4 5.0 5.2 4.2 4.5 5.1 3.3 5.3 4.7 3.2 4.4 3.6
##  [235] 5.4 4.4 4.9 3.9 5.2 3.9 4.8 5.2 4.7 3.4 5.0 5.7 4.7 5.5 3.5 4.0 4.0 4.1
##  [253] 5.6 2.3 4.8 4.7 5.0 3.6 5.7 5.9 4.7 6.1 4.8 3.8 5.3 4.7 3.8 4.8 5.3 5.0
##  [271] 4.6 4.5 4.9 3.9 4.7 5.0 5.6 5.4 4.8 4.3 4.0 5.9 3.2 6.2 3.2 3.8 6.3 3.6
##  [289] 5.7 4.3 5.3 5.6 4.6 4.0 5.5 4.1 4.8 4.3 4.3 4.9 5.4 3.9 4.4 5.3 4.0 3.9
##  [307] 4.9 4.4 3.6 4.2 3.8 3.9 5.2 5.0 4.1 4.2 4.1 3.6 3.9 5.3 4.8 4.5 6.9 2.8
##  [325] 3.5 4.5 6.1 3.6 4.6 4.1 3.0 4.4 5.8 5.3 4.7 5.3 5.0 4.2 2.9 5.9 3.6 3.2
##  [343] 6.1 4.0 3.1 4.6 4.3 3.4 5.0 3.8 4.4 4.4 3.3 4.1 4.4 5.4 4.0 4.9 6.5 5.1
##  [361] 3.7 4.7 4.5 4.8 5.2 4.6 3.9 3.9 5.1 3.4 5.0 5.7 3.6 3.7 5.2 5.6 4.3 4.4
##  [379] 6.0 4.9 5.4 5.2 5.0 4.0 3.7 4.3 5.0 4.6 4.4 3.2 4.5 5.4 5.8 5.1 4.9 3.4
##  [397] 3.2 4.1 5.1 5.4 3.6 3.9 5.1 5.5 4.3 4.4 5.5 4.7 4.8 6.3 3.2 4.1 3.2 5.5
##  [415] 3.6 3.5 3.9 3.2 4.7 4.3 3.1 5.0 3.1 3.5 2.5 4.8 3.6 4.2 3.6 5.8 4.8 4.4
##  [433] 6.4 4.8 3.4 5.4 5.3 4.0 5.6 4.1 5.1 3.9 4.4 6.2 5.8 4.3 5.2 4.7 4.9 4.9
##  [451] 5.5 3.6 2.4 3.7 4.8 6.6 4.8 3.2 4.9 5.4 3.1 4.5 5.0 5.4 5.1 4.5 4.3 3.7
##  [469] 3.9 4.8 4.6 4.9 5.3 4.2 3.9 5.3 4.8 4.1 4.4 3.2 4.7 5.4 4.4 5.3 4.2 6.9
##  [487] 5.1 3.7 6.3 3.5 4.0 5.2 6.1 4.7 4.5 4.4 4.2 5.3 4.8 3.9 5.4 6.2 6.1 4.6
##  [505] 4.0 4.3 5.7 3.7 3.8 4.8 4.5 4.1 4.0 4.2 4.3 2.3 3.8 5.2 5.0 5.1 4.0 4.0
##  [523] 3.9 4.7 6.0 4.7 3.5 3.6 6.0 4.3 4.5 2.9 5.5 3.9 4.2 5.3 4.4 6.0 5.2 4.8
##  [541] 6.6 4.0 4.1 4.4 4.2 4.5 4.6 4.5 3.4 3.0 3.2 4.9 3.8 6.3 3.0 3.5 5.4 5.2
##  [559] 3.0 3.1 4.8 5.0 5.6 4.2 4.1 4.3 6.1 6.0 3.1 3.7 5.5 3.3 3.8 3.4 4.9 4.4
##  [577] 3.4 4.4 4.6 3.8 4.2 5.8 4.5 4.4 3.7 6.5 4.9 4.7 2.2 4.9 4.7 5.5 4.5 5.3
##  [595] 5.3 5.7 6.6 5.0 3.7 3.9 4.8 4.2 5.7 5.3 3.5 5.5 4.9 5.5 4.2 5.4 4.6 5.4
##  [613] 3.9 6.5 4.4 4.3 4.6 4.0 4.5 3.5 4.3 3.6 3.9 4.2 5.0 4.8 4.7 4.1 5.2 4.0
##  [631] 6.0 4.8 6.6 5.1 4.3 3.9 3.8 6.1 5.4 4.5 4.8 5.0 5.2 5.8 6.8 4.3 3.9 5.4
##  [649] 5.4 5.5 4.2 5.1 4.2 2.3 5.9 5.1 4.8 4.4 4.4 4.8 4.7 4.4 3.5 4.1 5.4 5.3
##  [667] 5.0 4.3 4.9 5.0 2.7 3.4 6.0 4.0 4.3 4.7 5.9 5.1 4.9 4.4 4.4 5.1 5.1 4.7
##  [685] 4.8 5.0 5.9 6.8 5.1 4.6 4.8 5.5 3.1 4.3 5.1 3.2 5.4 5.8 3.0 4.9 4.7 5.3
##  [703] 4.8 4.8 3.8 5.0 5.8 4.5 4.2 4.4 4.0 4.4 4.7 4.0 4.4 4.1 6.9 4.4 4.0 4.1
##  [721] 4.8 4.9 5.3 3.2 4.8 5.6 4.3 4.3 4.2 5.0 5.9 4.7 5.1 2.8 4.7 3.3 5.9 4.1
##  [739] 4.0 5.0 4.3 4.3 4.8 8.2 3.4 4.8 4.3 4.1 5.1 4.8 5.4 2.8 4.4 3.3 2.5 5.8
##  [757] 4.2 3.7 5.3 5.8 3.2 4.0 5.4 6.2 3.6 4.6 4.2 4.5 4.3 3.8 5.5 2.3 3.6 5.8
##  [775] 2.9 5.6 4.0 4.4 3.7 5.4 4.1 3.3 4.4 5.3 4.0 5.5 4.0 4.9 5.6 3.5 4.5 4.5
##  [793] 4.7 3.7 4.1 5.4 5.2 3.2 4.3 5.2 4.2 3.7 4.5 4.9 4.6 4.2 5.5 4.7 4.1 2.8
##  [811] 5.3 2.6 4.2 3.7 5.0 6.5 5.5 5.0 5.8 5.8 4.3 3.8 5.3 5.0 3.7 3.3 4.9 3.8
##  [829] 5.0 5.3 3.7 5.6 4.6 5.2 4.7 5.2 3.7 6.0 4.9 4.1 4.7 5.2 3.7 4.3 4.3 4.3
##  [847] 4.6 2.6 4.2 6.3 5.8 5.4 3.7 3.2 2.5 4.1 3.6 4.9 3.9 3.7 6.5 4.7 4.7 3.5
##  [865] 4.5 4.2 4.2 5.5 4.2 3.6 3.9 3.8 3.6 4.4 4.3 3.4 4.9 3.4 5.0 3.8 4.8 4.4
##  [883] 5.9 5.1 4.9 5.5 5.3 4.7 3.9 3.8 3.9 6.5 4.8 5.4 4.0 4.3 4.3 5.4 2.9 4.4
##  [901] 5.4 3.5 3.8 5.5 4.7 4.8 3.2 3.7 4.7 3.7 3.5 4.3 4.7 3.1 4.9 4.0 5.3 4.5
##  [919] 3.5 5.0 4.1 4.9 3.5 4.6 6.3 5.2 5.0 6.1 4.4 4.7 4.1 4.3 5.0 5.1 3.9 3.5
##  [937] 4.0 3.7 3.3 4.7 3.4 6.2 5.4 5.1 4.6 5.5 4.5 6.7 4.9 3.9 5.6 3.1 2.9 4.0
##  [955] 3.8 3.8 4.9 5.2 4.7 3.8 4.5 5.1 4.6 4.3 5.0 5.3 5.2 6.3 6.2 4.4 4.7 5.0
##  [973] 3.4 3.5 4.0 5.3 5.1 5.1 4.7 3.1 4.4 5.3 5.7 5.2 5.5 4.0 4.5 4.7 3.8 4.6
##  [991] 4.4 5.8 3.7 5.6 6.7 5.7 3.5 2.1 4.2 4.1
## 
## $func.thetastar
## 72% 
## 5.1 
## 
## $jack.boot.val
##  [1] 5.500 5.400 5.356 5.300 5.100 5.000 4.900 4.800 4.700 4.500
## 
## $jack.boot.se
## [1] 0.9466996
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
## [1] 0.6117542
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
##   2.0878184   2.3883387 
##  (0.8694457) (1.1235780)
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
## [1]  1.0938322  0.7363768  0.1858618 -0.2263001  0.6325541  1.4719044
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
##    [1]  0.139104081  0.737221090  0.953016743  0.091538711  1.475412090
##    [6]  0.743937636  1.132488565  0.173909494 -0.299613683  0.677432043
##   [11]  0.906511493  0.327032652  0.344617686  0.292336370  0.111181350
##   [16]  0.020221684  0.458099580  0.814177937  0.494253972 -0.032331686
##   [21]  0.172684347  0.531475560  1.162998674  0.284975288 -0.092182647
##   [26]  0.810610753  0.762724371  0.194785310  0.119890716  0.975557534
##   [31]  0.469330455  0.271411688  0.070184514  1.078811745  0.226687977
##   [36] -0.388740900 -0.121473481  0.172988009  0.301410917 -0.466120381
##   [41]  0.123932230  0.393947453  0.333541970 -0.239499370  1.044596907
##   [46]  1.192556297  0.385121646  0.714095783 -0.477892520 -0.213662834
##   [51] -0.702217928  0.872805344 -0.388700382  0.151744042 -0.021225416
##   [56]  1.417741821 -0.169450175  0.095310790  0.773816749  1.143316431
##   [61]  0.069483586 -0.314127803  0.144869937 -0.179004784 -1.487589032
##   [66]  0.518416635  1.515931858  0.211470833  0.680349123  1.484532759
##   [71]  0.559332540  0.847175205  0.627169138 -0.040321741 -0.532460740
##   [76] -0.461635934  0.118688110  1.108693600  0.636969185  0.875150112
##   [81]  1.065961506  1.330812778  0.291206457  0.969340214  0.545563250
##   [86]  0.205166991 -0.758267051  0.279014883  1.311547735  0.504236067
##   [91] -0.204643110  0.417462870  0.077407648  0.192543913  0.577556359
##   [96] -0.019154138  0.661601119  0.299203541  0.428294810 -0.129819005
##  [101]  0.678442740  0.112730437  1.796741826  0.021040831  0.777418318
##  [106]  1.072931790  1.201317342  0.721603021  0.820834451  0.667297710
##  [111]  1.117167166 -0.174981104  0.627281414 -0.245156722  0.672973538
##  [116]  0.236139613  1.092923633  0.843174267  0.972031364  1.484532759
##  [121]  0.474146368  1.065984646 -0.308469238  0.605053733  0.285744904
##  [126]  0.344245317 -0.816038150  0.075203859  0.164957762  0.880495207
##  [131] -0.124902112  0.136350235  0.581596315  0.755644458  0.723088071
##  [136]  0.443056149 -0.082774299  1.296645740 -0.371950032  0.747860779
##  [141] -0.266603602  0.838772977  0.886813408  0.500440670  1.080995628
##  [146] -0.278409803  0.121750314 -0.058742434  0.642475789  0.768789573
##  [151]  0.178302824  0.452382288 -0.875900441  0.403349332  0.766755785
##  [156]  0.333200380  0.749101416  0.809631211  0.397835065  0.271426930
##  [161]  0.092317507  0.098963728  0.947981781  0.032997042  0.445874015
##  [166]  0.586963231  0.221442672  0.393719569  0.588804907  0.108963711
##  [171]  0.422924988  0.932622856  1.009407893 -0.745640008  0.016663472
##  [176]  1.716818489  0.044559087  0.781236860  1.211941624  0.498251026
##  [181]  1.018842610  0.626632592  0.302912699  0.939342013  1.231358616
##  [186]  0.688065321  0.815373133  0.025295147  0.172687898  0.372323682
##  [191]  1.145306938  0.420095254  1.085304623  0.487243421  0.613778104
##  [196]  0.759089606  0.414283480  1.190864748  0.581760615  1.371233101
##  [201]  0.654678955 -0.732467896  0.615248554  1.012940845  0.669830797
##  [206]  0.125210265  0.235717656  1.632566805  0.341315141  0.583624191
##  [211]  0.384930570  0.975615862  0.707317928  0.112220250  0.590170997
##  [216] -0.455723358  0.692713735  0.444067063  0.925155860  0.885836005
##  [221]  0.706112981  0.569647817  0.079626503 -0.312690766  0.626696791
##  [226]  0.692196218  0.533081911 -0.303983445  0.383478951  0.442402765
##  [231]  1.023777875  0.294291574 -0.413697670  0.076476731  0.382474327
##  [236] -0.909597294  0.971937482  0.485719280  0.769647799  1.368698625
##  [241]  0.939342013  0.388755851  0.827539942  0.519526932 -0.759566444
##  [246]  0.478897000  0.817409569  0.521527197  0.173646531  0.941702481
##  [251] -1.317040290  0.175974081 -0.135919599  0.434288840  1.386979944
##  [256]  0.590586195  0.075616334 -0.144306801  0.411686472  0.349587638
##  [261]  0.177255722 -0.392252565  0.209242127  2.059017591  1.313907123
##  [266] -0.268783304 -0.736623613  0.796181721 -0.164555732 -1.140274030
##  [271]  0.180863596  0.006915427  0.089903935  0.486761333  0.311097116
##  [276]  0.167371060 -0.162445380  0.862225575 -0.410922680  0.423469377
##  [281] -0.227913748  0.437021668  0.153145054  0.664798937  0.182095537
##  [286]  0.784715174  0.731771846  0.457325765  1.324080714  1.433658229
##  [291]  0.941317608 -0.447586477  1.002128472  0.213547334 -0.063661946
##  [296] -0.247580843  1.303439462  0.559214887 -0.112659352  0.003702752
##  [301]  0.169045798  0.082393565  0.551792268 -0.070296095  0.528280268
##  [306]  0.304827652 -0.975223536  0.533443887  0.953458109  0.594409028
##  [311]  0.277675073 -0.136799667  0.914282219  0.599018237  1.068644579
##  [316]  0.538696615 -0.490458128  0.401629557  0.627256963  1.133405481
##  [321]  0.516233231  0.449104304  0.701396357  1.396529580 -0.520472833
##  [326]  0.325620812 -0.408971183  0.649450553  0.202028342  0.947450693
##  [331]  1.651423225  0.114054006  0.658568462  0.909973051  0.938042739
##  [336]  0.556645225 -0.104120377  0.264030795 -0.697318325  0.439246947
##  [341]  0.385718617  0.662088467  1.394311632  0.622294421 -0.332083514
##  [346] -0.106500811  0.962892474  0.733677965  0.944059267 -0.909597294
##  [351]  0.307462206  0.716740572  0.337318541  0.584847160  0.402440158
##  [356] -0.243331636  0.625513673  0.745146279  0.515337650  0.612838239
##  [361]  1.188090217  0.947145796  0.265304847  0.326353090  0.274671819
##  [366]  1.297226481  0.224860579  1.438361814  0.427706749  1.227559355
##  [371]  0.636972988 -0.085044607  0.311478463  0.542341472  0.252884925
##  [376]  0.534476982  1.004208496  0.301215575  0.420185975  0.804074420
##  [381]  0.249821793  0.122493082  0.641588954  0.157118042 -0.432983822
##  [386]  0.402119570  0.167183806  0.869447673 -0.590940522  0.465259534
##  [391]  0.194840519 -1.452709775  0.578815918  1.122634446  0.425730852
##  [396] -0.082757522 -0.513923697  0.892207655  1.409307208  0.460716627
##  [401]  0.963655329  0.191428592  0.106774223  0.029971286 -0.148263283
##  [406]  0.597604239  0.652183357  0.512993742  0.443041549  0.888333843
##  [411]  0.537020321 -0.185501778 -0.255024224  0.189738708  1.810960451
##  [416]  0.920969765  0.535364570  0.427194846  0.832613225  0.427183611
##  [421]  0.104783789  0.696941816  0.657492754 -0.807203248 -0.484352746
##  [426]  0.825816880  0.379977828 -0.480099504  0.049332133  0.017943065
##  [431] -0.266752339 -0.457709689  0.501914204  0.400657177  0.569015583
##  [436]  1.074745649  0.509474796  0.250577491  0.504476074  1.040711624
##  [441]  0.113774389 -0.009093367  0.198173894 -0.487465635  0.871442026
##  [446] -0.286809331  0.561397783 -0.180324790  0.682521602  0.860550847
##  [451]  1.052379248 -0.228779459  0.523883104  0.214081114  0.336273220
##  [456]  0.093595217 -0.541192299  0.544008286  1.121567537  0.740608300
##  [461]  0.133881945 -0.003618894  0.871019878 -0.763963920  0.736729276
##  [466]  0.299797645  0.662446878  0.466811472  0.694298242  0.269608125
##  [471]  1.106640097  0.298587194  0.286487625  0.708992123  0.590356128
##  [476] -0.064748265  0.712688225  0.499698071  0.903480094  0.231970591
##  [481]  0.049562212  0.725791253  1.481416045  0.495900989  1.462111719
##  [486]  0.824275444  0.450834817 -0.277531894  0.925126341 -0.548002736
##  [491]  0.444084651  0.173487335  0.330635118  0.345091388  0.503405010
##  [496]  1.053800759  0.791832640 -0.252439920  0.474080091 -0.665533223
##  [501]  0.780178396  0.227176970  0.097448696  1.025651246  1.148088594
##  [506]  0.420426101  0.354953347 -0.213663688  1.341968364  0.201931107
##  [511]  0.363028005 -0.349953758 -0.848190870  0.827628068  0.616298154
##  [516]  0.036747109  0.513597055  0.931214443  0.549148396  0.465681784
##  [521]  0.301007872 -0.728313225 -0.065284886 -1.215523848  0.164990882
##  [526]  0.780762430  1.480299906  0.280743154  0.312131347  0.119565178
##  [531]  0.394296979  0.442381190  0.545110119  0.776040338  0.870405896
##  [536]  0.386048928  0.112683211 -0.361936810  0.580876941 -0.294964167
##  [541]  0.052621840  0.872416044  0.490475830  0.876970096 -0.093080465
##  [546]  0.373170647  0.667796104  0.297924919  0.214598919  0.187529980
##  [551]  0.146681661  0.349279129 -0.204646546  0.410415984 -0.090640414
##  [556]  0.630241525  0.441190925  0.115901604  0.756191805  0.620274718
##  [561] -0.040181443 -0.595826129  0.379573773  0.899823169  0.460663593
##  [566] -0.415746477  0.152812622  0.865701956  0.476479373  0.018395343
##  [571]  0.059591426  0.395338273  0.190641946  0.602735494 -0.996118616
##  [576] -0.367695530 -1.515299637  0.316422223  0.641790459  0.441810151
##  [581]  0.043479228 -0.077641425  0.406413180  1.322462621  0.751001835
##  [586]  0.244498328  0.817323057  1.016231625  0.596728816  0.652842548
##  [591]  0.820162555  0.679592967  0.197688475  0.716343392  0.302157937
##  [596]  0.405653460 -0.465049408  0.229003304  0.300764407 -0.355799361
##  [601]  1.112585512  0.090708770  0.991674815  0.484463866  0.514933820
##  [606]  0.896592299 -0.036739676  0.767711648 -1.296290279  1.635535907
##  [611] -0.377149158  0.409639232  0.530688800  1.247951040 -0.116247117
##  [616]  0.293735006  0.281078847  0.860557302  0.032667693  1.653562188
##  [621]  1.451788845  0.349572587 -0.110341021  0.718307149 -0.319743376
##  [626]  0.561427456  0.288705068 -0.316474491  0.271077157  0.800937220
##  [631]  0.113703685 -0.706938594  0.352882413  0.874541453  0.674927536
##  [636] -0.115662858 -0.284007724  0.739315977  0.700397828  0.602076138
##  [641]  0.561208197  0.342734962  0.764807153  0.358473685  0.639240909
##  [646]  0.669142600 -0.007946896  0.434826546  0.715421611  0.098508337
##  [651]  0.594457151  0.112866939  0.756242021  0.963758123  0.341958080
##  [656]  0.763278879  0.155461565 -0.069157425  0.179407348  0.803495423
##  [661]  0.029179241  0.823979418  0.386373852  0.659323940  0.356419644
##  [666]  0.491756934  0.060612651  0.710851836  0.900556523  0.646870514
##  [671]  0.964332047  0.153926585  0.457556233  0.530752541  1.130994489
##  [676]  0.203812718  0.410946083  1.103528034  0.516336964 -0.216399785
##  [681]  0.553632874  0.275469860 -0.097494053  0.171548717 -0.395378594
##  [686]  0.283790540  0.870242613  0.403349332  0.581760615  0.398074147
##  [691]  0.396074582  0.409701798  1.215052295  0.476272660  0.522018489
##  [696]  0.191235842 -0.421298584 -0.270772759  1.639464281  0.636545232
##  [701]  1.005576547  0.853792625  1.111278427  0.386747334  0.508873741
##  [706]  1.254592040 -0.244810887  0.619320934  1.394311632  0.758260306
##  [711] -0.019621213 -0.728764343 -0.780252263  0.960137763  0.298347274
##  [716]  0.396742975  0.742349529  1.332958042  0.487119907 -0.125036439
##  [721]  0.464777842  0.298957594  0.473666615  0.622472414  0.316949659
##  [726]  1.124645595 -0.367820296 -0.008199561 -0.508862514  0.611754153
##  [731] -0.038790071 -0.740270844  0.850692398 -0.110972848  0.772606270
##  [736]  0.389027595  0.884877256  0.649622557  1.549008204  1.454199472
##  [741]  0.440654219  0.663762772 -0.487279404  1.682296829  0.257911914
##  [746]  0.324339722  1.299481640  0.832833941  0.319845390  0.101901856
##  [751]  0.959623921  0.277122832  0.831504248  0.500432485  0.485941418
##  [756] -0.696470138  0.566865076  0.243900113 -0.050982395  0.887780340
##  [761]  0.662182359  0.117719935  0.656152343  0.621757238  0.554366996
##  [766]  0.176174409  0.127082798  0.460716627  0.937847807  0.291053516
##  [771]  1.158030941  0.513999991 -0.067616302  0.160950085  0.742563865
##  [776]  0.206544148  1.720490352  0.149870706  0.505606857  0.181662757
##  [781]  0.083275951  0.276210271  1.081560640  0.921175880 -0.046625401
##  [786]  0.767266868 -0.258204196 -1.008063693  1.213700036  0.653863944
##  [791]  0.020833812  0.114957504  0.639934416  1.093765436  0.158085508
##  [796]  0.459045941 -0.066549162  0.027162252  0.632712223  0.610015402
##  [801]  0.742497084  0.181187070  0.583388766  0.474712099  0.740790141
##  [806]  0.707609010  0.912713291  0.559098953  0.682583584 -0.327938467
##  [811]  0.632662698  0.120367960 -0.423733886  0.824275444  0.754945476
##  [816] -0.351581138  1.132276243  0.324232873  0.362544080  1.251434122
##  [821] -0.204643110  0.447867992  0.113544915 -0.787648286 -0.019096389
##  [826]  1.087860910  0.827610630  1.929744511  0.851533284 -0.731977437
##  [831]  0.446315596  0.804030683 -0.320553690  0.547762219  0.138227809
##  [836] -0.376343834  0.586114167  0.798259116  0.180879553  0.456541953
##  [841] -0.005658176  0.531832647  1.489576097  0.693765043  0.170417616
##  [846]  0.420208643  0.011959009 -0.380168647  0.131955327  0.675488043
##  [851]  1.021563592  0.254369221  1.158548579  0.782264708  0.099338817
##  [856]  0.802900085  0.180990749  0.972469311  0.503465284 -0.019275189
##  [861]  1.214465081  0.268203087  0.534984268  0.224167267 -0.275514297
##  [866]  0.660157783  0.490597646  0.792760328  0.582659982  0.200382474
##  [871]  0.420218329  0.220058808  0.512478211  0.484034590  0.568358193
##  [876]  1.215675759  1.094261599  1.082883571  0.069314218  0.132177520
##  [881]  0.872440979  0.336091583  0.158119058  0.569373299  0.176613752
##  [886] -0.008889769 -0.051774147  0.449753289  0.257535754 -0.108816259
##  [891]  0.888333843  0.863518931 -0.519691437 -0.377164136 -0.091570040
##  [896]  0.710973347  0.729116412  1.014344849  0.696478936  0.726850232
##  [901]  0.638447219  0.614209605 -0.822266129  0.491679216  0.602251355
##  [906]  0.505428167  0.223868474  0.123317920  0.699179714 -0.677169357
##  [911]  0.425865109  0.712082963  0.717436847  1.105695267  1.042370045
##  [916]  1.078983985 -0.075436667  0.452337696  0.665344347  0.157230813
##  [921]  1.063992472 -0.184356125  0.267989407 -0.013638574  0.464106810
##  [926]  0.289075280  0.001883019  0.035793006  1.149085538  0.312113627
##  [931] -0.715499365 -0.172055332  0.392895877  0.466819840  0.044352457
##  [936]  0.429985163 -0.564972796  0.194109499  0.138724605 -0.036005422
##  [941]  1.646853872  0.226576880  0.087525388  0.217340647 -0.698190742
##  [946]  0.796267339  0.421062056  1.020600459 -1.229510513  0.503511859
##  [951]  0.939325028  0.543510364 -0.856291519  0.615665173  0.399273336
##  [956] -0.199761038  0.264209699 -0.331346170 -0.594435445  0.334856597
##  [961]  0.544008286  2.260418615  0.811985769 -0.809422488  0.728871718
##  [966]  0.322193624  0.224256028 -0.253857113  0.504331646  0.733718508
##  [971]  0.218262108  0.000573694  0.641576321  0.422330696 -0.272128753
##  [976]  0.068517203  0.660682654  0.856481395 -0.202396740  1.455513761
##  [981]  0.581723763  0.409701798 -0.801141594  0.666118731  0.695559435
##  [986]  0.173727200  1.043228149  0.829047063  0.341425631 -0.107510060
##  [991]  0.467487155  0.779768397 -0.571433637  0.096182619  0.595608677
##  [996]  0.672147284  0.186587939 -0.004136920  0.489361994  0.632214615
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
##      mean         sd    
##   0.8741719   0.5786098 
##  (0.1829725) (0.1293775)
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
## [1]  0.3081841  0.3041269 -0.3656239 -0.4109135  0.8599588  0.7284039
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
## [1] 0.0425
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8927378
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
## t1*      4.5 -0.04674675   0.9156329
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 2 4 7 8 9 
## 1 1 2 2 2 2
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
## [1] -0.0044
```

```r
se.boot
```

```
## [1] 0.8992493
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

