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
## 1 4 5 6 7 8 9 
## 1 1 2 1 2 1 2
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
## [1] 0.0097
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
## [1] 2.734567
```

```r
UL.boot
```

```
## [1] 6.284833
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.8   6.3
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
##    [1] 4.6 3.9 3.2 3.1 4.5 2.7 4.9 5.5 5.5 6.2 7.5 2.8 4.5 3.6 4.7 3.2 3.9 4.4
##   [19] 3.0 3.4 3.6 3.6 4.5 3.3 3.1 3.7 4.7 5.2 5.3 5.1 4.1 4.8 5.0 4.7 4.7 5.2
##   [37] 5.5 4.5 3.2 4.4 4.6 4.9 3.5 6.8 3.3 5.1 4.7 4.0 2.7 3.7 5.8 4.1 4.7 4.7
##   [55] 4.3 3.2 4.4 6.2 4.7 4.0 4.7 2.2 5.3 3.9 5.0 4.5 3.4 3.7 3.2 5.3 4.0 3.8
##   [73] 4.4 4.3 2.9 4.6 4.7 5.4 4.2 6.1 5.2 4.5 4.5 5.5 5.1 4.9 6.8 6.4 3.5 5.4
##   [91] 4.8 3.6 5.1 4.5 3.4 2.9 3.2 4.5 5.6 4.7 4.4 5.3 4.4 2.2 4.2 3.3 5.2 6.4
##  [109] 3.1 5.0 4.7 6.2 4.4 4.3 4.7 5.1 2.6 5.9 4.4 3.8 5.4 3.8 3.8 4.8 4.7 4.6
##  [127] 4.3 4.0 5.1 5.7 5.7 5.0 3.9 4.2 5.2 4.1 3.7 4.5 4.8 4.8 4.0 3.2 3.5 4.4
##  [145] 4.1 3.5 3.2 3.8 4.1 4.6 5.1 4.7 4.0 4.9 4.6 4.7 2.4 5.5 4.3 4.9 4.7 4.9
##  [163] 5.3 3.9 4.6 3.2 5.7 4.4 3.3 4.1 4.7 4.7 5.4 5.2 5.0 5.0 3.8 3.6 6.2 5.6
##  [181] 4.7 3.6 5.6 4.6 4.7 4.6 3.5 5.4 3.7 5.9 5.4 5.0 5.3 3.7 3.0 4.3 5.1 3.8
##  [199] 5.5 4.6 4.1 5.5 4.1 4.0 4.4 3.6 4.4 2.6 4.1 4.8 4.6 5.5 5.1 4.1 5.1 5.9
##  [217] 5.7 4.8 3.0 4.5 5.3 5.9 5.7 2.8 4.2 5.0 5.5 4.3 4.7 5.1 5.6 5.1 3.0 2.9
##  [235] 2.6 4.0 3.5 4.4 4.1 4.8 4.1 4.1 4.4 4.3 2.3 5.6 3.6 6.1 4.3 3.8 4.2 3.3
##  [253] 4.4 5.5 2.8 4.1 3.2 3.9 5.0 3.8 4.6 4.0 3.1 3.9 3.7 4.9 3.9 4.0 4.9 4.1
##  [271] 6.0 4.1 4.6 5.8 5.3 3.5 5.3 5.0 5.5 4.3 2.5 4.4 2.8 4.6 4.8 3.1 3.4 5.0
##  [289] 4.2 4.0 4.5 4.1 4.2 3.4 5.2 5.4 2.1 5.5 3.5 4.9 5.6 5.2 3.9 3.6 3.9 2.9
##  [307] 4.4 5.2 4.1 4.5 3.8 4.5 4.8 6.7 5.4 3.5 4.8 4.9 3.0 3.8 4.5 4.7 5.5 5.2
##  [325] 5.6 2.5 3.0 5.2 4.2 5.2 4.5 5.1 3.4 3.6 4.9 2.9 4.3 5.4 3.7 5.3 4.8 4.1
##  [343] 5.3 3.9 4.0 3.9 4.6 4.3 6.5 4.4 4.6 6.7 4.3 4.5 5.7 3.9 4.5 4.1 2.2 5.1
##  [361] 4.3 4.2 5.1 5.7 5.0 5.4 3.1 2.9 5.0 4.1 4.4 3.4 3.6 4.4 3.0 5.2 5.5 3.5
##  [379] 4.2 3.6 3.9 5.1 3.8 3.7 4.9 3.9 6.0 4.3 3.3 3.0 5.2 5.2 4.3 3.8 3.3 5.1
##  [397] 4.9 5.2 5.4 5.0 4.7 4.4 5.6 5.3 3.1 4.1 5.0 3.9 4.3 3.4 6.1 3.4 3.7 5.6
##  [415] 5.8 3.7 3.7 3.9 3.8 6.4 4.1 4.8 3.9 4.1 5.1 4.8 5.3 5.0 5.8 4.9 4.4 3.7
##  [433] 4.1 4.0 4.3 3.6 3.6 5.0 3.6 5.3 4.2 3.7 4.4 4.8 5.0 4.5 4.7 5.0 5.4 5.6
##  [451] 4.5 4.9 5.0 4.1 5.3 5.0 4.2 3.7 4.8 3.9 5.8 4.2 2.6 4.7 4.0 3.2 4.2 3.4
##  [469] 4.6 4.0 2.4 5.0 4.9 4.9 4.1 5.4 5.0 4.9 5.7 4.3 4.9 2.2 4.3 4.9 5.0 4.4
##  [487] 4.3 4.3 4.5 3.8 5.6 3.3 4.3 4.9 3.4 3.0 5.9 4.1 4.2 3.4 4.4 3.6 3.4 4.1
##  [505] 4.1 2.3 4.3 3.8 5.6 5.7 3.8 5.0 4.8 6.1 4.7 4.0 4.2 3.9 4.3 4.5 4.4 5.1
##  [523] 3.6 5.2 5.5 4.5 5.3 2.9 3.1 4.7 4.8 3.8 4.7 5.0 4.5 5.5 3.9 5.3 4.3 5.6
##  [541] 4.1 3.5 4.9 5.4 3.5 4.6 6.3 5.3 2.9 5.4 4.7 5.4 4.8 5.1 4.7 3.4 6.3 3.9
##  [559] 4.0 5.6 4.4 5.3 3.4 2.6 5.7 3.9 3.5 5.3 4.8 5.1 5.9 3.2 4.4 5.6 4.6 3.7
##  [577] 3.5 4.8 3.2 4.6 4.5 3.6 3.4 3.8 4.4 5.1 5.2 4.4 2.8 4.0 4.9 5.1 4.2 5.5
##  [595] 3.2 4.1 3.2 4.2 5.3 4.8 3.8 5.3 3.5 3.6 4.5 5.2 5.4 4.3 4.7 5.0 4.7 4.9
##  [613] 4.0 4.7 3.9 4.3 3.9 4.6 4.6 4.9 5.5 4.1 6.3 4.1 4.6 5.1 5.3 3.9 4.3 5.5
##  [631] 5.2 4.4 4.4 4.9 4.4 5.6 2.8 5.9 3.1 5.1 5.2 4.6 3.7 4.4 2.9 5.8 4.7 2.3
##  [649] 4.5 5.5 3.5 3.0 4.9 3.6 4.3 4.3 4.8 4.4 3.7 5.3 4.6 3.3 4.3 3.2 6.2 5.4
##  [667] 3.9 4.3 4.3 3.8 4.7 4.6 3.9 5.6 3.8 4.4 4.3 4.8 4.3 5.8 5.0 5.3 2.6 4.9
##  [685] 4.1 5.1 5.1 5.9 4.4 4.5 4.8 3.2 4.9 5.1 2.5 2.9 4.8 4.5 4.5 3.8 5.1 4.3
##  [703] 6.0 3.5 4.9 4.3 3.2 4.7 4.2 5.4 4.9 5.6 4.7 6.2 4.6 3.6 4.1 5.0 3.7 3.9
##  [721] 4.2 4.6 5.6 5.9 6.0 4.4 4.8 4.9 3.5 5.8 4.6 4.0 4.3 5.1 4.4 4.6 5.1 5.4
##  [739] 5.5 4.8 3.3 5.2 4.7 4.3 4.3 4.7 4.4 5.1 5.4 4.1 3.0 3.8 4.6 6.1 4.2 4.9
##  [757] 4.0 3.6 2.7 5.0 3.5 4.5 4.7 3.4 2.7 4.3 4.9 3.5 4.8 4.3 4.9 5.3 4.5 4.5
##  [775] 4.6 4.0 4.6 4.7 4.8 5.8 3.8 5.5 6.3 4.2 4.3 2.9 6.2 5.5 4.8 4.5 3.5 4.8
##  [793] 3.9 4.5 4.5 5.6 6.4 4.3 3.3 3.2 5.9 5.6 4.8 2.8 6.2 4.6 4.3 3.9 3.9 4.1
##  [811] 5.0 4.8 4.9 6.7 4.7 5.1 3.3 4.7 3.8 5.3 5.6 4.1 5.1 4.1 3.2 5.5 2.5 4.8
##  [829] 4.4 3.1 5.2 4.0 3.7 4.0 3.7 6.2 3.5 4.5 2.6 4.8 4.0 5.8 4.0 5.3 4.1 3.7
##  [847] 5.8 4.4 3.8 5.7 5.0 4.8 2.3 4.0 5.2 4.9 4.2 3.0 3.5 5.0 3.2 5.2 4.7 6.8
##  [865] 4.1 3.9 4.6 4.9 4.6 4.2 3.3 2.6 3.9 5.6 4.3 3.8 4.2 4.5 4.7 3.5 4.3 4.5
##  [883] 4.0 6.0 5.1 4.4 4.0 3.7 4.1 5.5 5.2 5.1 3.8 3.2 3.1 2.0 5.1 5.2 3.7 5.0
##  [901] 4.0 4.8 6.3 5.3 4.3 3.7 6.0 4.5 3.8 3.9 4.8 4.5 5.0 4.7 4.0 4.7 5.6 5.3
##  [919] 4.8 3.4 4.2 5.3 4.5 4.5 5.6 3.4 5.7 5.1 4.7 3.6 4.7 5.2 2.7 3.5 4.7 2.9
##  [937] 4.1 3.4 4.4 4.1 4.2 5.7 4.6 4.3 5.1 4.4 4.0 5.9 4.4 5.4 4.2 3.7 3.4 4.8
##  [955] 7.2 3.9 4.4 3.2 5.0 3.8 5.5 4.7 3.2 4.3 3.1 3.5 5.5 5.2 6.6 4.4 5.0 4.2
##  [973] 4.0 4.7 4.4 4.6 5.7 6.0 3.2 6.0 4.4 2.5 5.6 5.0 5.1 3.0 3.5 3.7 4.5 5.0
##  [991] 4.7 4.3 3.6 4.6 4.9 4.4 4.5 4.8 4.9 6.2
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
## 2.6975 6.2000
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
##    [1] 5.0 3.9 5.9 7.1 5.1 5.7 4.6 6.1 3.4 4.1 3.4 2.8 3.0 3.8 3.8 2.7 4.7 4.2
##   [19] 6.0 5.2 4.3 3.8 3.3 3.5 4.6 4.5 3.8 2.9 4.9 3.7 3.1 5.1 2.0 4.5 4.7 3.2
##   [37] 4.2 5.2 4.0 3.8 5.5 4.6 5.1 5.0 4.7 5.5 6.6 4.3 6.5 6.2 3.5 2.5 4.0 3.7
##   [55] 5.5 4.6 5.7 4.4 5.5 5.2 5.4 4.1 4.3 3.5 5.7 5.2 4.5 5.5 4.8 4.5 5.1 4.3
##   [73] 5.6 4.5 5.1 3.6 3.5 4.3 4.2 4.1 4.9 6.7 4.9 4.5 2.9 3.6 4.5 6.2 7.1 4.2
##   [91] 3.0 5.4 5.0 3.9 5.5 5.6 4.3 4.4 5.5 4.5 4.3 4.2 4.6 6.8 4.6 3.4 3.4 4.1
##  [109] 4.0 5.2 3.9 4.8 3.5 5.3 3.3 4.2 4.5 5.0 4.6 4.8 4.6 3.4 4.8 5.3 4.8 3.1
##  [127] 4.9 5.4 4.2 5.8 3.1 5.3 3.5 5.8 4.6 4.3 5.9 3.5 5.0 4.4 3.2 4.3 3.8 3.0
##  [145] 2.6 4.2 5.0 5.2 4.1 3.7 6.6 4.1 3.6 4.8 5.3 4.9 3.7 4.9 3.5 3.7 3.6 4.3
##  [163] 5.8 4.0 5.2 4.8 3.7 3.4 4.7 3.6 3.8 3.9 3.1 4.2 4.8 5.5 6.0 5.2 4.8 4.1
##  [181] 5.4 4.8 5.6 4.6 3.7 3.7 3.9 4.5 4.2 4.6 3.8 3.5 5.7 3.6 3.9 5.7 4.9 5.2
##  [199] 3.5 4.3 5.7 5.2 5.5 6.3 5.7 4.6 4.8 6.1 4.0 3.6 3.3 3.7 4.1 3.4 3.8 3.3
##  [217] 4.4 3.8 3.6 5.0 5.0 4.4 3.8 5.6 3.2 3.5 4.4 4.5 4.1 4.8 4.1 3.0 4.5 3.5
##  [235] 5.5 3.0 4.4 4.9 4.6 4.3 3.6 4.6 3.6 3.6 3.2 3.9 4.7 4.8 6.6 4.0 4.5 4.7
##  [253] 4.8 5.1 3.1 5.7 4.3 5.4 5.0 5.5 4.9 2.8 3.7 4.1 2.5 4.3 3.4 4.7 3.8 5.3
##  [271] 4.7 5.6 5.8 3.3 5.8 5.9 4.5 6.5 4.6 4.7 4.2 4.9 4.6 4.5 3.7 5.4 4.2 4.3
##  [289] 4.7 3.3 4.4 5.5 3.1 4.5 5.9 4.6 4.1 5.3 5.4 5.3 3.6 3.9 4.0 5.2 4.6 4.3
##  [307] 4.5 3.7 6.0 4.4 4.8 5.4 4.6 4.1 3.0 6.1 5.6 5.0 4.3 4.2 5.6 4.2 3.6 5.5
##  [325] 4.1 3.4 5.0 4.5 6.3 4.7 5.1 3.2 5.0 3.6 3.7 5.2 4.0 4.6 5.5 4.0 4.7 3.1
##  [343] 5.5 5.2 4.0 4.0 3.1 4.1 4.6 5.4 5.1 4.7 3.1 2.5 5.3 5.6 3.3 2.7 4.5 4.3
##  [361] 3.8 5.5 5.2 4.6 4.8 2.7 4.4 4.7 3.1 4.9 5.0 3.5 3.9 4.8 4.9 4.8 4.4 4.5
##  [379] 4.8 3.9 4.5 4.5 6.0 3.1 5.8 5.8 6.2 5.0 4.6 4.0 3.9 5.5 4.8 5.1 4.0 4.5
##  [397] 3.7 4.5 4.6 5.1 4.4 4.8 4.4 5.7 5.2 3.2 3.5 3.9 4.8 4.7 4.6 3.5 4.5 4.0
##  [415] 5.7 5.2 3.7 3.3 5.6 4.3 4.9 3.4 5.7 3.5 3.3 2.4 4.2 5.7 5.0 5.5 3.9 4.3
##  [433] 4.3 4.5 5.1 5.9 5.2 3.1 5.4 5.5 3.7 2.7 4.2 6.2 4.2 4.0 4.2 4.9 4.8 2.4
##  [451] 3.2 6.2 2.1 5.4 5.5 4.5 4.9 2.6 4.3 4.8 5.0 3.2 4.1 4.6 5.3 3.9 4.1 4.0
##  [469] 5.3 4.1 4.5 4.6 4.5 3.8 4.2 4.8 4.4 4.7 5.0 4.9 3.6 3.6 4.3 6.1 4.6 2.3
##  [487] 4.9 5.4 5.4 3.5 4.7 4.5 5.3 3.3 5.9 4.6 5.7 3.4 4.1 6.8 4.6 2.9 5.8 2.0
##  [505] 3.6 5.2 4.0 5.4 3.1 5.0 4.6 5.3 3.9 4.4 5.4 4.8 4.6 4.5 4.9 4.2 1.8 4.2
##  [523] 3.1 4.3 5.2 3.3 5.1 3.1 3.7 5.1 3.8 3.7 3.2 4.5 4.5 5.2 5.3 5.6 3.8 4.2
##  [541] 4.1 4.2 4.4 4.1 4.2 4.8 4.6 4.6 4.8 1.9 5.4 5.3 4.7 5.7 5.1 4.5 4.2 4.6
##  [559] 5.0 4.2 5.2 5.2 6.0 3.4 3.9 3.9 4.3 3.5 4.8 4.3 4.4 2.9 5.0 5.4 4.4 3.7
##  [577] 4.8 3.8 4.2 4.8 5.6 4.5 4.1 3.3 3.3 5.0 5.1 3.9 4.9 5.1 5.9 3.3 4.2 4.7
##  [595] 4.4 4.7 4.1 3.8 3.8 6.1 4.5 4.7 5.0 4.0 4.4 5.0 4.1 4.6 3.7 7.3 3.6 3.2
##  [613] 5.7 4.8 4.4 5.2 4.0 5.3 4.6 3.7 5.1 5.3 3.9 5.0 4.1 4.8 2.6 6.4 5.0 5.2
##  [631] 6.1 3.3 5.9 5.0 4.9 5.4 3.9 4.9 4.6 5.2 6.0 4.2 5.4 3.6 3.1 5.6 5.6 4.3
##  [649] 5.0 4.0 3.2 4.1 4.7 4.8 4.7 5.1 3.8 3.9 5.6 3.2 3.7 4.8 4.9 4.2 6.2 4.8
##  [667] 4.8 3.7 5.8 3.8 5.2 4.7 4.9 5.3 2.6 6.2 5.1 5.9 6.1 4.0 4.8 4.3 4.0 3.3
##  [685] 5.5 4.1 4.0 4.4 3.0 5.0 5.8 4.0 5.7 4.8 5.6 4.3 3.0 5.7 4.6 4.9 4.2 3.5
##  [703] 4.6 4.3 3.8 3.9 6.2 4.0 4.2 5.6 4.3 4.6 4.3 4.1 5.4 3.7 4.1 5.1 3.0 4.6
##  [721] 3.3 3.3 3.8 4.8 5.0 3.6 4.8 4.3 4.5 3.2 4.2 5.9 5.8 3.1 3.7 5.5 7.0 2.3
##  [739] 4.0 3.4 5.0 2.6 3.9 5.9 4.3 4.9 4.4 5.0 5.5 4.2 4.1 5.2 2.9 4.0 4.9 5.4
##  [757] 3.6 5.0 5.8 3.9 5.1 4.7 3.1 6.2 3.9 4.8 4.6 5.6 4.6 4.2 3.9 3.7 5.5 4.3
##  [775] 4.8 3.4 5.7 5.8 5.0 5.3 4.9 3.6 4.4 3.9 4.3 5.1 4.3 4.1 4.5 5.7 4.7 3.8
##  [793] 3.7 4.4 4.9 4.3 4.8 4.1 4.1 4.8 4.2 5.4 3.5 5.4 5.3 5.7 3.8 6.9 3.8 5.4
##  [811] 3.9 3.8 3.5 4.5 4.2 5.2 4.4 5.0 4.6 3.8 4.5 5.3 3.5 5.6 2.2 5.8 4.6 3.8
##  [829] 5.2 6.0 3.6 4.2 3.9 5.6 2.7 5.6 5.7 5.3 4.7 3.8 4.1 4.5 5.0 3.3 4.0 4.5
##  [847] 4.1 5.0 6.2 4.1 4.9 5.0 5.0 3.6 5.7 5.2 3.3 3.1 3.4 4.5 2.8 4.9 4.9 2.3
##  [865] 5.2 3.3 5.0 4.9 3.8 6.0 5.7 4.5 4.9 3.9 4.9 5.6 5.2 6.0 4.7 5.0 6.3 3.8
##  [883] 4.5 4.5 3.8 4.2 4.4 4.8 3.6 6.6 4.9 3.1 2.6 3.7 5.1 3.9 5.0 6.4 3.8 3.8
##  [901] 5.0 5.9 4.6 6.2 4.5 5.3 5.7 4.8 5.3 5.5 6.0 4.6 4.7 3.4 4.5 5.2 4.2 6.8
##  [919] 5.0 4.5 6.4 5.2 4.9 4.9 5.8 5.3 3.8 4.6 4.8 3.6 4.2 6.3 3.3 2.6 4.0 5.1
##  [937] 3.5 6.1 4.5 3.8 4.7 3.6 4.3 4.4 3.3 5.2 4.4 4.7 5.6 5.0 1.9 5.4 6.1 4.1
##  [955] 4.0 3.4 2.7 5.0 3.3 4.5 5.5 4.4 5.1 3.8 3.3 4.8 3.7 3.2 5.4 3.7 4.1 4.3
##  [973] 3.4 5.1 4.7 4.1 4.9 2.3 3.4 4.2 3.5 3.9 4.7 5.5 2.4 4.9 4.5 5.6 3.5 4.6
##  [991] 5.4 4.6 5.4 5.6 4.6 3.3 5.9 3.9 4.6 3.8
## 
## $func.thetastar
## [1] 0.0022
## 
## $jack.boot.val
##  [1]  0.46806723  0.36985075  0.25365169  0.20372493  0.06565934 -0.05425220
##  [7] -0.14710983 -0.26303725 -0.47041420 -0.52146893
## 
## $jack.boot.se
## [1] 0.9730894
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
##    [1] 6.7 5.3 5.3 4.3 2.7 4.0 4.2 5.2 4.2 4.2 2.9 3.7 4.5 3.5 3.9 3.4 4.9 4.1
##   [19] 3.6 3.5 5.2 5.0 2.9 3.3 3.2 4.4 5.2 4.1 6.2 3.2 4.0 3.8 5.6 4.9 3.7 3.9
##   [37] 4.9 5.8 6.4 4.8 4.8 5.3 4.1 4.1 4.2 5.1 4.5 4.4 4.7 3.8 4.6 5.6 4.3 4.4
##   [55] 4.8 4.3 3.4 3.4 4.1 4.8 6.5 6.7 3.8 5.2 3.9 4.1 3.7 5.5 3.8 6.0 5.6 4.0
##   [73] 4.2 5.9 4.9 4.2 4.8 4.1 4.6 4.0 3.7 3.7 4.3 4.9 3.2 5.2 4.8 5.3 3.4 5.0
##   [91] 2.9 3.5 3.7 4.4 5.1 4.0 3.2 5.4 3.5 4.1 3.4 6.1 3.4 4.9 6.0 5.3 4.9 4.4
##  [109] 4.1 5.7 4.8 3.4 2.9 6.3 5.0 5.0 4.7 3.2 4.2 3.5 5.5 5.0 4.3 5.9 4.7 3.4
##  [127] 5.3 5.2 4.7 2.5 3.9 5.0 3.3 4.1 3.8 4.7 3.9 5.7 2.7 2.4 4.3 2.4 4.9 3.6
##  [145] 2.4 6.3 5.0 4.1 4.3 6.4 6.3 5.2 3.2 5.1 4.7 4.5 5.7 3.4 6.7 4.6 4.6 4.1
##  [163] 4.3 4.8 4.0 5.7 4.5 5.0 5.3 5.5 5.5 5.4 5.4 4.5 5.1 4.5 5.4 5.1 5.4 5.5
##  [181] 4.3 3.4 5.0 4.0 4.6 6.4 3.9 4.4 3.1 4.6 4.5 5.2 5.2 4.2 5.3 4.7 5.8 5.5
##  [199] 5.4 4.6 6.0 5.6 4.4 5.3 4.8 4.9 4.6 5.7 5.0 4.4 3.8 5.6 4.3 3.8 3.5 4.7
##  [217] 5.8 5.2 4.7 4.8 4.3 5.3 5.6 4.4 2.7 4.5 4.3 4.1 5.9 2.5 4.8 4.7 5.2 4.8
##  [235] 5.8 4.3 5.6 5.9 5.3 3.4 3.6 3.7 2.8 2.6 5.8 4.4 5.3 3.5 6.2 5.8 3.7 4.7
##  [253] 4.4 4.8 4.6 5.3 5.5 3.2 4.6 4.7 5.4 3.7 4.5 4.9 3.7 3.7 5.9 3.7 3.7 4.2
##  [271] 3.4 3.7 3.5 4.9 4.4 5.7 5.3 5.5 4.4 6.7 3.0 5.1 3.0 3.4 3.3 4.1 5.4 5.4
##  [289] 5.8 4.7 5.2 4.4 4.6 3.5 5.4 4.7 6.2 5.0 3.7 4.5 5.0 4.9 1.9 3.6 4.3 4.9
##  [307] 3.2 5.2 2.7 5.0 4.5 4.4 5.0 5.4 3.3 5.5 2.7 4.4 6.0 4.7 3.2 4.3 3.8 6.0
##  [325] 5.0 4.6 3.9 4.9 5.5 4.6 3.2 3.5 5.2 4.8 4.1 4.4 5.2 4.6 4.0 5.1 4.5 4.7
##  [343] 5.0 4.5 4.2 3.6 4.7 4.6 4.4 2.4 2.8 3.6 5.6 5.3 3.6 3.6 3.9 2.7 5.2 5.7
##  [361] 4.8 4.8 4.8 4.2 4.0 3.0 2.6 3.5 4.9 3.7 3.7 3.7 4.0 5.6 4.0 3.6 4.5 5.3
##  [379] 5.8 6.3 3.3 3.6 6.3 4.4 3.9 4.3 4.3 4.6 6.3 2.5 5.1 3.6 4.8 3.1 3.9 4.8
##  [397] 4.5 4.9 5.3 4.7 4.2 2.1 3.7 3.7 4.7 6.4 4.2 4.4 4.0 5.0 4.0 2.8 5.4 5.7
##  [415] 5.2 3.6 4.1 4.4 5.3 4.3 5.1 4.8 4.7 2.5 4.5 4.6 4.6 4.1 3.4 6.1 4.1 4.1
##  [433] 4.5 4.1 3.5 4.8 4.6 4.8 2.4 5.2 3.7 1.8 4.6 5.6 6.1 4.4 4.6 4.0 5.6 3.9
##  [451] 5.1 4.4 4.8 5.3 4.3 4.1 4.1 4.6 4.9 5.4 5.1 4.8 4.8 4.8 3.8 6.3 3.4 3.3
##  [469] 5.1 5.8 5.2 5.5 5.4 4.6 2.9 5.5 3.2 5.8 6.4 3.1 4.8 3.4 4.3 5.3 5.7 2.9
##  [487] 3.8 5.8 5.1 5.6 6.1 4.4 3.9 5.3 6.0 4.1 6.1 3.1 4.2 4.4 3.9 2.9 4.5 4.6
##  [505] 3.8 4.0 6.4 6.0 4.5 4.4 2.8 5.3 2.8 5.4 3.4 6.2 4.9 2.9 3.7 5.0 5.7 4.0
##  [523] 4.8 4.3 6.1 5.7 4.5 3.4 4.8 4.0 3.7 5.8 4.9 3.8 3.3 4.5 6.5 4.8 5.0 4.7
##  [541] 6.0 6.1 3.9 4.4 3.9 5.4 3.9 3.5 4.1 4.4 4.6 4.4 3.8 3.9 4.4 5.1 5.1 4.0
##  [559] 3.5 4.9 5.4 3.9 3.9 4.1 3.9 3.6 4.5 4.2 4.6 5.7 4.8 4.5 5.8 3.8 5.0 5.2
##  [577] 4.9 4.9 5.8 2.7 5.7 3.7 4.3 4.9 3.2 4.4 3.8 4.4 5.5 6.7 5.6 5.1 4.3 4.6
##  [595] 4.7 4.5 6.2 4.0 4.8 5.4 4.5 3.9 4.5 4.3 4.5 3.7 5.8 5.1 3.7 4.5 4.8 4.7
##  [613] 3.5 4.2 5.5 4.4 4.4 3.8 3.7 4.5 6.2 3.0 4.9 5.0 3.9 4.8 5.5 1.6 3.0 4.5
##  [631] 4.5 3.7 3.2 6.4 5.8 3.3 4.7 3.5 4.2 3.9 3.3 5.2 3.5 6.0 4.2 3.7 4.9 5.4
##  [649] 5.5 5.2 3.9 4.8 4.0 4.6 4.0 4.8 4.0 3.8 4.4 3.5 4.1 4.4 2.1 4.3 4.1 5.1
##  [667] 6.4 4.3 3.3 3.8 4.2 2.9 3.4 5.4 3.8 5.9 4.2 4.9 3.0 5.7 4.4 5.1 4.2 3.7
##  [685] 4.2 5.4 2.8 4.4 3.8 4.7 4.9 5.6 3.5 3.6 4.7 5.9 5.6 5.4 4.8 4.0 4.0 3.0
##  [703] 4.8 3.0 5.0 4.5 3.7 5.7 5.2 4.5 4.1 4.2 5.3 3.7 3.5 3.9 5.5 4.9 6.0 2.7
##  [721] 4.5 5.3 4.2 4.0 4.7 4.7 3.9 2.9 3.8 5.1 5.1 3.2 4.3 5.5 5.6 2.9 4.1 4.0
##  [739] 4.6 3.1 3.8 3.2 4.7 5.8 4.4 4.5 4.2 4.6 5.8 5.7 4.1 4.7 6.0 5.6 5.4 4.6
##  [757] 6.3 5.2 2.8 3.9 3.4 3.0 3.4 4.5 3.2 5.6 4.7 5.0 4.5 5.2 3.8 6.1 4.7 2.5
##  [775] 4.0 3.0 5.4 3.8 4.9 3.2 3.5 5.7 3.8 5.4 4.6 4.9 3.9 4.0 5.1 5.9 6.1 5.5
##  [793] 3.6 4.9 3.7 4.3 4.0 4.6 4.3 5.7 3.8 5.8 5.9 5.4 5.5 4.9 5.3 5.6 4.3 5.1
##  [811] 3.8 5.0 4.3 3.5 5.2 4.7 4.4 4.8 4.5 5.2 5.6 5.4 2.9 4.2 4.0 4.8 3.5 5.7
##  [829] 4.5 4.2 3.2 2.8 4.7 4.7 3.5 4.1 5.2 4.5 6.0 4.9 5.2 4.7 5.2 3.6 3.3 6.1
##  [847] 4.7 4.9 5.2 4.3 3.1 4.0 2.4 4.7 4.9 2.8 4.4 4.3 4.1 5.5 3.4 4.0 4.0 4.7
##  [865] 4.6 3.1 5.5 3.4 5.4 5.1 4.5 4.5 5.6 2.5 4.0 6.2 4.6 4.6 5.6 4.6 5.6 5.7
##  [883] 4.4 4.3 5.6 3.9 4.3 3.9 4.7 5.4 5.5 4.4 4.8 4.3 4.7 4.0 4.2 4.5 3.4 4.9
##  [901] 4.8 5.0 4.0 4.2 3.4 5.2 6.1 4.5 5.0 2.9 3.9 3.1 4.3 3.7 6.3 6.6 4.2 3.0
##  [919] 4.5 4.9 3.5 5.4 4.7 5.7 5.0 4.8 5.0 4.7 4.8 3.7 3.1 6.3 4.1 2.6 4.2 4.6
##  [937] 2.9 4.3 4.2 3.2 4.0 3.6 4.0 5.1 4.9 5.4 4.9 4.9 3.8 4.3 5.5 3.8 6.2 3.9
##  [955] 4.4 5.2 3.7 3.3 3.9 3.7 6.7 5.0 3.6 3.5 4.5 6.3 5.9 5.3 3.2 5.0 3.6 3.6
##  [973] 4.9 6.8 3.6 3.2 3.1 4.5 5.1 2.9 3.7 4.5 4.8 4.0 6.0 3.6 4.5 4.5 6.2 3.4
##  [991] 4.6 5.3 3.7 4.9 3.2 3.7 3.0 3.9 5.1 3.4
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.6 5.5 5.3 5.3 5.2 4.9 5.1 4.8 4.6 4.5
## 
## $jack.boot.se
## [1] 1.054704
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
## [1] 0.8073545
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
##   2.3473963   5.0176914 
##  (0.9843425) (2.3451530)
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
## [1] 0.07525576 0.32216788 0.13865065 0.39542531 1.09050839 0.83128571
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
##    [1]  0.3320519960  0.8573506038  0.7424228322  2.2238056279  0.8147858286
##    [6] -0.8445698773  0.7449661966  0.8802532007  0.8346052820  0.7123460076
##   [11]  0.4075143051  0.7987594046 -0.4065329638  0.2968230105  0.8452136700
##   [16]  0.3770451038  0.7401300210  1.4413154530  0.3542926953 -0.0599964576
##   [21]  1.9082602964  0.6861178596  2.0048484986  1.8257508913  0.4009742582
##   [26]  0.0173564315 -0.0334741797  0.4074428300  0.2930341408  0.5926290608
##   [31]  0.8241272520  1.2328044651  1.3629476765  1.2833583170  0.0336159434
##   [36]  0.8361910376  1.1496545757  1.2285391318  2.2192790819  1.4027073662
##   [41]  1.2733996384  0.7180861766  1.0096785799 -0.3724375961  0.8698145674
##   [46]  0.7565741099  1.4833973255  0.8382161825  0.0558662333  0.8722834957
##   [51]  0.8063076044 -0.3421582480  0.6540769035  0.8146671100  0.7922893369
##   [56] -0.3824918711  1.2862084250  1.1165536335  0.3276568966 -0.0221187376
##   [61] -0.0724711281  0.6184935937  1.3629476765  0.7137435106  0.8119492204
##   [66]  1.1401751874  2.2853952433  0.6720938487  2.1582420208  0.2626455095
##   [71] -0.0052654444  1.2151439705  2.3608738615  0.0363484632  0.3136067385
##   [76]  0.7269013600  1.3952597258  1.5043375098  1.8635896378  0.7979113559
##   [81]  0.8201644843  0.8650619313 -0.0020773599  0.4177007172  0.7949813538
##   [86]  0.4280758942  0.7572616588 -0.4222003144  0.8148684503  0.3786609274
##   [91]  0.0233672264  1.2391947397  0.3154074140  1.3102707984  1.4016569068
##   [96]  1.4300231733  0.3720744333  0.4357502294  2.3847658293  1.3252791837
##  [101]  0.7848563712  1.2285917185  0.7859181207  1.4284928207  0.3497058268
##  [106]  2.2315853843  0.7794046724  1.4155384836  0.0317073142  0.8665629620
##  [111]  1.2142898443  1.2929689488  0.7976972591  0.4433096715  1.2438528080
##  [116]  0.8437082333  0.0319460069  0.7630026268  1.5187577854  0.7024214624
##  [121]  1.2289484367  0.8057196124  0.8709691747  0.3793080413  0.0231079257
##  [126]  0.4282520813  1.4731537877  0.8560975893  1.4387631961  2.3572924102
##  [131]  0.3643269037  0.8148684503  0.7025048637  1.3719210841  0.3235934938
##  [136]  0.4255810134  0.3969037593  0.7843680676  0.6554571019  0.3022469595
##  [141]  1.3287115665  1.3795284535  0.3912377010  0.3867017625  0.0532529734
##  [146]  0.6528773853  2.3095786431  0.8000066387  0.3452237666  0.4123671782
##  [151]  0.8470755048  0.8366454085  1.4466005882  0.7973194444 -0.0065444946
##  [156]  1.3055888856  1.2969000432  1.3722860079  0.2676617243  0.0187694363
##  [161]  0.3723592608  0.8563965534 -1.4932510801 -0.0119627391  0.3044626958
##  [166]  0.7807471712  0.3533559439  1.2282521429  1.4392404307 -0.0059725011
##  [171]  0.7808251552  0.4076563888  1.0457872042  2.3771875467  0.3765309269
##  [176]  0.8519335995  0.4292429585  0.4234052518  0.7682726830  0.3347248583
##  [181]  0.4129019886  1.2744072956  0.3936040229  2.2937214134  1.4899850230
##  [186]  1.1013823629  0.7876947319  0.2680669246  0.4539834357  0.2081645354
##  [191]  0.8063845348  2.4757929925  1.3641002840  0.6864164407  1.2007211421
##  [196]  1.4524825064  0.8152465932  0.8981480586  1.4873109224  0.8420735456
##  [201]  1.4882591314  1.4178452161 -0.4205930894  0.7364224803  2.3250350364
##  [206]  0.3982585366  2.3597112022  1.4258057705 -0.3644189846  1.3303373724
##  [211]  2.4806027074  1.3436275832  0.3562444256  0.8725480245  0.8074288201
##  [216]  0.6343376371  1.4965801890  1.9652917511  1.2750915025  1.8590635208
##  [221]  0.4650639835  1.4053511171  0.0241270902  0.4349063370  0.3632753974
##  [226]  0.3887015809  0.8946966703  1.3256534860  0.4148416451  0.8311426282
##  [231]  0.8224760497  0.7145602664  0.4395545981  1.5079941334  0.2956858848
##  [236]  1.1978956964  2.2246945036  1.2882403030  0.3455156126  0.3828584419
##  [241]  1.2948033753  0.3540302598  2.3517857315  1.4388175586  2.3855863677
##  [246]  1.4921222147  0.8230069158  0.3158856223  1.3658941290  0.3386479352
##  [251]  1.2862279351  0.3348309652  1.3872624055  1.4892175072  0.7863039499
##  [256]  0.7257449298  0.7549147296  0.2783277576  1.2899639323  0.3405641256
##  [261]  1.4653626075  0.0177361201  2.1944847292  0.8191648290  1.1667408925
##  [266]  0.7734419289  0.2588479837  0.8022486282  0.3164054940  0.8594652256
##  [271] -0.0677204324  0.7871226797  1.1762774892  0.4370795578  1.4632331588
##  [276]  0.8619415416  1.4412301430  0.3308826896  0.8017205288  1.3396055802
##  [281]  0.6654245808  0.8608415178  1.3482498727  1.3787426870  1.0100792367
##  [286]  2.2152936445  0.4341034731 -0.0498523380  0.8130111725  0.8140713383
##  [291] -0.4182081169  0.8298724244  2.3233064042  1.3293667369  2.5708681312
##  [296]  0.4025280821  1.3390965776 -0.1245443651  0.7072384069  1.2782640383
##  [301]  0.8464071968  0.0046677994  1.4338245723  0.4002211938  1.1917122700
##  [306]  0.8142889827  0.4222980213  0.8636563663  0.8088845321  0.7897179045
##  [311]  1.3639112839  1.4076265447  0.3641350965  1.4374675539  0.9131872921
##  [316]  0.3879689680  0.8392908230  0.7345962653  1.4471096493  1.3252791837
##  [321]  0.3269059202  1.3657826102  0.8719083983  1.4753491806  0.7813829269
##  [326]  2.3249475067  2.3677545666  1.2129633369  1.4699650591  1.2911588005
##  [331]  1.3538148607  2.2393331067  2.0549927916  1.4944423015  0.3434139219
##  [336]  0.3594389688  0.8703567753  0.7918071933  1.4260931687  1.1504821834
##  [341]  0.3400792647 -0.8049883629  1.3341141978  0.0120587228  2.4226851458
##  [346] -0.0151694234  1.1715060013  0.3416962112  0.8487210901  0.3911858408
##  [351]  0.8380683502  1.1750577126  0.8768323432  0.4093011915  2.4419486439
##  [356]  2.5808510371  0.4471526073  1.2300810730  0.2756159255  0.4387767139
##  [361]  0.3333867085 -0.3869605501  2.0033221898  0.7470920977  0.3945620671
##  [366]  1.4760521524  0.3531914718  1.0999707164  0.8725151516  0.8504528441
##  [371]  0.3598820203  1.4042889691  0.8093881350  2.3095786431  2.2246945036
##  [376]  0.7698531312  2.5237326938  0.8478931233  0.8143166460  0.8499056876
##  [381]  2.1722477136  2.3517857315  1.4626021860  2.3688237919  0.8562573857
##  [386]  0.8079755803  2.1200748909  0.7432396032  1.3718611818  0.7423606764
##  [391]  0.8381975300  1.4282931171  0.8238923001  0.4352721892  0.4164625507
##  [396]  0.0346090891  0.3970001658  1.5030487832  0.7733870841  0.4169920275
##  [401]  2.0235544554  2.3514258452  2.6323193361  0.7874956160  1.2074407222
##  [406]  0.8001964590  1.0172011257 -0.3946725732  1.4504441065  1.3799221896
##  [411]  0.7871226797  0.8771308988  0.2738742234  2.2190849820  0.8541985544
##  [416]  0.3209503457  0.8093881350  1.3751856229  0.7573496793  2.3002461549
##  [421]  1.3895925875  0.2341956431  1.4121475714  1.2749605077 -0.0001524350
##  [426]  2.3178029608  1.4402464156  0.0167941415  2.4350376315  2.3427277378
##  [431] -0.0506161550  1.1510866292  0.3672540300  1.4726097356  0.0166365748
##  [436]  0.3649624083  0.5901156827  0.7174188420  1.4724073258  0.7854950690
##  [441]  0.8173580918  0.4511159580  0.2806654451  2.5305001349  0.4168014846
##  [446]  0.8161542710  0.8808885371  1.2015737762  0.8240236997  1.3559136894
##  [451]  0.0339966059  1.8587846747  0.8119492204  1.4395309199  1.5409866537
##  [456] -1.4247511327  0.4217805873  0.0420315537  0.3911027414  1.3239492059
##  [461] -0.0859183626  1.4821433999  0.7423606764  0.3160389061  0.0381994549
##  [466]  0.8506167983  0.4280571615  2.1439328934  1.2716813940  2.1748694139
##  [471]  1.3065326741  1.4024586380  0.7782503365 -0.0324793030  0.8803263477
##  [476]  0.4073933263  2.5272723072  0.8088867368  0.8778754276  1.2334558473
##  [481]  0.8207050878  0.8211938853  1.3820734405  0.2958789980 -0.3438297932
##  [486]  0.4265144669  1.4071676524  1.3289173429 -0.4162159537  2.2819848546
##  [491]  1.1855476098  0.8440576567  1.2956566083  1.3297866710  0.7680380763
##  [496]  1.3458000666  2.1306237713  1.3177814878  0.8734816346 -0.7970664118
##  [501]  0.4389917579  2.0750114391  1.2283385641 -0.1303044173  0.3502898366
##  [506]  1.3445036311  1.3702349033  0.7603543391  2.5138856888 -0.3951122495
##  [511]  0.3232628340  1.4777456005  0.8083094802  0.4475900112  0.4341346628
##  [516]  1.4093396551  1.3897221418  0.3916760851  0.3153882986  0.0241476688
##  [521]  0.8849154157  0.8316865229  2.2130186465  0.7614477430  0.8866149354
##  [526]  0.2981669179  0.7884004971  1.4183379990  2.2294811026  0.2326464846
##  [531]  2.5138856888  0.4038487726  0.4359555901  1.4429200133  1.5079693637
##  [536]  2.1484126018  0.3248016553  0.3741072317  0.3561222562  0.3268728554
##  [541]  0.8123992738  1.5238203351  0.4417461289  1.3126150291  1.5089808032
##  [546] -0.0149094290  1.2439048943  0.7447820031  0.4396565246 -0.4546374503
##  [551]  0.7451517861  0.3817364518  0.0121503804  0.7168812185  1.4588757542
##  [556]  0.4052118232 -0.3651367227  2.2646278911  1.3598163290  0.4015212244
##  [561]  0.4047148924  0.7262896159  0.7911288169  2.0918650297  2.5130196459
##  [566]  0.4287020283  1.4658428428 -0.0720523229  1.2726147301  1.3035342405
##  [571]  1.7818719314  1.2278830751  1.4445505776  0.8438078966  1.2995256487
##  [576]  0.8070158141  1.3702349033  0.8978827255  1.4855286062  1.1848050590
##  [581]  0.7094389804  0.8959884223  0.0092891166  0.4096013240  1.3672121082
##  [586]  0.3881068036  0.0900199189  2.2222220608  0.7813829269  0.4392459803
##  [591]  1.1917122700  0.8263634061  0.2905842438  0.4035553996  0.7383159220
##  [596]  0.7331258959  1.2685791729  0.3936188209  1.3870161413  0.8124892925
##  [601]  1.3037741714  1.4815565648  0.7909057532  0.7596200276  2.2954118865
##  [606] -0.0102754291  0.0102298258  0.7910157000  0.8379333321  0.3056524863
##  [611]  2.0132402965  1.4713510784  0.0681894974  1.3867379526  0.0502909610
##  [616]  0.0154957087  1.2423170123  0.3143930137  1.2590482914 -0.3442381989
##  [621]  0.7909057532 -0.0282429264  1.2391947397  2.1489727080  0.4397266334
##  [626]  0.3795438753  0.3667745354  0.8844378522  1.4769793825  0.3715363447
##  [631]  1.4416496777  1.3766022574  0.7752296224  0.8870002755 -0.3890574763
##  [636]  2.2067324646  0.8132914481  0.3888793353  0.8846042110  0.7693923479
##  [641]  0.7278242265  2.5347404558  1.2169764230 -0.0934360658  0.2869446880
##  [646]  1.3770319586  0.7386770770  0.3515897059  1.4965968649  0.8725699957
##  [651]  0.7964572178  2.3713455100  0.7212023504  0.4176914532 -0.3115368698
##  [656]  2.2713430800  0.8024088278  1.3559136894  0.3379529266  2.2965595179
##  [661]  1.4203637410  0.4120172780  1.3832663950  0.0539006471  0.8525032131
##  [666]  0.4027866292  0.8522148991  0.7653662829  1.0973366984  0.3283223162
##  [671]  1.3921120367  0.8091007763  1.8173798706  1.5452234112  0.8657901553
##  [676]  1.2624699651  0.7868138236  0.9052672460  1.4990701855  0.3730170111
##  [681]  0.3415406814 -0.0013522108  0.4495674247  0.3238170292  0.4169865143
##  [686]  0.4515722189  2.3021618710  1.9038663063  0.8787870586  0.3432791875
##  [691]  2.1537524924  0.8970451766  1.2027329175  0.7620537144  0.4150672427
##  [696]  0.4245250534  0.8978336563  0.4199174883  0.9001397801  1.4418965969
##  [701]  1.1440503813  0.7877574784  2.3332385736  2.4864457069  0.4271596412
##  [706] -0.0289533509 -0.0046729629  0.9120884443  0.8784040413  1.6667178814
##  [711]  0.0019809839  0.6621050784  0.8430870234  1.3563383428  1.3430785894
##  [716]  0.9022013178  1.8738718540  1.3113534749  1.3437512217  0.7521861003
##  [721]  1.2081350523  0.9571185540  0.3735693253  0.7455482039 -0.4007944041
##  [726]  0.3186006561  2.2168181572  1.3285956963  0.4305544162  1.3810849967
##  [731]  1.4574370104  1.2130964055 -0.0328617192  0.0236225150  0.3341475763
##  [736]  0.8487952995  2.1200037903  0.7185362116  0.0220816907 -0.0685950626
##  [741]  1.4643433906  0.4454483390  0.0467930390  0.4236404307  0.3341475763
##  [746]  0.0271580401  2.4935760695  1.4950576738  0.7816088722  0.7741488637
##  [751]  0.3735392627  0.8408233228  0.3542340815  0.7850331554  1.3634809717
##  [756]  0.3068819986  0.4125369815  1.1854446892 -0.4427498784  0.3448481755
##  [761]  2.5330530326  0.4381025821  0.7425733440  1.4688587446  1.4526311368
##  [766]  0.6978336973  1.5034947946  0.8200184412  0.7952899751  0.0020155566
##  [771]  0.3518385202  0.3769156077  1.5159776018  1.3607092772  2.1471882817
##  [776]  1.4669599390  0.6560160728  2.4857251698  0.8891152078  1.5144985834
##  [781]  1.4294426844  0.8927986009 -0.0099324604  0.0006363645  0.3697295998
##  [786]  1.2263062070  1.5126935118  1.3037741714 -0.0333970613  0.7875656525
##  [791] -0.8535125053 -0.4556729183  1.3034717738  2.4729344568 -0.0733899897
##  [796]  0.4413246695  1.2685078673  1.4472120040  1.3720094210  0.8884759533
##  [801]  2.0180406480 -0.8544183125  0.7994929906  0.0576587749  0.3212501063
##  [806]  0.4466712833  0.4076479437  0.8581233203  0.8214350929  0.8404428878
##  [811] -0.0282429264  0.3626006680  2.1969797434  0.4238109466  1.4293344078
##  [816]  1.2831419231  0.7336915335  0.8073544992  0.4367596211  1.4844300318
##  [821]  0.7100791946  2.2879007817 -0.0425072810  0.7650358439  0.6116392176
##  [826]  1.8507967549  0.3569256329  1.4069312724  1.9878746671  1.2263062070
##  [831]  0.3472282861  0.8075044678  0.8236437718 -0.0150427251  0.2853187967
##  [836]  0.0660472010  0.8033015838  0.8833478047  0.3899886647  0.4341986693
##  [841] -0.1051536574  1.2307152962  0.8068466684  2.5263888152  2.2828974316
##  [846]  1.2824278646  1.9784330906  0.8484201636  1.4182226943  1.4554820683
##  [851]  0.3034231556  1.3086983061  0.8238571584  0.9024469388  0.0380506142
##  [856]  1.3778545068  0.4187147664  2.3549129239  0.8648161561  0.0576587749
##  [861]  2.0759712621  1.5234339058  1.3795986994 -0.0087530876  2.2214208074
##  [866]  0.0187082472  1.4143563333  0.4345564817  0.8057734497 -0.0657716497
##  [871]  0.3902257557  0.7546611898  1.3913832473  0.7460452318  1.3024181198
##  [876]  2.5588198566  0.7515341064  0.7540987826  2.0532079027  2.2262183619
##  [881]  1.4044092322  0.3713630741  1.4659560324  0.8589183942  0.3242497591
##  [886]  0.0107995020  0.7834414335 -0.0772343017  0.7681643266  0.7450927475
##  [891]  0.7266413392  0.4213076909  2.3327855073  2.5577640621 -0.3430810470
##  [896]  0.8083104310  1.3759863293  0.7891306967  0.3671136995  2.0396279079
##  [901]  0.7780655296  1.3091446315  0.0478527923  2.1507715401  1.3691469239
##  [906]  1.2737215624  2.1489727080  1.2263817260  0.6649881017  0.7933257706
##  [911]  1.2726147301  0.8926011457  0.8489872259  0.3572158773  0.4035004398
##  [916]  0.4161273849  2.2430234643  0.0002758712  1.9803716281 -0.0059813936
##  [921]  0.4025280821  0.8024181397  0.8422887588  0.3972346702  0.8389499515
##  [926] -0.1551870825  1.1441795471  1.4934841999  0.4680543568  1.3177814878
##  [931]  2.4993426781  0.0095111494  0.7852374244 -0.0184614162  0.8265381778
##  [936]  0.8637555096  0.4079568740  1.2810575113  0.7946311449  1.3373066265
##  [941]  1.3460227325  0.8596808762  0.4373696839  1.1643709416  1.5058546706
##  [946]  0.8528999683  0.8685415585  1.3289020409 -0.4224897702  0.7875827963
##  [951]  2.2141273975  0.3273155011  0.7170006404  2.2173431566  0.7165083955
##  [956]  1.3225077226  0.8375248458 -0.0277476831  0.4290699737  1.3774535352
##  [961]  0.4233962581  0.8602987713  0.8552741501  0.2341111235  0.4210794175
##  [966]  0.2641680924 -0.0248822870  0.3485963188  0.3358715166 -0.3599826175
##  [971] -0.4857538527  0.3335502326  1.4784586474  0.4186636660 -0.3057795325
##  [976]  0.8062939090  0.3583247364  2.1824570163  1.2323708092  0.8361984148
##  [981]  1.1946496981  2.1790775637  2.5479970242 -0.4245331546  0.0328807925
##  [986]  0.8532356508  1.4488746820  0.9085969499 -0.0067898430 -0.0148019303
##  [991]  0.0001943736  0.3772317434  0.8489113328  0.3621306178 -0.0150198885
##  [996]  0.8439835211  2.2919910213  1.1515923739  2.5544482362  0.8280299529
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
```

```r
fit2
```

```
##      mean         sd    
##   0.4678238   0.3277465 
##  (0.1036425) (0.0732837)
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
## [1] -0.31826751  0.15828036  0.09302347 -0.56375736  0.06139986  1.01480553
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
## [1] 0.0478
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8955765
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
## t1*      4.5 -0.02092092   0.9087709
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 1 2 4 6 9 
## 2 1 1 3 1 2
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
## [1] -0.0272
```

```r
se.boot
```

```
## [1] 0.9127614
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

