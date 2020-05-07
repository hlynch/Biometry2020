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
## 0 1 2 8 9 
## 2 1 2 2 3
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
## [1] 0.0411
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
## [1] 2.723301
```

```r
UL.boot
```

```
## [1] 6.358899
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.7   6.2
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
##    [1] 4.2 4.8 4.1 5.6 4.6 3.9 4.0 3.7 4.1 5.7 5.7 4.1 3.1 6.1 3.0 3.8 4.4 4.2
##   [19] 3.8 4.3 5.6 5.7 3.2 3.5 4.6 6.0 3.0 4.5 4.1 5.6 3.4 3.5 3.8 3.3 5.1 5.1
##   [37] 3.3 6.1 4.1 4.0 5.0 5.4 4.5 4.6 5.9 4.8 3.4 5.8 2.9 4.9 2.9 5.8 4.6 4.0
##   [55] 3.8 3.2 3.1 4.0 6.0 5.1 4.5 4.6 6.3 5.5 6.1 2.3 4.9 4.1 3.4 4.5 5.1 5.1
##   [73] 4.9 4.6 4.8 5.8 5.9 4.6 4.4 4.1 4.4 4.4 3.8 5.5 3.0 4.6 4.5 3.9 3.1 3.7
##   [91] 4.9 3.3 4.0 4.1 4.7 4.9 5.3 6.1 4.5 4.0 5.9 5.3 4.4 4.1 3.3 4.9 4.1 4.0
##  [109] 3.1 4.9 4.6 2.6 2.9 6.3 5.5 6.4 4.6 4.1 4.3 4.6 4.1 3.2 3.4 4.4 4.7 4.5
##  [127] 4.6 4.2 4.0 4.6 5.2 4.5 3.5 5.1 4.0 5.6 4.0 5.6 4.8 5.4 4.1 4.2 2.9 5.5
##  [145] 7.0 4.3 3.9 5.8 3.2 3.9 5.3 4.2 4.5 4.1 5.3 4.2 5.0 5.1 5.1 4.0 3.4 4.9
##  [163] 4.7 4.9 4.2 5.4 4.8 4.2 5.1 5.0 5.0 5.3 4.5 5.3 2.2 3.3 4.9 5.6 5.8 2.8
##  [181] 4.0 4.7 5.1 5.4 4.9 4.7 4.5 3.8 4.7 4.9 3.0 4.6 4.8 3.1 6.1 3.5 6.0 4.1
##  [199] 3.7 4.9 5.3 3.3 4.2 3.6 5.7 3.8 5.5 3.1 4.2 4.7 2.9 4.3 4.8 4.8 4.7 4.8
##  [217] 3.0 4.1 6.2 4.6 4.6 4.4 4.3 3.7 5.7 5.1 5.3 4.8 4.0 3.7 4.6 3.9 5.0 5.1
##  [235] 4.9 4.1 4.2 4.0 4.8 4.6 5.7 3.0 4.2 6.8 4.9 4.1 4.3 3.8 5.1 5.3 5.1 4.6
##  [253] 4.2 3.1 4.1 4.6 5.2 5.0 3.9 5.0 4.4 3.3 2.7 4.1 4.1 5.6 4.1 4.4 4.3 4.1
##  [271] 6.2 5.6 3.9 4.6 2.8 3.0 2.0 6.6 5.8 4.0 5.2 3.9 4.9 5.2 3.7 3.6 5.5 4.5
##  [289] 5.2 5.1 5.5 4.1 5.5 3.4 4.6 5.4 5.6 4.8 5.2 4.6 3.7 5.5 4.1 6.0 3.4 5.0
##  [307] 5.2 3.2 3.3 3.5 5.7 5.2 3.5 5.0 4.3 6.2 5.2 3.9 4.9 5.1 4.3 3.6 5.2 4.9
##  [325] 5.9 3.0 5.3 3.5 3.2 3.1 3.7 3.4 6.5 4.0 2.4 5.6 2.3 4.7 4.5 4.6 3.1 5.2
##  [343] 4.4 4.9 4.9 6.0 4.4 4.0 5.0 3.3 3.6 3.8 4.4 3.7 3.5 5.1 5.7 5.0 3.3 4.5
##  [361] 3.8 3.6 5.0 5.2 4.9 4.8 4.8 4.5 3.0 5.0 4.3 4.2 5.6 3.3 4.6 3.5 3.8 3.6
##  [379] 5.6 4.5 6.8 4.4 4.3 4.7 5.3 3.9 6.4 4.9 5.7 4.8 5.4 4.3 3.7 5.4 6.0 3.5
##  [397] 3.3 3.5 3.0 4.4 4.4 4.7 4.0 4.5 3.8 5.4 5.6 2.7 4.3 4.9 3.2 5.9 2.4 6.0
##  [415] 3.7 6.1 4.3 2.3 3.8 5.4 3.1 4.7 5.0 4.8 4.0 5.5 4.3 5.3 3.7 4.5 5.0 4.7
##  [433] 6.3 5.3 5.1 3.9 5.7 5.1 5.5 3.5 4.8 3.5 5.5 3.6 4.0 4.3 3.8 4.2 4.4 5.2
##  [451] 5.2 4.7 3.7 4.3 4.5 3.7 5.2 5.6 6.7 3.8 3.2 4.7 5.5 5.2 4.5 5.9 4.9 3.2
##  [469] 4.3 3.3 5.6 7.1 6.6 4.4 4.6 3.5 3.7 3.7 4.1 5.2 5.3 3.2 4.9 3.2 5.1 5.9
##  [487] 1.5 4.8 3.7 6.0 3.9 3.9 4.1 4.6 4.8 4.0 5.3 3.9 4.8 5.2 3.1 4.9 5.9 2.9
##  [505] 4.5 5.0 3.5 4.9 5.6 3.9 6.1 3.4 5.0 4.0 4.5 3.9 4.7 5.6 4.2 5.1 5.2 4.6
##  [523] 4.8 6.2 3.3 4.7 2.6 4.5 4.2 4.5 3.0 5.5 4.6 3.4 3.9 4.5 4.8 4.2 6.4 4.0
##  [541] 5.6 4.8 4.1 4.1 4.4 4.9 3.6 5.2 3.3 5.7 3.1 3.8 3.6 5.4 4.5 3.8 4.5 3.2
##  [559] 4.5 4.1 3.8 3.4 6.1 4.9 5.8 2.6 4.8 2.2 4.7 4.5 4.8 4.1 5.0 4.8 7.2 6.5
##  [577] 5.3 5.5 5.6 4.9 4.6 4.2 3.9 4.5 3.5 3.2 3.8 3.9 4.9 5.8 4.0 6.1 4.3 5.4
##  [595] 4.3 5.3 4.3 4.4 4.5 4.1 5.3 2.8 3.0 4.7 5.1 4.7 6.1 2.8 6.5 5.9 4.8 5.5
##  [613] 4.5 5.8 3.0 4.2 4.3 3.7 3.7 5.0 3.2 5.5 3.1 3.3 4.3 4.2 4.2 2.9 5.3 3.1
##  [631] 5.8 4.3 4.0 4.7 7.2 6.8 3.6 3.5 5.2 6.1 5.4 5.7 5.3 5.3 4.5 4.6 4.6 4.6
##  [649] 4.8 3.9 5.0 4.4 3.7 5.9 5.7 4.1 3.9 4.5 4.4 5.0 5.3 3.2 6.8 5.7 6.0 6.3
##  [667] 3.7 4.4 2.6 3.7 3.8 3.7 5.1 4.8 5.0 3.5 4.8 4.1 4.7 4.1 2.6 3.8 3.6 4.2
##  [685] 5.3 3.5 4.3 5.1 4.9 4.1 4.0 4.2 4.5 5.1 4.8 5.3 5.3 3.2 6.2 3.4 5.4 3.4
##  [703] 5.3 5.4 5.3 4.0 4.0 4.1 3.8 2.9 4.0 4.2 3.3 5.8 3.7 5.7 4.1 4.1 4.6 6.1
##  [721] 5.9 5.5 4.1 4.5 4.0 6.1 4.8 4.5 4.4 4.7 2.4 5.8 4.2 5.0 4.4 4.4 4.6 2.5
##  [739] 4.9 5.6 5.0 3.6 4.5 4.5 5.6 4.6 4.5 3.4 4.8 4.9 5.6 4.7 4.3 4.9 5.2 6.8
##  [757] 3.9 4.2 4.5 3.8 5.0 5.5 3.9 4.0 4.6 3.5 3.8 4.4 4.6 5.1 5.4 3.8 4.0 4.7
##  [775] 4.2 5.1 4.9 5.0 3.3 3.5 4.6 4.5 5.6 3.3 4.4 5.2 4.7 5.4 6.6 5.2 2.4 4.6
##  [793] 6.0 6.8 4.4 3.7 5.4 3.7 5.3 3.6 4.8 2.8 5.0 4.1 5.0 4.9 4.6 2.8 3.9 4.1
##  [811] 4.8 5.1 5.9 5.6 6.7 4.2 5.8 4.8 6.0 4.4 3.9 4.7 5.0 3.9 4.7 5.3 2.3 5.9
##  [829] 4.2 5.1 4.7 3.9 4.1 5.8 3.5 4.5 4.2 4.9 4.5 4.1 4.8 3.5 4.4 4.0 4.9 4.1
##  [847] 3.4 4.1 4.5 5.5 5.3 5.4 4.4 4.7 3.1 4.9 4.4 3.7 3.4 2.5 5.4 3.7 4.8 2.3
##  [865] 4.0 3.7 4.7 4.6 4.2 4.2 4.5 4.6 4.1 6.5 3.8 4.2 5.1 5.4 3.5 3.6 3.5 4.5
##  [883] 3.7 4.1 5.3 3.2 4.4 3.6 3.6 5.0 4.0 3.7 3.8 4.8 3.9 4.5 4.1 4.5 5.0 3.8
##  [901] 4.5 4.8 4.4 5.1 4.4 3.7 4.2 6.1 5.9 4.4 3.6 3.8 6.0 4.6 5.8 3.1 6.4 5.1
##  [919] 4.8 6.1 5.4 2.9 4.6 4.6 4.4 4.5 4.2 5.1 4.3 5.9 3.8 5.9 4.1 4.3 3.8 3.8
##  [937] 3.8 5.3 3.6 4.1 5.1 5.6 4.3 3.7 5.5 5.5 3.9 4.8 5.0 3.3 4.8 4.8 4.4 3.3
##  [955] 3.4 4.0 4.8 5.1 4.2 4.7 3.5 5.2 3.5 4.0 5.9 3.9 4.7 3.5 3.1 5.1 4.8 5.6
##  [973] 4.2 4.0 5.0 3.5 4.5 4.6 5.1 4.5 3.8 2.2 4.6 4.1 6.0 3.7 4.6 4.0 3.7 4.3
##  [991] 4.2 4.7 3.8 4.6 4.0 4.3 3.7 4.4 4.8 4.2
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
##   2.8   6.3
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
##    [1] 4.9 5.8 5.6 3.7 4.9 3.3 5.3 5.1 4.0 5.5 4.2 3.6 4.2 4.6 3.7 4.4 5.1 4.2
##   [19] 4.8 4.9 4.5 5.7 3.9 4.4 5.1 5.4 3.7 5.2 4.7 4.4 5.9 4.4 5.0 4.7 3.7 4.4
##   [37] 5.3 5.4 4.3 3.6 2.7 5.1 5.5 4.6 5.3 4.2 3.0 4.4 3.9 5.2 5.9 4.1 4.1 3.5
##   [55] 3.5 3.1 2.4 4.0 5.9 4.7 4.9 4.9 4.5 4.8 4.3 5.1 5.0 3.1 3.7 4.7 4.8 3.8
##   [73] 4.6 4.6 5.1 3.4 3.9 5.4 4.8 3.9 3.3 4.8 4.8 3.7 3.6 3.3 4.3 3.6 3.4 4.3
##   [91] 5.6 4.0 4.3 2.2 3.7 5.7 3.9 4.8 3.8 5.2 3.7 5.0 3.6 6.2 4.2 4.6 4.3 4.6
##  [109] 4.3 5.2 6.0 4.4 4.9 4.6 3.2 5.2 4.1 4.8 4.5 4.6 6.4 5.3 3.9 3.7 2.9 4.5
##  [127] 4.8 3.2 4.8 3.9 5.5 5.9 4.9 4.0 6.3 5.2 5.3 4.9 3.5 3.6 5.1 3.2 4.9 3.9
##  [145] 2.4 4.1 4.8 6.1 3.5 3.3 3.7 5.2 3.6 3.8 3.6 3.1 3.6 4.9 4.6 4.4 4.6 4.8
##  [163] 4.9 3.5 4.9 2.8 4.8 5.2 4.6 3.1 3.5 4.5 4.0 4.6 4.1 4.2 2.2 5.4 5.8 4.9
##  [181] 4.7 2.9 3.9 4.7 5.9 3.6 4.5 3.8 3.2 5.5 4.0 3.3 4.3 5.4 5.2 5.0 3.9 5.0
##  [199] 5.9 4.8 3.2 4.3 4.6 5.4 4.3 5.3 4.5 3.4 5.8 5.3 4.7 4.8 4.6 5.8 3.7 3.4
##  [217] 4.0 3.9 2.8 4.1 3.6 4.9 3.4 3.8 5.3 1.7 4.3 5.1 5.5 5.2 4.3 5.7 5.2 4.6
##  [235] 4.3 4.6 3.9 4.3 3.2 3.9 3.2 3.0 4.0 4.2 2.6 4.6 3.3 3.9 4.4 4.7 3.5 3.3
##  [253] 3.2 5.2 4.1 4.7 5.1 3.2 2.7 6.1 3.9 4.4 4.5 5.9 5.9 5.6 3.2 3.9 4.8 4.3
##  [271] 5.0 3.5 5.7 3.9 5.9 3.7 3.4 3.8 4.8 5.4 2.3 5.0 3.7 5.8 4.2 4.8 3.6 4.1
##  [289] 4.7 4.9 5.1 4.5 5.8 5.8 3.0 4.1 3.7 3.3 5.5 3.3 3.1 3.7 4.0 4.6 4.5 5.4
##  [307] 4.7 4.1 4.2 5.9 4.9 5.1 3.7 5.7 5.5 4.6 4.1 5.5 5.2 5.6 3.5 4.7 5.0 3.8
##  [325] 4.9 3.7 4.3 4.3 6.1 5.3 4.2 4.8 5.2 5.0 4.9 4.4 5.3 3.7 4.2 5.1 4.8 4.0
##  [343] 4.2 4.9 3.1 2.9 5.7 6.2 2.9 3.8 6.1 4.0 4.1 4.7 4.5 6.1 4.7 4.6 3.7 5.3
##  [361] 2.8 3.9 3.8 3.9 4.0 3.6 5.2 5.0 4.7 4.8 4.5 4.5 5.1 5.8 5.3 5.5 4.8 4.6
##  [379] 4.4 4.1 4.5 4.5 4.4 5.2 3.3 4.7 4.8 5.0 4.1 5.2 3.2 4.8 5.9 5.1 5.1 3.8
##  [397] 4.5 4.6 4.7 4.2 3.7 4.5 4.7 4.7 4.1 3.6 3.9 2.5 4.1 4.1 4.7 4.3 5.3 4.4
##  [415] 4.9 5.5 3.4 4.7 4.3 6.0 4.7 3.0 3.3 5.7 3.3 4.1 3.2 4.8 4.2 4.6 4.7 4.1
##  [433] 2.8 3.9 5.1 4.8 3.4 4.3 4.0 3.5 5.0 3.8 5.1 3.9 4.0 5.5 5.7 3.0 5.8 5.0
##  [451] 5.2 4.8 7.1 5.1 4.2 4.4 6.2 3.7 5.7 5.3 4.4 5.5 5.4 4.8 7.0 6.2 2.9 3.5
##  [469] 5.4 7.4 5.9 4.0 3.0 4.7 3.5 4.4 5.8 3.4 4.3 4.6 2.4 4.2 4.2 5.2 4.5 5.3
##  [487] 4.6 4.8 4.5 5.2 3.9 4.8 4.9 4.7 4.9 5.9 2.6 5.3 3.8 4.3 4.8 4.9 3.5 4.8
##  [505] 3.0 6.0 6.0 4.9 3.3 5.1 6.1 4.0 4.0 4.0 3.4 4.7 4.3 6.6 4.2 5.3 5.1 4.0
##  [523] 3.6 4.0 3.5 3.6 4.6 4.7 4.2 3.8 4.0 5.0 4.3 4.5 4.0 6.0 5.0 2.6 2.6 3.5
##  [541] 5.0 3.3 3.3 4.1 5.0 5.1 3.7 3.0 4.3 5.5 4.0 3.3 3.8 4.7 5.9 5.8 5.2 4.6
##  [559] 5.2 1.9 5.1 3.3 5.3 5.7 4.5 3.9 5.8 5.4 2.8 3.4 5.1 3.5 3.4 5.2 4.6 4.3
##  [577] 4.9 3.6 5.4 6.4 5.2 5.4 5.4 2.9 4.7 4.1 4.7 3.9 4.2 2.3 2.9 3.3 3.5 5.5
##  [595] 3.3 7.2 5.0 4.8 3.4 3.6 4.4 5.1 4.0 5.4 4.4 5.5 5.6 4.0 4.5 3.6 4.4 5.5
##  [613] 4.6 3.5 5.5 6.0 3.1 4.3 4.9 6.4 3.5 4.5 3.5 4.7 2.3 5.7 4.8 3.4 4.8 3.9
##  [631] 4.2 5.4 4.7 4.8 4.1 4.8 4.2 3.9 5.1 4.0 5.0 4.9 3.2 5.6 5.1 4.8 4.7 3.5
##  [649] 4.8 3.7 6.0 4.6 4.4 3.6 4.4 5.0 3.8 4.4 5.9 4.3 4.3 4.5 3.8 6.1 6.6 5.6
##  [667] 5.8 5.2 6.0 4.4 4.2 4.3 3.6 3.2 3.5 4.6 4.2 4.6 4.0 4.1 3.9 4.7 5.8 4.1
##  [685] 3.9 2.7 4.7 4.6 3.4 5.0 3.0 5.3 3.9 3.3 3.4 5.3 3.3 4.0 5.4 3.5 5.3 4.4
##  [703] 5.3 3.5 5.0 3.8 4.4 4.3 4.2 4.8 4.3 4.9 3.8 5.5 5.4 5.4 5.0 4.5 5.2 5.1
##  [721] 3.8 2.3 4.1 5.7 5.8 5.2 3.6 4.6 5.7 5.6 5.6 5.7 5.5 5.0 4.8 5.8 2.5 5.1
##  [739] 4.0 5.4 5.7 5.1 4.8 3.2 3.7 5.0 4.1 4.4 4.5 6.0 3.9 5.5 3.8 5.3 4.7 5.5
##  [757] 4.5 5.0 4.8 3.9 4.3 2.7 4.6 4.3 5.3 3.6 2.8 4.8 5.8 4.5 3.2 4.0 5.6 4.9
##  [775] 6.3 4.3 5.2 2.6 3.4 3.8 4.9 4.9 4.1 5.3 4.8 3.7 4.1 5.5 4.7 5.4 4.6 4.3
##  [793] 3.2 3.2 2.5 4.4 3.2 5.9 5.5 3.9 3.9 3.2 5.5 4.6 5.0 5.1 4.6 4.0 3.0 4.7
##  [811] 4.5 3.9 4.3 3.8 5.2 5.8 3.7 5.2 4.9 4.1 4.9 5.2 3.3 4.6 4.6 5.4 3.5 4.9
##  [829] 2.9 5.9 4.4 3.8 3.4 4.3 3.1 5.7 5.6 4.7 4.4 5.1 3.5 4.7 4.7 5.7 3.0 5.5
##  [847] 4.2 4.2 3.7 4.9 4.2 5.2 3.4 2.8 4.9 2.7 3.9 4.6 6.4 2.8 4.4 4.6 4.4 4.7
##  [865] 5.0 5.9 4.5 4.3 5.0 5.9 3.8 4.5 5.2 4.3 5.0 3.9 4.9 3.2 5.8 3.4 4.8 4.4
##  [883] 4.5 5.1 5.0 3.4 4.7 4.2 5.7 6.2 3.3 4.3 4.9 4.1 4.6 4.3 4.4 4.2 4.6 4.9
##  [901] 5.9 3.8 3.6 4.8 4.1 4.7 3.9 3.3 6.0 4.7 4.1 3.6 3.4 4.6 5.1 4.4 4.5 3.5
##  [919] 5.2 4.9 4.3 2.7 4.2 6.0 2.5 4.5 4.8 4.0 2.6 4.6 4.6 3.5 4.9 3.3 3.8 4.9
##  [937] 5.6 5.6 5.6 2.5 4.3 6.0 5.0 3.8 4.1 5.4 4.4 5.6 2.9 5.2 4.3 3.8 3.5 6.0
##  [955] 4.5 6.1 5.3 4.1 3.2 4.4 3.6 4.1 4.1 4.8 4.6 7.0 5.0 3.4 5.1 4.4 4.6 5.2
##  [973] 3.3 4.8 3.4 5.5 2.1 6.4 3.7 4.8 5.4 3.9 5.1 4.5 5.2 4.7 4.1 3.9 4.9 4.1
##  [991] 5.4 3.1 3.9 4.3 3.4 4.9 5.0 5.4 5.8 5.0
## 
## $func.thetastar
## [1] -0.0316
## 
## $jack.boot.val
##  [1]  0.45866261  0.30720461  0.35718157  0.13343195 -0.01988473 -0.15683060
##  [7] -0.16006006 -0.32080925 -0.38685714 -0.52471910
## 
## $jack.boot.se
## [1] 0.9575921
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
##    [1] 5.6 4.3 3.7 5.7 5.1 2.6 4.7 5.0 4.6 5.8 4.0 4.6 4.7 4.6 4.5 5.5 4.5 5.3
##   [19] 5.7 4.7 4.5 4.0 4.7 3.1 5.0 4.6 5.1 5.9 3.5 2.5 5.3 4.9 4.8 4.4 4.7 5.8
##   [37] 4.6 4.5 6.8 4.4 4.5 3.7 4.5 4.5 4.7 2.8 3.8 4.3 4.5 4.8 3.9 2.5 3.7 5.0
##   [55] 4.2 4.3 3.5 3.6 4.6 7.1 3.7 6.3 3.7 5.0 5.3 6.0 4.5 4.4 5.2 6.0 5.4 5.2
##   [73] 4.1 5.0 5.2 6.5 4.4 3.8 5.9 3.6 3.0 5.2 5.4 5.7 4.5 3.4 5.0 5.7 4.8 5.1
##   [91] 2.9 5.2 3.1 4.0 4.6 5.6 4.6 3.6 4.4 4.1 4.8 5.0 3.7 4.1 4.1 3.4 3.6 5.8
##  [109] 3.6 4.6 5.2 5.8 3.7 3.6 3.1 4.5 6.1 5.6 4.6 6.3 5.8 4.1 4.1 3.0 4.9 6.0
##  [127] 6.4 5.0 3.1 3.9 4.5 5.4 5.9 2.9 4.3 5.2 4.6 4.8 4.1 3.5 5.2 3.4 3.8 5.2
##  [145] 4.8 3.8 4.4 3.7 6.0 5.4 5.8 4.9 4.5 4.6 5.3 4.4 3.5 6.2 3.3 4.8 4.2 3.5
##  [163] 5.9 3.6 5.2 5.5 5.3 3.9 4.8 6.5 5.7 5.4 5.0 5.1 4.6 5.5 3.7 3.4 4.2 2.1
##  [181] 4.3 6.6 4.4 3.5 4.4 4.3 6.0 6.1 4.0 4.5 5.6 4.5 5.3 2.9 4.8 6.5 4.1 5.5
##  [199] 5.3 4.9 4.1 4.9 4.6 4.2 4.1 5.5 3.8 5.4 5.1 5.4 5.2 3.4 5.9 5.3 4.7 2.1
##  [217] 4.1 3.7 2.7 4.7 4.6 5.0 5.2 2.1 5.2 3.4 2.7 5.1 3.2 5.2 5.1 4.1 6.4 3.3
##  [235] 4.9 5.0 4.2 6.1 5.1 4.0 3.7 4.2 3.1 3.1 3.3 3.3 3.1 3.1 5.6 3.4 4.7 4.5
##  [253] 4.0 4.2 4.4 5.0 5.0 4.5 5.5 5.8 4.3 4.8 4.3 4.9 3.9 4.4 5.0 4.7 4.6 4.8
##  [271] 4.3 4.8 5.1 5.6 2.6 5.3 4.3 5.1 4.4 4.0 5.4 3.9 2.6 4.0 4.0 3.8 5.1 6.0
##  [289] 4.8 4.0 4.7 4.0 4.9 3.5 5.4 3.7 4.2 4.5 3.8 3.5 5.5 5.0 4.0 2.7 4.3 3.7
##  [307] 5.3 5.5 3.9 5.6 3.9 3.5 6.5 5.3 6.5 3.9 4.8 4.5 4.0 4.8 3.8 6.3 3.4 5.9
##  [325] 6.3 2.9 4.8 5.1 4.4 4.3 4.5 4.9 4.7 4.2 2.3 4.9 3.5 4.7 3.9 5.4 2.9 4.2
##  [343] 3.9 6.0 3.6 3.6 5.7 5.1 5.6 4.4 4.2 3.8 4.5 5.7 5.9 4.9 5.3 4.6 5.3 3.1
##  [361] 4.7 3.6 2.9 4.8 2.9 3.7 4.7 5.2 5.8 4.7 4.5 6.4 4.5 3.5 2.7 3.7 4.8 2.7
##  [379] 5.1 4.2 5.0 4.1 5.7 5.8 4.2 5.8 5.2 4.9 4.7 5.0 4.7 4.9 5.2 3.1 4.9 5.0
##  [397] 3.1 5.1 3.8 4.0 3.9 4.3 5.0 4.3 4.4 4.2 3.2 4.5 4.4 3.7 5.3 4.4 5.3 3.9
##  [415] 4.0 5.8 4.9 4.6 3.8 5.0 5.1 4.2 4.9 5.9 3.7 4.0 3.5 4.0 3.9 3.8 5.1 2.3
##  [433] 4.9 2.5 3.4 3.2 4.2 2.6 3.7 1.8 3.9 3.6 4.1 4.8 5.0 3.6 3.9 3.0 4.1 4.6
##  [451] 5.1 3.7 5.4 3.2 3.8 5.0 2.9 3.5 3.8 6.3 7.0 3.2 5.9 4.1 4.0 4.8 4.1 4.9
##  [469] 3.8 4.2 3.9 4.5 6.5 4.6 4.5 4.3 4.3 4.5 3.8 4.3 4.9 4.3 1.3 4.9 3.3 4.5
##  [487] 3.9 4.0 5.1 3.5 5.2 2.7 6.3 4.8 4.5 5.5 4.8 5.3 5.8 3.6 4.8 2.4 2.9 5.6
##  [505] 4.0 4.1 6.4 4.6 3.3 4.4 4.4 4.7 3.6 4.9 4.9 6.5 5.1 3.5 4.4 4.5 3.5 4.2
##  [523] 6.0 3.7 4.5 6.1 4.6 5.3 4.7 4.1 3.5 5.1 5.7 4.9 4.1 4.0 5.2 4.0 4.7 4.2
##  [541] 5.1 4.1 4.8 4.2 4.5 5.1 3.8 4.2 4.4 4.0 3.8 6.6 4.3 4.1 4.6 5.7 4.5 4.4
##  [559] 3.1 3.0 6.6 3.7 5.2 5.8 4.8 4.5 3.2 2.8 4.7 6.1 4.0 6.0 4.3 6.2 3.7 4.8
##  [577] 5.0 4.2 5.1 4.4 4.2 5.0 5.3 4.1 4.6 6.1 4.9 3.5 4.4 4.8 6.9 5.2 4.3 4.3
##  [595] 3.7 6.3 5.6 4.6 3.7 4.9 5.2 4.3 3.4 6.4 6.2 4.9 3.2 2.9 3.6 5.4 5.2 4.4
##  [613] 3.9 3.7 4.0 3.7 4.5 3.7 4.0 6.3 4.3 4.8 4.1 3.7 4.1 5.9 5.3 4.7 3.7 4.1
##  [631] 3.8 4.2 5.0 3.8 5.7 4.7 4.8 5.1 1.9 4.9 5.0 3.8 3.1 5.6 4.3 4.1 3.6 3.6
##  [649] 5.4 4.1 5.1 3.9 3.4 4.0 5.7 4.2 4.6 3.0 3.8 4.6 6.4 4.4 5.3 4.8 4.5 3.9
##  [667] 4.6 5.0 7.0 4.5 5.0 3.2 3.9 5.3 4.8 4.3 4.6 4.2 4.7 3.8 3.5 5.7 5.9 5.2
##  [685] 5.3 4.3 2.5 4.2 6.0 3.5 6.6 5.0 4.2 5.0 4.0 6.1 4.3 5.1 3.6 4.2 4.8 3.6
##  [703] 5.2 3.6 5.3 2.6 5.5 5.3 4.2 4.4 5.2 4.9 6.3 5.1 3.1 5.5 5.1 3.4 4.6 4.7
##  [721] 4.4 3.6 2.7 3.5 5.1 3.8 5.2 4.0 6.0 3.6 4.3 4.0 3.0 5.9 3.2 3.3 4.9 3.5
##  [739] 4.8 5.0 4.4 2.8 4.4 6.3 5.2 4.4 4.5 3.6 4.1 3.9 4.1 3.0 3.5 3.5 5.6 3.5
##  [757] 3.8 5.6 5.5 5.8 4.1 5.9 4.0 4.7 6.6 3.9 3.4 5.0 6.5 4.9 3.6 3.4 4.9 5.0
##  [775] 3.5 5.9 5.0 5.6 4.1 3.7 5.5 4.6 4.2 5.1 5.1 3.8 5.0 4.7 3.8 5.2 3.8 4.6
##  [793] 4.1 4.1 5.5 4.1 5.1 5.2 4.8 4.4 4.3 4.2 2.7 5.1 4.0 4.1 4.9 2.9 6.1 3.8
##  [811] 4.5 4.3 3.0 3.8 3.8 4.4 6.3 4.2 4.8 4.6 5.2 3.9 2.4 3.8 5.3 4.4 3.9 5.9
##  [829] 5.2 4.3 4.8 5.1 3.8 5.0 4.5 4.3 4.6 4.0 5.4 4.2 4.6 5.9 4.7 4.5 3.5 4.4
##  [847] 2.8 3.7 4.0 3.9 5.0 4.5 4.2 5.8 3.8 5.0 6.0 5.8 3.5 4.4 3.9 4.6 3.4 4.2
##  [865] 3.3 5.1 5.1 3.9 5.5 4.0 4.0 3.7 2.8 5.0 5.1 4.4 4.0 3.8 4.1 3.9 3.1 4.6
##  [883] 5.0 4.0 3.7 4.6 5.5 4.4 4.6 4.5 3.6 3.6 4.6 2.9 3.9 4.4 5.2 3.5 4.9 5.2
##  [901] 4.0 3.0 6.4 5.2 4.8 5.2 4.7 2.7 4.7 5.5 4.6 3.2 6.0 4.1 4.1 4.2 3.4 4.8
##  [919] 3.7 3.2 4.3 5.4 4.5 5.1 4.6 4.5 4.5 5.7 3.4 4.2 5.2 3.7 3.2 5.4 4.4 3.9
##  [937] 5.9 5.7 3.7 4.7 4.3 6.0 5.2 3.9 3.3 4.5 5.3 4.0 4.5 5.2 5.1 5.1 2.3 4.0
##  [955] 3.2 3.4 6.1 4.5 5.2 4.7 4.7 4.9 6.6 6.4 4.8 4.5 3.8 6.1 5.1 4.9 3.6 6.4
##  [973] 4.4 3.7 5.7 2.9 3.4 4.7 6.0 4.1 5.0 2.7 5.9 4.6 3.8 5.6 2.9 4.4 6.4 4.7
##  [991] 4.2 4.2 5.5 5.2 4.6 5.4 4.8 4.7 5.4 3.8
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.460 5.500 5.288 5.200 5.100 5.000 4.900 4.800 4.600 4.500
## 
## $jack.boot.se
## [1] 0.9686022
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
## [1] -0.5062332
```

What is skew? Skew describes how assymetric a distribution is. A distribution with a positive skew is a distribution that is "slumped over" to the right, with a right tail that is longer than the left tail. Alternatively, a distribution with negative skew has a longer left tail. Here we are just using it for illustration, as a property of a distribution that you may want to estimate using your data.

Lets use 'fitdistr' to fit a gamma distribution to these data. This function is an extremely handy function that takes in your data, the name of the distribution you are fitting, and some starting values (for the estimation optimizer under the hood), and it will return the parameter values (and their standard errors). We will learn in a couple weeks how R is doing this, but for now we will just use it out of the box. (Because we generated the data, we happen to know that the data are gamma distributed. In general we wouldn't know that, and we will see in a second that our assumption about the shape of the data really does make a difference.)


```r
library(MASS)
fit<-fitdistr(original.data,dgamma,list(shape=1,rate=1))
```

```
## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced
```

```r
# fit<-fitdistr(original.data,"gamma")
# The second version would also work.
fit
```

```
##      shape       rate   
##    6.464232   17.845911 
##  ( 2.819446) ( 8.093985)
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
## [1]  0.2450764 -0.1613810  0.2885382  1.3546794  2.0024097  0.1606588
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
##    [1] -0.9945847462 -0.8752931097 -0.3446082388 -1.2397416762 -0.7644068768
##    [6]  0.0593807407 -0.7709767537 -0.7040122649 -0.4630146937 -0.2274579494
##   [11] -0.9551869239 -0.6700431445 -1.2417851685 -0.6735331874 -0.4758563717
##   [16]  0.6681867393 -0.7340689080 -0.2530391259  0.1785486320 -0.3003406205
##   [21] -0.7830391162 -0.4636027300 -0.6190613840 -0.6426655958 -1.1436465263
##   [26] -0.4811898455 -1.4024416949  0.1011732796 -0.6741204220 -1.0966663475
##   [31] -0.6633298790 -1.2119760972 -0.7840049460 -0.1134185275 -0.2649714057
##   [36] -0.6459955881 -0.2052429456 -0.6067308277 -0.7187661783 -0.4090671457
##   [41] -0.3871385358 -0.4232604347 -0.3308811336  0.4435908062 -1.3187832765
##   [46] -0.5435955569 -1.2597373938 -0.3126065981 -0.4649507600 -0.9686759682
##   [51]  0.1878958202 -1.2431874819 -0.1905193714 -0.6789280966 -0.5484015198
##   [56] -0.7845181457 -0.8483726633  0.1154592440  1.0411691925 -1.0588684609
##   [61] -0.5525238920  0.0376619074 -2.0771896534 -0.9911042218 -0.1808989386
##   [66] -1.1198685228 -0.0096410250 -0.2484012205 -0.0021391422 -0.0411503875
##   [71] -0.4268972494 -0.9303424821 -0.6811714218 -0.9001905592 -0.5849387835
##   [76] -0.3938051169 -0.4878612483 -0.7062203522 -0.9743021464 -1.1253049556
##   [81] -0.8516456239 -0.8815966002 -0.8536286279 -0.2454702775  0.1928589373
##   [86] -0.8288963544 -0.1850611930 -0.2864624248 -0.5290852489 -0.6565447303
##   [91] -0.0222732815  0.1485034569 -0.8256373775 -1.2404965796 -0.6057753700
##   [96] -0.1206137579 -1.0795942102  0.5298411965 -0.1847509896 -0.9564743030
##  [101] -0.3914560586 -0.1893248399  0.0923906550 -0.9858467035 -0.9065465376
##  [106]  0.0262851705 -0.4827868469 -0.5190524193 -0.2712422247 -0.3939852763
##  [111]  0.4890915341 -0.0546536647 -0.3569100393  0.5815581152 -0.1368182740
##  [116] -0.0171058314 -1.1430792522 -0.5345704827  0.1327250591 -1.5999306530
##  [121] -1.3977637273 -0.3557461095 -0.9662697040  0.0937994459 -0.7126709421
##  [126] -0.0690689978 -0.3379913394 -0.5425294379  0.0720147708 -0.0912270350
##  [131] -0.7402651069 -0.4021755006  0.2872157799  0.4586486370 -0.7522356327
##  [136] -0.0549066525 -0.9631695212 -0.6585643668 -0.6424462860 -0.3515521282
##  [141] -0.6243206452  0.1130521889 -0.3550052297 -0.7960080986 -0.6478795310
##  [146]  0.6678156318  0.4620980521 -0.5921580759 -0.2615307465 -0.7769867546
##  [151] -0.6527305957 -0.7520920133 -0.3023057817 -0.8148727939 -0.6661301028
##  [156] -0.9324898248 -0.2156670478 -0.0185426190  0.2097390525 -0.4454970932
##  [161] -0.4341431499  0.2231893713 -0.7408898856 -0.5233176451 -0.6364857870
##  [166] -0.2729991844  0.1034173477 -1.1653651680 -0.2122266662 -0.0147290212
##  [171] -1.2134395175 -0.2421263126 -0.3596177250  0.0503425901 -0.9123662560
##  [176] -0.4630752483 -0.3710862003  0.0302416746 -0.2317101420  0.0009160696
##  [181] -0.5035732360 -0.5890303929  0.2854132970 -0.1897499474 -0.7977929116
##  [186] -0.3760392410  0.1323327823 -0.3991427442  0.3528973131  0.3918998130
##  [191] -0.1178853693  0.1437284796 -0.4503590982 -0.0855103732 -0.0886665908
##  [196] -1.0099580953  0.4512451269 -0.1267329790 -0.4023490567 -0.3601284472
##  [201]  0.1827997338  0.3540435559  0.0559609977 -0.5870109691 -0.0799719203
##  [206] -0.1869422705 -1.0007782833 -0.1242323420 -0.4638666929  0.3995273636
##  [211] -0.1932380796  0.1204679877 -0.6262048207 -0.5039474418 -0.9216409472
##  [216] -0.3670039608  0.4859107094 -0.0244747651 -0.4589307696  0.2974905408
##  [221] -1.1274651215 -0.3693489221 -0.3481406613  0.0596534194 -0.9533052257
##  [226] -0.5952112579 -0.1883268919 -0.7966918661 -0.5859105749 -0.6425275741
##  [231] -0.8152024145 -0.2723837722 -0.6692279589 -0.8070631089 -0.1290973493
##  [236] -0.7559557381 -0.1420454498 -0.8990061674  0.0467170779 -0.9736902367
##  [241] -0.4315634850 -0.2067807981 -0.7437607469  0.0669794513 -0.2412582588
##  [246] -0.4476421362 -0.3681397505 -1.2975214175 -1.1734964342 -0.6126119679
##  [251] -0.2083678037 -0.0218335259 -0.8542219858 -0.1877024862  0.0341809397
##  [256] -0.2412159343 -0.3138028184  0.0449558345 -0.2135733216 -0.2212347268
##  [261] -0.5066185884 -1.4428144664 -1.2495641038  0.1829457465 -0.8581234145
##  [266] -1.0744319349 -0.2306336586 -1.1423055761 -0.8312577225 -0.2694526068
##  [271]  0.2385345695  0.1543444478  0.4134030128 -0.4117074258 -0.9506610183
##  [276] -0.5954333237 -0.4814720105 -0.5809313438  0.2526843835 -0.1122946625
##  [281] -0.2336953234  0.0142415720 -0.2750979565 -0.7959547083 -0.2794263946
##  [286] -1.2673156288 -0.1976653805  0.0063894612 -0.0613774370  0.2613890460
##  [291] -0.3386250211 -0.3732418068 -1.0919204963 -0.4761276792  0.0436296045
##  [296] -0.3312106650 -0.7277716242 -0.6561715554 -0.5843320512 -0.9487181090
##  [301] -0.0732173527 -0.5817237397 -0.3753958121 -0.9862151895 -0.7021344151
##  [306]  0.2091930116 -0.9869152064 -0.1599443301 -0.2468427801 -0.8625443371
##  [311] -0.0799519978 -0.9324898248 -0.0137439727 -1.2256405277 -0.6527991065
##  [316] -1.9866521533  0.2798996628  0.0253867956 -0.6515013962 -1.0255209934
##  [321] -0.3550052297 -1.0487677619 -0.7577267422 -1.0514174863 -1.2245211407
##  [326]  0.4578343984 -0.4129551464 -0.4842039026 -1.4689843532 -1.4839394809
##  [331]  0.9674455335 -0.8256373775 -0.3424200470 -0.8331373840 -1.3029191416
##  [336]  0.0715104518 -0.2337099284 -0.2129017816 -1.1487026725 -2.1521659023
##  [341] -0.4684443400 -0.6259081989 -0.7709767537 -0.9417299038  0.0083024728
##  [346] -0.5168071495  0.0850479946 -0.7245969433 -0.5918709994 -0.5833177784
##  [351] -0.3068264629 -0.7842416591 -0.8210506182  0.2137874950  0.8376779980
##  [356] -1.6961120058 -0.9880364801 -0.1983263050 -0.9008539823 -1.0045576290
##  [361] -0.3353701919 -1.4136294673 -0.4636986244 -0.4755073840 -0.8493147487
##  [366] -0.7464889195  0.3913867500 -0.4436490606 -0.4452957383 -0.3094761074
##  [371] -0.3902039761 -1.1765276005 -0.5157540086 -0.7634453974 -0.7102308062
##  [376] -0.2439909283 -0.9724125690 -0.0266171448 -0.2854163382 -0.9931218056
##  [381] -1.1272266569 -0.3515184973 -1.0573454289 -0.2833915783 -0.6364827843
##  [386]  0.6312185898  0.2768155632 -0.6074003902 -0.1409228776 -0.5539736683
##  [391] -0.5127889468 -0.7391464027 -0.5902840885  0.1255860749 -0.8122559172
##  [396] -0.6797175522 -0.1726977589  0.1557651319 -0.0865479609 -0.1911596944
##  [401] -0.2478689419  0.0925049705 -0.4366272761 -0.5733364006 -1.0168889057
##  [406] -1.3087797108 -1.4789292451 -0.6669515355 -0.8097395624 -0.8248428488
##  [411] -0.9688405358  0.1843880480 -0.6850292722  0.7949034513 -0.5818658437
##  [416] -0.8871311123 -0.1598568802 -0.3136049752 -0.5164470219 -1.9579459529
##  [421] -1.2354867944 -0.4175321738 -0.0104612645 -0.6899108130  0.2917195837
##  [426] -0.6148819885 -0.4243267667 -0.7664248470 -0.2294815262 -0.8232293392
##  [431] -1.0546783818 -0.4831836432 -0.3099154593 -0.0926184831 -0.6721078472
##  [436] -0.7681700126  0.0207665413 -0.4031750465 -0.5435215828 -0.1818104491
##  [441]  0.0840446112 -0.4687994108 -0.0238587666 -0.0346185008  0.3630339389
##  [446] -0.7279466087 -0.5175354279  0.4051159755  0.0433565616 -0.6064379675
##  [451]  0.2626727038 -0.5189445888 -0.2273801101 -1.3029278539 -0.5281504522
##  [456] -0.0733530452 -1.1094221083 -0.4512680634 -1.6208843296 -1.1774137991
##  [461] -0.2985650452 -0.6360721785 -0.9127969561 -0.7916493575 -0.9243690037
##  [466]  1.3697089867 -0.3988184093 -0.9143771244 -0.0206206123 -0.9950845454
##  [471] -0.3290760864 -1.0068041461 -0.5472267654 -0.5867354445 -0.5625927295
##  [476] -0.2485549038 -0.4204914939 -0.7873821112 -1.0434886291  0.1304643401
##  [481]  1.1202427196 -1.2386920146 -1.4456124215 -0.4819406755 -0.6424345105
##  [486] -0.5360345815 -0.4484063868 -0.2855951084 -0.8749025158 -1.5258261553
##  [491]  0.2066045696 -0.1615091476 -0.3429898404 -0.2362490045 -0.3303856837
##  [496]  0.0119178457 -0.8372177074 -1.0317484864  0.8445076659 -0.5242079345
##  [501] -0.1615091476 -0.9586124987 -0.4300414487 -0.3982462042 -0.1193067028
##  [506]  0.3696332191  0.3136408261 -0.9934061391 -0.1638744269 -0.7559112864
##  [511] -0.3528670640 -0.3337271421 -0.6495561233 -0.0368232530 -0.3790816458
##  [516] -1.2511489165  0.1323525358 -0.4049177307  0.8118233057 -0.7493917081
##  [521] -0.8672216942 -0.3067538002 -0.7826899400 -0.4922156048 -0.3369368397
##  [526] -0.1752291087 -0.5903541055 -0.6578691258 -1.0253858533 -0.7312829802
##  [531] -0.4517013254 -0.5728850843 -1.4083977557 -0.9419334002  0.1003609037
##  [536]  0.3506665931 -0.5271069061 -0.5612671117 -0.4877319862 -0.3422786719
##  [541] -0.6228861875 -0.5872832480 -0.7973172910 -1.2604720952 -0.6959067837
##  [546] -0.3430106612 -0.5225371842  0.1711445311 -0.9806541809 -0.1958223672
##  [551]  0.1761109969 -0.8606482171 -0.6020886748 -0.3729029817 -0.3811340260
##  [556] -0.2111630254 -0.2376610513  0.0029293280 -0.7211192005 -0.9864568877
##  [561] -0.2702637536 -0.1410000161 -0.6267986698 -0.2564053502  0.0222225992
##  [566] -0.7002640222  0.4289136766 -0.2940161346 -0.9831231533 -1.3516680738
##  [571] -0.2038038331  0.1277199841 -0.3004355239 -0.6048539151 -0.2919776281
##  [576] -0.9391411592  0.8413315314 -0.1423629777 -0.5521810227 -0.6301501182
##  [581] -0.6004894992 -1.4817858930  0.2203811078  0.1017696839 -0.4073972208
##  [586] -0.7617715330 -0.2487335769  0.0964560418 -0.0211709465 -0.3488434768
##  [591]  0.0553916848 -1.4636151406 -0.4367203370 -0.7671362014 -0.9212275788
##  [596] -0.6432022941 -1.2450592779 -0.4185923383  0.0502195274 -0.1039959118
##  [601] -0.1889966466  0.4138808612 -0.5242079345  0.0885499802  0.1536762991
##  [606] -0.5768669813 -1.4357274097 -0.6147278649 -0.5774653664 -1.3053564501
##  [611] -0.5759514444  0.0902288418 -0.9803228760 -0.0643241847 -0.3432593035
##  [616] -0.3367097228  0.1880602440 -0.6126027069 -0.7277772293 -0.6478795310
##  [621] -0.8874695429 -0.9290909741 -0.1352327380  0.1675073370 -0.0405728805
##  [626] -0.1718412293 -0.7537420098 -0.6480125679 -1.1819179034 -0.2669080555
##  [631] -0.8485879518 -0.8900204246 -0.0824399858  0.0933246247 -0.2530918220
##  [636] -0.4182413576  0.1679310259 -0.5533488839 -0.4116531699 -0.7222174157
##  [641]  0.0452757411 -0.6854132278  0.1065340353  0.9342063487 -0.2231221360
##  [646]  0.5769749521 -0.6417237499 -0.3239461712 -0.2840450638 -0.8874384928
##  [651] -2.0178129452  0.1592346874 -0.3007905454 -0.6193420554 -0.7556954171
##  [656] -0.5190524193 -1.5538920620 -0.8595153572 -1.6554790197  0.2933004566
##  [661] -1.0099195012 -0.5990839825 -0.1700568861 -1.4026491177 -0.6819663217
##  [666] -0.6369453090 -1.0215968803 -0.0081320710 -0.0467874033  0.0720147708
##  [671] -0.0856286686 -0.6121629717 -0.4809160304 -0.7002640222 -1.3105074716
##  [676] -0.5358770453  0.5179184985 -0.5110982171 -0.6153372463 -0.2997082470
##  [681] -0.1895335388  0.2025159656 -0.1649447844  0.1269748321 -0.4445899920
##  [686] -1.0118962125  1.6938895401 -0.5289389445 -0.7709492642 -1.1835196234
##  [691] -0.8546429374  0.3028290392  0.5156119710 -0.6618761523  0.3744226163
##  [696] -0.7399143424 -0.0286798097 -0.6745360520 -0.7355299870  0.1939305648
##  [701] -0.8225292594 -1.0178615236 -0.1103102659 -1.2545858122  0.0411295813
##  [706] -0.5192657618 -0.5886824016 -0.5262973319 -1.4853359115 -0.7086989247
##  [711] -0.5977883737 -0.4208469827  0.0882746379  0.1053955522 -0.8867380654
##  [716] -0.5539472023 -1.2796364257 -1.5258261553 -0.5498302284 -0.1456463426
##  [721] -0.2503319336  0.4757121716 -0.6746571953 -0.6922962175  0.0933246247
##  [726] -0.1813724772 -0.1058788505 -0.5430077407 -0.4348327155 -0.9762247187
##  [731] -0.0997179958  0.0598036506 -1.3936013922  0.1257985776  0.0169270574
##  [736] -0.9219367097 -0.8526801295  0.2749102013 -1.2862598485 -0.7569104677
##  [741]  0.0544310585 -0.9584606681 -0.6008244572 -0.4746332907 -0.1947537863
##  [746]  1.2387728822  0.5968976844 -0.1734634922 -0.0085376788 -0.1291938039
##  [751] -0.4315045098 -0.3808162098 -0.1863972566  0.1510934328 -0.5746559893
##  [756]  0.2938811407 -0.8237060892 -0.2862345289 -0.6966698835 -1.3246047219
##  [761] -0.1977814868  0.4152499267 -1.3094524842 -0.8116291175 -0.5000271155
##  [766] -0.6760787043 -0.1192338230  0.0580903627 -0.7570537238  0.5761045320
##  [771] -0.7026495412  0.1351325527 -0.4726465627 -0.5845914691 -0.1555046569
##  [776] -0.2774515956 -1.2332746599 -0.0337121993 -0.7467395080 -0.1317221817
##  [781]  0.4484632657 -0.1257378775 -0.3577086146  0.0542917956 -0.6777997963
##  [786] -0.6073767923 -0.3622262849  0.6233183942  0.0151125993 -0.2050618696
##  [791]  0.3759833126  0.0856052826 -0.6556998673 -0.7553134717 -0.2268479025
##  [796] -0.3071828405 -0.4810920248 -0.3942671364 -0.8517424332 -0.3936162206
##  [801] -0.7168961582 -0.3754586102 -0.9189637338 -0.9645181325 -1.5717445819
##  [806] -0.5795460079 -0.9888222967 -0.1980161972 -0.5753051381 -0.0372244449
##  [811] -0.9869152064 -1.1775322228 -0.2649153296 -0.5147504414 -0.1228679524
##  [816] -0.8973864346 -1.0523857376 -0.8976402694 -0.6168523069 -0.6363076369
##  [821] -0.1622416605  0.4729990045 -0.7840829784 -0.5160195356 -0.4422098154
##  [826] -0.9862151895  0.4309892018 -0.5056074608  0.0331225028  0.2361211587
##  [831] -0.2554730326 -1.2880521590 -0.2052631664 -0.0649313646 -0.1875420090
##  [836]  0.1417056484 -0.2463184330  0.7287191879 -0.5729036211 -1.0903178835
##  [841] -0.4051777197 -1.0239838590 -1.0056260288 -0.3010508870 -0.6760941682
##  [846] -0.4211466262 -0.3658757571 -0.0824399858 -0.8444199450 -1.2236341421
##  [851] -0.4869558601 -1.2417377928  0.0860289861 -0.0224862795 -0.6472954103
##  [856]  0.1430856241 -0.3745372891 -0.1995174781 -0.7190679969 -0.0403280367
##  [861] -0.1721970754 -0.4222451552 -0.6923441819  0.0487412495 -1.3082256472
##  [866] -0.5741353760 -0.5341550313 -0.7129592822 -0.3878886856 -0.2298057707
##  [871] -1.1909681775  0.2320179784 -1.4293423404 -0.5916503326 -0.4856418792
##  [876] -0.3067497851 -1.4612083485 -0.3872396555  0.0321123215 -0.3548971443
##  [881]  0.2568608922 -0.9753294001 -0.4542146287 -0.3031944644 -1.5658627004
##  [886] -0.0571165613 -0.7549022002 -0.6666815908 -0.9562592161 -1.6565854004
##  [891] -0.2739936297 -2.0089346875 -0.1152137746 -0.7702916124 -0.4443078862
##  [896] -0.3028526160 -0.6836653981 -0.8090922089 -1.2825328843 -0.6680024566
##  [901] -0.9108082524 -0.2230270909 -0.9025098372 -0.0941738802 -0.9984097348
##  [906] -0.7994356469 -0.4497762981 -0.5291211568  0.0278512155 -0.7098009561
##  [911]  0.0173449152 -0.1527157019 -0.2991061769 -1.1847354858 -0.4163933609
##  [916] -1.0169340402 -1.0259153553 -0.6748164820 -0.8878035619  0.2558593847
##  [921] -1.2731956482 -0.6107635516 -0.5793912945 -0.1556280048 -0.0408258320
##  [926] -0.3491495817 -0.4876174637 -0.3024211614  0.2609328377 -0.9645181325
##  [931] -0.1825668534 -0.0581960773 -0.7604359825  0.1436795397 -0.7455469711
##  [936] -0.1191249807 -0.9482685441 -1.3668002058 -0.3248717228 -0.2246423075
##  [941] -1.1799613413 -0.8799232174 -0.6066232525  0.3441562857 -0.0468017906
##  [946] -0.3505679816 -0.8310575486 -0.6663600896 -0.9758975786 -0.5903148689
##  [951]  0.1866928500 -0.0207934453  0.0056178445 -0.0272673636 -1.2834941196
##  [956] -0.6460801886 -0.4336016019 -0.8037901907  0.3339293162  0.2215834506
##  [961] -0.3208039572  0.3581874483 -0.3379459393 -0.9177567759 -1.2253631443
##  [966] -0.4231431616 -1.2344035125 -0.1530149626 -0.1610995583  0.3554739512
##  [971] -0.2802946922 -0.7528325711 -0.3597926634 -0.9262047486  0.3810982744
##  [976] -0.9011290646 -0.4307486065 -0.2622120053  0.2210379989 -0.5621220433
##  [981] -0.7723266700  0.5110024464 -0.8542375194  0.6857389816 -0.8019417541
##  [986] -0.7176764911 -0.6050672726  0.1354756301 -0.3071299178 -0.7208379104
##  [991] -0.2155586372 -1.3229362970 -0.5205235567 -1.2660296381  0.2840789310
##  [996] -0.8655937850 -0.1029413038  0.2327847107 -1.0261562341 -1.1354070909
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
##   0.36221429   0.12527983 
##  (0.03961696) (0.02800757)
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
## [1]  0.7928999  0.4216782  0.1539632 -0.4440067  0.7133714  0.1246946
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
## [1] 0.0356
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9054571
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
## t1*      4.5 -0.05235235   0.9213521
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 2 3 4 5 8 9 
## 1 2 1 1 2 1 2
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
## [1] 0.0165
```

```r
se.boot
```

```
## [1] 0.9223603
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

