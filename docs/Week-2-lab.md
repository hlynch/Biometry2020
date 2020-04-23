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
## 0 1 2 4 6 8 9 
## 1 2 1 1 1 2 2
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
## [1] 0.0167
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
## [1] 2.766486
```

```r
UL.boot
```

```
## [1] 6.266914
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##   2.5%  97.5% 
## 2.7975 6.2025
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
##    [1] 4.1 4.0 5.1 4.4 3.8 4.3 5.2 3.7 4.9 6.1 4.1 3.6 6.1 4.0 3.8 4.6 3.4 4.6
##   [19] 3.7 5.5 5.4 4.1 4.9 4.3 3.7 4.7 5.0 4.9 4.4 5.9 4.3 3.9 5.5 4.6 4.5 6.1
##   [37] 4.9 5.4 4.8 4.8 3.3 4.7 4.4 2.8 3.8 4.3 3.9 4.1 4.5 4.4 3.6 4.8 5.1 2.0
##   [55] 5.7 6.1 4.8 5.9 4.9 4.3 5.6 3.3 4.2 4.9 5.2 4.4 5.2 5.3 4.6 5.0 6.1 5.2
##   [73] 4.0 4.9 4.1 4.1 4.0 2.8 5.4 3.1 3.8 3.7 5.1 2.2 5.3 3.3 4.3 5.4 3.6 5.1
##   [91] 3.3 3.8 3.3 4.6 2.0 4.5 4.7 5.2 6.5 4.3 6.3 3.6 3.7 4.2 4.9 6.2 4.5 4.8
##  [109] 4.4 4.7 4.4 3.9 5.5 4.9 4.0 5.6 3.2 4.0 3.9 2.8 5.1 3.6 5.9 4.8 2.6 5.2
##  [127] 3.6 4.3 5.8 3.9 3.9 3.3 4.9 5.2 4.7 6.1 5.1 5.6 3.4 4.8 5.5 4.4 4.9 5.3
##  [145] 3.8 3.5 5.0 4.8 3.9 4.4 4.5 5.3 5.2 4.6 5.0 5.6 2.4 3.8 5.5 3.0 2.2 4.9
##  [163] 4.0 3.0 5.0 3.7 5.8 3.2 4.6 4.6 3.7 3.6 3.9 2.8 5.3 3.9 4.1 4.9 4.6 3.2
##  [181] 4.2 4.3 5.4 3.7 4.1 5.3 3.8 4.0 4.3 4.2 3.7 3.0 3.1 5.0 3.1 5.2 3.5 3.8
##  [199] 3.1 4.7 4.0 3.9 3.8 4.8 4.0 5.0 4.3 3.0 4.9 3.9 4.3 4.2 4.9 4.0 4.9 3.8
##  [217] 4.1 3.4 3.6 3.4 3.9 3.1 4.3 5.4 5.1 4.2 3.5 4.1 2.4 4.3 5.3 4.2 3.4 6.7
##  [235] 5.7 3.0 4.5 5.8 3.6 5.3 5.7 4.0 4.4 4.9 3.0 4.1 3.0 5.0 3.1 4.2 3.8 4.5
##  [253] 3.9 5.9 4.9 3.8 4.5 4.7 5.0 4.1 5.9 4.2 4.6 3.6 4.4 5.4 5.3 3.7 6.3 5.1
##  [271] 4.8 4.2 5.0 3.9 5.4 3.3 5.9 4.4 2.6 4.9 5.3 4.4 4.3 4.3 5.7 5.0 2.6 5.1
##  [289] 4.7 5.8 4.6 5.2 5.7 4.3 4.6 4.2 4.5 3.2 4.6 5.2 4.8 5.1 5.4 2.1 4.3 4.0
##  [307] 4.5 4.5 3.5 4.4 3.0 3.7 5.6 4.8 4.7 4.2 3.2 4.7 5.3 5.4 5.5 3.0 4.5 4.2
##  [325] 5.8 4.3 5.4 4.6 4.3 4.8 5.7 5.2 5.2 5.2 4.6 4.7 5.1 4.5 4.3 4.3 4.0 5.5
##  [343] 4.4 4.7 4.0 4.8 4.9 3.9 5.1 4.1 4.5 4.3 4.5 5.0 4.9 3.2 4.9 3.8 3.9 3.9
##  [361] 2.9 6.9 4.7 6.0 4.6 4.8 5.1 4.5 2.5 4.1 4.4 4.5 3.8 3.5 4.3 3.6 6.1 4.2
##  [379] 5.6 3.9 4.4 5.2 3.8 3.3 5.2 2.8 4.5 4.9 4.2 4.1 4.9 3.8 4.2 3.4 5.2 3.1
##  [397] 5.4 3.6 4.5 5.3 3.9 3.9 4.1 3.4 4.7 3.5 2.9 4.7 4.5 3.4 5.8 4.1 3.5 4.7
##  [415] 4.5 4.0 3.7 6.0 4.7 3.0 3.4 4.9 4.3 3.2 3.7 4.0 4.7 5.9 4.8 4.5 6.1 3.1
##  [433] 3.7 3.4 4.4 5.1 3.7 5.2 4.8 4.8 2.4 3.8 4.7 5.1 5.3 3.0 3.8 4.0 4.1 5.5
##  [451] 4.1 4.5 3.6 4.1 3.7 3.7 4.8 5.0 4.4 5.1 4.3 4.3 3.9 5.1 5.7 4.9 2.8 5.3
##  [469] 5.2 5.8 5.4 3.4 5.3 3.6 4.4 6.6 6.0 5.9 5.8 5.2 4.3 2.9 5.4 5.7 4.7 2.6
##  [487] 3.2 3.2 3.9 3.2 5.7 4.6 4.4 5.7 5.3 5.7 4.5 4.4 3.4 4.0 4.9 6.9 6.1 4.2
##  [505] 2.8 3.5 4.4 4.6 5.7 5.0 4.2 4.6 3.8 2.8 5.3 3.5 4.1 4.3 5.6 4.9 4.2 3.7
##  [523] 4.8 6.1 3.5 4.1 3.3 5.3 5.2 4.6 5.1 3.0 5.2 4.2 4.4 3.2 5.4 3.5 6.0 4.4
##  [541] 4.7 5.8 5.0 3.1 5.0 2.6 5.1 4.3 5.1 3.1 5.0 4.5 6.2 4.5 3.9 5.6 4.3 4.2
##  [559] 4.1 4.4 3.8 5.1 4.7 3.9 2.1 3.8 4.3 4.1 4.5 4.9 4.8 4.8 4.6 4.6 5.5 7.0
##  [577] 5.1 5.0 5.2 3.9 4.2 3.3 5.4 3.8 4.8 6.0 4.5 5.7 4.1 4.6 3.4 4.6 6.4 4.5
##  [595] 4.9 3.6 4.0 3.4 4.8 5.6 4.8 3.9 5.1 2.7 5.6 4.3 6.0 6.7 3.4 5.8 4.2 4.9
##  [613] 6.6 3.8 5.4 5.5 3.1 4.3 4.6 4.6 4.1 4.2 5.0 3.6 6.0 3.3 4.8 4.3 5.0 4.1
##  [631] 4.7 4.2 3.7 3.8 5.2 4.9 4.7 5.6 4.3 3.9 3.9 4.6 3.6 6.3 4.1 4.7 6.1 5.5
##  [649] 6.2 4.2 4.6 5.0 3.6 3.9 3.4 4.2 5.0 4.3 4.6 4.5 5.4 3.3 5.3 5.1 5.5 3.9
##  [667] 3.7 4.8 4.9 4.3 3.5 6.0 2.6 4.2 5.1 3.4 4.7 3.7 4.5 6.5 4.2 3.9 4.0 3.6
##  [685] 3.4 3.1 2.6 5.0 4.6 3.0 2.1 4.2 3.7 4.7 5.5 4.5 6.2 4.2 2.9 3.9 4.9 2.8
##  [703] 3.2 6.0 5.8 5.5 4.0 5.3 4.0 4.3 6.1 4.5 4.9 6.4 4.5 5.1 5.0 4.8 4.4 4.6
##  [721] 4.5 3.2 5.8 4.7 3.2 4.8 3.6 4.7 4.2 4.1 5.0 5.0 4.2 3.9 4.7 4.9 4.7 4.0
##  [739] 5.8 4.7 5.4 4.2 3.8 4.2 5.7 4.1 4.2 4.8 4.5 3.9 4.1 5.2 4.7 5.6 6.0 4.0
##  [757] 4.1 3.8 6.4 6.3 6.9 3.6 2.9 3.2 5.5 4.3 5.5 4.8 5.8 3.7 4.5 5.0 4.3 4.4
##  [775] 4.9 5.3 3.3 5.0 2.7 3.1 3.6 5.1 4.4 2.6 4.4 4.4 3.5 5.4 5.0 5.7 4.4 3.1
##  [793] 2.2 3.4 5.1 6.0 4.8 4.9 3.8 4.0 4.4 4.4 4.8 3.9 5.1 4.4 4.9 4.4 2.9 6.0
##  [811] 2.9 5.7 3.5 3.4 4.4 4.1 3.9 5.1 4.3 4.8 4.2 5.2 3.5 4.1 4.6 3.5 5.3 3.7
##  [829] 5.2 2.1 4.7 4.3 5.1 4.2 3.5 4.4 4.6 5.3 5.2 5.1 4.4 5.8 4.8 3.9 4.8 4.7
##  [847] 4.4 5.7 5.1 4.6 4.2 3.6 3.3 6.5 4.0 3.0 4.6 4.0 5.1 5.4 4.6 3.2 5.1 5.0
##  [865] 3.1 5.1 4.5 4.7 4.1 3.4 4.2 4.8 3.3 4.0 3.8 6.0 3.5 5.2 4.3 5.7 5.0 7.0
##  [883] 4.5 4.1 4.5 3.6 4.1 4.3 3.1 5.1 5.3 4.9 4.2 4.6 3.6 5.4 3.3 4.6 4.4 3.7
##  [901] 5.5 5.2 4.4 3.0 4.5 4.4 4.7 3.0 4.0 4.5 4.9 6.2 6.6 4.3 3.4 4.8 5.4 3.9
##  [919] 5.4 5.0 3.9 3.5 4.5 3.4 5.2 3.9 3.8 4.3 4.7 5.4 4.2 4.8 4.9 6.4 4.9 2.7
##  [937] 4.6 5.6 4.4 4.5 4.4 4.0 6.5 6.2 4.1 3.0 4.0 4.1 5.2 2.2 4.9 3.7 3.7 4.4
##  [955] 5.6 5.1 4.9 5.3 4.7 4.6 4.8 4.4 5.1 5.5 5.2 3.9 4.5 4.6 4.6 5.9 3.5 3.6
##  [973] 4.2 4.4 4.6 4.0 4.1 4.0 4.6 5.5 4.5 3.9 3.5 5.0 2.0 4.4 5.1 4.4 3.8 4.4
##  [991] 5.1 5.6 5.9 4.1 3.4 3.6 5.5 6.1 4.8 5.6
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
##    [1] 5.5 4.8 6.8 3.5 6.2 3.3 5.6 3.6 5.3 3.4 2.8 5.6 5.4 4.5 3.6 5.3 6.4 5.2
##   [19] 3.8 5.4 3.3 4.9 3.4 5.6 4.1 4.3 3.6 3.5 4.5 5.6 4.7 4.4 5.3 5.9 3.4 4.4
##   [37] 4.6 4.9 4.5 5.1 4.1 4.4 4.5 4.2 4.0 4.0 4.9 3.2 4.9 6.0 4.0 4.4 3.4 4.2
##   [55] 4.0 4.2 4.7 5.0 5.1 4.0 4.4 5.0 4.8 3.6 4.1 4.3 4.5 4.2 5.2 5.6 4.7 5.8
##   [73] 6.5 4.7 2.7 6.1 4.5 4.5 5.0 4.4 5.0 4.2 4.7 4.9 4.6 5.5 4.3 6.1 4.6 3.6
##   [91] 4.7 3.3 4.3 4.4 5.5 3.6 4.1 3.9 3.5 3.4 6.5 3.8 5.3 3.2 3.8 4.6 5.3 4.2
##  [109] 5.1 4.4 2.7 5.3 4.1 5.0 4.1 3.9 5.2 3.4 5.7 3.4 2.6 4.9 5.1 4.7 6.0 5.1
##  [127] 5.3 4.3 4.2 5.1 4.2 3.8 4.0 4.0 4.5 4.2 3.2 5.2 4.9 6.7 3.8 3.4 2.2 5.0
##  [145] 6.5 5.9 5.2 5.6 5.1 5.9 4.7 4.3 4.5 4.9 3.9 4.2 4.7 5.1 3.8 5.5 3.5 4.3
##  [163] 5.7 4.9 4.1 2.8 2.4 6.0 4.8 5.2 6.2 6.1 2.9 3.6 3.9 4.4 5.2 4.0 4.0 4.9
##  [181] 4.7 5.0 5.4 5.0 5.6 5.0 4.7 5.3 4.1 4.6 5.3 4.3 2.5 6.2 4.3 4.8 4.1 5.3
##  [199] 3.3 5.7 4.3 5.8 3.4 4.6 3.6 3.0 4.1 4.3 3.2 5.0 5.0 5.2 4.8 6.3 3.8 4.7
##  [217] 5.4 5.1 4.4 3.2 5.3 4.2 6.9 6.7 6.5 5.2 4.1 4.8 4.5 4.3 4.3 1.9 5.2 4.7
##  [235] 5.3 4.6 3.6 3.4 5.3 5.4 3.4 4.1 4.5 4.1 5.6 5.0 6.2 5.6 4.1 4.2 3.3 5.0
##  [253] 4.5 4.8 2.8 4.4 5.9 5.6 5.0 3.7 5.4 4.6 6.3 5.7 4.8 4.0 5.0 4.9 5.5 5.4
##  [271] 4.2 3.8 5.5 2.7 3.3 6.5 3.9 6.2 4.1 3.4 5.0 4.5 4.9 6.2 4.2 3.9 6.2 5.5
##  [289] 6.0 4.3 3.7 4.7 4.0 3.4 6.6 5.1 5.0 6.0 5.8 4.5 5.3 4.2 3.5 5.1 5.1 3.6
##  [307] 5.2 5.5 3.4 3.5 5.4 3.4 5.0 4.6 4.1 4.3 4.6 6.0 3.7 5.9 4.9 4.6 3.8 5.0
##  [325] 5.6 4.7 3.5 4.7 5.1 4.5 4.6 3.3 3.9 4.2 2.6 5.4 5.3 5.8 4.2 4.7 3.6 5.0
##  [343] 5.4 3.6 3.8 3.8 4.7 3.6 3.5 5.8 4.7 4.9 3.0 5.2 4.3 4.7 3.7 4.5 4.0 3.6
##  [361] 5.5 4.9 5.7 4.9 5.4 5.8 3.1 6.0 5.7 5.1 5.0 5.3 3.9 5.6 4.5 6.2 4.1 6.6
##  [379] 4.6 5.4 5.4 4.9 4.2 5.3 3.9 3.5 5.0 5.0 5.4 4.6 5.0 3.8 4.8 3.7 4.3 4.3
##  [397] 6.6 3.7 4.5 4.9 3.3 3.2 5.5 4.8 4.9 3.8 4.0 4.7 4.4 5.5 3.1 4.5 5.2 3.9
##  [415] 3.8 3.6 4.4 3.9 4.6 2.8 4.6 4.0 3.6 4.3 4.1 4.4 4.6 3.8 3.7 5.2 5.6 3.8
##  [433] 3.9 4.3 1.6 3.5 3.7 4.8 6.0 4.2 4.9 5.2 4.8 4.5 4.9 6.0 5.4 4.3 2.8 4.6
##  [451] 3.1 3.8 4.8 5.1 3.9 4.2 4.5 4.3 3.6 5.6 6.0 4.7 5.2 4.0 4.6 4.3 4.8 6.8
##  [469] 5.7 4.8 4.9 6.3 5.1 3.8 3.9 3.8 4.2 2.8 4.5 3.3 3.3 4.2 4.9 5.5 4.6 3.8
##  [487] 4.4 5.5 4.1 5.2 4.6 4.3 5.5 2.5 5.8 7.3 5.3 3.4 5.3 4.1 3.7 4.3 4.0 4.3
##  [505] 4.0 3.5 4.6 4.0 5.1 4.8 5.2 4.2 5.9 3.8 3.2 5.2 5.2 4.4 4.6 4.8 5.3 5.7
##  [523] 5.2 5.5 4.3 4.4 4.6 4.5 3.6 3.1 4.1 5.4 5.1 5.0 4.0 4.0 5.5 4.0 4.5 3.9
##  [541] 3.5 4.5 4.8 5.7 3.9 2.8 3.4 3.0 4.9 5.1 3.5 4.9 4.4 5.1 3.8 5.1 3.6 4.5
##  [559] 3.6 4.5 3.7 5.5 3.5 4.7 6.3 5.4 4.2 3.2 4.2 6.2 4.0 3.4 5.4 4.1 3.3 4.0
##  [577] 4.5 4.7 2.3 4.9 5.1 5.2 3.9 4.1 4.4 4.3 4.9 5.5 4.5 4.2 4.2 3.6 4.1 4.8
##  [595] 5.1 3.6 3.0 4.4 3.1 3.6 3.5 6.1 5.5 3.9 4.5 4.9 2.4 5.0 5.2 4.2 3.2 4.2
##  [613] 6.2 5.7 4.0 5.3 6.1 4.6 4.9 3.2 6.1 5.7 5.3 6.2 3.5 5.4 3.8 4.7 4.3 5.3
##  [631] 3.5 5.7 5.3 2.0 6.5 4.3 6.4 3.8 5.3 5.4 3.9 5.2 4.3 4.4 5.8 3.9 4.3 4.1
##  [649] 4.2 4.4 4.9 5.2 4.1 4.3 4.5 2.9 4.2 4.2 5.8 4.6 4.8 3.4 4.5 3.5 4.6 5.3
##  [667] 3.0 5.2 3.4 3.9 5.1 5.4 4.7 4.2 3.4 4.6 5.0 4.7 3.0 4.1 3.4 5.3 4.4 2.7
##  [685] 5.6 5.1 5.0 4.7 3.4 4.0 5.1 4.0 4.9 4.3 5.6 4.7 5.3 5.6 4.7 4.7 5.5 5.8
##  [703] 4.7 2.1 4.1 5.6 4.7 3.9 4.9 3.9 4.0 4.4 3.3 4.6 2.6 4.8 5.2 4.6 5.3 5.0
##  [721] 5.2 4.4 5.9 4.7 5.7 4.2 3.4 4.9 3.0 4.4 6.1 5.2 3.6 4.2 4.0 4.1 5.3 6.0
##  [739] 4.1 5.0 4.7 6.1 3.1 4.9 3.7 4.9 4.9 5.6 6.4 3.7 4.1 4.5 3.6 4.2 5.2 5.6
##  [757] 4.4 3.6 4.0 3.6 5.3 3.8 5.3 5.8 3.7 3.9 4.2 6.4 3.5 4.9 4.8 3.4 5.3 5.1
##  [775] 3.8 3.8 3.4 3.7 4.5 4.3 2.5 6.4 3.8 5.4 4.6 4.5 4.8 3.7 4.3 5.2 5.7 5.4
##  [793] 3.5 2.5 4.1 5.1 4.2 4.1 5.2 4.7 3.5 5.3 4.8 5.1 5.0 6.0 3.3 6.4 5.1 5.7
##  [811] 5.3 5.5 3.4 4.4 6.4 3.2 5.0 5.0 4.4 4.8 4.3 4.7 3.5 5.3 3.4 6.5 4.2 5.3
##  [829] 5.5 4.2 4.6 4.2 5.4 5.6 4.8 4.2 5.8 3.9 6.0 3.1 4.8 3.5 3.6 6.0 4.0 3.0
##  [847] 5.4 3.7 4.8 5.1 5.2 4.2 4.0 4.4 4.5 5.2 4.3 4.6 6.6 4.1 4.8 4.7 4.1 4.9
##  [865] 4.7 4.7 4.9 5.8 4.3 2.8 3.3 5.0 4.4 3.5 3.5 5.1 5.8 4.1 5.4 3.0 6.1 5.0
##  [883] 3.9 3.8 3.6 4.1 3.7 4.0 4.9 3.7 4.8 6.0 3.5 5.8 4.6 4.3 3.5 4.2 4.4 4.1
##  [901] 3.0 5.5 4.7 3.1 4.0 3.4 4.3 4.3 5.7 3.2 4.2 2.6 5.9 5.5 4.6 3.9 5.0 3.8
##  [919] 5.1 3.1 4.8 4.9 4.8 4.5 4.4 2.5 4.1 5.2 5.8 4.9 4.6 5.0 4.2 6.0 5.5 5.1
##  [937] 5.3 5.3 5.1 4.8 3.3 5.4 4.0 4.9 4.1 4.4 4.0 6.3 3.6 4.5 5.7 4.1 5.5 3.0
##  [955] 4.0 4.7 5.5 3.6 4.4 4.6 4.0 4.9 5.0 5.1 4.7 3.5 3.2 5.2 4.5 4.8 3.6 4.9
##  [973] 5.8 3.9 4.4 4.7 4.8 5.3 5.1 5.5 5.2 3.1 4.4 3.2 4.2 4.6 3.2 5.3 5.5 5.0
##  [991] 4.2 4.1 4.8 4.6 4.5 5.1 5.1 3.0 5.0 5.2
## 
## $func.thetastar
## [1] 0.0575
## 
## $jack.boot.val
##  [1]  0.521138211  0.404956268  0.367500000  0.216111111  0.082561308
##  [6]  0.006303725 -0.052173913 -0.176811594 -0.332653061 -0.421282799
## 
## $jack.boot.se
## [1] 0.9041954
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
##    [1] 4.0 4.4 5.1 4.9 5.2 6.2 4.3 4.6 4.6 3.9 5.7 3.9 5.1 3.9 3.0 3.3 3.4 5.9
##   [19] 6.0 3.9 3.5 5.1 6.6 2.8 5.7 5.9 4.6 4.3 4.4 3.4 3.9 4.2 5.2 4.6 3.6 4.4
##   [37] 3.5 4.8 3.8 4.6 4.4 3.6 4.2 4.7 4.0 4.4 4.7 6.2 4.1 2.1 4.5 4.1 4.9 5.1
##   [55] 4.4 5.9 5.7 6.2 3.5 2.9 5.6 5.7 5.3 4.6 3.9 6.6 5.6 5.1 6.4 4.4 4.6 5.4
##   [73] 4.4 4.9 4.4 3.0 4.6 3.3 4.8 4.4 5.2 4.3 6.0 4.7 6.4 4.1 6.9 3.3 3.6 3.5
##   [91] 2.4 4.0 5.2 5.1 6.6 4.2 4.6 1.5 3.7 5.8 6.0 6.3 4.4 3.7 4.6 4.6 3.0 4.6
##  [109] 5.0 3.5 3.6 4.0 4.4 4.2 4.2 3.3 3.7 4.3 4.2 3.6 4.7 4.7 5.3 4.3 4.6 5.7
##  [127] 4.2 2.9 3.8 3.9 4.8 3.2 4.0 4.3 5.2 5.5 5.1 4.3 3.3 5.0 4.4 5.3 4.7 4.3
##  [145] 3.3 3.6 3.9 5.1 3.7 6.0 3.9 4.8 6.6 4.9 4.5 4.8 4.6 4.2 5.6 2.1 2.5 6.0
##  [163] 4.1 6.0 5.3 3.6 3.7 5.3 3.6 4.6 5.3 4.5 4.6 5.6 5.0 5.1 4.1 4.7 4.2 2.5
##  [181] 4.9 5.4 5.1 4.9 4.3 7.1 3.4 4.8 4.5 5.0 3.9 3.4 4.8 3.8 6.6 4.1 3.2 3.8
##  [199] 6.4 2.6 5.1 4.3 5.3 4.9 3.9 3.9 2.9 3.2 4.3 5.5 4.1 5.4 3.7 4.0 6.4 5.9
##  [217] 4.2 4.1 4.9 5.6 3.1 4.1 4.3 4.7 4.8 5.8 4.2 4.0 5.2 3.5 5.4 6.3 4.5 5.2
##  [235] 4.5 3.9 4.4 6.2 3.6 3.6 3.9 4.8 3.7 4.7 3.3 3.5 2.3 4.1 3.7 3.5 2.9 4.5
##  [253] 3.8 4.5 3.3 4.8 4.4 3.8 5.0 4.9 5.3 4.2 3.3 3.9 5.9 4.7 4.6 3.4 6.0 3.6
##  [271] 4.7 4.7 4.3 4.3 2.8 4.5 3.6 4.6 3.0 3.6 6.3 4.0 6.8 3.8 4.8 4.6 5.2 4.6
##  [289] 5.0 4.1 4.1 3.5 4.2 4.0 5.5 4.5 4.1 5.2 5.3 5.4 4.5 4.9 4.0 5.3 5.0 6.2
##  [307] 3.6 3.3 4.7 5.6 3.3 5.1 3.8 4.9 4.8 5.4 5.2 4.7 3.5 3.7 4.9 4.4 3.5 4.5
##  [325] 4.0 5.0 3.9 4.4 4.2 3.4 3.9 4.3 5.9 4.4 4.3 4.3 5.1 3.8 4.4 5.0 4.7 3.4
##  [343] 3.6 4.5 5.7 5.4 4.0 4.4 3.4 4.6 4.8 5.2 3.9 5.2 6.6 4.5 4.7 4.8 5.0 2.4
##  [361] 3.4 5.5 3.0 3.8 3.3 4.5 5.4 5.8 6.5 4.6 4.7 4.7 4.6 4.3 4.5 3.6 5.4 3.5
##  [379] 4.5 4.0 5.1 3.7 3.7 2.6 5.2 4.4 5.6 6.2 6.2 4.2 4.8 3.3 4.8 4.9 3.4 5.2
##  [397] 4.9 3.2 2.8 3.9 6.2 4.8 3.6 4.4 5.3 6.0 4.5 6.1 4.7 3.3 4.6 4.4 2.8 4.2
##  [415] 4.6 2.5 5.5 5.3 4.1 6.1 4.8 5.3 3.4 4.7 5.2 3.6 4.9 4.3 4.3 3.3 4.4 4.8
##  [433] 5.2 5.6 4.4 4.4 5.3 5.1 4.3 4.1 5.7 3.4 4.6 4.6 5.3 3.3 5.0 4.8 4.6 3.9
##  [451] 5.5 6.5 4.1 4.4 4.3 3.8 3.7 4.2 3.7 6.8 3.7 4.9 4.2 5.3 2.8 4.2 5.0 4.0
##  [469] 4.8 4.9 4.7 2.9 3.7 3.4 3.8 4.8 2.8 5.2 3.7 3.9 4.2 4.3 3.2 2.9 4.9 3.9
##  [487] 3.3 3.4 5.4 3.9 5.5 4.0 2.4 2.8 3.2 3.6 4.7 3.7 4.2 2.9 6.1 4.1 4.1 4.9
##  [505] 5.0 5.1 4.7 4.8 4.3 3.8 5.8 3.8 4.4 4.9 4.9 5.8 2.7 3.8 3.4 5.1 5.1 3.9
##  [523] 4.8 6.0 5.5 4.8 4.4 3.8 4.5 5.6 5.0 4.5 5.3 5.4 4.6 5.8 3.5 6.6 3.4 5.8
##  [541] 3.5 2.4 5.8 6.0 4.2 4.5 6.1 2.4 5.0 5.9 4.9 5.7 4.3 4.2 3.8 4.0 4.7 3.7
##  [559] 4.3 5.5 5.3 6.0 5.3 3.7 4.6 4.5 3.4 4.7 2.7 4.3 4.3 3.6 5.1 5.4 3.5 5.4
##  [577] 6.8 4.5 4.8 3.9 4.3 4.1 5.1 4.7 4.2 4.3 4.0 2.9 5.0 5.0 4.0 3.7 4.1 3.8
##  [595] 3.7 4.0 4.2 4.7 3.2 3.8 5.4 6.2 4.3 4.9 5.1 5.8 4.5 4.2 4.5 4.1 4.6 6.5
##  [613] 5.4 4.1 3.4 4.0 3.2 3.2 4.5 5.0 3.3 3.1 6.3 4.5 4.0 5.1 4.4 5.6 5.2 3.9
##  [631] 4.3 3.2 4.7 4.3 4.4 5.6 5.6 3.6 4.1 4.1 5.3 5.2 4.1 5.1 3.1 4.7 3.8 4.2
##  [649] 3.7 4.4 4.0 3.8 3.3 5.2 4.4 5.3 6.2 3.6 6.0 5.2 5.2 5.8 5.5 4.5 5.6 4.9
##  [667] 4.4 4.5 6.4 2.9 3.2 5.0 6.3 5.1 5.0 5.3 2.9 4.4 5.9 3.0 6.1 2.6 4.2 4.5
##  [685] 3.4 5.5 4.4 4.6 4.0 4.1 4.4 3.2 4.3 5.6 5.5 5.5 3.6 3.5 4.3 5.3 5.9 3.8
##  [703] 3.5 3.0 5.1 5.2 3.6 5.9 5.4 4.3 4.7 3.2 5.1 4.1 4.6 4.0 2.7 3.7 3.7 3.2
##  [721] 5.1 4.6 5.1 3.2 5.3 4.3 2.8 3.0 5.3 5.7 3.0 5.8 3.4 4.7 4.5 2.4 3.7 5.5
##  [739] 4.0 4.7 4.0 5.2 6.0 3.7 4.6 6.2 4.1 3.2 5.1 4.8 5.9 6.6 5.5 3.6 5.5 3.4
##  [757] 5.8 5.3 3.5 5.1 3.9 5.4 4.7 3.3 3.6 4.6 4.6 4.6 4.8 6.0 5.2 4.2 4.9 4.7
##  [775] 4.4 3.6 3.4 4.3 4.8 6.3 3.4 4.3 5.3 4.2 3.4 4.9 3.6 6.5 4.2 4.3 3.7 4.9
##  [793] 5.0 3.8 4.7 4.1 3.0 5.4 4.0 6.1 3.4 5.3 5.7 3.0 4.0 4.7 4.4 4.3 4.8 4.6
##  [811] 4.8 5.0 3.6 4.2 3.7 5.9 5.2 4.2 4.1 3.3 4.6 6.0 3.9 5.1 5.4 4.2 5.9 4.1
##  [829] 3.8 5.0 4.7 4.1 2.8 4.8 4.9 2.8 6.7 2.9 3.1 4.5 3.1 3.9 3.7 6.0 5.7 3.4
##  [847] 5.2 4.7 5.3 5.3 4.0 5.1 4.2 3.2 4.4 5.0 4.3 6.3 3.8 4.8 5.0 3.2 5.2 4.3
##  [865] 4.6 4.4 4.6 4.6 4.0 4.8 4.4 4.4 4.3 5.0 3.6 4.1 4.0 4.4 3.8 5.4 5.1 4.1
##  [883] 4.6 4.6 4.1 4.9 5.0 2.8 4.6 4.6 3.9 4.8 4.9 4.8 5.1 3.0 4.3 5.0 3.9 3.8
##  [901] 5.7 5.4 2.9 4.0 5.2 5.6 3.7 4.3 3.7 4.8 5.6 4.4 4.3 3.7 6.2 4.3 4.7 5.1
##  [919] 4.6 5.7 4.7 7.2 4.1 4.7 5.0 4.5 3.5 4.2 5.8 3.3 4.4 5.0 3.4 5.4 4.9 2.8
##  [937] 4.6 4.9 3.9 4.1 6.4 5.8 5.0 5.1 4.1 3.9 5.6 4.5 5.1 4.0 5.0 6.1 4.4 4.0
##  [955] 4.5 5.4 3.7 3.5 4.1 4.8 4.9 3.8 5.3 3.7 5.6 5.4 3.5 5.1 5.4 4.0 5.1 5.1
##  [973] 4.0 5.4 5.6 4.9 4.5 6.0 4.3 4.7 5.4 4.6 4.6 5.3 4.8 3.4 4.2 4.6 4.1 3.7
##  [991] 4.3 4.3 4.7 4.4 3.5 4.2 5.2 3.6 4.0 4.4
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.400 5.348 5.300 5.200 5.200 5.000 5.000 4.700 4.600 4.400
## 
## $jack.boot.se
## [1] 0.9762286
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
## [1] 1.999068
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
##   1.9022469   2.8351829 
##  (0.7874994) (1.3417094)
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
## [1] 0.07572954 0.15378948 1.21300888 1.09425744 1.61295011 0.89292409
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
##    [1]  0.2173485169  0.7574728493 -0.1544538586  0.6934694221  1.2701221561
##    [6]  1.1941064631 -0.8847716858  0.5966539534  0.1173296531  1.4125028635
##   [11]  1.2159414383  0.5917044103  0.0467397593  1.0246152326  0.6418289434
##   [16]  1.8178038243  1.0574448497  1.7971442413  0.1036083065  2.3567121479
##   [21] -0.1112617297  2.0049117674  2.1601227293  1.1847815850  2.2265983017
##   [26]  1.1952057517  0.7844095444  1.1350611094  1.3646680656  1.7467154195
##   [31]  1.2282391955  2.0607790688  1.3252594913  1.3441900409  0.0425335798
##   [36]  2.5251456731  0.7766279844  2.0834287608  2.2053133563  1.8767699214
##   [41] -0.6268326921  0.9392903516 -0.1555112686  1.2024578931  0.3692545648
##   [46]  2.3472644451  1.9738691538  0.0666202639  0.9692406715  2.5717829104
##   [51]  1.0896162629  2.3865946941  1.2906656549  0.6785488235  1.9366416778
##   [56] -0.6821730340  1.6848500096  0.7213747106  2.0372816325  0.2327128928
##   [61]  1.1359953426  2.0066074777  1.2659286164  2.0363029039  2.4156172370
##   [66]  0.5790690263  1.6873691945  0.0746690589  1.4158395001  1.3102494375
##   [71]  1.2506177836  1.6380820292  1.9496870904  1.4581482277  1.3847480177
##   [76]  0.2732571680  0.6036071631  0.5029722366  1.0110740616  1.1940850428
##   [81]  1.2497608969  1.6723812303  1.2250941863  2.0975272605  0.7931183112
##   [86]  1.4381357306  1.2190734425  0.0739557751  2.0380368803  0.7220340959
##   [91]  2.0337091785  0.2018764397  0.8986217617  1.2447160795  0.4020421376
##   [96]  0.9794524439  1.4075406511  0.4617331714  1.8242530217  2.2291814604
##  [101]  2.3190548027  0.3356523407  0.8551923489  2.4420785840  0.2609134080
##  [106]  0.5064704273  0.8418110753  0.3638437691  2.2215820038  0.6158621281
##  [111]  2.0925951061  0.7772203522  0.2461057497 -0.1756442736  1.0999489569
##  [116]  2.0432077686  1.2545015339  0.6215425741  0.4420263275  0.6895753922
##  [121]  0.3542927239 -1.4540725669  1.3284960777 -1.0755416344  2.2381684948
##  [126]  1.1841248434  1.9994182398  0.3022918070  1.6380820292  0.2436844986
##  [131]  2.3124762941  1.1834352034  1.1677303594  0.1730980038  2.3768788085
##  [136]  0.0515896291  2.3318697869  1.2012641085  2.1896241097  1.2218035794
##  [141]  1.0730472314  1.7490497479  2.3768204655  0.6747497434  1.9913245508
##  [146]  2.3191290014  1.9103110632  1.3374577694  2.3497877452  2.0523838395
##  [151]  0.4940244471  0.4663530310  1.2767994868  2.0116418171  0.4359833601
##  [156]  2.5097277747  1.9729515036  1.7573986259  1.2542817124  0.0158913484
##  [161]  0.1374347343  1.2431960484  0.3999492277  1.9332623249  1.2017639951
##  [166]  1.0410702683  1.8173373488  0.7833826101  2.4146142229  1.9893976068
##  [171]  2.1043249417  1.4527040968  1.0153470221  0.7898333559  2.3431234089
##  [176]  2.5063289494  1.2909785643  2.0389642660  2.1027066566  0.5472292724
##  [181]  2.1119864594  0.1499864565  1.3000212790  1.9839154141  2.4024612303
##  [186]  1.8272827261  1.0986630362  2.4066173474  1.3660790262  0.7039034150
##  [191]  2.2616164738  1.9540795217  2.0480154772  2.4047812610  0.0797487791
##  [196] -0.9509495205  1.3466775778  1.5619120260  1.1940850428  2.3247113768
##  [201]  1.4593788512  0.2042407331  0.2804521920  0.6389952925  0.3139187400
##  [206]  1.0027477681  1.4627023745  2.3154694128  0.7974380096  1.6882220390
##  [211]  1.7413267989  1.3579579506  2.3266707610  1.4133226173  1.9806275155
##  [216]  1.9537517368  0.7871570533  2.3687265679  1.2840624552  2.4734267113
##  [221]  1.8583937255 -0.4312204576  1.2605379682  1.9000552528  0.4613690157
##  [226]  1.2221751057  2.0363029039  2.3194593459  1.8058245239  2.3966909874
##  [231]  2.0377200072  0.3335084530 -0.8497176158  2.1614746122  0.4977565779
##  [236]  0.6987947646  1.1617148394  1.5002971273  0.3435693903  0.4166470337
##  [241]  2.2209908225  1.7666909018  1.9399923727  0.2845704143  1.0164151858
##  [246] -1.3367608703  1.2251966132  2.4317756851  1.1400181585  0.5469775097
##  [251] -0.3917376701  0.6762676113  0.7477988015  0.3559190694  1.7227779407
##  [256]  2.3076080365  1.2962991151  1.9255735566  1.1086471116  2.3454361382
##  [261]  0.3501127849  2.4029780347  1.8852140096  1.8267178062  2.1195934766
##  [266]  0.2663827307  0.4746592356  1.3653663346  1.9742381432  1.0607780319
##  [271]  0.6610518099  0.0839993073  1.0719365739  2.1402175366  1.3241911143
##  [276]  2.4048010933  1.6280207129  1.9868265620 -0.2695087872  0.8523445959
##  [281]  0.2881853451  1.7162284091  1.7067974722  0.0019635139  0.9841549706
##  [286]  2.0329802653  2.4018081257  2.1601029969  2.0602030165  0.3015091962
##  [291]  2.4387547134  1.7185759640  1.2619153048  0.1148779105 -0.0295300765
##  [296]  1.3113937842  2.2611766381  2.1840979187  1.7959884049  1.9140034504
##  [301]  1.8703204751  2.1061121390  0.2263907475  1.6303503676  1.3436458826
##  [306]  1.3789749441  0.5454609774  2.0354029745  1.7125260998  1.1679371822
##  [311]  2.0389250124  1.8265021566 -1.3160371965  0.1911502973  1.4179313272
##  [316]  0.9321162619  0.2617730983  1.4141447661  1.8511763835  1.4452823368
##  [321] -0.9760823792  1.7196204775  1.9903480133  2.3072897219  1.2080233322
##  [326]  0.9383002516  1.2210854720  2.4803261929  1.1346269404  1.4071486135
##  [331]  2.2141117653  2.0045452345  0.7631961474 -0.0094619835  2.3310637947
##  [336]  1.2694721863  0.4950049688  0.2610957241  1.2541687605  0.9683606487
##  [341]  1.7790470082 -0.7042376489  2.0092966218  0.6555325373  2.2329394583
##  [346]  1.4073452117  1.9924387779  2.3516174584  0.2050581302  2.0067462137
##  [351]  1.6799435284  2.3635751640 -0.9640156069  1.0101241326  2.3149253983
##  [356]  0.9023331704  1.0742861747  2.2590587432  2.2523647151  2.0181573777
##  [361]  0.3875614146  1.8124985172  2.2204372833 -0.6684515468  2.2893394286
##  [366]  2.2155281399  0.7228111789  1.9856819263  1.2730496845  2.1331597010
##  [371]  0.6784646623  1.1089261370  1.3416080084  0.5933480015  0.7519763595
##  [376]  2.4127032723  1.4620633210 -0.4524607992  1.2410302076  0.8305991130
##  [381] -0.0929771399  2.2470890363  0.1664148464  2.2398067442  2.4591666368
##  [386]  0.7075348048  1.2502307818  0.6202731242  1.7296257060  1.0387812714
##  [391]  1.0203036022  2.1004681406  0.5489424554  0.3373233463  1.0312775568
##  [396]  2.3267966337  1.3781492308  1.7018931633  1.7383760542 -0.2721877664
##  [401]  0.2716326644  0.7366345493  1.7666909018  1.6378213218  2.0963782031
##  [406]  1.0766292840  1.0461961655  0.9147741423  2.5111012730  1.8248000320
##  [411]  1.8044372392  1.3673956218  2.4543452841  0.7626932877 -0.0949185248
##  [416] -1.0502013770  2.5050047983  0.3989345952  1.3476322285  0.1581407482
##  [421]  0.8260545320  1.4224639305  1.5684200171  1.3445967094  2.3580757011
##  [426]  2.0793527736 -0.7829777832  1.1884022936  0.3596164996  2.1058943381
##  [431]  2.3904135313  2.1776929433  2.2196357790  1.2741700239  1.3385066149
##  [436] -0.2370863084  2.3971348816 -0.1268141071  1.6577519308  1.8481140107
##  [441]  2.0574501350  1.5217711288  2.1885415153  2.1981645022  1.8712911402
##  [446]  1.6722584727  1.6004385832  1.1844447054  2.3047353033  1.7592439496
##  [451]  0.5723924019  0.6385775786  0.2752221232  1.8999034064  1.2641791983
##  [456]  2.2045172374  1.9868232996  2.4239382333  2.3623644185  1.5784665146
##  [461]  1.4269327985  1.2630807367 -0.2433382328  1.2924860535  2.1885415153
##  [466]  2.5321071381  1.4171053988  1.3342798938 -0.0370411519  2.2104582367
##  [471]  0.4776021091 -0.0299772418  0.9887348457  2.3072146579  2.5296474467
##  [476]  1.4319345047  1.0031546815 -0.2936248740  0.6284942168  2.2014099252
##  [481]  0.8002447472  1.3935214485  0.5401873306  1.3762575054  0.0873284921
##  [486]  0.3402398048  2.3540243704  0.6688818696  2.1696553010  2.4548451383
##  [491]  0.6733863829 -0.0402756463  2.4780438041  1.6168965389  2.3417356486
##  [496]  2.2011938973  2.2737844401  1.4105759817  2.3963926692  0.8569612819
##  [501]  0.5987738569  2.1878081790  0.6912051924  2.2928678339  2.0538331320
##  [506]  0.7887193360  2.3644741882 -0.4301590936  1.3902177219  1.2430061452
##  [511]  2.1085356826  2.0071087812  1.1504739240  0.8225476408  1.2171287046
##  [516]  2.4768223284  0.3695845161  0.4944148482  0.5626616856  2.4449317814
##  [521]  2.0793214524  1.1178990629  2.1024411273  0.9839042090  1.7797189080
##  [526]  1.3460355915  2.0222088213  2.4787574856  0.5528609058  0.8624671994
##  [531]  1.2837107040  2.3687175980  2.1058943381  0.5380073395  0.7931057360
##  [536]  2.0398413937  2.0572975776  2.1048663091  0.1997881292  1.9678659908
##  [541]  2.3957686185  0.5036830646  1.9292789800  1.9476575814  1.2089989534
##  [546]  2.1713779034  0.0010740625  0.8010980955  0.1241758737  2.2355932120
##  [551]  1.3761579705  2.2907541899  1.1140130368  2.3188357969  2.2119327466
##  [556]  0.1612755265  2.0150436214  0.6046723313  2.3807283046  2.1040788307
##  [561]  0.6675129742  2.3135440705  1.0635266517  2.0880630510  1.3005744706
##  [566]  0.6976087026  1.2801898964  2.3159243921  2.1398332241  2.3501275808
##  [571]  2.2918079726  0.5395092811  0.7161446615  1.3829616199  1.2289675391
##  [576]  2.4116601395  0.6789211828  1.6783339928  2.1301401141  0.4454994366
##  [581] -0.0124969616  1.2730435213  1.2749655952  0.3635790141  0.8747323230
##  [586]  0.3357809787  2.2750351688  1.7004267589  2.0181573777  0.8280098474
##  [591]  2.1516755600 -0.2832113838  2.2668003739  0.9270924138  1.2144719655
##  [596]  1.2711066261  1.0348905948  2.4354880606  0.7174140380  1.9894815794
##  [601]  0.6604496271  1.7280347272  0.4458276724  2.0260235756  1.7934213839
##  [606] -0.1464720957  2.0037859710  1.0667977853  1.6441819908  0.7326352691
##  [611]  2.0544046667  0.0629052096  2.0464359717 -0.7645323613  1.0793190039
##  [616]  2.3181041761  1.3937541075  1.2173331257 -0.2954790477 -0.2494123087
##  [621]  2.2638968332  0.5536068584  2.0427633382  1.1664431552 -0.0856706996
##  [626]  0.6626808646  2.2529694217  1.2444107756  0.5509798348 -0.3109702558
##  [631]  0.6415271422  0.5459523129  2.2604167034  2.2217170353  1.3735102297
##  [636]  2.0890339733  0.4609163157  2.3738469700  1.3372268614  1.3877215360
##  [641]  2.1166692983  2.4387762356  1.3068020275  1.9235483529  2.4853366856
##  [646]  1.6962038294  0.6823780246  2.1038140630 -0.2021245491  1.7335131628
##  [651]  1.3868571657  2.2887624038  2.4342273330 -0.2675275522  2.1453081919
##  [656]  1.4307319515  1.2482484971  1.2476624896  0.3496843981  1.7614768160
##  [661]  1.7186013892  0.8605403383  2.0403382353 -0.0008190266  1.3545711942
##  [666]  0.4135659705  2.2715932607  1.0620997133  0.9793131555  0.6155541307
##  [671]  1.0696367429 -0.2991812099  0.8008962675  1.3849979106  1.6991883655
##  [676]  1.2234569075  2.1092764946  2.5356105246  2.2820553562  0.7469034453
##  [681]  0.5324754199  1.7218379757  1.9971522504 -0.1846423839  0.4667625065
##  [686]  1.8450967462  2.2764995483  1.2478034316  2.1556912698  1.3759854017
##  [691] -0.7891468260  1.0230239873  1.2045242226  0.1088873560  1.0734227367
##  [696]  0.6385851821  1.3079403179  1.7344847965  2.1395400702  2.0650862954
##  [701] -0.7025465946  0.3680642906  1.1986265198  1.3496287547  0.6732110405
##  [706]  1.3130284764  1.7761024309  1.8374084444  2.0742692544  1.3622239087
##  [711]  1.6984646837  1.0093469629  2.4342273330  2.0507286654  1.4893509695
##  [716]  0.6275146361  0.2673832169  1.4106521818  1.2877265389  2.4631873514
##  [721]  1.8914221899  1.0193922660 -0.6011160077  0.2052185927  1.1003398371
##  [726]  0.8789937027  2.4851350127  1.9934958896  1.5351627241 -0.0195051508
##  [731]  1.3639022833  2.2685725397  1.3720596681  2.3627341245  1.0091461192
##  [736]  0.3487333160  1.7499150914  0.1529320296  0.7992714806  1.5470113149
##  [741]  2.0023544963  0.6108170400  2.4403341995  1.9813804365  1.2336020322
##  [746]  0.8124595282  0.4974619577  2.3423037793  1.6563406912  0.8258418650
##  [751]  2.0284004656  0.5723962110  0.0765688972  0.0362946777  0.3013731982
##  [756] -0.4348081846  0.3139187400  2.2347734589  1.0188223548  0.5309995126
##  [761] -0.1783791343  0.4111871754  0.5062074936  2.3484560300  0.6690164055
##  [766]  2.1154327939  2.4493598961  1.8919606117 -0.3455922388  0.6924753995
##  [771]  2.0469946011  1.1758955668  2.2869066342  0.2838705092 -0.2146767277
##  [776]  1.2035598395  2.4691527922  0.6545313465  0.7387971344 -0.4700374327
##  [781]  1.7790976144  0.5105411550  2.0020591409  0.8213847300  2.5136819267
##  [786]  2.3662965917  2.3082723970  2.0257850412  1.3652045460  2.4731522379
##  [791]  2.2701121942  2.2501488636  2.3095215460  2.0329748073 -0.5607186567
##  [796]  0.6732717084  0.9120178818  1.4125300845  1.3637060911  2.3969655846
##  [801]  2.0087497589  2.1453921078  1.9738541099  0.8818524050  1.7345754971
##  [806]  0.7618845227  2.4402666756  2.0171402910 -0.2005873683  2.3216151554
##  [811]  0.6912051924  0.0541367770  0.4625438084  1.4268421778  1.9676477343
##  [816]  1.0638427272  2.1266379948  2.1020627436  2.2556794584  0.7900534176
##  [821]  1.2224158983  2.3754857365  0.8885515351  0.8149858559  0.3124797149
##  [826]  2.0039535848  1.9443453895  0.8411741216  2.3371015012  0.6732578643
##  [831]  1.6394958645  0.0390270147  0.6944169100  0.1397178554  0.5789466807
##  [836]  1.3180736752 -0.2922457129  1.2105217081  1.3869950844  2.0397736801
##  [841]  1.6394675673  2.5117105978  2.4268278113  0.6031030060  1.9083397180
##  [846]  2.2827884838  0.4655503492  1.7616001093  1.6873973190  0.4897341084
##  [851]  1.6565398305  1.6673938621  0.5301933647  0.2771853892 -0.7257213114
##  [856]  1.3497815446  1.9480942071  1.7748657327  2.1363820348  2.4725515899
##  [861]  0.9690135586  2.2339536086  2.0465705220  0.5669030841  0.1356496799
##  [866]  2.1307419459 -0.0278275193  0.6750085066  0.3338820437  0.9485287906
##  [871]  2.2505302753  2.4844833137  1.2350820476  1.2535972378  2.0308465549
##  [876]  1.4142108060  1.2005625592  0.8198809086  1.3453878241  1.2813262812
##  [881]  2.3199940002  0.1760393742  1.3015454935  0.4055821397  2.1604050093
##  [886]  0.3563063883  0.4990473462  2.4897196584 -0.5379351116  0.3510287546
##  [891]  1.4108933945  1.9035403638  2.3247858307  2.2985953148  0.5028632330
##  [896]  1.0081726978 -0.4913030003  2.3820906927  1.3870474658  1.8840756706
##  [901]  2.3812750494  1.2265362444  0.9545604822  1.3482021879  0.0197436195
##  [906] -0.3167372641  2.4387762356  0.9264304660  0.1051795448  1.6818473751
##  [911]  1.4434964106  1.1870314573  0.0162780246  1.2064933073  0.8270380520
##  [916]  0.9833005749  1.9600142710  2.3361348751  1.0472740445  0.8396297768
##  [921] -0.2264039847 -0.4107694861  0.1559032419  0.6737409614  1.8327926047
##  [926]  1.3927295235  1.1797606926  2.3052591407  0.5109145214  0.8460386691
##  [931]  2.1095253816  2.3559840658  0.7833598936  1.1572170267 -0.2430496205
##  [936] -0.8299484493  1.3810622538  1.9678659908  0.7800191614  1.8468292354
##  [941]  0.9643334956  2.3362884432 -0.0407862327  0.4043176018  0.4901246482
##  [946]  2.5481238217  0.5902487966  1.1039747829  0.7129669221  2.0677131522
##  [951]  0.5137309632 -0.4453542678  2.0248700622 -0.4965384751  0.7231478159
##  [956]  1.2248290338  1.0117413485  0.3384675230  1.1802494034  0.6904844431
##  [961]  2.0479575189  0.9900168079  2.0354029745  0.2613725970  1.8644430669
##  [966]  0.2540874155  1.9965247163 -0.8933672356  0.6545981216  0.8260545320
##  [971]  0.7718998935  2.1900763192  2.4047221423  0.2201036823  2.1748619179
##  [976]  1.2514960415  0.9030682927  0.8449416550  1.2101191261  0.1838958682
##  [981]  2.0559404287  1.2689127026  0.2636058629  0.3644599069  0.0429520972
##  [986]  1.0472740445  2.1653955205  1.9434632914  1.9107819498  2.1231063445
##  [991]  1.0800667055  0.7307418836 -0.4992255309  0.3001402435  2.0087497589
##  [996]  2.4491401419  0.7591295004  2.1306360711  2.4236349985  1.6149912857
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
##   0.6709436   0.5958446 
##  (0.1884226) (0.1332329)
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
## [1]  0.08282930 -0.94374956  0.49940693  0.04092662  0.01479075  0.60415601
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
## [1] -0.0287
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9023084
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
##     original       bias    std. error
## t1*      4.5 -0.006106106   0.9241433
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 1 2 6 7 8 9 
## 2 1 2 1 1 1 2
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
## [1] 0.0276
```

```r
se.boot
```

```
## [1] 0.9322593
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

