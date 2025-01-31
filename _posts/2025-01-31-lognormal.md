---
title: The Lognormal Distribution
author: Tarik Benaddi
category: Statistics
layout: post
published: true
---

## Introduction

The Normal distribution is often presented as a universal approximation for a wide range of phenomena, from the height of the population in a country to the weight of apples in a bag, and even income distribution or asset prices in stock markets. In wireless communications, the Normal distribution serves as the first approximation for various noise sources.

However, a closer look reveals a fundamental inconsistency: the domain of a Normal random variable spans from $$-\infty$$ to $$+\infty$$. Sure the probabilties at the boundaries get very small, still the mentionned real-world quantities, such as height or stocks prices, cannot never be negative…

XXX figure normal

<!-- excerpt_separator -->

In statistics, there exists a less popular but equally powerful distribution that addresses this issue: the Lognormal distribution. Unlike the Normal distribution, the Lognormal ensures that all values are non-negative, making it a better first approximation for data where negative values are physically impossible.

Back to wireless communications, while additive Gaussian noise is widely used due to its simplicity and mathematical tractability (thanks to nice properties like the independence of mean and variance, symmetry, the Gaussian integral, Gaussian stability with respect to  convolution/product/Fourier transform…, …), it doesn’t always align with observations. And the use of Additive Gaussian noise to model this noise is often justified by the Central Limit Theorem (CLT)[^1], which states that the sum of many independent, identically distributed random variables tends to follow a Normal distribution.

However, this additive assumption doesn’t hold in all systems. In some cases, noise is not additive but rather multiplicative. Take free-space optical communications as an example. Here, atmospheric turbulence causes by scintillation results in a multiplicative effects on the signal intensity. The received signal power, say $$q(t)$$, is the product of the transmitted power $$p(t)$$ and a multiplicative noise term $$n(t)$$, which itself may be the product of several independent noise factors: $$n(t) = n_1(t) \times n_2(t) \times \dots$$

In such cases, we would like to have an equivalent of the Central Limit Theorem (CLT) that applies to the product of a large number of random variables instead of their sum. It turns out that deriving such theorem is not straightforward and is only addressed in specific cases in the literature. The Lognormal distribution hence emerges as an interesting workaround for our case.


## Construction
Let us consider the signal model $$q = y \cdot p$$, where the random variable $$Y = \prod_{i=1}^n Y_i$$, and $$Y_i$$ are independent, identically distributed (i.i.d.) strictly positive random variables (time dependency is omitted for notational simplicity). 

The key idea behind deriving the Lognormal distribution lies in transforming the multiplicative noise model into an additive one by moving to the logarithmic domain. By taking the natural logarithm of both sides, we obtain:

$$\ln(q) = \sum_{i=1}^n \ln(Y_i) + \ln(p)$$

Of course this holds for strictly positive values, think of signal power levels for example or stock prices. This transformation gives another signal model where the new signal of interest, $$\ln(p)$$, is corrupted by an additive noise equal to $$\sum_{i=1}^n \ln(Y_i)$$.

Since $$\{\ln(Y_i)\}_i$$ are i.i.d. random variables, by classical CLT again, the summation tends to follow a Normal distribution as $$n$$ goes to infinity. Hence we say that $$X$$, $$(\triangleq e^{Y})$$, is Lognormally distributed.

With an abuse of notation,  we can write $$\log\mathcal{N}(\mu, \sigma^2) \triangleq e^{\mathcal{N}(\mu, \sigma^2)}$$, which is quite misleading: the Lognormal distribution is actually taking the exponential of a normal distribution, "Lognormal" should be understood as “Normal in the sense of the log” and not “the log of a normal distirbution”.

## Probability Density Function
Let $$Y \sim \mathcal{N}(\mu, \sigma^2)$$, and we seek the probability density function (PDF) of $$X = e^Y$$. The PDF of $$Y$$ is given by

$$f_Y(y) = \frac{1}{\sigma\sqrt{2\pi}} \exp\left(-\frac{(y - \mu)^2}{2\sigma^2}\right)$$

Now consider the transformation $$h: x \mapsto e^x$$ and its inverse $$h^{-1}: x \mapsto \ln(x)$$. This gives

$$X = h(Y) = e^Y \text{ and }Y = h^{-1}(X) = \ln(X)$$

Since $$X = h(Y)$$ is a continuous random variable and $$h^{-1}$$ is strictly monotone, the transformation theorem[^2]  yields:

$$f_X(x) = f_Y(h^{-1}(x)) \cdot \left|\frac{dh^{-1}(x)}{dx}\right| \quad \text{for } x \in \mathbb{R}^{+*}$$ 

After simplification, we obtain the PDF of the Lognormal distribution as:

$$f_X(x) = \frac{1}{x\sigma\sqrt{2\pi}} \exp\left(-\frac{(\ln(x) - \mu)^2}{2\sigma^2}\right), \quad x > 0$$

Some observations here:
1. it has the same form as the PDF of the Normal distribution, but with $$x$$ replaced by $$\ln(x)$$ and the whole scaled by $$x$$. 
2. the PDF of the Lognormal distribution is not symmetric, unlike the Normal distribution: since $$X$$ is always positive, it cannot exhibit the symmetry characteristic of the Normal distribution.
3. Finally, notice that the support of a Lognormal distribution is $$]0,+\infty[$$. Hence $$c+Y, c\in\mathbb{R}$$ can't be a lognormal.

## Cumulative distribution function

The cumulative distribution function (CDF) of the Normal distribution is given by:

$$F_Y(y) = \int_{-\infty}^y \mathcal{N}(t; \mu, \sigma^2) \, dt = \frac{1}{2} \left[1 + \text{erf}\left(\frac{y - \mu}{\sqrt{2}\sigma}\right)\right]$$

where $$\text{erf}(x) = \frac{2}{\sqrt{\pi}} \int_0^x \exp(-t^2) \, dt$$

We can compute its CDF, following a similar derivation and variable substitution for the Lognormal distribution. Given the PDF of the Lognormal distribution:

$$f_X(x) = \frac{1}{x\sigma\sqrt{2\pi}} \exp\left(-\frac{(\ln(x) - \mu)^2}{2\sigma^2}\right), \quad x > 0$$

Using the substitution $$t' = \ln(t)$$, we have $$dt' = \frac{1}{t} dt$$, and the CDF integral transforms as follows:

$$F_X(x) = \int_{-\infty}^{\ln(x)} \frac{1}{\sigma\sqrt{2\pi}} \exp\left(-\frac{(y - \mu)^2}{2\sigma^2}\right) \, dy$$

This integral is identical in form to the CDF of the Normal distribution, but with the upper limit $$\ln(x)$$ instead of $$x$$. Therefore, the CDF of the Lognormal distribution can be simplified to:

$$F_X(x) = \frac{1}{2} \left[1 + \text{erf}\left(\frac{\ln(x) - \mu}{\sqrt{2}\sigma}\right)\right]$$

At first glance, this result may seem counterintuitive, as it appears to simply substitute $$x$$ with $$\ln(x)$$ in the Normal CDF. But actually becomes clear since $$P(X \leq x)$$ is equivalent to $$P(Y \leq \ln(x))$$. This connection makes another more straightforward derivation of the CDF without resorting to the PDF.

## Some Moments
Using the same integration techniques as in the Normal distribution case, we can derive key properties of the Lognormal distribution:

- Expectation:
    - For $$Y \sim \mathcal{N}(\mu, \sigma^2)$$, the expectation is $$E[Y] = \mu$$

    - For $$X = e^Y$$, the expectation is:
    
    $$E[X] = \int_0^\infty x f_X(x) \, dx = \exp\left(\mu + \frac{1}{2}\sigma^2\right)$$
- Variance:
    - For $$Y$$, the variance is $$\text{Var}(Y) = \sigma^2$$.
    
    - For $$X$$, the variance is:
    
    $$\text{Var}(X) = E[(X - E(X))^2] = e^{2\mu}\left(e^{2\sigma^2} - e^{\sigma^2}\right)$$

- Median:
    - For $$Y$$, the median is $$\text{med}(Y) = \mu$$.
    
    - For $$X$$, the median is $$\text{med}(X) = e^\mu$$
    
- Mode:
    - For $$Y$$, the mode is  $$\text{mod}(Y) = \mu$$.
    
    - For $$X$$, the mode is $$\text{mod}(X) = e^{\mu} e^{-\sigma^2}$$

We see that the mean and variance of the Lognormal distribution depend on both $$\mu$$ and $$\sigma^2$$ of the underlying Normal distribution. Also, unlike the Normal distribution where the mean and the median are equaln the Lognormal distirbution is asymmetry. Which causes a skew to the right, with a long tail extending to the left toward larger values.

XXX graph lognormal with mode median etc

## Interpretation of $$e^\mu$$ and $$e^{\sigma^2}$$

We saw that where $$\mu$$ is the median $$\mathcal{N}(\mu, \sigma^2)$$, $$e^\mu$$ is the median of $$Y \sim \ln\mathcal{N}(\mu,\sigma^2)$$. However, the variance of $$\ln\mathcal{N}(\mu,\sigma^2)$$ is not $$e^{\sigma^2}$$. What may this quantity represent? It turns out that $$e^{\sigma^2}$$ does not have an easy interpretation, but still it is related to the scale/dispersion of the distribution. Actually, the Lognormal distribution has a natural connection to geometric measures due to its multiplicative nature. Specifically:

* The geometric expectation of $$X$$ is nothing but the exponential of the arithmetic expectation of $$Y$$ since:

$$\text{GE}(X) = \left(\prod_i x_i\right)^{\frac{1}{n}} = \left(\prod_i e^{y_i}\right)^{\frac{1}{n}} = \exp\left(\frac{1}{n}\sum_i y_i\right) = e^\mu$$

* Same for the geometric variance of $$X$$ where $$\text{GV}(X) = e^{\sigma^2}$$.

The geometric expectation is particularly useful in contexts where the underlying process is multiplicative rather than additive. For example, in finance, the geometric mean is used to calculate average rates of return over time, since investment returns compound exponentially.

Another interesting observation that exhibits the quantity $$e^{\sigma^2}$$ is the ratio between the arithmetic mean and the geometric mean (AM-GM ratio) of $$X$$:

$$\left(\frac{\text{AM}(X)}{\text{GM}(X)}\right)^2 = \left(\frac{\exp(\mu + \frac{\sigma^2}{2})}{\exp(\mu)}\right)^2 = e^{\sigma^2}$$

A slightly different ratio is interesting when computing the so-called coefficient of variation[^3] ($$C$$):

$$C=\frac{\sqrt{Var(X)}}{E[X]}=\sqrt{e^{\sigma^2}-1}$$ 

$$C$$ is a valuable metric as it is more context agnostic in comparison to the standard deviation: unlike the standard deviation, the CV is a dimensionless number and this property makes it particularly useful for comparing variability across different data sets with varying units or significantly different means.

## Moment generating function
The moments of a distribution can be derived using the moment generating function (MGF), defined as:

$$M_Y(t) = E[e^{tY}] = E[X^t] = \exp\left( \mu t + \frac{1}{2}\sigma^2 t^2\right) , \quad t \in \mathbb{N}$$

An interesting propoerty of the Lognormal distribution is that while it has finite moments of all orders $$(M_Y(t) \in \mathbb{R}, \forall t \in \mathbb{N})$$, its MGF is not defined elsewhere.


## Entropy

Computing the entropy, or “surprise meter”, is of interest in statistics as it quantifies the amount of information we get from sampling a PDF. The differential entropy of a random variable $$X$$ is defined as:

$$h(X) = -E[\ln f_X]$$

As we have seen, for a Lognormal distribution, the PDF is given by:

$$f_X(x) = \frac{1}{x\sigma\sqrt{2\pi}} \exp\left(-\frac{(\ln x - \mu)^2}{2\sigma^2}\right)$$

Taking the logarithm of $$f_X(x)$$ and after simplification:

$$\ln f_X(x) = -\frac{1}{2}\ln(2\pi\sigma^2) - \ln(x) - \frac{(\ln x - \mu)^2}{2\sigma^2}$$

Using this and the fact that $$\log(X) \sim \mathcal{N}(\mu, \sigma^2)$$ we obtain:

$$h(X) = -E[\ln f_X(x)] = E\left[\frac{1}{2}\ln(2\pi\sigma^2) + \ln(x) + \frac{(\ln x - \mu)^2}{2\sigma^2}\right]=\frac{1}{2}\ln(2e\pi\sigma^2) + \mu$$

Observe that the entropy of the Lognormal distribution is a translated version of the entropy of the Normal distribution by the mean $$\mu$$.

The entropy of the Lognormal distribution depends on $$\mu$$ since the widening of its PDF, and hence the amount of "surprise" of the Lognormal distribution, varies with $$\mu$$. At the contrary, for the normal distribution, the $$\mu$$ just translates the distribution as it is without any change in probabilities.

On the other hand, while the entropy of the Normal distribution depends only on the variance $$\sigma^2$$, for the Lognormal distribution, the entropy depends on both the mean $$\mu$$ and the variance $$\sigma^2$$ of the underlying Normal distribution.


[^1]: [https://en.wikipedia.org/wiki/Central_limit_theorem](https://en.wikipedia.org/wiki/Central_limit_theorem)
[^2]: [https://www.cl.cam.ac.uk/teaching/2002/Probability/prob11.pdf](https://www.cl.cam.ac.uk/teaching/2002/Probability/prob11.pdf)
[^3]: [https://en.wikipedia.org/wiki/Coefficient_of_variation](https://en.wikipedia.org/wiki/Coefficient_of_variation)