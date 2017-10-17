## these are function to analyze the power to detect QTL and to estimate the true number of loci in QTL mapping studies using intercross data
## lynch_walsh_sample calculates the minimum F2 sample size needed to detect a QTL of a given effect size under different type 1 and type 2 errors (equation 15.37, Lynch & Wlash, 1998, Genetics and Analysis of Quantitative Traits)
## lynch_walsh_detect approximates the smallest effect size one can detect in a QTL experiment given the F2 sample size based on Lynch & Walsh, 1998
## otto_jones_detect approximates the smalles effect size one can detect (threshold theta) based on equation 11 in Otto & Jones, 2000, Detecting the undetected: estimating the total number of loci underlying a quantitative trait, Genetics 156.
## otto_jones_trueQTLnumber calculates the expected true number of loci underlying a trait given the approximated theta, the minimum and average effect size, the parental difference delta z, and the number of detected QTL. Note that these values can be both standardized (delta Z = 0.5, and minimum effect = a/parental_difference) and on the observed scale (latter is preferred)
## I recommend trying the examples in Example 11 in Chapeter 15 of Lynch and Walsh to get a feel for the first two functions and the example in Table 5 in the Otto & Jones paper to get a feel for how the latter two functions work.

lynch_walsh_sample<-function(effect.size=0.5, alpha=0.05, beta=0.1, k=0) {
#effect.size = expected effect size of QTL
# alpha = type I error, false positive rate
# beta = typ II error, 1-expected power to detect a QTL
# k = dominant/recessive effect; 0 means fully additive, positive implies dominance effects, negative recessive effects
 
	 prob.beta=qnorm(1-beta, mean = 0, sd = 1, log = FALSE)
	 prob.alpha=qnorm(1-(0.5*alpha), mean = 0, sd = 1, log = FALSE)
	 sample.size=((1-effect.size)/effect.size)*(((prob.alpha/(sqrt(1-effect.size)))+prob.beta)^2)*(1+((k^2)/2))
	 sample.size
	 }
	 
lynch_walsh_detect<-function(sample.size=100, alpha=0.05, beta=0.1, k=0,res=0.001) {
#sample.size=sample size of QTL experiment
#res is resolution at which to estimate the minimum effect size of detectable QTL

	prob.beta=qnorm(1-beta, mean = 0, sd = 1, log = FALSE)
	prob.alpha=qnorm(1-(0.5*alpha), mean = 0, sd = 1, log = FALSE)
	vector.sample.size<-data.frame(effectsize=1, samplesize=lynch_walsh_sample(1,alpha,beta,k))
	for(i in seq(from=1-res, to=0+res, by=-res)) {
		effect.size=i
		vector.sample.size=rbind(vector.sample.size,data.frame(effectsize=i, samplesize=lynch_walsh_sample(i,alpha,beta,k)))
		}
	detect=min(vector.sample.size[which(vector.sample.size$samplesize<sample.size),"effectsize"])
	detect
	}



otto_detect<-function(amin=0.03,M=0.09,nd=7, res=4) {
#theta = detectability threshold	
#D = sum of additive effects across all possible QTL, ie half the parental difference on the original scale (in which case amin should also be on the scale of the measured variable) or standardized ( in which case, always 0.5)
#amin = minimum effect size that was detected in a QTL experiment
#nd = number of detected QTL
#M = average effect size of detected QTL
#res = the resolution used to approximate zero, the larger the more precise, but also the more computation time involved. 

equation11=function(theta){amin-theta-((M-theta)/nd)+((amin*exp((-1*amin*nd)/(M-theta)))/(1-exp((-1*amin*nd)/(M-theta))))}
	repeat{
		theta=sample(seq(from=0, to=amin, by=1*10^(-1*res)),1)
		solution=round(equation11(theta),res)
		if(solution == 0) break
		}
	return(theta)
	}
	
otto_trueQTLnumber<-function(D=0.83,amin=0.03,M=0.09,nd=7, res=4) {
	nqtl=D/(M-otto_detect(amin,M,nd,res))
	return(nqtl)
	}
	