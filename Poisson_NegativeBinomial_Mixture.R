#OK let us start with some minor observations of a gaussian distribution with different means
a = rpois(10,1)
b = rnbinom(100,size = 2,mu = 15)
data = c(a[1:5],b[1:30],a[6:10],b[40:100])
plot(data, type="l")
#we can see that we have mixture of two distributions
#Let us initialise our Mixture EM data for every state
s11m = s12m = mean(data)-sd(data)
s21m = s22m = mean(data)+sd(data)
s11s = s12s = 1
s21s = s22s = 2

#Let us now initialise the state probability i.e. P(Si)
S1 = 0.5
S2 = 0.5
#Transition probability self, then the other
TS1 = c(0.5, 0.5)
TS2 = c(0.5, 0.5)
S11weights = runif(1)
S12weights = 1-S11weights
S21weights = runif(1)
S22weights = 1-S21weights

#Emission probability is calculated as in the EM approach, this makes sense that we have a probability for every observation, since they are not the same observation
#Except that we don't have a mixture in the states, essentially we want every state to map to one component of the mixture
bs1w1 = dpois(data,s11m)
bs1w2 = dnbinom(data,size = s12s , mu = s12m)
bs2w1 = dpois(data,s21m)
bs2w2 = dnbinom(data,size = s22s , mu = s22m)

bs1 = bs1w1*S11weights + bs1w2*S12weights
bs2 = bs2w1*S21weights + bs2w2*S22weights
  
for(x in 1:1000)
{
#OK, now we are set to do the first HMM iteration
#Forward probability
#As mentioned above of the observation that we have, we will calculate the forward probability 
alpha = matrix(ncol=2,nrow=length(data))
#Forward for state 1 and first of the observation
alpha[1,1] = S1 * bs1[1] #Probability of State1 * Probability of observing the first number
#Forward for state 2 and first of the observation
alpha[1,2] = S2 * bs2[1] #Probability of State1 * Probability of observing the first number

#Second is just conditional on the previous one (could be either s1 or s2) and ENDING in state 1
alpha[2,1] = alpha[1,1] * TS1[1] *  bs1[2] + # ending in s1, P(S1) * Self transition * P(O2|S1)
            alpha[1,2] * TS2[2] *  bs1[2] #ending in s1, P(S2) * T[S2->S1] * P(O2|S1)

alpha[2,2] = alpha[1,2] * TS2[1] *  bs2[2] + # ending in s2, P(S2) * Self transition * P(O2|S2)
             alpha[1,1] * TS1[2] *  bs2[2] #ending in s2, P(S1) * T[S2->S1] * P(O2|S1)

for(i in 2:length(data))
{
  #Ok, let us see the pattern
  #Probability of Observation given state 1 * sum of ((probability of previous state1 * self probability)+
  #(Probability of previous state2 * transition from S2 to S1))
  alpha[i,1] = bs1[i] * (alpha[i-1,1] * TS1[1] + alpha[i-1,2] * TS2[2])
  #Similarly
  alpha[i,2] = bs2[i] * (alpha[i-1,2] * TS2[1] + alpha[i-1,1] * TS1[2])
}

#Probability of the whole observation across the two states
po = alpha[length(data),1]+alpha[length(data),2]

#OK now the backward probability

beta = matrix(ncol=2,nrow=length(data))
beta[length(data),1] = 1
beta[length(data),2] = 1
#Stand at the previous step and then look ahead
#Probility at length(data)-1 at state 1 = T[S1->S1] * P(O|S1) * B(length(data),S1) + T[S1->S2] * P(O|S2) * B(length(data),S2)
beta[length(data)-1,1] = TS1[1] * bs1[length(data)-1] * beta[length(data),1] + TS1[2] * bs2[length(data)-1] * beta[length(data),2]
beta[length(data)-1,2] = TS2[1] * bs2[length(data)-1] * beta[length(data),2] + TS2[2] * bs1[length(data)-1] * beta[length(data),1]

#pattern
for(j in (length(data)-1):1)
{
  beta[j,1] = TS1[1] * bs1[j+1] * beta[j+1,1] + TS1[2] * bs2[j+1] * beta[j+1,2]
  beta[j,2] = TS2[1] * bs2[j+1] * beta[j+1,2] + TS2[2] * bs1[j+1] * beta[j+1,1]
}

#Gamma, which is the probability of the state given the observation
gamma = matrix(ncol=2, nrow=length(data))
gamma[,1] = alpha[,1]*beta[,1]
gamma[,2] = alpha[,2]*beta[,2]
gamma = t(apply(gamma,1,function(x){x/sum(x)}))

#Eta, which is the transition probability from state i to j at data length t
eta1 = matrix(ncol=2,nrow=(length(data)-1)) #for every data observation
eta2 = matrix(ncol=2,nrow=(length(data)-1)) #for every data observation
#going from 1 to 1 across the first observation
eta1[1,1] = (gamma[1,1] * TS1[1] * bs1[2] * beta[2,1])/beta[1,1]
#going from 1 to 2 across the first observation
eta1[1,2] = (gamma[1,1] * TS1[2] * bs2[2] * beta[2,2])/beta[1,1]
#going from 2 to 2 across the first observation
eta2[1,1] = (gamma[1,2] * TS2[1] * bs2[2] * beta[2,2])/beta[1,2]
#going from 1 to 2 across the first observation
eta2[1,2] = (gamma[1,2] * TS2[2] * bs1[2] * beta[2,1])/beta[1,2]

#Pattern

for(t in 1:(length(data)-1))
{
  #going from 1 to 1 across the first observation
  eta1[t,1] = (gamma[t,1] * TS1[1] * bs1[(t+1)] * beta[(t+1),1])/beta[t,1]
  #going from 1 to 2 across the first observation
  eta1[t,2] = (gamma[t,1] * TS1[2] * bs2[(t+1)] * beta[(t+1),2])/beta[t,1]
  #going from 2 to 2 across the first observation
  eta2[t,1] = (gamma[t,2] * TS2[1] * bs2[(t+1)] * beta[(t+1),2])/beta[t,2]
  #going from 1 to 2 across the first observation
  eta2[t,2] = (gamma[t,2] * TS2[2] * bs1[(t+1)] * beta[(t+1),1])/beta[t,2]
}

TS1[1] = sum(eta1[,1])/sum(gamma[1:(length(data)-1),1])
TS1[2] = sum(eta1[,2])/sum(gamma[1:(length(data)-1),1])
TS2[1] = sum(eta2[,1])/sum(gamma[1:(length(data)-1),2])
TS2[2] = sum(eta2[,2])/sum(gamma[1:(length(data)-1),2])

#if(gamma[1,1]==S1&&gamma[1,2]==S2) break;

S1 = gamma[1,1]
S2 = gamma[1,2]

#Mixture gamma
gammaM = matrix(ncol=4, nrow=length(data))
gammaM[,1] = gamma[,1]*S11weights*bs1w1/bs1
gammaM[,2] = gamma[,1]*S12weights*bs1w2/bs1
gammaM[,3] = gamma[,2]*S21weights*bs2w1/bs2
gammaM[,4] = gamma[,2]*S22weights*bs2w2/bs2

S11weights = sum(gammaM[,1])/sum(gamma[,1])
S12weights = sum(gammaM[,2])/sum(gamma[,1])
S21weights = sum(gammaM[,3])/sum(gamma[,2])
S22weights = sum(gammaM[,4])/sum(gamma[,2])

#Update the distribution parameters, mean and standard deviation similarly
s11m = sum(gammaM[,1]*data)/sum(gammaM[,1])
s12m = sum(gammaM[,2]*data)/sum(gammaM[,2])
s21m = sum(gammaM[,3]*data)/sum(gammaM[,3])
s22m = sum(gammaM[,4]*data)/sum(gammaM[,4])

#v11s = sum(gammaM[,1]*((data-s11m)*t(data-s11m)))/sum(gammaM[,1])
v12s = sum(gammaM[,2]*((data-s12m)*t(data-s12m)))/sum(gammaM[,2])
#v21s = sum(gammaM[,3]*((data-s21m)*t(data-s21m)))/sum(gammaM[,3])
v22s = sum(gammaM[,4]*((data-s22m)*t(data-s22m)))/sum(gammaM[,4])

be = c(v12s/s12m,v22s/s22m)
for(x in 1:length(be)) if(be[x] < 1) be[x]=2 

#s11s = s11m/(be[1]-1)
s12s = s12m/(be[1]-1)
#s21s = s21m/(be[3]-1)
s22s = s22m/(be[2]-1)

bs1w1 = dpois(data,s11m)
bs1w2 = dnbinom(data,size = s12s , mu = s12m)
bs2w1 = dpois(data,s21m)
bs2w2 = dnbinom(data,size = s22s , mu = s22m)

bs1 = bs1w1*S11weights + bs1w2*S12weights
bs2 = bs2w1*S21weights + bs2w2*S22weights

print(c(s11m,s12m,s21m,s22m))
print(c(s11s,s12s,s21s,s22s))
print(c(S11weights,S12weights,S21weights,S22weights))
print(c(S1,S2))
print(c(TS1,TS2))
}

#Now we can get the state model for the data
StateModel = c("H","L")
print(StateModel[as.integer(apply(gamma,1,function(x){x[1]>x[2]}))+1])

colors = c("red","blue")
plot(data,type="b", col=colors[as.integer(apply(gamma,1,function(x){x[1]>x[2]}))+1])

