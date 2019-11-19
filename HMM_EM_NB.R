#The main aim is to obtimise the allocation of the states and the transition between them
#So we use an EM approach. Essentially estimate the transition probability and emission probability
#Should the current parameters be the most optimum and then use those probabilities to find
#The most likely parameters of the model

#OK we are going to initialise so that we have what we need for the first E-Step
#We are then going to calculate the E-Step and use it recalculate the parameters and so on
#Until convergence.

#The E-Step
#Forward Step, is the probability of obtaining a certain sequence from model H
#Calculate the probability of seeing the sequence y1..t and being in state i at the final observation yt 
#is only dependent on the previous probability of being in state i-1 at t-1, 
#The probability of seeing the first observation y1 of the sequecne and be in state i 
#will then be P(y1 & being in state 1) = P(S1)*P(y1)
#But we have multiple possibilities for the first state, either state 1, 2, .. N possible states that we could start from
#The second observation then I could have come from any state and then ended in state j
#So I need to sum for all the possible states I could have come from ending in state 2
#Hence, the probability of being in state 2 is P(S2)*P(y2|S2). The probability of State 2 is then
#The sum over all the previous state probabilities * the probability of transition to this state
#In general then, we could write this as follows 
#P(y,j,s) = Current state probability * sum(Previous State probabilities*transition)
#This will have to be done recursively for every sequence of observation and you will then end up
#with a probability of ending in every state for a sequence of observation. If we add them up, this will give us the probability
#of the observation given our model
#This gave us the probability of certain state given the starting model parameters


#Backward Step, this is the probability of looking at the model backwards rather forward.
#We start from the end state and then go backward asking what is the probability of observation given the state
#The probability of being in previous state j at time t is a function of being in state i at t+1
#P(End State T) = 1, just start with one
#P(Observing j of the observation in the Current State) = P(Current State)*P(J Observation|Current state i)*P(Transition from Current State to Previous state)
#Since you could have gone to multiple states, then you need to sum over all the states that you could have gone to
#P(O,length-1,Ending State) = sum(P(Observation at length-1|current state)*P(current->next)*P(Next)) for all the next states
#You can do this iteratively.
#At the End you will get probability for every starting state
#This will then give us the probability of Observation up to a certain length given that I am in a certain time and certain state

#Ok, now the first aim is to find the probability of being in a certain state given the observation and model parameters
#P(S|O,M) = P(O,S,M)/P(O,M), since P(O,S,M) = P(S|O,M)P(O|M)P(M)
#P(O,S,M) = P(O,S|M)P(M)
#So P(S|O,M) = P(O,S|M)/P(O|M)
#Ok, the P(O,S|M) = #forward P(O|S,M) #backward P(S|M), the forward probability gives is the P(S|M), the backward gives the probability of observation given states
#Now, the denominator is P(O|M), which is simply the marginal over all the states of P(O,S|M)
#So the probability of ending up in state i is then forward*backward of state i/sum over all states

#The second aim is to find the probability of a sequence of states as this will help in the final transition matrix estimation
#The probability of being in state i and then going to state j globally
#P(Si,Si+1|O,M) similar to the above = P(Si,Si+1,O|M)/P(O|M)
#Numerator, forward probability P(Si|M)*P(Si+1|M)*transition probability from Si to Si+1
#Deonominator is simply summing over all possible si and si+1, which is actually Backward probability

#Ok from the first aim, for every state, summing over all obervations Sum(P(S|O,M)) for all lengths of observations, 
#will give us the probability of being in state S|M, which is what we started with as our assumption
#Similarly summing over all the length of observations, will give us the transition matrix, again one of our initial assumptions
#Remaining, is the probability of our observation given a certain state

#To covert them to probabilities, You essentially, want to sum P(S|O,M) for the first observation across all the sequences, for the first state,
#Then the second state and so on. The final probability will be just the normalised quantity

#The transition from Si to Si+1 will be the transition function above for all the observations irrespective of the position of the observation
#The ratio between the transitions from Si to Si+1 relative to all transitions away from Si

#The probability of the output for a gaussian mixture is as follows will be the sum over all components with the associated parameters for this component
#Updating this probability will be similar to the normal EM approach.
#Ok, let us remember, we have calculated the weights of every component, as the normalised likelihood for this observation
#In this case, we also have the states, so we split things further. 
#First, the probability that the lth component in state i have generated a particular observation t is defined as:
#Probability of state i * probability of observation given mixture l * weight of mixture l / sum over all components at state i
#We have the initial weights to calculate this.
#We can now get the updates as weight as before, but instead of summing for one observation, we get it for two states
#So on


#OK let us start with some minor observations of a gaussian distribution with different means
library(lubridate)
a = rnbinom(350,size=1,mu=9)
b = rnbinom(50,size=1,mu=10)
data = c(a[1:50],b[1:10],a[50:100])


traffic = fread("~/Downloads/TEST_TRAFFIC.csv")
traffic$date = anytime(traffic$date)
traffic = traffic %>% mutate(day=day(date), year=year(date),month=month(date), hour=hour(date), minute=minute(date))
daily = traffic %>% group_by(year,month,day,hour) %>% summarise(value=sum(value))
data = daily$value
plot(data, type="l")
#we can see that we have mixture of two distributions
#Let us initialise our Mixture EM data for every state

s1m = m[1]
s2m = m[2]

s1s = s[1]
s2s = s[2]

#Let us now initialise the state probability i.e. P(Si)
S1 = 0.5
S2 = 0.5
#Transition probability self, then the other
TS1 = c(0.1, 0.9)
TS2 = c(0.1, 0.9)

#Emission probability is calculated as in the EM approach, this makes sense that we have a probability for every observation, since they are not the same observation
#Except that we don't have a mixture in the states, essentially we want every state to map to one component of the mixture
bs1 = dnbinom(data,size = s1s , mu = s1m)
bs2 = dnbinom(data,size = s2s , mu = s2m)

for(x in 1:10)
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
  alpha[i,1] = exp(log(bs1[i]) + log(alpha[i-1,1] * TS1[1] + alpha[i-1,2] * TS2[2]))+1e-300
  #Similarly
  alpha[i,2] = exp(log(bs2[i]) + log(alpha[i-1,2] * TS2[1] + alpha[i-1,1] * TS1[2]))+1e-300
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
  beta[j,1] = exp(log(TS1[1]) + log(bs1[j+1]) + log(beta[j+1,1])) + exp(log(TS1[2]) + log(bs2[j+1]) + log(beta[j+1,2]))+1e-300
  beta[j,2] = exp(log(TS2[1]) + log(bs2[j+1]) + log(beta[j+1,2])) + exp(log(TS2[2]) + log(bs1[j+1]) + log(beta[j+1,1]))+1e-300
}

#Gamma, which is the probability of the state given the observation
gamma = matrix(ncol=2, nrow=length(data))
gamma[,1] = exp(log(alpha[,1])+log(beta[,1]))+1e-300
gamma[,2] = exp(log(alpha[,2])+log(beta[,2]))+1e-300
gamma = t(apply(gamma,1,function(x){x/sum(x)}))

#Eta, which is the transition probability from state i to j at data length t
eta1 = matrix(ncol=2,nrow=(length(data)-1)) #for every data observation
eta2 = matrix(ncol=2,nrow=(length(data)-1)) #for every data observation
#going from 1 to 1 across the first observation
eta1[1,1] = exp(log(gamma[1,1]) + log(TS1[1]) + log(bs1[2]) + log(beta[2,1]) - log(beta[1,1]))+1e-300
#going from 1 to 2 across the first observation
eta1[1,2] = exp(log(gamma[1,1]) + log(TS1[2]) + log(bs2[2]) + log(beta[2,2]) - log(beta[1,1]))+1e-300
#going from 2 to 2 across the first observation
eta2[1,1] = exp(log(gamma[1,2]) + log(TS2[1]) + log(bs2[2]) + log(beta[2,2]) - log(beta[1,2]))+1e-300
#going from 1 to 2 across the first observation
eta2[1,2] = exp(log(gamma[1,2]) + log(TS2[2]) + log(bs1[2]) + log(beta[2,1]) - log(beta[1,2]))+1e-300

#Pattern

for(t in 1:(length(data)-1))
{
  #going from 1 to 1 across the first observation
  eta1[t,1] = exp(log(gamma[t,1]) + log(TS1[1]) + log(bs1[(t+1)]) + log(beta[(t+1),1])-log(beta[t,1]))+1e-300
  #going from 1 to 2 across the first observation
  eta1[t,2] = exp(log(gamma[t,1]) + log(TS1[2]) + log(bs2[(t+1)]) + log(beta[(t+1),2])-log(beta[t,1]))+1e-300
  #going from 2 to 2 across the first observation
  eta2[t,1] = exp(log(gamma[t,2]) + log(TS2[1]) + log(bs2[(t+1)]) + log(beta[(t+1),2])-log(beta[t,2]))+1e-300
  #going from 1 to 2 across the first observation
  eta2[t,2] = exp(log(gamma[t,2]) + log(TS2[2]) + log(bs1[(t+1)]) + log(beta[(t+1),1])-log(beta[t,2]))+1e-300
}

TS1[1] = sum(eta1[,1])/sum(gamma[1:(length(data)-1),1])
TS1[2] = sum(eta1[,2])/sum(gamma[1:(length(data)-1),1])
TS2[1] = sum(eta2[,1])/sum(gamma[1:(length(data)-1),2])
TS2[2] = sum(eta2[,2])/sum(gamma[1:(length(data)-1),2])

#if(gamma[1,1]==S1&&gamma[1,2]==S2) break;

S1 = gamma[1,1]
S2 = gamma[1,2]

#Update the distribution parameters, mean and standard deviation similarly
s1m = sum(gamma[,1]*data)/sum(gamma[,1])
s2m = sum(gamma[,2]*data)/sum(gamma[,2])

v1s = sum(gamma[,1]*((data-s1m)*t(data-s1m)))/sum(gamma[,1])
v2s = sum(gamma[,2]*((data-s2m)*t(data-s2m)))/sum(gamma[,2])

be = c(v1s/s1m,v2s/s2m)
for(x in 1:length(be)) if(be[x] < 1) be[x]=2 

s1s = s1m/(be[1]-1)
s2s = s2m/(be[2]-1)

bs1 = dnbinom(data,size = s1s , mu = s1m)
bs2 = dnbinom(data,size = s2s , mu = s2m)
print(c(s1m,s2m))
print(c(s1s,s2s))
print(c(S1,S2))
print(c(TS1,TS2))
}

#Now we can get the state model for the data
StateModel = c("H","L")
print(StateModel[as.integer(apply(gamma,1,function(x){x[1]>x[2]}))+1])

colors = c("red","blue")
plot(data,type="b", col=colors[as.integer(apply(gamma,1,function(x){x[1]>x[2]}))+1])

