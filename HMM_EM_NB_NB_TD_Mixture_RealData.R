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
B = 5.4
G = -0.47
X = 1:200
l = exp(B)*((X)^G)
beta = (l*0.2313)+1
a = rnbinom(200,mu=exp(B)*((X)^G), size=exp(B)*((X)^G)/((exp(B)*((X)^G)*0.2313)))
b = rnbinom(80,size=.5,mu=400)
data = c(a[1:60],b[1:20],a[61:100],b[21:40],a[101:160])+1
plot(data, type="l")

data = read.table("~/Documents/ResearchProjects/CC_Capture_Project/786O_pools/786O_pool5_1_2016-05-24/786O_pool5_1_2016-05-24_REdig_sorted_CC2_chr8:128888906-128889026_MYC_1_rs6470589_c.gff")
data = round(data[data[,1]=="chr8",6])
max = which(data==max(data))
X = abs(c(1:length(data))-max)
range = c((max-1000):(max+1000))
X = X[range]+1
data = data[range]
plot(data, type="l")
I = cbind(rep(1,(length(data))),log1p(X))
Y = log1p(data)
b = solve(t(I)%*%I)%*%t(I)%*%Y
B = log(max(data))
G = sum(((log(data)-B)/log1p(X)))/length(data)
#we can see that we have mixture of two distributions
#Let us initialise our Mixture EM data for every state
s21m = s22m = mean(data)+sd(data)
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
B1 = B2 = B
G1 = G2 = G

s12s = s22s = 2
s12m = s22m = mean(data)

Mean = (exp(B1)*(X^G1))
beta = ((Mean*0.46)+15.59)^2/Mean
r = Mean/(beta-1)
bs1w1 = dnbinom(data,mu=Mean,size=r)+1e-300
bs1w2 = dnbinom(round(data),size = s12s , mu = s12m)+1e-300

Mean = (exp(B2)*(X^G2))
beta = ((Mean*0.46)+15.59)^2/Mean
r = Mean/(beta-1)
bs2w1 = dnbinom(data,mu=Mean,size=r)+1e-300
bs2w2 = dnbinom(round(data),size = s22s , mu = s22m)+1e-300

bs1 = bs1w1*S11weights + bs1w2*S12weights
bs2 = bs2w1*S21weights + bs2w2*S22weights
  
for(x in 1:100)
{
#OK, now we are set to do the first HMM iteration
#Forward probability
#As mentioned above of the observation that we have, we will calculate the forward probability 
alpha = matrix(ncol=2,nrow=length(data))
#Forward for state 1 and first of the observation
alpha[1,1] = log(S1) + log(bs1[1]) #Probability of State1 * Probability of observing the first number
#Forward for state 2 and first of the observation
alpha[1,2] = log(S2) + log(bs2[1]) #Probability of State1 * Probability of observing the first number

#Second is just conditional on the previous one (could be either s1 or s2) and ENDING in state 1
a = alpha[1,1] + log(TS1[1]) +  log(bs1[2])
b = alpha[1,2] + log(TS2[2]) +  log(bs1[2])
alpha[2,1] =  log(sum(exp(c(a,b)-max(c(a,b)))))+max(c(a,b)) # ending in s1, P(S1) * Self transition * P(O2|S1)
             #ending in s1, P(S2) * T[S2->S1] * P(O2|S1)

a = alpha[1,2] + log(TS2[1]) +  log(bs2[2])
b = alpha[1,1] + log(TS1[2]) +  log(bs2[2])
alpha[2,2] =  log(sum(exp(c(a,b)-max(c(a,b)))))+max(c(a,b)) # ending in s2, P(S2) * Self transition * P(O2|S2)
              #ending in s2, P(S1) * T[S2->S1] * P(O2|S1)


for(i in 2:length(data))
{
  #Ok, let us see the pattern
  #Probability of Observation given state 1 * sum of ((probability of previous state1 * self probability)+
  #(Probability of previous state2 * transition from S2 to S1))
  a = alpha[i-1,1] + log(TS1[1]) + log(bs1[i])
  b = alpha[i-1,2] + log(TS2[2]) + log(bs1[i])
  alpha[i,1] = log(sum(exp(c(a,b)-max(c(a,b)))))+max(c(a,b))
  #Similarly
  a = alpha[i-1,2] + log(TS2[1]) + log(bs2[i])
  b = alpha[i-1,1] + log(TS1[2]) + log(bs2[i])
  alpha[i,2] = log(sum(exp(c(a,b)-max(c(a,b)))))+max(c(a,b))
}

#alpha = alpha+1e-300
#Probability of the whole observation across the two states
a = alpha[length(data),1]
b = alpha[length(data),2]
po = log(sum(exp(c(a,b)-max(c(a,b)))))+max(c(a,b))

print(a+b)
#OK now the backward probability

beta = matrix(ncol=2,nrow=length(data))
beta[length(data),1] = log(1)
beta[length(data),2] = log(1)
#Stand at the previous step and then look ahead
#Probility at length(data)-1 at state 1 = T[S1->S1] * P(O|S1) * B(length(data),S1) + T[S1->S2] * P(O|S2) * B(length(data),S2)
a = log(TS1[1]) + log(bs1[length(data)-1]) + beta[length(data),1]
b = log(TS1[2]) + log(bs2[length(data)-1]) + beta[length(data),2]
beta[length(data)-1,1] = log(sum(exp(c(a,b)-max(c(a,b)))))+max(c(a,b)) 
a = log(TS2[1]) + log(bs2[length(data)-1]) + beta[length(data),2]
b = log(TS2[2]) + log(bs1[length(data)-1]) + beta[length(data),1]
beta[length(data)-1,2] = log(sum(exp(c(a,b)-max(c(a,b)))))+max(c(a,b)) 

#pattern
for(j in (length(data)-1):1)
{
  a = log(TS1[1]) + log(bs1[j+1]) + beta[j+1,1]
  b = log(TS1[2]) + log(bs2[j+1]) + beta[j+1,2]
  beta[j,1] = log(sum(exp(c(a,b)-max(c(a,b)))))+max(c(a,b)) 
  a = log(TS2[1]) + log(bs2[j+1]) + beta[j+1,2]
  b = log(TS2[2]) + log(bs1[j+1]) + beta[j+1,1]
  beta[j,2] = log(sum(exp(c(a,b)-max(c(a,b)))))+max(c(a,b)) 
}

#Gamma, which is the probability of the state given the observation
gamma = matrix(ncol=2, nrow=length(data))
gamma[,1] = alpha[,1]+beta[,1]
gamma[,2] = alpha[,2]+beta[,2]
gamma = t(apply(gamma,1,function(b){exp(b-log(sum(exp(b-max(b))))-max(b))}))

#Eta, which is the transition probability from state i to j at data length t
eta1 = matrix(ncol=2,nrow=(length(data)-1)) #for every data observation
eta2 = matrix(ncol=2,nrow=(length(data)-1)) #for every data observation
#going from 1 to 1 across the first observation
eta1[1,1] = (log(gamma[1,1]) + log(TS1[1]) + log(bs1[2]) + beta[2,1])-beta[1,1]
#going from 1 to 2 across the first observation
eta1[1,2] = (log(gamma[1,1]) + log(TS1[2]) + log(bs2[2]) + beta[2,2])-beta[1,1]
#going from 2 to 2 across the first observation
eta2[1,1] = (log(gamma[1,2]) + log(TS2[1]) + log(bs2[2]) + beta[2,2])-beta[1,2]
#going from 1 to 2 across the first observation
eta2[1,2] = (log(gamma[1,2]) + log(TS2[2]) + log(bs1[2]) + beta[2,1])-beta[1,2]

#Pattern

for(t in 1:(length(data)-1))
{
  eta1[t,1] = (log(gamma[t,1]) + log(TS1[1]) + log(bs1[t+1]) + beta[(t+1),1])-beta[t,1]
  #going from 1 to 2 across the first observation
  eta1[t,2] = (log(gamma[t,1]) + log(TS1[2]) + log(bs2[t+1]) + beta[(t+1),2])-beta[t,1]
  #going from 2 to 2 across the first observation
  eta2[t,1] = (log(gamma[t,2]) + log(TS2[1]) + log(bs2[t+1]) + beta[(t+1),2])-beta[t,2]
  #going from 1 to 2 across the first observation
  eta2[t,2] = (log(gamma[t,2]) + log(TS2[2]) + log(bs1[t+1]) + beta[(t+1),1])-beta[t,2]
}

TS1[1] = sum(exp(eta1[,1]))/sum(gamma[1:(length(data)-1),1])
TS1[2] = sum(exp(eta1[,2]))/sum(gamma[1:(length(data)-1),1])
TS2[1] = sum(exp(eta2[,1]))/sum(gamma[1:(length(data)-1),2])
TS2[2] = sum(exp(eta2[,2]))/sum(gamma[1:(length(data)-1),2])

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
#s11m = sum(gammaM[,1]*data)/sum(gammaM[,1])
s12m = sum(gammaM[,2]*data)/sum(gammaM[,2])
#s21m = sum(gammaM[,3]*data)/sum(gammaM[,3])
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

#Estimate Decay parameter
Y = log1p(data)
I = cbind(rep(1,(length(data))),log1p(X))
W = gammaM[,1]
b<-solve(t(I)%*%diag(W)%*%I)%*%t(I)%*%diag(W)%*%Y
B1 = b[1]
G1 = b[2]
#G1 = sum(((log1p(data)-B1)/log1p(X))*gammaM[,1])/sum(gammaM[,1])

W = gammaM[,3]
b<-solve(t(I)%*%diag(W)%*%I)%*%t(I)%*%diag(W)%*%Y
B2 = b[1]
G2 = b[2]
#G2 = sum(((log1p(data)-B2)/log1p(X))*gammaM[,3])/sum(gammaM[,3])



# y2 = (data*gammaM[,3])
# #y2[y2<1] = exp(B2)*(X[y2<1]^G2)
# filter= gammaM[,3]>gammaM[,4]
# if(sum(filter)>2)model2 = lm(log((y2[filter]))~log(X[filter]))
#B1 = coef(model1)[1]
#G2 = sum(((log(data[-1])-B2)/log(X[-1]))*gammaM[-1,3])/sum(gammaM[-1,3])
Mean = (exp(B1)*(X^G1))
beta = ((Mean*0.5)+15.59)^2/Mean
r = Mean/(beta-1)
bs1w1 = dnbinom(data,mu=Mean,size=r)+1e-300
bs1w2 = dnbinom(round(data),size = s12s , mu = s12m)+1e-300

Mean = (exp(B2)*(X^G2))
beta = ((Mean*0.5)+15.59)^2/Mean
r = Mean/(beta-1)
bs2w1 = dnbinom(data,mu=Mean,size=r)+1e-300
bs2w2 = dnbinom(round(data),size = s22s , mu = s22m)+1e-300

bs1 = bs1w1*S11weights + bs1w2*S12weights
bs2 = bs2w1*S21weights + bs2w2*S22weights

print(c(B1,G1,s12m,s12s,B2,G2,s22m,s22s))
#print(c(s11s,s12s,s21s,s22s))
#print(c(S11weights,S12weights,S21weights,S22weights))
#print(c(S1,S2))
#print(c(TS1,TS2))
}

#Now we can get the state model for the data
StateModel = c("H","L")
print(StateModel[as.integer(apply(gammaM,1,function(x){(x[1]+x[3])>(x[2]+x[4])}))+1])

colors = c("red","blue")
plot(data,type="b", col=colors[as.integer(apply(gamma,1,function(x){(x[1])>(x[2])}))+1])

