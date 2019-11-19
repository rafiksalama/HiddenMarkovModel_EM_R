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
#Y = A*exp(-lx) => log(Y) = Log(A) - l*x

data = read.table("~/Documents/ResearchProjects/CC_Capture_Project/786O/CapC-786O/CapC-786O_REdig_sorted_CC1_chr16:30076905-30077296_ALDOA.gff")
data = round(data[data[,1]=="chr16",6]/((data[data[,1]=="chr16",5]-data[data[,1]=="chr16",4])/1000))+1
max = which(data==max(data))
X = abs(c(1:length(data))-max)
range = c((max-750):(max+750))
X = X[range]+1
data = data[range]
pdata = data
I = cbind(rep(1,(length(data))),log1p(X))
Y = log1p(pdata)
b = solve(t(I)%*%I)%*%t(I)%*%Y
norm_data = log(pdata/(exp(b[1])*X^b[2]))
plot(pdata, type="l")
lines((exp(b[1])*X^b[2]), col="red")
plot(norm_data, type="l")
#Bayesian Fit for the exponenial decay
#P(theta|data) = likelihood*priorjoint(theta)/sum(likelihood*priorjoint(theta)) for all possible theta
#Let us take theta1 as a uniform distribution over 0,1, theta2 as a gamma distribution centered around the max
#likelihood
loglik <- function(y,x,theta){dnorm(log(y/(theta[3]+theta[1]*exp(-theta[2]*x))+(max(data)*exp(-theta[4]*x))),mean=0,sd=2, log=T)}
prior <- c(function(theta){log(dgamma(theta,max(data)))},function(theta){log(1/100)})
#Posterior
posterior <- function(theta,prior,loglik){loglik+prior[[1]](theta[1])+prior[[2]](theta[2])}
#Grid approximation
m = matrix(ncol=100,nrow=100)
theta1_prior = seq(100,1750,length=100)
theta2_prior = seq(0,0.5,length=100)
theta3_prior = seq(0,100,length=10)
theta4_prior = seq(0,0.5,length=100)
mx = -1e09
curr_theta = c(theta1_prior[1],theta2_prior[1],theta3_prior[1])
for(i in 1:length(theta1_prior))
  for(j in 1:length(theta2_prior))
    for(x in 1:length(theta3_prior))
      for(n in 1:length(theta4_prior))
      {
      l = sum(loglik(data,X,c(theta1_prior[i],theta2_prior[j],theta3_prior[x],theta4_prior[n])))
      if(l>mx)
      {
        mx = l
        curr_theta = c(theta1_prior[i],theta2_prior[j],theta3_prior[x],theta4_prior[n])
      }
      #m[i,j] = posterior(loglik=sum(loglik(data,X,c(theta1_prior[i],theta2_prior[j]))),prior=prior,theta=c(theta1_prior[i],theta2_prior[j]))
      #if(j==100&&i%%10==0)print(i)
    }

mr = apply(m,1,max)
mc = apply(m,2,max)
print(theta1_prior[which(mr==max(mr))])
print(theta2_prior[which(mc==max(mr))])
#What you are aiming at is simply constructing a matrix of all thetas together
#
#data = norm_data
#we can see that we have mixture of two distributions
#Let us initialise our Mixture EM data for every state
NumberOfStates = 2
NumberOfMixture = 2
NumberOfParameters = 2
likelihoodfunction <- list(
  #Gaussian Fit
  function(data,theta){
    #theta is the fit parameters
    #norm_data = log(data/(exp(theta[2])*(X^theta[1])))
    norm_data = log(data/(11+133*exp(-.02*X)+(1200*exp(-0.5*X))))
    dnorm(norm_data,0,1)+1e-300
  },
  #Negative Binomial
  function(data,theta){
    dnbinom(round(data),mu = theta[1], size = theta[2])+1e-300
  }
)

MixtureParameters = rep(list(matrix(ncol=NumberOfStates, nrow=NumberOfParameters)),NumberOfMixture)
#for(i in 1:NumberOfMixture)
#{
#Gaussian fit
MixtureParameters[[1]][1,] = b[1]+runif(NumberOfStates)
MixtureParameters[[1]][2,] = b[2]+runif(NumberOfStates)

#Negative Binomial
MixtureParameters[[2]][1,] = mean(data)+runif(NumberOfStates)
MixtureParameters[[2]][2,] = 2+runif(NumberOfStates)
#}
#Let us now initialise the state probability i.e. P(Si)
S = rep(0.5, NumberOfStates)
#Transition probability self, then the other
TS = matrix(ncol=NumberOfStates,nrow=NumberOfStates, 0.5)

Sweights = matrix(nrow=NumberOfStates, ncol=NumberOfMixture, 0.5)

#Emission probability is calculated as in the EM approach, this makes sense that we have a probability for every observation, since they are not the same observation
#Except that we don't have a mixture in the states, essentially we want every state to map to one component of the mixture

#Likelihood
Mixturelikelihood = matrix(nrow=length(data), ncol=NumberOfMixture,0)
StateLikelihood = rep(list(Mixturelikelihood),NumberOfStates)
bs = matrix(ncol=NumberOfStates,nrow=length(data))

for(i in 1:NumberOfStates)
  for(j in 1:NumberOfMixture)
  {
    StateLikelihood[[i]][,j] = likelihoodfunction[[j]](data,MixtureParameters[[j]][,i])+1e-300
  }

for(i in 1:NumberOfStates){
  bs[,i] = StateLikelihood[[i]] %*% Sweights[i,]
}
diff=1
prev = -10000
continue = T
while(continue)
{
#OK, now we are set to do the first HMM iteration
#Forward probability
#As mentioned above of the observation that we have, we will calculate the forward probability 
alpha = matrix(ncol=NumberOfStates,nrow=length(data))
#Forward for state 1 and first of the observation
for(i in 1:NumberOfStates) alpha[1,i] = log(S[i]) + log(bs[1,i])

for(i in 2:length(data))
{
  for(j in 1:NumberOfStates)
  {
    v = rep(0, NumberOfStates)
    for(l in 1:NumberOfStates)
    {
      #Ok, let us see the pattern
      #Probability of Observation given state 1 * sum of ((probability of previous state1 * self probability)+
      #(Probability of previous state2 * transition from S2 to S1))
      v[l] = alpha[i-1,l] + log(TS[l,j]) + log(bs[i,j])
    }
    alpha[i,j] = log(sum(exp(v-max(v))))+max(v)
  }
}

#OK now the backward probability

beta = matrix(ncol=NumberOfStates,nrow=length(data))
for(i in 1:NumberOfStates) beta[length(data),i] = log(1)

#Stand at the previous step and then look ahead
#Probility at length(data)-1 at state 1 = T[S1->S1] * P(O|S1) * B(length(data),S1) + T[S1->S2] * P(O|S2) * B(length(data),S2)
#pattern
for(i in (length(data)-1):1)
{
  for(j in 1:NumberOfStates)
  {
    for(l in 1:NumberOfStates)
    {
      v[l] = log(TS[j,l]) + log(bs[i+1,j]) + beta[i+1,l]
    }
    beta[i,j] = log(sum(exp(v-max(v))))+max(v)
  }
}

#Gamma, which is the probability of the state given the observation
gamma = matrix(ncol=NumberOfStates, nrow=length(data))
for(i in 1:NumberOfStates) gamma[,i] = alpha[,i]+beta[,i]
gamma = t(apply(gamma,1,function(b){exp(b-log(sum(exp(b-max(b))))-max(b))}))

#Eta, which is the transition probability from state i to j at data length t
eta = rep(list(matrix(ncol=NumberOfStates,nrow=(length(data)-1))), NumberOfStates)

#Pattern

for(t in 1:(length(data)-1))
  for(j in 1:NumberOfStates)
    for(l in 1:NumberOfStates)
      eta[[j]][t,l] = (log(gamma[t,j]) + log(TS[j,l]) + log(bs[t+1,l]) + beta[(t+1),l])-beta[t,j]

for(j in 1:NumberOfStates)
  for(l in 1:NumberOfStates)
    TS[j,l] = sum(exp(eta[[j]][,l]))/sum(gamma[1:(length(data)-1),j]) 

#if(gamma[1,1]==S1&&gamma[1,2]==S2) break;

for(i in 1:NumberOfStates) S[i] = gamma[1,i]

#Mixture gamma
gammaM = rep(list(matrix(ncol=NumberOfMixture, nrow=length(data))), NumberOfStates)
for(i in 1:NumberOfStates)
  for(j in 1:NumberOfMixture)
    gammaM[[i]][,j] = gamma[,i]*Sweights[i,j]*StateLikelihood[[i]][,j]/bs[,i]

  
for(i in 1:NumberOfStates)
  for(j in 1:NumberOfMixture)
    Sweights[i,j] = sum(gammaM[[i]][,j])/sum(gamma[,i])    
  
#Normal Distribution Fit
for(i in 1:NumberOfStates)
    MixtureParameters[[1]][,i] = solve(t(I)%*%diag(gammaM[[i]][,1])%*%I)%*%t(I)%*%diag(gammaM[[i]][,1])%*%Y

#Negative Binomial
for(i in 1:NumberOfStates)
{
  MixtureParameters[[2]][1,i] = sum(gammaM[[i]][,2]*data)/sum(gammaM[[i]][,2])
  MixtureParameters[[2]][2,i] = sqrt(sum(gammaM[[i]][,2]*((data-MixtureParameters[[2]][1,i])*t(data-MixtureParameters[[2]][1,i])))/sum(gammaM[[i]][,2]))
  vv = sum(gammaM[[i]][,2]*((data-MixtureParameters[[2]][1,i])*t(data-MixtureParameters[[2]][1,i])))/sum(gammaM[[i]][,2])
  be = vv/MixtureParameters[[2]][1,i]
  MixtureParameters[[2]][2,i] = MixtureParameters[[2]][1,i]/(be-1)
}

#Update the distribution parameters, mean and standard deviation similarly

for(i in 1:NumberOfStates)
  for(j in 1:NumberOfMixture)
    StateLikelihood[[i]][,j] = likelihoodfunction[[j]](data,MixtureParameters[[j]][,i])+1e-300


for(i in 1:NumberOfStates)
  bs[,i] = StateLikelihood[[i]] %*% Sweights[i,]


pd = t(gamma %*% S) %*% t(data %*% diag(rep(1/length(data), length(data))))
diff = pd - prev
prev = pd

if(diff<1E-06) stop = stop+1
else stop=0

if(stop>40) continue = F

print(paste(diff,stop))

}

#Now we can get the state model for the data
#StateModel = c("H","L")
#print(StateModel[as.integer(apply(gammaM,1,function(x){(x[1]+x[3])>(x[2]+x[4])}))+1])

colors = c("red","blue", "green", "red")
plot(data,type="b", col=colors[apply(gamma,1,function(x){which(x==max(x))})])

