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
