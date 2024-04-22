# Changing the shape of the distribution of w
# import w

library(R.matlab)
w_vec <- readMat("...w_vec_6Perso.mat")
w_vec <- w_vec$w.vec

plot(w_vec, type="l")

g=6
b=6
w=31

# Info
age <- 1:w
perso <- 1:g
state <- 1:b


positions<- data.frame(
  age = rep(1:w, b*g),
  state = rep(rep(1:b, each=w), g),
  perso = rep(1:g, each = w*b),
  line = c(1:(b*g*w))
)


# Get the proportion in each age class
prop_age <-
  data.frame(
    age=age,
    prop=rep(0,w)
  )

for(i in 1:w){
  pos <- positions[which(positions$age==i),"line"]
  prop_age[i, "prop"] <- sum(w_vec[pos])
}


# Get the proportion in each breeding state
prop_state <-
  data.frame(
    state=state,
    prop=rep(0,b)
  )

for(i in 1:b){
  pos <- positions[which(positions$state==i),"line"]
  prop_state[i, "prop"] <- sum(w_vec[pos])
}

# Get the proportion in each personality class
prop_perso <-
  data.frame(
    perso=perso,
    prop=rep(0,g)
  )

for(i in 1:g){
  pos <- positions[which(positions$perso==i),"line"]
  prop_perso[i, "prop"] <- sum(w_vec[pos])
}


new_w_vec <- w_vec

for(j in 1:b){
  for(i in 1:w){
    pos <- positions[which(positions$age==i & positions$state==j), "line"]
    new_w_vec[pos] <- rep(1/g*sum(w_vec[pos]), g)
  }
}

#test
test <- positions[which(positions$age==20 & positions$state==2), "line"]
plot(new_w_vec[test])
plot(w_vec[test])
# it works

new_w_vec <- as.matrix(new_w_vec)

# export the new w vec
writeMat(new_w_vec = new_w_vec, "C:/Users/joani/Documents/PostDoc/Manuscripts/Personality_Demography/Matlab code/Hyperstate/Adult2008_standardized_correctfor20192020/new_w_vec_6perso_uniform.mat")
