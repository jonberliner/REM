# legend:
#   W: vis-hid layer weight matrix (nVis x nHid)
#   p*: prob of *
#   v/V: visible layer
#   h/H: hidden layer
#   a*: activation of *
#   _ : sometimes means "given" (e.g. a_b can mean "a given b" or "a | b")

logistic <- function(a){
    return (1 / (1+exp(-a)))
}

rlu <- function(a){
    # rectified linear unit
    a[a<0] <- 0
    return(a)
}

pH_aV <- function(W,aV){
# prob_hidden_activation given visible_activation
# get hidden layer activaiton probs given visible activation
    return(logistic(aV %*% W))
}


pV_aH <- function(W,aH){
# get hidden layer activaiton probs given visible activation
    return(logistic(aH %*% t(W)) )
}

get_associations <- function(v,h){
    return( t(v) %*% h )
}

fire <- function(p){
    # stochastically fire length(p) units, each unit i with prob(p[i])
    return( (p > matrix(runif(length(p)),nrow=nrow(p)))+0)
}

rbm_train <- function(W,input,lrate=0.05){
    stopifnot(sum(is.infinite(W))==0)
    # perform CD1 update on W, training on input vector
    v <- input
    pH <- pH_aV(W,v) # get prob of hidden firing given current vis activation (i.e. the input)
    # print(max(pH))
    h <- fire(pH) # stochasitically fire each hidden unit i with p(pH[i])
    pos_associations <- get_associations(v,pH) # store update gradient from posphase

    pV <- pV_aH(W,h)
    # print(max(pH))
    v <- fire(pV) # get reconstruction from neural activation
    pH <- pH_aV(W,v)
    # print(max(pH))
    h <- fire(pH) # stochasitically fire each hidden unit i with p(pH[i])
    neg_associations <- get_associations(v,pH)

    # perform CD update
    W <- W + lrate*(pos_associations-neg_associations)

    return(W)
}

flip_p <- function(a,p){
    # a is a binary array    
    # will flip bits with prob p
    flipmat <- matrix(rbinom(length(a),1,p),nrow=nrow(a))
    return(xor(a,flipmat))
}