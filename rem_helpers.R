fire <- function(probs){
    # function fires each unit i in vec probs with prob probs[i]
    return( (probs > runif(length(probs)))+0 )
}

get_recency_weights <- function(n,rec_mag){
    # rec_mag adds more to weight to things seen more recently when making firing combinations
    #   - 0=uniform weight
    #   - as goes towards inf, probs will approach last seen item
    w <- exp( rec_mag*seq(n))
    return( w / max(w) )
}

swirl_probs <- function(a,rec_weights){
    # a is data array to be warped into training examples
    # rec_weights are recency weights (make with get_recency_weights)
    # returns prob vector to sample from to make REM samples
    
    # weight a by recency
    weighted_a <- apply(a,2,function(c){rec_weights*c})
    
    probs <- colSums(weighted_a) / sum(rec_weights)
    stopifnot(sum(probs>1)==0) # assert valid prob for every unit (i.e. between 0 and 1)
    return(probs)
}

add_noise <- function(a,noise_sd){
    # adds uniform noise to array a with sd noise_sd
    return( a+matrix(rnorm(length(a),sd=noise_sd),nrow=nrow(a)))
}
