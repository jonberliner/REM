source('~/Dropbox/Sleep/datastream_fcns.R') # load fcns

cats <- sample(c(1,2,3,4,5,6,7,8,9,0),ncat) # randomly choose ncat categories
if(!exists('train') || !exists('test')){
    load_mnist('~/Documents/Datasets')
}

prep_single_waking <- function(catdata,l.catdata,n.train,n.test){
    # takes in samples and labels and returns a randomly partitioned train and test set
    n.tot <- nrow(catdata)
    i.tot <- seq(n.tot)
    i.train <- sample(i.tot,n.train) # get train inds
    i.test <- sample(setdiff(i.tot,i.train),n.test) # get non-overlapping test inds
    
    set.train <- catdata[i.train,]
    l.set.train <- l.catdata[i.train]
    set.test <- catdata[i.test,]
    l.set.test <- l.catdata[i.test]
    
    return(list('train'=set.train,'l.train'=l.set.train,
                'test'=set.test,'l.test'=l.set.test))
}




prep_all_waking <- function(mnist,n.train,n.test,cats,N){
    # mnist is the mnist$train list made with load_mnist()
    # n.train is the number of UNIQUE training samples, NOT total number of training presentations
    # n.test is the size of the testing set
    # cats is a vec with cats to use
    # N is total number of training samples to be made
    #   - will go through random cycles of the n.train unique samples until N presentations made
    
    n.cat <- length(cats)
    nPerEpi <- n.cat*n.train
    nTest <- n.cat*n.test
    
    # init arrays
    trainpool <- matrix(NA,nrow=nPerEpi,ncol=784)
    l.trainpool <- matrix(NA,nrow=nPerEpi,ncol=1)
    
    set.test <- matrix(NA,nrow=nTest,ncol=784)
    l.set.test <- matrix(NA,nrow=nTest,ncol=1)
    
    set.train <- matrix(NA,nrow=N,ncol=784)
    l.set.train <- matrix(NA,nrow=N,ncol=1)
    
    
    ci <- 1
    for(c in cats){
        catdata <-  (mnist$x[mnist$y==c,]/255 > .5)+0 # +0 is a trick to turn a boolean matrix numeric
        l.catdata <-  (mnist$x[mnist$y==c,]/255 > .5)+0 # also scaling to 0 and 1 making binary
        
        catsorted <- prep_single_waking(catdata,l.catdata,n.train,n.test)
        
        inds <- seq((ci-1)*n.train,ci*n.train-1)+1
        trainpool[inds,] <- catsorted$train
        l.trainpool[inds] <- catsorted$l.train
        
        inds <- seq((ci-1)*n.test,ci*n.test-1)+1
        set.test[inds,] <- catsorted$test
        l.set.test[inds] <- catsorted$l.test
        
        ci <- ci+1
    }
    
    # epicycles of samples to be presented until N samples
    n2go <- N
    ci <- 1
    shuffle <- seq(nPerEpi) # pool for shuffled inds
    while(n2go > 0){
        if(n2go >= n.train){
            shuffle <- sample(shuffle) # reshuffle
            inds <- seq((ci-1)*nPerEpi,ci*nPerEpi-1)+1
        } else {
            shuffle <- sample(shuffle,n2go)
            inds <- seq((ci-1)*nPerEpi,N-1)+1
        }
        
        set.train[inds,] <- trainpool[shuffle,]
        l.set.train[inds] <- l.trainpool[shuffle]
        
        n2go <- n2go - nPerEpi
        
        ci <- ci+1
    }
    
    return( list('set.train'=set.train,'l.set.train'=l.set.train,'set.test'=set.test,'l.set.test'=l.set.test))
}


make_waking_rem_stream <- function (waking,n.wake,n.sleep,rec_mag=1,p.noise=.2){
    # waking is a dataset made with prep_all_waking
    # n.wake is the number of waking trials in a cycle
    # n.sleep is the number of sleep samples in a trial
    #
    # returns a dataset with 5 REM conditions:  synthesis, repeat, repeat+noise,
    #                                           sparsity-matched noise, and random noise
    
    
    #waking <- prep_all_waking(train,n.train,n.test,cats,N)
    n.totwake <- nrow(waking)
    n.vis <- ncol(waking) # number visible units
    n.sleep_cycles <- floor(n.totwake / n.wake) # number of rem cycles to go through
    
    # note that we clip and dangling waking data
    n.tot <- (n.wake+n.sleep)*n.sleep_cycles
    
    # init full arrays
    rem.syn <- rem.rep <- rem.repn <- rem.noises <- rem.noiser <- matrix(NA,nrow=n.tot,ncol=n.vis)
    
    # inds.w is the waking indices in the full dataset.  inds.wo is the inds in the waking-only dataset
    inds.w <- inds.wo <- seq(n.waking)
    inds.r <- seq(n.sleep) + tail(inds.w,1)
    # let us begin!
    for(i in seq(n.sleep_cycle)){
        #         inds.w <- seq((i-1)*n.waking,i*n.waking-1)+1 # get this cycle's waking inds
        #         inds.r <- (seq((i-1)*n.sleep,i*n.sleep-1)+1)+tail(inds.w,1) # get sleep inds
        
        w0 <- waking[inds.wo,] # get wakings samples
        
        # add waking to full dataset
        rem.syn[inds.w,] <- rem.rep[inds.w,] <- rem.repn[inds.w,] <-
            rem.noises[inds.w,] <- rem.noiser[inds.w,] <- w0
        
        ##### ADD REM SAMPLES! #####
        # using a lot of functions found in rem_helpers.R
        #  ^ see these for detailed explanation of what's going on
        rec_weights <- get_recency_weights(n.waking,rec_mag) # get recency weights
        rempty <- matrix(NA,nrow=n.sleep,ncol=n.vis) # for re-initing
        ## synthesis rem
        r0 <- rempty # init empty rem data matrix
        
        fire_probs <- swirl_probs(w0,rec_weights) # get prob of each neuron firing
        
        # fill rem data with binary presentations made from fire_probs
        r0[] <- t(replicate(n.sleep,fire_probs)) # fire with prob fire_probs n.sleep times
        rem.syn[inds.r,] <- r0 # add to datastream
        
        ## replay rem
        r0 <- rempty # reinit rem matrix
        # get indices of samples to replay
        replay_inds <- sample(seq(n.wake),n.sleep,replace=TRUE,prob=rec_weights)
        r0 <- w0[replay_inds,] # add replayed samples to rem matrix
        rem.rep[inds.r,] <- r0 # add to datastream
        
        ## replay rem+noise
        # p.noise chance for a bit to be flipped
        flip <- matrix(rbinom(n.sleep*n.vis,1,p.noise),nrow=n.sleep) # get inds of bits to be flipped
        for(ii in seq(n.sleep)){
            r0[flip[ii,]] <- !r0[flip[ii,]] # flip
        }
        rem.repn[inds.r,] <- r0 # add to datastream
        
        ## sparsity-matched noise
        # matched to total sparsity, not neuron by neuron
        p <- mean(m0) # get total mean activation
        r0 <- matrix(rbinom(n.sleep*n.vis,1,p),nrow=n.sleep) # fire with prob p
        rem.noises <- r0
        
        ## random noise 
        # p changes for each sleep cycle
        # not constant p=0.5 because want to ensure not corrlated with actual data sparsity
        p <- runif(1) # get random p
        r0 <- matrix(rbinom(n.sleep*n.vis,1,p),nrow=n.sleep) # fire
        rem.noiser <- r0
        
        
        # increment indices
        inds.w <- inds.w + n.waking + n.sleep
        inds.r <- inds.r + n.waking + n.sleep
        inds.wo <- inds.wo + n.waking        
    }
    
    return( list('rem.syn'=rem.syn,'rem.rep'=rem.rep,'rem.repn'=rem.repn,
                 'rem.noises'=rem.noises,'rem.noiser'=rem.noiser))
}