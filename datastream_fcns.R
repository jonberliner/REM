
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




make_wakingstream_and_testset <- function(mnist,n.train,n.test,cats,N){
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
        l.catdata <-  mnist$y[mnist$y==c] # also scaling to 0 and 1 making binary
        
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
        
        if(n2go < nPerEpi){
            inds <- inds[1:n2go]
            shuffle <- shuffle[1:n2go]
        }
            
        set.train[inds,] <- trainpool[shuffle,]
        l.set.train[inds] <- l.trainpool[shuffle]
                    
        n2go <- n2go - nPerEpi
        
        ci <- ci+1
    }
    
    return( list('set.train'=set.train,'l.set.train'=l.set.train,'set.test'=set.test,'l.set.test'=l.set.test))
}


make_waking_rem_streams <- function (waking,l.waking,n.wake,n.sleep,rec_mag,p.noise){
    # waking is a dataset made with prep_all_waking
    # n.wake is the number of waking trials in a cycle
    # n.sleep is the number of sleep samples in a trial
    #
    # returns a dataset with 5 REM conditions:  synthesis, repeat, repeat+noise,
    #                                           sparsity-matched noise, and random noise
    
    
    #waking <- prep_all_waking(train,n.train,n.test,cats,N)
    n.totwake <- nrow(waking)
    n.vis <- ncol(waking) # number visible units
    n.sleep_cycle <- floor(n.totwake / n.wake) # number of rem cycles to go through
    
    # note that we clip and dangling waking data
    n.tot <- (n.wake+n.sleep)*n.sleep_cycle
    
    # init full arrays
    rem.syn <- rem.rep <- rem.repn <- rem.noises <- rem.noiser <- matrix(NA,nrow=n.tot,ncol=n.vis)
    l.train <- matrix(NA,nrow=n.tot,ncol=1)
    
    # inds.w is the waking indices in the full dataset.  inds.wo is the inds in the waking-only dataset
    inds.w <- inds.wo <- seq(n.wake)
    inds.r <- seq(n.sleep) + tail(inds.w,1)
    # let us begin!
    for(i in seq(n.sleep_cycle)){
        
        w0 <- waking[inds.wo,] # get wakings samples
        l.w0 <- l.waking[inds.wo] # get labels
        
        # add waking to full dataset
        rem.syn[inds.w,] <- rem.rep[inds.w,] <- rem.repn[inds.w,] <-
            rem.noises[inds.w,] <- rem.noiser[inds.w,] <- w0
        
        # add waking labels
        l.train[inds.w] <- l.w0
        l.train[inds.r] <- NA
        
        
        ##### ADD REM SAMPLES! #####
        # using a lot of functions found in rem_helpers.R
        #  ^ see these for detailed explanation of what's going on
        rec_weights <- get_recency_weights(n.wake,rec_mag) # get recency weights
        rempty <- matrix(NA,nrow=n.sleep,ncol=n.vis) # for re-initing
        
        ## synthesis rem
        r0 <- rempty # init empty rem data matrix
        
        fire_probs <- swirl_probs(w0,rec_weights) # get prob of each neuron firing
        
        # fill rem data with binary presentations made from fire_probs
        r0[] <- t(replicate(n.sleep,fire(fire_probs))) # fire with prob fire_probs n.sleep times
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
        r0 <- xor(r0,flip)
        rem.repn[inds.r,] <- r0 # add to datastream
        
        ## sparsity-matched noise
        # matched to total sparsity, not neuron by neuron
        r0 <- rempty
        p <- mean(fire_probs) # get recency-weighted total mean activation
        r0 <- matrix(rbinom(n.sleep*n.vis,1,p),nrow=n.sleep) # fire with prob p
        rem.noises[inds.r,] <- r0
        
        ## random noise 
        # p changes for each sleep cycle
        # not constant p=0.5 because want to ensure not corrlated with actual data sparsity
        r0 <- rempty
        p <- runif(1) # get random p
        r0 <- matrix(rbinom(n.sleep*n.vis,1,p),nrow=n.sleep) # fire
        rem.noiser[inds.r,] <- r0
        
        # increment indices
        inds.w <- inds.w + n.wake + n.sleep
        inds.r <- inds.r + n.wake + n.sleep
        inds.wo <- inds.wo + n.wake   
    }
    
    return( list('rem.syn'=rem.syn,'rem.rep'=rem.rep,'rem.repn'=rem.repn,
                 'rem.noises'=rem.noises,'rem.noiser'=rem.noiser,'l.waking'=l.train))
}


make_datastreams <- function(mnist,n.cat,n.train,n.test,N,n.wake,n.sleep,rec_mag,p.noise){
    #     mnist # mnist dataset made with load_mnist
    #     n.cat  <- 2 # number of categories
    #     n.train # number of unique samples per catergory in the training set
    #     n.test # number of unique samples per catergory in the testing set
    #     N # total size of the waking training set
    #     
    #     n.wake # number of waking trials per cycle
    #     n.sleep # number of REM trials per cycle
    #     rec_mag # recency weighting curve
    #       #   0 means all waking samples have same weight
    #       #   approaching inf means only most recent has weight
    #       #   effects rem.syn, rem.rep, rem.repn 
    #      p.noise <- 0.2 # p.noise bits will be flipped per replay
    #       #   effects rem.repn
    
    source('~/Dropbox/Sleep/rem_helpers.R') # load helpers for making data
    
    cats <- sample(c(1,2,3,4,5,6,7,8,9,0),n.cat) # randomly choose ncat categories
    waking_and_test <- make_wakingstream_and_testset(mnist,n.train,n.test,cats,N) # make waking data stream AND test set
    
    waking_rem_streams <- make_waking_rem_streams(waking_and_test$set.train,waking_and_test$l.set.train,n.wake,n.sleep,rec_mag,p.noise) # make full waking+REM streams!
    return(list('test'=waking_and_test$set.test,'l.test'=waking_and_test$l.set.test,
                'train'=waking_rem_streams,'l.train'=waking_rem_streams$l.waking))
}