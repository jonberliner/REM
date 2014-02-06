    library('stats')
    library('plyr')
    library('reshape2')
    library('ggplot2')
    library('LiblineaR')

    source('/mnt/cd/taylor/berliner/Sleep_12414/load_mnist.R')
    source('/mnt/cd/taylor/berliner/Sleep_12414/rbm.R')


######### NETWORK AND EXPERIMENT FREE PARAMS #########
    nHid <- 100 # number of hidden nodes
    lrate <- 0.05 # network learning rate
    sd_init_weights <- 0.05 # sd of normally distd random net weights

    # experiment free parameters
    maxI <- 10000 # number total number of examples presented for training
    nInEpoch <- 200 # an epoch here is number of samples before testing, not a run through every example
    # set training set sizes
    train_size <- 500
    test_size <- 2000


######### PREP TRAINING AND TESTING DATA #########
    if(!exists('fives') || !exists('twos') |
           !exists('l.fives') || !exists('l.twos')){
        load_mnist('~/Documents/Datasets')
        fives <-  (train$x[train$y==5,]/255 > .5)+0 # +0 is a trick to turn a boolean matrix numeric
        twos <-  (train$x[train$y==2,]/255 > .5)+0
        l.fives <- (train$y[train$y==5])
        l.twos <- (train$y[train$y==2])
    }

    iTwos.all <- seq(length(twos[,1]))
    iFives.all <- seq(length(fives[,1]))

    # get training and testing indices
    # get training indices
    iTwos.train <- sample(iTwos.all,train_size)
    iFives.train <- sample(iFives.all,train_size)

    # get test indices, ensuring no overlap with train indices
    iTwos.test <- sample(setdiff(iTwos.all,iTwos.train),test_size)
    iFives.test <- sample(setdiff(iFives.all,iFives.train),test_size)

    # get training and testing samples
    twos.train <- twos[iTwos.train,]
    fives.train <- fives[iFives.train,]
    twos.test <- twos[iTwos.test,]
    fives.test <- fives[iFives.test,]
    # get labels
    l.twos.train <- l.twos[iTwos.train]
    l.fives.train <- l.fives[iFives.train]
    l.twos.test <- l.twos[iTwos.test]
    l.fives.test <- l.fives[iFives.test]

    # combine training pool
    train_set <- rbind(twos.train,fives.train)
    l.train_set <- c(l.twos.train,l.fives.train)

    # shuffle training pool
    iShuffle <- sample(dim(train_set)[1])
    train_set <- train_set[iShuffle,]
    l.train_set <- l.train_set[iShuffle]
    nTrain <- dim(train_set)[1]



######### INITIATE NETWORK AND EXPERIMENT PARAMS #########
    nEpoch <- maxI / nInEpoch # will probe net performance nEpoch times

    # init network
    nVis <- dim(train_set)[2] # number of visible units
    W <- matrix(rnorm(nHid*nVis,mean=0,sd=sd_init_weights),
                nrow=nVis,ncol=nHid) # vis-hid weight matrix

    # init stats collection datafra
    dimslist <- c('epoch','m_mdist_22','m_mdist_25','m_mdist_52','m_mdist_55','m_err_recon2','m_err_recon5',
                       'sd_mdist_22','sd_mdist_25','sd_mdist_52','sd_mdist_55','sd_err_recon2','sd_err_recon5')
    progress <- data.frame(matrix(NaN,nrow=nEpoch,ncol=length(dimslist)))
    names(progress) <- dimslist


######### CONDUCT EXPERIMENT #########
    i <- 0 # i = n samples trained on so far
    trained = FALSE # end when trained = TRUE
    while(!trained){
        i <- i+1

        if(i%%nInEpoch == 0){
        ## TEST NETWORK PERFORMANCE
            image(W) # display weight matrix

            # get neural probs for samples in test set
            test.hidact.2 <- pH_aV(W,twos.test)
            test.hidact.5 <- pH_aV(W,fives.test)

            # get reconstructions of test data
            temp <- fire(test.hidact.2)
            test.recon.2 <- pV_aH(W,temp)
            temp <- fire(test.hidact.5)
            test.recon.5 <- pV_aH(W,temp)

            # calc avg activation of hidden layer
            mu.2 <- colMeans(test.hidact.2)
            mu.5 <- colMeans(test.hidact.5)
            mu.all <- colMeans(rbind(test.hidact.2,test.hidact.5))
            # calc avg co-activation of hidden layer
            cov.5 <- cov(test.hidact.5)
            cov.2 <- cov(test.hidact.2)
            cov.all <- cov(rbind(test.hidact.2,test.hidact.5))

            # get mahalonobis distance between and within cats
            m2dist.2 <- mahalanobis(test.hidact.2,center=mu.2,cov=cov.2)
            m2dist.5 <- mahalanobis(test.hidact.5,center=mu.2,cov=cov.2)
            m5dist.2 <- mahalanobis(test.hidact.2,center=mu.5,cov=cov.5)
            m5dist.5 <- mahalanobis(test.hidact.5,center=mu.5,cov=cov.5)

            mdist.2 <- mahalanobis(test.hidact.2,center=mu.all,cov=cov.all)
            mdist.5 <- mahalanobis(test.hidact.5,center=mu.all,cov=cov.all)
            # get euc dists between and within



            # get avg euclidean reconstruction error
            test_sse_2 <- rowSums( sqrt((twos.test-test.recon.2)^2) )
            test_sse_5 <- rowSums( sqrt((fives.test-test.recon.5)^2) )


            # STORE DATA
            progress$m_mdist_22[i/nInEpoch] <- mean(m2dist.2) # avg mahal dist from test '2' neural activation to mean test '2' neural act
            progress$m_mdist_25[i/nInEpoch] <- mean(m2dist.5)
            progress$m_mdist_52[i/nInEpoch] <- mean(m5dist.2)
            progress$m_mdist_55[i/nInEpoch] <- mean(m5dist.5)

            progress$sd_mdist_22[i/nInEpoch] <- sd(m2dist.2)
            progress$sd_mdist_25[i/nInEpoch] <- sd(m2dist.5)
            progress$sd_mdist_52[i/nInEpoch] <- sd(m5dist.2)
            progress$sd_mdist_55[i/nInEpoch] <- sd(m5dist.5)

            progress$m_err_recon2[i/nInEpoch] <- mean(test_sse_2)
            progress$m_err_recon5[i/nInEpoch] <- mean(test_sse_5)

            progress$sd_err_recon2[i/nInEpoch] <- sd(test_sse_2)
            progress$sd_err_recon5[i/nInEpoch] <- sd(test_sse_5)

            progress$epoch[i/nInEpoch] <- i # store how many examples seen so far
        }

        ## TRAIN NETWORK
        sample0 <- matrix(train_set[(i-1)%%nTrain+1,],nrow=1)
        W <- rbm_train(W,sample0,lrate=lrate)
        print(i)

        # stop if fully trained
        if(i>maxI){
            trained <- TRUE
        }
    }
    
    write.csv(progress,"jsb_test1.txt")