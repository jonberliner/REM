library('plyr')
library('ggplot2')
library('gtools')
library('LiblineaR')

source('datastream_fcns.R')
source('nnet_metrics_fcns.R')
source('rbm.R')
source('~/Dropbox/Sleep/mnist_helpers.R') # load fcns for loading and viewing mnist raw data from yann lecunn's site

perm <- function(n,k){choose(n,k) * factorial(k)}

# load mnist dataset
if(!exists('mnist')){
    load_mnist('~/Documents/Datasets')
}

######### EXPERIMENT FREE PARAMS #########
    n.cat  <- 3 # number of categories
    n.train <- 200 # number of unique samples per catergory in the training set
    n.test <- 200 # number of unique samples per catergory in the testing set
    N.waking <- 1000 # total size of the WAKING training set
    
    n.wake <- 100 # number of waking trials per cycle
    n.sleep <- 10 # number of REM trials per cycle
    rec_mag <- 1 # recency weighting curve
      #   0 means all waking samples have same weight
      #   approaching inf means only most recent has weight
      #   effects rem.syn, rem.rep, rem.repn 
     p.noise <- 0.2 # p.noise bits will be flipped per replay
      #   effects rem.repn

######### NETWORK FREE PARAMS #########
n.hid <- 100 # number of hidden nodes
lrate <- 0.05 # network learning rate
sd_init_weights <- 0.05 # sd of normally distd random net weights


######### CREATE DATA STREAMS #########
sets <- make_datastreams(mnist=mnist,n.cat=n.cat,n.train=n.train,n.test=n.test,N=N.waking,
                         n.wake=n.wake,n.sleep=n.sleep,rec_mag=rec_mag,p.noise=p.noise)

######### INITIATE NETWORK AND EXPERIMENT PARAMS #########
# init datastreams
test_set <- sets$test
l.test_set <- sets$l.test
train_set <- sets$train$rem.syn
l.train_set <- sets$train$l.train
N.total <- nrow(train_set)

# exp params
n.cycle <- n.wake + n.sleep # a cycle is a wake-sleep cycle
N.cycle <- N.total / n.cycle

# init network
n.vis <- ncol(train_set) # number of visible units
W <- matrix(rnorm(n.hid*n.vis,mean=0,sd=sd_init_weights),
            nrow=n.vis,ncol=n.hid) # vis-hid weight matrix

# init stats collection lists
nrow.mahal <- perm(n.cat,2)
nrow.rec_err <- n.cat
stats_out.rec_err <- as.data.frame(matrix(NA,nrow=nrow.rec_err*N.cycle,ncol=5))
names(stats_out.rec_err) <- c('cat','mu','std','cycle','time_of_day')
stats_out.mahal <- as.data.frame(matrix(NA,nrow=nrow.mahal*N.cycle,ncol=6))
names(stats_out.mahal) <- c('to','from','mu','std','cycle','time_of_day')

######### CONDUCT EXPERIMENT #########
i <- 1 # i = ith training example
i.cycle <- 0
inds.mahal <- seq(nrow.mahal)
inds.rec_err <- seq(nrow.rec_err)

for(i in seq(N.total)){
    
    ## TRAIN NETWORK
    W <- rbm_train(W,matrix(train_set[i,],nrow=1),lrate=lrate)
        
    ## TEST NETWORK PERFORMANCE
    if( i%%n.cycle == 0 | # (waking up)
           i%%n.cycle == n.wake){ # (going to sleep)
        image(W) # display weight matrix
        
        if(i%%n.cycle == 0){
            i.cycle <- i.cycle + 1
            time_of_day <- 'morning'
            print(cat('cycle',as.character(i.cycle),'of',as.character(N.cycle)))
        } else{
            time_of_day <- 'evening'
        }
        
        cats <- unique(c(l.test_set))
        n.cat <- length(cats)
        catnames <- as.character(cats)
        
        # get probs of hidden activation on testing set
        hidprobs.all <- pH_aV(W,test_set) # get hidden layer activation for all samples
        
        # get mahalanobis metrics
        mahal <- mahal_metrics(W,hidprobs.all,l.test_set)
        mahal$cycle <- i.cycle
        mahal$'time_of_day' <- time_of_day
        
        # get reconstructions of all samples
        temp <- fire(hidprobs.all)
        recon.all <- pV_aH(W,temp)
        
        # get euclidean reconstruction error
        rec_err.indiv <- llply(cats,function(c){
            rowSums( sqrt((test_set[l.test_set==c,] - recon.all[l.test_set==c,])^2) )
        })
        names(rec_err.indiv) <- catnames
        
        # get avg recon err
        rec_err <- ldply(catnames,function(c){
            data.frame('cat'=c,
                       'mu'=mean(rec_err.indiv[[c]]),
                       'std'=sd(rec_err.indiv[[c]])
            )
        })
        rec_err$cycle <- i.cycle
        rec_err$time_of_day <- time_of_day
        
        # STORE DATA
        stats_out.rec_err[inds.rec_err,] <- rec_err
        stats_out.mahal[inds.mahal,] <- mahal
        
        inds.rec_err <- inds.rec_err + nrow.rec_err
        inds.mahal <- inds.mahal + nrow.mahal
        
    }
}

# 
# # visualize
g <- ggplot(stats_out.mahal)
g + geom_point(aes(x=cycle,y=mu,colour=to,shape=time_of_day)) + scale_y_log10()
