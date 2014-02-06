setwd('~/Dropbox/Sleep/')

library('plyr')
library('ggplot2')
library('gtools')
library('LiblineaR')
library('corpcor')

source('datastream_fcns.R')
source('nnet_metrics_combined.R')
source('rbm.R')
source('~/Dropbox/Sleep/mnist_helpers.R') # load fcns for loading and viewing mnist raw data from yann lecunn's site


output_name <- '~/Dropbox/Sleep/Outputs/test.csv'
perm <- function(n,k){choose(n,k) * factorial(k)}

# load mnist dataset
if(!exists('mnist')){
    load_mnist('~/Documents/Datasets')
}

######### EXPERIMENT FREE PARAMS #########
n.cat  <- 2 # number of categories
n.train <- 100 # number of unique samples per catergory in the training set
n.test <- 200 # number of unique samples per catergory in the testing set
N.waking <- 3000 # total size of the WAKING training set

n.wake <- 200 # number of waking trials per cycle
n.sleep <- 40 # number of REM trials per cycle
rec_mag <- 1 # recency weighting curve
#   0 means all waking samples have same weight
#   approaching inf means only most recent has weight
#   effects rem.syn, rem.rep, rem.repn 
p.noise <- 0.2 # p.noise bits will be flipped per replay
#   effects rem.repn
p.nrem_reduce <- 0.95 # forget rate for nrem portion

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

# set later
train_set <- sets$train$rem.syn # rem.syn arbitarily chosen at this point

l.train_set <- sets$train$l.waking
N.total <- nrow(train_set)

# exp params
n.cycle <- n.wake + n.sleep # a cycle is a wake-sleep cycle
N.cycle <- N.total / n.cycle

# init network
n.vis <- ncol(train_set) # number of visible units
W0 <- matrix(rnorm(n.hid*n.vis,mean=0,sd=sd_init_weights),
            nrow=n.vis,ncol=n.hid) # vis-hid weight matrix

condnames <- setdiff(names(sets$train),'l.waking')
n.cond <- length(condnames)

# init stats collection dataframe
nrow.metrics <- (perm(n.cat,2) + n.cat) * 4 * 2# four is 2 recon (mahal and euc) * 2 hidden (mahal and euc). 2 is times a day
names_metrics <- c('condition','cycle','time_of_day','from','to','metric','mu','std')
ncol.metrics <- length(names_metrics)
metrics <- as.data.frame(matrix(NA,nrow=nrow.metrics*N.cycle*n.cond,ncol=ncol.metrics))
nrow.metrics.percond <- nrow(metrics) / n.cond
names(metrics) <- names_metrics


######### CONDUCT EXPERIMENT #########
i <- 1 # i = ith training example
inds.mahal <- seq(nrow.metrics)
i.nrem <- round(seq(0,3/4,1/4)*n.sleep)+n.wake

inds.metrics <- seq(nrow.metrics)

base <- 0 # increases after every condition has fully run (for adding to metrics dataframe properly)
for(COND in condnames){
    
    base <- base + nrow.metrics.percond
    
    print(cat('condition',as.character(COND),'of',as.character(n.cond)))
    
    # get training set for this condition
    train_set <- sets[['train']][[COND]]
        
    i.cycle <- 0
    
    W <- W0
    
    for(i in seq(N.total)){
        
        ## TRAIN NETWORK
        W[] <- rbm_train(W,matrix(train_set[i,],nrow=1),lrate=lrate)
        
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
            
            # get all the metrics!
            metrics0 <- nnet_metrics_combined(W,test_set,l.test_set)
            
            # STORE DATA
            metrics[inds.metrics,4:ncol.metrics] <- metrics0
            metrics[inds.metrics,'condition'] <- COND
            metrics[inds.metrics,'cycle'] <- i.cycle
            metrics[inds.metrics,'time_of_day'] <- time_of_day
            
            
            inds.metrics[] <- inds.metrics + nrow.metrics + base
        }
        
        ## NREM (weight reduction)
        if(i%%n.cycle %in% i.nrem){
            W[] <- p.nrem_reduce*W
        }
    }
}
# 

write.csv(W,'W.csv')
write.csv(metrics,'metrics.csv')
# # # visualize
# g <- ggplot(subset(metrics,metric='hm'))
# 
# g + geom_point(aes(x=cycle,y=mu,colour=condition)) + scale_y_log10()
# g + geom_point(aes(x=to,y=from,size=mu))
