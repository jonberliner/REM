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
    source('~/Dropbox/Sleep/datastream_fcns.R') # load fcns for waking waking+rem datastreams
    
    cats <- sample(c(1,2,3,4,5,6,7,8,9,0),n.cat) # randomly choose ncat categories
    waking_and_test <- make_wakingstream_and_testset(mnist,n.train,n.test,cats,N) # make waking data stream AND test set
    
    waking_rem_streams <- make_waking_rem_streams(waking_and_test$set.train,n.wake,n.sleep,rec_mag,p.noise) # make full waking+REM streams!
    return(list('test'=waking_and_test$set.test,'l.test'=waking_and_test$l.set.test,
                'train'=waking_rem_streams$set,'l.train'=waking_rem_streams$l.waking))
}