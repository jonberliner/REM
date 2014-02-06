# repmat = function(X,m,n){
#     ##R equivalent of repmat (matlab)
#     mx = dim(X)[1]
#     nx = dim(X)[2]
#     matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
# }

n.times <- 10
# sweep of params
n.cat  <- c(2) # number of categories
n.train <- c(seq(10,100,10),seq(200,1000,100),seq(2000,3000,1000)) # number of unique samples per catergory in the training set
n.test <- 1000 # number of unique samples per catergory in the testing set
N.waking <- 10000 # total size of the WAKING training set

n.wake <- 500 # number of waking trials per cycle
n.sleep <- 50 # number of REM trials per cycle
rec_mag <- 1 # recency weighting curve
#   0 means all waking samples have same weight
#   approaching inf means only most recent has weight
#   effects rem.syn, rem.rep, rem.repn 
p.noise <- 0.1 # p.noise bits will be flipped per replay
#   effects rem.repn

######### NETWORK FREE PARAMS #########
n.hid <- c(5,seq(10,80,20),seq(100,100,1000)) # number of hidden nodes
l.rate <- c(0.05,0.1) # network learning rate
sd_init_weights <- 0.05 # sd of normally distd random net weights

n <- length(n.cat)*length(n.train)*length(p.noise)*length(n.hid)*length(l.rate)*n.times
colnames <- c('ncat','ntrain','ntest','ntotalwaking','nwake','nsleep','rec_mag','pnoise','nhid','lrate','sd_init_weights')
out <- as.data.frame(matrix(NA,nrow=n,ncol=length(colnames)))
names(out) <- colnames

i <- 1
for(c in n.cat){
    for(tr in n.train){
        for(n in p.noise){
            for(h in n.hid){
                for(l in l.rate){
                    for(ii in seq(n.times)){
                        out[i,'ncat'] <- c
                        out[i,'ntrain'] <- tr
                        out[i,'ntest'] <- n.cat
                        out[i,'ntotalwaking'] <- N.waking
                        out[i,'nwake'] <- n.wake
                        out[i,'nsleep'] <- n.sleep
                        out[i,'rec_mag'] <- rec_mag
                        out[i,'pnoise'] <- p.noise
                        out[i,'nhid'] <- h
                        out[i,'lrate'] <- l
                        out[i,'sd_init_weights'] <- sd_init_weights
                        
                        i <- i+1
                    }
                }
            }
        }
    }
}

#out <- repmat(out,n.per_paramset,1)

write.csv(out,'~/Desktop/test.csv')
