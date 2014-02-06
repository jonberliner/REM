###### BETWEEN MAHAL-DISTS, EUC DISTS, AND RECON ERR ######
nnet_metrics_combined <- function(W,data,l.data){
    
    cats <- unique(c(l.data))
    n.cat <- length(cats)
    catnames <- as.character(cats)
    
    # get hidden neuron probs and covariance
    hidprobs.all <- pH_aV(W,data) # get hidden layer activation for all samples
    
    mu.hid.all <- colMeans(hidprobs.all)
    cv.hid.all <- cov(hidprobs.all)
    icv.hid.all <- pseudoinverse(cv.hid.all)
    
    mu.hid <- llply(cats,function(c){
        colMeans(hidprobs.all[l.data==c,])
    })
    names(mu.hid) <- catnames
    
    icv.hid <- llply(cats,function(c){
        cv <- cov(hidprobs.all[l.data==c,])
        return(pseudoinverse(cv))
    })
    names(icv.hid) <- catnames
    
    # get reconstructions
    temp <- fire(hidprobs.all)
    recon.all <- pV_aH(W,temp)
    
    # get recon probs and covariance
    mu.rec <- llply(cats,function(c){
        colMeans(recon.all[l.data==c,])
    })
    names(mu.rec) <- catnames
    
    icv.rec <- llply(cats,function(c){
        cv <- cov(recon.all[l.data==c,])
        return(pseudoinverse(cv))
    })
    names(icv.rec) <- catnames
    
    euc <- function(a1,a2){
        return( rowSums( sqrt((a1-a2)^2) ) )
    }
    
    
#     mahal_else_na <- function(a,mu,cv){
#         out <- tryCatch({
#             mahalanobis(a,center=mu,cov=cv) # try to get mahal dist
#         }, error=function(cond){
#             NA
#         }, finally = {
#             
#         })
#         return(out)
#     }
    
    # mahal dist for every between pair, and to whole dataset.  no use doing within because, by def, will be 1*nunits
    dists_a2b <- llply(catnames,function(c1){
        llply(catnames,function(c2){
            list(
            #if(c1!=c2){    
                'hm' = mahalanobis(hidprobs.all[l.data==as.numeric(c1),], center=mu.hid[[c2]], cov=icv.hid[[c2]],inverted=TRUE), # hid mahal dist
                'rm' = mahalanobis(recon.all[l.data==as.numeric(c1),], center=mu.rec[[c2]], cov=icv.rec[[c2]],inverted=TRUE), # recon mahal dist
            #}
                'he' = euc(hidprobs.all[l.data==as.numeric(c1),], mu.hid[[c2]]), # hid euc dist
                're' = euc(recon.all[l.data==as.numeric(c1),], mu.rec[[c2]]) # recon euc dist
            )
        })
    })
    
    # label
    names(dists_a2b) <- catnames
    for(c in catnames){
        names(dists_a2b[[c]]) <- catnames
    }
    
    df <- ldply(catnames,function(c1){
        ldply(catnames,function(c2){
            a2b <- dists_a2b[[c1]][[c2]]
            ldply(names(a2b),function(m){
                data.frame('from'=c1,
                           'to'=c2,
                           'metric'=m,
                           'mu'=mean(a2b[[m]]),
                           'std'=sd(a2b[[m]])
            )}
        )}
    )})
    
    return(df)
}