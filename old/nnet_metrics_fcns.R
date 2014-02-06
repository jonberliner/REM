
###### BETWEEN CAT MAHAL METRICS ######
mahal_metrics <- function(W,hidprobs.all,l.data){
    
    cats <- unique(c(l.data))
    n.cat <- length(cats)
    catnames <- as.character(cats)
    
    
    mu.hid.all <- colMeans(hidprobs.all)
    cv.hid.all <- cov(hidprobs.all)
    
    mu.hid <- llply(cats,function(c){
        colMeans(hidprobs.all[l.data==c,])
    })
    names(mu.hid) <- catnames
    
    cv.hid <- llply(cats,function(c){
        cov(hidprobs.all[l.data==c,])
    })
    names(cv.hid) <- catnames
    
    # mahal dist for every between pair, and to whole dataset.  no use doing within because, by def, will be 1*nunits
    mahal_b2a <- llply(catnames,function(c1){
        llply(catnames,function(c2){
            if(c1!=c2){
                mahalanobis(hidprobs.all[l.data==as.numeric(c2),],center=mu.hid[[c1]],cov=cv.hid[[c1]])
            }
        })
    })
    
    # label
    names(mahal_b2a) <- catnames
    for(c in catnames){
        names(mahal_b2a[[c]]) <- catnames
    }
    
    
    # get summary stats
    perm <- function(n,k){choose(n,k) * factorial(k)}
    nperm <- perm(n.cat,2) # minue n.cat because don't have from all to things
    df <- as.data.frame(matrix(nrow=nperm,ncol=4))
    names(df) <- c('to','from','mu','std')
    ir <- 1
    for(c1 in catnames){
        for(c2 in catnames){
            if(c1!=c2){
                df[ir,'to'] <- c1
                df[ir,'from'] <- c2
                df[ir,'mu'] <- mean(mahal_b2a[[c1]][[c2]])
                df[ir,'std'] <- sd(mahal_b2a[[c1]][[c2]])
                
                ir <- ir+1
            }
        }
    }
    return(df)   
}


###### BETWEEN AND WITHIN CAT EUC METRICS ######
euc_metrics <- function(W,hidprobs.all,l.data){
    
    cats <- unique(c(l.data))
    n.cat <- length(cats)
    catnames <- as.character(cats)
    
    mu.hid.all <- colMeans(hidprobs.all)
    
    mu.hid <- llply(cats,function(c){
        colMeans(hidprobs.all[l.data==c,])
    })
    names(mu.hid) <- catnames
    
    
    # euc dist for every between AND within pair
    euc_b2a <- llply(catnames,function(c1){
        llply(catnames,function(c2){
            sum(sqrt((hidprobs.all[l.data==as.numeric(c2),] - mu.hid[[c1]])^2))
        })
    })
    
    # label
    names(euc_b2a) <- catnames
    for(c in catnames){
        names(euc_b2a[[c]]) <- catnames
    }
    
    
    # get summary stats
    perm <- function(n,k){choose(n,k) * factorial(k)}
    nperm <- perm(n.cat,2) # minue n.cat because don't have from all to things
    df <- as.data.frame(matrix(nrow=nperm,ncol=4))
    names(df) <- c('to','from','mu','std')
    ir <- 1
    for(c1 in catnames){
        for(c2 in catnames){
            df[ir,'to'] <- c1
            df[ir,'from'] <- c2
            df[ir,'mu'] <- mean(euc_b2a[[c1]][[c2]])
            df[ir,'std'] <- sd(euc_b2a[[c1]][[c2]])        
            
            ir <- ir+1
        }
    }
    return(df)   
}




###### RECONSTRUCTION ERROR METRICS ######
recon_metrics <- function(W,hidprobs.all,l.data){
    temp <- fire(hidprobs.all)
    recon.all <- pV_aH(W,temp)
    
    # get euclidean reconstruction error
    rec_err.indiv <- llply(cats,function(c){
        rowSums( sqrt((test_set[l.test_set==c,] - recon.all[l.test_set==c,])^2) )
    })
    names(rec_err.indiv) <- catnames
    
    # get avg recon err
    return( ldply(catnames,function(c){
        data.frame('cat'=c,
                   'mu'=mean(rec_err.indiv[[c]]),
                   'std'=sd(rec_err.indiv[[c]]),
                   'cycle'=i.cycle,
                   'time_of_day'=time_of_day
        )
    }) )
}