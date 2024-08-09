filter_ds = function(ds,ls){
    ds %>% filter(
            virus   %in% ls$virus   &
            antigen %in% ls$antigen &
            vaccine %in% ls$vaccine &
            group   %in% ls$group &
            ls$DPFD_fn(DPFD) &
            ls$DPSD_fn(DPSD) &
            ls$DPTD_fn(DPTD)
        )
}

get_same<-function(cnt,exp){
    paste(na.omit(sapply(names(cnt[1:4]),function(ni){
        if( cnt[[ni]] == exp[[ni]] ) {
            return(paste0(ni,':',cnt[[ni]]))
        }
        return(NA)
    })),collapse=', ')
}

get_diff<-function(cnt,exp){
    paste(na.omit(sapply(names(cnt[1:4]),function(ni){
        if( cnt[[ni]] != exp[[ni]] ) {
            return(paste0(cnt[[ni]],'_',exp[[ni]]))
        }
        return(NA)
    })),collapse=', ')
}

fit_lmer<-function(data_i,cnt,exp,compare_var,time_var,subj_var='Lab.ID',n='ngrps'){

    formula = as.formula(paste('log(value) ~',compare_var,'*',time_var,' +(1|',subj_var,')'))
#     formula = as.formula(paste('log(value) ~',compare_var,'*',time_var,' +(1+',time_var,'|',subj_var,')'))

    ds=as.data.frame(coef(summary(fitI<-lme4::lmer(formula,data = data_i))))%>% 
            tibble::rownames_to_column('parameter') %>%
            mutate(
                vaccine=paste(unique(c(cnt$vaccine,exp$vaccine)),collapse='_'),
                antigen=paste(unique(c(cnt$antigen,exp$antigen)),collapse='_'),
                virus=paste(unique(c(cnt$virus,exp$virus)),collapse='_'),
                time=time_var,
                ngrps = summary(fitI)$ngrps,
                n=length(summary(fitI)$residuals))

    colnames(ds) = gsub(' ','',colnames(ds))
    
    if(!n%in%c('n','ngrps')){stop('n must be "n" or "ngrps"')}
    ds$Pr_t = 2*pt(-abs(ds$tvalue),df = ds[[n]],lower.tail = TRUE)
        
    list(ds,fitI)
}



constructCIRibbon <- function(model, newdata=NULL,z=0.95,formula=NULL) {

    if(is.null(newdata)){
      newdata = model@frame  
    }

    df = newdata
    if(is.null(df$fit)){
        df$fit = predict(model, newdata = newdata)
    }

    if(is.null(formula)){
    formula = as.formula( paste( '~', gsub( '\\+.*','', as.character(model@call$formula)[3] )))
    }
    
    mm <- model.matrix(formula, data = df)

    vars <- mm %*% vcov(model) %*% t(mm)
#      print(vars)
    sds <- sqrt(diag(vars))
#      print(sds)
    
    df=df %>% 
        mutate(
            fit.lw = fit - qnorm(z) * sds,
            fit.hi = fit + qnorm(z) * sds#,
#             across(fit:fit.hi, exp)
        )
#     print(head(df))
    df$fit = exp(df$fit)
    df$fit.hi = exp(df$fit.hi)
    df$fit.lw = exp(df$fit.lw)
    df
}

get_projection_uni <- function(cf,time){
    exp(cf[1] + cf[2]*time )
}


get_projection <- function(cf,comp_num,time){
    exp(cf[1]      +cf[2]*comp_num +
        cf[3]*time +cf[4]*comp_num*time)
}



#########################################3

### main function
compare_decay = function(data_i,cnt,exp,compare_var,time_var,subj_var='Lab.ID',
                         null_ls=list(vaccine='Naive',
                            DPFD_fn=function(x) is.na(x),
                            DPSD_fn=function(x) is.na(x)),
                         trend_type = 'variance_based', # or 'simple'
                         plot=F,run_dbg=F,trim_incs=T,z=0.95){
                             
    ### Init null distribution
    nld = cnt
    for(inull in names(null_ls)){
        nld[[inull]] = null_ls[[inull]]
    }
                             
    null_data = filter_ds(data_i,nld)

    # get quantiles and IQR multiples for the null distribution
    qt1 = quantile( null_data$value , na.rm=T )
    iqr1 = qt1[4] - qt1[2]
    null_thresh = list(median=qt1[3]*1,
               'M+IQR*5' =qt1[3]+iqr1*5,
               'M+IQR*10'=qt1[3]+iqr1*10)

    ### Compile and trim data
    data_i = rbind(
        filter_ds(data_i,cnt),
        filter_ds(data_i,exp)
    )  %>% droplevels()
   

    ### Init stratification
    data_i$strat = paste(data_i[[subj_var]],data_i[[compare_var]])
    
    ### Drop increasing values
    if(trim_incs){
        
        data_i = data_i[order(data_i[[time_var]]),] # sort by time variable
        
        for(subj in unique(data_i[[subj_var]])){
            idx_subj = which(data_i[[subj_var]]==subj)
            v_si = data_i$value[idx_subj] # subj-specific values
            t_si = data_i[[time_var]][idx_subj] # subj-specific values
            is_inc = ifelse(t_si>100,
                            c(FALSE,diff(v_si) > 0),
                            c(diff(v_si) > 0,FALSE)
                            )
            data_i$value[idx_subj[is_inc]] = NA
        }
    }
    
#     ### drop single observations 
#     tab = table(data_i[[subj_var]])
#     rm_sbj = names(tab)[tab<2]
#     data_i = data_i[ !data_i[[subj_var]]%in%rm_sbj , ]
    
    if(run_dbg){
        var_dbg<<-list(data_i,cnt,exp,compare_var,time_var)
    }
    
    ###Fit model
    model_out <-  fit_lmer(data_i,cnt,exp,compare_var,time_var,subj_var,'ngrps')
    fitI = model_out[[2]]
                             
    # simple regression projection
    cf = coef(summary(fitI))[,1]
    cf_hi = coef(summary(fitI))[,1] + coef(summary(fitI))[,2]*qnorm(z)
    cf_lw = coef(summary(fitI))[,1] - coef(summary(fitI))[,2]*qnorm(z)
    
    ### Get factors
    if(!is.factor(data_i[[compare_var]])){
        data_i[[compare_var]] = factor( data_i[[compare_var]] )
    }
    
    ### Get trend lines
    comp_num = as.numeric(data_i[[compare_var]])-1
    if(!any(comp_num%in%c(0,1))){stop("comp_num should be 0 or 1. Are there more than 2 variables in comp_var?")}
                             
    if(trend_type=='variance_based'){
        
        data_i$fit =    log(get_projection(cf,    comp_num,data_i[[time_var]]))
        data_i = constructCIRibbon(fitI,data_i,z=z)
        data_i$fit =    (get_projection(cf,    comp_num,data_i[[time_var]]))

                
#         thresh=120
#         data_i$fit[data_i[[time_var]]>thresh] = NA
#         data_i$fit.lw[data_i[[time_var]]>thresh] = NA
#         data_i$fit.hi[data_i[[time_var]]>thresh] = NA
        
    }else if(trend_type=='simple'){

        data_i$fit    = get_projection(cf,    comp_num,data_i[[time_var]])
        data_i$fit.hi = get_projection(cf_hi, comp_num,data_i[[time_var]])  
        data_i$fit.lw = get_projection(cf_lw, comp_num,data_i[[time_var]]) 
        
    }else{
        stop('trend_type can be simple or variance_based')
    }
    
    # remove fit values beyond observed range
    data_i$fit[data_i[[time_var]]>max(data_i[[time_var]][comp_num==0]) & comp_num==0] = NA
    data_i$fit[data_i[[time_var]]>max(data_i[[time_var]][comp_num==1]) & comp_num==1] = NA

    # spagetti plot
    g1=data_i %>%
        ggplot(aes_string(x=time_var,y='value',shape=compare_var,lty=compare_var))+
            geom_point()+
            geom_line(aes_string(group='strat',color='strat'))+ # patient-specific lines
    
            geom_line(aes(y=fit),size=3,color='black')+         # trend lines
#         geom_ribbon(aes_string(ymin='lowCI',ymax='HiCI',y='Predict',fill=compare_var),alpha=.2) +
            geom_ribbon(aes_string(ymin='fit.lw',ymax='fit.hi',y='fit',fill=compare_var),alpha=.25) +
    
            theme_bw(base_size = 15)+
            theme(legend.key.width = unit(2.5,"cm"))+
            guides(color = 'none')+ ylab('Antigen-specific IgG1 Titer')+
            ggtitle(title<-paste0('Compare ',compare_var,': ',get_diff(cnt,exp),' over ',time_var),
                    subtitle = (sub<-get_same(cnt,exp)))+
#             geom_hline(aes(yintercept=null_thresh[[1]]))+
#                 geom_text(aes(0,null_thresh[[1]],label = names(null_thresh)[1], vjust = -1))+
                geom_hline(aes(yintercept=null_thresh[[2]]))+
                geom_text(aes(0,null_thresh[[2]],label = names(null_thresh)[2], vjust = -.25))+
                geom_hline(aes(yintercept=null_thresh[[3]]))+
                geom_text(aes(0,null_thresh[[3]],label = names(null_thresh)[3], vjust = -.25))

#     }

    if(plot){
        print(g1)
        anyl = paste0('Compare_',compare_var,'_',get_diff(cnt,exp),'.over_',time_var,'/')
        dir.create(paste0('../results//03-Generalize_differential_dose/',anyl))
        spec = paste0(gsub(',','',gsub(' ','.',gsub(':','_',get_same(cnt,exp)))),'.pdf')
        spec_lg = paste0(gsub(',','',gsub(' ','.',gsub(':','_',get_same(cnt,exp)))),'.log10.pdf')
        ggsave(g1,filename=paste0('../results//03-Generalize_differential_dose/',anyl,spec),
               height=10,width=10)
        g2 = g1+scale_y_log10()
        print(g2)
        ggsave(g2,filename=paste0('../results//03-Generalize_differential_dose/',anyl,spec_lg),
               height=10,width=10)
    }

    if(run_dbg){
        var_dbg<<-list(data_i,cnt,exp,compare_var,time_var)
    }
    
    return(model_out[[1]])
}