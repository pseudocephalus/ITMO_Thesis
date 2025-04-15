library(dplyr)
library(stringr)
library(SVDFunctions)

#https://github.com/na89/weightedLR
linearReg <- function(cases, controls) {
  table <- data.frame(geno = c(0, 1, 2, 0, 1, 2), phe = c(0, 0, 0, 1, 1, 1))
  ws <- c(controls, cases)
  m <- lm(phe ~ geno, data = table, weights = ws)
  m$df.residual <- sum(ws) - 2
  ks <- summary(m)$coefficients
  list(beta = ks[2, 1], pval = ks[2, 4])
}

#https://www.rdocumentation.org/packages/GenABEL/versions/1.8-0/topics/estlambda
"estlambda" <- function(data, plot=FALSE, proportion=1.0,
                        method="regression", filter=TRUE, df=1,... ) {
  data <- data[which(!is.na(data))]
  if (proportion>1.0 || proportion<=0)
    stop("proportion argument should be greater then zero and less than or equal to one")
  
  ntp <- round( proportion * length(data) )
  if ( ntp<1 ) stop("no valid measurements")
  if ( ntp==1 ) {
    warning(paste("One measurement, lambda = 1 returned"))
    return(list(estimate=1.0, se=999.99))
  }
  if ( ntp<10 ) warning(paste("number of points is too small:", ntp))
  if ( min(data)<0 ) stop("data argument has values <0")
  if ( max(data)<=1 ) {
    #		lt16 <- (data < 1.e-16)
    #		if (any(lt16)) {
    #			warning(paste("Some probabilities < 1.e-16; set to 1.e-16"))
    #			data[lt16] <- 1.e-16
    #		}
    data <- qchisq(data, 1, lower.tail=FALSE)
  }
  if (filter)
  {
    data[which(abs(data)<1e-8)] <- NA
  }
  data <- sort(data)
  ppoi <- ppoints(data)
  ppoi <- sort(qchisq(ppoi, df=df, lower.tail=FALSE))
  data <- data[1:ntp]
  ppoi <- ppoi[1:ntp]
  #	s <- summary(lm(data~offset(ppoi)))$coeff
  #       bug fix thanks to Franz Quehenberger
  
  out <- list()
  if (method=="regression") {
    s <- summary( lm(data~0+ppoi) )$coeff
    out$estimate <- s[1,1]
    out$se <- s[1,2]
  } else if (method=="median") {
    out$estimate <- median(data, na.rm=TRUE)/qchisq(0.5, df)
    out$se <- NA
  } else if (method=="KS") {
    limits <- c(0.5, 100)
    out$estimate <- estLambdaKS(data, limits=limits, df=df)
    if ( abs(out$estimate-limits[1])<1e-4 || abs(out$estimate-limits[2])<1e-4 )
      warning("using method='KS' lambda too close to limits, use other method")
    out$se <- NA
  } else {
    stop("'method' should be either 'regression' or 'median'!")
  }
  
  if (plot) {
    lim <- c(0, max(data, ppoi,na.rm=TRUE))
    #		plot(ppoi,data,xlim=lim,ylim=lim,xlab="Expected",ylab="Observed", ...)
    oldmargins <- par()$mar
    par(mar=oldmargins + 0.2)
    plot(ppoi, data,
         xlab=expression("Expected " ~ chi^2),
         ylab=expression("Observed " ~ chi^2),
         ...)
    abline(a=0, b=1)
    abline(a=0, b=out$estimate, col="red")
    par(mar=oldmargins)
  }
  
  out
}

plot_lam <- function(pvals, title = 'title'){
  p_values <- pvals
  
  observed <- -log10(sort(p_values))
  n <- length(p_values)
  expected <- -log10(ppoints(n))  
  
  ci <- 0.95
  lower <- -log10(qbeta((1 - ci)/2, 1:n, n:1))
  upper <- -log10(qbeta((1 + ci)/2, 1:n, n:1))
  
  plot(expected, observed, 
       xlab = "Expected (-log₁₀ p)", 
       ylab = "Observed (-log₁₀ p)", 
       main = title,
       pch = 19, col = "blue", cex = 0.6)
  
  abline(0, 1, col = "red", lwd = 2)
  
  lines(expected, lower, lty = 2, col = "gray")
  lines(expected, upper, lty = 2, col = "gray")
  
  chi_sq <- qchisq(p_values, df = 1, lower.tail = FALSE)
  lambda <-estlambda(p_values)$estimate
  text(x = max(expected) - 1, y = 1, 
       labels = paste0("λ = ", round(lambda, 3)), 
       col = "red", cex = 0.8)
}

read_controls <- function(path){
  read.csv(path, sep = '\t', skip = 3) %>%
    mutate(chr_pos = str_extract(X, '(chr[0-9X]*.*?)\t', group = 1),
           ref_gt  = str_extract(X, 'chr[0-9X]*:[0-9]*\t(.*)\t.*$', group = 1),
           alt_gt  = str_extract(X, 'chr[0-9X]*:[0-9]*\t.*\t(.*)$', group = 1)
    ) %>% 
    mutate(chr_pos_gt =paste0(chr_pos, ':', ref_gt, ' ', alt_gt)) %>% 
    group_by(chr_pos_gt) %>% 
    summarise(hom_ref = sum(hom_ref),
              het = sum(het),
              hom_alt = sum(hom_alt)) %>% 
    filter(checkAlleleCounts(matrix(c(hom_ref, het, hom_alt),ncol = 3), maf = 0.05, mac=10)) 
}


# Imputed
## Cluster 1
controls_imp1 <- read_controls('controls_counts/counts_imputed_1.tsv') 
cases_imp1 <- read_cases('case_counts/counts_imputed_1') 
res_imp1 <- do_regr(cases_imp1, controls_imp1)
#plot_lam(res_imp1$p_regr, 'Cluster 1; imputed variants')
## Cluster 2
controls_imp2 <- read_controls('controls_counts/counts_imputed_2.tsv') 
cases_imp2<- read_cases('case_counts/counts_imputed2') 
res_imp2 <- do_regr(cases_imp2, controls_imp2)
## Cluster 3
controls_imp3 <- read_controls('controls_counts/counts_imputed_3.tsv') 
cases_imp3 <- read_cases('case_counts/counts_imputed3') 
res_imp3 <- do_regr(cases_imp3, controls_imp3)

# Typed
## Cluster 1
controls_typed1 <- read_controls('controls_counts/counts_typed_1.tsv') 
cases_typed1<- read_cases('case_counts/counts_typed_1') 
res_typed1 <- do_regr(cases_typed1, controls_typed1)
## Cluster 2
controls_typed2 <- read_controls('controls_counts/counts_typed_2.tsv') 
cases_typed2<- read_cases('case_counts/counts_typed_2') 
res_typed2 <- do_regr(cases_typed2, controls_typed2)
## Cluster 3
controls_typed3 <- read_controls('controls_counts/counts_typed_3.tsv') 
cases_typed3 <- read_cases('case_counts/counts_typed_3') 
res_typed3 <- do_regr(cases_typed3, controls_typed3)

# output
res_imp1 <- res_imp1 %>% 
  mutate(type = 'imputed',
         cluster = 1)
res_imp2 <- res_imp2 %>% 
  mutate(type = 'imputed',
         cluster = 2)
res_imp3 <- res_imp3 %>% 
  mutate(type = 'imputed',
         cluster = 3)
res_typed1 <- res_typed1 %>% 
  mutate(type = 'typed',
         cluster = 1)
res_typed2 <- res_typed2 %>% 
  mutate(type = 'typed',
         cluster = 2)
res_typed3 <- res_typed3 %>% 
  mutate(type = 'typed',
         cluster = 3)

res_full<- rbind(res_imp1, res_imp2, res_imp3, 
                 res_typed1, res_typed2, res_typed3) 


anns <- read.csv('variant_annotations.csv', header = F)
colnames(anns) <- c('snp', 'chr', 'pos', 'ref',  'alt', 'R2', 'MAF', 'CSQ')


anns %>% 
  mutate(chr_pos_gt = paste0(chr,':', pos, ':',
                             ref, ' ', alt),
         CSQ = str_extract(CSQ, '.*?\\|(.*?)\\|', group = 1))  %>% 
  select(snp,chr_pos_gt, R2, MAF, CSQ) -> anns


df <- merge(res_full, anns, by = 'chr_pos_gt')
write.csv(df, 'results_full_annotated.csv')
