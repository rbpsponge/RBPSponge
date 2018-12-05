rm(list=ls(all=TRUE))
library(ggplot2)
library(DAAG)
library(lmtest)
# devtools::install_github('gokceneraslan/DAAG')
corr_RBP_lncRNA <- function(existing_list, rbp_data, lncRNA_data){
  num_genes = length(existing_list)
  list_corrs_RBP_CV = numeric(length = num_genes)
  list_corrs_RBPlncRNA_CV = numeric(length = num_genes)
  lr_test_pvalues = numeric(length = num_genes)
  i = 1
  for (gene in existing_list){
    gene_data = exp_data[gene,skip_col:ncol(exp_data)];

    if (RBP == "PUM2")
    {
      data_RBP = data.frame(gene_data, rbp_data, PUM1_data);
      cvres_RBP = cv.lm(data_RBP, formula(gene_data ~ rbp_data + PUM1_data), m=10, plotit = FALSE, printit = FALSE);
      data_RBP_lncRNA = data.frame(gene_data, rbp_data, lncRNA_data, PUM1_data);
      cvres_RBP_lncRNA = cv.lm(data_RBP_lncRNA, formula(gene_data ~ rbp_data + PUM1_data + lncRNA_data), m=10,plotit = FALSE, printit = FALSE);
      RBP_lm = lm(gene_data ~ rbp_data + PUM1_data, data = data_RBP)
      RBP_lncRNA_lm = lm(gene_data ~ rbp_data + PUM1_data + lncRNA_data, data = data_RBP_lncRNA);
      
      lr_test_result = lrtest(RBP_lm, RBP_lncRNA_lm)$`Pr(>Chisq)`[2]
    }
    else
    {
      data_RBP = data.frame(gene_data, rbp_data)  ;
      cvres_RBP = cv.lm(data_RBP, formula(gene_data ~ rbp_data), m=10, plotit = FALSE, printit = FALSE);
      data_RBP_lncRNA = data.frame(gene_data, rbp_data, lncRNA_data);
      cvres_RBP_lncRNA = cv.lm(data_RBP_lncRNA, formula(gene_data ~ rbp_data + lncRNA_data), m=10, plotit = FALSE, printit = FALSE);
      RBP_lm = lm(gene_data ~ rbp_data, data = data_RBP)
      RBP_lncRNA_lm = lm(gene_data ~ rbp_data + lncRNA_data, data = data_RBP_lncRNA);
      # likelihood ratio test
      lr_test_result = lrtest(RBP_lm, RBP_lncRNA_lm)$`Pr(>Chisq)`[2]
    }
    corr_RBP_CV = cor(cvres_RBP$cvpred, gene_data, method = "spearman");
    list_corrs_RBP_CV[i] <- round(corr_RBP_CV, digits = 5);
    corr_RBPlncRNA_CV = cor(cvres_RBP_lncRNA$cvpred, gene_data, method = "spearman");
    list_corrs_RBPlncRNA_CV[i] <-  round(corr_RBPlncRNA_CV,5);
    lr_test_pvalues[i] = lr_test_result;
    i = i + 1;
  }

  return (list(list_corrs_RBP_CV, list_corrs_RBPlncRNA_CV, lr_test_pvalues))
}


args = commandArgs(trailingOnly=TRUE)
lncRNA = args[1]
RBP = args[2]
dataset = args[3]
job_id = args[4]

if(length(args) > 4) # user supplied target  file
{
  targetfile = args[5]
}

infofileh <- file(paste(job_id, '_analysis2_info.txt', sep = ''), open="wt")
errorfileh <- file(paste(job_id, '_analysis2_error.txt', sep = ''), open="wt")
sink(infofileh, type="output", append=TRUE) # if this is empty, we should display 'unexpected error occurred'
sink(errorfileh, type="message", append=TRUE)

targetdir = '/results/'
datadir = '/datasets/'
outdir = '/results/'

present = TRUE
if(startsWith(dataset, 'EMTAB2706')){
  load(paste(datadir, 'EMTAB2706_mapping_frame.RData', sep = ''))
  if(!(lncRNA %in%  mapping_frame$ENSG_id) || !(RBP %in%  mapping_frame$gene_symbol)){
    print('ERROR: this gene do not exist in this expression data')
    present = FALSE
  }else{
    skip_col = 1
    filename = paste(datadir, dataset, '.RData', sep = '')
    if (dataset == "EMTAB2706"){
      filename = paste(datadir, dataset, '_tpm.RData', sep = '')
    }
    
    if (file.exists(filename)){
      load(filename)
    }else{
      print("ERROR: cannot open input data file")
    }
  }
  if(startsWith(dataset, 'EMTAB2770')){
  load(paste(datadir, 'EMTAB2770_mapping_frame.RData', sep = ''))
  if(!(lncRNA %in%  mapping_frame$ENSG_id) || !(RBP %in%  mapping_frame$gene_symbol)){
    print('ERROR: this gene do not exist in this expression data')
    present = FALSE
  }else{        
    skip_col = 1
    filename = paste(datadir, dataset, '.RData', sep = '')
    if (dataset == "EMTAB2770"){
      filename = paste(datadir, dataset, '_tpm.RData', sep = '')
    }
    if (file.exists(filename)){
      load(filename)
    }else{
      print("ERROR: cannot open input data file")
    }
  }
} else if(startsWith(dataset, 'GTEX')){
  load(paste(datadir, 'GTEX_mapping_frame.RData', sep = ''))
  if(!(lncRNA %in%  mapping_frame$ENSG_id) || !(RBP %in%  mapping_frame$gene_symbol)){
    print('ERROR: this gene do not exist in this expression data')
    present = FALSE
  }else{ 
    skip_col = 1
    filename = paste(datadir, dataset, '.RData', sep = '')    
    if (file.exists(filename)){
      load(filename)
      exp_data = data.matrix(exp_data)
    }else{
      print("ERROR: cannot open input data file")
    }
  }
} else if (startsWith(dataset, 'User'))# user uploaded data, dataset should refer to the name of the file
{
  # this will already have ensg ids  
  skip_col = 1
  filename = paste(datadir, dataset, '.RData', sep = '')
  if (file.exists(filename)){
    load(filename)
    if(!(lncRNA %in%  rownames(exp_data))){
      present = FALSE
    }
    
  }else{
    print("ERROR: cannot open input data file")
  }      
}

if(present){

  RBP_ENSG_id = toString(mapping_frame[mapping_frame$gene_symbol == RBP,]$ENSG_id)

  if (RBP == 'PUM2'){

    PUM1_ENSG_id = 'ENSG00000134644'
    PUM1_data = exp_data[PUM1_ENSG_id,skip_col:ncol(exp_data)];
  }


  if(length(args) > 4){
	  path_target = c("/var/www/laravel/storage/app/results/",args[5])
    path_target = paste(path_target, collapse="")
    if (file.exists(path_target)){

        target_list = readLines(path_target)
    }
    else{
      print("ERROR: user supplied target file do not exist")
    }
  } else {

    target_list = readLines(paste(targetdir, RBP ,'_target_genes.txt', sep = ''))
    nontarget_list =  readLines(paste(targetdir, RBP ,'_nontarget_genes.txt', sep = ''))
    }


  existing_target_list = intersect(target_list, row.names(exp_data))
  num_targets = length(existing_target_list)
  
  existing_nontarget_list = intersect(nontarget_list, row.names(exp_data))
  num_nontargets = length(existing_nontarget_list)

  #existing_nontarget_list = intersect(nontarget_list, row.names(exp_data))
  #num_nontargets = length(existing_nontarget_list)
  lncRNA_data = exp_data[lncRNA,skip_col:ncol(exp_data)]
  rbp_data = exp_data[RBP_ENSG_id,skip_col:ncol(exp_data)];
  result_target = corr_RBP_lncRNA(existing_target_list, rbp_data, lncRNA_data);
  target_corrs_RBP = result_target[[1]];
  target_corrs_RBPlncRNA = result_target[[2]];
  target_lr_pvalues = result_target[[3]]

  result_nontarget = corr_RBP_lncRNA(existing_nontarget_list, rbp_data, lncRNA_data);
  nontarget_lr_pvalues = result_nontarget[[3]]
  #nontarget_corrs_RBP = result_nontarget[[1]];
  #nontarget_corss_RBPlncRNA = result_nontarget[[2]];

  print(length(which(target_lr_pvalues < 0.01)) / length(target_list))
  print(length(which(nontarget_lr_pvalues < 0.01)) / length(nontarget_list))
  
  print(length(which(target_lr_pvalues < 0.05)) / length(target_list))
  print(length(which(nontarget_lr_pvalues < 0.05)) / length(nontarget_list))
  
  # convert these to percentages
  target_lr_01 =  signif(length(which(target_lr_pvalues < 0.01)) / length(target_list), digits = 2) * 100
  nontarget_lr_01 = signif(length(which(nontarget_lr_pvalues < 0.01)) / length(nontarget_list), digits = 2) * 100
  
  
  #targets VS RBPlncRNA
  dataframe_a1 <- data.frame(group = "RBP only" , value = target_corrs_RBP);
  dataframe_b1 <- data.frame(group = "RBP + lncRNA" , value = target_corrs_RBPlncRNA);
  plot.data1 <- rbind(dataframe_a1,dataframe_b1);

  data.print = data.frame('gene_id' = existing_target_list,'RBP_only_corr' = target_corrs_RBP,'RBP_and_lncRNA_corr' = target_corrs_RBPlncRNA)

   if(length(args) > 4){
        outfilekey = paste(outdir,job_id,  '_',RBP, '_', lncRNA, '_analysis2_', dataset, sep = '')

   }
    else
    {
        outfilekey = paste(outdir, RBP, '_', lncRNA, '_analysis2_', dataset, sep = '')
    }
  

  write.table(data.print, file = paste(outfilekey, '.txt', sep = ''),  quote = FALSE, sep = "\t", row.names = FALSE)

  result = wilcox.test(target_corrs_RBP, target_corrs_RBPlncRNA);
  pvalue <- result$p.value;


  ggplot(plot.data1, aes(x=group, y=value, fill= group)) + scale_fill_brewer(palette="Dark2") + geom_boxplot() + theme(axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size = rel(1.75)), axis.title.y=element_text(size = rel(1.5), angle = 90)) + labs(y = "Spearman correlation coefficient") +  annotate("text", x=1.5, y=-0.2, label= paste('Wilcox. p-val: ', signif(pvalue, digits = 3))) + annotate("text", x=1.5, y=-0.1, label= paste('LR test: ',target_lr_01, '% vs ', nontarget_lr_01, '%', sep = ''))

  
  filename1 = paste(outfilekey, '.jpeg' , sep= '');
  ggsave(filename1, plot = last_plot())


sink(type="message")
sink(type="output")
}



