rm(list=ls(all=TRUE))
library(ggplot2)
library(data.table)


args = commandArgs(trailingOnly=TRUE)
lncRNA = args[1]
RBP = args[2]
dataset = args[3]
job_id = args[4]


infofileh <- file(paste(job_id, '_analysis1_info.txt', sep = ''), open="wt")
errorfileh <- file(paste(job_id, '_analysis1_error.txt', sep = ''), open="wt")
sink(infofileh, type="output", append=TRUE) # if this is empty, we should display 'unexpected error occurred'
sink(errorfileh, type="message", append=TRUE)

targetdir = '/results/'
datadir = '/datasets/'
outdir = '/results/'

present = TRUE
if(startsWith(dataset, 'EMTAB2706')){
  load(paste(datadir, 'EMTAB2706_mapping_frame.RData', sep = ''))
  if(!(lncRNA %in%  mapping_frame$ENSG_id) && !(lncRNA %in%  mapping_frame$gene_symbol)){
    print('ERROR: this gene do not exist in this expression data')
    present = FALSE
  }else{
    if(!(lncRNA %in%  mapping_frame$ENSG_id) && (lncRNA %in%  mapping_frame$gene_symbol) )
      lncRNA = toString(mapping_frame[mapping_frame$gene_symbol == lncRNA,]$ENSG_id)
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
} else if(startsWith(dataset, 'EMTAB2770')){
  load(paste(datadir, 'EMTAB2770_mapping_frame.RData', sep = ''))
  if(!(lncRNA %in%  mapping_frame$ENSG_id) && !(lncRNA %in%  mapping_frame$gene_symbol)){
    print('ERROR: this gene do not exist in this expression data')
    present = FALSE
  }else{        
    if(!(lncRNA %in%  mapping_frame$ENSG_id) && (lncRNA %in%  mapping_frame$gene_symbol) )
      lncRNA = toString(mapping_frame[mapping_frame$gene_symbol == lncRNA,]$ENSG_id)
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
  if(!(lncRNA %in%  mapping_frame$ENSG_id) && !(lncRNA %in%  mapping_frame$gene_symbol)){
    print('ERROR: this gene do not exist in this expression data')
    present = FALSE
  }else{ 
    if(!(lncRNA %in%  mapping_frame$ENSG_id) && (lncRNA %in%  mapping_frame$gene_symbol) )
      lncRNA = toString(mapping_frame[mapping_frame$gene_symbol == lncRNA,]$ENSG_id)      
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
    if(length(args) == 6){
        path_target = c("/var/www/laravel/storage/app/results/",args[5])
        path_target = paste(path_target, collapse="")

        path_nontarget = c("/var/www/laravel/storage/app/results/",args[6])
        path_nontarget = paste(path_nontarget, collapse="")

        target_list = readLines(path_target)
        nontarget_list = readLines(path_nontarget)
        
    } else {
        target_list = readLines(paste(targetdir, RBP ,'_target_genes.txt', sep = ''))
        nontarget_list = readLines(paste(targetdir, RBP ,'_nontarget_genes.txt', sep = ''))
    }

    existing_target_list = intersect(target_list, row.names(exp_data))
    num_targets = length(existing_target_list)

    existing_nontarget_list = intersect(nontarget_list, row.names(exp_data))
    num_nontargets = length(existing_nontarget_list)
    lncRNA_data = exp_data[lncRNA,skip_col:ncol(exp_data)]

    target_corrs = numeric(length = num_targets)
    i = 1
    for (target in existing_target_list){
        target_data = exp_data[target,skip_col:ncol(exp_data)]
        corr_value = cor(target_data, lncRNA_data, method = 'spearman')
        target_corrs[i] <- round(corr_value, digits = 5)
        i = i + 1
    }
    nontarget_corrs = numeric(length = num_nontargets)
    i = 1
    for (nottarget in existing_nontarget_list){
        nontarget_data = exp_data[nottarget,skip_col:ncol(exp_data)]
        corr_value = cor(nontarget_data, lncRNA_data, method = 'spearman')
        nontarget_corrs[i] <- round(corr_value, digits = 5)
        i = i + 1
    }

    targets_df <- data.frame(group = "Targets" , corr_value = target_corrs)
    nontargets_df <- data.frame(group = "Background" , corr_value = nontarget_corrs)
    plot.data <- rbind(targets_df,nontargets_df)

    plot.data$gene_id = c(existing_target_list, existing_nontarget_list);
    if(length(args) > 4){
        outfilekey = paste(outdir,job_id,  '_',RBP, '_', lncRNA, '_analysis1_', dataset, sep = '')

   }
    else
    {
        outfilekey = paste(outdir, RBP, '_', lncRNA, '_analysis1_', dataset, sep = '')
    }
    write.table(plot.data, file = paste(outfilekey, '.txt', sep = ''),  quote = FALSE, sep = "\t", row.names = FALSE)

    result = wilcox.test(target_corrs, nontarget_corrs)
    pvalue <- result$p.value


    ggplot(plot.data, aes(x=group, y=corr_value, fill= group)) + geom_boxplot() + theme(axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size = rel(1.75)), axis.title.y=element_text(size = rel(1.5), angle = 90)) + labs(y = "Spearman correlation coefficient") +
        annotate("text", x=1.5, y=-0.4, label= paste('p-val: ', signif(pvalue, digits = 3)))



    filename = paste(outfilekey, '.jpeg' , sep= '')
    print(filename)

    ggsave(filename, plot = last_plot())



sink(type="message")
sink(type="output")
}
