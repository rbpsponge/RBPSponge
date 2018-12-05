rm(list=ls(all=TRUE))
library(ggplot2)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
lncRNA = args[1]
RBP = args[2]
job_id = args[3]

if(length(args) > 3){ # user supplied target background file
  targetfile = args[4]
  nontargetfile = args[5]
}

infofileh <- file(paste(job_id, '_analysis3_info.txt', sep = ''), open="wt")
errorfileh <- file(paste(job_id, '_analysis3_error.txt', sep = ''), open="wt")
sink(infofileh, type="output", append=TRUE) # if this is empty, we should display 'unexpected error occurred'
sink(errorfileh, type="message", append=TRUE)

targetdir = '/results/'
datadir = '/datasets/lncRNA_KD_datasets/'
outdir = '/results/'


lncRNA_KD_mapping = data.table(read.table(paste(datadir, "ENSG_id_KD_data_mapping.txt", sep = ''), header = TRUE,  sep = "\t"))

if(!(lncRNA %in%  lncRNA_KD_mapping$gene_id)){

  print("ERROR: There is no knockdown data for this lncRNA")

} else {

  KD_files = lncRNA_KD_mapping[gene_id == lncRNA,]$KD_data_key
  count = 1

  for (LFC_filekey in KD_files){
	  # print(LFC_filekey)
  
    LFC_filename = paste(datadir, LFC_filekey, '_KD_LFCs.txt', sep = '')
    if(file.exists(LFC_filename)) {
         LFC_data = read.table(LFC_filename, header = TRUE, row.names =1, sep = "\t")
    } else {
         print("ERROR: cannot open the knockdown file")
    }
    
    if(length(args) > 4){
      path_target = c("/var/www/laravel/storage/app/results/",args[5])
      path_target = paste(path_target, collapse="")

      path_nontarget = c("/var/www/laravel/storage/app/results/",args[6])
      path_nontarget = paste(path_nontarget, collapse="")
      
      outfilekey = paste(outdir,job_id,  '_',RBP, '_', lncRNA, '_analysis3','_KDdata', count, sep = '')
      target_list = readLines(path_target)
      nontarget_list = readLines(path_nontarget)
     
    } else {
      outfilekey =  paste(outdir,  RBP, '_', lncRNA, '_analysis3','_KDdata', count, sep = '')

      target_list = readLines(paste(targetdir, RBP ,'_target_genes.txt', sep = ''))
      nontarget_list = readLines(paste(targetdir, RBP ,'_nontarget_genes.txt', sep = ''))
    }

    existing_target_list = intersect(target_list, row.names(LFC_data));
    num_targets = length(existing_target_list);

    existing_nontarget_list = intersect(nontarget_list, row.names(LFC_data));
    num_nontargets = length(existing_nontarget_list);

    target_data = LFC_data[existing_target_list, 1]
    nontarget_data = LFC_data[existing_nontarget_list, 1]
    
    if (num_targets > 0 && num_nontargets > 0) {
      
      targets_df <- data.frame(group = "Targets" , value = target_data, gene_id = existing_target_list );
      background_df <- data.frame(group = "Background" , value = nontarget_data, gene_id = existing_nontarget_list);
      save.data <- rbind(targets_df,background_df);
      write.table(save.data, file = paste(outfilekey, '.txt', sep = ''),  quote = FALSE, sep = "\t", row.names = FALSE)
    
      # target_data and nontarget_data can be empty. 
      result = wilcox.test(target_data, nontarget_data)
      pvalue <- result$p.value
      #print(pvalue)
      target_legend = paste('Targets (',length(existing_target_list), ')', sep = '')
      background_legend = paste('Background (',length(existing_nontarget_list), ')', sep = '')
      df <- data.frame(x = c(target_data, nontarget_data), ggg=factor(rep(1:2, c(length(target_data),length(nontarget_data)))))
      ggplot(df, aes(x, colour = ggg)) + ggtitle(LFC_filekey) + scale_size_manual(values=c(3,3)) + stat_ecdf(size=1.5)+ scale_colour_hue(name="", labels=c(paste('Targets (',length(existing_target_list), ')', sep = ''),paste('Background (',length(existing_nontarget_list), ')', sep =''))) +  xlim(-0.5, 0.5) + labs(x = "log fold change") + labs(y = "Cumulative Fraction") + theme(text=element_text(size=16)) + annotate("text", x=0.1, y=-0.02, label= paste('p-val: ', signif(pvalue, digits = 3)))
      

      filename = paste(outfilekey,'.jpeg', sep = '')
      ggsave(filename, plot = last_plot())
      count = count + 1

    } else { # write to log file
      print('ERROR: cannot perform the analysis as the number of target or background genes with expression data is 0.')
    }
  }
  
}

sink(type="message")
sink(type="output")
