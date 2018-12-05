rm(list=ls(all=TRUE))
library(ggplot2)


args = commandArgs(trailingOnly=TRUE)
gene_id_input = args[1]
gene_id = gene_id_input
dataset = args[2]
job_id = args[3]

infofileh <- file(paste(job_id, '_analysis0_info.txt', sep = ''), open="wt")
errorfileh <- file(paste(job_id, '_analysis0_error.txt', sep = ''), open="wt")
sink(infofileh, type="output", append=TRUE) # if this is empty, we should display 'unexpected error occurred'
sink(errorfileh, type="message", append=TRUE)

datadir = '/datasets/'
outdir = '/results/'

present = TRUE
    
if(dataset == 'EMTAB2706'){
  load(paste(datadir, 'EMTAB2706_mapping_frame.RData', sep = ''))
  if(!(gene_id %in%  mapping_frame$ENSG_id) && !(gene_id %in%  mapping_frame$gene_symbol)){
    print('ERROR: this gene do not exist in this expression data')
    present = FALSE
  }else{
    if(!(gene_id %in%  mapping_frame$ENSG_id) && (gene_id %in%  mapping_frame$gene_symbol) )
      gene_id = toString(mapping_frame[mapping_frame$gene_symbol == gene_id,]$ENSG_id)
    skip_col = 1
    filename = paste(datadir, 'EMTAB2706_median_tpm.RData', sep = '')
    if (file.exists(filename)){
      load(filename)
    }else{
      print("ERROR: cannot open input data file")
    }
  }
} else if(dataset == 'EMTAB2770'){
  load(paste(datadir, 'EMTAB2770_mapping_frame.RData', sep = ''))
  if(!(gene_id %in%  mapping_frame$ENSG_id) && !(gene_id %in%  mapping_frame$gene_symbol)){
    print('ERROR: this gene do not exist in this expression data')
    present = FALSE
  }else{        
    if(!(gene_id %in%  mapping_frame$ENSG_id) && (gene_id %in%  mapping_frame$gene_symbol) )
      gene_id = toString(mapping_frame[mapping_frame$gene_symbol == gene_id,]$ENSG_id)
    skip_col = 1
    filename = paste(datadir, 'EMTAB2770_median_tpm.RData', sep = '')
    if (file.exists(filename)){
      load(filename)
    }else{
      print("ERROR: cannot open input data file")
    }
  }
} else if(dataset == 'GTEX'){
  load(paste(datadir, 'GTEX_mapping_frame.RData', sep = ''))
  if(!(gene_id %in%  mapping_frame$ENSG_id) && !(gene_id %in%  mapping_frame$gene_symbol)){
    print('ERROR: this gene do not exist in this expression data')
    present = FALSE
  }else{ 
    if(!(gene_id %in%  mapping_frame$ENSG_id) && (gene_id %in%  mapping_frame$gene_symbol) )
      gene_id = toString(mapping_frame[mapping_frame$gene_symbol == gene_id,]$ENSG_id)      
    skip_col = 1
    filename = paste(datadir, 'GTEX_median_tpm.RData', sep = '')
    if (file.exists(filename)){
      load(filename)
    }else{
      print("ERROR: cannot open input data file")
    }
  }
} else if (startsWith(dataset, 'User'))# user uploaded data, dataset should refer to the name of the file
{
  # this will already have ensg ids  
  
  skip_col = 1
  filename = paste(datadir, dataset, '_median.RData', sep = '')
  if (file.exists(filename)){
    load(filename)
    if(!(gene_id %in%  rownames(exp_data))){
      present = FALSE
    }
    
  }else{
    print("ERROR: cannot open input data file")
  }      
}
if (present)
{
    outfilekey = paste(outdir, gene_id_input, '_analysis0_', dataset, sep = '')
    #write.csv(args[1], "SalesData.csv", row.names = FALSE)

    data = exp_data[gene_id, skip_col:ncol(exp_data)]
    sorted_data = sort(names(data), index.return = TRUE)
    data = data[sorted_data$ix]
    xlim_max = as.numeric(max(data) + 5)
    if((dataset == 'EMTAB2706' || dataset == "GTEX") || ncol(exp_data) < 120  ){
        df <- data.frame(tissues = names(data), expression = data) #coord_flip() +
        ggplot(df, aes(tissues, expression)) + geom_col() +  theme(text=element_text(size=16)) + coord_flip() #theme(axis.text.x=element_text(angle=90, hjust=1))
        filename = paste(outfilekey, '_.jpeg' , sep= '')
        ggsave(filename, plot = last_plot())
    } else if(dataset == 'EMTAB2770' || ncol(exp_data) >= 120){
        print("here")
        df <- data.frame(tissues = names(data[1:70]), expression = data[1:70]) #coord_flip() +
        ggplot(df, aes(tissues, expression)) + geom_col() +  theme(text=element_text(size=8)) + coord_flip() #theme(axis.text.x=element_text(angle=90, hjust=1))
        filename = paste(outfilekey, '_part1.jpeg' , sep= '')
        ggsave(filename, plot = last_plot())
        print("here 2")
        df <- data.frame(tissues = names(data[71:length(data)]), expression = data[71:length(data)]) #coord_flip() +
        ggplot(df, aes(tissues, expression)) + geom_col() +  theme(text=element_text(size=8)) + coord_flip()
        filename = paste(outfilekey, '_part2.jpeg' , sep= '')
        ggsave(filename, plot = last_plot())

    }

sink(type="message")
sink(type="output")
}
