#!/usr/bin/env Rscript
# rm(list=ls())

#for DEBUGGING
# setwd("D:/Dropbox/_git/progression/")
# #FOR DEBUGGING
# options<-NULL
# ###INPUT OPTIONS
# filename<-"sample1vdm_q200_subtracted.vcf";vcf_format<-"GATK"
# filename<-"MMD_ref955.vcf";vcf_format<-"MiModD"
# filename<-"nor22.vcf";vcf_format<-"Freebayes"
# #input file
# options$inf<-filename
# #vcf format
# options$vcff<-vcf_format

# ###BASE OPTIONS
# #interval type
# interval_type<-"C.elegans"
# options$itype<-interval_type
# #quality filter
# read_qual<-200
# options$qual<-read_qual
# #for variants to exclude, bottom threshold for which to be used, i.e. 0=ALL, 1=HOM only, 0.8=near HOM)
# freq_threshold<-"0-1"
# options$freqthr<-freq_threshold
# #type of allelic ratio (AB/ratio) only applies to FB
# allele_ratio<-"ratio"
# options$allr<-allele_ratio
# #include only variants of type SNP or all
# snp_only<-FALSE
# options$snp<-snp_only

# ###PLOT OPTIONS
# #loess span
# loess_span<-0.4
# options$lsp<-loess_span
# #point colour
# plot_color<-"black"
# options$pcol<-plot_color
# #loess plot colour
# loess_color<-"red"
# options$lcol<-loess_color
# xaxis_standard<-TRUE
# #standardize x-axis interval (e.g. 1Mb interval)
# options$xstand<-xaxis_standard
# #bin size for barplot
# bin_size<-1000000
# options$bsize<-bin_size
# #normalization for barplot
# bfreq_norm<-FALSE
# options$bnorm<-bfreq_norm

# ###VARIANT EXCLUSION OPTIONS
# #files with variants to exclude
# # exclusion_list$exclf<-c("nor22chunk1.txt","nor22chunk5.txt")
# exclusion_list<-"FALSE"
# options$exclf<-exclusion_list
# excl_threshold<-"0-1"
# options$excl_thr<-excl_threshold
# excl_loess_color<-"green"
# #additional colour option for pre-subtraction line
# options$exclcol<-excl_loess_color

# ###OUTPUT OPTIONS
# #fixed file names
# vcfoutput_filename<-"vcf_filted_output.txt"
# options$outn<-vcfoutput_filename
# vcfparsed_filename<-"vcf_parsed_output.txt"
# varbin_filename<-"varlist_binned.txt"
# options$vbinn<-varbin_filename
# options$outpn<-vcfparsed_filename
# pdf_filename<-"vdm_mapping_plot.pdf"
# options$pdfn<-pdf_filename

#if TRUE, kept simple with only the relevant table and plot output. if FALSE, additional output tables generated.
#only parsing the INFO column of the QUAL filtered variant list greatly speeds up runs
# simplemode<-FALSE
# options$smode<-simplemode

library("getopt")
#input from trailing line arguments
args <- commandArgs(trailingOnly = TRUE)
#read the options from input commandArgs
option_specification = matrix(c(
  'inf','in',1,'character',
  'itype','i',2,'character',
  'vcff','vf',2,'character',  
  'qual','q',2,'double', 
  'allr','a',2,'character',
  'snp','s',2,'logical',
  'lsp','l',2,'double',
  'pcol','pc',2,'character',
  'lcol','lc',2,'character',
  'xstand','x',2,'logical',
  'bsize','b',2,'integer',
  'freqthr','ft',2,'character',
  'bnorm','n',2,'logical',
  'exclf','ef',2,'character',
  'exclcol','ec',2,'character',
  'outn','o',2,'character',
  'outpn','op',2,'character',
  'vbinn','vb',2,'character',
  'pdfn','p',2,'character',
  'smode','sm',2,'logical'
), byrow=TRUE, ncol=4)
options = getopt(option_specification)


myfunction<-function(inf,itype,vcff,qual,allr,snp,lsp,pcol,lcol,xstand,bsize,freqthr,bnorm,exclf,exclcol,outn,outpn,vbinn,pdfn,smode){
  #input file
  filename<-inf
  #choose species for chrom config
  interval_type<-itype
  #vcf format
  vcf_format<-vcff
  #quality filter
  read_qual<-as.numeric(qual)
  #allele ratio from AB or AO/(AO+RO)
  allele_ratio<-allr
  #snp only or include all variant types
  snp_only<-snp
  #loess span
  loess_span<-as.numeric(lsp)
  #scatterplot point colour
  plot_color<-pcol
  #loess curve colour
  loess_color<-lcol
  #uniform y-axis scaling
  xaxis_standard<-xstand
  #binsize
  bin_size<-as.numeric(bsize)
  #freq thresholds input as range ###-###, e.g. "0.2-1.0"
  freq_threshold<-freqthr
  #allele ratio greater or equal to this value considered hom ALT
  freq_threshold<-freqthr
  threshold_upper<-as.numeric(gsub(".*-","",freq_threshold))
  #allele ratio less than or equal to this value considered hom REF
  threshold_lower<-as.numeric(gsub("-.*","",freq_threshold))
  #normalize frequency barplot
  bfreq_norm<-bnorm
  #filenames for variant subtraction
  exclusion_list<-exclf
  #pre-subtraction loess curve colour
  excl_loess_color<-exclcol
 

 #filename for ouput table
  vcfoutput_filename<-outn
  #filename for ouput table with parsed full variant list
  vcfparsed_filename<-outpn
  #filename for output table with binned variant list
  binvar_filename<-vbinn
  #filename for pdf
  pdf_filename<-pdfn

  simplemode<-smode

  ###ADDITIONAL OPTIONAL PARAMETERS
  #only use variants with allele ratio in this range for subtraction (dependent on AB/ratio option)
  # exclusion_threshold<-exclthr
	exclusion_threshold<-"0-1"
  excl_threshold_upper<-as.numeric(gsub(".*-","",exclusion_threshold))
  excl_threshold_lower<-as.numeric(gsub("-.*","",exclusion_threshold))
   
  ###FIXED PARAMETERS
  #chromosome intervals in Mb rather than custom
  interval_unit<-1000000 
  #linkage scatter plot yaxis max value=1
  sp_yaxis<-1 

#informative files names overrides any specified in arguments
autoname<-TRUE
if(autoname==TRUE){
	if(options$exclf=="FALSE"){exclname<-"nosubtract"}else{exclname<-options$exclf}
	# vcfoutput_filename<-paste(gsub("vcf","",options$inf),"_output_q",options$qual,"-",paste(options$exclf,sep="",collapse=""),".txt",sep="")
	vcfoutput_filename<-paste(gsub("\\.vcf","",options$inf),"_output_q",options$qual,"-",paste(exclname,sep="",collapse=""),".txt",sep="")
	vcfparsed_filename<-paste(gsub("\\.vcf","",options$inf),"_parsed_q",options$qual,"-",paste(exclname,sep="",collapse=""),".txt",sep="")
	pdf_filename<-paste(gsub("\\.vcf","",options$inf),"_plot_q",options$qual,"-",paste(exclname,sep="",collapse=""),".pdf",sep="")
	varbin_filename<-paste(gsub("\\.vcf","",options$inf),"_varlist_binned",options$qual,"-",paste(exclname,sep="",collapse=""),".txt",sep="")
}
  
  #transparency for selected colour (to see plot points underneath)
  plot_color<-rgb(c(col2rgb(plot_color)[1]),c(col2rgb(plot_color)[2]),c(col2rgb(plot_color)[3]),max=255,alpha=150)
  loess_color<-rgb(c(col2rgb(loess_color)[1]),c(col2rgb(loess_color)[2]),c(col2rgb(loess_color)[3]),max=255,alpha=150)
  excl_loess_color<-rgb(c(col2rgb(excl_loess_color)[1]),c(col2rgb(excl_loess_color)[2]),c(col2rgb(excl_loess_color)[3]),max=255,alpha=150)
  
  
######################
###READ IN VCF FILE
#extract column names
vcf_readin<-readLines(filename)   
#find header line, i.e. last line to begin with #
for(l in 1:length(vcf_readin)){
  vcf_readinl<-vcf_readin[l]
  if(substr(vcf_readinl,1,1)=="#"){next}
  else if(substr(vcf_readinl,1,1)!="#"){rowline<-l-1;break}
}
vcf_header<-vcf_readin[rowline]
#e.g. CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\trgSM"
vcf_header<-gsub("#","",vcf_header)
vcf_colnames<-unlist(strsplit(vcf_header,"\t"))

#extract data (hashed vcf header skipped with read.table)
vcf_rtable<-read.table(filename,sep="\t",stringsAsFactors=FALSE)
names(vcf_rtable)<-vcf_colnames

if(simplemode=="TRUE"){
  #to speed up runtime- quality filter here by skipping parsing of full vcf (parsed output table will be incomplete) 
  vcf_rtable<-subset(vcf_rtable,QUAL>=read_qual)
}


######################
###PREPARE DATA
#list all INFO columns
vcfinfo_splitall<-strsplit(vcf_rtable$INFO,split=";")	
vcfinfo_colnall<-gsub("=.*","",unlist(vcfinfo_splitall))
vcfinfo_allcol<-unique(vcfinfo_colnall)

#parse data
vcfinfo_dat<-NULL
multiallele_counter<-0
diviserror_counter<-0
for(i in c(1:nrow(vcf_rtable))){
  vcf_line<-vcf_rtable[i,]
  #remove chrom or chr prefix from chromosome value
  if(grepl("chrom",vcf_line$CHROM,ignore.case=TRUE)==TRUE){
    vcf_line$CHROM<-gsub("chrom","",vcf_line$CHROM,ignore.case=TRUE)
  }else if(grepl("chr",vcf_line$CHROM,ignore.case=TRUE)==TRUE){
    vcf_line$CHROM<-gsub("chr","",vcf_line$CHROM,ignore.case=TRUE)
  }
#PARSE INFO
  vcfinfo_split<-strsplit(vcf_line$INFO,split=";")	
  vcfinfo_coln<-gsub("=.*","",unlist(vcfinfo_split))
  vcfinfo_cold<-gsub(".*=","",unlist(vcfinfo_split))
  vcfinfo_ldat<-data.frame(t(vcfinfo_cold),stringsAsFactors=FALSE)
  names(vcfinfo_ldat)<-vcfinfo_coln
   
	#ensure all columns are populated (blanks filled with "NA")
	for(allcoli in vcfinfo_allcol){
		if((allcoli %in% names(vcfinfo_ldat))==TRUE){
			assign(allcoli,vcfinfo_ldat[,names(vcfinfo_ldat)==allcoli])

			vcfinfo_ldat<-vcfinfo_ldat[,!names(vcfinfo_ldat) %in% allcoli]
			vcfinfo_coln<-c(names(vcfinfo_ldat),allcoli)
			vcfinfo_ldat<-cbind(vcfinfo_ldat,get(allcoli))
			names(vcfinfo_ldat)<-vcfinfo_coln
			}else{
			assign(allcoli,"NA")
			vcfinfo_coln<-c(names(vcfinfo_ldat),allcoli)
			vcfinfo_ldat<-cbind(vcfinfo_ldat,get(allcoli))
			names(vcfinfo_ldat)<-vcfinfo_coln
		}
	}

	if(vcf_format=="MiModD"){
	  #skip if commas in values to avoid returning errors
	  if(grepl(",",vcfinfo_ldat$AC)==TRUE){
		 multiallele_counter<-multiallele_counter+1
		 next
	  }
		AO<-as.numeric(vcfinfo_ldat$AC)
		RO<-as.numeric(vcfinfo_ldat$DP)-as.numeric(vcfinfo_ldat$AC)
		vcfinfo_ldat<-cbind(AO,RO,vcfinfo_ldat)
	}
	
	if(vcf_format=="GATK"){
		#ref+alt count information in rgSM column
		if(("rgSM" %in% names(vcf_line))==TRUE){
		vcf_line_rg<-vcf_line$rgSM
		}else{
		rgcol_pos<-grep("FORMAT",names(vcf_line))+1
		vcf_line_rg<-vcf_line[,rgcol_pos]
		}
		###alt method of identifying rg column (ie may be assigned as sample name)	
		vcfrgsm_split<-strsplit(vcf_line_rg,split=":")
		refalt_pos<-grep("[[:digit:]]/[[:digit:]]",unlist(vcfrgsm_split))+1
		refalt_val<-unlist(vcfrgsm_split)[refalt_pos]
		#ref value
		RO<-gsub(",.*","",refalt_val)
		#alt value
		AO<-gsub(".*,","",refalt_val)
		DP<-unlist(vcfrgsm_split)[refalt_pos+1]
		vcfinfo_ldat<-cbind(AO,RO,DP,vcfinfo_ldat)
	}  
	
  #skip if commas in values to avoid returning errors
  if(grepl(",",vcfinfo_ldat$AO)==TRUE){
    multiallele_counter<-multiallele_counter+1
    next
  }
  #skip divide by zero errors (under "ratio" setting for ratio calculation)
  if(as.numeric(vcfinfo_ldat$AO)+as.numeric(vcfinfo_ldat$RO)=="0"){
    diviserror_counter<-diviserror_counter+1
    next
  }		 
  
  #PARSE ANNOTATION
  #check if annotation column present
	if("ANN" %in% names(vcfinfo_ldat)=="TRUE"){
		  ann_rparsed<-unlist(strsplit(vcfinfo_ldat$ANN[1],split="\\|"))[1:20]
		  ann_rparsed[ann_rparsed==""]<-"novalue"
		  ann_parsed<-data.frame(t(ann_rparsed),stringsAsFactors=FALSE)
		  names(ann_parsed)<-paste("ANN",c(1:dim(ann_parsed)[2]),sep="")
		  #remove duplicate redundant INFO column (fully parsed)
		  vcf_line<-vcf_line[,names(vcf_line)!="INFO"]
		  
		  #dataset keeping parsed annotation (partial parsed-for relevant)
		  vcfinfo_ldat<-vcfinfo_ldat[,names(vcfinfo_ldat)!="ANN"]
		  #append copy of original data to parsed data
		  vcfinfo_lldat<-cbind(vcf_line,vcfinfo_ldat,ann_parsed)
		  vcfinfo_dat<-rbind(vcfinfo_dat,vcfinfo_lldat)
	}else{
		#remove duplicate redundant INFO column (fully parsed)
		vcf_line<-vcf_line[,names(vcf_line)!="INFO"]
		#append copy of original data to parsed data
		vcfinfo_lldat<-cbind(vcf_line,vcfinfo_ldat)
		vcfinfo_dat<-rbind(vcfinfo_dat,vcfinfo_lldat)
	}
}

print(paste("rows with multiple alleles skipped: ",multiallele_counter,sep=""))
# print(paste("rows with AO+RO=0 (not multiple alleles) skipped: ",diviserror_counter,sep=""))

#ENSURE CORRECT DATATYPES
#convert dataframe columns of factor type to character type
vcfinfo_dat<-data.frame(lapply(vcfinfo_dat,as.character),stringsAsFactors=FALSE)
#convert to numeric if column is numeric and not string
for(n in c(1:dim(vcfinfo_dat)[2])){
  #suppress warnings when columns with strings encountered (not converted)
  suppressWarnings(
    colnum_index<-!is.na(as.numeric(vcfinfo_dat[,n]))
  )
  if(all(colnum_index)==TRUE){
    vcfinfo_dat[,n]<-as.numeric(vcfinfo_dat[,n])
  }
}


######################
#RATIO CALCULATION
#ratio calculation from AC and DP
if(vcf_format=="Freebayes"){
	if(allele_ratio=="AB"){
		adj_ratio<-replace(vcfinfo_dat$AB,vcfinfo_dat$AB==0,1)
		RATIO<-adj_ratio
		}else if(allele_ratio=="ratio"){
			RATIO<-c(vcfinfo_dat$AO/(vcfinfo_dat$AO+vcfinfo_dat$RO))
		}
	}else if(vcf_format=="GATK"){
			RATIO<-c(vcfinfo_dat$AO/(vcfinfo_dat$AO+vcfinfo_dat$RO))
	}else if(vcf_format=="MiModD"){
		#ratio calculation from AC and DP
		RATIO<-c(vcfinfo_dat$AC/vcfinfo_dat$DP)
		vcfinfo_dat<-cbind(vcfinfo_dat,RATIO)
		#reorder to ensure chrom and pos as first columns
		vcfinfo_dat<-vcfinfo_dat[with(vcfinfo_dat,order(CHROM,POS)),]
		vcfinfo_pdat<-vcfinfo_dat
		vcfinfo_dat<-subset(vcfinfo_dat,QUAL>=read_qual)
	}
	
	vcfinfo_dat<-cbind(vcfinfo_dat,RATIO)
	vcfinfo_dat<-vcfinfo_dat[with(vcfinfo_dat,order(CHROM,POS)),]
	#keep a copy of unfiltered before restricting by qual
	vcfinfo_pdat<-vcfinfo_dat
	vcfinfo_dat<-subset(vcfinfo_dat,QUAL>=read_qual)
	vcfinfo_dat<-vcfinfo_dat[with(vcfinfo_dat,order(CHROM,POS)),]

# vcfinfo_dat<-vcfinfo_pdat

#CONSIDER ONLY SNP VARIANTS
if(snp_only==TRUE){
	if(vcf_format=="Freebayes"){
		vcfinfo_dat<-subset(vcfinfo_dat,CIGAR=="1X")
	}else{
		snpindex<-NULL
		for(i in 1:nrow(vcfinfo_dat)){
			vcfinfo_dati<-vcfinfo_dat[i,]
			if(nchar(vcfinfo_dati$REF)==1 && nchar(vcfinfo_dati$ALT)==1){
				snpindex<-c(snpindex,i)
			}
		}	
	vcfinfo_dat<-vcfinfo_dat[snpindex,]
	}
}


######################
#SUBTRACT VARIANTS FROM EXCLUSION LIST
if(exclusion_list!="FALSE"){
  #keep copy of pre-subtraction data for later plotting
  vcfinfo_origdat<-vcfinfo_dat
  
  #identifiers for exclusion list based on CHROM/POS/REF/ALT
  index1<-paste(vcfinfo_dat$CHROM,vcfinfo_dat$POS,sep="_")
  index2<-paste(vcfinfo_dat$CHROM,vcfinfo_dat$POS,vcfinfo_dat$REF,sep="_")
  index3<-paste(vcfinfo_dat$CHROM,vcfinfo_dat$POS,vcfinfo_dat$REF,vcfinfo_dat$ALT,sep="_")
  vcfinfo_dat<-cbind(vcfinfo_dat,index1,index2,index3)
  
  print(paste("before subtraction: ",nrow(vcfinfo_dat),sep=""))
  #loop and subtract through exclusion lists (if multiple files)
  for(exclusion_ind in exclusion_list){
    exclin<-read.table(exclusion_ind,header=TRUE)
  
  #THRESHOLD FILTER ON EXCLUSION LIST VARIANTS  
    if(allele_ratio=="AB"){
      exclin<-subset(exclin,adj_AB>=excl_threshold_lower)
      exclin<-subset(exclin,adj_AB<=excl_threshold_upper)
    }
    if(allele_ratio=="ratio"){
      exclin<-subset(exclin,ratio>=excl_threshold_lower)
      exclin<-subset(exclin,ratio>=excl_threshold_upper)
    }
  
    #identifiers for vcf data based on CHROM/POS/REF/ALT  
    index1<-paste(exclin$CHROM,exclin$POS,sep="_")
    index2<-paste(exclin$CHROM,exclin$POS,exclin$REF,sep="_")
    index3<-paste(exclin$CHROM,exclin$POS,exclin$REF,exclin$ALT,sep="_")
    exclin<-cbind(exclin,index1,index2,index3)
    #exclude based on CHROM+POS+REF+ALT
    vcfinfo_dat<-subset(vcfinfo_dat,!(index3 %in% exclin$index3))
  }
  print(paste("after subtraction: ",nrow(vcfinfo_dat),sep=""))
}


######################
#WRITE TO OUTPUT
#select relevant columns 2 variant type; 4 gene; !5 wormbase ID; 8 type change; 10 nucleotide change; 11 amino acid change; 16 warning message
#check if annotation present in input vcf
if("ANN" %in% names(vcfinfo_dat)=="TRUE"){
	vcfinfo_simp<-subset(vcfinfo_dat,select=c("CHROM","POS","QUAL","DP","REF","ALT","AO","RO","RATIO","ANN2","ANN4","ANN8","ANN10","ANN11","ANN16"))
	names(vcfinfo_simp)<-c("CHROM","POS","QUAL","DP","REF","ALT","AO","RO","RATIO","VARTYPE","GENE","TYPE","NTCHANGE","PRCHANGE","WARNINGS")
	vcfsimp_dat<-vcfinfo_simp[with(vcfinfo_simp,order(CHROM,POS)),]
}else{
	vcfinfo_simp<-subset(vcfinfo_dat,select=c("CHROM","POS","QUAL","DP","REF","ALT","AO","RO","RATIO"))
	names(vcfinfo_simp)<-c("CHROM","POS","QUAL","DP","REF","ALT","AO","RO","RATIO")
	vcfsimp_dat<-vcfinfo_simp[with(vcfinfo_simp,order(CHROM,POS)),]
}
#write table with quality filtered variants for VDM plotting and relevant columns
try(
  write.table(vcfsimp_dat,vcfoutput_filename,sep="\t",quote=FALSE,row.names=FALSE)
  ,silent=TRUE)

#write table with all unfiltered variants and all columns including parsed INFO

if(simplemode!="TRUE"){
	try(
	write.table(vcfinfo_pdat,vcfparsed_filename,sep="\t",quote=FALSE,row.names=FALSE)
	,silent=TRUE)
}


######################
###CHROMOSOME (INTERVAL) ARRANGEMENT
#define chromosome and chromosome size in Mb
if(interval_type == 'C.elegans'){
  chrom_n<-c('I','II','III','IV','V','X')
  chrom_mb<-c(16,16,14,18,21,18)
  interval_frame<-data.frame(chrom_n,chrom_mb)
} else if(interval_type == 'Zebrafish'){
  chrom_n<-c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25')
  chrom_mb<-c(61,61,64,63,76,60,78,57,59,47,47,51,55,54,48,59,54,50,51,56,45,43,47,44,39)
  interval_frame<-data.frame(chrom_n,chrom_mb)
} else if(interval_type == 'Brachypodium'){
  chrom_n<-c('1','2','3','4','5')
  chrom_mb<-c(75,60,60,50,30)
  interval_frame<-data.frame(chrom_n,chrom_mb)
} else if(interval_type == 'Arabidopsis'){
  chrom_n<-c('1','2','3','4','5')
  chrom_mb<-c(31, 20,24,19,27 )
  interval_frame<-data.frame(chrom_n,chrom_mb)
} else{
  #user interval file- no headers, with chromosome in column 1 (format CHR# or CHROM#) and size in Mb (rounded up) in column 2
  user_interval_type<-read.table(interval_type)
  if(grepl("chrom",user_interval_type[1,1],ignore.case=TRUE)==TRUE){
    user_interval_type[,1]<-gsub("chrom","",user_interval_type[,1],ignore.case=TRUE)
  }else if(grep("chr",user_interval_type[1,1],ignore.case=TRUE)==TRUE){
    user_interval_type[,1]<-gsub("chr","",user_interval_type[,1],ignore.case=TRUE)
  }
  chrom_n<-user_interval_type[,1]
  chrom_mb<-user_interval_type[,2]
  interval_frame<-data.frame(chrom_n,chrom_mb)
}
names(interval_frame)<-c("CHROM","INTERVAL")


######################
###PLOTTING
#VDM SCATTER PLOT
#save to pdf
pdf(file=pdf_filename,width=9,height=8)
#par(mfrow=c(2,3))
for(chromind in interval_frame$CHROM){
  #subset by data by chromosome for plotting
  intervalind<-interval_frame$INTERVAL[interval_frame$CHROM==chromind]
  chr_dat<-subset(vcfsimp_dat,CHROM==chromind,silent=TRUE)
  
  #same subsetting by chromosome for pre-subtraction data
  if(exclusion_list!="FALSE"){
    chr_origdat<-subset(vcfinfo_origdat,CHROM==chromind,silent=TRUE)
  }
  
  #define x-axis upper limit
  if(xaxis_standard==TRUE){
    #for standardized x-axis (max x-axis chromosome length)
    scupper_xaxis<-max(interval_frame$INTERVAL)
    scupper_xval<-scupper_xaxis*interval_unit
  } else if(xaxis_standard==FALSE){
    scupper_xaxis<-intervalind
    scupper_xval<-intervalind*interval_unit
  }
	 plot(chr_dat$POS,chr_dat$RATIO,cex=0.60,xlim=c(0,scupper_xval),ylim=c(0,sp_yaxis),main=paste("Chr",chromind," Variant Discovery Mapping",sep=""),xlab="Position along Chromosome (in Mb)",ylab='Allelic Frequency',pch=10, col=plot_color,xaxt='n')				
	 try(lines(loess.smooth(chr_dat$POS,chr_dat$RATIO,span=loess_span),lwd=5,col=loess_color))
	 #plot loess curve for data without subtraction of exclusion variants
	 if(exclusion_list!="FALSE"){
		try(lines(loess.smooth(chr_origdat$POS,chr_origdat$RATIO,span=loess_span),lty="longdash",lwd=4,col=excl_loess_color))
	 }
	 axis(1,at=seq(0,scupper_xval,by=interval_unit),labels=c(0:scupper_xaxis))     		  
	 abline(h=seq(0,sp_yaxis,by=0.1),v=c(1:scupper_xaxis)*interval_unit,col="gray")

}


######################
#graph barplots
location_index<-NULL
meanSNP_dat<-NULL
# prepare table of counts and calculations
for(chromind in interval_frame$CHROM){
  #for standardized x-axis
  if(xaxis_standard==TRUE){
    intervalind<-max(interval_frame$INTERVAL)*1000000/bin_size
  } else if(xaxis_standard==FALSE){
    intervalind<-interval_frame$INTERVAL[interval_frame$CHROM==chromind]*1000000/bin_size
  }
  #start intervals with **1 and end with **0
  interval_begin<-c(((0:(intervalind-1))*bin_size)+1)
  interval_end<-c((1:intervalind)*bin_size)
  
  #define x-axis upper limit
  if(xaxis_standard==TRUE){
    upper_xaxis<-max(interval_frame$INTERVAL)
  } else if(xaxis_standard==FALSE){
    upper_xaxis<-interval_frame$INTERVAL[interval_frame$CHROM==chromind]
  }
  #prepare columns
  snp_counter<-0
  purealt_counter<-0
  pureref_counter<-0
  het_counter<-0
  chr_mean<-0
  normed_freq<-0
  ratio<-0
  
  interval_index<-data.frame(chromind,interval_begin,interval_end,snp_counter,purealt_counter,pureref_counter,het_counter,chr_mean,normed_freq)
  chr_dat<-subset(vcfinfo_dat,CHROM==chromind)
  #ratio calculation
  ratio<-chr_dat$AO/(chr_dat$AO+chr_dat$RO)
  
  #if counter based on adj_AB or ratio
	chr_purealtdat<-subset(chr_dat,ratio>=threshold_upper)
	chr_purerefdat<-subset(chr_dat,ratio<=threshold_lower)
	chr_hetdat<-subset(chr_dat,ratio>threshold_lower & ratio<threshold_upper)

  #if chromosome with data, count number of snps within each bin (positions rounded up to nearest bin), else skip to next chromosome
  if(dim(chr_dat)[1]>0){
    for(i in 1:dim(chr_dat)[1]){
      chr_datind<-chr_dat[i,]
      #round up to nearest bin-size interval			
      chr_datind_upper<-ceiling(chr_datind$POS/bin_size)*bin_size
      interval_coln<-NULL;interval_rown<-NULL
      #identify row and and counter column to increment
      interval_coln<-which(names(interval_index)=="snp_counter")
      interval_rown<-match(chr_datind_upper,interval_index$interval_end)		
      interval_index[interval_rown,interval_coln]<-c(interval_index$snp_counter[interval_rown]+1)
    }
  }else{
    next
  }
  #if chromosome with pure AO, count number of snps with each bin (positions rounded up to nearest bin)
  if(dim(chr_purealtdat)[1]>0){
    for(i in 1:dim(chr_purealtdat)[1]){
      chr_purealtind<-chr_purealtdat[i,]
      chr_purealtind_upper<-ceiling(chr_purealtind$POS/bin_size)*bin_size
      interval_coln<-NULL;interval_rown<-NULL
      interval_coln<-which(names(interval_index)=="purealt_counter")
      interval_rown<-match(chr_purealtind_upper,interval_index$interval_end)		
      interval_index[interval_rown,interval_coln]<-c(interval_index$purealt_counter[interval_rown]+1)
    }
  }
  #if chromosome with pure RO, count number of snps with each bin (positions rounded up to nearest bin)
  if(dim(chr_purerefdat)[1]>0){
    for(i in 1:dim(chr_purerefdat)[1]){
      chr_purerefind<-chr_purerefdat[i,]
      chr_purerefind_upper<-ceiling(chr_purerefind$POS/bin_size)*bin_size
      interval_coln<-NULL;interval_rown<-NULL
      interval_coln<-which(names(interval_index)=="pureref_counter")
      interval_rown<-match(chr_purerefind_upper,interval_index$interval_end)			
      interval_index[interval_rown,interval_coln]<-c(interval_index$pureref_counter[interval_rown]+1)
    }
  }
  #if chromosome with hets, count number of snps with each bin (positions rounded up to nearest bin)
  if(dim(chr_hetdat)[1]>0){
    for(i in 1:dim(chr_hetdat)[1]){
      chr_hetind<-chr_hetdat[i,]
      chr_hetind_upper<-ceiling(chr_hetind$POS/bin_size)*bin_size
      interval_coln<-NULL;interval_rown<-NULL
      interval_coln<-which(names(interval_index)=="het_counter")
      interval_rown<-match(chr_hetind_upper,interval_index$interval_end)			
      interval_index[interval_rown,interval_coln]<-c(interval_index$het_counter[interval_rown]+1)
    }
  }
  #irrespective of standardized x-axis, mean should be calculated from actual interval range of chromosome
  chr_mean<-sum(interval_index$purealt_counter)/(interval_frame$INTERVAL[interval_frame$CHROM==chromind]*1000000/bin_size)
  interval_index$chr_mean<-chr_mean
  meanSNP_dat<-rbind(meanSNP_dat,data.frame(chromind,chr_mean))
  #normalization treatment for if SNPs are AO=0, AO=total SNPs, or AO and RO in bin
  for(i in 1:dim(interval_index)[1]){
    chr_intind<-interval_index[i,]
    #hom definition based on specified upper and lower thresholds
    if(chr_intind$purealt_counter<=threshold_lower){
      interval_coln<-NULL
      interval_coln<-which(names(interval_index)=="normed_freq")
      interval_index[i,interval_coln]=0
    } else if (chr_intind$purealt_counter==chr_intind$snp_counter){
      interval_coln<-NULL
      interval_coln<-which(names(interval_index)=="normed_freq")
      interval_index[i,interval_coln]=(chr_intind$purealt_counter)^2/chr_mean
    } else {
      interval_coln<-NULL
      interval_coln<-which(names(interval_index)=="normed_freq")
      interval_index[i,interval_coln]=chr_mean*(chr_intind$purealt_counter)^2/(chr_intind$snp_counter-chr_intind$purealt_counter)
    }
  }
  location_index<-rbind(location_index,interval_index)
}

###WRITE TO OUTPUT
if(smode==FALSE){
#write table with quality filtered variants for VDM plotting and relevant columns
try(
  write.table(location_index,binvar_filename,sep="\t",quote=FALSE,row.names=FALSE)
  ,silent=TRUE)
}

for(chromind in interval_frame$CHROM){
	interval_index<-location_index[location_index$chromind==chromind,]
	#assign 0 values to avoid empty datatable error
	if(dim(interval_index)[1]==0){
		interval_index[1,]<-rep(0,dim(interval_index)[2])		
	}
	#set up x_axis 
	if(xaxis_standard==TRUE){
		#for standardized x-axis (max x-axis chromosome length)
		bpupper_xaxis<-max(interval_frame$INTERVAL)
		bpupper_xval<-bpupper_xaxis*interval_unit
	}else if(xaxis_standard==FALSE){
		bpupper_xaxis<-intervalind
		bpupper_xval<-intervalind*interval_unit
	}
	#set up y_axis range for barplots
	if(bfreq_norm==TRUE){
		bp_yaxis<-5*ceiling(max(location_index$normed_freq)/5)
	#assign non-0 value to yaxis to avoid error
		if(bp_yaxis==0){
			bp_yaxis<-10
		}
		if(xaxis_standard==TRUE){
			bplot<-barplot(interval_index$normed_freq,space=0,ylim=c(0,bp_yaxis),main=paste("Chr",chromind," Variant Only",sep=""),xlab="Position along Chromosome (in Mb)",ylab='Normalised Frequency')	
		}else if(xaxis_standard==FALSE){
			bplot<-barplot(interval_index$normed_freq,space=0,xlim=c(0,bpupper_xaxis),ylim=c(0,bp_yaxis),main=paste("Chr",chromind," Variant Only",sep=""),xlab="Position along Chromosome (in Mb)",ylab='Normalised Frequency')	
		}
		
	}else if(bfreq_norm==FALSE){
		bp_yaxis<-5*ceiling(max(location_index$purealt_counter)/5)
		#assign non-0 value to yaxis to avoid error
		if(bp_yaxis==0){
		bp_yaxis<-10
		}
		if(xaxis_standard==TRUE){
			bplot<-barplot(interval_index$purealt_counter,space=0,ylim=c(0,bp_yaxis),main=paste("Chr",chromind," Variant Only",sep=""),xlab="Position along Chromosome (in Mb)",ylab='Frequency')	
		}else if(xaxis_standard==FALSE){
			bplot<-barplot(interval_index$purealt_counter,space=0,xlim=c(0,bpupper_xaxis),ylim=c(0,bp_yaxis),main=paste("Chr",chromind," Variant Only",sep=""),xlab="Position along Chromosome (in Mb)",ylab='Frequency')	
		}
	}
	bp_xaxis1<-as.numeric(bplot)
	bp_xaxis2<-c(bp_xaxis1,tail(bp_xaxis1,1)+bp_xaxis1[2]-bp_xaxis1[1])
	bp_xaxis<-bp_xaxis2-bp_xaxis1[1]

	axis(1,at=bp_xaxis,labels=seq(0,bpupper_xaxis,by=c(bin_size/1000000)))
	}
	dev.off()
}

myfunction(inf=options$inf,itype=options$itype,vcff=options$vcff,qual=options$qual,allr=options$allr,snp=options$snp,lsp=options$lsp,pcol=options$pcol,lcol=options$lcol,
           xstand=options$xstand,bsize=options$bsize,bnorm=options$bnorm,freqthr=options$freqthr,exclf=options$exclf,exclcol=options$exclcol,outn=options$outn,pdfn=options$pdfn,outpn=options$outpn,vbinn=options$vbinn,smode=options$smode)

#Rscript vdm_tool.R --inf "nor22.vcf" --vcff "Freebayes" --itype "C.elegans" --qual 200 --allr "ratio" --snp TRUE --lsp 0.4 --pcol "black" --lcol "red" --xstand TRUE --bsize 1000000 --freqthr "0-1" --bnorm FALSE --exclf "FALSE" --exclcol "green" --outn "test-output.txt" --pdfn "test-plot.pdf" --outpn "test-parsed.txt" --sm "TRUE"
#Rscript vdm_tool.R --inf "iz1905gatk.vcf" --vcff "GATK" --itype "C.elegans" --qual 200 --allr "ratio" --snp TRUE --lsp 0.4 --pcol "black" --lcol "red" --xstand TRUE --bsize 1000000 --freqthr "0-1" --bnorm FALSE --exclf "FALSE" --exclcol "green" --outn "test-output.txt" --pdfn "test-plot.pdf" --outpn "test-parsed.txt" --sm "TRUE"
#Rscript vdm_tool.R --inf "MMD_ref955.vcf" --vcff "MiModD" --itype "C.elegans" --qual 200 --allr "ratio" --snp TRUE --lsp 0.4 --pcol "black" --lcol "red" --xstand TRUE --bsize 1000000 --freqthr "0-1" --bnorm FALSE --exclf "FALSE" --exclcol "green" --outn "test-output.txt" --pdfn "test-plot.pdf" --outpn "test-parsed.txt" --sm "TRUE"
