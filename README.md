
#Aug 31st, 2016
---
###Working on MinION data analysis


#August 30th, 2016
---
###Performed Lamda control DNA MinION sequencing using the Rapid sequencing kit (SQK-RAD001)
###TAKEHOMES:
```
- Take photographs of all flowcells before and after sequencing for recording purposes
- The flowcells themselves are bubble traps!
- Try to prevent ever removing liquid from the sensor array on flowcell
- Sequencing did appear to work. Will work on analysis later.
```
<img src=./bioinformatics_notebook_images/20160830_080324.jpg height=410x>

Started earling in the morning to perform the lamda sequencing run in which all MinION customers are to do for general QC and working with the wet lab component for sequencing.

The protocol itself for creating the library is dead simple. Following the CompanION for 'Rapid Lamda control experiment' had no issues. Only really need a thermocycler for a 1min 30C and a 1min 75C temp.

Flowcell FAD23939 had air bubbles across sensor pore array (the part that matters). Attempted to remove by pipetting out small volume of buffer from sample port to no avail. QC run indicated basically no useable pores (!). Ben Schwessinger was present to observe and tried pulling buffer completely off array and re-adding. This did remove the observable bubbles. QC after this still indicated zero active pores. Compare these results from today against initla QC check when flowcell was received:

###FAD23939 QC pore counts to date:
| Pore group | Pore count - Aug 11th | Pore count - Aug 30 (bubbles) | Pore count - Aug 30 (buffer movement, no bubbles) |
|------------|------------|---|---|
| 1          |         361|4 | 0|
| 2          |         317|1 | 0|
| 3          |         210|0 | 0|
| 4          |          74|0 | 0|
| Total      |         962|5 | 0|

So something went seriously wrong with chip after storing it in our fridge (at the desired 2-8C temp). I do not remember there being any bubbles present when package was opened, so perhase they developed over time in fridge?

Here is a picture of the flowcell after Ben re-applied the buffer across the sensor array. Note that there still appears to be a bubble-ish points near the ends of the sensor array:

<img src=./bioinformatics_notebook_images/20160830_134749.jpg width=800x>

Given this mess, I used the other flowcell that was provided (FAD24036) for the lamda control experiment. This flowcell was QC'd for pore counts twice. Once specifically, and a second time during the sequencing run protocol (so pore counts are always determined right before seuqencing). The results are as follows:

###FAD24036 QC pore counts to date:
| Pore group | Pore count - Aug 11th | Pore count - Aug 30 QC | Pore count - Aug 30 in sequencing |
|------------|------------|---|---|
| 1          |         506|505|498|
| 2          |         448|448|439|
| 3          |         308|313|275|
| 4          |          104|111|77|
| Total      |      1366|1377|1289|

So the pore counts appeared just fine and fairly consistant going into the sequencing itself. No bubbles were visable prior to sequencing.

I seledted to have MinKNOW perform local basecalling rather than using Metrichor to call bases in the cloud. Therefore, sequencing commenced using the MinKNOW protocol ```NC_6Hr_Lambda_Burn_In_Run_FLO_MIN104_plus_1D_Basecalling.py```. Everything appeared to work properly, however the tab for local basecalling metrics (quality and read length) never filled in. The overall read size histogram did though. 

<img src=./bioinformatics_notebook_images/lamda_seq_stats.png height=650x>

Will check with Nanopore if there is some issue with MinKNOW currently.

This ran for 6 hours. At the end of 6 hours the program did appear to complete on its own and return me to the MinION status page. However the computer was still quite busy using its processors for something with MinKNOW.

The flowcell was then washed using the flowcell wash kit. After the 6 hour run, there were multiple bubbles now present on the flowcell:
<img src=./bioinformatics_notebook_images/20160830_132454.jpg width=800>

These bubbles remained after washing and returning to fridge:
<img src=./bioinformatics_notebook_images/20160830_133349.jpg width=800x>



#Aug 9th, 2016
---

Getting B stacei annotation files created for Beth

```{r}
getAttributeField <- function (x, field, attrsep = ";") {
     s = strsplit(x, split = attrsep, fixed = TRUE)
     sapply(s, function(atts) {
         a = strsplit(atts, split = "=", fixed = TRUE)
         m = match(field, sapply(a, "[", 1))
         if (!is.na(m)) {
             rv = a[[m]][2]
         }
         else {
             rv = as.character(NA)
         }
         return(rv)
     })
}
##########
gffRead <- function(gffFile, nrows = -1) {
     cat("Reading ", gffFile, ": ", sep="")
     gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
     header=FALSE, comment.char="#", nrows = nrows,
     colClasses=c("character", "character", "character", "integer",  
"integer",
     "character", "character", "character", "character"))
     colnames(gff) = c("seqname", "source", "feature", "start", "end",
             "score", "strand", "frame", "attributes")
     cat("found", nrow(gff), "rows with classes:",
         paste(sapply(gff, class), collapse=", "), "\n")
     stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
     return(gff)
}
```

```{r}
data=gffRead('Bstacei_316_v1.1.gene.gff3.gz')
geneID=getAttributeField(data$attributes, "Name")
data2=cbind(data,geneID)
data2.gene=subset(data2,data$feature=='gene')
out=data2.gene[,c(1,4,5,10,6,7)]

write.table(out,'Bsta.gene.bed',sep='\t',row.names=F,quote=F,col.names=F)
```

#Aug 8th, 2016
---

MinION arrived

attempting configuration test cell

- 19089 MinION Identifcaion
- flowcell 285382121
- sample sre_ctc1
- NC\_CTC_Run

This is the configuration test cell that was present in the MinION upon arrival. Was run through the MinKNOW software. I need to determine if all was successful (I think it was)

Also going through the test data upload of Metrichor to confirm that bass calling will work well.



#Aug 4th, 2016
---

Attempting to run DSS on C24 and Ler CG methylation as I had to remap Ler yesterday.

So you are required to use data smoothing when you do not have biological replicates in DSS. Therefore the output I created was:

- CG DMRS
- p threshold 0.01
- data smoothed

This results in 10224 CG DMRs between C24 and Ler. I have passed them to Ian in case he wants to compare with our tile approach.

#Aug 3rd, 2016
---

####Preparing Ian's samples for DSS analysis

We may need to attempt using DSS (bioconductor.org) to call DMRs across the C24 and Ler parental inbreds for Ian's TCM and TCdM reviews. I am attempting to prepare samples and data formats for this possibility.

There is a specific file type that is wanted for input into DSS. From the documentation:

---
DSS requires data from each BS-seq experiment to be summarized into following information for each CG position:
- chromosome number, 
- genomic coordinate, 
- total number of reads, 
- and number of reads showing methylation. 

For a sample, this information are saved in a simple text file, with each row representing a CpG site. Below shows an example
of a small part of such a file:

| chr   | pos     | N  | X  |
|-------|---------|----|----|
| chr18 | 3014904 | 26 | 2  |
| chr18 | 3031032 | 33 | 12 |
| chr18 | 3031044 | 33 | 13 |

---

So, we can make this data from our bismark cov files

```{bash}
#from bismark alignments output folder...
awk '{print $1"\t"$2"\t"$5+$6"\t"$5}' C24_CHH.bed.bismark.cov > C24_CHH.dss.txt
awk '{print $1"\t"$2"\t"$5+$6"\t"$5}' C24_CHG.bed.bismark.cov > C24_CHG.dss.txt
awk '{print $1"\t"$2"\t"$5+$6"\t"$5}' C24_CpG.bed.bismark.cov > C24_CG.dss.txt

awk '{print $1"\t"$2"\t"$5+$6"\t"$5}' Ler_CHH.bed.bismark.cov > Ler_CHH.dss.txt
awk '{print $1"\t"$2"\t"$5+$6"\t"$5}' Ler_CHG.bed.bismark.cov > Ler_CHG.dss.txt
awk '{print $1"\t"$2"\t"$5+$6"\t"$5}' Ler_CpG.bed.bismark.cov > Ler_CG.dss.txt
```

Throw on 'chr', 'pos', 'N', and 'X' as column names using nano.

Can then start attempting DSS:

```{r}
library(DSS)
library(bsseq)
c24=read.delim('C24_CG.dss.txt',head=T)
ler=read.delim('Ler_CG.dss.txt',head=T)
BSobj = makeBSseqData( list(c24, ler),c('C24','Ler'))
BSobj

dmlTest = DMLtest(BSobj, group1=c("C24"), group2=c("Ler"))
dmlTest.smoothed = DMLtest(BSobj, group1=c("C24"), group2=c("Ler"),smoothing=T)
head(dmlTest)
head(dmlTest.smoothed)

dmrs=callDMR(dmlTest,p.threshold=0.01)
dmrs.smoothed=callDMR(dmlTest.smoothed,p.threshold=0.01)
head(dmrs)
head(dmrs.smoothed)

#showOneDMR(dmrs[1,],BSobj)
write.table(dmrs.smoothed,'C24vsLer_CGDMRS.dss.txt',sep='\t',row.names=F)

```

#July 29th, 2016
---
###Justin Borevitz - 'omics for selection

__Plant climate__

- specificity (selecting zones and seasons for what to plant where)
- sensitivity (diversification under variable climates and environments; require general hardiness)

__'Pre-breeding for Adaptation'__

Geno - Pheno - Enviro triangle

```{}
Better understanding statistical models that can association and predict phenotype from genotype. 

A whole lot of talk, little substance so far
```

__Breeding paradigm__

collect lines
cycle in with phenotyping, selection, and breeding

Requires global phenotypic screens across a large number of conditions

Make predictions across huge swaths of accessions and populations

```
So many buzz words
```

DivSeek ```http://www.divseek.org/```

```
It's all good ideas and concepts, but requires serious manpower and $$$$ to actually perform in any meaningfull way
```

###Steven Swain - CSIRO

data61

Digiscape FSP

FSP Environomics Collaboration-Hub

CSIRO is global impact and revenue focused

transforming nitrogen fixation into crops and plants

Topical RNAi application (a la Monsanto Bioactive)

Association expression level with phenotypes across multiple environments


#July 19th, 2016
---

Tim Stuart has provided new files to attempt the trans TE-DMR analysis for his manuscript. Previously, both TEPID insertions and deletions were combined into a 'TE_poly' file containing both types of elements. However, further analysis by Tim has indicated that it may be worthwhile to look at these two types of variants seperatly. Beyond this, there was an issue with the methylation DMR files indicating no coverage as a value of 0, which could seriously screw things up.

Therefore, I am revisiting the TE-DMR trans association analysis and performing it four times for these datasets:

- C DMRs (13484) vs TE insertions (15077)
- C DMRs (13484) vs TE deletions (5856)
- CG DMRs (40268) vs TE insertions (15077)
- CG DMRs (40268) vs TE deletions (5856)

Tim has provided new files for methylation ```c_dmr_allC.tsv``` and ```cg_dmrs_allC.tsv``` as well as new dataframes of TE insertions ```TE_insertions_matrix.tsv``` and deletions ```TE_deletions_matrix.tsv```

All code and results for each analysis are found in their respective folder and ```TIM_TRANS_ANALYSIS_JULY2016.html``` file.

To begin, dataset accessions are confirmed identical resulting in 136 accessions. We trim this down to 124 accessions that were identified in the TE-to-SNP LD analysis to eliminate some of the populations structure found across these accessions.

TE variants are subset to only those 'common' variants with a minor-allele-frequence > 3% (4 accessions). This substantially cuts down on the number of TE variants that we can test:

| TE set     | raw variants | >3% MAF | proportion kept |
|------------|--------------|---------|-----------------|
| insertions | 15077        | 2782    | 0.185           |
| deletions  | 5856         | 2859    | 0.488           |

We can also calculate the FDR for each of these association tests by counting the number of associations above our 1% permuted data threshold as so:

(count above real - count above perm) / (count above real)

Which gives us....not the best FDR in the world:

| FDR %  | individual r2 values | sum r2 values | individual binary states | sum binary states |
|--------|----------------------|---------------|--------------------------|-------------------|
| C_ins  | 14.38                | 87.77         | 14.38                    | 77.77             |
| C_del  | 17.17                | 87.65         | 17.16                    | 81.29             |
| CG_ins | 7.55                 | 52.54         | 7.55                     | 21.62             |
| CG_del | 10.27                | 47.27         | 10.26                    | 25                |

Even so, we can get a count of our punative trans-associated TE variants for each contrast:

| possible trans Tes | count | max DMRs it hits | max %DMRs it hits |
|--------------------|-------|------------------|-------------------|
| C_ins              | 126   | 526              | 3.901             |
| C_del              | 155   | 457              | 3.389             |
| CG_ins             | 37    | 2021             | 5.019             |
| CG_del             | 40    | 2395             | 5.948             |

So if I'm honest, I think we don't really have any evidence of trans bands. In fact, if we look at the plotted heatmaps it is much harder to identify the cis-band going up the diagonals of the C DMR set.

If anything, we see lots of horizontal bands, which may indicate that we would need to do some heavy filtering of the Schmitz DMR calls in order to try and clear things up. I'm not sure how to best approach that however.




#June 17th, 2016
____


There are some revisions to be done in regards to the Brachypodium distachyon reference methylomes manuscript. I am trying to organize my thoughts and scripts to have a clearer picture of exactly what is being done and where any new analyses will fit in:

In the end what I, and everyone else, should be doing is identifying the time in which a manuscript will be created and directly begin organizational steps in order to prepare for publication and reproducibility.

There are key steps to making this a reality

- __Start by uploading any created sequencing data to a public repo EARLY ON.__ By doing so you have completed one of the most annoying steps of any genomic paper prep. A specific benifit of doing this early on is that you can develop your analysis scripts to start directly from pulling sequence data (e.g. from the SRA). If you are smart you can create scripts that can flag off this step for subsequent analyses when you have the sequencing data in hand (the same holds true for any time-consuming steps such as read alignment).

- __Keep an up-to-date tab on any third-party data that is used within your analyses.__ Everythign we do often builds on other labs data, or other public annotation information that is associated with your experimental system. This includes reference genome(s), gene annotations, and other public annotations which may be used within your analysis. By knowing what _version_ as well as the public _location_ of such data (TAIR, Phytozome, genome connsortium, etc) you will make your own life, and others reproducing your results, easier.

- __Timestamp everything.__ Analysis often beings in a 'dirty phase' of analysis with rapid testing and code/plots/data tables that are good, bad, and ugly. Often you begin building off of these for downstream steps. You often never know that a specific step or output will make it out of this phase. Devise a clear timestamp method for any output of any script. This allows you to keep track of when certain files were made, and forces a clear naming convention for your own sanity.

- __Try to get as close to final publication figures from first output.__ From my experience, the effort you put into a 'final' figure is often left behind when the 'final' figure is no longer actually final. There are often changes, edits, tweaks that may require a new primary figure or plot to be developed which is not always easy to slot into the edited manuscript figure. I have often spent a fair amount of time in Illustrator getting things just right only to have a _plot.v3USETHISONE.pdf_ get spit out right at the end requiring an additional round of finalization for publication. This is inevitable as far as I can tell. Therefore, taking the time to make your plotting software (ggplot is my go to) get as close as possible to what you want in the end. This includes colors, legends, labels, font sizes, etc. The more you can code into your scripts, the easier it will be to deal with any last minute changes.

- __When you have completed all of this, make sure the first steps of the script highlight the required software, libraries, and disk space needed to make things work.__ If you make it this far, you may just have a script which will generate everything you made. Well done! However, take the time to clearly identify the specific software and libraries required for your scripts to work. Best practices would also include the specific versions of all items that you used as they seem to change and break things quite often.


My personal goal would be to develop a series of scripts in which with the scripts alone it would:

- Download all primary data
- Download all annotation third-part data required
- Perform all base processing of data
- Perform all biological analyses performed in the paper
- Output semi-finished versions of all primary and supplemental figures and tables (noting that they may be cleaned further in Illustrator for publication)

This is much easier said than done. I have yet to personally meet this goal. I imagine that this level of organization requires one to accept that you will get so far with first-pass analyses and developing your results, then have to _start from the beginning_ with this organization in mind. Things should move faster at this point because you know where you are headed (rather than in discovery phase).

I have also been continuing to try and find the best method for pulling scripts together to include code for all aspects of the required software. I often move from bash scripts to R, back to bash, and sometimes others which would require a more unified method. Right now I often develop R scripts with commented out bash code which I manually execute while walking through the script. I don't like doing this and either would any person trying to re-develop your results and analyses.

I imagine that beyond this would be a magical world of VMs that someone could spin up to do all of this in the cloud, however I'm still focused on these baby steps.


#May 9th, 2016
----

Trying to understand _Probability Theory - The Logic of Science_ by E.T. Jaynes

###Probability defined for:
- Frequency & Fair odds 
(dice rolls, gambling)

vs 

- credible, sound judgements
(finding expectations from incomplete information)
- Kolmogorov
- 
###Statistics
- from the state for mathematical analysis of affairs
- defining measurement error
- sampling variation

Laplace

- Linear regression
- Distribution of errors 

Inference (Fisher & Pearson)

- How do we infer the true 'state of being'?
- p-values
- confidence intervals
- hypothesis testing

Bayesian statistics

- How do prob estimates get updated given new data?
- conditional probabilities

Given all of this, what exactly is probability in the real world?

###Freqentist

- counting, _measurements_ of a physical prop, sampling, 
- This only really makes sense when we can make the measurements
vs

###Subjectionist 

__'A degree of belief'__ - something we all have

- Given we cannot simulate worlds to test outcomes, what do we do when we cannot measure directly
- It is naturally uncomfortable to include our subjective degrees of belief
- What if we have different degrees of belief when we can properly test (roll the dice)?
- Kolmogorov axioms don't necessesarily fit well here

Is probability forwarding looking and statistics backward looking?

in this book, Jaynes frames probability as an extention of LOGIC specifically.
We are reasoning with incomplete information. We can get Kolmogorov axioms as a derivitive of logic. From this, the separation of probability vs statistics dissapears. Also destroys frequentist ad hoceries and methods as a fundamental misunderstanding of probability.

In doing so, Jaynes shows that a subjectionist view allows for more powerful statements.

Hypothesis to Data (test hypothesis to explain our data) is not good. We want __Data to hypothesis__ (the thing in the real world)

###Resolves about 400 years of prob theory by extending logic!

----

#April 18th, 2016
----

Created folder structure on local machine for IGV sessions for Joanne's methylc-seq data visualizations.

bis snp snp-calling from bisulfite sequencing data for Diep's work.

----
#April 14th, 2016
----
The alignments of Joanne's Methyl-seq data has been completed on edmund. The summary information is provided below:

| Input fastq | Sample | ref genome | method | bismark ver | samtools ver | total reads | flt reads | % unique aligned | unique map | multi-map | no-map | %CG | %CHG | %CHH |
|---------------------|------------------------|---------------|--------------------------------|----|---------|-----------------|----------|----------|------|----------|---------|---------|------|-----|-----|
| Col_mock_C_R1.fastq.gz | Col_mock_3DAI | ../../genomes/TAIR10/assembly/ | se | v0.13.0 | 1.1-26-g29b0367 | 32705946 | 32652653 | 78.1 | 25507247 | 1958063 | 5187343 | 23.2 | 7.5 | 2.2 |
| Col_Fox_D_R1.fastq.gz  | Col_Fox_3DPI  | ../../genomes/TAIR10/assembly/ | se | v0.13.0 | 1.1-26-g29b0367 | 29501104 | 29460765 | 67.6 | 19925070 | 5170754 | 4364941 | 21.8 | 7.0 | 2.1 |
| rdd_Fox_B_R1.fastq.gz  | rdd_Fox_3DPI  | ../../genomes/TAIR10/assembly/ | se | v0.13.0 | 1.1-26-g29b0367 | 29645007 | 29612191 | 67.7 | 20055251 | 5865639 | 3691301 | 23.0 | 7.0 | 2.0 |
| rdd_mock_A_R1.fastq.gz | rdd_mock_3DPI | ../../genomes/TAIR10/assembly/ | se | v0.13.0 | 1.1-26-g29b0367 | 36957919 | 36915555 | 80.1 | 29574380 | 1810380 | 5530795 | 23.6 | 7.3 | 2.1 |

Another important note is calculating the bisulfite conversion rate by looking for methylated cytosines in the chloroplast (which is unmethylated). This can be done via the output files and a bash one-liner:

~~~bash
grep "ChrC" *CHH.bed.bismark.cov | awk '{ met+= $5} { unmet += $6} { total = met + unmet } END {print 100-((met / total)*100)}'
~~~

This results in the folowing conversion rates (higher the better):

| Sample | Conversion Rate (%) |
|--------|-----------------|
| Col_mock_3DPI | 99.6488 |
| Col_Fox_3DPI | 99.6395 |
| rdd_mock_3DPI | 99.6367 |
| rdd_Fox_3DPI | 99.6164 |

So overall, it looks like conversion rate will not be a problem with this dataset. Although one can correct for it, it would make minimal difference with these rates.

We can also look at global methylation levels to determine if there are any clear differences across samples.

Start by pulling together base-pair resolution data for the four samples into a single file in R:

~~~r

context=c('CpG','CHG','CHH')
for(i in 1:3){
  out=NULL
  files=dir(pattern=paste('*',context[i],'.bed.bismark.cov',sep=''))
  input=read.delim(files[1],head=F)
  colnames(input)[4:6]=paste(files[1],names(input)[4:6],sep='')
  
  for(q in 2:length(files)){
    add=read.delim(files[q],head=F)
    colnames(add)[4:6]=paste(files[q],names(add)[4:6],sep='')
    input=merge(input,add,by=c('V1','V2','V3'))
  }
 write.table(input,paste(context[i],'_allsamples.cov',sep='\t',row.names=F,quote=F)
}
~~~

Ran 100bp DMR caller for CG, 60% difference, 10coverage and at least 3 sites in the window.

Also did CHG 40%; CHH 10%

Have to do some serious filtering

----
#April 7th, 2016
----

Today I am attempting to make an annotation file of TAIR10 genes that have a flanking TE. This is to better define genes that fall on a euchromatin / heterochromatin boundary within the genome. I will also add information regarding orientation of chromatin in relation to the geen (5prime / 3prime / strand).

My first attempts have used the `TAIR10_gene.bed` and the `TAIR10_transposable_element.bed` files, however I kept getting a Segmentation Fault when running on _edmund_. So, I will attempt it on my own laptop to see if results change:

I still get segmentation faults after 12,957 lines when running:

~~~bash
closestBed -iu -D ref -a TAIN10_transposable_elements.bed -b TAIR10_gene.sorted.bed
~~~

Perhaps bedops will get the job done by giving the closest upstream _and_ downstream elements. I can further require that it does not report overlaps via ```--no-overlaps``` and also report back the distance for each ``--dist``. So, I could take all TAIR10 genes and map them against a file that contains all TAIR10 genes as well as all transposable elements:

~~~bash
closest-features --no-overlaps --dist --delim '\t' TAIR10_gene.sorted.bed TAIR10_gene_te.sorted.bed > test.output
~~~

This reports out a bed-like file that contains:

- initial query sequence
- upstream hit (which would be itself if not for --no-overlaps)
- distance to upstream feature
- downstream hit
- distance to downstream feature

These sections are all seperated by tabs as specified in the ``--delim '\t' `` flag

Need to make sure that the TE and gene files that are combined contain the same number of columns. Added two columns of NA to the gene file to match the number found in the TE file.

~~~r
gene=read.delim('TAIR10_gene.sorted.bed',head=F)
out=cbind(gene[,1:7],rep(NA,nrow(gene)),rep(NA,nrow(gene)))
write.table(out,'TAIR10_gene.sorted.forcomb.bed',sep='\t',row.names=F,col.names=F,quote=F)

~~~

~~~bash
cat TAIR10_transposable_element.bed TAIR10_gene.sorted.forcomb.bed > TAIR10_gene_te.bed
sort-bed TAIR10_gene_te.bed > TAIR10_gene_te.sorted.bed
closest-features --no-overlaps --dist --delim '\t' TAIR10_gene.sorted.bed TAIR10_gene_te.sorted.bed > test.output
~~~

We can now go look into this bedops output to see the results:

~~~r
data = read.delim('test.output',head=F,sep='\t')

#genes with nothing upstream are writted out with NA 
#however it does not match full formatting. 
#Therefore, we can pull these columns seperatly.

no.upstream=subset(data,is.na(data$V9)==T)
both=subset(data,is.na(data$V9)==F)
no.downstream=subset(both,is.na(both$V19)==T)
both=subset(both,is.na(both$V19)==F)
class=paste(both$V13,both$V23,sep='-')
both=cbind(both,class)
out=subset(both,both$class=='transposable_element-transposable_element')
oneside=subset(both,both$class=='transposable_element-gene' | both$class=='gene-transposable_element')
write.table(out,'TAIR10_tete_gene.bed',sep='\t',row.names=F,quote=F,col.names=F)
write.table(oneside,'TAIR10_teoneside_gene.bed',sep='\t',row.names=F,quote=F,col.names=F)
~~~

Of the __27,404__ genes with calls on both sides, just over half (56%) are flanked by other genes:
![barplot of gene neighbor classes](./bioinformatics_notebook_images/p1.png)

So I think that the file contains only protein coding genes. No TE genes, no pseudogenes, tRNAs, etc.

From this, we see that although overal gene density decreases near the centromers (TE-rich regions), that is largely where we find genes with a TE flanking both sides:
![barplot of gene neighbor classes](./bioinformatics_notebook_images/p2.png)

A local view to highlight it:
![barplot of gene neighbor classes](./bioinformatics_notebook_images/p3.png)

The output file has been creaed as ```TAIR10_tete_gene.bed``` containing 3134 genes along with a file of genes with TEs flanking one side in ```TAIR10_teoneside_gene.bed``` containing 8762 genes.

----

#April 6th, 2016
----

Today I have attempted to align Joanne's bisulfite sequencing data of 3DPI Col-0 and _rdd_ mutants for both mock and control samples. This data has been passed to myself and Peter Crisp, along with sRNA sequencing data (1/3/6 DPI) and mRNA-seq (3 DPI) data.

The fastq files have been passed via my processing script as single-end sequencing reads:

```
/home/steve/scripts/wgbs_pipelinev0.4.sh -se <input fastq> ../../../genomes/TAIR10/assembly/ <outname>
```

