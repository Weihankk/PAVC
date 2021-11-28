# **PAVC**
## **Introduction**
PAVC (PAV Classifier) is developed for presence/absence variation (PAV) identification and easily obtain vcf format results. It is based on the results of [SyRI](https://github.com/schneebergerlab/syri), which is a accurate structural variation detect tools. Use PAVC, you can performs accurate classification and get results files in vcf format, which can be very conveniently used for downstream analysis like graph pan-genome construction, GWAS, population genetic analysis, etc.
## **Requirement**
 - [R](https://www.r-project.org/)
 - R packages - [data.table](https://github.com/Rdatatable/data.table) and [Biostrings](https://www.bioconductor.org/packages/release/bioc/html/Biostrings.html)
   - `install.packages("data.table")`
   - ```
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")  
      BiocManager::install("Biostrings")```
 - [bcftools](http://www.htslib.org/download/)

## **Usuage**
**First of all, please download PAVC.**

`git clone https://github.com/Weihankk/PAVC.git`

--------------------------------------------------
1. **Now suppose you have prepared two genomes for SV calling and would like to classify them into PAV.**
   - `refgenome`: Reference genome
   - `qrygenome`: Query genome
2. **SV calling by SyRI.**  
   `nucmer --maxmatch -l 50 -c 100 -t 10 refgenome qrygenome`  
   `delta-filter -m -i 90 -l 100 out.delta > out.filtered.delta`  
   `show-coords -THrd out.filtered.delta > out.filtered.coords`  
   `syri -c out.filtered.coords -d out.filtered.delta --allow-offset 0 -r refgenome -q qrygenome`
   > **syri.out** is raw results of syri, we will use this file for PAV classification.
3. **Run PAV_Classifier**  
   `Rscript GetPAV.R syri.out`  
   `Rscript GetVCF.R pavc.txt query_name`  
   `bgzip pavc.vcf -@ 6`  
   `bcftools index pavc.vcf.gz --threads 6`  
   `bcftools norm -d all -cx -f refgenome pavc.vcf.gz > pavc.norm.vcf` 
   `Rscript FilterPAV.R pavc.norm.vcf 50`
   > For **FilterPAV.R**, _50_ indicated only keep PAV with length >= 50bp.


## **Additional Information**
## Methods
The result of SyRI include two hierarchy: Differences in structure (abbreviated as **DSTR**) & Differences in sequence (abbreviated as **DSEQ**). Note that **DSEQs** are located in **DSTRs**.
- Differences in structure in include five types:
  1. **SYN**, syntenic region
  2. **INV**, inverted region
  3. **TRANS/INVTR**, translocated region or inverted translocated region
  4. **DUP/INVDP**, duplicated region or inverted duplicated region
  5. **NOTAL**, un-aligned region
- Differences in sequence include seven types:
  1. **SNP**, single nucleotide polymorphism
  2. **CPG**, copy gain in query
  3. **CPL**, copy loss in query
  4. **HDR**, highly diverged regions
  5. **TDM**, tandem repeat
  6. **INS**, insertion in query
  7. **DEL**, deletion in query

Therefore, concat all DSTRs we can got the refgenome/qrygenome. However, the coordinates of some DSTRs on the genome are not clear, like NOTAL. 

#### Step 1. Extract DSTRs and sort them according to the refgenome and qrygenome, respectively.
#### Step 2. Classify some DSTRs into presence/absence.
- **SYN** : Neither belongs to presence or absence
- **INV** : Neither belongs to presence or absence
- **TRANS/INVTR** : The part located on refgenome is absence, while the part located on qrygenome is presence.
- **DUP/INVDP** : SyRI classified them into **copygain** and **copyloss**, so we can easily classified them. **DUP/INVDP-copyloss** is absence, while **DUP/INVDP-copygain** is presence.
- **NOTAL** : Refgenome sequence can not align on qrygenome is absence, qrygenome seqeunce can not align on refgenome is presence.
#### Step 3. Detect the breakpoints of these PAVs on refgenome.
#### Step 4. Classify the DSEQs into presence/absence.
- **SNP** : Neither belongs to presence or absence
- **CPG** : Can be regarded as presence
- **CPL** : Can be regarded as absence
- **HDR** : The part located on refgenome is absence, while the part located on qrygenome is presence.
- **TDM** : If refgenome segment length > qry segment length, then this is a absence, otherwise this is a presence
- **INS** : Can be regarded as presnece
- **DEL** : Can be regarded as absence
#### Step 5. Combine PAV


## **Issues**
If you have any questions with installation and usage, please open a new issue in [Issues](https://github.com/Weihankk/PAVC/issues).
