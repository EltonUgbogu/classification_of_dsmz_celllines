clear
ls
cd colDataAnnData/
ls
cd logs/
ls
cat pan_cancer_umap_27421300.err 
cat pan_cancer_umap_27421300.out 
lcear
l
clear
ls
cat pan_cancer_umap_27426895.err 
cat pan_cancer_umap_27426895.out 
squeue -u chu25
cat pan_cancer_umap_27426895.out 
cat pan_cancer_umap_27426895.err 
ls
rm pan_cancer_umap_274*
clear
ls
cat pan_cancer_umap_27428321.err 
cat pan_cancer_umap_27428321.out 
squeue -u chu25
cat pan_cancer_umap_27428321.out 
cat pan_cancer_umap_27428321.err 
clear
ls
rm pan_cancer_umap_27428321.*
ls
cat pan_cancer_umap_27428322.err 
cat pan_cancer_umap_27428322.out 
ls
rm pan_cancer_umap_27428322.*
clear
ls
cat pan_cancer_umap_27429748.err 
cat pan_cancer_umap_27429748.out 
cat pan_cancer_umap_27429748.err 
cat pan_cancer_umap_27429748.out 
clear
cat pan_cancer_umap_27429748.err 
cat pan_cancer_umap_27429748.out 
cat pan_cancer_umap_27429748.err 
cat pan_cancer_umap_27429748.out 
clear
cat pan_cancer_umap_27429748.out 
clear
cat pan_cancer_umap_27429748.out 
clear
cat pan_cancer_umap_27429748.out 
cat pan_cancer_umap_27429748.err 
cat pan_cancer_umap_27429748.out 
cat pan_cancer_umap_27429748.err 
cat pan_cancer_umap_27429748.out 
clear
cat pan_cancer_umap_27429748.out 
clear
cat pan_cancer_umap_27429748.out 
clear
cat pan_cancer_umap_27429748.out 
cat pan_cancer_umap_27429748.err 
cat pan_cancer_umap_27429748.out 
clear
cat pan_cancer_umap_27429748.out 
cat pan_cancer_umap_27429748.err 
clear
cd ..
ls
cd TCGA/
ls
cd tcga
ls
cd ,,
cd ..
cd  ..
cd data/
clear
ls
cd icgc_mmml/
ls
cd logs/
ls
cd ..
cd tcga/
ls
cd brca_processed/
ls
cd ..
cd tcga_dlbc/
ls
tree
cd ..
cd .
cd ..
cd tcga/brca_processed/
clear
ls
cd ..
cd TCGA/
ls
cd logs/
ls
ls -lh
cat tumour_brca_purity_analysis_tumour_brca_purity_27764978.err 
cat tumour_brca_purity_analysis_tumour_brca_purity_27764978.out 
squeue -u chu25
cat tumour_brca_purity_analysis_tumour_brca_purity_27764978.out 
cat tumour_brca_purity_analysis_tumour_brca_purity_27764978.err 
squeue -u chu25
cat tumour_brca_purity_analysis_tumour_brca_purity_27764978.err 
cat tumour_brca_purity_analysis_tumour_brca_purity_27764978.out 
clear
ls
rm tumour_laml_purity_analysis_tumour_laml_purity_27766902.*
rm tumour_brca_purity_analysis_tumour_brca_purity_27766903.*
cd ..
cd data/
ls
cd tcga/
ls
cd tcga_brca/
ls
cd ~/dsmz/leu_lym_dsmz_analysis/
ls
cd preprocessing/
ls
cat io.
cat io.R
clear
cd ~/data/
ls
cd hema_aggregated/
ls
cdcd tumour_purity_results/
ls
cd tumour_purity_results/
ls
pwd
cd ..
clear
ls
cd mpn/
ls
cd GSE228758/
ls
cd results/
ls
cd aggregated/
l
sls
ls
cd tumour_purity/
ls
pwd
clear
conda env list
conda activate tcga-r.env
conda activate tcga-r-env
clear
Rscript /home/chu25/TCGA/hema_pipeline/hema_pipeline/scripts/hema_tumour_purity_estimation.R   --expr   /home/chu25/data/cll/GSE237907/results/aggregated/raw_count_matrix.rds   --map    /home/chu25/data/cll/GSE237907/results/aggregated/ensembl_to_hgnc.tsv   --outdir /home/chu25/data/cll/GSE237907/results/aggregated/tumour_purity   --threshold 0.7   --meta   /home/chu25/data/cll/GSE237907/results/aggregated/sample_metadata.csv
conda activate tcga-r-env
conda env list
conda activate haema-r-env
clear
Rscript /home/chu25/TCGA/hema_pipeline/hema_pipeline/scripts/hema_tumour_purity_estimation.R   --expr   /home/chu25/data/cll/GSE237907/results/aggregated/raw_count_matrix.rds   --map    /home/chu25/data/cll/GSE237907/results/aggregated/ensembl_to_hgnc.tsv   --outdir /home/chu25/data/cll/GSE237907/results/aggregated/tumour_purity   --threshold 0.7   --meta   /home/chu25/data/cll/GSE237907/results/aggregated/sample_metadata.csv
conda install tidyestimate
clear
conda env list
conda activate tidyestimate
clear
Rscript /home/chu25/TCGA/hema_pipeline/hema_pipeline/scripts/hema_tumour_purity_estimation.R   --expr   /home/chu25/data/cll/GSE237907/results/aggregated/raw_count_matrix.rds   --map    /home/chu25/data/cll/GSE237907/results/aggregated/ensembl_to_hgnc.tsv   --outdir /home/chu25/data/cll/GSE237907/results/aggregated/tumour_purity   --threshold 0.7   --meta   /home/chu25/data/cll/GSE237907/results/aggregated/sample_metadata.csv
Rscript /home/chu25/TCGA/hema_pipeline/hema_pipeline/scripts/hema_tumour_purity_estimation.R   --expr      /home/chu25/data/mpn/GSE228758/results/aggregated/raw_count_matrix.rds   --map       /home/chu25/data/mpn/GSE228758/results/aggregated/ensembl_to_hgnc.tsv   --outdir    /home/chu25/data/mpn/GSE228758/results/aggregated/tumour_purity   --threshold 0.7   --meta      /home/chu25/data/mpn/GSE228758/results/aggregated/sample_metadata.csv
clear
Rscript /home/chu25/TCGA/hema_pipeline/hema_pipeline/scripts/hema_tumour_purity_estimation.R   --expr      /home/chu25/data/bl/GSE199413/results/aggregated/raw_count_matrix.rds   --map       /home/chu25/data/bl/GSE199413/results/aggregated/ensembl_to_hgnc.tsv   --outdir    /home/chu25/data/bl/GSE199413/results/aggregated/tumour_purity   --threshold 0.7   --meta      /home/chu25/data/bl/GSE199413/results/aggregated/sample_metadata.csv
tail -n 120 /home/chu25/data/bl/GSE199413/results/aggregated/tumour_purity/hema_purity_pipeline_20260216_163142.log
tail -n 120 /home/chu25/data/bl/GSE199413/results/aggregated/tumour_purity/hema_purity_pipeline_20260216_163142.logclear
clear
tail -n 120 /home/chu25/data/bl/GSE199413/results/aggregated/tumour_purity/hema_purity_pipeline_20260216_163142.log
clear
tail -n 120 /home/chu25/data/bl/GSE199413/results/aggregated/tumour_purity/hema_purity_pipeline_20260216_163142.log
ls -lh /home/chu25/data/bl/GSE199413/results/aggregated/
find /home/chu25/data/bl/GSE199413/results/aggregated -maxdepth 2 -type f -name "*.rds"
clear
Rscript /home/chu25/TCGA/hema_pipeline/hema_pipeline/scripts/hema_tumour_purity_estimation.R   --expr      /home/chu25/data/pri_bl/GSE199413/results/aggregated/raw_count_matrix.rds   --map       /home/chu25/data/pri_bl/GSE199413/results/aggregated/ensembl_to_hgnc.tsv   --outdir    /home/chu25/data/pri_bl/GSE199413/results/aggregated/tumour_purity   --threshold 0.7   --meta      /home/chu25/data/pri_bl/GSE199413/results/aggregated/sample_metadata.csv
clear
cd ~/TCGA/
ls
cd tcga
cd ..
ls
cd tcga_dlbc/
ls
cat tumour_purity_estimation.
cat tumour_purity_estimation.sh 
clear
ls
sbatch tumour_purity_estimation.sh 
cd logs/
ls
cat dlbc_purity_analysis_28760435.out 
cat dlbc_purity_analysis_28760435.err 
cat dlbc_purity_analysis_28760435.out 
ls -lh
cat dlbc_purity_analysis_28760435.out 
cat dlbc_purity_analysis_28760435.err 
squeue -u chu25
cat dlbc_purity_analysis_28760435.err 
cat dlbc_purity_analysis_28760435.out 
clear
ls
cat dlbc_purity_analysis_28760435.out 
clear
cat dlbc_purity_analysis_28760435.err 
cat dlbc_purity_analysis_28760435.out 
clear
cd /home/chu25/data/tcga/tcga_dlbc/tumour_purity_results0.7
ls
cd ..
ls
cd ..
ls
cd ..
ls
cd tcga/
ls
cd ..
clear
cd ..
l
ls
cd TCGA/
ls
cd tcga
ls
cd ..
cd hema_pipeline/
ls
cd hema_pipeline/
ls
cd logs/
ls
clear
cd ..
ls
cd logs/
ls
cd ..
ls
cd ..
ls
cd ..
ls
cd ..
ls
cd dsmz/
cd ..
clear
ls
cd data/
l
ls
cd tcga/
cd ..
cd hema_aggregated/
ls
cd tumour_purity_results/
ls
clear
cd ..
ls
cd tcga/
ls
ls -lh
clear
cd ~/colDataAnnData/
ls
cd logs/
ls
cat pan_cancer_umap_27429748.out 
clear
ls
cd colDataAnnData/
ls
cd results/
ls
c pan_cancer_umap_main/
ls
cd pan_cancer_umap_main/
ls
ls -lh
cd ..
ls
mkdir pan_cancer_umapp_archieve2
mv pan_cancer_umap pan_cancer_umapp_archieve2/
mv pan_cancer_umap_main/ pan_cancer_umapp_archieve2/
claer
ls
clear
ls
cd ../scripts/
ls
sbatch pan_cancer_main.sh 
sbatch pan_cancer_umap.sh 
cd ..
ls
cd dsmz/
clear
ls
cd leu_lym_dsmz_analysis/
ls
cd results/
ls
cd ..
cd preprocessing/
ls
cd logs/
ls
cat leu_lym_dsmz_io_25735347.err 
cat leu_lym_dsmz_io_25735347.out 
clear
cd ..
cd ~
ls
cd TCGA/
clear
ls
cd icgc_mmml_bl_dlbcl_fl/
ls
cd logs/
ls
cat download_mmml_rnaseq_25389100.err 
cat download_mmml_rnaseq_25389100.out 
clear
cd ..
ls
cd ..
ls
cd LUSC_LUAD/
ls
cd ..
ls
cd tcga
ls
cd ..
ls
clear
cd ..
ls
cd TCGA/
ls
cd all_tcga_count_data/
ls
cd tumour_purity/
ls
cd results/
ls
cd ..
ls
sbatch tumour_brca_purity.sh 
ls
sbatch tumour_tcga_laml_purity.sh 
sbatch tumour_brca_purity.sh 
clear
ls
cd results/
ls
cd ..
ls
cd ..
ls
cd tumour_purity/
ls
cd results/
ls
cd ..
ls
cd logs/
ls
clear
ls
cd ..
ls
cd tumour_purity/
ls
cd results/
ls
pwd
cd ..
ls
sbatch tumour_brca_purity.sh 
sbatch tumour_tcga_laml_purity.sh 
sbatch tumour_brca_purity.sh 
squeue -u chu25
sbatch tumour_brca_purity.sh 
sbatch tumour_tcga_laml_purity.sh 
sbatch tumour_brca_purity.sh 
cd ~/data/
clear
ls
cd cll/
ls
tree
clear
ls
cd ..
ls
cd tcga/
ls
cd ..
cd hema_aggregated/
ls
cd tumour_purity_results/
ls
pwd
cd ..
ls
clear
cd ~/TCGA/hema_pipeline/
ls
cd hema_pipeline/
ls
cd scripts/
ls
cat hema_tumour_purity_estimation.R 
clear
ls
cd ..
ls
cd ..
ls
cd nbl
ls
cd nbl_scripts/
ls
cat fastq_2_countdata.R
clear
ls
cat fastq_2_countdata.R
cd ../../dsmz
clear
cd ..
cd ...
cd ..
cd dsmz/
clear
ls
cd fastqc_preprocessing/
ls
cd fastqc/
ls
cd .
cd ..
clear
ls
cat Snakefile 
clear
ls
cd config/
ls
car config.yaml
cart config.yaml
clear
cat config.yaml 
cd ~/TCGA/logs/
clear
ls
cat tumour_brca_purity_analysis_tumour_brca_purity_27781050.out 
cat tumour_brca_purity_analysis_tumour_brca_purity_27781050.err 
cd /home/chu25/TCGA/all_tcga_count_data/tumour_purity/results/TCGA_BRCA/TCGA-BRCA/
ls
clear
cd =/home/chu25/TCGA/logs
cd /home/chu25/TCGA/logs
ls
cat tumour_brca_purity_analysis_tumour_brca_purity_28756529.err 
cat tumour_brca_purity_analysis_tumour_brca_purity_28756529.out 
squeue -u chu25
cat tumour_brca_purity_analysis_tumour_brca_purity_28756529.out 
squeue -u chu25
cat tumour_brca_purity_analysis_tumour_brca_purity_28756529.out 
cat tumour_brca_purity_analysis_tumour_brca_purity_28756529.err 
cat /home/chu25/TCGA/all_tcga_count_data/tumour_purity/results/TCGA_BRCA/tcga_purity_pipeline_20260216_145147.log
clear
cd /home/chu25/TCGA/all_tcga_count_data/tumour_purity/results/TCGA_BRCA/TCGA-BRCA
ls
cd ~/data/
ls
cd pri_bl/
ls
cd GSE199413/
ls
cd results/
ls
cd agre
cd aggregated/
ls
cd ~/data/bl/
ls
cd GSE199413/
ls
cd results/
ls
cd aggregated/
ls
cd tumour_purity/
ls
rm -rf ~/data/bl/GSE199413
clear
Rscript /home/chu25/TCGA/hema_pipeline/hema_pipeline/scripts/hema_tumour_purity_estimation.R   --expr      /home/chu25/data/pri_bl/GSE199413/results/aggregated/raw_count_matrix.rds   --map       /home/chu25/data/pri_bl/GSE199413/results/aggregated/ensembl_to_hgnc.tsv   --outdir    /home/chu25/data/pri_bl/GSE199413/results/aggregated/tumour_purity   --threshold 0.7   --meta      /home/chu25/data/pri_bl/GSE199413/results/aggregated/sample_metadata.csv
ls
cd ..
cd ~/data/
ls
cd bl
ls
cd ..
rm -r  bl
clear
ls
cd pri_bl/
ls
cd GSE199413/
ls
cd results/
ls
cd aggregated/
ls
cd tumour_purity/
ls
pwd
cd ..
ls
cd tcga/
clear
ls
cd tcga_dlbc/
ls
cd dlbc_tumour_purity_results/
ls
cd ..
cd tcga_dlbc_tumour_purity_results/
ls
cd ..
cd tumour_purity_results0.5/
ls
cd ..
ls
cd ~/TCGA/
clear
ls
cd tcga_dlbc/
ls
cd tumour_purity/
lsl
ls
pwd
ls
cd ..
ls
sbatch tumour_purity_estimation.sh 
cd ../logs/
ls
cat tcga_dlbc_count_data_25105842.out 
ls
clear
ls
cd ..
ls
cd tcga_dlbc/
ls
cd logs/
ls
cat dlbc_purity_analysis_28760412.out 
ls -lh
cat dlbc_purity_analysis_28760412.out 
cat dlbc_purity_analysis_28760412.err 
clear
ls
rm dlbc_purity_analysis_2*
ls
ls -lh
squeue -u chu25
cat dlbc_purity_analysis_28760435.
cat dlbc_purity_analysis_28760435.err 
cat  /home/chu25/data/tcga/tcga_dlbc/tumour_purity_results0.7/tcga_dlbc_purity_pipeline_20260216_180116.log
cd ..
clear
ls
cd tcga
ls
cd .
cd ..
cd all_tcga_count_data/
ls
cd tumour_purity/
ls
cat TCGA_LAML.R 
cd results/
ls
cd TCGA_LAML/
ls
cd TCGA-LAML/
ls
pwd
ls
cd ..
clear
ls
cd ..
ls
cd TCGA_LAML/
ls
cd TCGA-LAML/
ls
cd ..
ls
cd ..
ls
cd ..
ls
cd TCGA/
ls
cd results/
ls
cd tumour_purity_results/
ls
cd --
cd TCGA/all_tcga_count_data/
ls
cd tumour_purity/
ls
cd results/
l
ls
cd TCGA_BRCA/
ls
cat tcga_purity_pipeline_20260216_145147.log 
clear
cd ..
ls
cd ..
ls
cd rbl/
ls
cd logs/
ls
cat rbl_purity_analysis_28755885.out 
clear
ls
cat rbl_purity_analysis_28755885.out 
clear
ls
cat rbl_purity_analysis_25117236.out 
cd ..
cd nbl/
ls
cd logs/
ls
cat nbl_purity_analysis_28756345.out 
Rscript --vanilla -e '
x<-readRDS("/home/chu25/data/tcga/ALL_TCGA_STAR_Counts_SummarizedExperiment_filtered_TCGA-LAML.rds");
m<-if(inherits(x,"SummarizedExperiment")) SummarizedExperiment::assay(x) else as.matrix(x);
cat("TCGA rowname examples:\n"); print(head(rownames(m), 5));
cat("Looks like Ensembl?", all(grepl("^ENSG", head(rownames(m), 50))), "\n");
'
clear
conda env list
conda env lists
conda actiavte tcga-r-env
conda env list
conda activate tcga-r-env
clear
Rscript --vanilla -e '
x<-readRDS("/home/chu25/data/tcga/ALL_TCGA_STAR_Counts_SummarizedExperiment_filtered_TCGA-LAML.rds");
m<-if(inherits(x,"SummarizedExperiment")) SummarizedExperiment::assay(x) else as.matrix(x);
cat("TCGA rowname examples:\n"); print(head(rownames(m), 5));
cat("Looks like Ensembl?", all(grepl("^ENSG", head(rownames(m), 50))), "\n");
'
clear
ls /home/chu25/dsmz/leu_lym_dsmz_analysis/results
clear
ls
cd colDataAnnData/
ls
cd results/
ls
tree
clear
cd /home/chu25/colDataAnnData/scripts/
sbatch /home/chu25/colDataAnnData/scripts/pan_cancer_umap.sh
cd ..
cd logs/
ls
cat pan_cancer_umap_27278185.out 
cat pan_cancer_umap_27278185.err 
squeue -u chu25
clear
ls
cat pan_cancer_umap_27278185.err 
clear
cd ..
ls
cd results/
ls
cd pan_cancer_umap_main/
ls
pwd
cd ..
rm -r pan_cancer_umap_main/
clear
cd ..
cd /home/chu25/colDataAnnData/scripts/
sbatch /home/chu25/colDataAnnData/scripts/pan_cancer_umap.sh
._Kpäü#+
cd ../logs/
ls
rm pan_cancer_umap_27278185.*
ls
cat pan_cancer_umap_27282383.err 
-öä+
.-ö#+
cat pan_cancer_umap_27282383.err 
clear
cat pan_cancer_umap_27282383.err 
squeue -u chu25
cat pan_cancer_umap_27282383.err 
clear
cat pan_cancer_umap_27282383.err 
clear
cat pan_cancer_umap_27282383.err 
clear
cat pan_cancer_umap_27282383.err 
clear
squeue -u chu25
cat pan_cancer_umap_27282383.err 
ls
rm pan_cancer_umap_27282383.*
clear
cd ..
cd results/
ls
rm -r pan_cancer_umap_main/
cd ..
cd /home/chu25/colDataAnnData/scripts/
sbatch /home/chu25/colDataAnnData/scripts/pan_cancer_umap.sh
cd log
cd ..
ls
cd logs/
ls
cat pan_cancer_umap_27290667.err 
cat pan_cancer_umap_27290667.out 
clear
cat pan_cancer_umap_27290667.err 
clear
cat pan_cancer_umap_27290667.err 
ls
rm pan_cancer_umap_27290667.*
clear
cd ..
ls
cd results/
ls
rm -r pan_cancer_umap_main/
cd ../scripts/
clear
sbatch /home/chu25/colDataAnnData/scripts/pan_cancer_umap.sh
cd ..
cd logs/
ls
cat pan_cancer_umap_27295064.err 
cat pan_cancer_umap_27295064.out 
squeue -u chu25
cat pan_cancer_umap_27295064.out 
squeue -u chu25
cat pan_cancer_umap_27295064.out 
squeue -u chu25
clear
squeue -u chu25
cat pan_cancer_umap_27295064.out 
clear
cat pan_cancer_umap_27295064.out 
squeue -u chu25
cat pan_cancer_umap_27295064.out 
squeue -u chu25
clear
squeue -u chu25
cat pan_cancer_umap_27295064.out 
clear
cat pan_cancer_umap_27295064.out 
clear
cat pan_cancer_umap_27295064.out 
clear
ls
rm pan_cancer_umap_27295064.*
cd ,,
cd ..
ls
cd results/
ls
cd pan_cancer_umap_main/
ls
cdc ..
cd ..
mkdir pan_cancer_umapp_archieve
mv pan_cancer_umap_main/ pan_cancer_umapp_archieve/
ls
cd ..
cd scripts/
ls
cat pan_cancer_main.sh 
cat pan_cancer_umap.sh 
sbatch pan_cancer_umap.sh 
cd ..
ls
cd logs/
clear
ls
cat pan_cancer_umap_27417093.err 
cat pan_cancer_umap_27417093.out 
cat pan_cancer_umap_27417093.err 
cat pan_cancer_umap_27417093.out 
clear
cat pan_cancer_umap_27417093.out 
cat pan_cancer_umap_27417093.err 
cat pan_cancer_umap_27417093.out 
clear
cat pan_cancer_umap_27417093.err 
clear
cd ../results/
ls
cd pan_cancer_umapp_archieve2/
pwd
ls
cd ..
ls
cd pan_cancer_umap_main/
pwd
clear
cd ..
ls
cd dsmz/
ls
cd leu_lym_dsmz_analysis/
ls
cd preprocessing/
ls
cd results/
ls
cd ..
clear
ls
cd ..
ls
cd results/
ls
pwd
cd ..
ls
cd preprocessing/
ls
cat tcga_dsmz_base_functions.R 
cd logs/
ls
cat leu_lym_dsmz_io_25750900.out 
sbatch /home/chu25/TCGA/all_tcga_count_data/tumour_purity/tumour_tcga_laml_purity.sh
cd home/chu25/TCGA/logs
cd ~/home/chu25/TCGA/logs
cd ~/home/chu25/TCGA/logs/
cd ~/home/chu25/TCGA/
cd ..
cd TCGA/
clar
clear
ls
cd logs/
ls
clear
ls
squeue -u chu25
scancel  27764978
scancel  27766897
squeue -u chu25
rm tumour_brca_purity_analysis_tumour_brca_purity_2776*
ls -lh
clear
cd ..
ls
cd logs/
ls
cat tumour_brca_purity_analysis_tumour_brca_purity_27766903.out 
cat tumour_brca_purity_analysis_tumour_brca_purity_27766903.err 
cat tumour_laml_purity_analysis_tumour_brca_purity_27766903.err
cat tumour_laml_purity_analysis_tumour_brca_purity_27766903.out
cat tumour_laml_purity_analysis_tumour_laml_purity_27766902.out 
cat tumour_laml_purity_analysis_tumour_laml_purity_27766902.err 
cat tumour_brca_purity_analysis_tumour_brca_purity_27766903.err 
cat tumour_brca_purity_analysis_tumour_brca_purity_27766903.out 
cat tumour_laml_purity_analysis_tumour_laml_purity_27766902.out 
squeue -u chu25
cat tumour_laml_purity_analysis_tumour_laml_purity_27766902.out 
squeue -u chu25
cat tumour_laml_purity_analysis_tumour_laml_purity_27766902.err 
cd /home/chu25/dsmz/results/tumor_purity_results/TCGA-LAML
ls
pls
cd ..
tree
cd ..
clear
ls
rm -r tumor_purity_results/
ls
cd ..
ls
cd ~/TCGA/
ls
cd logs/
clear
ls
cat tumour_brca_purity_analysis_tumour_brca_purity_27770658.out 
cat tumour_laml_purity_analysis_tumour_laml_purity_27770659.out 
cat tumour_laml_purity_analysis_tumour_laml_purity_27770659.err 
squeue -u chu25
cd ..
ls
cd tcga_dlbc/
ls
cd logs/
ls
cat dlbc_purity_analysis_25121562.err 
clear
ls
cd ..
ls
cd all_tcga_count_data/
ls
cd tumour_purity/
ls
squeue -u chu25
cd ..
ls
cd logs/
ls
cd ..
ls
cd logs/
ls
cd ..
cd tumour_purity/
ls
cd ../..
cd logs/
ls
cat tumour_laml_purity_analysis_tumour_laml_purity_27770659.err 
cat tumour_brca_purity_analysis_tumour_brca_purity_27770658.err 
rm tumour_*
cat tumour_brca_purity_analysis_tumour_brca_purity_27770658.err 
cat tumour_brca_purity_analysis_tumour_brca_purity_27772608.err 
clear
ls
cat tumour_laml_purity_analysis_tumour_laml_purity_27772607.err 
cat tumour_laml_purity_analysis_tumour_laml_purity_27772607.out 
cat tumour_laml_purity_analysis_tumour_laml_puri
cat tumour_laml_purity_analysis_tumour_laml_purity_27772607.err 
rm tumour_*
clear
ls
cat tumour_laml_purity_analysis_tumour_laml_purity_27776389.err 
rm tumour_*
clear
ls
cat tumour_brca_purity_analysis_tumour_brca_purity_27781050.err 
cat tumour_laml_purity_analysis_tumour_laml_purity_27781049.err 
cat tumour_laml_purity_analysis_tumour_laml_purity_27781049.out 
cat tumour_brca_purity_analysis_tumour_brca_purity_27781050.out 
squeue -u chu25
cat tumour_brca_purity_analysis_tumour_brca_purity_27781050.out 
cat tumour_laml_purity_analysis_tumour_laml_purity_27781049.out 
cat tumour_laml_purity_analysis_tumour_laml_purity_27781049.err 
cat tumour_brca_purity_analysis_tumour_brca_purity_27781050.out 
squeue -u chu25
cat tumour_brca_purity_analysis_tumour_brca_purity_27781050.out 
cat tumour_laml_purity_analysis_tumour_laml_purity_27781049.err 
clar
clear
ls
cd ..
ls
cd nbl
ls
cd results/
ls
cd tumour_purity_results/
s
clear
ls
ls -lh
pwd
ls
cd ..
cd .
cd ..
clear
ls
cd mpn
ls
cd logs/
ls
cd ..
cd data/
ls
cd mpn/
ls
cd GSE228758/
ls
cd results/
ls
cd aggregated/
ls
cd tumour_purity/
ls
pwd
clear
cd ..
cd .
cd ~/colDataAnnData/
ls
pwd
cd ..
clear
ls
cd data/
ls
cd GTF
ls
ls -lh
pwd
clear
ls
cd ..
l
cd ..
cler
clear
ls
cd dsmz/
ls
cd ~/data/GTF/
ls
cd ~/TCGA/
ls
cd tcga
ls
cd ../all_tcga_count_data/
ls
cd tumour_purity/
ls
cat tumour_brca_purity.sh 
ls
pwd
clear
ls
cat TCGA_BRCA.R 
clear
ls
cat tumour_purity.sh 
cat tumour_brca_purity.sh 
Rscript tumour_brca_purity.sh 
sbatch tumour_brca_purity.sh 
clear
ls
cd ..
ls
cd ..
ls
cd rbl
ls
cd results/
ls
cd tumour_purity_results/
ls
cd ~/data/
ls
cd rbl
clear
ls
cd results/
ls
cd tumour_purity_results/
ls
conda env list
cd ~/data
clear
ls
cd tcga/g
cd tcga/
ls
cd ~/data/cll/
ls
cd GSE237907/
ls
cd results/
ls
cl aggregated/
ls
cd aggregated/
ls
cd tumour_purity/
ls
pwd
cd colDataAnnData/
ls
cd scripts/
ls
cd leu_lym/
ls
cd ..
clear
ls
cd ..
ls
cd logs/
ls
cat pan_cancer_umap_27429748.err 
cat pan_cancer_umap_27429748.out 
clear
cat pan_cancer_umap_27429748.out 
clear
cd ..
mkdir logs_new
mn logs logs_new/
mv logs logs_new/
ls
mkdir logs_new/
mkdir logs
mv logs_new logs_old
mkdir results_old
mv results results_old/
mkdir results_old/
mkdir results
ls
clear
ls
cd ..
ls
cd colDataAnnData/
ls
cd scripts/
clear
ls
sbatch pan_cancer_main.sh 
cursor pan_cancer_main.sh 
cd ..
ls
cd logs
ls
cat pan_cancer_umap_29422872.err 
squeue -u chu25
cat pan_cancer_umap_29422872.out 
squeue -u chu25
cat pan_cancer_umap_29422872.out 
squeue -u chu25
ls
cat pan_cancer_umap_29422872.err 
cat pan_cancer_umap_29422872.out 
ls -lhg
rm pan_cancer_umap_29422872.*
clear
cd ..
ls
cd scripts/
ls
sbatch pan_cancer_main.sh 
cd ../logs/
ls
cat pan_cancer_umap_29437128.out 
cat pan_cancer_umap_29437128.err 
squeue -u chu25
cat pan_cancer_umap_29437128.err 
cat pan_cancer_umap_29437128.out 
squeue -u chu25
cat pan_cancer_umap_29437128.out 
cat pan_cancer_umap_29437128.err 
cat pan_cancer_umap_29437128.out 
squeue -u chu25
cat pan_cancer_umap_29437128.out 
squeue -u chu25
cat pan_cancer_umap_29437128.out 
clear
ls
cat pan_cancer_umap_29437128.err 
cat pan_cancer_umap_29437128.out 
cd ..
cd results
ls
cd pan_cancer_umap_main/
ls
ls -lh
cursor Fig_umap_BATCH_CORR_VST_INTEGRATED_cosine.pdf 
ls
cd dsmz/results/
ls
cd tumor_purity_results/
ls
cd ..
ls
cd dsmz_tcga_analysis/
ls
cd results/
ls
clear
cd ..
cd data/
clear
ls
cd tcga/
ls
pwd
cd tcga_dlbc/
ls
cd ..
clear
ls
cd tcga/
ls
mkdir tcga_brca
pwd
ls
cd tcga_brca/
pwd
ls
cd /home/chu25/dsmz/leu_lym_dsmz_analysis/preprocessing
ls
cursor leu_lym_dsmz_main.R 
cd  /home/chu25/data/hema_aggregated/tumour_purity_result
cd /home/chu25/data/hema_aggregated/tumour_purity_result
cd /home/chu25/data/hema_aggregated/tumour_purity_results
ls
cursor purity_tumour_retrieval.R 
clear
cd ..
clear
ls
cd mpn
ls
cd GSE228758/
ls
cd results/
ls
cd aggregated/
ls
cd logs/
ls
cat aggregate_counts.log 
cd ..
cd .
cd ..
ls
cd pri_bl
ls
cd GSE199413/
ls
cd results/
ls
cd aggregated/
ls
curspr final_hema_sample_metadata.csv 
cursor final_hema_sample_metadata.csv 
cd ..
ls
cd icgc_mmml/
ld
ls
tree
cd EGAD00001001672/
ls
cd fastq/
cursor pyega3_output.log 
clear
cd .
cd ..
cd d..
cd ..
cd ~/colDataAnnData/
ls
clear
ls
cd scripts/
ls
cursor pan_cacer_featureset.R
clear
ls
cd ../results/
ls
cd ..
cd /home/chu25/data/cll/GSE237907/results/aggregated/
ls
xs tumour_purity/
ls
clear
ls
cd tumour_purity/
ls
clear
cd ..
clear
ls
cd rbl/
ls
cd results/
ls
cd tumour_purity_results/
ls
cd ~/TCGA/
clear
ls
cd rbl/
ls
cd results/
ls
cd tumour_purity_results/
ls
pwd
clear
cd ..
ls
cd ..
ls
cursor tumour_purity_estimation.R
clear
pwd
ls
cd ..
cd nbl/
ls
cd results/
l
ls
cd tumour_purity_results/
clear
ls
pwd
cd ..
ls
cd ..
ls
cd preprocessing/
ls
cursor tumour.R 
ls
clear
ls
pwd
cd /home/chu25/TCGA/all_tcga_count_data/tumour_purity
ls
cursor TCGA_BRCA.R 
cd ..
cd ../hema_pipeline/
ls
clear
s
ls
cd hema_pipeline/
ls
cd scripts/
ls
cd ..
cursor hema_pipeline.sh 
cd ~/dsmz/
clear
ls
cd dsmz_rbl_analysis/
ls
cd preprocessing/
ls
curosr rbl_main.R 
cursor rbl_main.R 
cursor rbl_io.R
cd data/
ls
cd ..
cd TCGA/
ls
cd all_tcga_count_data/
ls
cd tumour_purity/
ls
cd results/
ls
cursor tcga_estimate_purity.Rmd 
ls
clear
cursor check_purity.R 
ls
cd ..
ls
cursor tcga_estimate_purity_20602524.err 
cursor tcga_estimate_purity_20602524.out 
cursor tcga_estimate_purity_20656149.err 
cursor tcga_estimate_purity_20656149.out 
ls
rm tcga_estimate_purity_206*
clear
ls
cursot tumour_purity.sh 
cursor tumour_purity.sh 
cd ..
ls
cd tcga
ls
cd ..
cd all_tcga_count_data/
clear
ls
cd logs/
ls
cd ..
cd tumour_purity/
ls
cd results/
ls
cursor tumour_purity.sh 
cd ..
clear
ls
cursor tumour_purity.R
ls
cd results/
ls
cd ..
ls
cursor tumour_brca_purity.sh
ls
cursor tumour_purity.sh 
clear
cd ..
cd .
cd ..
ls
cd rbl
ls
pwd
cursor tumour_purity_estimation.R
clear
cd ..
ls
cd cll
ls
cd ../hema_pipeline/
ls
cd hema_pipeline/
ls
cd scripts/
ls
cursor hema_tumour_purity_estimation.R 
cd ..
cursor Snakefile 
clear
ls
cd ..
ls
cd ..
ls
cd mpn
cd ..
cd mpn/
ls
cd ~/data/mpn/
ls
cd GSE228758/
l
clear
ls
cd results/
ls
cd aggregated/
ls
cd tumour_purity/
ls
pwd
clear
ls
cd ..
ls
cd aggregated/
ls
cd ~/TCGA/
ls
cd tcga_dlbc/
clear
ls
cd tumour_purity/
ls
cd ..
ls
cursor tumour_purity_estimation.sh 
cd tumour_purity/
ls
rm *r
rm *
clear
ls
cd ..
ls
pwd
cd ~/data/cll/GSE237907/results/aggregated/tumour_purity/
ls
clear
cd ..
cd ~/colDataAnnData/
ls
cd scripts/
ls
cursor pan_cancer_umap_main.R 
clear
cd ..
cd dsmz/
ls
cd leu_lym_dsmz_analysis/
clear
ls
cd preprocessing/
ls
cd results/
ls
cd ..
clear
ls
cursor io.R
cd logs/
ls
cat leu_lym_dsmz_io_25768257.out 
clear
cat leu_lym_dsmz_io_25768257.out 
clear
cd ..
ls
cat io.sh
ls
cursor leu_lym_dsmz_main.R 
clear
cd ..
cd preprocessing/
ls
sbatch io.sh
cd logs/
ls
cat leu_lym_dsmz_io_29375044.out 
cat leu_lym_dsmz_io_29375044.err 
cat leu_lym_dsmz_io_29375044.out 
squeue -u chu25
cat leu_lym_dsmz_io_29375044.out 
squeue -u chu25
clear
squeue -u chu25
cat leu_lym_dsmz_io_29375044.out 
cat leu_lym_dsmz_io_29375044.err 
rm leu_lym_dsmz_io_29375044.*
cd ..
clear
sbatch io.sh
cd logs/
ls
cat leu_lym_dsmz_io_29375668.err 
cat leu_lym_dsmz_io_29375668.out 
cat leu_lym_dsmz_io_29375668.err 
cd ..
ls
mv normalization.R normalisation.R
cd logs/
ls
rm leu_lym_dsmz_io_29375668.err 
clear
ls -lh
rm leu_lym_dsmz_io_29375668.out 
cd ..
clear
sbatch io.sh 
cd logs/
ls
cat leu_lym_dsmz_io_29376579.err 
cd ..
ls
mv visualize.R visualise.R
cd logs/
rm leu_lym_dsmz_io_29376579.*
cd ..
sbatch io.sh 
squeue -u chu25
cd lgs
ls
cat leu_lym_dsmz_io_29377991.out 
cat leu_lym_dsmz_io_29377991.err 
squeue -u chu25
cat leu_lym_dsmz_io_29377991.err 
cat leu_lym_dsmz_io_29377991.out 
cat leu_lym_dsmz_io_29377991.err 
cat leu_lym_dsmz_io_29377991.out 
cd ..
cd results/
ls
clear
cd ..
ls
cd results/
ls
ls -lh
# See if job is active and using CPU
sstat -j 29377991 --format=JobID,AveCPU,MaxRSS
# Watch the output file in real-time
tail -f logs/leu_lym_dsmz_io_29377991.out
clear
cd ..
ls
cd preprocessing/
ls
cd logs/
cat leu_lym_dsmz_io_29377991.err 
cat leu_lym_dsmz_io_29377991.out 
rm leu_lym_dsmz_io_29377991.*
cd ..
ls
sbatch io.sh 
cd logs/
cat leu_lym_dsmz_io_29378342.out 
cat leu_lym_dsmz_io_29378342.err 
squeue -u chu25
cat leu_lym_dsmz_io_29378342.err 
cat leu_lym_dsmz_io_29378342.out 
squeue -u chu25
cat leu_lym_dsmz_io_29378342.out 
cat leu_lym_dsmz_io_29378342.err 
squeue -u chu25
cat leu_lym_dsmz_io_29378342.err 
cat leu_lym_dsmz_io_29378342.out 
cat leu_lym_dsmz_io_29378342.err 
squeue -u chu25
cd ..
ls
cursor io.sh 
pwd
sbatch io.sh
cd logs/
l
ls
cd ..
cat logs/leu_lym_dsmz_io_29380881.out 
cat logs/leu_lym_dsmz_io_29380881.err 
squeue -u chu25
cat logs/leu_lym_dsmz_io_29380881.err 
cat logs/leu_lym_dsmz_io_29380881.out 
squeue -u chu25
cat logs/leu_lym_dsmz_io_29380881.out 
cat logs/leu_lym_dsmz_io_29380881.err 
cd /home/chu25/dsmz/leu_lym_dsmz_analysis/preprocessing
Rscript --vanilla -e 'source("io.R"); cat("has config? ", exists("config"), "\n"); if (exists("config")) print(config)'
clear
rm logs/leu_lym_dsmz_io_29380881.*
sbatch io.sh
cd logs/
ls
cst leu_lym_dsmz_io_29382583.err 
cat leu_lym_dsmz_io_29382583.err 
cat leu_lym_dsmz_io_29382583.out 
squeue -u chu25
cat leu_lym_dsmz_io_29382583.out 
cat leu_lym_dsmz_io_29382583.err 
cat leu_lym_dsmz_io_29382583.out 
clear
squeue -u chu25
cat leu_lym_dsmz_io_29382583.out 
squeue -u chu25
cat leu_lym_dsmz_io_29382583.out 
clear
cat leu_lym_dsmz_io_29382583.out 
clear
echo "[DEBUG] checking mounts"
ls -lah /home/chu25/data/tcga || true
ls -lah /mnt/beegfs/mnt/work/bioinfo/home/chu25/data/tcga || true
cd ~/data/
ls
cd ~/dsmz/leu_lym_dsmz_analysis/
clear
ls
cd preprocessing/
ls
sbatch io.sh 
cd logs/
ls
cat leu_lym_dsmz_io_29383873.out 
cat leu_lym_dsmz_io_29383873.err 
clear
squeue -u chu25
cat leu_lym_dsmz_io_29383873.out 
cat leu_lym_dsmz_io_29383873.err 
squeue -u chu25
cat leu_lym_dsmz_io_29383873.err 
cat leu_lym_dsmz_io_29383873.out 
clear
cat leu_lym_dsmz_io_29383873.out 
cat leu_lym_dsmz_io_29383873.err 
ls -lh /home/chu25/data/tcga/ALL_TCGA_STAR_Counts_SummarizedExperiment_filtered_TCGA-LAML.rds 
rm leu_lym_dsmz_io_29383873.*
cd ..
clear
sbatch io.sh 
cd logs/
ls
cat leu_lym_dsmz_io_29384451.out 
cat leu_lym_dsmz_io_29384451.err 
squeue -u chu25
cat leu_lym_dsmz_io_29384451.err 
cat leu_lym_dsmz_io_29384451.out 
squeue -u chu25
clear
squeue -u chu25
cat leu_lym_dsmz_io_29384451.err 
rm  rerls
rm leu_lym_dsmz_io_29384451.*
cd ..
clear
sbatch io.sh 
cat logs/leu_lym_dsmz_io_29385707.err 
cat logs/leu_lym_dsmz_io_29385707.out 
squeue -u chu25
cat logs/leu_lym_dsmz_io_29385707.out 
cat logs/leu_lym_dsmz_io_29385707.err 
cd logs/
ls
rm leu_lym_dsmz_io_29385707.*
clear
cd ..
sbatch io.sh 
ls
cd logs/
ls
cat leu_lym_dsmz_io_29386306.out 
cat leu_lym_dsmz_io_29386306.err 
squeue -u chu25
clear
squeue -u chu25
cat leu_lym_dsmz_io_29386306.err 
cat leu_lym_dsmz_io_29386306.out 
cat leu_lym_dsmz_io_29386306.err 
squeue -u chu25
cat leu_lym_dsmz_io_29386306.err 
cat leu_lym_dsmz_io_29386306.out 
clear
cat leu_lym_dsmz_io_29386306.out 
cat leu_lym_dsmz_io_29386306.err 
cd ..
ls
cat tcga_dsmz_base_functions.R 
ls
cat normalisation.R 
clear
ls
cd logs/
ls
rm leu_lym_dsmz_io_29382583.*
cd ..
sbatch io.sh 
cd logs/
ls
cat leu_lym_dsmz_io_29411851.out 
cat leu_lym_dsmz_io_29411851.err 
sbatch io.sh 
sqeueu -u chu25
squeue -u chu25
cat leu_lym_dsmz_io_29411851.err 
cat leu_lym_dsmz_io_29411851.out 
cat leu_lym_dsmz_io_29411851.err 
clear
cat leu_lym_dsmz_io_29411851.err 
cat leu_lym_dsmz_io_29411851.out 
cd lo
ls
rm leu_lym_dsmz_io_257*
ls
rm leu_lym_dsmz_io_293*
ls
cd ..
ls
cd results/
ls
ls -lh
clear
cd ..
ls
cd preprocessing/
ls
cd results/
ls
cd ..
ls
cd results/
ls
clear
pwd
ls -lh
rm HEMA*
ls
ls -lh
rm dispersion_DSMZ_LL100_RAW*
ls
ls -lh
rm batch_evaluation_results.rds 
ls -lh
rm dispersion_TCGA_HEMA_LL_DLBC_LAML_ST_RAW
rm dispersion_TCGA_HEMA_LL_DLBC_LAML_ST_RAW*
rm vst_DSMZ_LL100_RAW.rds 
rm vst_TCGA_HEMA_LL_DLBC_LAML_ST_RAW.rds 
rm tcga_vst_dlbc_laml_batch_corrected.rds 
rm tcga_hema-ll_dlbc_laml_counts.rds 
ls -lh
rm dsmz_leu_lym_counts.rds 
rm dsmz_vst_ll100_batch_corrected.rds 
ls -lh
rm mean_sd_DSMZ_LL100_RAW.pdf 
rm mean_sd_TCGA_HEMA_LL_DLBC_LAML_ST_RAW.pdf 
rm hema-ll_dlbc_laml_batch_removal_metrics.rds 
ls -lh
clear
cd ..
ls
cd preprocessing/
ls
cd logs/
ls
rm leu_lym_dsmz_io_29411851.*
cd ..
sbatch leu_lym_dsmz_main.R 
sbatch io.sh
cd logs/
ls
clear
ls
cat leu_lym_dsmz_io_29413393.err 
cat leu_lym_dsmz_io_29413393.out 
squeue -u chu25
cat leu_lym_dsmz_io_29413393.out 
cd ../../results/
ls
squeue -u chu25
clear
squeue -u chu25
cat leu_lym_dsmz_io_29413393.out 
clear
cat leu_lym_dsmz_io_29413393.out 
squeue -u chu25
clear
squeue -u chu25
ls
cd ..
ls
cd preprocessing/
ls
cd logs/
ls
cat leu_lym_dsmz_io_29413393.err 
cat leu_lym_dsmz_io_29413393.out 
ls
rm leu_lym_dsmz_io_29413393.*
clear
cd ..
ls
sbatch io.sh
cd ..
ls
cd results/
ls
squeue -u chu25
clear
cd ../log
ls
cd ..
ls
cd preprocessing/
ls
cd logs/
clear
ls
cat leu_lym_dsmz_io_29416202.err 
cat leu_lym_dsmz_io_29416202.out 
clear
ls
cd ..
ls
cd results/
ls
cd ..
ls
cd results/
ls
clear
ls -lh
clclear
clear
ls
cd colDataAnnData/
ls
cd scripts/
ls
cursor pan_cancer_umap_main.R 
cursor pan_cancer_umap.sh 
clear
cd ..
ls
cd dsmz/
clear
ls
cd leu_lym_dsmz_analysis/
ls
cd ..
clear
ls
cd TCGA/
clear
ls
cd tcga_dlbc/
ls
cd logs/
ls
cursor dlbc_purity_analysis_25121562.err 
cd ..
cd sample_types/
clear
ls
cd ..
clear
ls
cd logs/
ls
clear
cd ..
ls
cd mpn/
ls
cd logs/
ls
cursor download_mpn_primary_fastq_25393833.err 
cursor download_mpn_primary_fastq_25393833.out 
cd ..
clear
cd ..
ls
cd pri_bl/
ls
cd logs/
ls
cursor download_bl_primary_fastq_25391108.out 
clear
ls
cd ..
ls
cd ..
ls
cd ..
cd /data/tcga
cd data/
clear
ls
cd tcga/
ls
cd ../results/
ls
cd ..
cd hema_aggregated/
ls
cd ../../TCGA/results/
clear
ls
cd tumour_purity_results/
ls
tree
pwd
clear
cd /home/chu25/TCGA/rbl
ls
cursor tumour_purity_estimation.R
cd ../../dsmz/
ls
cd leu_lym_dsmz_analysis/
clear
ls
cd preprocessing/
ls
cursor leu_lym_dsmz_main.R 
cursor io.R
cd ~/colDataAnnData/
clear
ls
cd scripts/
ls
mkdir leu_lym
cd leu_lym/
cursor leu_lym_io.R
clear
