#  ShapeIt2.sh -- read aware phasing
#citation: O. Delaneau, B. Howie, A. Cox, J-F. Zagury, J. Marchini (2013) Haplotype estimation using sequence reads. American Journal of Human Genetics 93 (4) 787-696

#  Created by Nina Tombers on 30.09.21.



mkdir /gss/work/abal8898/IBD/phasing/logs/
mkdir /gss/work/abal8898/IBD/phasing/outputs/
mkdir /gss/work/abal8898/IBD/phasing/outputs/1_BAM/
mkdir /gss/work/abal8898/IBD/phasing/outputs/2_singlChrom/
mkdir /gss/work/abal8898/IBD/phasing/outputs/3_PIRs/
mkdir /gss/work/abal8898/IBD/phasing/outputs/4_Phased/
--------------------------------------------------------------------------------
#!/bin/bash
#SBATCH --job-name=1_bam
#SBATCH --partition=carl.p
#SBATCH --output=1_bam__%A_%a.out
#SBATCH --error=1_bam_%A_%a.err
#SBATCH --time=1-14:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=75G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nina.tombers@bluewin.ch

#SBATCH --array=1-24

printf -v i "%02d" $SLURM_ARRAY_TASK_ID

awk -v lg="$i" '{ print $0, "LG"lg }' BAM-all.txt | sort -k 3 > /gss/work/abal8898/IBD/phasing/outputs/1_BAM/bamlist-LG"$i".txt
--------------------------------------------------------------------------------
#!/bin/bash
#SBATCH --job-name=2_splitChr
#SBATCH --partition=carl.p
#SBATCH --output=2_splitChr_%A_%a.out
#SBATCH --error=2_splitChr_%A_%a.err
#SBATCH --time=01-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=64G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nina.tombers@bluewin.ch

#SBATCH --array=1-24

printf -v i "%02d" $SLURM_ARRAY_TASK_ID

  vcftools \
  --gzvcf /gss/work/abal8898/IBD/phasing/filterd_bi-allelic.vcf.gz \
  --chr LG"$i" \
  --recode \
  --stdout |
  bgzip > /gss/work/abal8898/IBD/phasing/outputs/2_singlChrom/filterd_bi-allelic.LG"$i".vcf.gz
  tabix -p vcf /gss/work/abal8898/IBD/phasing/outputs/2_singlChrom/filterd_bi-allelic.LG"$i".vcf.gz

--------------------------------------------------------------------------------
#!/bin/bash
#SBATCH --job-name=3_runPIRs
#SBATCH --partition=carl.p
#SBATCH --time=01-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=64G
#SBATCH --output=3_runPIRs_%A_%a.out
#SBATCH --error=3_runPIRs_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nina.tombers@bluewin.ch

#SBATCH --array=1-24

printf -v i "%02d" $SLURM_ARRAY_TASK_ID

extractPIRs --bam /gss/work/abal8898/IBD/phasing/outputs/1_BAM/bamlist-LG"$i".txt \
            --vcf /gss/work/abal8898/IBD/phasing/outputs/2_singlChrom/filterd_bi-allelic.LG"$i".vcf.gz \
            --out /gss/work/abal8898/IBD/phasing/outputs/3_PIRs/PIRsList-LG"$i".txt \
            --base-quality 20 \
            --read-quality 15

--------------------------------------------------------------------------------
##4.    Phase the VCF using extracted PIRs

#!/bin/bash
#SBATCH --job-name=4_phase
#SBATCH --partition=carl.p
#SBATCH --output=4_phase__%A_%a.out
#SBATCH --error=4_phase_%A_%a.err
#SBATCH --time=01-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=64G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nina.tombers@bluewin.ch

#SBATCH --array=1-24

printf -v i "%02d" $SLURM_ARRAY_TASK_ID

shapeit -assemble --force\
         --input-vcf /gss/work/abal8898/IBD/phasing/outputs/2_singlChrom/filterd_bi-allelic.LG"$i".vcf.gz \
         --input-pir /gss/work/abal8898/IBD/phasing/outputs/3_PIRs/PIRsList-LG"$i".txt \
         -O /gss/work/abal8898/IBD/phasing/outputs/4_Phased/HapTyData-LG"$i" 


shapeit -convert \
        --input-hap /gss/work/abal8898/IBD/phasing/outputs/4_Phased/HapTyData-LG"$i" \
        --output-vcf /gss/work/abal8898/IBD/phasing/outputs/4_Phased/HapTyData-LG"$i".vcf 

   bgzip /gss/work/abal8898/IBD/phasing/outputs/4_Phased/HapTyData-LG"$i".vcf 

--------------------------------------------------------------------------------
##5.    Merge LGs back together

#!/bin/bash
#SBATCH --job-name=5_merge
#SBATCH --partition=carl.p
#SBATCH --output=5_merge_%j.out
#SBATCH --error=5_merge_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nina.tombers@bluewin.ch

ml VCFtools/0.1.14

  vcf-concat \
           /gss/work/abal8898/IBD/phasing/outputs/4_Phased/HapTyData-LG*.vcf.gz | \
grep -v ^\$ | \
tee phased.vcf | \
   bgzip phased.vcf
 tabix -p vcf phased.vcf.gz

