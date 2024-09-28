#----------------------- ADMIXTURE: global ancestry --------------------------------#
# A. Config file
configfile:"config_admixture.yaml"

# B. Define wildcard
wildcard_constraints:
    panel = config["panel"], 
    folder = config["dir"],
    dataset = config["dataset"],
    #K = range(config["Kstart"], config["Kend"]) 
    

# C. Rule all
rule all:
     input:
        expand("{folder}/results/{folder}.{panel}.{dataset}.merged.{K}.{ext}", folder = config["dir"], panel = config["panel"], dataset = config["dataset"], K = range(config["Kstart"], config["Kend"]), ext =["Q", "P", "log"]),
        expand("{folder}/files/{folder}.{panel}.{dataset}.merged.{ext}", folder = config["dir"], panel = config["panel"], dataset = config["dataset"], ext =["bed", "bim", "fam", "log"]),
        expand("{folder}/files/{panel}_all_QCed.{ext}", folder = config["dir"], panel = config["panel"], ext =["bed", "bim", "fam"]),
        expand("{folder}/files/{panel}_N_SNPs_variants_final.txt", folder = config["dir"], panel = config["panel"]),
        expand("{folder}/population/{panel}_{dataset}.population.txt", folder = config["dir"], panel = config["panel"], dataset = config["dataset"]),
        expand("{folder}/files/{panel}_{dataset}_N_SNPs_variants.txt", folder = config["dir"], panel = config["panel"], dataset = config["dataset"]),
        expand("{folder}/results/K/{panel}.{dataset}.cv.error", folder = config["dir"], panel = config["panel"], dataset = config["dataset"]),
        expand("{folder}/results/K/{panel}.{dataset}.Kplot.pdf", folder = config["dir"], panel = config["panel"], dataset = config["dataset"]),
        expand("{folder}/results/K/{panel}.{dataset}.{K}.perID.txt", folder = config["dir"], panel = config["panel"], dataset = config["dataset"], K = range(config["Kstart"], config["Kend"])),
        expand("{folder}/results/K/{panel}.{dataset}.{K}.popmean.txt", folder = config["dir"], panel = config["panel"], dataset = config["dataset"], K = range(config["Kstart"], config["Kend"])),
        expand("{folder}/results/K/{panel}.{dataset}.{K}.admixtureplot.pdf", folder = config["dir"], panel = config["panel"], dataset = config["dataset"], K = range(config["Kstart"], config["Kend"])),
        expand( "{folder}/config/config_admixture.yaml", folder = config["dir"])

#----------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------ Reference Sample (1KGenome) -----------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------------#

rule select_samples:
    input:
        "../1KGenomes/1000GP_Phase3/1000GP_Phase3.sample"
    output:
        ref_pop = "{folder}/population/{panel}.sample",
        ref_list = "{folder}/population/{panel}.ids"
    params:
        ref_sample = config["ref_samples"]
    shell:
        """
        mkdir -p {wildcards.folder}/population/
        egrep "{params.ref_sample}" {input} > {output.ref_pop}
        awk '{{print $1}}' {output.ref_pop} > {output.ref_list}
        """
# Select populations to be used as reference from 1KGenomes

rule select_ref:
    input:
        ref_list = "{folder}/population/{panel}.ids",
        vcf = "../1KGenomes/Reference_vcf/ALL.{ch}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    output:
        subset_vcf = temp("{folder}/temp/{panel}_{ch}.vcf.gz")
    conda: "./envs/bcftools.yaml"
    shell:
        """
        bcftools view -S {input.ref_list} {input.vcf} -m2 -M2 -v snps -Oz -o {output.subset_vcf}
        """
# Select the samples from the vcf file (one per chr) and only keep biallelic variants

rule vcf_plink:
    input:
        subset_vcf = "{folder}/temp/{panel}_{ch}.vcf.gz"
    output:
        subset_bim = temp("{folder}/temp/{panel}_{ch}.bim"),
        subset_bed = temp("{folder}/temp/{panel}_{ch}.bed"),
        subset_fam = temp("{folder}/temp/{panel}_{ch}.fam")
    conda: "./envs/plink.yaml"
    params:
        maf = config["maf_ref"],
        out = temp("{folder}/temp/{panel}_{ch}")
    shell:
        """
        plink --vcf {input.subset_vcf} --const-fid 0 --maf {params.maf} --make-bed --out {params.out}
        """
# Convert each vcf file to plink binaries (--const-fid sets FID for all individuals) (Filtering out MAF thr)

rule change_bim:
    input:
        subset_bim = "{folder}/temp/{panel}_{ch}.bim",
        subset_bed = "{folder}/temp/{panel}_{ch}.bed",
        subset_fam = "{folder}/temp/{panel}_{ch}.fam"
    output:
        changed_bim = temp("{folder}/temp/{panel}_{ch}_changed.bim"),
        changed_bed = temp("{folder}/temp/{panel}_{ch}_changed.bed"),
        changed_fam = temp("{folder}/temp/{panel}_{ch}_changed.fam")

    shell:
        """
        awk '{{print $1 , $1":"$4 , $3 , $4 , $5 , $6}}' {input.subset_bim} > {output.changed_bim}
        cp {input.subset_bed} {output.changed_bed}
        cp {input.subset_fam} {output.changed_fam}
        """
# VCFs don't have rsIDs, so we change bim file 

rule remove_dup:
    input:
        changed_bim = "{folder}/temp/{panel}_{ch}_changed.bim",
        changed_bed = "{folder}/temp/{panel}_{ch}_changed.bed",
        changed_fam = "{folder}/temp/{panel}_{ch}_changed.fam"
    output:
        dups = temp("{folder}/temp/{panel}_{ch}.dupvar"),
        nodup_bim = temp("{folder}/temp/{panel}_{ch}_nodup.bim"),
        nodup_bed = temp("{folder}/temp/{panel}_{ch}_nodup.bed"),
        nodup_fam = temp("{folder}/temp/{panel}_{ch}_nodup.fam"),
        nodup_nosex = temp("{folder}/temp/{panel}_{ch}_nodup.nosex")
    conda: "./envs/plink.yaml"
    params: 
        inp = "{folder}/temp/{panel}_{ch}_changed",
        out = "{folder}/temp/{panel}_{ch}_nodup",
        dupvar = "{folder}/temp/{panel}_{ch}"
    shell:
        """
        plink --bfile {params.inp} --list-duplicate-vars --out {params.dupvar}
        plink --bfile {params.inp} --exclude {output.dups} --make-bed --out {params.out}
        """
# Removed duplicated variants to avoid issues in merging (recommended by plink)

rule merge_list:
    input:
        expand("{folder}/temp/{panel}_{ch}_nodup.{ext}", folder = config["dir"], panel = config["panel"], ext = ["bed", "bim", "fam"], ch=[f"chr{i}" for i in range(1, 23)])
    output:
        merge_list_1 = temp("{folder}/temp/{panel}.merge_list_temp.txt"),
        merge_list = temp("{folder}/temp/{panel}.merge_list.txt")
    shell:
        r"""
        ls {input} > {output.merge_list_1}
        sed 's/\(.*\)\..*/\1/' {output.merge_list_1} | sort -u > {output.merge_list}
        """
# r""" “raw string” leaves the interpretation of the code to bash shell. It's better to make sure that the ser command will work
# Make list of all chr files to be merged

rule check_merge:
    input:
        binaries = expand("{folder}/temp/{panel}_{ch}_nodup.{ext}", folder = config["dir"], panel = config["panel"], ext = ["bed", "bim", "fam"], ch=[f"chr{i}" for i in range(1, 23)]),
        mergelist = "{folder}/temp/{panel}.merge_list.txt"
    output:
        temp("{folder}/temp/{panel}_all_check.missnp")
    params: 
        out = "{folder}/temp/{panel}_all_check"
    conda: "./envs/plink.yaml"
    shell:
        """
        set +e
        plink --merge-list {input.mergelist} --out {params.out}
        set -e
        """
# Create list of SNPs that will fail merge. The command would automatically generate binaries. But since it fails and we are only interested in the 
    # misssnp file, I added set +e and set -e so that snakemake doesn't crash (disable the exit-on-error behavior)

rule exclude_snps:
    input:
        nodup_bim = "{folder}/temp/{panel}_{ch}_nodup.bim",
        nodup_bed = "{folder}/temp/{panel}_{ch}_nodup.bed",
        nodup_fam = "{folder}/temp/{panel}_{ch}_nodup.fam",
        exclude = "{folder}/temp/{panel}_all_check.missnp"
    output:
        rm_bim = temp("{folder}/temp/{panel}_{ch}_rm.bim"),
        rm_bed = temp("{folder}/temp/{panel}_{ch}_rm.bed"),
        rm_fam = temp("{folder}/temp/{panel}_{ch}_rm.fam"),
        rm_nonex = temp("{folder}/temp/{panel}_{ch}_rm.nosex")
    params:
        ipt = "{folder}/temp/{panel}_{ch}_nodup",
        out = "{folder}/temp/{panel}_{ch}_rm"
    conda: "./envs/plink.yaml"
    shell:
        """
        plink --bfile {params.ipt} --exclude {input.exclude} --make-bed --out {params.out}
        """
# Exclude problematic snps that will keep files from merging

rule merge_list_final:
    input:
        expand("{folder}/temp/{panel}_{ch}_rm.{ext}", folder = config["dir"], panel = config["panel"], ext = ["bed", "bim", "fam"], ch=[f"chr{i}" for i in range(1, 23)])
    output:
        merge_list_1 = temp("{folder}/temp/{panel}.merge_list_final_temp.txt"),
        merge_list = temp("{folder}/temp/{panel}.merge_list_final.txt")
    shell:
        r"""
        ls {input} > {output.merge_list_1}
        sed 's/\(.*\)\..*/\1/' {output.merge_list_1} | sort -u > {output.merge_list}
        """
# Create a list for merge after excluding problematic snps

rule merge_chr:
    input:
        binaries = expand("{folder}/temp/{panel}_{ch}_rm.{ext}", folder = config["dir"], panel = config["panel"], ext = ["bed", "bim", "fam"], ch=[f"chr{i}" for i in range(1, 23)]),
        mergelist = "{folder}/temp/{panel}.merge_list_final.txt"     
    output:
        bim = "{folder}/files/{panel}_all.bim",
        bed = "{folder}/files/{panel}_all.bed",
        fam = "{folder}/files/{panel}_all.fam"
    params: 
        out = "{folder}/files/{panel}_all"
    conda: "./envs/plink.yaml"
    shell:
        """
        plink --merge-list {input.mergelist} --make-bed --out {params.out}
        """
# Merge all chromosomes

rule qc_panel:
    input:
        multiext("{folder}/files/{panel}_all", ".bed", ".bim", ".fam")
    output:
        multiext("{folder}/files/{panel}_all_QCed", ".bed", ".bim", ".fam")
    params:
        ipt = "{folder}/files/{panel}_all",
        out = "{folder}/files/{panel}_all_QCed",
        maf = config["qcmaf"],
        mind = config["qcmind"],
        hwe = config["qchwe"],
        geno = config["qcgeno"]
    conda: "./envs/plink.yaml"
    shell:
        """
        plink --bfile {params.ipt} --geno {params.geno} --mind {params.mind} --maf {params.maf} --hwe {params.hwe} --make-bed --out {params.out}
        """
    # Standart QC for reference panel

rule final_report:
    input:
        initial_bim = "{folder}/files/{panel}_all.bim",
        initial_fam = "{folder}/files/{panel}_all.fam",
        final_bim = "{folder}/files/{panel}_all_QCed.bim",
        final_fam = "{folder}/files/{panel}_all_QCed.fam"
    conda: "./envs/R.yaml"
    output:
        report = "{folder}/files/{panel}_N_SNPs_variants_final.txt"
    script:
        "scripts/variants_individuals.R"
#Report N individuals and SNPs after QC

#----------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------ Target Sample--------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------------#

rule prunning:
    input:
        multiext("../QC/Binary_final/{dataset}_QCed_final", ".bed", ".bim", ".fam")
    output:
        prune_in = "{folder}/files/{dataset}.prune.in",
        pruned = multiext("{folder}/files/{dataset}_prunned", ".bed", ".bim", ".fam")     
    params:
        ipt = "../QC/Binary_final/{dataset}_QCed_final",
        prune = "{folder}/files/{dataset}",
        out = "{folder}/files/{dataset}_prunned",
        window = config["prune_window"],
        step = config["prune_step"],
        r2 = config["prune_r2"]
    conda: "./envs/plink.yaml"
    shell:
        """
        plink \
        --bfile {params.ipt} \
        --indep-pairwise {params.window} {params.step} {params.r2} \
        --out {params.prune}

        plink \
        --bfile {params.ipt} \
        --extract {output.prune_in} \
        --make-bed \
        --out {params.out}
        """
# Prune variants from target sample (recommended by admixture)
rule change_target_bim:
    input:
        bim = "{folder}/files/{dataset}_prunned.bim",
        bed = "{folder}/files/{dataset}_prunned.bed",
        fam = "{folder}/files/{dataset}_prunned.fam"
    output:
        bim = "{folder}/files/{dataset}.bim",
        bed = "{folder}/files/{dataset}.bed",
        fam = "{folder}/files/{dataset}.fam"
    shell:
        """
        awk '{{print $1 , $1":"$4 , $3 , $4 , $5 , $6}}' {input.bim} > {output.bim}
        cp {input.bed} {output.bed}
        cp {input.fam} {output.fam}
        """

rule extract_snps_A:
    input:
        target = multiext("{folder}/files/{dataset}", ".bed", ".bim", ".fam"),
        ref = multiext("{folder}/files/{panel}_all_QCed", ".bed", ".bim", ".fam")
    output:
        snplist1 = "{folder}/temp/{panel}.{dataset}.snplist1.txt",
        snplist2 = "{folder}/temp/{panel}.{dataset}.snplist2.txt",
        target = multiext("{folder}/files/{dataset}.{panel}.target", ".bed", ".bim", ".fam"),
        ref = multiext("{folder}/files/{panel}.{dataset}_all.ref", ".bed", ".bim", ".fam")
    conda: "./envs/plink.yaml"
    params:
        target = "{folder}/files/{dataset}.bim",
        target_plink = "{folder}/files/{dataset}",
        ref = "{folder}/files/{panel}_all_QCed",
        reffinal = "{folder}/files/{panel}.{dataset}_all.ref.bim",
        reffinal_plink = "{folder}/files/{panel}.{dataset}_all.ref",
        targetfinal = "{folder}/files/{dataset}.{panel}.target"
    shell:
        """
        awk '{{print $2}}' {params.target} > {output.snplist1}
        plink --bfile {params.ref} --extract {output.snplist1} --make-bed --out {params.reffinal_plink}
        awk '{{print $2}}' {params.reffinal} > {output.snplist2}
        plink --bfile {params.target_plink} --extract {output.snplist2} --make-bed --out {params.targetfinal}    
        """
# Extract snps in common in both datasets and merge plink files

rule check_flip_merge:
    input:
        target = multiext("{folder}/files/{dataset}.{panel}.target", ".bed", ".bim", ".fam"),
        ref = multiext("{folder}/files/{panel}.{dataset}_all.ref", ".bed", ".bim", ".fam")
    output:
        snp = temp("{folder}/temp/{panel}.{dataset}.testing.missnp"),
        ref_flipped =  multiext("{folder}/files/{panel}.{dataset}_all_flipped.ref", ".bed", ".bim", ".fam")
    conda: "./envs/plink.yaml"
    params:
        target = "{folder}/files/{dataset}.{panel}.target",
        ref = "{folder}/files/{panel}.{dataset}_all.ref",
        out = "{folder}/temp/{panel}.{dataset}.testing",
        ref_flipped = "{folder}/files/{panel}.{dataset}_all_flipped.ref"
    shell:
        """
        set +e
        plink --bfile {params.target} --bmerge {params.ref} --out {params.out}
        set -e
        plink --bfile {params.ref} --flip {output.snp} --make-bed --out {params.ref_flipped}
        """
# Merge two datasets that now have intersecting snps

rule target_ref:
    input:
        target = multiext("{folder}/files/{dataset}.{panel}.target", ".bed", ".bim", ".fam"),
        ref = multiext("{folder}/files/{panel}.{dataset}_all_flipped.ref", ".bed", ".bim", ".fam")
    output:
        multiext("{folder}.{panel}.{dataset}.merged", ".bed", ".bim", ".fam")
    conda: "./envs/plink.yaml"
    params:
        target = "{folder}/files/{dataset}.{panel}.target",
        ref = "{folder}/files/{panel}.{dataset}_all_flipped.ref",
        out = "{folder}.{panel}.{dataset}.merged"
    shell:
        """
        plink --bfile {params.target} --bmerge {params.ref} --out {params.out}
        """
# Merge two datasets that now have intersecting snps

rule report:
    input:
        final_bim = "{folder}.{panel}.{dataset}.merged.bim",
        final_fam = "{folder}.{panel}.{dataset}.merged.fam"
    conda: "./envs/R.yaml"
    output:
        report = "{folder}/files/{panel}_{dataset}_N_SNPs_variants.txt"
    script:
        "scripts/variants_individuals_merge.R"
# Report how many variants and individuals are in the final merged file

rule populations:
    input:
        fam = "{folder}.{panel}.{dataset}.merged.fam",
        fambr = "{folder}/files/{dataset}.fam",
        pop1k = "../1KGenomes/1000GP_Phase3/1000GP_Phase3.sample"
    output:
        file = "{folder}/population/{panel}_{dataset}.population.txt"
    conda: "./envs/R.yaml"
    script:
        "scripts/population_file.R"
# Making a file with the population from all Ids

#----------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------ ADMIXTURE --------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------------#

rule admixture:
    input:
        bed = "{folder}.{panel}.{dataset}.merged.bed"
    output:
        q = "{folder}.{panel}.{dataset}.merged.{K}.Q",
        p = "{folder}.{panel}.{dataset}.merged.{K}.P",
        log = "{folder}.{panel}.{dataset}.merged.{K}.log",
    conda: "./envs/admixture.yaml"
    benchmark:
        "{folder}/benchmarks/{dataset}.{panel}.{K}.admixture.benchmark.txt"
    params:
        K = lambda wildcards: [wildcards.K],
        j = 12
    shell:
        """
        admixture --cv {input.bed} {params.K} -j{params.j} > {output.log}
        """
# Run admixture for several Ks

rule mv_results_admixture:
    input:
        q = "{folder}.{panel}.{dataset}.merged.{K}.Q",
        p = "{folder}.{panel}.{dataset}.merged.{K}.P",
        log = "{folder}.{panel}.{dataset}.merged.{K}.log"
    output:
        q = "{folder}/results/{folder}.{panel}.{dataset}.merged.{K}.Q",
        p = "{folder}/results/{folder}.{panel}.{dataset}.merged.{K}.P",
        log = "{folder}/results/{folder}.{panel}.{dataset}.merged.{K}.log"
    params:
        path = "{folder}/results/"
    shell:
        """
        mv {input.q} {input.p} {input.log} {params.path}
        """
# Move admixture results to final directory (output location can't be specified in admixture)

rule mv_results_plink:
    input:
        bed = "{folder}.{panel}.{dataset}.merged.bed",
        bim = "{folder}.{panel}.{dataset}.merged.bim",
        fam = "{folder}.{panel}.{dataset}.merged.fam",
        log = "{folder}.{panel}.{dataset}.merged.log"     
    output:
        bed = "{folder}/files/{folder}.{panel}.{dataset}.merged.bed",
        bim = "{folder}/files/{folder}.{panel}.{dataset}.merged.bim",
        fam = "{folder}/files/{folder}.{panel}.{dataset}.merged.fam",
        log = "{folder}/files/{folder}.{panel}.{dataset}.merged.log"
    params:
        path = "{folder}/files/"
    shell:
        """
        mv {input.bed} {input.bim} {input.fam} {input.log} {params.path}
        """
# Move plink results to final directory 

rule CVerror:
    input:
        expand("{folder}/results/{folder}.{panel}.{dataset}.merged.{K}.log", folder = config["dir"], panel = config["panel"], dataset = config["dataset"], K = range(config["Kstart"], config["Kend"]))
    output:
        "{folder}/results/K/{panel}.{dataset}.cv.error"
    params:
        pattern = "{folder}/results/{folder}.{panel}.{dataset}.merged.*.log"
    shell:
        """
        awk '/CV/ {{print $3,$4}}' {params.pattern} > {output}
        """
# Select CV from all Ks in their log files and combine them in one file

rule plot_K:
    input:
        cv = "{folder}/results/K/{panel}.{dataset}.cv.error"
    output:
        plot = "{folder}/results/K/{panel}.{dataset}.Kplot.pdf"
    conda: "./envs/Rplotting.yaml"
    script: 
        "scripts/plotK.R"
# Plot all Ks

rule plot_admixture:
    input:
        population = "{folder}/"popula"tion/{panel}_{dataset}.population.txt",
        admixture = "{folder}/results/{folder}.{panel}.{dataset}.merged.{K}.Q",
        fam = "{folder}/files/{folder}.{panel}.{dataset}.merged.fam"
    output:
        table = "{folder}/results/K/{panel}.{dataset}.{K}.perID.txt",
        poptable = "{folder}/results/K/{panel}.{dataset}.{K}.popmean.txt",
        plot = "{folder}/results/K/{panel}.{dataset}.{K}.admixtureplot.pdf"
    conda: "./envs/Rplotting.yaml"
    script: 
        "scripts/plotAdmixture.R"
# Plot admixture results per K and generate tables

rule config:
    input:
        "config_admixture.yaml"
    output:
        "{folder}/config/config_admixture.yaml"
    params:
        path = "{folder}/config/"
    shell:
        """
        cp {input} {params.path}
        """
# Copy config file to corresponding folder, so config informations are saved to each run

