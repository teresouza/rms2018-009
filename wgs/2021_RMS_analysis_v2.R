library(tidyverse)
library(VariantAnnotation)
library(ComplexHeatmap)
library(circlize)
library(ccfindR)
library(devtools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(RColorBrewer)
library(ggsci)
library(broom)
library(ggpubr)
load_all("~/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Freek/MutationalPatterns/")


#Function to get the vafs from a vcf
get_vaf = function(vcf, sample){
    ad = geno(vcf)$AD[,sample]
    vaf = sapply(ad, function(x) x[[2]] / sum(x))
    vaf[is.na(vaf)] = 0 #For sites with 0 reads
    return(vaf)
}

# Read and filter vcf
read_filter_vcf = function(vcf_fname, sample_name, min_vaf){
    # Read vcf
    vcf = readVcf(vcf_fname)
    
    # Determine name of normal
 
    both_samples = samples(header(vcf))
    normal_name = both_samples[both_samples != sample_name]
    
    # Determine quality parameters
    vaf_sample = get_vaf(vcf, sample_name)
    vaf_norm = get_vaf(vcf, normal_name)
    
    dp_sample = lapply(geno(vcf)$AD[, sample_name], sum)
    dp_norm = lapply(geno(vcf)$AD[, normal_name], sum)

    
    # Filter variants
    variant_f = vaf_sample > min_vaf & vaf_norm == 0 & dp_sample >= 20 & dp_norm >= 20
    vcf = vcf[variant_f]
    vaf_sample = vaf_sample[variant_f]
    
    # Add variant type
    csq = info(vcf)$CSQ
    csq = purrr::map(as.list(csq), ~str_split(.x, "\\|")) %>% 
        purrr::map(1)
    impact = purrr::map_chr(csq, 9)
    type = purrr::map_chr(csq, 10)
    
    # Create GRanges
    gr = granges(vcf)
    gr$vaf = vaf_sample
    gr$impact = impact
    gr$type = type
    
    return(gr)
}

do_signature_analysis = function(mut_mat, signatures, signatures_exposure, meta_tb, sigs_used, refit_method, type){
    
    #profile_fig = plot_96_profile(mut_mat)
    
    # Create cosine similarity heatmaps
    
    # cos_sim_m = cos_sim_matrix(mut_mat, mut_mat)
    # cos_sim_heat_fig = plot_cosine_heatmap(cos_sim_m, cluster_cols = F, cluster_rows = F, plot_values = F)
    
    # Create complex heatmap figure
    ha = columnAnnotation(Sample = meta_tb$Sample,
                          Subtype = meta_tb$Subtype,
                          Type = meta_tb$Type,
                          col = list(Type = c("Tumor" = "#f4a71d", "Tumoroid" = "#1760a9", "TumoroidL" = "brown", "Tumoroid2" = "red"),
                                     Subtype = c(c("eRMS" = "#139174", "FN aRMS" = "#f19f89", "P3F aRMS" = "#7ec7db", "P7F aRMS" = "#2f4175"))))
    
    
    cos_sim_m = cos_sim_matrix(mut_mat, mut_mat)
    row_clust = cos_sim_m %>% 
        dist() %>%
        hclust()
    
    cos_sim_heat_fig = Heatmap(cos_sim_m,
                               name = "Cosine \nsimilarity",
                               col = colorRamp2(seq(0,1, length.out = 9),brewer.pal(9, "YlGnBu")),
                               cluster_columns = row_clust,
                               cluster_rows = row_clust,
                               row_names_gp = gpar(fontsize = 6), 
                               column_names_gp = gpar(fontsize = 6),
                               top_annotation = ha)
    pdf(file.path(paste0(type, "_results"), "cosine_similarity_heatmap.pdf"))
    print(cos_sim_heat_fig)
    dev.off()
    
    
    cos_sim_m = cos_sim_matrix(mut_mat, signatures)
    cos_sim_heat_fig = plot_cosine_heatmap(cos_sim_m, cluster_cols = T, cluster_rows = T, plot_values = F)
    ggsave(file.path(paste0(type, "_results"), "cosine_similarity_sigs_heatmap.pdf"), cos_sim_heat_fig)
    
    if (!.is_na(signatures_exposure)){
        cos_sim_m = cos_sim_matrix(mut_mat, signatures_exposure)
        cos_sim_heat_fig = plot_cosine_heatmap(cos_sim_m, cluster_cols = T, cluster_rows = T, plot_values = F)
        ggsave(file.path(paste0(type, "_results"), "cosine_similarity_exposigs_heatmap.pdf"), cos_sim_heat_fig)
    }
    
    
    
    # Perform refitting with signatures found in nmf
    strict_refit = fit_to_signatures_strict(mut_mat, sigs_used, method = refit_method, max_delta = 0.004)
    
    # Plot refitting results
    tb = strict_refit$fit_res$contribution %>%
        as.data.frame() %>%
        tibble::rownames_to_column("Signature") %>%
        tidyr::pivot_longer(-Signature, names_to = "Sample_name", values_to = "Contribution") %>%
        left_join(meta_tb, by = "Sample_name") %>% 
        dplyr::mutate(
            BM_nr = factor(BM_nr, levels = unique(BM_nr)),
            Signature = factor(Signature, levels = unique(Signature)),
            Type = factor(Type, levels = c("Tumor", "Tumoroid", "Tumoroid2", "TumoroidL"))
        )
    
    contri_abs_fig = ggplot(tb, aes(x = Type, y = Contribution, fill = Signature)) +
        geom_bar(stat = "identity", colour = "black") +
        facet_grid(. ~ Sample, scales = "free", space = "free") +
        scale_fill_npg() +
        theme_bw() +
        theme(
            panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.spacing = unit(0, "lines"),
            axis.text.x = element_text(angle = 90)
        )
    ggsave(file.path(paste0(type, "_results"), "refit_contribution_absolute.pdf"), contri_abs_fig, width = 18)
    
    
    contri_rel_fig = ggplot(tb, aes(x = Type, y = Contribution, fill = Signature)) +
        geom_bar(position = "fill", stat = "identity", colour = "black") +
        facet_grid(. ~ Sample, scales = "free", space = "free") +
        scale_fill_npg() +
        theme_bw() +
        theme(
            panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.spacing = unit(0, "lines"),
            axis.text.x = element_text(angle = 90)
        )
    ggsave(file.path(paste0(type, "_results"), "refit_contribution_relative.pdf"), contri_rel_fig, width = 18)
    
    
    recon_fig = plot_original_vs_reconstructed(mut_mat, strict_refit$fit_res$reconstructed)
    ggsave(file.path(paste0(type, "_results"), "refit_reconstruction.pdf"), recon_fig)
    
    
    
    
    
    # Perform bootstrapped refitting
    contri_boots = fit_to_signatures_bootstrapped(mut_mat, sigs_used, n_boots = 100, method = "strict_backwards", max_delta = 0.004)
    boot_fig = plot_bootstrapped_contribution(contri_boots, plot_type = "dotplot")
    ggsave(file.path(paste0(type, "_results"), "bootstrapped_refit_contribution.pdf"), boot_fig, width = 14, height = 14)
    
    
    # Pooled by type analysis
    pooling_factor = factor(factor(meta_tb$Type, levels = c("Tumor", "Tumoroid", "Tumoroid2", "TumoroidL"))) # double factor to drop unused levels
    mut_mat_type = pool_mut_mat(mut_mat, pooling_factor)
    
    # Create profile
    if (str_detect(type, "snv")){
        type_profile_fig = plot_96_profile(mut_mat_type, condensed = T)
    } else{
        type_profile_fig = plot_indel_contexts(mut_mat_type, condensed = T)
    }
    ggsave(file.path(paste0(type, "_results"), "pooled_profile.pdf"), type_profile_fig)
    
    # Calculate if profiles are significantly different.
    
    if (ncol(mut_mat_type) == 2){
        chisq_res = broom::tidy(chisq.test(as.matrix(mut_mat_type)))
        write_tsv(chisq_res, file.path(paste0(type, "_results"), "pooled_profile_chisq.txt"))
    } else{
        chisq_res1 = broom::tidy(chisq.test(as.matrix(mut_mat_type[,c(1,2)]), simulate.p.value = TRUE, B = 10000))
        chisq_res2 = broom::tidy(chisq.test(as.matrix(mut_mat_type[,c(1,3)]), simulate.p.value = TRUE, B = 10000))
        chisq_res = rbind(chisq_res1, chisq_res2)
        chisq_res$group = colnames(mut_mat_type)[2:3]
        write_tsv(chisq_res, file.path(paste0(type, "_results"), "pooled_profile_chisq.txt"))
    }
    
    # Create similarity matrix
    cos_sim_m = cos_sim_matrix(mut_mat_type, mut_mat_type)
    cos_sim_heat_fig = plot_cosine_heatmap(cos_sim_m, cluster_cols = T, cluster_rows = T, plot_values = T)
    ggsave(file.path(paste0(type, "_results"), "pooled_cosine_similarity_heatmap.pdf"), cos_sim_heat_fig)
    
    # Refit
    strict_refit = fit_to_signatures_strict(mut_mat_type, sigs_used, max_delta = 0.004, method = refit_method)
    contri_fig = plot_contribution(strict_refit$fit_res$contribution)
    ggsave(file.path(paste0(type, "_results"), "pooled_refit.pdf"), contri_fig)
    
    recon_fig = plot_original_vs_reconstructed(mut_mat_type, strict_refit$fit_res$reconstructed)
    ggsave(file.path(paste0(type, "_results"), "pooled_refit_reconstruction.pdf"), recon_fig)
    
    invisible(0)
}


setwd("~/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Freek/Freek_Michael/2021_analysis/")


# Read meta data
full_meta_tb = read_tsv("BM_WGS_Apr2021_V3.txt") %>% 
    dplyr::mutate(Type = case_when(Type == "tumor" ~ "Tumor",
                                   Type == "organoid" ~ "Tumoroid",
                                   Type == "organoid2" ~ "Tumoroid2",
                                   Type == "organoidL" ~ "TumoroidL",
                                   Type == "normal" ~ "normal"))
age_tbl = read_tsv("220208_RMS_tumoroid_samples_age.txt")
full_meta_tb = left_join(full_meta_tb, age_tbl, by = "Sample")

# Determine names
vcf_fnames = list.files("vcf_filter_RMS/", full.names = TRUE)
sample_names = vcf_fnames %>% 
    basename() %>% 
    stringr::str_remove("_filter.vcf")

# Read and filter vcfs
grl = purrr::map2(vcf_fnames, sample_names, read_filter_vcf, min_vaf = 0.4) %>% 
    GRangesList()
names(grl) = sample_names

GenomeInfoDb::genome(grl) = 'hg38'
seqlevels(grl, pruning.mode = "coarse") = paste0("chr", c(1:22, "X", "Y"))

# Get snvs and remove samples with too few muts
grl_snv = get_mut_type(grl, "snv")
grl_snv = grl_snv[elementNROWS(grl_snv) >= 50]


# Order samples based on meta data. (This ensures they are sorted based on the patient.)
snv_meta_tb = inner_join(full_meta_tb, tibble("BM_nr" = names(grl_snv)), by = "BM_nr")
grl_snv = grl_snv[snv_meta_tb$BM_nr]
names(grl_snv) = snv_meta_tb$Sample_name

# Get signatures
signatures = get_known_signatures(genome = "GRCh38")
signatures_exposure = get_known_signatures(source = "SIGNAL", sig_type = "exposure")

# Get mut mat
ref_genome = "BSgenome.Hsapiens.UCSC.hg38"
mut_mat = mut_matrix(grl_snv, ref_genome)


# Perform nmf
sc <- scNMFSet(count = mut_mat)
estimate_bayes <- vb_factorize(sc, ranks = 1:8, nrun = 200, 
                               progress.bar = FALSE, verbose = 0)
plot(estimate_bayes)

set.seed(44)
nmf_res = extract_signatures(mut_mat, 4, nmf_type = "variational_bayes", nrun = 200)
sigs_combi = cbind(signatures, signatures_exposure)
nmf_res <- rename_nmf_signatures(nmf_res, sigs_combi, cutoff = 0.7, suffix = "")
nmf_res_sigs_fig = plot_96_profile(nmf_res$signatures)
sigs_used = sigs_combi[,colnames(sigs_combi) %in% colnames(nmf_res$signatures)]


# Determine which samples are organoid and tumoroid
tumor_comp_samples_f = snv_meta_tb$Type %in% c("Tumor", "Tumoroid")
organoid_comp_samples_f = snv_meta_tb$Sample %in% c("RMS006", "RMS007", "RMS012", "RMS335", "RMS444") & !snv_meta_tb$Type %in% c("Tumor")

# Run signature analysis
do_signature_analysis(mut_mat[,tumor_comp_samples_f, drop = FALSE], signatures, signatures_exposure, snv_meta_tb[tumor_comp_samples_f,], sigs_used, "best_subset", "snv_tumorvstumoroid")
do_signature_analysis(mut_mat[,organoid_comp_samples_f, drop = FALSE], signatures, signatures_exposure, snv_meta_tb[organoid_comp_samples_f,], sigs_used, "best_subset", "snv_tumoroidvsextratumoroid")



# Look at rainfall plot of samples with many dbs
grl_snv2 = get_mut_type(grl, "snv", predefined_dbs_mbs = TRUE)
chromosomes = paste0("chr", c(1:22, "X"))
rain_fig1 = plot_rainfall(grl_snv2[["PMOBM000AAJ"]], chromosomes)
rain_fig2 = plot_rainfall(grl_snv2[["PMRBM000AAB"]], chromosomes)
ggsave("snv_results/rainfall_PMOBM000AAJ.pdf", rain_fig1)
ggsave("snv_results/rainfall_PMRBM000AAB.pdf", rain_fig2)









###### Use a lower min vaf
grl = purrr::map2(vcf_fnames, sample_names, read_filter_vcf, min_vaf = 0.3) %>%
    GRangesList()
names(grl) = sample_names

GenomeInfoDb::genome(grl) = 'hg38'
seqlevels(grl, pruning.mode = "coarse") = paste0("chr", c(1:22, "X", "Y"))

# Get snvs and remove samples with too few muts
grl_snv = get_mut_type(grl, "snv")


# Count mutation burden
nr_muts = elementNROWS(grl_snv)

grl_snv_nonsynonymous = purrr::map(as.list(grl_snv), function(gr){gr[gr$impact == "MODERATE" | gr$impact == "HIGH"]})
nr_muts_nonysnonymous = elementNROWS(grl_snv_nonsynonymous)
nr_muts_tbl = tibble("BM_nr" = names(grl_snv), "nr_muts" = nr_muts, "nr_muts_nonsynonymous" = nr_muts_nonysnonymous)
write_tsv(nr_muts_tbl, "nr_muts_vaf0.3.txt")

nr_muts_tbl = nr_muts_tbl %>% 
    left_join(full_meta_tb, by = "BM_nr") %>% 
    dplyr::mutate(fusion = ifelse(Sample %in% c("RMS000ETY", "RMS007", "RMS012", "RMS444"), "Negative", "Positive")) %>% 
    dplyr::filter(Type != "TumoroidL" & Type != "Tumoroid2")

mean_mut_comparison_fig = ggplot(nr_muts_tbl, aes(x = fusion, y = nr_muts)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point() +
    facet_grid(. ~ Type) +
    stat_compare_means() +
    labs(y = "Nr mutations", x = "Fusion status") +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.spacing = unit(0, "lines"),
          axis.text.x = element_text(angle = 90)
    )
ggsave("mean_mut_comparison.pdf", mean_mut_comparison_fig)


grl_snv = grl_snv[elementNROWS(grl_snv) >= 50]


# Order samples based on meta data. (This ensures they are sorted based on the patient.)
snv_meta_tb = inner_join(full_meta_tb, tibble("BM_nr" = names(grl_snv)), by = "BM_nr")
grl_snv = grl_snv[snv_meta_tb$BM_nr]
names(grl_snv) = snv_meta_tb$Sample_name

# Get signatures
signatures = get_known_signatures(genome = "GRCh38")
signatures_exposure = get_known_signatures(source = "SIGNAL", sig_type = "exposure")

# Get mut mat
ref_genome = "BSgenome.Hsapiens.UCSC.hg38"
mut_mat = mut_matrix(grl_snv, ref_genome)

# Perform nmf
sc <- scNMFSet(count = mut_mat)
estimate_bayes <- vb_factorize(sc, ranks = 1:8, nrun = 200,
                               progress.bar = FALSE, verbose = 0)
plot(estimate_bayes)

set.seed(44)
nmf_res = extract_signatures(mut_mat, 4, nmf_type = "variational_bayes", nrun = 200)
sigs_combi = cbind(signatures, signatures_exposure)
nmf_res <- rename_nmf_signatures(nmf_res, sigs_combi, cutoff = 0.75, suffix = "")
nmf_res_sigs_fig = plot_96_profile(nmf_res$signatures)
## The extracted SBSA is a combi of SBS18 and SBS1. SBS18 was also extracted separately.
## Therefore both SBS1 and SBS18 are used for further refitting.
#sigs_used = sigs_combi[,c("SBS1", "SBS5", "SBS18", "Temozolomide..200.uM..1")]
sigs_used = sigs_combi[,colnames(sigs_combi) %in% colnames(nmf_res$signatures)]

# Determine which samples are organoid and tumoroid
tumor_comp_samples_f = snv_meta_tb$Type %in% c("Tumor", "Tumoroid")
organoid_comp_samples_f = snv_meta_tb$Sample %in% c("RMS006", "RMS007", "RMS012", "RMS335", "RMS444") & !snv_meta_tb$Type %in% c("Tumor")

# Run signature analysis
do_signature_analysis(mut_mat[,tumor_comp_samples_f, drop = FALSE], signatures, signatures_exposure, snv_meta_tb[tumor_comp_samples_f,], sigs_used, "best_subset", "snv_vaf0.3_tumorvstumoroid")
do_signature_analysis(mut_mat[,organoid_comp_samples_f, drop = FALSE], signatures, signatures_exposure, snv_meta_tb[organoid_comp_samples_f,], sigs_used, "best_subset", "snv_vaf0.3_tumoroidvsextratumoroid")


# Look at relation signatures with age
strict_refit = fit_to_signatures_strict(mut_mat, sigs_used, method = "best_subset", max_delta = 0.004)
contri_tbl = strict_refit$fit_res$contribution %>% 
    as.data.frame() %>% 
    rownames_to_column("Signature") %>% 
    tidyr::gather(key = "Sample_name", value = "Contribution", -Signature) %>% 
    left_join(full_meta_tb, by = "Sample_name") %>% 
    dplyr::filter(Type != "TumoroidL" & Type != "Tumoroid2")

lm_tbl = contri_tbl %>% 
    split(., paste0(contri_tbl$Signature, "_", contri_tbl$Type)) %>% 
    purrr::map(~lm(Contribution ~ age_at_sample, data = .x)) %>% 
    purrr::map(~broom::tidy(.x)) %>% 
    bind_rows(.id = "Signature_Type")
write_tsv(lm_tbl, "linear_models_age.txt")

age_fig = ggplot(contri_tbl, aes(x = age_at_sample, y = Contribution)) +
    geom_point(aes(colour = Sample)) +
    geom_smooth(method = "lm") +
    facet_grid(Type ~ Signature) +
    theme_bw() +
    theme(
        panel.spacing = unit(0, "lines"),
        axis.text.x = element_text(angle = 90)
    )
ggsave("age_fig.pdf", age_fig)


# Look at indels
# Get indels
grl_indel = get_mut_type(grl, "indel")
grl_indel = grl_indel[elementNROWS(grl_indel) >= 20]

# Order samples based on meta data. (This ensures they are sorted based on the patient.)
indel_meta_tb = inner_join(full_meta_tb, tibble("BM_nr" = names(grl_indel)), by = "BM_nr")
grl_indel = grl_indel[indel_meta_tb$BM_nr]

grl_indel = get_indel_context(grl_indel, ref_genome)
indel_counts = count_indel_contexts(grl_indel)

indel_signatures = get_known_signatures("indel")


# Perform nmf
sc <- scNMFSet(count = indel_counts)
set.seed(4)
estimate_bayes <- vb_factorize(sc, ranks = 1:8, nrun = 200, 
                               progress.bar = FALSE, verbose = 0)
plot(estimate_bayes)

set.seed(44)
nmf_res = extract_signatures(indel_counts, 3, nmf_type = "variational_bayes", nrun = 200)
nmf_res <- rename_nmf_signatures(nmf_res, indel_signatures, cutoff = 0.7, suffix = "")
nmf_res_sigs_fig = plot_indel_contexts(nmf_res$signatures)
cossim_m = cos_sim_matrix(nmf_res$signatures, indel_signatures)
plot_cosine_heatmap(cossim_m, cluster_cols = F, cluster_rows = F)

sigs_used = indel_signatures[,colnames(indel_signatures) %in% colnames(nmf_res$signatures)]

# Use all indel signatures, because there are so few that overfitting is less of a problem.
# Also the signature extraction is not working great with so few mutations.
# Also the indels seem to contain more artefacts.
do_signature_analysis(indel_counts, indel_signatures, NA, indel_meta_tb, indel_signatures, "backwards", "indel")

# Look at dbs
# All samples have 0 dbs, except for two samples which have a lot.
dbs_grl = get_mut_type(grl, "dbs")
dbs_grl = dbs_grl[elementNROWS(dbs_grl) > 10]
dbs_grl = get_dbs_context(dbs_grl)
dbs_counts = count_dbs_contexts(dbs_grl)
dbs_profile_fig = plot_dbs_contexts(dbs_counts)
ggsave(file.path("dbs_results", "dbs_profile.pdf"), dbs_profile_fig)


# Look at the profile between different VAFs.
grl = purrr::map2(vcf_fnames[1], sample_names[1], read_filter_vcf, min_vaf = 0.0000001) %>% 
    GRangesList()
names(grl) = sample_names[1]

gr = grl[[1]]
bins = cut(gr$vaf, breaks = c(0, 0.1, 0.2, 0.3, 0.4, 1))
grl_vaf = split(gr, bins)
GenomeInfoDb::genome(grl_vaf) = 'hg38'
seqlevels(grl_vaf, pruning.mode = "coarse") = paste0("chr", c(1:22, "X", "Y"))
grl_vaf = get_mut_type(grl_vaf, "snv")
mut_mat = mut_matrix(grl_vaf, ref_genome)
vaf_profile_fig = plot_96_profile(mut_mat, condensed = TRUE)
ggsave("vaf_profile.pdf", vaf_profile_fig)
