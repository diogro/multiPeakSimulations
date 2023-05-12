if(!require(optparse)){install.packages("optparse"); library(optparse)}
if(!require(Hmisc)){install.packages("Hmisc"); library(Hmisc)}
if(!require(doMC)){install.packages("doMC"); library(doMC)}

if (sys.nframe() == 0L) {
  option_list <- list(
    make_option("--label", 
                help = ("Integrated matrix label"),
                metavar = "label"),
    make_option("--type", 
                 help = ("Type of landscape:random or enriched."),
                 metavar = "type"),
    make_option("--diag", default = "low",
                help = ("Type of diagonal matrix: [diag] or low"),
                metavar = "diag"),
    make_option("--min_dist", default = 3, 
                 help = ("Minimum distance of peaks from origin. Default is 3."),
                 metavar = "min_dist"),
    make_option("--max_dist", default = 10,
                 help = ("Maximum distance of peaks from origin. Default is 10."),
                 metavar = "max_dist"),
    make_option("--max_gens", default = 20000, 
                 help = ("Maximum number of generations. Default is 20k."),
                 metavar = "max_gens"),
    make_option("--diff_cut_off", default = 1e-5, 
                 help = ("Cuf_off for considering population has not moved. Default is 1e-5."),
                 metavar = "diff_cut_off"),
    make_option("--max_stand_still", default = 1000,
                 help = ("Maximum number of generations for population with movement bellow diff_cut_off. Default is 1k."),
                 metavar = "max_stand_still"),
    make_option("--n_sims", default = 1000,
                 help = ("Number of simulations in each scenario. Default is 1k."),
                 metavar = "n_sims"),
    make_option("--n_peaks", default = 50,
                help = ("Number of peaks in the multi peak landscape. Default is 50."),
                metavar = "n_peaks"),
    make_option("--trim_trajectory", default = TRUE,
                help = ("Trim trajectory? If FALSE, will produce very large output. Default is TRUE."),
                metavar = "trim_trajectory"),
    make_option("--random_seed", default = 42,
                help = ("Random seed. Default is 42."),
                metavar = "random_seed"),
    make_option("--n_cores", default = min(detectCores()-1, 64),
                help = ("Number of cores to use in computation. Default is up to 64."),
                metavar = "n_cores")
  )
  parser_object <- OptionParser(usage = "Rscript %prog --label MATRIX_LABEL --type [random|enriched] [aditional options]\n",
                                option_list = option_list,
                                description = "Run multi peak simulations.")
  
  ## aliases
  opt <- parse_args(parser_object, args = commandArgs(trailingOnly = TRUE), positional_arguments = TRUE)
  label <- opt$options$label
  type <- opt$options$type
  diag_mat_type <- opt$options$diag
  min_dist = opt$options$min_dist
  space_size = opt$options$max_dist
  max_gens = opt$options$max_gens
  diff_cut_off = opt$options$diff_cut_off
  max_stand_still = opt$options$max_stand_still
  n_sims = opt$options$n_sims
  n_peaks = opt$options$n_peaks
  trim_trajectory = opt$options$trim_trajectory
  random_seed = opt$options$random_seed
  n_cores = opt$options$n_cores
}
#type = "random"
#n_sims = 10
#n_peaks = 2
#n_cores = 1

set.seed(random_seed)
if(n_cores < 1) n_cores = 1
registerDoMC(cores=n_cores)

source("./prepareHighDim.R")

if(is.null(type)) stop("Please specify landscape type.")
if(type == "random"){
    peakPool_diag = peakPool_corr = peakPool_random
} else if(type == "enriched"){
    peakPool_diag = peakPool_G_diag_enriched
    peakPool_corr = peakPool_G_corr_enriched
} else 
    stop("Unknown landscape type. Please use random or enriched.")

say(paste("Running ", n_sims, " simulations for ", type," landscape with ", n_peaks, " peaks, using ", n_cores, " cores."), by = "spider")

# Test runs
#########################

#x = runSimulation("Integrated", G_corr, n_peaks = 1, n_traits, scale = 40, peakPool = peakPool_G_corr_enriched)
#sapply(x, object.size)
#runSimulation("Diagonal", G_diag, n_peaks = 1, n_traits, scale = 40, peakPool = peakPool_G_diag_enriched)
#runSimulation("Integrated", G_corr, n_peaks = 50, n_traits, scale = 40, peakPool = peakPool_random)
#runSimulation("Diagonal", G_diag, n_peaks = 50, n_traits, scale = 40, peakPool = peakPool_random)

# Simulations
#########################

results = runTrypitch(G_diag, peakPool_diag,          
                      G_corr, peakPool_corr,          
                      n = n_sims, n_peaks = n_peaks, scale = 40)

output.dir <- "./output"
if (!file.exists(file.path(output.dir, "plots")))
  dir.create(file.path(output.dir, "plots"), showWarnings = TRUE, recursive = TRUE)
if (!file.exists(file.path(output.dir, "Rds")))
  dir.create(file.path(output.dir, "Rds"), showWarnings = TRUE, recursive = TRUE)

output_name = paste0(label, "-", type,
                    "-diag_mat_type-", diag_mat_type,
                    "_minDist-", min_dist,
                    "_dcoff-", diff_cut_off,
                    "_maxGens-", max_gens,
                    "_nPeaks-", n_peaks,
                    "_nSims-", n_sims,
                    "_seed-", random_seed)
saveRDS(results, file = file.path(output.dir, "Rds", paste0(output_name, ".Rds")))
 
plots = plot_grid(plotDzgmax_normdz(results$DS, ylim = c(2, space_size), main = "Diagonal G - Single Peak"),
                  plotDzgmax_normdz(results$CS, ylim = c(2, space_size), main = "Integrated G - Single Peak"),
                  plotDzgmax_normdz(results$DM, ylim = c(2, space_size), main = "Diagonal G - Multiple Peaks"),
                  plotDzgmax_normdz(results$CM, ylim = c(2, space_size), main = "Integrated G - Multiple Peaks"),
                  ncol = 2, labels = LETTERS[1:4])
save_plot(file.path(output.dir, "plots", paste0(output_name, ".png")), plots, base_height = 5, base_asp = 1.3, ncol = 2, nrow = 2)

