library(ape)
tre = read.tree("/Users/Xin/Documents/Ratsnake_project/manuscript/plylo+div/tree_oldcal.txt")
tips = c("Coelognathus_helena", "Coelognathis_subradiatus", "Coelognathus_erythrurus", "Coelognathus_flavolineatus", "Coelognathus_radiatus",
                "Gonyosoma_oxycephalum", "Nerodia_sipedon", "Heterodon_platirhinos", "Gonyosoma_frenatum", "Gonyosoma_prasinum", "Gonyosoma_boulengeri", 
                  "Coluber_constrictor", "Drymobius_margaritiferus", "Gyalopion_canum", "Hapsidophrys_lineatus", "Hemorrhois_ravergieri","Tantilla_coronata")
       
tre = drop.tip(tre, tips)
write.tree(tre, file = "/Users/Xin/Documents/Ratsnake_project/manuscript/plylo+div/tree_oldcal3.txt")
# Load the package (after installation, see above).
library(optimx)         # You need to have some version of optimx available
# as it is a BioGeoBEARS dependency; however, if you
# don't want to use optimx, and use optim() (from R core) 
# you can set:
# BioGeoBEARS_run_object$use_optimx = FALSE
# ...everything should work either way -- NJM 2014-01-08
library(FD)       # for FD::maxent() (make sure this is up-to-date)
library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
library(parallel)
library(BioGeoBEARS)
library(BioGeoBEARS)
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_basics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_classes_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_models_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_plots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_readwrite_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_simulate_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stratified_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_univ_model_v1.R")
source("http://phylo.wikidot.com/local--files/biogeobears/calc_loglike_sp_v01.R")
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)
# slight speedup hopefully
# You will need to set your working directory to match your local system
# You can find the input files at:
path="/Users/Xin/Documents/Ratsnake_project/manuscript/plylo+div/"  
temp_dir = np(path)
temp_dir
list.files(temp_dir)

# Note these very handy functions!
# Command "setwd(x)" sets your working directory
# Command "getwd()" gets your working directory and tells you what it is.
# Command "list.files()" lists the files in your working directory
# To get help on any command, use "?".  E.g., "?list.files"

# Set your working directory for output files
# default here is your home directory ("~")
# Change this as you like
#wd = np("~")
#setwd(wd)

# Double-check your working directory with getwd()
getwd()
# This is the example Newick file for Hawaiian Psychotria
# (from Ree & Smith 2008)
# "trfn" = "tree file name"
trfn = np(paste(addslash(path), "tree_oldcal2.txt", sep=""))
#trfn = np(paste(addslash(path), "tree_oldcal3.txt", sep=""))

# Look at the raw Newick file:
moref(trfn)

# Look at your phylogeny:
tr = read.tree(trfn)
tr<-ladderize(tr)
plot(tr)
#title("Example Psychotria phylogeny from Ree & Smith (2008)")
axisPhylo() # plots timescale
# This is the example geography file for Hawaiian Psychotria
# (from Ree & Smith 2008)
geogfn = np(paste(addslash(path), "ratgeo.txt", sep=""))
#geogfn = np(paste(addslash(path), "ratgeo2.txt", sep=""))
# Look at the raw geography text file:
moref(geogfn)

# Look at your geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

# Set the maximum number of areas any species may occupy; this cannot be larger 
# than the number of areas you set up, but it can be smaller.
max_range_size = 3
# DEC AND DEC+J ANALYSIS
#######################################################
####################################################
####################################################
# KEY HINT: The number of states (= number of different possible geographic ranges)
# depends on (a) the number of areas and (b) max_range_size.
# If you have more than about 500-600 states, the calculations will get REALLY slow,
# since the program has to exponentiate a matrix of e.g. 600x600.  Often the computer
# will just sit there and crunch, and never get through the calculation of the first
# likelihood.
# 
# (this is also what is usually happening when LAGRANGE hangs: you have too many states!)
#
 
# Set up empty tables to hold the statistical results
restable = NULL
teststable = NULL

#######################################################
# Run DEC
#######################################################

BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$force_sparse=FALSE    # sparse=FALSE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object$speedup=TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$calc_ancprobs=TRUE    # get ancestral states from optim run

# Set up a time-stratified analysis 
# (un-comment to use; see example files in extdata_dir, 
#  and BioGeoBEARS google group posts for further hints)
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
read_dispersal_multipliers_fn(inputs = NULL,dispersal_multipliers_fn = "manual_dispersal_multipliers.txt")
dispfn = np(paste(addslash(path), "manual_dispersal_multipliers.txt", sep=""))
moref(dispfn)

BioGeoBEARS_run_object$dispersal_multipliers_fn = dispfn

#read_areas_allowed_fn(inputs = NULL,areas_allowed_fn = "area_allowed.txt")
#areafn = np(paste(addslash(path), "area_allowed.txt", sep=""))
#moref(areafn)
#BioGeoBEARS_run_object$areas_allowed = areafn
 
#BioGeoBEARS_run_object$areas_allowed_fn = areafn
# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

# Multicore processing if desired
BioGeoBEARS_run_object$num_cores_to_use=1
# (use more cores to speed it up; this requires
# library(parallel), which is default on Macs in
# R 3.0+, but apparently still has to be typed
# on Windows machines. Note: apparently parallel
# works on Mac command-line R, but not R.app.
# BioGeoBEARS checks for this and resets to 1
# core with R.app)

# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
# I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
# but the results are imprecise and so I haven't explored it further.
# In a Bayesian analysis, it might work OK, but the ML point estimates are
# not identical.
# Also, I have not implemented all functions to work with force_sparse=TRUE.
# Volunteers are welcome to work on it!!
#BioGeoBEARS_run_object$force_sparse=FALSE

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object

# Read in the input files. Read the error messages if you get them!
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
# Divide the tree up by strata (stratified analysis)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$master_table

# Set up DEC model
# (nothing to do; defaults)

# Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
BioGeoBEARS_run_object

# This contains the model object
BioGeoBEARS_run_object$BioGeoBEARS_model_object

# This table contains the parameters of the model 
para<-BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
para[1,3:4]<-c(1e-6,1)
para[2,3:4]<-c(1e-6,1)
#para[7,1]<-c("free")
 
#para
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table<-para
# Run this to check inputs. Read the error messages if you get them!
check_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object
# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = "Ratgeo_DEC.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDEC = res
} else {
  # Loads to "res"
  load(resfn)
  resDEC = res
}

#######################################################
# Run DEC+J
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size

# Set up the stratified part
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"

BioGeoBEARS_run_object$use_optimx = TRUE
BioGeoBEARS_run_object$num_cores_to_use=1
BioGeoBEARS_run_object$force_sparse=FALSE    # sparse=FALSE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object$speedup=TRUE        # seems to work OK
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
BioGeoBEARS_run_object$dispersal_multipliers_fn = dispfn
#BioGeoBEARS_run_object$areas_allowed = areafn
# Run to check inputs
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
# Divide the tree up by strata
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$master_table

# Set up DEC+J model
# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 1e-6
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 1
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 1e-6
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 1
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "ratgeo_DEC+J.Rdata"
runslow = TRUE
if (runslow)
{
  #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
  
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDECj = res
} else {
  # Loads to "res"
  load(resfn)
  resDECj = res
}

#######################################################
# PDF plots
#######################################################
pdffn = "ratgeo2.pdf"
pdf(pdffn, width=15, height=15)

#######################################################
# Plot ancestral states - DEC
#######################################################
analysis_titletxt ="BioGeoBEARS DEC on ratsnake"

# Setup
results_object = resDEC
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=1.2, tipcex=1, statecex=0.7, splitcex=0.7, titlecex=1, plotsplits=F, cornercoords_loc=scriptdir,plotlegend=T,include_null_range=FALSE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=1.2, tipcex=1, statecex=0.7, splitcex=0.7, titlecex=1, plotsplits=F, cornercoords_loc=scriptdir,plotlegend=T, include_null_range=FALSE, tr=tr, tipranges=tipranges)

#######################################################
# Plot ancestral states - DECJ
#######################################################
analysis_titletxt ="BioGeoBEARS DEC+J on ratsnake "

# Setup
results_object = resDECj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=1.2, tipcex=1, statecex=0.7, splitcex=0.7, titlecex=1, plotsplits=F, cornercoords_loc=scriptdir, include_null_range=FALSE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=1.2, tipcex=1, statecex=0.7, splitcex=0.7, titlecex=1, plotsplits=F, cornercoords_loc=scriptdir, include_null_range=FALSE, tr=tr, tipranges=tipranges)

dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

#######################################################
# Statistics
#######################################################
# We have to extract the log-likelihood differently, depending on the 
# version of optim/optimx
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDEC)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECj)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

#alt null            LnLalt           LnLnull DFalt DFnull DF            Dstatistic
#1          -74.7713706143037 -74.7709418033916     3      2  1 -0.000857621824252419
#pval        test       tail     AIC1     AIC2    AICwt1    AICwt2
#1    1 chi-squared one-tailed 155.5427 153.5419 0.2688571 0.7311429
#AICweight_ratio_model1 AICweight_ratio_model2
#1              0.3677217               2.719448
#############model 2, which is DEC model is better than DEC+J for ratsnake tree.

res2    # DEC, null model for Likelihood Ratio Test (LRT)
res1    # DEC+J, alternative model for Likelihood Ratio Test (LRT)

# The null hypothesis for a Likelihood Ratio Test (LRT) is that two models
# confer the same likelihood on the data. See: Brian O'Meara's webpage:
# http://www.brianomeara.info/tutorials/aic
# ...for an intro to LRT, AIC, and AICc

rbind(res2, res1)
tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)
teststable = rbind(teststable, tmp_tests)