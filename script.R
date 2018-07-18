# Metagene Projection Methodology R script -- Leukemia 1 Example
#
# Pablo Tamayo 3/29/2007 (tamayo@broad.mit.edu)
#
# This script runs the first Leukemia example (ALL/AML) from the paper:
#
# Metagene projection for cross-platform, cross-species characterization of global transcriptional states
# P. Tamayo, D. Scanfeld, B. L. Ebert, M. A. Gillette, C. W. M. Roberts, and J.P. Mesirov  
# Proc. Natl. Acad. Sci. USA, 104: 5959-5964 2007. http://www.pnas.org/cgi/content/abstract/0701068104v1
#
# It uses the main function "MetaGene.Projection(...)" which implement most of the method described in the
# paper (except for the "k" number-of-metagenes selection).
# The MP library is loaded (sourced) with other functions from library: MP.Library.R"
# The GSEA analysis of the resulting metagenes and analysis of different choices for k (model selection)
# are implemented using separate functions and scripts.
#
# To run: cut and paste (or source) this code inside the R GUI console. The plots will be produced in
# the R GUI screen and also saved in files. Before running on new datasets try to reproduce the Leukemia1 example
# from the paper. Then to run on different datasets e.g. modify the pathnames and parameters below accordingly.
# This script takes about 1:40 minutes of CPU to run using R 2.2.1 on a Windows-XP Dell Inspiron 630m laptop.
#
# While running this script calls "MetaGene.Projection" (see below) which will produce the following output files
# under the directory specified by "output.dir":
#
# Main output files:
# <identifier>.<date>_<time>.log.txt = File containing the parameters used in the run and the data ans time
# <identifier>.model_dataset.H.gct = projected model dataset
# <identifier>.all.H.cls =  projection of model + test datasets (cls phenotypes)
# <identifier>.all.H.gct = projection of model + test datasets (gct dataset)
# <identifier>.heatmap.jpeg = heat map of projection
# <identifier>.heatmap.sorted.jpeg = heat map of projection sorted inside each phenotype
# <identifier>.2D.proj.jpeg = 2D biplot of projected model and test datasets
# <identifier>.H.htree.jpeg = hierarchical tree built on the projected model and test datasets
# <identifier>.pred.gct = projection-based SVM prediction results (gct dataset)
# <identifier>.pred.txt = projection-based SVM prediction results (text file)
# <identifier>.H.mem.txt = clustering membership based on metagene with largest amplitude
#
# Other complementary output files:
# <identifier>.model.H.gct = H matrix from the NMF decomposition
# <identifier>.model.W.gct = W matrix from the NMF decomposition
# <identifier>.model_set.2.cls = model dataset after pre-preprocessing and refinement (cls phenotypes)
# <identifier>.model_set.2.gct = model dataset after pre-preprocessing and refinement (gct file)
# <identifier>.model_set.1.cls = model dataset after pre-preprocessing and before refinement (cls phenotypes)
# <identifier>.model_set.1.gct = model dataset after pre-preprocessing and before refinement (gct files)
# <identifier>.model_set.0.cls = model dataset before pre-preprocessing but containing samples after refinement (cls phenotypes)
# <identifier>.model_set.0.gct = model dataset before pre-preprocessing but containing samples after refinement (gct file)
# <identifier>.htree.jpeg = hierarchical tree on original pre-projection dataset
# <identifier>.all.cls = consolidated model + test dataset in the space of common genes (cls phenotypes)
# <identifier>.all.gct = consolidated model + test dataset in the space of common genes (gct dataset)
# <identifier>.prelim.pred.txt = preliminary projection-based SVM prediction results (used in refinement) (text file)

MP.library.location  <-  "/Users/aaronabraham/Documents/HYRS/mouse_project/script/MP.Library.R"
source(MP.library.location, verbose=T, max.deparse.length=9999)   # Load Metagene Projection library 

# Define model & test datasets and parameters

model.dataset.table <-  # Defines the input model dataset and pre-processing options for it
  list( # Subset of samples from Ross et al 2003 (PMID: 12730115) and Ross et al 2004 (PMID: 15226186)
    gct.file = "/Users/aaronabraham/Documents/HYRS/mouse_project/data/model_data.gct",  # Gene expression data
    cls.file = "/Users/aaronabraham/Documents/HYRS/mouse_project/data/model_classes.cls",  # Annotation of classes
    column.subset = "ALL",        # Which subset (samples or phenotypes) to include (default: "ALL" : all of them)
    column.sel.type = "samples",  # Selection type: "sample": or "phenotypes"
    thres = "NULL",                   # Threshold to apply to dataset before projection
    ceil = "NULL",                # Ceiling to apply to dataset before projection
    fold = 2,                     # Fold change (max/min) for variation filter before projection
    delta = "NULL",                  # Absolute difference (max - min) for variation filter before projection
    norm = 6)                     # Normalization before projection (default 6 column-rank and rescaling normalization)

test.datasets.table <- # Defines one or more input test datasets and pre-processing options for each of them
  list( 
    list( # Using same model dataset as the test dataset
      gct.file = "/Users/aaronabraham/Documents/HYRS/mouse_project/data/test_data.gct", # Test gene expression data
      cls.file = "/Users/aaronabraham/Documents/HYRS/mouse_project/data/test_classes.cls",  # Test class annotation
      column.subset = "ALL",       # Which subset (sample numbers( to include (default: "ALL" : all of them)
      column.sel.type = "samples", # Selection type: "sample": or "phenotypes"
      thres = "NULL",                  # Threshold to apply to dataset before projection
      ceil = "NULL",               # Ceiling to apply to dataset before projection
      norm = 6)                    # Normalization before projection (default 6 column-rank and rescaling normalization)
  )

# Define parameters for this specific run (see detailed definitions below)

identifier           <-    "mouse_model1"   
k.proj               <-    2     #K value for NMF
alg                  <-    "NMF.div"
niter                <-    2000
seed                 <-    1234
nchar.phen           <-    9
postprojnorm         <-    TRUE
heatmap.row.norm     <-    FALSE
heatmap.cmap.type    <-    6
output.dir           <-    "/Users/aaronabraham/Documents/HYRS/mouse_project/output/"
high.conf.thres      <-    0.3
kernel               <-    "radial"
cost                 <-    1
gamma                <-    5
model.set.refinement <-    TRUE #Can be true, but for test run it will be false
symbol.scaling       <-    0.55

# These are the symbols and colors to use for each phenotype in the model and test sets 
#          model samples:   square symbols
#                  color         symbol      phenotype
legend.list <- c("steelblue2",    21,        # D5 
                 "red",           21,        # D0
                 #          test samples:    cicle symbols
                 #                  color         symbol      phenotype                 
                 "steelblue2",    22,        # D5 
                 "red",           22        # D0
)

col <- legend.list[seq(1, length(legend.list), 2)]
symbs <- as.numeric(legend.list[seq(2, length(legend.list), 2)])

# This is the call to the Metagene Projection function:

MetaGene.Projection(                           # Runs the entire methodology
  #   (except for the GSEA analysis of metagenes and the model (k) selection)
  model.dataset.table = model.dataset.table,   # R list with model dataset parameters (see model.dataset.table above)
  test.datasets.table = test.datasets.table,   # R list with test dataset(s) parameters (see model.dataset.table above)
  identifier = identifier,                     # Prefix to name output files
  k.proj = k.proj,                             # Number of metagenes in projection
  alg = alg,                                   # Algorithm for Metagene Projection (default NMF.div):
  #    "NMF.div" : Non-Negative Matrix Factorization using the divergence cost
  #  (other algorithms for projection are internally supported but have not being tested)
  niter = niter,                               # Number of algorithm iterations (default: 2000)
  seed = seed,                                 # Random seed to initialize metagene matrices (default: 1234)
  nchar.phen =  nchar.phen,                    # Number of characters to use to identify classes from the CLS files
  postprojnorm = postprojnorm,                 # TRUE or FALSE: apply post-projection normalization (i.e. scale points to unit
  #     hypersphere, default: T)
  heatmap.row.norm = heatmap.row.norm,         # TRUE or FALSE: row-normalize (standardize) the rows in the heat map (default F)
  heatmap.cmap.type = heatmap.cmap.type,       # 1 = vintage pinkogram, 2 = scale of grays, 4 = high-resolution pinkogram,
  #   6 = redish color map for metagene factors (default: 6)
  high.conf.thres = high.conf.thres,           # Confidence threshold (Brier score) to separate call from no-calls (default 0.3)
  output.dir = output.dir,                     # Output directory where the resulting output files will be produced
  col = col,                                   # Colors for the legend symbols for phenotypes: first model and then test dataset(s)
  symbs = symbs,                               # Plotting symbols for phenotypes: first model and then test dataset(s)
  symbol.scaling = symbol.scaling,             # Graphical scaling for symbols in plots and plot legends (default: 1)
  kernel = kernel,                             # Kernel function for SVM: "radial" or "linear" (default: "radial")
  cost = cost,                                 # Cost parameter for SVM (default: 1)
  gamma = gamma,                               # Gamma coefficient for radial base function kernel (default:  0.05 )
  model.set.refinement = model.set.refinement) # TRUE or FALSE: perform model set refinement (default: T)

# end of metagene projection script example
