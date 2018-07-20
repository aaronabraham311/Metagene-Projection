########################################################################################################
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2007 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.
########################################################################################################
# Library of R functions implementing:
#
# Metagene projection for cross-platform, cross-species characterization of global transcriptional states
# P. Tamayo, D. Scanfeld, B. L. Ebert, M. A. Gillette, C. W. M. Roberts, and J.P. Mesirov  
# Proc. Natl. Acad. Sci. USA, 104: 5959-5964 2007. http://www.pnas.org/cgi/content/abstract/0701068104v1
#
# Author: Pablo Tamayo  -  April 12, 2007
# The function "MetaGene.Projection" below implement most of the method described in the
# paper.
# The GSEA analysis of the resulting metagenes and analysis of different choices for k (model selection)
# are implemented using separate functions and scripts.
#
# While running "MetaGene.Projection" will produce the following output files:
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

MetaGene.Projection <- function(
  model.dataset.table,         # Defines the input model dataset and pre-processing options for it
  test.datasets.table,         # Defines one or more input test datasets and pre-processing options for each of them
  identifier,                  # Prefix to name output files
  k.proj,                      # Number of metagenes in projection
  alg = "NMF.div",             # Algorithm for Metagene Projection (default NMF.div):
  #    "NMF.div" : Non-Negative Matrix Factorization using the divergence cost
  #  (other algorithms for projection are internally supported but have not being tested)
  niter = 2000,                # Number of algorithm iterations (default: 2000)
  seed = 1234,                 # Random seed to initialize metagene matrices (default: 1234)
  nchar.phen = 9,              # Number of characters to use to identify classes from the CLS files
  postprojnorm = TRUE,         # TRUE or FALSE: apply post-projection normalization (i.e. scale points to unit
  #     hypersphere, default: T)
  heatmap.row.norm = FALSE,    # TRUE or FALSE: row-normalize (standardize) the rows in the heat map (default F)
  heatmap.cmap.type = 6,       # 1 = vintage pinkogram, 2 = scale of grays, 4 = high-resolution pinkogram,
  #   6 = redish color map for metagene factors (default: 6)
  high.conf.thres = 0.3,       # Confidence threshold (Brier score) to separate call from no-calls (default 0.3)
  output.dir,                  # Output directory where the resulting output files will be produced
  col = c("green", "blue", "pink", "red", "orange", "red4", "steelblue2", "violet"), # Colors for the legend symbols for phenotypes: first model and then test dataset(s)  
  symbs = c(22, 21, 20, 23, 24, 25, 21, 20, 23, 24, 25), # Plotting symbols for phenotypes: first model and then test dataset(s)
  symbol.scaling = 1,          # Graphical scaling for symbols in plots and plot legends (default: 1)
  kernel = "radial",           # Kernel function for SVM: "radial" or "linear" (default: "radial")
  cost = 1,                    # Cost parameter for SVM (default: 1)
  gamma = 0.05,                # Gamma coefficient for radial base function kernel (default:  0.05 )
  theta = 0,                   # Smoothing parameter only used for the NSNMF.div algorithm (default 0)
  lambda = 1,                  # Sparse parameter used only for the SNMF algorithm (default 0)
  model.set.refinement = T) {  # TRUE or FALSE: perform model set refinement (default: T)
  
  print(c(model.dataset.table))
  print(c(test.datasets.table))
  
  date.string <- date()
  date.string2 <- paste(unlist(strsplit(date.string, " ")), sep="_", collapse="_")
  date.string3 <- paste(unlist(strsplit(date.string2, ":")), sep="_", collapse="_")
  params.file <- paste(output.dir, identifier, ".", date.string3, ".params.txt", sep="")
  
  write(c("Metagene Projection Run on ", date.string2), file= params.file, ncolumns=100, append=F)
  write(c("  "), file= params.file, ncolumns=100, append=T)
  
  write("Model dataset: ", file= params.file, ncolumns=100, append=T)
  model.frame <- as.matrix(cbind(model.dataset.table))
  write.table(model.frame, quote = F, file= params.file, append=T)
  write("  ", file= params.file, ncolumns=100, append=T)
  for (i in 1:length(test.datasets.table)) {
    write(c("Test dataset:", i), file= params.file, ncolumns=100, append=T)
    test.frame <- as.matrix(cbind(test.datasets.table[[i]]))
    write.table(test.frame, quote = F, file= params.file, append=T)
    write("   ", file= params.file, ncolumns=100, append=T)
  }
  
  write("Parameters:", file= params.file, ncolumns=100, append=T)
  write(c("identifier=",identifier), file=params.file, ncolumns=100, append=T)
  write(c("k.proj=",k.proj), file= params.file, ncolumns=100, append=T)
  write(c("alg=",alg), file= params.file, ncolumns=100, append=T)
  write(c("niter=",niter), file= params.file, ncolumns=100, append=T)
  write(c("seed=",seed), file= params.file, ncolumns=100, append=T)
  write(c("nchar.phen=",nchar.phen), file= params.file, ncolumns=100, append=T)
  write(c("postprojnorm=",postprojnorm), file= params.file, ncolumns=100, append=T)
  write(c("heatmap.row.norm=",heatmap.row.norm), file= params.file, ncolumns=100, append=T)
  write(c("heatmap.cmap.type=",heatmap.cmap.type), file= params.file, ncolumns=100, append=T)
  write(c("high.conf.thres=",high.conf.thres), file= params.file, ncolumns=100, append=T)
  write(c("output.dir=",output.dir), file= params.file, ncolumns=100, append=T)
  write(c("col=",col), file= params.file, ncolumns=100, append=T)
  write(c("symbs=",symbs), file= params.file, ncolumns=100, append=T)  
  write(c("symbol.scaling=",symbol.scaling), file= params.file, ncolumns=100, append=T)  
  write(c("kernel=",kernel), file= params.file, ncolumns=100, append=T)  
  write(c("cost=",cost), file= params.file, ncolumns=100, append=T)  
  write(c("gamma=",gamma), file= params.file, ncolumns=100, append=T)  
  write(c("theta=",theta), file= params.file, ncolumns=100, append=T)  
  write(c("lambda=",lambda), file= params.file, ncolumns=100, append=T)  
  write(c("model.set.refinement=",model.set.refinement), file= params.file, ncolumns=100, append=T)  
  
  # -------------------------------------------------- Projection Methodology
  
  set.seed(seed=seed, kind = NULL)
  
  O <- MP.Subset.Dataset(
    input.ds =         model.dataset.table$gct.file,
    input.cls =        model.dataset.table$cls.file,
    column.subset =    model.dataset.table$column.subset,
    column.sel.type =  model.dataset.table$column.sel.type,
    row.subset =       "ALL", 
    output.ds =        paste(output.dir, "temp1.gct", sep=""),
    output.cls =       paste(output.dir, identifier, ".model_set.1.cls", sep=""))
  
  O <- MP.Preprocess.Dataset(
    input.ds =       paste(output.dir, "temp1.gct", sep=""),
    output.ds =      paste(output.dir, identifier, ".model_set.1.gct", sep=""),          
    thres =          model.dataset.table$thres,
    ceil =           model.dataset.table$ceil,
    fold =           model.dataset.table$fold,
    delta =          model.dataset.table$delta,
    normalization =  model.dataset.table$norm) 
  
  O <- MP.Extract.Factors(
    input.ds =        paste(output.dir, identifier, ".model_set.1.gct", sep=""),          
    input.cls =       paste(output.dir, identifier, ".model_set.1.cls", sep=""),
    output.W.file =   paste(output.dir, identifier, ".model.W.gct", sep=""),
    output.H.file =   paste(output.dir, identifier, ".model.H.gct", sep=""),
    k.proj =          k.proj,
    alg =             alg,
    niter =           niter,
    seed =            seed,
    theta =           theta,
    lambda =          lambda,
    sort.factors =    T) 
  
  O <- MP.Factors.Project(
    input.ds =          paste(output.dir, identifier, ".model_set.1.gct", sep=""),          
    factors.ds =        paste(output.dir, identifier, ".model.W.gct", sep=""),
    postprojnorm =      postprojnorm,
    output.file =       paste(output.dir, identifier, ".model_dataset.H.gct", sep=""))
  
  CLS <- MP.ReadClsFile(file =  paste(output.dir, identifier, ".model_set.1.cls", sep=""))
  model.size <- length(CLS$class.v)
  
  O <- MP.Evaluate.Projection(
    input.ds =                    paste(output.dir, identifier, ".model_dataset.H.gct", sep=""),
    input.cls =                   paste(output.dir, identifier, ".model_set.1.cls", sep=""),
    model.set =                   seq(1, model.size),
    prediction.results.file =     paste(output.dir, identifier, ".prelim.pred.txt", sep=""),
    prediction.matrix.file =      paste(output.dir, identifier, ".prelim.pred.gct", sep=""),
    col =                         col,
    use.feature.names =           use.feature.names,
    nchar.phen =                  nchar.phen,
    high.conf.thres =             high.conf.thres,
    symbs          =              symbs,
    symbol.scaling =              symbol.scaling,
    levels =                      levels,
    nlevels =                     nlevels,
    kernel =                      kernel,
    cost =                        cost,
    gamma =                       gamma)
  
  input.txt <- paste(output.dir, identifier, ".prelim.pred.txt", sep="")
  pred.table <- read.table(file=input.txt, skip=2, nrow=model.size, sep="\t", header=T, comment.char="", as.is=T)
  conf.list <- pred.table["Conf..H.L."]
  actual.list <- pred.table["Actual"]
  predicted <- pred.table["Predicted"]
  sample.select <- ((conf.list == " H ") & (actual.list == predicted))
  
  if (model.set.refinement == T) {
    high.conf.set <- seq(1, model.size)[sample.select]
  } else {
    high.conf.set <- seq(1, model.size)
  }
  
  print(c("pred table from file=", pred.table))
  print(c("model sizet=", model.size))
  print(c("high.conf.set=", high.conf.set))
  print(c("sample.select", sample.select))
  print(paste("Original: ", model.size, " (samples); New: ", length(high.conf.set), " (samples); diff: ", model.size - length(high.conf.set), sep= " "))
  
  O <- MP.Subset.Dataset(
    input.ds =         paste(output.dir, identifier, ".model_set.1.gct", sep=""),
    input.cls =        paste(output.dir, identifier, ".model_set.1.cls", sep=""),
    column.subset =    high.conf.set,
    row.subset =       "ALL", 
    output.ds =        paste(output.dir, identifier, ".model_set.2.gct", sep=""),
    output.cls =       paste(output.dir, identifier, ".model_set.2.cls", sep=""))
  
  O <- MP.Subset.Dataset(
    input.ds =         model.dataset.table$gct.file,
    input.cls =        model.dataset.table$cls.file,
    column.subset =    high.conf.set,
    row.subset =       "ALL", 
    output.ds =        paste(output.dir, identifier, ".model_set.0.gct", sep=""),
    output.cls =       paste(output.dir, identifier, ".model_set.0.cls", sep=""))
  
  O <- MP.Extract.Factors(
    input.ds =          paste(output.dir, identifier, ".model_set.2.gct", sep=""),          
    input.cls =         paste(output.dir, identifier, ".model_set.2.cls", sep=""),
    output.W.file =     paste(output.dir, identifier, ".model.W.gct", sep=""),
    output.H.file =     paste(output.dir, identifier, ".model.H.gct", sep=""),
    k.proj =            k.proj,
    alg =               alg,
    niter =             niter,
    seed =              seed,
    theta =             theta,
    lambda =            lambda,
    sort.factors =    T) 
  
  O <- MP.Factors.Project(
    input.ds =          paste(output.dir, identifier, ".model_set.2.gct", sep=""),          
    factors.ds =        paste(output.dir, identifier, ".model.W.gct", sep=""),
    postprojnorm =      postprojnorm,
    output.file =       paste(output.dir, identifier, ".model_dataset.H.gct", sep=""))
  
  O <- MP.Subset.Dataset(
    input.ds =        paste(output.dir, identifier, ".model_dataset.H.gct", sep=""),
    input.cls =       paste(output.dir, identifier, ".model_set.2.cls", sep=""),
    column.subset =   "ALL",
    column.sel.type = "samples",
    row.subset =      "ALL", 
    output.ds =       paste(output.dir, identifier, ".all.H.gct", sep=""),
    output.cls =      paste(output.dir, identifier, ".all.H.cls", sep=""))
  
  O <- MP.Subset.Dataset(
    input.ds =          paste(output.dir, identifier, ".model_set.2.gct", sep=""),          
    input.cls =         paste(output.dir, identifier, ".model_set.2.cls", sep=""),
    column.subset =   "ALL",
    column.sel.type = "samples",
    row.subset =      "ALL", 
    output.ds =       paste(output.dir, identifier, ".all.gct", sep=""),
    output.cls =      paste(output.dir, identifier, ".all.cls", sep=""))
  
  # Pre-process and project all the test datasets
  
  if (! is.null(test.datasets.table)) {
    max.files <- length(test.datasets.table)
    for (ds in 1:max.files) {
      
      print(c("Processing test file: ", test.datasets.table[[ds]]$gct.file))
      
      O <- MP.Subset.Dataset(
        input.ds =         test.datasets.table[[ds]]$gct.file,
        input.cls =        test.datasets.table[[ds]]$cls.file,
        column.subset =    test.datasets.table[[ds]]$column.subset,
        column.sel.type =  test.datasets.table[[ds]]$column.sel.type,
        row.subset =       "ALL", 
        output.ds =        paste(output.dir, "temp1.gct", sep=""),
        output.cls =       paste(output.dir, "temp1.cls", sep=""))
      
      O <- MP.Preprocess.Dataset(
        input.ds =         paste(output.dir, "temp1.gct", sep=""),
        output.ds =        paste(output.dir, "temp2.gct", sep=""),
        thres =            test.datasets.table[[ds]]$thres,
        ceil =             test.datasets.table[[ds]]$ceil,
        normalization =    "NULL") 
      
      O <- MP.Match.and.Select(
        input1.ds =      paste(output.dir, identifier, ".model_set.2.gct", sep=""),          
        input2.ds =      paste(output.dir, "temp2.gct", sep=""),
        output.ds =      paste(output.dir, "temp3.gct", sep=""))
      
      O <- MP.Preprocess.Dataset(
        input.ds =         paste(output.dir, "temp3.gct", sep=""),
        output.ds =        paste(output.dir, "temp4.gct", sep=""),
        normalization =    test.datasets.table[[ds]]$norm)
      
      O <- MP.Factors.Project(
        input.ds =         paste(output.dir, "temp4.gct", sep=""),
        factors.ds =       paste(output.dir, identifier, ".model.W.gct", sep=""),
        postprojnorm =     postprojnorm,
        output.file =      paste(output.dir, "temp5.gct", sep=""))
      
      O <- MP.Match.and.Merge(
        input1.ds =        paste(output.dir, identifier, ".all.H.gct", sep=""),
        input1.cls =       paste(output.dir, identifier, ".all.H.cls", sep=""),
        input2.ds =        paste(output.dir, "temp5.gct", sep=""),
        input2.cls =       paste(output.dir, "temp1.cls", sep=""),
        output.ds =        paste(output.dir, identifier, ".all.H.gct", sep=""),
        output.cls =       paste(output.dir, identifier, ".all.H.cls", sep=""))
      
      O <- MP.Match.and.Merge(
        input1.ds =         paste(output.dir, identifier, ".all.gct", sep=""),
        input1.cls =        paste(output.dir, identifier, ".all.cls", sep=""),
        input2.ds =         paste(output.dir, "temp4.gct", sep=""),
        input2.cls =        paste(output.dir, "temp1.cls", sep=""),
        output.ds =         paste(output.dir, identifier, ".all.gct", sep=""),
        output.cls =        paste(output.dir, identifier, ".all.cls", sep=""))
      
      print(c("Done processing test file: ", test.datasets.table[[ds]]$gct.file))
    }
    
    unlink(paste(output.dir, "temp1.gct", sep=""))
    unlink(paste(output.dir, "temp2.gct", sep=""))
    unlink(paste(output.dir, "temp3.gct", sep=""))
    unlink(paste(output.dir, "temp4.gct", sep=""))
    unlink(paste(output.dir, "temp5.gct", sep=""))
    unlink(paste(output.dir, "temp1.cls", sep=""))
    
  }
  
  CLS <- MP.ReadClsFile(file =  paste(output.dir, identifier, ".model_set.2.cls", sep=""))
  model.size <- length(CLS$class.v)
  
  O <- MP.Projection.Plots(
    input.ds =                    paste(output.dir, identifier, ".all.H.gct", sep=""),
    input.cls =                   paste(output.dir, identifier, ".all.H.cls", sep=""),
    model.set =                   seq(1, model.size),
    output.2D.proj.plot =         paste(output.dir, identifier, ".2D.proj", sep=""),
    output.heatmap.plot =         paste(output.dir, identifier, ".heatmap", sep=""),
    output.heatmap.sorted.plot =  paste(output.dir, identifier, ".heatmap.sorted", sep=""),
    title =                       identifier,
    seed =                        seed, 
    heatmap.row.norm =            heatmap.row.norm,
    heatmap.cmap.type =           heatmap.cmap.type,
    symbol.scaling =              symbol.scaling,
    col =                         col,
    symbs =                       symbs)
  
  # Evaluate projection 
  
  O <- MP.Evaluate.Projection(
    input.ds =                    paste(output.dir, identifier, ".all.H.gct", sep=""),
    input.cls =                   paste(output.dir, identifier, ".all.H.cls", sep=""),
    model.set =                   seq(1, model.size),
    prediction.results.file =     paste(output.dir, identifier, ".pred.txt", sep=""),
    prediction.matrix.file =      paste(output.dir, identifier, ".pred.gct", sep=""),
    col =                         col,
    use.feature.names =           use.feature.names,
    nchar.phen =                  nchar.phen,
    high.conf.thres =             high.conf.thres,
    symbol.scaling =              symbol.scaling,
    symbs          =              symbs,
    levels =                      levels,
    nlevels =                     nlevels,
    kernel =                      kernel,
    cost =                        cost,
    gamma =                       gamma)
  
  # Compute hierarchical clustering
  
  input.ds <- paste(output.dir, identifier, ".all.gct", sep="")
  input.cls <- paste(output.dir, identifier, ".all.cls", sep="")
  dataset <- MP.Gct2Frame(filename = input.ds)
  m.ds <- data.matrix(dataset$ds)
  N.ds <- length(m.ds[,1])
  M.ds <- length(m.ds[1,])
  ds.names <- dataset$row.names
  ds.descs <- dataset$descs
  ds.sample.names <- dataset$names
  
  # Read class labels
  
  CLS <- MP.ReadClsFile(file=input.cls)
  class.labels <- CLS$class.v
  class.phen <- CLS$phen
  class.list <- CLS$class.list 
  class.labels <- match(class.list, class.phen)
  
  # Compute hierarchical tree clustering
  
  dist.matrix <- dist(t(m.ds))
  HC <- hclust(dist.matrix, method="complete")
  
  quartz(height = 20, width = 30)
  nf <- layout(matrix(c(1, 2), 1, 2, byrow=T), widths = c(7, 1), respect = FALSE)
  
  HC$labels <- class.list
  dhc <- as.dendrogram(HC, hang = 0.05, edge.root = T, dLeaf = 4, edgePar = list(col = c("blue", "green"), lty = c(1, 1), lwd = c(2, 2), t.col = c(1, 1)))
  local({
    colLab <<- function(n) {
      if(is.leaf(n)) {
        a <- attributes(n)
        i <<- i+1
        attr(n, "nodePar") <-
          #                    c(a$nodePar, list(lab.col = mycols[i], pch = c(21, 21), col = c(1, 1), bg = c(mycols[i], mycols[i]), cex = c(0.8, 0.8),
          c(a$nodePar, list(lab.col = 0, pch = c(mysymbs[i], mysymbs[i]), col = c(1, 1), bg = c(mycols[i], mycols[i]), cex = c(1.5, 1.5),
                            lab.font= i%%1))
      }
      n
    }
    mycols <- col[class.labels[HC$order]]
    mysymbs <- symbs[class.labels[HC$order]]
    i <- 0
  })
  dL <- dendrapply(dhc, colLab)
  plot(dL, cex=1, edge.root = T, main = " Hierarchical Clustering (original data)", xlab = "samples") ## --> colored labels!
  
  leg.txt <- class.phen
  n.phen <- length(class.phen)
  p.vec <- symbs[1:n.phen]
  c.vec <- col[1:n.phen]
  par(mar = c(0, 0, 0, 0))
  plot(c(0,0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes=F, type="n", xlab = "", ylab="")
  legend(x=0, y=1, legend=leg.txt, bty="n", xjust=0, yjust= 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = symbol.scaling*2.00, pt.cex=symbol.scaling*3, x.intersp = 2, y.intersp = 5)
  
  plot.filename <- paste(output.dir, identifier, ".htree", sep="")
  jpeg(filename = plot.filename, type ="quartz")
  
  # Read projected dataset 
  
  input.ds <- paste(output.dir, identifier, ".all.H.gct", sep="")
  
  input.cls <- paste(output.dir, identifier, ".all.H.cls", sep="")
  dataset <- MP.Gct2Frame(filename = input.ds)
  m.ds <- data.matrix(dataset$ds)
  N.ds <- length(m.ds[,1])
  M.ds <- length(m.ds[1,])
  ds.names <- dataset$row.names
  ds.descs <- dataset$descs
  ds.sample.names <- dataset$names
  
  # Compute hierarchical tree clustering
  
  dist.matrix <- dist(t(m.ds))
  HC <- hclust(dist.matrix, method="complete")
  
  quartz(height = 20, width = 30)
  nf <- layout(matrix(c(1, 2), 1, 2, byrow=T), widths = c(7, 1), respect = FALSE)
  
  HC$labels <- class.list
  dhc <- as.dendrogram(HC, hang = 0.05, edge.root = T, dLeaf = 4, edgePar = list(col = c("blue", "green"), lty = c(1, 1), lwd = c(2, 2), t.col = c(1, 1)))
  local({
    colLab <<- function(n) {
      if(is.leaf(n)) {
        a <- attributes(n)
        i <<- i+1
        attr(n, "nodePar") <-
          #                    c(a$nodePar, list(lab.col = mycols[i], pch = c(21, 21), col = c(1, 1), bg = c(mycols[i], mycols[i]), cex = c(0.8, 0.8),
          c(a$nodePar, list(lab.col = 0, pch = c(mysymbs[i], mysymbs[i]), col = c(1, 1), bg = c(mycols[i], mycols[i]), cex = c(1.5, 1.5),
                            lab.font= i%%1))
      }
      n
    }
    mycols <- col[class.labels[HC$order]]
    mysymbs <- symbs[class.labels[HC$order]]
    i <- 0
  })
  dL <- dendrapply(dhc, colLab)
  plot(dL, cex=1, edge.root = T, main = " Hierarchical Clustering (projected data)", xlab = "samples") ## --> colored labels!
  
  leg.txt <- class.phen
  n.phen <- length(class.phen)
  p.vec <- symbs[1:n.phen]
  c.vec <- col[1:n.phen]
  par(mar = c(0, 0, 0, 0))
  plot(c(0,0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes=F, type="n", xlab = "", ylab="")
  legend(x=0, y=1, legend=leg.txt, bty="n", xjust=0, yjust= 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = symbol.scaling*1.35, pt.cex=symbol.scaling*3, x.intersp = 2, y.intersp = 5)
  
  plot.filename <- paste(output.dir, identifier, ".H.htree", sep="")
  jpeg(filename = plot.filename, type ="quartz")
  
  # Compute class membership
  
  membership <- vector(length=M.ds, mode="numeric")
  
  for (j in 1:M.ds) { # Find membership
    membership[j] <- order(m.ds[,j], decreasing=T)
  }
  
  mem.order <- order(membership, decreasing=F)
  membership.sorted <- membership[mem.order]
  ds.sample.names <- paste(class.list, ds.sample.names, sep="_")
  ds.sample.names.sorted <- ds.sample.names[mem.order]
  class.list.sorted <- class.list[mem.order]
  
  mem.table <- data.frame(cbind(class.list, ds.sample.names, membership, rep(" ", M.ds), class.list.sorted, ds.sample.names.sorted, membership.sorted))
  row.names(mem.table) <- seq(1, M.ds)
  names(mem.table) <- c("Phen", "Sample Names", "Membership", " ", "Phen Sorted", "Sample Names Sorted", "Membership Sorted")
  
  mem.filename <- paste(output.dir, identifier, ".H.mem.txt", sep="")
  
  write.table(file = mem.filename, mem.table, quote=F, sep="\t")
  
  table(class.list.sorted, membership.sorted)
  
}

MSIG.Match.and.Merge.2 <- function(
  input1.ds,
  input1.cls = "",
  input2.ds,
  input2.cls = "",
  output.ds,
  output.cls,
  mode = "intersection") {  # gene set mode = intersection (default) or union
  
  # start of methodology
  
  print(c("Running MSIG.Match.and.Merge... on: ", input1.ds, " ", input2.ds))
  
  # Read input datasets
  
  dataset1 <- MP.Gct2Frame(filename = input1.ds)
  m1 <- data.matrix(dataset1$ds)
  gs.names1 <- dataset1$row.names
  gs.descs1 <- dataset1$descs
  sample.names1 <- dataset1$names
  N1 <- length(m1[1,])
  
  dataset2 <- MP.Gct2Frame(filename = input2.ds)
  m2 <- data.matrix(dataset2$ds)
  gs.names2 <- dataset2$row.names
  gs.descs2 <- dataset2$descs
  sample.names2 <- dataset2$names
  N2 <- length(m2[1,])
  
  # Read CLS files 
  
  if ((input1.cls != "") & (input2.cls != "")) {
    CLS1 <- MP.ReadClsFile(file=input1.cls)
    class.labels1 <- CLS1$class.v
    class.phen1 <- CLS1$phen
    
    CLS2 <- MP.ReadClsFile(file=input2.cls)
    class.labels2 <- CLS2$class.v
    class.phen2 <- CLS2$phen
  }
  
  # Match features to first dataset and create matching m2 dataset
  
  sample.names3 <- c(sample.names1, sample.names2)
  N3 <- N1 + N2
  
  if (mode == "intersection") {
    gs.names3 <- intersect(gs.names1, gs.names2)
    
    locations1 <- match(gs.names3, gs.names1, nomatch=0)
    m1 <- m1[locations1, ]
    gs.descs1 <- gs.descs1[locations1]
    
    locations2 <- match(gs.names3, gs.names2, nomatch=0)
    m2 <- m2[locations2, ]
    gs.descs2 <- gs.descs2[locations2]
    m3 <- cbind(m1, m2)
    
  } else if (mode == "union") {
    
    gs.names3 <- union(gs.names1, gs.names2)
    M3 <- length(gs.names3)
    
    m3 <- matrix(0, nrow = M3, ncol=N3)
    
    locations1 <- match(gs.names3, gs.names1, nomatch=0)
    locations2 <- match(gs.names3, gs.names2, nomatch=0)
    
    for (i in 1:M3) {
      if (locations1[i] != 0) {
        m3[i, 1:N1] <- m1[locations1[i],]
      }
      if (locations2[i] != 0) {
        m3[i, (N1+1):N3] <- m2[locations2[i],]
      }
    }
  } else {
    stop(c("unknown mode", mode))
  }
  
  
  if ((input1.cls != "") & (input2.cls != "")) {
    class.labels3 <- c(class.labels1, class.labels2 + length(class.phen1))
    class.phen3 <- c(class.phen1, class.phen2)
  }
  
  # Save datasets
  
  V <- data.frame(m3)
  names(V) <- sample.names3
  row.names(V) <- gs.names3
  gs.descs1 <- gs.names3
  MP.WriteGct(gct.data.frame = V, descs = gs.descs1, filename = output.ds)  
  
  if ((input1.cls != "") & (input2.cls != "")) {
    write.cls(class.v = class.labels3, phen = class.phen3, filename = output.cls) 
  }
}

MP.Reconstruct.Dataset <- function(
  input.H.ds, 
  input.W.ds,
  output.file) {
  
  library(MASS)
  
  # start of methodology
  
  print(c("Running MP.Reconstruct.Dataset... on: ", input.H.ds, input.W.ds))
  
  # Read input datasets
  
  dataset <- MP.Gct2Frame(filename = input.W.ds)
  W <- data.matrix(dataset$ds)
  W.row.names <- dataset$row.names
  W.row.descs <- dataset$descs
  W.names <- dataset$names
  
  dataset <- MP.Gct2Frame(filename = input.H.ds)
  H <- data.matrix(dataset$ds)
  H.row.names <- dataset$row.names
  H.row.descs <- dataset$descs
  H.names <- dataset$names
  
  # Project input dataset using factors input
  
  A <- W %*% H
  
  # Save reconstructed dataset
  
  V <- data.frame(A)
  names(V) <- H.names
  row.names(V) <- W.row.names
  MP.WriteGct(gct.data.frame = V, filename = output.file)  
  
}

MP.HeatMapPlot <- function(
  V, 
  row.names = "NA", 
  col.labels = "NA", 
  col.classes = "NA", 
  col.names = "NA", 
  main = " ", 
  sub = " ", 
  xlab=" ", 
  ylab=" ",
  row.norm = TRUE,
  char.rescale = 1.0,                               
  cmap.type = 1,   # 1 = vintage pinkogram, 2 = scale of grays, 3 = high-resolution pinkogram for probabilities [0, 1], 4 = high-resolution pinkogram for general values, 5 = color map for normalized enrichment scores, 6 = color map for metagene factors
  max.v = "NA",
  rotated.col.labels = F)
{
  #
  # Plots a heatmap "pinkogram" of a gene expression matrix including phenotype vector and gene, sample and phenotype labels
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  
  n.rows <- length(V[,1])
  n.cols <- length(V[1,])
  n.phen <- length(col.classes)
  V1 <- matrix(0, nrow=n.rows, ncol=n.cols)
  
  if ((cmap.type == 3) | (cmap.type == 5)) {
    row.norm <- F
  }
  
  if (row.norm == TRUE) {
    row.mean <- apply(V, MARGIN=1, FUN=mean)
    row.sd <- apply(V, MARGIN=1, FUN=sd)
    row.n <- length(V[,1])
    for (i in 1:n.rows) {
      if (row.sd[i] == 0) {
        V1[i,] <- 0
      } else {
        V1[i,] <- (V[i,] - row.mean[i])/(0.333 * row.sd[i])
      }
      V1[i,] <- ifelse(V1[i,] < -6, -6, V1[i,])
      V1[i,] <- ifelse(V1[i,] > 6, 6, V1[i,])
    }
  } else {
    V1 <- V
  }
  
  if (cmap.type == 1) { 
    mycol <- c("#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0", "#FF7080", 
               "#FF5A5A", "#FF4040", "#FF0D1D") # blue-pinkogram colors. This is the 1998-vintage, pre-gene cluster, original pinkogram color map
  } else if (cmap.type == 2) {
    mycol <- vector(length=256, mode = "numeric")
    for (k in 1:256) {
      red <-  (k - 1)*0.80 + 50
      green <- (k - 1)*0.68 + 80
      blue <- (k - 1)*0.60 + 100
      mycol[k] <- rgb(red, green, blue, maxColorValue=255)
    }
    mycol <- rev(mycol)
    
  } else if ((cmap.type == 3) | (cmap.type == 4) | (cmap.type == 5)) {
    mycol <- vector(length=512, mode = "numeric")
    
    for (k in 1:256) {
      mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
    }
    for (k in 257:512) {
      mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
    }
    mycol <- rev(mycol)
  } else if (cmap.type == 6) {
    
    mycol <- vector(length=256, mode = "numeric")
    max.pos.V1 <- max(V1)
    if (min(V1) < 0) {
      min.neg.V1 <- min(V1)
    } else {
      min.neg.V1 <- 0
    }
    neg.k <- ceiling(256*(abs(min.neg.V1)/(max.pos.V1 - min.neg.V1)))
    print(c("neg.k=", neg.k, " max.pos.V1=", max.pos.V1, "min.neg.V1=", min.neg.V1)) 
    neg.frac <- ifelse(abs(min.neg.V1/max.pos.V1) < 1, abs(min.neg.V1/max.pos.V1), 1)
    for (k in 1:neg.k) {
      max.red <- 255 -  (255 - 50)*neg.frac
      min.red <- 255
      red <-  max.red + (min.red - max.red) * (k - 1)/(neg.k - 1)
      
      max.green <- 255 - (255 - 50)*neg.frac
      min.green <- 255
      green <-  max.green + (min.green - max.green) * (k - 1)/(neg.k - 1)
      
      max.blue <- 255 - (255 - 205)*neg.frac
      min.blue <- 255
      blue <-  max.blue + (min.blue - max.blue) * (k - 1)/(neg.k - 1)
      
      print(c(k, red, green, blue))
      mycol[k] <- rgb(red, green, blue, maxColorValue=255)
    }
    for (k in (neg.k + 1):256) {
      max.red <- 205
      min.red <- 255
      red <-  min.red - (min.red - max.red) * (k - (neg.k + 1))/(256 - (neg.k + 1))
      max.green <- 50
      min.green <- 255
      green <-  min.green - (min.green - max.green) * (k - (neg.k + 1))/(256 - (neg.k + 1))
      max.blue <- 50
      min.blue <- 255
      blue <-  min.blue - (min.blue - max.blue) * (k - (neg.k + 1))/(256 - (neg.k + 1))
      mycol[k] <- rgb(red, green, blue, maxColorValue=255)
    }
  }
  
  ncolors <- length(mycol)
  
  if (cmap.type == 5) {
    if (max.v == "NA") {
      max.v <- max(max(V1), -min(V1))
    }
    V2 <- ceiling(ncolors * (V1 - (- max.v))/(1.001*(max.v - (- max.v))))
    
  } else if (cmap.type == 3) {
    V2 <- ceiling(ncolors * (V1/1.001))
  } else {
    V2 <- ceiling(ncolors * (V1 - min(V1))/(1.001*(max(V1) - min(V1))))
  }
  
  heatm <- matrix(0, nrow = n.rows, ncol = n.cols)
  heatm[1:n.rows,] <- V2[, seq(n.cols, 1, -1)]
  
  col.labels <- col.labels[seq(n.cols, 1, -1)]
  if (col.names[1] != "NA") {
    col.names <-  col.names[seq(n.cols, 1, -1)]
  }        
  height <- ifelse(n.rows >= 25, 25, n.rows*0.8 + 2)
  
  quartz(width=31, height=19)
  
  nf <- layout(matrix(c(1, 2), 1, 2, byrow=T), widths = c(7, 1), respect = FALSE)
  
  n.rows2 <- ifelse(n.rows < 4, 4, n.rows)
  a <- -12/16
  b <- 27
  margin <- a*n.rows2 + b
  margin <- ifelse(margin < 2, 2, ifelse(margin > 24, 24, margin))
  
  if (rotated.col.labels == F) {            
    par(mar = c(4, margin, 4, 10))
  } else {
    par(mar = c(4, margin, 10, 10))
  }
  
  #       print(heatm)
  
  image(1:n.rows, 1:n.cols, heatm, zlim = c(0, ncolors), col=mycol, axes=FALSE, main=main, sub = sub, xlab= xlab, ylab=ylab)
  
  if (row.names[1] != "NA") {
    numC <- nchar(row.names)
    size.row.char <- char.rescale*15/(n.rows + 15)
    size.col.char <- char.rescale*30/(n.cols + 15)
    size.lab.char <- char.rescale*30/(n.phen + 15)
    for (i in 1:n.rows) {
      row.names[i] <- substr(row.names[i], 1, 20)
    }
    if (rotated.col.labels == F) {            
      axis(3, at=1:n.rows, labels=row.names, adj= 1, tick=FALSE, las = 1, cex.axis=size.row.char, font.axis=1, line=-0.75)
    } else {
      axis(3, at=1:n.rows, labels=row.names, adj= 1, tick=FALSE, las = 3, cex.axis=size.row.char, font.axis=1, line=-0.75)
    }
  }
  
  #       print(c(" col.names=", col.names))
  
  if (col.names[1] != "NA") {
    axis(2, at=1:n.cols, labels=col.names, tick=FALSE, las = 1, cex.axis=size.col.char, font.axis=2.5, line=-1)
  }
  
  for (i in 1:(n.rows)) {
    lines(x = c(i + 0.5, i + 0.5), y = c(0.5, n.cols + 0.5), type = "l", lwd =1, col = "black")
  }
  
  boundaries <- cumsum(sapply(split(rep(1, n.cols), col.labels), sum))
  boundaries <- n.cols - boundaries
  
  loc.phen <- vector(length=n.phen, mode="numeric")
  for (i in 1:(n.phen)) {
    lines(x = c(1 - 0.5, n.rows + 0.5), y = c(boundaries[i] + 0.5, boundaries[i] + 0.5), type = "l", lwd =1, col = "black")
    if (i > 1) {
      loc.phen[i] <- mean(c(boundaries[i - 1], boundaries[i]))
    } else {
      loc.phen[i] <- mean(c(n.cols, boundaries[i]))
    }
    #          print(c(i, " line=", boundaries[i], " label=", loc.phen[i]))
  }
  axis(4, at=loc.phen, labels=col.classes, tick=FALSE, las = 1, cex.axis=size.lab.char, font.axis=1, line=-0.75)
  
  
  lines(x = c(0.50, n.rows + 0.50), y = c(0.50, 0.50), type = "l", lwd =1, col = "black")
  lines(x = c(0.50, n.rows + 0.50), y = c(n.cols + 0.50, n.cols + 0.50), type = "l", lwd =1, col = "black")
  lines(x = c(0.50, 0.50), y = c(0.50, n.cols + 0.50), type = "l", lwd =1, col = "black")
  lines(x = c(n.rows + 0.50, n.rows + 0.50), y = c(0.50, n.cols + 0.50), type = "l", lwd =1, col = "black")
  
  # Color map legend
  
  #       print(c("range V=", range(V)))
  #       print(c("range V1=", range(V1)))
  #       print(c("range V2=", range(V2)))
  
  par(mar = c(10,2,10,2))
  num.v <- 20
  if (cmap.type == 3) {
    range.v <- c(0, ncolors)
  } else {
    range.v <- range(V2)
  }
  incr <-  (range.v[1] - range.v[2])/(num.v - 1)
  heatm.v <- matrix(rev(seq(range.v[2], range.v[1], incr)), nrow=num.v, ncol=1)
  image(1:1, 1:num.v, t(heatm.v), zlim = c(0, ncolors + max(col.labels)), col=mycol, axes=FALSE, sub="Color \n Legend ", main= " ", xlab= xlab, ylab=ylab)
  if (cmap.type == 3) {
    range.v <- c(0, 1)
  } else {
    range.v <- range(V1)
  }
  incr <-  (range.v[1] - range.v[2])/(num.v - 1)
  heatm.v2 <- matrix(signif(rev(seq(range.v[2], range.v[1], incr)), digits=3), nrow=num.v, ncol=1)
  #          print(c("heatm.v2=", heatm.v2))
  axis(2, at=1:num.v, labels=heatm.v2, adj= 0.5, tick=FALSE, las = 1, cex.axis=char.rescale*0.6, font.axis=1.25, line=-0.8)
  
  lines(x = c(0.5, 1.5), y = c(0.5, 0.5), type = "l", lwd =1, col = "black")
  lines(x = c(0.5, 1.5), y = c(num.v + 0.5, num.v + 0.5), type = "l", lwd =1, col = "black")
  
  lines(x = c(0.6, 0.6), y = c(0.5, num.v + 0.5), type = "l", lwd =1, col = "black")
  lines(x = c(1.4, 1.4), y = c(0.5, num.v + 0.5), type = "l", lwd =1, col = "black")
  
  
  return()
  
}


MP.Evaluate.Projection <- function(
  input.ds,
  input.cls,
  model.set,
  prediction.results.file,
  prediction.matrix.file,
  col = c("grey3", "cadetblue", "darkgreen", "chartreuse2", "red1", "darkred", "orange", "blue2", "lightblue", "pink", "coral"),
  symbs = c(22, 21, 20, 23, 24, 25, 21, 20, 23, 24, 25),
  use.feature.names = F,
  nchar.phen = 3,
  high.conf.thres = 0.75,
  symbol.scaling = symbol.scaling,
  levels = NULL,
  nlevels = 10,
  kernel = "radial",
  cost = 5,
  gamma = 0.05) {
  
  print(c("Running MP.Evaluate.Projection... on:", input.ds))
  
  library(e1071)
  library(tree)
  
  # Read dataset
  
  dataset <- MP.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  max.m <- max(m)
  m <- m/max.m
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names
  
  dim(m) 
  
  Ns <- length(m[1,])
  N <- length(m[,1])
  
  CLS <- MP.ReadClsFile(file=input.cls)
  class.labels <- CLS$class.v
  class.list <- CLS$class.list
  class.phen <- CLS$phen
  
  num.classes <- length(class.phen)
  
  print("Reading dataset completed...")
  
  # Use first nchar.phen characters of phenotype class to define new phenotypes
  
  class.list2 <- vector(length = Ns, mode = "character")
  for (i in 1:Ns) {
    class.list2[i] <- substr(class.list[i], 1, nchar.phen)
  }
  class.phen2 <- vector(length = num.classes, mode = "character")
  for (i in 1:num.classes) {
    class.phen2[i] <- substr(class.phen[i], 1, nchar.phen)
  }
  true.num.classes <- length(table(class.phen2))
  
  class.labels2 <- match(class.list2, class.phen2)
  
  # Separate data into train and test pieces
  
  m.train <- m[,model.set]
  n.train <- length(model.set)
  num.samples.train <- n.train
  sample.names.train <- as.factor(sample.names[model.set])
  class.list.train <- class.list2[model.set]
  class.phen.train <- unique(class.list.train)
  class.labels.train <- class.labels2[model.set]
  orig.class.labels.train <- class.labels[model.set]
  
  if (Ns - length(model.set) > 0) { 
    m.test <- as.matrix(m[, - model.set])
    n.test <- length(m.test[1,])
    sample.names.test <- as.factor(sample.names[- model.set])
    class.list.test <- class.list2[- model.set]
    class.phen.test <- unique(class.list.test)
    class.labels.test <- class.labels2[- model.set]
  }
  
  # Build SVM and tree models
  
  print("Building SVM model...")
  
  one.over <- function(x) { return(100/length(x)) }
  class.number.list <- split(rep(1, length(class.list.train)) , class.list.train)
  class.weights  <- sapply(class.number.list, one.over)
  print(c("class.weights=", class.weights))
  
  svm.model <- svm(x = t(m.train), y = class.list.train, scale = F, tolerance = 0.00001, type = "C-classification", kernel = kernel, cost = cost, gamma = gamma, probability = T)
  
  print("Computing train set predictions...")
  
  train.pred <- predict(object = svm.model, newdata = t(m.train), decision.values = T, probability = T)  
  
  dec.vals.train <- attr(train.pred, "decision.values")
  prob.train <- signif(attr(train.pred, "probabilities"), digits=2)
  confidence.vector <- vector(length=n.train, mode="numeric")
  bscore <- vector(length=n.train, mode = "numeric")
  max.k <- length(prob.train[1,])
  random.pred.conf <- ((max.k - 1)/max.k)^2 + (max.k - 1)*(1/max.k)^2
  for (ii in 1:n.train) {
    probs <- sort(prob.train[ii,], decreasing=T)
    confidence.vector[ii] <-  1 - ((1 - probs[1])^2 + sum(probs[2:max.k]^2))/random.pred.conf
    confidence.vector[ii] <- signif(confidence.vector[ii], digits=3)
    if (class.list.train[ii] == as.character(train.pred[ii])) {
      bscore[ii] <- signif((1 - probs[1])^2, digits=2)
    } else {
      bscore[ii] <- signif(probs[1]^2, digits=2)
    }
  }
  confidence.call <- ifelse(confidence.vector >= high.conf.thres, " H ", " L ")
  error.call <- ifelse(class.list.train == as.character(train.pred), "   ", " * ")
  no.call <- ifelse(confidence.vector >= high.conf.thres, 0, 1)
  real.error <- ifelse(((no.call == 0) & (error.call == " * ")), 1, 0)
  correct.call <- ifelse(((no.call == 0) & (error.call == "   ")), 1, 0)
  
  col.symbols.train <- paste(confidence.call, error.call)
  class.names <- names(data.frame(prob.train))
  Brier.train <- signif(mean(bscore), digits=2)
  
  train.results <- data.frame(cbind(as.character(sample.names.train), class.list.train, as.character(train.pred), error.call, confidence.call, confidence.vector, no.call, real.error, correct.call, prob.train, bscore))
  names(train.results)[1] <- "Train Sample Name"
  names(train.results)[2] <- "Actual"
  names(train.results)[3] <- "Predicted"
  names(train.results)[4] <- "Error (*)"
  names(train.results)[5] <- "Conf (H/L)"
  names(train.results)[6] <- "Conf"
  names(train.results)[7] <- "No Call"
  names(train.results)[8] <- "Real Error"
  names(train.results)[9] <- "Correct Call"
  
  names(train.results)[10 + length(class.phen.train)] <- "Brier score"
  #   print(train.results)
  print(c("Brier score (Train) = ", Brier.train))
  
  write("Training Results \n", file = prediction.results.file, append = F)
  write.table(train.results, file = prediction.results.file, append = T, quote=F, row.names=F, sep = "\t")
  
  write(c("\n\n Brier score (Train) = ", Brier.train), file = prediction.results.file, append = T)
  
  no.call.list <- split(no.call, class.list.train)
  real.error.list <- split(real.error, class.list.train)
  correct.call.list <- split(correct.call, class.list.train)
  count.class <- c(sapply(no.call.list, length), length(no.call))
  no.call.class <- c(sapply(no.call.list, sum), sum(no.call))
  real.error.class <- c(sapply(real.error.list, sum), sum(real.error))
  correct.call.class <- c(sapply(correct.call.list, sum), sum(correct.call))
  train.pred.high.conf <- ifelse(no.call == 0,  as.character(train.pred), "-- no call")
  #   print(c("train.pred.high.conf =", train.pred.high.conf))
  
  no.call.class.pct <- no.call.class/count.class
  real.error.class.pct <- real.error.class/count.class
  correct.call.class.pct <- correct.call.class/count.class
  perf.table.train <- data.frame(cbind(c(names(no.call.list), "Total"), count.class, no.call.class, no.call.class.pct, real.error.class, real.error.class.pct, correct.call.class, correct.call.class.pct))
  names(perf.table.train) <-  c("Class", "Count", "No Call", "No Call (%)", "Real Error", "Real Error (%)", "Correct Call", "Correct Call (%)")
  write.table(perf.table.train, file = prediction.results.file, append = T, quote=F, row.names=F, sep = "\t")
  print(perf.table.train)
  
  conf.table.train <- table(class.list.train, train.pred.high.conf)
  conf.table.train <- data.frame(cbind(row.names(conf.table.train), conf.table.train))
  print(conf.table.train)
  write("\n\n Confusion Matrix (Train) \n", file = prediction.results.file, append = T)
  write.table(conf.table.train, file = prediction.results.file, append = T, quote=F, row.names=F, sep = "\t")
  
  print("Building SVM model completed. Predicting test data...")
  
  if (Ns - length(model.set) > 0) { 
    
    test.pred <- predict(object = svm.model, newdata = t(m.test), decision.values = T, probability = T)  
    dec.vals.test <- attr(test.pred, "decision.values")
    prob.test <- signif(attr(test.pred, "probabilities"), digits=2)
    confidence.vector <- vector(length=n.test, mode="numeric")
    bscore <- vector(length=n.test, mode = "numeric")
    max.k <- length(prob.train[1,])
    random.pred.conf <- ((max.k - 1)/max.k)^2 + (max.k - 1)*(1/max.k)^2
    for (ii in 1:n.test) {
      probs <- sort(prob.test[ii,], decreasing=T)
      confidence.vector[ii] <-  1 - ((1 - probs[1])^2 + sum(probs[2:max.k]^2))/random.pred.conf
      confidence.vector[ii] <- signif(confidence.vector[ii], digits=3)
      if (class.list.test[ii] == as.character(test.pred[ii])) {
        bscore[ii] <- signif((1 - probs[1])^2, digits=2)
      } else {
        bscore[ii] <- signif(probs[1]^2, digits=2)
      }
      
    }
    
    confidence.call <- ifelse(confidence.vector >= high.conf.thres, " H ", " L ")
    error.call <- ifelse(class.list.test == as.character(test.pred), "   ", " * ")
    no.call <- ifelse(confidence.vector >= high.conf.thres, 0, 1)
    real.error <- ifelse(((no.call == 0) & (error.call == " * ")), 1, 0)
    correct.call <- ifelse(((no.call == 0) & (error.call == "   ")), 1, 0)
    col.symbols.test <- paste(confidence.call, error.call)
    class.names <- names(data.frame(prob.test))
    Brier.test <- signif(mean(bscore), digits=2)
    
    test.results <- data.frame(cbind(as.character(sample.names.test), class.list.test, as.character(test.pred), error.call, confidence.call, confidence.vector, no.call, real.error, correct.call, prob.test, bscore))
    names(test.results)[1] <- "Test Sample Name"
    names(test.results)[2] <- "Actual"
    names(test.results)[3] <- "Predicted"
    names(test.results)[4] <- "Error (*)"
    names(test.results)[5] <- "Conf (H/L)"
    names(test.results)[6] <- "Conf"
    names(test.results)[7] <- "No Call"
    names(test.results)[8] <- "Real Error"
    names(test.results)[9] <- "Correct Call"
    
    names(test.results)[10 + length(class.phen.train)] <- "Brier score"
    #      print(test.results)
    print(c("Brier score (Test) = ", Brier.test))
    
    write("\n Test Results \n", file = prediction.results.file, append = T)
    write.table(test.results, file = prediction.results.file, append = T, quote=F, row.names=F, sep = "\t")
    
    write(c("\n\n Brier score (Test) = ", Brier.test), file = prediction.results.file, append = T)
    
    no.call.list <- split(no.call, class.list.test)
    real.error.list <- split(real.error, class.list.test)
    correct.call.list <- split(correct.call, class.list.test)
    count.class <- c(sapply(no.call.list, length), length(no.call))
    no.call.class <- c(sapply(no.call.list, sum), sum(no.call))
    real.error.class <- c(sapply(real.error.list, sum), sum(real.error))
    correct.call.class <- c(sapply(correct.call.list, sum), sum(correct.call))
    test.pred.high.conf <- ifelse(no.call == 0,  as.character(test.pred), "-- no call")
    #   print(c("test.pred.high.conf =", test.pred.high.conf))
    
    no.call.class.pct <- no.call.class/count.class
    real.error.class.pct <- real.error.class/count.class
    correct.call.class.pct <- correct.call.class/count.class
    
    perf.table.test <- data.frame(cbind(c(names(no.call.list), "Total"), count.class, no.call.class, no.call.class.pct, real.error.class, real.error.class.pct, correct.call.class, correct.call.class.pct))
    names(perf.table.test) <-  c("Class", "Count", "No Call", "No Call (%)", "Real Error", "Real Error (%)", "Correct Call", "Correct Call (%)")
    write.table(perf.table.test, file = prediction.results.file, append = T, quote=F, row.names=F, sep = "\t")
    print(perf.table.test)
    
    conf.table.test <- table(class.list.test, test.pred.high.conf)
    conf.table.test <- data.frame(cbind(row.names(conf.table.test), conf.table.test))
    print(conf.table.test)
    write("\n\n Confusion Matrix (Test) \n", file = prediction.results.file, append = T)
    write.table(conf.table.test, file = prediction.results.file, append = T, quote=F, row.names=F, sep = "\t")
    
    # Save predictions for all classes in a gct and cls files
    
    #   print(c("dim train=", dim(prob.train)))
    #   print(c("dim test=", dim(prob.test)))
    
    if (Ns - length(model.set) > 0) { 
      V <- cbind(t(prob.train), t(prob.test))
      W <- data.frame(V)
      names(W) <- c(as.character(sample.names.train), as.character(sample.names.test))
      row.names(W) <- names(train.results)[seq(7, 7 + length(class.phen.train) - 1)]
    } else {
      V <- t(prob.train)
      W <- data.frame(V)
      names(W) <- as.character(sample.names.train)
      row.names(W) <- names(train.results)[seq(7, 7 + length(class.phen.train) - 1)]
    }
    #   print(W)
    MP.WriteGct(gct.data.frame = W, descs = row.names(W), filename = prediction.matrix.file)  
    
    print("Done predicting test data...")
  }
  
}     


MP.nnls.fit <- function(x,y,wsqrt=1,eps=0,rank.tol=1e-07) {
  ## Purpose: Nonnegative Least Squares (similar to the S-Plus function
  ## with the same name) with the help of the R-library quadprog
  ## ------------------------------------------------------------------------
  ## Attention:
  ## - weights are square roots of usual weights
  ## - the constraint is coefficient>=eps
  ## ------------------------------------------------------------------------
  ## Author: Marcel Wolbers, July 99
  ##
  ##========================================================================
  require ("quadprog")
  m <- NCOL(x)
  if (length(eps)==1) eps <- rep(eps,m)
  x <- x * wsqrt
  y <- y * wsqrt
  #  sometimes a rescaling of x and y helps (if solve.QP.compact fails otherwise)
  xscale <- apply(abs(x),2,mean)
  yscale <- mean(abs(y))
  x <- t(t(x)/xscale)
  y <- y/yscale
  Rinv <- backsolve(qr.R(qr(x)),diag(m))
  cf <- solve.QP.compact(Dmat=Rinv,dvec=t(x)%*%y,Amat=rbind(rep(1,m)),
                         Aind=rbind(rep(1,m),1:m),bvec=eps*xscale/yscale,
                         factorized=TRUE)$sol
  cf <- cf*yscale/xscale  #scale back
  cf
}

MP.replace.row.col <- function(input.ds, output.ds, mode = "row", number, values) {
  
  dataset <- MP.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  
  
  if (mode == "row") {
    m[number, ] <-  ifelse (length(values) == 1, rep(values, length(m[number,])), values)
  } else if (mode == "col") {
    m[, number] <-  ifelse (length(values) == 1, rep(values, length(m[, number])), values)
  } else {
    stop(c("unknown mode:", mode))
  }
  
  V <- data.frame(m)
  names(V) <- dataset$names
  row.names(V) <- dataset$row.names
  
  MP.WriteGgct(gct.data.frame = V, descs = dataset$descs, filename = output.ds)
}

MP.Subset.Dataset <- function(
  input.ds,
  input.cls = "",
  column.subset = "ALL",    # subset of column numbers or names (or phenotypes)
  column.sel.type = "samples",  # "samples" or "phenotype"
  row.subset = "ALL",       # subset of row numbers or names
  output.ds,
  output.cls) {
  
  # start of methodology
  
  print(c("Running MP.Subset.Dataset... on GCT file:", input.ds, " and CLS file:", input.cls))
  
  # Read input datasets
  
  dataset <- MP.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names
  
  # Read CLS file
  
  if (input.cls != "") {
    CLS <- MP.ReadClsFile(file=input.cls)
    class.labels <- CLS$class.v
    class.phen <- CLS$phen
    class.list <- CLS$class.list 
  }
  
  # Select desired column subset
  
  if (column.sel.type == "samples") {
    if (column.subset[1] == "ALL") {
      m2 <- m
      sample.names2 <- sample.names
      if (input.cls != "") {
        class.labels2 <- class.labels
      }
    } else {
      m2 <- m[,column.subset]
      if (is.numeric(column.subset[1])) {
        sample.names2 <- sample.names[column.subset]
        if (input.cls != "") {
          class.labels2 <- class.labels[column.subset]
        }
      } else {
        locations <- match(sample.names, column.subset)
        sample.names2 <- sample.names[locations]
        if (input.cls != "") {
          class.labels2 <- class.labels[locations]
        }
      }
    }
  } else if (column.sel.type == "phenotype") {
    locations <- !is.na(match(class.list, column.subset))
    sample.names2 <- sample.names[locations]
    m2 <- m[,locations]
    if (input.cls != "") {
      class.labels2 <- class.labels[locations]
    }
  }
  
  if (row.subset[1] == "ALL") {
    m3 <- m2
    gs.names2 <- gs.names
    gs.descs2 <- gs.descs
  } else {
    m3 <- m2[row.subset,]
    if (is.numeric(row.subset[1])) {
      gs.names2 <- gs.names[row.subset]
      gs.descs2 <- gs.descs[row.subset]
    } else {
      locations <- match(row.subset, gs.names)
      gs.names2 <- gs.names[locations]
      gs.descs2 <- gs.descs[locations]
    }
  }
  
  # Save datasets
  
  V <- data.frame(m3)
  names(V) <- sample.names2
  row.names(V) <- gs.names2
  MP.WriteGct(gct.data.frame = V, descs = gs.descs2, filename = output.ds)  
  
  if (input.cls != "") {
    MP.WriteCls(class.v = class.labels2, phen = class.phen, filename = output.cls) 
  }
}

MP.Preprocess.Dataset <- function(
  input.ds, 
  output.ds,
  thres = "NULL", 
  ceil = "NULL", 
  shift = "NULL",
  fold = "NULL", 
  delta = "NULL", 
  normalization = "NULL") {
  
  print(c("Running MP.Preprocess.Dataset... on:", input.ds))
  print(c("output file:", output.ds))
  print(c("normalization =", normalization))
  
  # Read dataset
  
  dataset <- MP.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names
  
  # threshold, ceiling and shift
  
  if (thres != "NULL") {
    m[m < thres] <- thres
  }
  if (ceil != "NULL") {
    m[m > ceil] <- ceil
  }
  if (shift != "NULL") {
    m <- m + shift
  }
  
  if (normalization == 9) {
    m2 <- m
    cols <- length(m2[1,])
    for (j in 1:cols) {  # column rank normalization from 0 to N - 1
      m2[,j] <- rank(m[,j], ties.method = "average") - 1
    }
    m2 <- 10000*m2/(length(m2[,1]) - 1)
    locations <- c(rep(TRUE, length(m2[,1])))
  }
  
  # variation filter
  
  if ((fold != "NULL") && (delta != "NULL")) {
    temp <- MP.VarFilter(V = m, fold = fold, delta = delta, gene.names = gs.names, gene.descs = gs.descs) 
    m <- temp$V
    locations <- temp$locations
    gs.names <- temp$new.gene.names
    gs.descs <- temp$new.gene.descs
    dim(m) 
  }
  
  # normalization
  
  if (normalization == 1) {
    m <- MP.NormalizeCols.Rank(m)
  } else if (normalization == 2) {
    m <- MP.NormalizeCols.Rank(m)/length(m[,1])
  } else if (normalization == 3) {
    m <- MP.NormalizeCols(m) + 3
    m <- MP.Threshold(m, 0.001, 100000) 
  } else if (normalization == 4) {
    m <- MP.NormalizeCols.Rank(m)/length(m[,1])
  } else if (normalization == 5) {
    m <- MP.NormalizeCols.Rescale(m)
  } else if (normalization == 6) {
    cols <- length(m[1,])
    for (j in 1:cols) {  # column rank normalization from 0 to N - 1
      m[,j] <- rank(m[,j], ties.method = "average") - 1
    }
    m <- 10000*m/(length(m[,1]) - 1)
  } else if (normalization == 7) {
    m <- ((100*MP.NormalizeCols.Rank(m))%/%length(m[,1]) + 1)
  } else if (normalization == 8) { 
    row.mean <- apply(m, MARGIN=1, FUN=mean)
    for (i in 1:length(m[,1])) {
      m[i,] <- m[i,] / row.mean[i]
    }
  } else if (normalization == 9) { 
    cols <- length(m[1,])
    size <- sum(locations)
    m <- matrix(0, nrow = size, ncol = cols)
    rows <- length(m2[,1])
    j <- 1
    for (i in 1:rows) {
      if (locations[i]) {
        m[j,] <- m2[i,]
        j <- j + 1
      }
    }
  }
  
  V <- data.frame(m)
  names(V) <- sample.names
  row.names(V) <- gs.names
  MP.WriteGct(gct.data.frame = V, descs = gs.descs, filename = output.ds)  
  
}

MP.Extract.Factors <- function(
  input.ds, 
  input.cls = "", 
  output.W.file, 
  output.H.file, 
  k.proj = 2, 
  alg = "NMF.div",         # decomposition algorithm: NMF.div, NMF, SNMF, NSNMF.div or PCA
  niter = 1000,
  seed = 1234,
  theta = 0,
  lambda = 1,
  sort.factors = F) {
  
  # start of methodology
  
  print(c("Running MP.Extract.Factors... on", input.ds))
  
  
  # Read dataset
  
  dataset <- MP.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names
  
  dim(m) 
  
  Ns <- length(m[1,])
  Ng <- length(m[,1])
  
  if (input.cls != "") {
    CLS <- MP.ReadClsFile(file=input.cls)
    class.labels <- CLS$class.v
    class.phen <- CLS$phen
  } else {
    class.labels <- rep(1, Ns)
    class.phen <- "Samples"
  }
  
  if (alg == "PCA") {    # PCA projection 
    svd.proj <- svd(m, nv = Ns, nu = Ns)
    H.full <- diag(x = svd.proj$d, nrow=Ns, ncol=Ns) %*% t(svd.proj$v)                
    H <- H.full[1:k.proj,]
    W <- svd.proj$u[,1:k.proj]
  } else if (alg == "NMF.div") {  # NMF divergence
    NMF.out <- NMF.div(V = m, k = k.proj, maxniter = niter, seed = seed, stopconv = 40, stopfreq = 10)
    H <- NMF.out$H
    W <- NMF.out$W
  } else if (alg == "NMF") { # NMF Euclidean
    NMF.out <- NMF(V = m, k = k.proj, maxniter = niter, seed = seed, stopconv = 40, stopfreq = 10)
    H <- NMF.out$H
    W <- NMF.out$W
  } else if (alg == "SNMF") { # Sparse NMF Euclidean
    NMF.out <- SNMF(V = m, k = k.proj, maxniter = niter, seed = seed, lambda = lambda, stopconv = 40, stopfreq = 10)
    H <- NMF.out$H
    W <- NMF.out$W
  } else if (alg == "NSNMF.div") { # non-smooth NMF
    NMF.out <- NSNMF.div(V = m, k = k.proj, theta = theta, maxniter = niter, seed = seed, stopconv = 40, stopfreq = 10)
    H <- NMF.out$H
    W <- NMF.out$W
  } else {
    stop (c("unknown algorithm:", alg))
  }
  
  # sort W columns and H rows to make similar projections get similar ordering
  
  if (sort.factors == T) {
    dist.matrix <- dist(t(W))
    HC <- hclust(dist.matrix, method="complete")
    W <- W[, HC$order]
    H <- H[HC$order, ]
  }
  
  factor.names <- paste("F", seq(1, k.proj), sep = "")
  factor.descs <- paste("NMF Extracted Factor Number ", seq(1, k.proj), sep = "")
  
  # save extracted factors datasets W and H
  
  V <- data.frame(W)
  names(V) <- factor.names
  row.names(V) <- gs.names
  MP.WriteGct(gct.data.frame = V, descs = gs.descs, filename = output.W.file)  
  
  V <- data.frame(H)
  names(V) <- sample.names
  row.names(V) <- factor.names
  MP.WriteGct(gct.data.frame = V, descs = factor.descs, filename = output.H.file)  
  
}

MP.Factors.Project <- function(
  input.ds, 
  factors.ds,
  postprojnorm = TRUE,
  output.file,
  method = "pseudo-inverse") {  # method: pseudo-inverse, nnls-solver
  
  library(MASS)
  
  # start of methodology
  
  print(c("Running MP.Factors.Project... on: ", input.ds))
  
  # Read input dataset
  
  dataset <- MP.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names
  
  # Read factors dataset
  
  dataset <- MP.Gct2Frame(filename = factors.ds)
  W <- data.matrix(dataset$ds)
  W.row.names <- dataset$row.names
  W.row.descs <- dataset$descs
  W.names <- dataset$names
  
  # Match features to first dataset and create matching m2 dataset
  
  overlap <- intersect(gs.names, W.row.names)
  
  print(c("Size of Input dataset=", length(gs.names), " genes"))
  print(c("Size of W matrix (rows)=", length(W.row.names), " genes"))
  print(c("Size of overlap=", length(overlap), " genes"))
  
  locations.m <- match(overlap, gs.names, nomatch=0)
  m2 <- m[locations.m, ]
  
  locations.W <- match(overlap, W.row.names, nomatch=0)
  W2 <- W[locations.W, ]
  
  if (method == "pseudo-inverse") {
    
    # Project input dataset using factors input
    
    H <- ginv(W2) %*% m2
    
    # three potential ways to deal with negative values created in the approximated projection
    #
    # I:
    #   max.H <- max(H)
    #   min.H <- min(H)
    #   H <- (H - min.H)/(max.H - min.H)
    #
    # II:
    #   H <- ifelse(H < 0, 0, H)
    #
    # III:
    #  n.col <- length(H[1,])
    #  for (i in 1:n.col) {
    #        max.H <- max(H[,i])
    #        min.H <- min(H[,i])
    #        H[,i] <- (H[,i] - min.H)/(max.H - min.H)
    #  }
    
    print(c("projecting using pseudo-inverse..."))
    
  } else if  (method == "nnls-solver") {  # using a non-negative least square solver 
    
    H <- matrix(0, nrow=length(W2[1,]), ncol= length(m2[1,]))
    
    for (i in 1:length(m2[1,])) {
      H[, i] <- MP.nnls.fit(W2, m2[, i], wsqrt=1, eps=0, rank.tol=1e-07)
    }
    
    print(c("projecting using NNLS solver..."))
    
    
  } else {
    stop("unknown method")
  }
  
  # Normalize projected dataset to the unit hypersphere
  
  if (postprojnorm == TRUE) {
    n.col <- length(H[1,])
    for (i in 1:n.col) {
      S.2 <- sqrt(sum(H[,i]*H[,i]))
      #        S.2 <- sum(H[,i])
      H[,i] <- H[,i]/S.2
    }
  }
  
  
  # Save projected dataset
  
  V <- data.frame(H)
  names(V) <- sample.names
  row.names(V) <- W.names
  MP.WriteGct(gct.data.frame = V, filename = output.file)  
  
}
MP.ReadClsFile <- function(file = "NULL") { 
  #
  # Reads a class vector CLS file and defines phenotype and class labels vectors (numeric and character) for the samples in a gene expression file (RES or GCT format)
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  
  cls.cont <- readLines(file)
  cls.cont <- cls.cont[1:3]
  num.lines <- length(cls.cont)
  class.list <- unlist(strsplit(cls.cont[[3]], " "))
  s <- length(class.list)
  t <- table(class.list)
  l <- length(t)
  phen <- vector(length=l, mode="character")
  class.v <- vector(length=s, mode="numeric")
  
  current.label <- class.list[1]
  current.number <- 1
  class.v[1] <- current.number
  phen[1] <- current.label
  phen.count <- 1
  
  if (length(class.list) > 1) {
    for (i in 2:s) {
      if (class.list[i] == current.label) {
        class.v[i] <- current.number
      } else {
        phen.count <- phen.count + 1
        current.number <- current.number + 1
        current.label <- class.list[i]
        phen[phen.count] <- current.label
        class.v[i] <- current.number
      }
    }
  }
  return(list(phen = phen, class.v = class.v, class.list = class.list))
}
MP.Match.and.Select <- function(
  input1.ds,
  input2.ds,
  output.ds) {
  
  # Match the genes of the first dataset on the second and select those rows from the second
  
  # start of methodology
  
  print(c("Running MP.Match.and.Select... on: ", input1.ds, " ", input2.ds))
  
  # Read input datasets
  
  dataset1 <- MP.Gct2Frame(filename = input1.ds)
  m1 <- data.matrix(dataset1$ds)
  gs.names1 <- dataset1$row.names
  gs.descs1 <- dataset1$descs
  sample.names1 <- dataset1$names
  
  dataset2 <- MP.Gct2Frame(filename = input2.ds)
  m2 <- data.matrix(dataset2$ds)
  gs.names2 <- dataset2$row.names
  gs.descs2 <- dataset2$descs
  sample.names2 <- dataset2$names
  
  # Match features to first dataset and create matching m2 dataset
  
  gs.names3 <- intersect(gs.names1, gs.names2)
  
  locations2 <- match(gs.names3, gs.names2, nomatch=0)
  gs.names2 <- gs.names2[locations2]
  gs.descs2 <- gs.descs2[locations2]
  m2 <- m2[locations2, ]
  
  # Save dataset
  
  V <- data.frame(m2)
  names(V) <- sample.names2
  row.names(V) <- gs.names2
  MP.WriteGct(gct.data.frame = V, descs = gs.descs2, filename = output.ds)  
  
}
MP.Match.and.Merge <- function(
  input1.ds,
  input1.cls = "",
  input2.ds,
  input2.cls = "",
  output.ds,
  output.cls) {
  
  # start of methodology
  
  print(c("Running MP.Match.and.Merge... on: ", input1.ds, " ", input2.ds))
  
  # Read input datasets
  
  dataset1 <- MP.Gct2Frame(filename = input1.ds)
  m1 <- data.matrix(dataset1$ds)
  gs.names1 <- dataset1$row.names
  gs.descs1 <- dataset1$descs
  sample.names1 <- dataset1$names
  
  dataset2 <- MP.Gct2Frame(filename = input2.ds)
  m2 <- data.matrix(dataset2$ds)
  gs.names2 <- dataset2$row.names
  gs.descs2 <- dataset2$descs
  sample.names2 <- dataset2$names
  
  # Read CLS files 
  
  if ((input1.cls != "") & (input2.cls != "")) {
    CLS1 <- MP.ReadClsFile(file=input1.cls)
    class.labels1 <- CLS1$class.v
    class.phen1 <- CLS1$phen
    
    CLS2 <- MP.ReadClsFile(file=input2.cls)
    class.labels2 <- CLS2$class.v
    class.phen2 <- CLS2$phen
  }
  
  # Match features to first dataset and create matching m2 dataset
  
  gs.names3 <- intersect(gs.names1, gs.names2)
  
  locations1 <- match(gs.names3, gs.names1, nomatch=0)
  m1 <- m1[locations1, ]
  gs.descs1 <- gs.descs1[locations1]
  
  locations2 <- match(gs.names3, gs.names2, nomatch=0)
  m2 <- m2[locations2, ]
  gs.descs2 <- gs.descs2[locations2]
  
  # Merge datasets
  
  m3 <- cbind(m1, m2)
  sample.names3 <- c(sample.names1, sample.names2)
  
  if ((input1.cls != "") & (input2.cls != "")) {
    class.labels3 <- c(class.labels1, class.labels2 + length(class.phen1))
    class.phen3 <- c(class.phen1, class.phen2)
  }
  
  # Save datasets
  
  V <- data.frame(m3)
  names(V) <- sample.names3
  row.names(V) <- gs.names3
  MP.WriteGct(gct.data.frame = V, descs = gs.descs1, filename = output.ds)  
  
  if ((input1.cls != "") & (input2.cls != "")) {
    MP.WriteCls(class.v = class.labels3, phen = class.phen3, filename = output.cls) 
  }
}

MP.Gct2Frame <- function(filename = "NULL") { 
  #
  # Reads a gene expression dataset in GCT format and converts it into an R data frame
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  
  ds <- read.delim(filename, header=T, sep="\t", skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T, na.strings = "") # row.names = 1 is supposed to be it
  descs <- ds[,1]
  ds <- ds[-1]
  row.names <- row.names(ds)
  names <- names(ds)
  return(list(ds = ds, row.names = row.names, descs = descs, names = names))
}

MP.VarFilter <- function(V, fold, delta, gene.names = "", gene.descs = "") { 
  
  # Variation filter pre-processing for gene expression matrix
  
  cols <- length(V[1,])
  rows <- length(V[,1])
  row.max <- apply(V, MARGIN=1, FUN=max)
  row.min <- apply(V, MARGIN=1, FUN=min)
  flag <- array(dim=rows)
  flag <- (row.max /row.min >= fold) & (row.max - row.min >= delta)
  size <- sum(flag)
  B <- matrix(0, nrow = size, ncol = cols)
  j <- 1
  if (length(gene.names) == 1) {
    for (i in 1:rows) {
      if (flag[i]) {
        B[j,] <- V[i,]
        j <- j + 1
      }
    }
    return(B)
  } else {
    new.gene.names <- vector(mode = "character", length = size)
    new.gene.descs <- vector(mode = "character", length = size)
    for (i in 1:rows) {
      if (flag[i]) {
        B[j,] <- V[i,]
        new.gene.names[j] <- gene.names[i]
        new.gene.descs[j] <- gene.descs[i]
        j <- j + 1
      }
    }
    return(list(V = B, new.gene.names = new.gene.names, new.gene.descs = new.gene.descs, locations = flag))
  }
}

MP.NormalizeCols.Rank <- function(V) { 
  
  cols <- length(V[1,])
  rows <- length(V[,1])
  for (j in 1:cols) {  # column rank normalization
    V[,j] <- rank(V[,j], ties.method = "average")
  }
  
  return(V)
}

MP.NormalizeCols <- function(V) { 
  #
  # Stardardize columns of a gene expression matrix
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  
  col.mean <- apply(V, MARGIN=2, FUN=mean)
  col.sd <- apply(V, MARGIN=2, FUN=sd)
  col.n <- length(V[1,])
  for (i in 1:col.n) {
    if (col.sd[i] == 0) {
      V[i,] <- 0
    } else {
      V[,i] <- (V[,i] - col.mean[i])/col.sd[i]
    }
  }
  return(V)
}

MP.NormalizeCols.Rank <- function(V) { 
  #
  cols <- length(V[1,])
  rows <- length(V[,1])
  for (j in 1:cols) {  # column rank normalization
    V[,j] <- rank(V[,j], ties.method = "average")
  }
  
  return(V)
}
MP.Threshold <- function(V, thres, ceil) { 
  #
  # Threshold and ceiling pre-processing for gene expression matrix
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  
  V[V < thres] <- thres
  V[V > ceil] <- ceil
  return(V)
}

MP.NormalizeCols.Rescale <- function(V) { 
  
  epsilon <- 0.00001
  cols <- length(V[1,])
  for (j in 1:cols) {  # column rank normalization
    max.v <- max(V[,j])
    min.v <- min(V[,j])
    V[,j] <- (V[,j] - min.v + epsilon)/(max.v - min.v)
  }
  
  return(V)
}

MP.WriteGct <- function (gct.data.frame, descs = "", filename) 
{
  f <- file(filename, "w")
  cat("#1.2", "\n", file = f, append = TRUE, sep = "")
  cat(dim(gct.data.frame)[1], "\t", dim(gct.data.frame)[2], "\n", file = f, append = TRUE, sep = "")
  cat("Name", "\t", file = f, append = TRUE, sep = "")
  cat("Description", file = f, append = TRUE, sep = "")
  
  names <- names(gct.data.frame)
  cat("\t", names[1], file = f, append = TRUE, sep = "")
  
  if (length(names) > 1) {
    for (j in 2:length(names)) {
      cat("\t", names[j], file = f, append = TRUE, sep = "")
    }
  }
  cat("\n", file = f, append = TRUE, sep = "\t")
  
  oldWarn <- options(warn = -1)
  m <- matrix(nrow = dim(gct.data.frame)[1], ncol = dim(gct.data.frame)[2] +  2)
  m[, 1] <- row.names(gct.data.frame)
  if (length(descs) > 1) {
    m[, 2] <- descs
  } else {
    m[, 2] <- row.names(gct.data.frame)
  }
  index <- 3
  for (i in 1:dim(gct.data.frame)[2]) {
    m[, index] <- gct.data.frame[, i]
    index <- index + 1
  }
  write.table(m, file = f, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", col.names = FALSE, row.names = FALSE)
  close(f)
  options(warn = 0)
  
}

MP.WriteCls <- function (class.v, phen, filename) 
{
  f <- file(filename, "w")
  n <- length(phen)
  l <- length(class.v)
  cat(l, n, "1", "\n", file = f, append = TRUE, sep = " ")
  cat("#", phen, "\n", file = f, append = TRUE, sep = " ")
  class.v <- phen[class.v]
  cat(class.v, "\n", file = f, append = TRUE, sep = " ")
  close(f)
}

NMF.div <- function(V, k, maxniter = 2000, seed = 123456, stopconv = 40, stopfreq = 10) {
  
  N <- length(V[,1])
  M <- length(V[1,])
  set.seed(seed)
  W <- matrix(runif(N*k), nrow = N, ncol = k)  # Initialize W and H with random numbers
  H <- matrix(runif(k*M), nrow = k, ncol = M)
  VP <- matrix(nrow = N, ncol = M)
  error.v <- vector(mode = "numeric", length = maxniter)
  new.membership <- vector(mode = "numeric", length = M)
  old.membership <- vector(mode = "numeric", length = M)
  no.change.count <- 0
  eps <- .Machine$double.eps
  for (t in 1:maxniter) {
    VP = W %*% H
    W.t <- t(W)
    H <- H * (W.t %*% (V/VP)) + eps
    norm <- apply(W, MARGIN=2, FUN=sum)
    for (i in 1:k) {
      H[i,] <- H[i,]/norm[i]
    }
    VP = W %*% H
    H.t <- t(H)
    W <- W * ((V/VP) %*% H.t) + eps
    norm <- apply(H, MARGIN=1, FUN=sum)
    for (i in 1:k) {
      W[,i] <- W[,i]/norm[i]
    }
    error.v[t] <- sum(V * log((V + eps)/(VP + eps)) - V + VP)/(M * N)
    if (t %% stopfreq == 0) {
      
      for (j in 1:M) {
        class <- order(H[,j], decreasing=T)
        new.membership[j] <- class[1]
      }
      if (sum(new.membership == old.membership) == M) {
        no.change.count <- no.change.count + 1
      } else {
        no.change.count <- 0
      }
      if (no.change.count == stopconv) break
      old.membership <- new.membership
    }
  }
  return(list(W = W, H = H, t = t, error.v = error.v))
}

NMF <- function(V, k, maxniter = 2000, seed = 123456, stopconv = 40, stopfreq = 10) {
  N <- length(V[,1])
  M <- length(V[1,])
  set.seed(seed)
  W <- matrix(runif(N*k), nrow = N, ncol = k)  # Initialize W and H with random numbers
  H <- matrix(runif(k*M), nrow = k, ncol = M)
  VP <- matrix(nrow = N, ncol = M)
  error.v <- vector(mode = "numeric", length = maxniter)
  new.membership <- vector(mode = "numeric", length = M)
  old.membership <- vector(mode = "numeric", length = M)
  eps <- .Machine$double.eps
  for (t in 1:maxniter) {
    VP = W %*% H
    H <- H * (crossprod(W, V)/crossprod(W, VP)) + eps
    VP = W %*% H
    H.t <- t(H)
    W <- W * (V %*% H.t)/(VP %*% H.t) + eps
    error.v[t] <- sqrt(sum((V - VP)^2))/(N * M)
    if (t %% stopfreq == 0) {
      for (j in 1:M) {
        class <- order(H[,j], decreasing=T)
        new.membership[j] <- class[1]
      }
      if (sum(new.membership == old.membership) == M) {
        no.change.count <- no.change.count + 1
      } else {
        no.change.count <- 0
      }
      if (no.change.count == stopconv) break
      old.membership <- new.membership
    }
  }
  return(list(W = W, H = H, t = t, error.v = error.v))
}

# Computes a Sparse NMF decomposition (Gao and Church 2005, Bioinformatics 2005 21(21):3970-3975)
# Adapted from an original C++ program kindly provided by Yuan Gao 

SNMF <- function(V, k, maxniter = 2000, seed = 123456, lambda = 1, stopconv = 40, stopfreq = 10) {
  N <- length(V[,1])
  M <- length(V[1,])
  set.seed(seed)
  W <- matrix(runif(N*k), nrow = N, ncol = k)  # Initialize W and H with random numbers
  H <- matrix(runif(k*M), nrow = k, ncol = M)
  VP <- matrix(nrow = N, ncol = M)
  error.v <- vector(mode = "numeric", length = maxniter)
  new.membership <- vector(mode = "numeric", length = M)
  old.membership <- vector(mode = "numeric", length = M)
  eps <- .Machine$double.eps
  for (t in 1:maxniter) {
    T <- crossprod(W, W) + lambda * diag(k)
    for (j in 1:M) {
      b <- crossprod(W, V[,j])
      H[,j] <- solve(T, b)
    }
    H[ H < 0] <- 0
    H.t <- t(H)
    VH <- V %*% H.t
    HH <- H %*% H.t
    WH <- W %*% HH
    W <- W * VH/WH + eps
    VP <- W %*% H
    error.v[t] <- sqrt(sum((V - VP)^2))/(N * M)
    if (t %% stopfreq == 0) {
      for (j in 1:M) {
        class <- order(H[,j], decreasing=T)
        new.membership[j] <- class[1]
      }
      if (sum(new.membership == old.membership) == M) {
        no.change.count <- no.change.count + 1
      } else {
        no.change.count <- 0
      }
      if (no.change.count == stopconv) break
      
      old.membership <- new.membership
    }
  }
  return(list(W = W, H = H, t = t, error.v = error.v))
}

NSNMF.div <- function(V, k, theta, maxniter = 2000, seed = 123456, stopconv = 40, stopfreq = 10) {
  
  # Non-smooth NMF from Carmona-Saez et al 2006 BMC Bioinformatics. 2006; 7: 78.
  
  N <- length(V[,1])
  M <- length(V[1,])
  set.seed(seed)
  W <- matrix(runif(N*k), nrow = N, ncol = k)  # Initialize W and H with random numbers
  H <- matrix(runif(k*M), nrow = k, ncol = M)
  VP <- matrix(nrow = N, ncol = M)
  error.v <- vector(mode = "numeric", length = maxniter)
  new.membership <- vector(mode = "numeric", length = M)
  old.membership <- vector(mode = "numeric", length = M)
  no.change.count <- 0
  eps <- .Machine$double.eps
  
  # smoothing matrix
  SS <- theta/k * (rep(1, k) %*% t(rep(1, k)))
  for (i in 1:k) {
    SS[i, i] <- SS[i, i] + (1 - theta)
  }
  
  for (t in 1:maxniter) {
    W <- W %*% SS
    VP = W %*% H
    W.t <- t(W)
    H <- H * (W.t %*% (V/VP)) + eps
    norm <- apply(W, MARGIN=2, FUN=sum)
    for (i in 1:k) {
      H[i,] <- H[i,]/norm[i]
    }
    H <- SS %*% H
    VP = W %*% H
    H.t <- t(H)
    W <- W * ((V/VP) %*% H.t) + eps
    norm <- apply(H, MARGIN=1, FUN=sum)
    for (i in 1:k) {
      W[,i] <- W[,i]/norm[i]
    }
    error.v[t] <- sum(V * log((V + eps)/(VP + eps)) - V + VP)/(M * N)
    if (t %% stopfreq == 0) {
      
      for (j in 1:M) {
        class <- order(H[,j], decreasing=T)
        new.membership[j] <- class[1]
      }
      if (sum(new.membership == old.membership) == M) {
        no.change.count <- no.change.count + 1
      } else {
        no.change.count <- 0
      }
      if (no.change.count == stopconv) break
      old.membership <- new.membership
    }
  }
  return(list(W = W, H = H, t = t, error.v = error.v))
}

MP.Projection.Plots <- function(
  input.ds, 
  input.cls = "", 
  model.set = "ALL",
  output.2D.proj.plot, 
  output.heatmap.plot,
  output.heatmap.sorted.plot,
  title = "",
  seed = 1234, 
  heatmap.row.norm = T,
  heatmap.cmap.type = 1,
  symbol.scaling = 1,
  col = c("greny3", "cadetblue", "darkgreen", "chartreuse2", "red1", "darkred", "orange", "blue2", "lightblue", "pink", "coral"),
  symbs = c(22, 21, 20, 23, 24, 25, 21, 20, 23, 24, 25)         
) {
  
  print(c("Running MP.Projection.Plots... on:", input.ds))
  
  library("scatterplot3d")
  library(MASS)
  
  set.seed(seed=seed, kind = NULL)
  
  # Read dataset and initialize plot filenames
  
  dataset <- MP.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names
  plot.filename <- paste(output.2D.proj.plot, ".pdf", sep="", collapse="")
  pdf(file=plot.filename, height = 20, width = 30)
  
  dim(m) 
  
  Ns <- length(m[1,])
  k.proj <- length(m[,1])
  
  if (input.cls != "") {
    CLS <- MP.ReadClsFile(file=input.cls)
    class.labels <- CLS$class.v
    class.phen <- CLS$phen
  } else {
    class.labels <- rep(1, Ns)
    class.phen <- "Samples"
  }
  
  # Separate data into train and test pieces
  
  if (model.set == "ALL") {
    model.set <- seq(1, Ns)
  }
  
  m.train <- as.matrix(m[, model.set])
  num.samples.train <- length(model.set)
  sample.names.train <- sample.names[model.set]
  if (input.cls != "") {
    class.labels.train <- class.labels[model.set]
  }
  m.test <- as.matrix(m[, - model.set])
  sample.names.test <- sample.names[- model.set]
  if (input.cls != "") {
    class.labels.test <- class.labels[- model.set]
  }
  
  pca <- prcomp(t(m.train), retx = TRUE, center = TRUE, scale. = TRUE)
  
  S1 <- pca$x[,1]
  S2 <- pca$x[,2]
  
  X1 <- pca$rotation[,1]
  X2 <- pca$rotation[,2]
  
  # 2D plots
  
  max.S <- max(sqrt(S1*S1 + S2*S2))
  max.X <- max(sqrt(X1*X1 + X2*X2))
  X1 <-  max.S * X1/max.X
  X2 <-  max.S * X2/max.X
  max.A <- max(max.S, max.X)
  
  
  c0 <- col
  c1 <- col
  #   c1 <- colors()[match(c0, colors())]
  color <- c1[class.labels]
  
  if (.Platform$OS.type == "windows") {
    plot.filename <- output.2D.proj.plot
    windows(height = 20, width = 30)
  }
  
  nf <- layout(matrix(c(1, 2, 3), 1, 3, byrow=T), widths = c(3, 3, 1), heights = 1, respect = FALSE)
  
  # 1st subplot 
  
  plot(S1, S2, xlim = c(-max.A, max.A), ylim = c(-max.A, max.A), type = "n", main = paste(title, " -- Model Samples Biplot", sep=""), sub = input.ds)
  
  for (j in 1:num.samples.train) {
    if (min(class.labels) == 0) {
      color.code <- c1[class.labels[j] + 1]
      symb <- symbs[class.labels[j] + 1]
    } else {
      color.code <- c1[class.labels[j]]
      symb <- symbs[class.labels[j]]
    }
    points(S1[j], S2[j], pch=symb, type="p", cex = symbol.scaling*3, bg = color.code, col = "black")   
  }
  
  
  for (j in 1:k.proj) {
    x.coor <- X1[j]*0.925
    y.coor <- X2[j]*0.925
    arrows(0, 0, x.coor, y.coor, lwd = 1, length = 0.15, angle = 20, col = "black")
    
    leg.txt <- paste("F", j, sep = "")
    
    text(X1[j], X2[j], labels = leg.txt, adj = NULL, pos = NULL, offset = 0.5, vfont = NULL, cex = symbol.scaling*2, col = "black")
  }
  
  ang <- vector(length = k.proj, mode = "numeric")
  for (j in 1:k.proj) {
    ang[j] <- ifelse(atan2(X2[j], X1[j]) > 0,  atan2(X2[j], X1[j]), 2*pi + atan2(X2[j], X1[j]))
  }
  
  ang.index <- order(ang, decreasing=F)
  ang2 <- ang[ang.index]
  
  for (j in 1:k.proj) {
    if (j == k.proj) {
      angle.in.between <- (ang2[1] - ang2[j] - 2*pi)/2 + ang2[j] - 2*pi
    } else {
      angle.in.between <- (ang2[j + 1] - ang2[j])/2 + ang2[j]
    }
    x <- max.S * cos(angle.in.between)
    y <- max.S * sin(angle.in.between)
    arrows(0, 0, x, y, lwd = 4, length = 0, lty = 3, col = "grey50")
  }
  
  # 2nd subplot Project all data
  
  test.scores <- predict(pca, t(m.test))
  
  S1 <- c(pca$x[,1], test.scores[,1])
  S2 <- c(pca$x[,2], test.scores[,2])
  
  max.S <- max(sqrt(S1*S1 + S2*S2))
  max.X <- max(sqrt(X1*X1 + X2*X2))
  X1 <-  max.S * X1/max.X
  X2 <-  max.S * X2/max.X
  num.samples <- length(S1)
  
  plot(S1, S2, xlim = c(-max.A, max.A), ylim = c(-max.A, max.A), type = "n", main = paste(title, " -- Model + Test Samples Biplot", sep=""), sub = input.ds)
  
  for (j in 1:num.samples) {
    if (min(class.labels) == 0) {
      symb <- symbs[class.labels[j] + 1]
      color.code <- c1[class.labels[j] + 1]
    } else {
      symb <- symbs[class.labels[j]]
      color.code <- c1[class.labels[j]]
    }
    points(S1[j], S2[j], pch=symb, type="p", cex = symbol.scaling*3, bg = color.code, col = "black")
  }
  
  for (j in 1:k.proj) {
    x.coor <- X1[j]*0.925
    y.coor <- X2[j]*0.925
    arrows(0, 0, x.coor, y.coor, lwd = 1, length = 0.15, angle = 20, col = "black")
    
    leg.txt <- paste("F", j, sep = "")
    
    text (X1[j], X2[j], labels = leg.txt, adj = NULL, pos = NULL, offset = 0.5, vfont = NULL, cex = symbol.scaling*2, col = "black")
  }
  
  ang <- vector(length = k.proj, mode = "numeric")
  for (j in 1:k.proj) {
    ang[j] <- ifelse(atan2(X2[j], X1[j]) > 0,  atan2(X2[j], X1[j]), 2*pi + atan2(X2[j], X1[j]))
  }
  
  ang.index <- order(ang, decreasing=F)
  ang2 <- ang[ang.index]
  
  for (j in 1:k.proj) {
    if (j == k.proj) {
      angle.in.between <- (ang2[1] - ang2[j] - 2*pi)/2 + ang2[j] - 2*pi
    } else {
      angle.in.between <- (ang2[j + 1] - ang2[j])/2 + ang2[j]
    }
    x <- max.S * cos(angle.in.between)
    y <- max.S * sin(angle.in.between)
    arrows(0, 0, x, y, lwd = 4, length = 0, lty = 3, col = "grey50")
  }
  
  
  # 3nd subplot: legend 
  
  class.phen.train <- unique(class.labels.train)
  leg.txt <- class.phen
  n.phen <- length(class.phen)
  p.vec <- symbs[1:n.phen]
  c.vec <- c1[1:n.phen]
  
  par(mar = c(0, 0, 0, 0))
  plot(c(0,0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes=F, type="n", xlab = "", ylab="")
  legend(x=0, y=1, legend=leg.txt, bty="n", xjust=0, yjust= 1, pch = p.vec, pt.bg = c.vec, col = "black", cex = symbol.scaling*1.35, pt.cex=symbol.scaling*3, x.intersp = 2, y.intersp = 5)
  
  if (.Platform$OS.type == "windows") {
    savePlot(filename = plot.filename, type ="jpeg", device = dev.cur())
  } else {
    dev.off()
  }
  
  # Heat map plot
  plot.filename <- paste(output.heatmap.plot, ".pdf", sep="", collapse="")
  height <- ifelse(k.proj > 50, 20, 0.50*k.proj + 7)
  pdf(file=plot.filename, height = height, width = 25)
  
  if (.Platform$OS.type == "windows") {
    plot.filename <- output.heatmap.plot
  }
  
  MP.HeatMapPlot(V = m, row.names = gs.names, col.labels = class.labels, col.classes = class.phen, col.names = sample.names, main= paste(title, " ", sep=""), sub = " ", xlab=" ", ylab=" ", row.norm = heatmap.row.norm,  cmap.type = heatmap.cmap.type, char.rescale = 1) 
  
  # legend
  
  if (.Platform$OS.type == "windows") {
    savePlot(filename = plot.filename, type ="jpeg", device = dev.cur())
  } else {
    dev.off()
  }
  
  # heat map sorted inside each phenotype 
  
  m2 <- m
  gs.names2 <- gs.names
  sample.names2 <- sample.names
  max.classes <- max(class.labels)
  plot.filename <- paste(output.heatmap.sorted.plot, ".pdf", sep="", collapse="")
  
  for (k in 1:max.classes) {
    if (sum(class.labels==k) > 1) {
      m3 <- m2[,class.labels==k]   
      sn <- sample.names2[class.labels==k]
      dist.matrix <- dist(t(m3))
      HC <- hclust(dist.matrix, method="complete")
      m3 <- m3[, HC$order]
      sn <- sn[HC$order]
      m2[,class.labels==k] <- m3
      sample.names2[class.labels==k] <- sn
    }
  }
  
  height <- ifelse(k.proj > 50, 20, 0.50*k.proj + 7)
  
  if (.Platform$OS.type == "windows") {
    plot.filename <- output.heatmap.sorted.plot
  } else {
    pdf(file=plot.filename, height = height, width = 25)
  }
  
  MP.HeatMapPlot(V = m2, row.names = gs.names2, col.labels = class.labels, col.classes = class.phen, col.names = sample.names2, main= paste(title, " (sorted inside each class)", sep=""), sub = " ", xlab=" ", ylab=" ", row.norm = heatmap.row.norm,  cmap.type = heatmap.cmap.type, char.rescale = 1) 
  
  if (.Platform$OS.type == "windows") {
    savePlot(filename = plot.filename, type ="jpeg", device = dev.cur())
  } else {
    dev.off()
  }
  
}

MP.replace.row.col <- function(input.ds, output.ds, mode = "row", number, values) {
  
  dataset <- MP.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  
  
  if (mode == "row") {
    m[number, ] <-  ifelse (length(values) == 1, rep(values, length(m[number,])), values)
  } else if (mode == "col") {
    m[, number] <-  ifelse (length(values) == 1, rep(values, length(m[, number])), values)
  } else {
    stop(c("unknown mode:", mode))
  }
  
  V <- data.frame(m)
  names(V) <- dataset$names
  row.names(V) <- dataset$row.names
  
  MP.WriteGct(gct.data.frame = V, descs = dataset$descs, filename = output.ds)
}


MP.Resort.Rows <- function(
  input.ds,
  new.order,
  new.row.labs = NULL,
  new.row.descs = NULL,                             
  output.ds) {
  
  # Read input datasets
  
  dataset <- MP.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names
  
  m <- m[new.order,] 
  
  if (length(new.row.labs) == 1) {
    gs.names <- gs.names[new.order]
    gs.descs <- gs.descs[new.order]
  } else {
    gs.names <- new.row.labs
    gs.names <- new.row.descs
  }
  
  # Save dataset
  
  V <- data.frame(m)
  names(V) <- sample.names
  row.names(V) <- gs.names
  MP.WriteGct(gct.data.frame = V, descs = gs.descs, filename = output.ds)  
  
}

MP.consensusNMF <- function(input.ds, k.init, k.final, num.clusterings, maxniter, error.function, rseed=123456789, directory = "", stopconv = 40, stopfreq = 10, non.interactive.run = F, doc.string = "", ...) {
  
  #
  #  GenePattern Methodology for:
  #
  #  Metagenes and Molecular Pattern Discovery using Matrix Factorization
  #  Jean-Philippe Brunet, Pablo Tamayo, Todd. R. Golub, and Jill P. Mesirov
  # 
  #  Author:  Pablo Tamayo (tamayo@genome.wi.mit.edu)
  #
  #  Based on the original matlab version written by Jean-Philippe Brunet (brunet@broad.mit.edu) and
  #  with additional contributions from: Ted Liefeld (liefeld@broad.mit.edu)   
  #  Date:  November 27, 2003
  #
  #  Last change March 3, 2005: modifications to make the output more readable.
  #
  #  Execute from an R console window with this command:
  #  source("<this file>", echo = TRUE)
  #  E.g. someoutput <- mynmf2(input.ds="c:\\nmf\\all_aml.res",k.init=2,k.final=5,num.clusterings=20,maxniter=500) 
  #
  #  For details on the method see:
  #
  #  Proc. Natl. Acad. Sci. USA 2004 101: 4164-4169
  #  http://www.broad.mit.edu/cgi-bin/cancer/publications/pub_paper.cgi?mode=view&paper_id=89
  #
  #  Input parameters
  #
  #   input.ds
  #                       input gene expression dataset in GCT or RES format
  #   k.init
  #                       initial value of k
  #   k.final
  #                       final value of k
  #   num.clusterings
  #                       number of NMF clusterings to build consensus matrix
  #   maxniter
  #                       maximum number of NMF iterations
  #   error.function
  #                       NMF error function: "divergence" of "euclidean"
  #   rseed
  #                       random number generator seed
  #   directory
  #                       file directory where to store the result files
  #   stopconv
  #                       how many no change checks are needed to stop NMF iterations (convergence)
  #   stopfreq
  #                       frequency (NMF iterations) of "no change" checks 
  #   non.interactive.run 
  #                       flag controlling if the plots are produced interatively (Rgui and saved) or only saved in files
  #   doc.string
  #                       prefix to be added to the output files
  #
  #  Output files are (prefix with "doc.string")
  #
  #   params.txt 
  #                       run parameters and time of execution
  #   membership.gct		
  #			membership results for samples at all values of K
  #   cophenetic.txt 
  #			cophenetic values for each K
  #   cophenetic.plot.jpeg
  #			plot of cophenetic for each value of K		
  #   consensus.k#.gct (for each value of K)
  #			consensus matrix for k=#
  #   consensus.plot.k#.jpeg (for each value of K)
  #			plot of consensus matrix for k=#
  #   graphs.k#.jpeg (for each value of K)
  
  # save input parameters
  
  filename <- paste(directory, doc.string, ".params.txt", sep="", collapse="")  
  
  time.string <- as.character(as.POSIXlt(Sys.time(),"GMT"))
  write(paste("Run of NMF on ", time.string), file=filename)
  
  write(paste("input.ds =", input.ds, sep=" "), file=filename, append=T) 
  write(paste("k.init = ", k.init, sep=" "), file=filename, append=T) 
  write(paste("k.final =", k.final, sep=" "), file=filename, append=T) 
  write(paste("num.clusterings =", num.clusterings, sep=" "), file=filename, append=T) 
  write(paste("maxniter =", maxniter, sep=" "), file=filename, append=T) 
  write(paste("error.function =", error.function, sep=" "), file=filename, append=T) 
  write(paste("rseed =", rseed, sep=" "), file=filename, append=T) 
  write(paste("directory =", directory, sep=" "), file=filename, append=T) 
  write(paste("stopconv =", stopconv, sep=" "), file=filename, append=T) 
  write(paste("stopfreq =", stopfreq, sep=" "), file=filename, append=T)
  write(paste("non.interctive.run =", non.interactive.run, sep=" "), file=filename, append=T) 
  write(paste("doc.string =", doc.string, sep=" "), file=filename, append=T) 
  
  
  k.init<-as.integer(k.init)
  k.final<-as.integer(k.final)
  num.clusterings<-as.integer(num.clusterings)
  n.iter<-as.integer(maxniter)
  if (!is.na(rseed)){
    seed <- as.integer(rseed)
  }
  
  
  # library(mva)
  # library(MASS)
  # library(GenePattern)
  
  D <- read.dataset(input.ds)
  A <- data.matrix(D)
  
  # Threshold negative values to small quantity 
  
  eps <- .Machine$double.eps
  A[A < 0] <- eps
  
  
  
  cols <- length(A[1,])
  rows <- length(A[,1])
  
  col.names <- names(D)
  
  num.k <- k.final - k.init + 1
  
  rho <- vector(mode = "numeric", length = num.k)
  k.vector <- vector(mode = "numeric", length = num.k)
  
  k.index <- 1
  
  connect.matrix.ordered <- array(0, c(num.k, cols, cols))
  
  for (k in k.init:k.final) { 
    
    if (non.interactive.run == F) {
      if (.Platform$OS.type == "windows") {
        filename <- paste(directory, doc.string, ".", "graphs.k", k, sep="", collapse="")
        #             windows(width = 9, height = 11)
        windows(width = 40, height = 22)
      } else if (.Platform$OS.type == "unix") {
        filename <- paste(directory, doc.string, ".", "graphs.k", k, ".pdf", sep="", collapse="")
        pdf(file=filename, width = 9, height = 11)
      }
    } else {
      if (.Platform$OS.type == "unix") {
        filename <- paste(directory, doc.string, ".", "graphs.k", k, ".pdf", sep="", collapse="")
        pdf(file=filename, width = 9, height = 11)
      } else if (.Platform$OS.type == "windows") {
        filename <- paste(directory, doc.string, ".", "graphs.k", k, ".pdf", sep="", collapse="")
        pdf(file=filename, width = 9, height = 11)
        dev.off()
      }
    }
    
    #   nf <- layout(matrix(c(1,2,3,4,5,6,7,8), 4, 2, byrow=T), c(1, 1, 1, 1), c(1, 1), TRUE)
    nf <- layout(matrix(c(1,2,3,4,5,6,7,8), 2, 4, byrow=T), c(1, 1), c(1, 1, 1, 1), TRUE)
    assign <- matrix(0, nrow = num.clusterings, ncol = cols)
    
    for (i in 1:num.clusterings) {
      
      print(paste("Computing clustering number=", i, " for k=", k, sep=""))
      
      if (error.function == "divergence"){
        NMF.out <- NMF.div(V = A, k = k, maxniter = n.iter, seed = seed + i, stopconv = stopconv, stopfreq = stopfreq)
      } else if (error.function == "euclidean"){
        NMF.out <- NMF(V = A, k = k, maxniter = n.iter, seed = seed + i, stopconv = stopconv, stopfreq = stopfreq)
      } else {
        stop(paste("Un-supported error function=", error.function, sep=""))
      }
      print(paste(NMF.out$t, " NMF iterations performed", sep=""))
      
      for (j in 1:cols) { # Find membership
        class <- order(NMF.out$H[,j], decreasing=T)
        assign[i, j] <- class[1]
      }
      
      if (i == 1) {  # Plot example for first clustering iteration
        H.saved <- NMF.out$H
        sub.string <- paste(doc.string, " k=", k, sep="")
        plot(1:NMF.out$t, NMF.out$error.v[1:NMF.out$t], pch = 20, cex = 1.5, col = 1, xlab="time", ylab="NMF error", sub=sub.string, main=paste("Convergence plot k=", k, " example", sep=""))
        
        
        if (rows < 1000) {
          W <- NMF.out$W
        } else {
          W <- NMF.out$W[sample(x = 1:rows, size = 1000),]
        }
        sub.string <- paste(doc.string, " k=", k, sep="")
        matrix.abs.plot(W, sub = sub.string, log = F, main = "Example W matrix (orig. order)", ylab = "genes", xlab ="metasamples")
        matrix.abs.plot(H.saved, sub = sub.string, log = F, main = "Example H matrix (orig. order)", ylab = "metagenes", xlab ="samples")
        metagene.plot(H = H.saved, main = "Metagenes Example (orig. order)", sub = sub.string, xlab = "samples", ylab = "metagenes")
        
      }
      
      rm(NMF.out)
      
    }  ## end  for (i in 1:num.clusterings)
    
    
    # compute consensus matrix
    connect.matrix <- matrix(0, nrow = cols, ncol = cols)
    
    for (i in 1:num.clusterings) {
      for (j in 1:cols) {
        for (p in 1:cols) {
          if (j != p) {
            if (assign[i, j] == assign[i, p]) {
              connect.matrix[j, p] <- connect.matrix[j, p] + 1
            } 
          } else {
            connect.matrix[j, p] <- connect.matrix[j, p] + 1
          }
        }
      }
    }
    
    connect.matrix <- connect.matrix / num.clusterings
    
    dist.matrix <- 1 - connect.matrix
    dist.matrix <- as.dist(dist.matrix)
    HC <- hclust(dist.matrix, method="average")
    
    dist.coph <- cophenetic(HC)
    k.vector[k.index] <- k
    rho[k.index] <- cor(dist.matrix, dist.coph)
    rho[k.index] <- signif(rho[k.index], digits = 4)
    
    #     connect.matrix.ordered <- matrix(0, nrow=cols, ncol = cols)
    
    for (i in 1:cols) {
      for (j in 1:cols) {
        connect.matrix.ordered[k.index, i, j] <- connect.matrix[HC$order[i], HC$order[j]]
      }
    }
    
    # compute consensus clustering membership
    
    membership <- cutree(HC, k = k)
    
    max.k <- max(membership)
    items.names.ordered <- col.names[HC$order]
    membership.ordered <- membership[HC$order]
    results <- data.frame(cbind(membership.ordered, items.names.ordered))
    
    if (k > k.init){
      all.membership <- cbind(all.membership, membership);
    } else {
      all.membership <- cbind(membership);
    }
    
    sub.string <- paste(doc.string, " k=", k, sep="")
    matrix.abs.plot(connect.matrix.ordered[k.index,,], sub=sub.string, log = F, main = "Ordered Consensus Matrix", ylab = "samples", xlab ="samples")
    plot(HC, xlab="samples", cex = 0.75, labels = col.names, sub = sub.string, col = "blue", main = paste("Ordered Linkage Tree. Coph=", rho[k.index]))
    
    resultsGct <- data.frame(membership.ordered)
    row.names(resultsGct) <- items.names.ordered
    filename <- paste(directory, doc.string, ".", "consensus.k.",k, ".gct", sep="", collapse="")
    MP.WriteGCt(resultsGct, filename)
    
    H.sorted <- H.saved[,HC$order]
    sub.string <- paste(doc.string, " k=", k, sep="")
    matrix.abs.plot(H.sorted, sub = sub.string, log = F, main = "Example H matrix (ordered)", ylab = "metagenes", xlab ="samples")
    metagene.plot(H = H.sorted, sub = sub.string, main = "Metagenes Example (ordered)", xlab = "samples", ylab = "metagenes")
    
    if (non.interactive.run == F) {  
      if (.Platform$OS.type == "windows") {
        savePlot(filename = filename, type ="jpeg", device = dev.cur())
      } else if (.Platform$OS.type == "unix") {
        dev.off()
      }
    } else {
      dev.off()
    }
    
    if (non.interactive.run == F) {
      if (.Platform$OS.type == "windows") {
        filename <- paste(directory, doc.string, ".", "consensus.plot.k", k, sep="", collapse="")
        #             windows(width = 8.5, height = 11)
        windows(width = 20, height = 20)
      } else if (.Platform$OS.type == "unix") {
        filename <- paste(directory, doc.string, ".", "consensus.plot.k", k, ".pdf", sep="", collapse="")
        pdf(file=filename, width = 8.5, height = 11)
        dev.off()
      }
    } else {
      if (.Platform$OS.type == "unix") {
        filename <- paste(directory, doc.string, ".", "consensus.plot.k", k, ".pdf", sep="", collapse="")
        pdf(file=filename, width = 8.5, height = 11)
      } else if (.Platform$OS.type == "windows") {
        filename <- paste(directory, doc.string, ".", "consensus.plot.k", k, ".pdf", sep="", collapse="")
        pdf(file=filename, width = 8.5, height = 11)
        dev.off()
      }
    }
    
    nf <- layout(matrix(c(1), 1, 1, byrow=T), c(1, 1), c(1, 1), TRUE)
    
    conlabel <- paste("Consensus k =", k, sep=" ", collapse="")
    
    sub.string <- paste("Consensus matrix k=", k, "; dataset= ", input.ds, sep="")
    ConsPlot(connect.matrix.ordered[k.index,,], col.labels = membership.ordered, col.names = items.names.ordered, main = " ", sub=sub.string, xlab=" ", ylab=" ")
    
    if (non.interactive.run == F) {  
      if (.Platform$OS.type == "windows") {
        savePlot(filename = filename, type ="jpeg", device = dev.cur())
      } else if (.Platform$OS.type == "unix") {
        dev.off()
      }
    } else {
      dev.off()
    }
    
    k.index <- k.index + 1
    
  } # end of loop over k
  
  
  # Save consensus matrices in one file
  
  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".", "consensus.all.k.plot", sep="")
      #             windows(width = 8.5, height = 11)
      windows(width = 30, height = 22)
    } else if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".", "consensus.all.k.plot.pdf", sep="")
      pdf(file=filename, width = 8.5, height = 11)
      dev.off()
    }
  } else {
    if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".", "consensus.all.k.plot.pdf", sep="")
      pdf(file=filename, width = 8.5, height = 11)
      dev.off()
    } else if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".", "consensus.all.k.plot.pdf", sep="")
      pdf(file=filename, width = 8.5, height = 11)
      dev.off()
    }
  }
  
  nf <- layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), 4, 4, byrow=T), c(1, 1, 1, 1), c(1, 1, 1, 1), TRUE)
  
  for (k in 1:num.k) { 
    matrix.abs.plot(connect.matrix.ordered[k,,], log = F, main = paste("k=", k.vector[k]), 
                    sub = paste("Cophenetic coef.=", rho[k]), ylab = "samples", xlab ="samples")
  }
  
  y.range <- c(1 - 2*(1 - min(rho)), 1)
  plot(k.vector, rho, main ="Cophenetic Coefficient", xlim=c(k.init, k.final), ylim=y.range, 
       xlab = "k", ylab="Cophenetic correlation", type = "n")
  lines(k.vector, rho, type = "l", col = "black")
  points(k.vector, rho, pch=22, type = "p", cex = 1.25, bg = "black", col = "black")
  
  if (non.interactive.run == F) {  
    if (.Platform$OS.type == "windows") {
      savePlot(filename = filename, type ="jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }
  
  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".", "cophenetic.plot", sep="")
      #             windows(width = 8.5, height = 11)
      windows(width = 20, height = 20)
    } else if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".", "cophenetic.plot.pdf", sep="")
      pdf(file=filename, width = 8.5, height = 11)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".", "cophenetic.plot.pdf", sep="")
      pdf(file=filename, width = 8.5, height = 11)
    } else if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".", "cophenetic.plot.pdf", sep="")
      pdf(file=filename, width = 8.5, height = 11)
      dev.off()
    }
  }
  
  
  # Write the membership matrix
  
  resultsmembership <- data.frame(all.membership)
  row.names(resultsmembership) <- col.names
  filename <- paste(directory, doc.string, ".", "membership", ".gct", sep="", collapse="")
  MP.WriteGct(resultsmembership , filename)
  
  y.range <- c(1 - 2*(1 - min(rho)), 1)
  plot(k.vector, rho, main ="Cophenetic Coefficient", xlim=c(k.init, k.final), ylim=y.range, xlab = "k", ylab="Cophenetic correlation", type = "n")
  lines(k.vector, rho, type = "l", col = "black")
  points(k.vector, rho, pch=22, type = "p", cex = 1.25, bg = "black", col = "black")
  
  if (non.interactive.run == F) {  
    if (.Platform$OS.type == "windows") {
      savePlot(filename = filename, type ="jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }
  
  xx <- cbind(k.vector, rho)
  write(xx, file= paste(directory, doc.string, ".", "cophenetic.txt", sep=""))
}


read.dataset <- function(file) {
  result <- regexpr(paste(".gct","$",sep=""), tolower(file))
  if(result[[1]] != -1)
    return(read.gct(file))
  result <- regexpr(paste(".res","$",sep=""), tolower(file))
  if(result[[1]] != -1)
    return(read.res(file))
  stop("Input is not a res or gct file.")	
}

matrix.abs.plot <- function(V, axes = F, log = F, norm = T, transpose = T, matrix.order = T, max.v = 1, min.v = 0, main = " ", sub = " ", xlab = " ", ylab = "  ") {
  rows <- length(V[,1])
  cols <- length(V[1,])
  if (log == T) {
    V <- log(V)
  }
  B <- matrix(0, nrow=rows, ncol=cols)
  for (i in 1:rows) {
    for (j in 1:cols) {
      if (matrix.order == T) {
        k <- rows - i + 1
      } else {
        k <- i
      }
      if (norm == T) {
        if ((max.v == 1) && (min.v == 0)) {
          max.val <- max(V)
          min.val <- min(V)
        } else {
          max.val = max.v
          min.val = min.v
        }
      }
      B[k, j] <-  max.val - V[i, j] + min.val
    }
  }
  if (transpose == T) {
    B <- t(B)
  }
  if (norm == T) {
    image(z = B, zlim = c(min.val, max.val), axes = axes, col = rainbow(100, s = 1.0, v = 0.75, start = 0.0, end = 0.75, gamma = 1.5), main = main, sub = sub, xlab = xlab, ylab = ylab) 
  } else {
    image(z = B, axes = axes, col = rainbow(100, s = 1, v = 0.6, start = 0.1, end = 0.9, gamma = 1), main = main, sub = sub, xlab = xlab, ylab = ylab) 
  }
  return(list(B, max.val, min.val))
}

metagene.plot <- function(H, main = " ", sub = " ", xlab = "samples ", ylab = "amplitude") {
  k <- length(H[,1])
  S <- length(H[1,])
  index <- 1:S
  maxval <- max(H)
  minval <- min(H)
  plot(index, H[1,], xlim=c(1, S), ylim=c(minval, maxval), main = main, sub = sub, ylab = ylab, xlab = xlab, type="n")
  for (i in 1:k) {
    lines(index, H[i,], type="l", col = i, lwd=2)
  }
}

ConsPlot <- function(V, col.labels, col.names, main = " ", sub = " ", xlab=" ", ylab=" ") {
  # Plots a heatmap plot of a consensus matrix
  
  cols <- length(V[1,])
  B <- matrix(0, nrow=cols, ncol=cols)
  max.val <- max(V)
  min.val <- min(V)
  for (i in 1:cols) {
    for (j in 1:cols) {
      k <- cols - i + 1
      B[k, j] <-  max.val - V[i, j] + min.val
    }
  }
  
  col.names2 <- rev(col.names)
  col.labels2 <- rev(col.labels)
  D <- matrix(0, nrow=(cols + 1), ncol=(cols + 1))
  
  col.tag <- vector(length=cols, mode="numeric")
  current.tag <- 0
  col.tag[1] <- current.tag
  for (i in 2:cols) {
    if (col.labels[i] != col.labels[i - 1]) {
      current.tag <- 1 - current.tag
    }
    col.tag[i] <- current.tag
  }
  col.tag2 <- rev(col.tag)
  D[(cols + 1), 2:(cols + 1)] <- ifelse(col.tag %% 2 == 0, 1.02, 1.01)
  D[1:cols, 1] <- ifelse(col.tag2 %% 2 == 0, 1.02, 1.01)
  D[(cols + 1), 1] <- 1.03
  D[1:cols, 2:(cols + 1)] <- B[1:cols, 1:cols]
  
  col.map <- c(rainbow(100, s = 1.0, v = 0.75, start = 0.0, end = 0.75, gamma = 1.5), "#BBBBBB", "#333333", "#FFFFFF")
  image(1:(cols + 1), 1:(cols + 1), t(D), col = col.map, axes=FALSE, main=main, sub=sub, xlab= xlab, ylab=ylab)
  for (i in 1:cols) {
    col.names[i]  <- paste("      ", substr(col.names[i], 1, 12), sep="")
    col.names2[i] <- paste(substr(col.names2[i], 1, 12), "     ", sep="")
  }
  
  axis(2, at=1:cols, labels=col.names2, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.50, font.axis=1, line=-1)
  axis(2, at=1:cols, labels=col.labels2, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.65, font.axis=1, line=-1)
  
  axis(3, at=2:(cols + 1), labels=col.names, adj= 1, tick=FALSE, las = 3, cex.axis=0.50, font.axis=1, line=-1)
  axis(3, at=2:(cols + 1), labels=as.character(col.labels), adj = 1, tick=FALSE, las = 1, cex.axis=0.65, font.axis=1, line=-1)
  
  return()
}

read.res <- function(filename = "NULL") { 
  #
  # Reads a gene expression dataset in RES format and converts it into an R data frame
  #
  header.cont <- readLines(filename, n = 1)
  temp <- unlist(strsplit(header.cont, "\t"))
  colst <- length(temp)
  header.labels <- temp[seq(3, colst, 2)]
  ds <- read.delim(filename, header=F, row.names = 2, sep="\t", skip=3, blank.lines.skip=T, comment.char="", as.is=T)
  colst <- length(ds[1,])
  cols <- (colst - 1)/2
  rows <- length(ds[,1])
  A <- matrix(nrow=rows - 1, ncol=cols)
  A <- ds[1:rows, seq(2, colst, 2)]
  table1 <- data.frame(A)
  names(table1) <- header.labels
  return(table1)
}

read.gct <- function(filename = "NULL") { 
  #
  # Reads a gene expression dataset in GCT format and converts it into an R data frame
  #
  ds <- read.delim(filename, header=T, sep="\t", skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T)
  ds <- ds[-1]
  return(ds)
}