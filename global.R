
data("Achilles2.0")
data("Achilles2.4")
data("Achilles2.4_Data")
data("Achilles2.4_annotation")

RNAiObject <- new(Class = "RNAi",
                  genes = Achilles2.4_Data$Gene,
                  sequences = Achilles2.4_Data$Name,
                  cancers = colnames(Achilles2.4_Data)[-(1:2)],
                  cancersP53 = Achilles2.4_annotation,
                  values = Achilles2.4_Data[,-(1:2)])

no_cancers <- c(3,30,10,35,2,7,21,6,4,17,13,20,29,2,1,3,1,10,2)
cancers <- c("PROSTATE", "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "KIDNEY", "CENTRAL_NERVOUS_SYSTEM",
             "SOFT_TISSUE", "SKIN", "LUNG", "BONE",
             "STOMACH", "PANCREAS", "BREAST", "LARGE_INTESTINE",
             "OVARY", "ENDOMETRIUM", "LIVER", "URINARY_TRACT",
             "SMALL_INTESTINE", "OESOPHAGUS", "PLEURA")



multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}