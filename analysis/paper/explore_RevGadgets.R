# explore RevGadgets

library(devtools)
devtools::install_github("cmt2/RevGadgets")
library(RevGadgets)

library(ggplot2)
library(ggtree)
library(grid)
library(gridExtra)


# download the example dataset to working directory
url <-
  "https://revbayes.github.io/tutorials/intro/data/primates_cytb_GTR.log"

dest_path <- "primates_cytb_GTR.log"

download.file(url, dest_path)


# to run on your own data, change this to the path to your data file
file <- dest_path

# read the trace and discard burnin
trace_quant <- readTrace(path = file, burnin = 0.1)

# or read the trace _then_ discard burnin
trace_quant <- readTrace(path = file, burnin = 0)
trace_quant <- removeBurnin(trace = trace_quant, burnin = 0.1)

# assess convergence with coda
library(coda)
trace_quant_MCMC <- as.mcmc(trace_quant[[1]])
effectiveSize(trace_quant_MCMC)
traceplot(trace_quant_MCMC)

summarizeTrace(trace = trace_quant, vars =  c("pi[1]","pi[2]","pi[3]","pi[4]"))

plotTrace(trace = trace_quant, vars = c("pi[1]","pi[2]","pi[3]","pi[4]"))[[1]]

##########example 1###########
# download the example dataset to working directory
url_rj <- "https://revbayes.github.io/tutorials/intro/data/freeK_RJ.log"
dest_path_rj <- "freeK_RJ.log"
download.file(url_rj, dest_path_rj)
file <- dest_path_rj
# read in trace
trace <- readTrace(path = file)
plots <- plotTrace(trace = trace,
                   vars = c("prob_rate_12", "prob_rate_13",
                            "prob_rate_31", "prob_rate_32"))
plots[[1]]
trace_qual <- readTrace(path = file)

# summarize parameters
summarizeTrace(trace_qual,
               vars = c("prob_rate_12", "prob_rate_13", "prob_rate_21",
                        "prob_rate_23", "prob_rate_31", "prob_rate_32"))

plots <- plotTrace(trace = trace_qual,
                   vars = c("prob_rate_12", "prob_rate_13",
                            "prob_rate_31", "prob_rate_32",
                            "rate_31", "rate_32"))

grid.newpage()
grid.draw( # draw the following matrix of plots
  rbind( # bind together the columns
    ggplotGrob(plots[[1]]),
    ggplotGrob(plots[[2]]))
)

############example 2 and Basic tree plots###########
# download the example dataset to working directory
url_nex <-
  "https://revbayes.github.io/tutorials/intro/data/primates_cytb_GTR_MAP.tre"
dest_path_nex <- "primates_cytb_GTR_MAP.tre"
download.file(url_nex, dest_path_nex)

# to run on your own data, change this to the path to your data file
file <- dest_path_nex

tree <- readTrees(paths = file)

tree_rooted <- rerootPhylo(tree = tree, outgroup = "Galeopterus_variegatus")

# create the plot of the rooted tree
plot <- plotTree(tree = tree_rooted,
                 # label nodes the with posterior probabilities
                 node_labels = "posterior",
                 # offset the node labels from the nodes
                 node_labels_offset = 0.005,
                 # make tree lines more narrow
                 line_width = 0.5,
                 # italicize tip labels
                 tip_labels_italics = TRUE)

# add scale bar to the tree and plot with ggtree
library(ggtree)
plot + geom_treescale(x = -0.35, y = -1)


############example 3 and Fossilized birth-death trees plots###########
# download the example dataset to working directory

file <- system.file("extdata", "fbd/bears.mcc.tre", package="RevGadgets")
tree <- readTrees(paths = file)
plotFBDTree(tree = tree, timeline = TRUE, tip_labels_italics = FALSE,
            tip_labels_remove_underscore = TRUE,
            node_age_bars = TRUE, age_bars_colored_by = "posterior",
            age_bars_color = rev(colFun(2))) +
  ggplot2::theme(legend.position=c(.25, .85))



# plot the FBD tree
plotFBDTree(tree = tree,
            timeline = T,
            geo_units = "epochs",
            tip_labels_italics = T,
            tip_labels_remove_underscore = T,
            tip_labels_size = 3,
            tip_age_bars = T,
            node_age_bars = T,
            age_bars_colored_by = "posterior",
            label_sampled_ancs = TRUE) +
  # use ggplot2 to move the legend and make
  # the legend background transparent
  theme(legend.position=c(.05, .6),
        legend.background = element_rect(fill="transparent"))



############example 4 and Coloring branches by variables###########
file <- system.file("extdata",
                    "relaxed_ou/relaxed_OU_MAP.tre",
                    package="RevGadgets")


# read in the tree
tree <- readTrees(paths = file)

# plot the tree with rates
plotTree(tree = tree,
         # italicize tip labels
         tip_labels_italics = FALSE,
         # specify variable to color branches
         color_branch_by = "branch_thetas",
         # thicken the tree lines
         line_width = 1.7) +
  # move the legend with ggplot2
  theme(legend.position=c(.1, .9))


############Standard (anagenetic) models###########

file <- system.file("extdata",
                    "comp_method_disc/ase_freeK.tree",
                    package="RevGadgets")

# process the ancestral states
freeK <- processAncStates(file,
                          # Specify state labels.
                          # These numbers correspond to
                          # your input data file.
                          state_labels = c("1" = "Epitheliochorial",
                                           "2" = "Endotheliochorial",
                                           "3" = "Hemochorial"))

# produce the plot object, showing MAP states at nodes.
# color corresponds to state, size to the state's posterior probability
plotAncStatesMAP(t = freeK,
                 tree_layout = "circular") +
  # modify legend location using ggplot2
  theme(legend.position = c(0.57,0.41))


###################Cladogenetic models##################

file <- system.file("extdata", "dec/simple.ase.tre", package="RevGadgets")

# Create the labels vector.
# This is a named vector where names correspond
# to the computer-readable numbers generated
# in the biogeographic analysis and the values
# are character strings of whatever you'd like
# as labels on the figure. The state.labels.txt
# file produced in the analysis links the
# computer-readable numbers with presence/ absence
# data for individual ranges.
labs <- c("1"  = "K",   "2"  = "O",
          "3"  = "M",   "4"  = "H",
          "5"  = "KO",  "6"  = "KM",
          "7"  = "OM",  "8"  = "KH",
          "9"  = "OH",  "10" = "MH",
          "11" = "KOM", "12" = "KOH",
          "13" = "KMH", "14" = "OMH",
          "15" = "KOMH")

# pass the labels vector and file name to the processing script
dec_example <- processAncStates(file, state_labels = labs)

# You can see the states sampled in the analysis in the
# dec_example@state_labels vector. This may be different
# from the `labs` vector you provided above if not all
# possible states are included in the annotated tree.
dec_example@state_labels

# We are going to generate colors for these states using
# a color palette, but you could also specify a color per
# state manually.

# Get the length of the dec_example$state_labels vector
# to know how many colors you need.
ncol <- length(dec_example@state_labels)

# We use colorRampPalette() to generate a function that will
# expand the RevGadgets color palette (colFun) to the necessary
# number of colors, but you can use any colors you like as long
# as each state_label has a color.
colors <- colorRampPalette(colFun(12))(ncol)

# Name the color vector with your state labels and then order
# it in the order you'd like the ranges to appear in your legend.
# Otherwise, they will appear alphabetically.
names(colors) <- dec_example@state_labels
colors <- colors[c(1,2,9,11,
                   3,4,6,10,12,13,
                   5,7,14,
                   8)]

# Plot the results with pies at nodes
pie <- plotAncStatesPie(t = dec_example,
                        # Include cladogenetic events
                        cladogenetic = TRUE,
                        # Add text labels to the tip pie symbols
                        tip_labels_states = TRUE,
                        # Offset those text labels slightly
                        tip_labels_states_offset = .05,
                        # Pass in your named and ordered color vector
                        pie_colors = colors,
                        # Offset the tip labels to make room for tip pies
                        tip_labels_offset = .2,
                        # Move tip pies right slightly
                        tip_pie_nudge_x = .07,
                        # Change the size of node and tip pies
                        tip_pie_size = 0.8,
                        node_pie_size = 1.5) +
  # Move the legend
  theme(legend.position = c(0.1, 0.75))

map <- plotAncStatesMAP(t = dec_example,
                        # Include cladogenetic events
                        cladogenetic = T,
                        # Pass in the same color vector
                        node_color = colors,
                        # adjust tip labels
                        tip_labels_offset = 0.1,
                        # increase tip states symbol size
                        tip_states_size = 3) +
  # adjust legend position and remove color guide
  theme(legend.position = c(0.2, 0.87)) +
  guides(color = FALSE)


grid.arrange(pie,map, ncol = 2)


############## State-Dependendent Diversification Analysis #########
url <-
  "https://revbayes.github.io/tutorials/intro/data/primates_BiSSE_activity_period.log"

bisse_file <- "primates_BiSSE_activity_period.log"
download.file(url, bisse_file)


pdata <- processSSE(bisse_file)

# plot the rates
plotMuSSE(pdata)

##############file doesn't exist #########3
# bisse_anc_states_file <- "data/anc_states_BiSSE.tree"
# p_anc <- processAncStates(path = bisse_anc_states_file)
#
# # plot the ancestral states
# plotAncStatesMAP(p_anc, tree_layout = "circular")

#################
url_rates <-
  "https://revbayes.github.io/tutorials/intro/data/primates_BDS_rates.log"



branch_specific_file <- "primates_BDS_rates.log"
branch_specific_tree_file <- "primates_tree.nex"
download.file(url_rates, branch_specific_file)
rates <- readTrace(branch_specific_file)
tree  <- readTrees(branch_specific_tree_file)

combined <- processBranchData(tree    = tree,
                              dat     = rates,
                              net_div = TRUE)

plotTree(combined, color_branch_by = "net_div",
         tip_labels_size = 2, tree_layout = "circular")

#########Episodic Diversification Analysis












