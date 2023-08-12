library("ggplot2");
library("ggpubr");
library("data.table");
library("viridis");

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
# SOC: Start Of Code
#  o Sec-01:
#    - Functions
#  o sec-02:
#    - Variable declerations
#  o sec-03:
#    - Set up for data analysis
#  o sec-04:
#    - Make the time graph
#  o sec-05:
#    - Make the memory graph
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-01:
#  - Functions
#  o sec-01 sub-01:
#    - Uses ggsave to save a graph as an svg
#    - Allows me to change graph saving options easily for
#      every graph this code makes
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*********************************************************\
# Sec-01 Sub-01:
#  - Output: till: Saves output graph as a svg
#  - This is mainly here so I can change graph saving
#    options easily
#*********************************************************/
saveGraph = function(
  nameStr           # name of graph to save
) # Use: saves a graph using ggsave function
{ # saveGraph
    ggsave(paste(nameStr, ".svg", sep = ""), 
           device = "svg", # save as tiff file
           dpi = 300,
    ); # ggsave (save the graph)
} # saveGraph

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-02:
#  - Variable declerations
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

dataDT =
  setDT(                        # convert to a data table
    read.csv(
      "alnSeq-bench-stats.tsv", # File with data
      sep = "\t",               # for tsv files
      header = TRUE             # Has header
)); # Read in my data

graphDT = NULL;       # modified data for graphing
graphObj = NULL;      # Holds the output graph

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
# Sec-03:
#  - Set up for data analysis
#  o sec-03 sub-01:
#    - Rename tests done and organize by ram size
#  o sec-03 sub-02:
#    - Remove the exta O compiled entries for alnSeq and
#      add in hb tag for Hirschberg from bio-align
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>/

#*********************************************************\
# Sec-03 Sub-01:
#  - Rename tests done and organize by ram size
#*********************************************************/

# Sizes for each genome
# query-small is 1761 bases
# reference-small is 1767 bases
# query-mid is 11141 bases
# refrence-mid is 10884 bases
# query-large is 28030 bases
# reference-large is 28030 bases
# query-ramKiller is 199984 bases
# reference-ramKiller is 199982 bases

dataDT$test=paste(dataDT$query,dataDT$reference,sep = "-");
dataDT$test = gsub("ramKiller", "huge", dataDT$test);

dataDT$test = gsub("mid-small", "small-mid", dataDT$test);
dataDT$test=gsub("large-small","small-large",dataDT$test);
dataDT$test = gsub("huge-small","small-huge",dataDT$test);
dataDT$test = gsub("large-mid", "mid-large", dataDT$test);
dataDT$test = gsub("huge-mid", "mid-huge", dataDT$test);
dataDT$test = gsub("huge-large","large-huge",dataDT$test);

dataDT$test =
  factor(
    dataDT$test,
    c(
      "small-small",
      "small-mid",
      "small-large",
      "mid-mid",
      "mid-large",
      "small-huge",
      "large-large",
      "mid-huge",
      "large-huge",
      "huge-huge"
)); # Order data for graphing

#*********************************************************\
# Sec-03 Sub-02:
#  - Remove the exta O compiled entries for alnSeq and add
#    in the hb tag for Hirschberg from bio-align
#*********************************************************/

dataDT =
  dataDT[
    dataDT$Program != "alnSeqOFast" &
    dataDT$Program != "alnSeqO2" &
    dataDT$Program != "alnSeqO0"
    , # Not filtering by columns
]; # Remove the various compiler options

dataDT$Program =
  paste(
    dataDT$Program,
    dataDT$algorithm,
    sep = "-"
); # Comibine program name and aligorithm

# Remove the "-water" and "-needle" program ids
dataDT$Program = gsub("-water", "", dataDT$Program);
dataDT$Program = gsub("-needle", "", dataDT$Program);

# Replacne hb with needle, so is kept same in graphs
dataDT$algorithm = gsub("hb", "needle", dataDT$algorithm);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
# Sec-04:
#  - Make the time graph
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>/

#dev.new();

graphObj = 
  ggplot(
    data =
      dataDT[
        #dataDT$algorithm == "water" &
        dataDT$Program != "seq-align" &
        (
          paste(
            dataDT$Program,
            dataDT$algorithm,
            sep = "-"
          ) != "bio-align-water"
        )
        ,
    ],
    aes(
      x = test,
      y = elapsedTime, # time in seconds
      color = Program,
      shape = Program    
    ),
);  

graphObj =
  graphObj +
  geom_point(
    size = 3,
    alpha = 0.5,
    position = position_jitter(height = 0, width = 0.3)
); # Add points to graph and jitter, to reduce overlap

# separate graphs by neelde/water
graphObj = graphObj + facet_grid(cols = vars(algorithm));

# Transform y-axis by log10, so trends are more clear
graphObj = graphObj + scale_y_sqrt();

graphObj =
  graphObj +
  scale_fill_viridis(
    discrete = FALSE,   # Using continous data
    direction = -1,     # Darkets color used first
    option = "D",       # Default color scheme
    name = "Program"    # name of label
); # Color scheme to use

# Add x and y axis labels
graphObj = graphObj + xlab("Size of aligned genomes");
graphObj = graphObj + ylab("Time used in seconds");

# Add in the theme I like to use
graphObj = graphObj + theme_pubr();

graphObj =
  graphObj +
  theme(axis.text.x = element_text(angle = 90)
); # Rotate x-axis text 90 degrees

saveGraph("alnSeq-time");

#dev.off(); # remove all plots

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
# Sec-05:
#  - Make the memory graph
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>/

#dev.new();

graphObj = 
  ggplot(
    data = dataDT[ dataDT$Program != "seq-align",],
    aes(
      x = test,
      y = maxResidentMemory / 1000, # convert from kb to mb
      color = Program,
      shape = Program    
    ),
);  

graphObj =
  graphObj +
  geom_point(
    size = 3,
    alpha = 0.5,
    position = position_jitter(height = 0, width = 0.3)
); # Add points to graph and jitter, to reduce overlap

# separate graphs by neelde/water
graphObj = graphObj + facet_grid(cols = vars(algorithm));

graphObj =
  graphObj +
  scale_fill_viridis(
    discrete = FALSE,   # Using continous data
    direction = -1,     # Darkets color used first
    option = "D",       # Default color scheme
    name = "Program"    # name of label
); # Color scheme to use

# Transform y-axis to log10
graphObj = graphObj + scale_y_log10();

# Add x and y axis labels
graphObj = graphObj + xlab("Size of aligned genomes");
graphObj = graphObj +ylab("Memory used in megabytes");

# Add in the theme I like to use
graphObj = graphObj + theme_pubr();

graphObj =
  graphObj +
  theme(axis.text.x = element_text(angle = 90)
); # Rotate x-axis text 90 degrees

saveGraph("alnSeq-memory");

#dev.off(); # remove all plots
