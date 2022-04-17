# drop tree tip

suppressMessages(library("argparse"))
suppressMessages(library("treeio"))
suppressMessages(library("dplyr"))

parser <- ArgumentParser()
parser$add_argument("--tree", required=TRUE, help="Input original tree file path")
parser$add_argument("--lineage", required=TRUE, help="Input lineage file path")
parser$add_argument("--output", required=TRUE, help="Output directory path")
args <- parser$parse_args()

# tree <- read.tree(args$tree)
tree <- read.tree(paste(args$tree, 'timetree.nwk', sep='/'))
sublineages <- read.csv(args$lineage)
tree_info <- as_tibble(tree)
target = sublineages[,1]
not_target <- filter(tree_info, !(label %in% target))
to_drop <- as.matrix(filter(not_target, !substring(label,1,4) == "NODE")['label'])
subtree <- drop.tip(tree, as.matrix(to_drop))
write.tree(subtree, file=paste(args$output, 'subtree.nwk', sep=''))
