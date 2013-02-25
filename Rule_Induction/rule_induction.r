#Run Command: R --slave --args -i otu_table.txt -c "confidence" -s "support" -t "sort" --source_dir ~/qiime_software/qiime-1.5.0-release/qiime/support_files/R < rule_induction.r


args <- commandArgs(trailingOnly=TRUE)
if(!is.element('--source_dir', args)){
    stop("\n\nPlease use '--source_dir' to specify the R source code directory.\n\n")
}
sourcedir <- args[which(args == '--source_dir') + 1]
source(sprintf('%s/loaddata.r',sourcedir))
source(sprintf('%s/util.r',sourcedir))

load.library('optparse')
load.library('arulesViz')


option_list = list(
    make_option(c("--source_dir"), type="character",
        help="Path to R source directory [required]."),
    make_option(c("-i", "--otutable"), type="character",
        help="Input otu table [required]."),
    make_option(c("-c", "--confidence"), type="double",
        help="Input confidence threshold [required]."),
    make_option(c("-s", "--support"), type="double",
        help="Input support threshold [required]."),
    make_option(c("-o", "--sort"), type="character",
        help="Input sort threshold [required].")
)


opts <- parse_args(OptionParser(option_list=option_list), args=args)

adata <- load.qiime.otu.table(opts$otutable)

save(adata,file="otu.rda");
mydata <- load("otu.rda");
newdata <- get(mydata);


finaldata <- as(newdata, "transactions");

conf <- as.double(opts$confidence)
sup <- as.double(opts$support)
sort_option <- c(opts$sort)

sink('logfile.txt')

## find all closed frequent itemsets
closed <- apriori(finaldata, parameter = list(target = "closed", confidence = conf,  support = sup));
## rule induction
rules <- ruleInduction(closed, finaldata, confidence = conf, control = NULL);
rules_apriori <- ruleInduction(closed, finaldata, confidence = conf, control = list(method = "apriori", verbose = TRUE));
sorted_rules <- sort(rules, by = sort_option)
sorted_rules_sub <- subset(rules, subset = rhs %pin% "HIV")
sorted_rules_apriori <- sort(rules_apriori, by = sort_option)
sorted_rules_apriori_sub <- subset(rules_apriori, subset = rhs %pin% "HIV")
plot(sorted_rules, method="grouped")
# sort rules by user input and print out two text files
sink('eclat_rules.txt')
inspect(sorted_rules[1:min(length(sorted_rules), 500)])
sink('apriori_rules.txt')
inspect(sorted_rules_apriori[1:min(length(sorted_rules_apriori), 500)])
