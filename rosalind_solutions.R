
# A Rapid Introduction to Molecular Biology -------------------------------
s1 <- "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"

composistionCount <- function(seq) {
  seq <- unlist(strsplit(seq, split = ""))
  out <- table(seq)
  return(out)
}

s2 <- 
  readLines("1_rosalind_dna.txt")
composistionCount(seq = s2)



# The Second Nucleic Acid -------------------------------------------------

s1 <- "GATGGAACTTGACTACGTAAATT"

transcribeDNA <- function(seq){
  seq <- unlist(strsplit(seq, split = ""))
  pos <- which(unlist(strsplit(seq, split = "")) == "T")
  seq[pos] <- "U"
  seq <- paste0(seq, collapse = "")
  return(seq)
}

s2 <- 
  readLines("2_rosalind_rna.txt")

transcribeDNA(seq = s2)


# The Secondary and Tertiary Structures of DNA ----------------------------
s1 <- "AAAACCCGGT"

reverse_complement <- function(seq, complement_only = FALSE) {
  
  seq <- rev(unlist(strsplit(seq, split = "")))
  
  a <- which(seq == "A")
  t <- which(seq == "T")
  c <- which(seq == "C")
  g <- which(seq == "G")
  
  seq[a] <- "T"
  seq[t] <- "A"
  seq[c] <- "G"
  seq[g] <- "C"
  
  out <- paste0(seq, collapse = "")
  
  if(complement_only) {
    out <- paste0(rev(unlist(strsplit(out, split = ""))), collapse = "")
  } else {
    out <- out
  }
  
  return(out)
}

s2 <- 
  readLines("3_rosalind_revc.txt")

reverse_complement(seq = s2)

# hammingDistance ---------------------------------------------------------

s1 <- "GAGCCTACTAACGGGAT"
s2 <- "CATCGTAATGACGGCCT"

hammingDistance <- function(seq1, 
                            seq2) {
  seq1 <- unlist(strsplit(seq1, split = ""))
  seq2 <- unlist(strsplit(seq2, split = ""))
  
  out <- sum(! seq1 == seq2)
  return(out)
}

s2 <- 
  readLines("4_rosalind_hamm.txt")

hammingDistance(seq1 = s2[[1]],
                seq2 = s2[[2]])



# rabbits -----------------------------------------------------------------
# needed to get help from forum
rF <- c()
n = 28
k = 2
i = 0
for(i in 1:n) {
  if(i == 1) {
    rF[i] <- 1
  } else if(i == 2) {
    rF[i] <- 1
  } else rF[i] = rF[i-1] + k*rF[i-2]
}
format(rF[length(rF)], scientific = F)

# Identifying Unknown DNA Quickly -----------------------------------------

s2 <-
  readLines("5_rosalind_gc.txt")

s1 <-
  readLines("5_test.txt")



gc_content <- function(seq) {
  
  fasta_names <- grep(x = seq, pattern = "^>.+[0-9]$", value = T)
  l <- lapply(strsplit(paste0(seq, collapse = ""), 
                       split = ">Rosalind_", 
                       fixed = F, 
                       perl = T), 
              substring, 5)
  l <- l[[1]][-1]
  l <- lapply(X = l, FUN = function(x) {table(strsplit(x, split = ""))})
  
  names(l) <- fasta_names
  out <- 
    lapply(X = l, 
           FUN = function(x) {
             round(sum(x[c(2,3)])/sum(x) * 100, digits = 4)
           })
  out <- sort(unlist(out), decreasing = T)[1]
  return(out)
}

gc_content(seq = s2)

# The Genetic Code --------------------------------------------------------

amino_acids <- read.delim2("~/practice/rosalind/aminoacids.txt", header = F)

get_aa_dict <- function(path) {
  amino_acids <- read.delim2(path, header = F)
  aa_split <- gsub(x = strsplit(amino_acids$V1, split = "\t"),
                   replacement = "",
                   pattern = " ")
  
  all <- unlist(strsplit(aa_split, split = ""))
  
  is_lower <- which(grepl("[[:lower:]]", all))
  lower_start <- is_lower[seq(1, length(is_lower), 3)]-1
  
  all[lower_start] <- "STOP"
  all <- all[which(grepl("[[:upper:]]", all))]
  
  amino_names <- all[seq(4, length(all), 4)]
  start <- seq(1, length(all), 4)
  stop <- seq(3, length(all), 4)
  
  seq <- c()
  for(i in 1:length(amino_names)) {
    seq[i] <- paste0(all[start[i]:stop[i]], collapse = "")
  }
  
  names(seq) <- amino_names
  
  return(seq)
}



seq <- get_aa_dict(path = "~/practice/rosalind/aminoacids.txt")
seq

get_amino <- function(sequence) {
  test_seq <- unlist(strsplit(sequence, split = ""))
  start <- seq(1, length(test_seq), 3)
  stop <- seq(3, length(test_seq), 3)
  
  out <- c()
  for(i in 1:(length(test_seq)/3)) {
    out[i] <- paste0(test_seq[start[i]:stop[i]], collapse = "")
  }
  
  aa <- c()
  for(i in 1:length(out)) {
    aa[i] <- names(seq[seq %in% out[i]])
  }
  
  final_out <- paste0(aa[1:length(out) -1], collapse = "")
  return(final_out)
}

s <- "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
get_amino(sequence = s)

task_string <-  paste(readLines("6_rosalind_prot.txt"), 
                      collapse="\n")
get_amino(sequence = task_string)


# Finding a Motif in DNA --------------------------------------------------

# help from https://stat.ethz.ch/pipermail/r-help/2009-December/415501.html
# since could not catch overlaps
s <- "GATATATGCATATACTT"
m <- "ATAT"

findall <- function(motif, sequence) {
  stopifnot(length(motif) == 1, length(sequence) == 1)
  pos <- c()  # positions of matches
  i <- 1; n <- nchar(sequence)
  found <- regexpr(motif, substr(sequence, i, n), perl=TRUE)
  while (found > 0) {
    pos <- c(pos, i + found - 1)
    i <- i + found
    found <- regexpr(motif, substr(sequence, i, n), perl=TRUE)
  }
  return(pos)
}

task_string <-  paste(readLines("7_rosalind_subs.txt"), 
                      collapse="\n")
s <- strsplit(task_string, split = "\n")

s1<- s[[1]][1]
m1 <- s[[1]][2]

findall(motif = m1, sequence = s1)


# Consensus and Profile ---------------------------------------------------

new_consensus <- function(file) {
  t <- readLines(file)
  headers <- grep(">", t)
  s <- data.frame(headers = t[headers], 
                  from = headers+1, 
                  to = c((headers-1)[-1], length(t))
  )
  # Process sequence lines
  seqs <- rep(NA, length(headers))
  for(i in 1:length(headers)) {
    seqs[i]<-paste(t[s$from[i]:s$to[i]], collapse="")
  }
  
  DF <- data.frame(name=gsub(">", "", t[headers]), sequence=seqs)
  
  string_list <- as.list(DF$sequence)
  #names(string_list) <- DF$name
  l <- lapply(string_list, strsplit, split = "")
  m <- matrix(unlist(l), ncol = length(l[[1]][[1]]), byrow = TRUE)
  
  mtx_out <- matrix(nrow = 4, 
                    ncol = ncol(m), 
                    dimnames = 
                      list(c("A", "C", "G", "T"),
                           NULL)
  )
  
  for(i in seq_along(rownames(mtx_out))){
    mtx_out[i, ] <-  
      apply(X = m,
            MARGIN = 2,
            FUN = function(x) {
              sum(x == rownames(mtx_out)[i])
            }
      )
  }
  
  cons <- 
    apply(X = mtx_out,
          MARGIN = 2,
          FUN = function(x) {
            #sum(x == rownames(mtx_out)[i])
            names(which(x == max(x)))
          }
    )
  
  out <- paste0(unlist(lapply(X = cons, FUN = function(x) {x[1] })), collapse = "")
  
  res <- list(mtx_out, out)
  return(res)
}

#t <- strsplit(t, split = t[grepl(pattern = '>Rosalind_([0-9]{4})+$', x = t)])

#solution
t <- new_consensus(file = "9_rosalind_cons.txt")
t[[2]]
paste0(t[[1]][1,], collapse = " ")
paste0(t[[1]][2,], collapse = " ")
paste0(t[[1]][3,], collapse = " ")
paste0(t[[1]][4,], collapse = " ")

# Introduction to Mendelian Inheritance -----------------------------------

mendel_first_law <- function(AA, Aa, aa){
  s = AA + Aa + aa
  P_aa = (aa/s)*((aa-1)/(s-1)+Aa/((s-1)*2))
  P_Aa = 0.5*(Aa/s)*(aa/(s-1)+(Aa-1)/((s-1)*2))
  out <- 1-(P_aa+P_Aa)
  return(out)
}

mendel_first_law(29, 27, 15)

# variables and some arithmetic --------------------------------------------
hypotenuse_square <- function(a,b) {
  a^2 + b^2
}

hypotenuse_square(929, 823)


# Conditions and Loops ----------------------------------------------------

sum_odd_numbers <- function(a,b) {
  sum(seq(a, b, 1)[seq(a, b, 1) %% 2 == 1])
}

a = 4877 
b = 9828

sum_odd_numbers(a, b)


# strings and lists -------------------------------------------------------

#examlple
s <- "HumptyDumptysatonawallHumptyDumptyhadagreatfallAlltheKingshorsesandalltheKingsmenCouldntputHumptyDumptyinhisplaceagain."
a = 22 
b = 27 
c = 97 
d = 102
ranges <- c("a" = a, "b" = b, "c" = c, "d" = d)

extract_string_in_R <- function(s, a, b, c, d) {
  ranges <- c("a" = a, "b" = b, "c" = c, "d" = d)
  ranges <- ranges + 1
  string <- unlist(strsplit(s, split = ""))
  
  out <- 
    paste(
      paste0(string[ranges["a"]:ranges["b"]], collapse = ""),
      paste0(string[ranges["c"]:ranges["d"]], collapse = ""),
      collapse = " ")
  return(out)
  
}

extract_string_in_R(s = s, 22, 27, 97, 102)

#solution
s <- paste(readLines("rosalind_ini3.txt"), 
           collapse="\n")

s <- strsplit(s, split = "\n")

extract_string_in_R(s = s[[1]][1], 3, 9, 160, 162)


# reading and writing -----------------------------------------------------
##example
s <- readLines("example_rosalind_reading_and_writing.txt")

## solution
s <- readLines("rosalind_ini5.txt")
cat(paste(s[seq(2, length(s), 2)], collapse = "\n"))           



# dictionaries -------------------------------------------------------------

s <- "We tried list and we tried dicts also we tried Zen"
# s <- readLines("rosalind_ini6.txt")
s <- strsplit(s, split = " ")

out <- 
  as.data.frame(table(s[[1]]))

out



# graph -------------------------------------------------------------------

read_fasta <- function(file) {
  t <- readLines(file, warn = FALSE)
  headers <- grep(">", t)
  s <- data.frame(headers = t[headers], 
                  from = headers+1, 
                  to = c((headers-1)[-1], length(t))
  )
  # Process sequence lines
  seqs <- rep(NA, length(headers))
  for(i in 1:length(headers)) {
    seqs[i]<-paste(t[s$from[i]:s$to[i]], collapse="")
  }
  
  DF <- data.frame(name = t[headers], 
                   sequence = seqs) 
  
}

# working solution
s <- read_fasta(file = "/home/jason/rosalind/ex_rosalind_graph.txt")
test <- read_fasta(file = "/home/jason/rosalind/rosalind_grph.txt")
s <- read_fasta(file = "/home/jason/rosalind/rosalind_grph_2.txt")

find_graph <- function(data) {
  s <- data
  l <- lapply(as.list(s$sequence), strsplit, split = "")
  string_length <- sapply(X = s$sequence, nchar)
  fasta_names <- gsub(x = unlist(strsplit(s$name, split = " ")),
                      pattern = ">",
                      replacement = "")
  
  t <- list()
  for(row1 in 1:length(l)) {
    for(row2 in 1:length(l)) {
      if(row1 == row2)
        next
      t[[paste(fasta_names[row1], fasta_names[row2])]] <-
        all(tail(l[[row1]][[1]], 3) == head(l[[row2]][[1]], 3))
    }
  }
  
  t <- t[sapply(t, FUN = function(x) { x[1] == TRUE})]
  
  cat(paste(names(t), collapse = "\n"))
}


find_graph(data = s)



# The Need for Averages ---------------------------------------------------


data <- c(17862, 19995, 16460, 19114, 19485, 16802)

chances <- c(4,4,4,3,2,0)

sum(data * (chances/2))

# from formula:
(min(data*chances) + sum(data * chances))/2


# Finding a Shared Motif --------------------------------------------------

e <- c("A" = "GATTACA", 
       "B" = "TAGACCA", 
       "C" = "ATACA")



# e <- sapply(e, strsplit, split = "")
# sapply(e, length)
# vector_max <- max(sapply(e, length))
# 
# z <- list()
# for (i in 1:length(e)) {
#   if (length(e[[i]]) < vector_max) {
#     z[[i]] <- append(e[[i]], rep("0", vector_max - length(e[[i]])))
#     #append(e[[x]], rep("0", max(sapply(e, length)) - length(e[[x]])))
#   } else{
#     z[[i]] <- e[[i]]
#   }
# }
# 
# t <- list()
# for(row1 in 1:length(z)) {
#   for(row2 in 1:length(z)) {
#     if(row1 == row2)
#       next
#     t[[row1]] <-
#       z[[row1]] == z[[row2]]
#   }
# }


# PTXQC
PTXQC::LCSn(e, min_LCS_length = 0)

#solution
# slow but works

e <- read_fasta("rosalind/rosalind_lcsm.txt")
vec <- e$sequence
names(vec) <- e$name
PTXQC::LCSn(vec, min_LCS_length = 2)

# faster
x <- PTXQC::LCSn(vec[1:2], min_LCS_length = 2)
PTXQC::LCSn(c(x, vec[3]), min_LCS_length = 2)



# Finding a Protein Motif -------------------------------------------------

ids1 <- c("A2Z669",
         "B5ZC00",
         "P07204_TRBM_HUMAN",
         "P20840_SAG1_YEAST")

ids1 <- c("Q4FZD7",
          "P06870_KLK1_HUMAN",
          "P02974_FMM1_NEIGO",
          "Q0AYI5",
          "P55067_PGCN_RAT",
          "P01588_EPO_HUMAN",
          "P07204_TRBM_HUMAN",
          "P01215_GLHA_HUMAN",
          "A1X8C2",
          "P10761_ZP3_MOUSE",
          "Q2GBA3",
          "Q00001_RHGA_ASPAC")

getMotifPosition <- function(ids, 
                             SearchPattern = "N(?=[^P](S|T)[^P])") {
  
  ids <- sub(x = ids, pattern = "_.*", replacement = "")
  protein_seq <-
    sapply(X = ids, 
           FUN = function(x) {
             protr::getUniProt(id = x)
             }
           )
  
  out <- 
    sapply(X = protein_seq, 
           FUN = function(x) {
             unlist(gregexpr(text = x, pattern = SearchPattern, perl = TRUE))
             }
           )
  
  names(out) <-ids1
  out <- out[sapply(out, FUN = function(x) { x[1] > 0 } )]
  
  for(i in seq_along(out)){
    cat(names(out)[i], "\n")
    cat(out[[i]], "\n")
    }
} 

getMotifPosition(ids = ids1)



# Inferring mRNA from Protein ---------------------------------------------

# inspired by https://github.com/Rgtemze/rosalind-solutions/blob/master/bioinformatics-stronghold/MRNA-Inferring-mRNA-from-Protein.py
# at the begining I did not codnisder 3 stop codons

mRNA_from_Protein <- function(amino_acid) {
  protein_to_possible_rna_dict <- 
    sort(
      c("F" = 2, "L" = 6, "S" = 6, "Y" = 2, "C" = 2, 
        "W" = 1, "P" = 4, "H" = 2, "Q" = 2, "R" = 6, 
        "I" = 3, "M" = 1, "T" = 4, "N" = 2, "K" = 2, 
        "V" = 4, "A" = 4, "D" = 2, "E" = 2, "G" = 4)
    )
  
  stop_codons = c("UAA","UGA","UGA")
  number_of_stop_codons = 3
  
  # get sequence
  aa_sequence <- sapply(amino_acid, strsplit, split = "")
  
  #Assuming input is provided properly
  possibility_amount <- 1
  for (aa in seq_along(aa_sequence[[1]])) {
    possibility_amount <- c(possibility_amount * protein_to_possible_rna_dict[aa_sequence[[1]][aa]]) %% 1E6
  }
  
  #Number of stop codons is taken into account. 
  #Despite being not explicity defined in the sequence they exist
  possibility_amount <- (possibility_amount * number_of_stop_codons)
  
  print(possibility_amount[[1]])
}


mRNA_from_Protein(amino_acid = "MA")

test_aa <- "MFAESLCDIIYNGPFPECHLRAMSREKALNMEKMVEKYQQWGYWHFMNISMDVYHPAASMHFWFYMFNMMKNSLMCHNVTGYQRCHQTCKRVCMIWHGAIPLCSPGEGAQKWESTANRHSNCNWADGDITSGWWCARISELETFISFSEKMHTGRWSYVHCVDWSGLGLTEFQWDNGQAILLRMVPYGRVTLHGMAMDPMTRGDKNEEFICIHREVPSAAAMLERLEHHMAPYCLRSHGDLSYKTHCIWFWPEAHFPNGAWRIQGLPQGGELHMCRIYCGDKWLGCFILYPHPWMNWDMHRRARPMYKRQFFMQPFGHWPTMVEVIKEQSIHHYTFRYCHPCDMCHWREPIGKYLMFMWETACLPTREDHASIASFGTFDRSTTMIHSTQMYADPGKFDGMTTCYFSMIQWARYSMLEPGLHPIPVKQACSVWSPVGSWDTVYPLLHDWQYPLEWKMGSCFHNVQLRCWIHNTSERMVERYNAFPSQYCYHGSQWKHHQACNWYLMHLTKNNDCWCPMWDDIFVFLYQNCAPPLVDNRPMIYPHQDLLFYHIHQSTQFECLCTERCMKSNTLEITHEIWWQCEACGVECCRWWLERHLLECNYTNLFFLRIVDALWMTEAGCMTRDNYPPNLWEYPWKERISPSHKCECQGALWGYCEPLMYKWWIFWISMITVTTSLKPLFKHYMHEMHTLFFDNRDACTLVQWPVWPWWCEVCKTMSAIRVWLSHKNFMFLWHQNFHNALARMVFRSHKSYVFIRVATWWQTSHMKWECQASSPGSVQRAPGKAQSFVTSKKASQVYTCNNLIYFWGTWCKTKLKLNKLTDCTKVQMNADYKIRLKSHLYWIANPDVNDVYTSTYSGNTVARHVNVQCKYTGFTGNEVFHRCEPTFHKNNYIDTMAGFVIQTLVAGGMHFKTICGACHRQVRAWWVYTCCPQKEMGQMKKYDWDIQMKIMQLITIYGIENNPKQFSQTWRESRKWGTHGHNFKLEQFMNRKMFYAYW"

mRNA_from_Protein(amino_acid = test_aa)

# Translating RNA into Protein --------------------------------------------

RNA_to_DNA <- function(seq) {
  seq <- unlist(strsplit(seq, split = ""))
  pos <- which(unlist(strsplit(seq, split = "")) == "U")
  seq[pos] <- "T"
  seq <- paste0(seq, collapse = "")
  return(seq)
}

alternative_aa_from_ORF <- function(seq) {
  # searching for codons
  codons <- strsplit(x = seq, split = "(?<=.{3})", perl = TRUE)[[1]]
  
  # codons positions
  start <- which(strsplit(x = seq, split = "(?<=.{3})", perl = TRUE)[[1]] == "ATG")
  stop <- which(strsplit(x = seq, split = "(?<=.{3})", perl = TRUE)[[1]] %in% c("TAA", "TAG", "TGA"))
  
  # cleaning 
  start <- start[start < max(stop)]
  stop <- c(stop[stop >= min(start) & stop <= max(start)], stop[which(stop > max(start))][1])
  
  out <- c()
  for(i in seq_along(start)) {
    my_stop <- stop[(stop > start[i])][1]
        out <- append(out, 
                      get_amino(sequence = gsub(x = paste0(codons[start[i]:my_stop], collapse = ""), 
                                                pattern = "T", 
                                                replacement = "U")
                      )
        )
  }
  
  out <- out[! grepl(out, pattern = "STOP")]
  return(out)
}

# get data
s <- paste(readLines("~/practice/rosalind/rosalind_orf_3.txt"), 
           collapse="")
s <- gsub(x = s, pattern = ">Rosalind_[0-9]+", replacement = "")
rev_s <- reverse_complement(seq = s, complement_only = FALSE)

#reverse
out <- 
  c(alternative_aa_from_ORF(seq = s),
    alternative_aa_from_ORF(seq = substring(s, 2)),
    alternative_aa_from_ORF(seq = substring(s, 3)),
    alternative_aa_from_ORF(seq = rev_s),
    alternative_aa_from_ORF(seq = substring(rev_s,2)),
    alternative_aa_from_ORF(seq = substring(rev_s,3))
    )

cat(paste(unique(out), collapse = "\n"))



# Enumerating Gene Orders -------------------------------------------------

perm <- function(v) {
  n <- length(v)
  if (n == 1) v
  else {
    X <- NULL
    for (i in 1:n) X <- rbind(X, cbind(v[i], perm(v[-i])))
    X
  }
}

all_permuatations <- function(n) {
  print(nrow(perm(v = 1:n)))
  cat(paste(apply(X = perm(v = 1:n), MARGIN = 1, paste, collapse = " "), collapse = "\n"))
}

all_permuatations(n = 5)



# Calculating Protein Mass ------------------------------------------------

weighted_alphabet <- function(file) {
  t <- paste(readLines(file), 
             collapse = "\n")
  
  t <- unlist(strsplit(t, split = "\n"))
  t <- unlist(strsplit(t, split = " "))
  t <- t[nzchar(t)]
  t_dict <- setNames(object = round(as.numeric(t[seq(2, length(t), 2)]), digits = 8),
                     nm = t[seq(1, length(t), 2)])
  return(t_dict)
}

# calculate mass from mass symbol
calculate_mass <- function(string) {
  s <- unlist(strsplit(string, split = ""))
  out <- sum(t_dict[names(t_dict) %in% names(table(s))] * table(s))
  return(out)
}


t_dict <- weighted_alphabet(file = "~/practice/rosalind/mass_table.txt")

res <-calculate_mass(string = "SKADYEK")
sprintf("%.3f", res)

ex <- readLines("~/practice/rosalind/rosalind_prtm (5).txt")
res2 <- calculate_mass(string = ex)
sprintf("%.3f", res2)


# Locating Restriction Sites ----------------------------------------------
# first solution
get_palindroms <- function(seq) {

  s <- seq
  all_strings <- c()
  start <- c()
  stop <- c()
  for(i in 1:nchar(s)) {
    for(j in 1:nchar(s)) {
      str <-  substring(s, i, j)
      if(nchar(str) >= 4 & nchar(str) <=12 & (nchar(str) %% 2) == 0) {
        if(i < j) {
          all_strings <- append(all_strings, substring(s, first = i,j))
          start <- append(start, i)
          stop <- append(stop, j)
        }
        
      }
    }
  }
  
  out_list <- list()
  is.palindrom <- c()
  for(i in seq_along(all_strings)){
    
    fwd <- substring(all_strings[i], 1, nchar(all_strings[i])/2)
    rev_comp <- reverse_complement(seq = fwd, complement_only = FALSE)
    to_compare_with <- substring(all_strings[i], nchar(all_strings[i])/2 +1, nchar(all_strings[i]))
    is.match <- rev_comp == to_compare_with
    
    if(is.match) {
      is.palindrom <- append(is.palindrom, i)
    }
  }
  
  out <- 
    data.frame("start" = start[is.palindrom], 
               "length" = stop[is.palindrom] - start[is.palindrom] +1)
  out <- out[order(out[,1],decreasing=FALSE),]
  return(out)
  
}

s <- "TCAATGCATGCGGGTCTATATGCAT"
my_out <- get_palindroms(seq = s)
cat(paste(paste(my_out[, 1], my_out[,2], sep = " "), collapse = "\n"))

my_out <- get_palindroms(seq = "TATATA")
cat(paste(paste(my_out[, 1], my_out[,2], sep = " "), collapse = "\n"))


m <- paste(readLines("~/practice/rosalind/rosalind_revp_final.txt"), collapse = "")
m <- gsub(x = m, pattern = ">Rosalind_[0-9]+", replacement = "")
my_out <- get_palindroms(seq = m)
cat(paste(paste(my_out[, 1], my_out[,2], sep = " "), collapse = "\n"))


# second cleaner solution
get_palindroms2 <- function(seq) {
  
  s <- seq
  for(i in 1:nchar(s)) {
    for(j in 1:nchar(s)) {
      str <-  substring(s, i, j)
      if(nchar(str) >= 4 & nchar(str) <=12 & (nchar(str) %% 2) == 0) {
        if(i < j) {
          
          fwd <- substring(str, 1, nchar(str)/2)
          rev_comp <- reverse_complement(seq = fwd, complement_only = FALSE)
          to_compare_with <- substring(str, nchar(str)/2 +1, nchar(str))
          is.match <- rev_comp == to_compare_with
          if(is.match) {
            cat(paste(i , j-i +1, "\n"))
            }
        }
        
      }
    }
  }
}

get_palindroms2(seq = m)



# RNA Splicing ------------------------------------------------------------

# require read_fasta()
# require transcribeDNA()
# require get_amino()
# get_aa_dict()
translateDNA2protein <- function(file) {
  
  data <- read_fasta(file)
  sequence <- data$sequence[1]
  introns <- data$sequence[-1]
  seq <- get_aa_dict(path = "~/practice/rosalind/aminoacids.txt")
  
  for(i in seq_along(introns)) {
    
    sequence <- gsub(x = sequence, pattern = introns[i], replacement = "")
  }
  
  sequence <- transcribeDNA(seq = sequence)
  sequence <- get_amino(sequence = sequence)
  return(sequence)
}


translateDNA2protein(file = "~/practice/rosalind/splice_RNA_example.txt")

translateDNA2protein(file = "~/practice/rosalind/rosalind_splc.txt")



         