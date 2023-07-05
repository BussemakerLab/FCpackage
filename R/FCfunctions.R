#' Write Bash
#'
#' Write bash file with UNIX commands
#' @param script Vector of string that are lines of UNIX commands
#' @param outfile Path and name of bash file to be output
#' @return successful transformation shows a reminder text: "bash script done, remember to dos2unix"
#' @export
write.bash <- function(script, outfile){
  cat("#!/bin/bash\n", file = outfile, sep = "\n", append = F)
  for(i in 1:length(script)){
    cat(script[i], file = outfile, sep = "\n", append = T)
  }
  return('bash script done, remember to dos2unix')
}

#' Write Fasta
#'
#' Write fasta files from data frame containing name and sequence
#' @param data data.frame that has two columns: name and seq, the colnames is crucial
#' @param filename Path and name of fasta file to be output
#' @return Nothing
#' @export
writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
  return()
}


#' Matrix to Sequence
#'
#' Output most preferred sequence of a motif matrix
#' @param scoringMatrix motif matrix with row names labeled as base types
#' @return A string that contains the most preferred sequence of the scoring matrix
#' @export
matrix2seq <- function(scoringMatrix){
  seq <- ''
  for(i in 1:ncol(scoringMatrix)){
    add <- rownames(scoringMatrix)[scoringMatrix[,i] == max(scoringMatrix[,i])]
    if(length(add) > 1){
      add <- 'N'
    }
    seq <- paste(seq, add, sep = '')
  }
  return(seq)
}

#' Pymol: open files
#'
#' Open structure files in pymol showing protiens in sticks and DNA in lines (.pdb/.cif)
#' @param RA.files vector of strings containing the path to the structure files
#' @param pymol.dir Path to pymol executable
#' @param pml Path to safe the pml file
#' @return pymol program with structures loaded
#' @export
pymolOpenFiles <- function(RA.files, pymol.dir, pml = '~/pymolBash.pml'){
  pymol <- c()
  for(i in 1:length(RA.files)){
    pymol[i] <- paste('load ', RA.files[i], sep = '')
  }
  pymol <- c(pymol,'hide all', 'select polymer.protein', 'show stick, sele', 'select polymer.nucleic', 'show line, sele', 'hide lines, hydrogen', 'hide stick, hydrogen')
  run.pymol(pymol.dir = pymol.dir, pml = pml, script = pymol)
}

#' Run pymol
#'
#' Run pymol commands from R
#' @param pymol.dir Path to pymol executable
#' @param pml Path to pml file
#' @param script Vector of String with each element containing a line of command
#' @return promt: pymol started
#' @export
run.pymol <- function(pymol.dir = 'E:/pymol/pymol_app/pyMOLWin.exe',
                      pml = '~/pymolScript.pml',
                      script){
  cat('#run pymol \n', file = pml, append = F)
  for(i in 1:length(script)){
    cat(script[i],' \n', sep = '', file = pml, append = T)
  }
  cmd <- paste(pymol.dir, ' ', pml, sep = '')
  system(cmd, wait = TRUE)
  return('Pymol started')
}

####function: JSON2Matrix (read JSON file of motif energy matrix and output as matrix)####
#' JSON to Matrix
#'
#' Parse JSON file containing binding motifs to a scoring matrix
#' @param JSON_text a string containing the raw JSON information by readLines()
#' @param expo Factor of exponential to be applied to the matrix, 0 if no exponential should be applied
#' @param mode Binding mode to be parsed
#' @return Binding matrix
#' @export
JSON2Matrix <- function(JSON_text, expo = exp(1), mode = 1){
  scoring <- RJSONIO::fromJSON(JSON_text)
  scoring_matrix <- scoring$coefficients$bindingModes[[mode + 1]]$mononucleotide
  DN_seq <- strsplit(scoring$modelSettings$letterOrder, "")[[1]]
  N_pos <- scoring$modelSettings$bindingModes[[mode+1]]$size
  if(expo < 0){
    out_matrix <- matrix(data = scoring_matrix, nrow = 4, ncol = N_pos, byrow = FALSE)
  }else{
    out_matrix <- matrix(data = expo^scoring_matrix, nrow = 4, ncol = N_pos, byrow = FALSE)
  }
  colnames(out_matrix) <- paste('P',c(1:N_pos),sep = '')
  rownames(out_matrix) <- DN_seq
  return(out_matrix)
}

#' Normalize Sum to 1
#'
#' Normalize the sum of each row or column of a matrix to 1
#' @param mat Matrix to be normalized
#' @param rows Logical that if true normalizes rows to 1, if false normalizes columns to 1
#' @return Normalized matrix
#' @export
normalize_sum1 <- function(mat, rows = TRUE){
  if(rows){
    mat <- t(mat)
  }
  for(i in 1:ncol(mat)){
    sum <- sum(mat[,i])
    for(j in 1:nrow(mat)){
      mat[j,i] <- mat[j,i] / sum
    }
  }
  if(rows){
    mat <- t(mat)
  }
  return(mat)
}

#' Concatenate Alignment table
#'
#' Concatenate .ali/.sto files from clustal omega and extract the alignment of the DNA-binding domain
#' @param Ali Data frame with 2 columns, column 1 has names as identifiers and column 2 has a string of aligned sequence
#' @param start Starting position of the DNA-binding domain
#' @param length length of the DNA-binding domain
#' @return Data frame with 2 columns: name and aligned sequence of DNA-binding domain
#' @export
concatAli <- function(Ali, start = 1, length = nchar(Ali[1,2])){
  names <- unique(Ali[,1])
  bHLH_pbAlignment <- matrix(nrow = length(names), ncol = 2)
  for(i in 1:length(names)){
    ali <- Ali[Ali[,1] == names[i],2]
    alignment <- paste(ali,collapse = '')
    add <- c(names[i], alignment)
    bHLH_pbAlignment[i,] <- add
  }
  colnames(bHLH_pbAlignment) <- c('name', 'alignment')
  bHLH_pbAlignment <- as.data.frame(bHLH_pbAlignment)
  bHLH_pbAlignment$alignment <- substr(bHLH_pbAlignment$alignment, start, (start + length-1))

  return(bHLH_pbAlignment)
}

#' Load mono_motifs
#'
#' Load binding motifs from dir containing ProBound results, the model with the highest
#' consensus sequence recognition score will be added to the resulting list
#' @param index data frame containing the identifiers of each experiment,
#' contains protein name as $gene_symbol and study name as $study in the first two column,
#' ProBound running folders should be name as $gene_symbol_$study.
#' @param modelFile_Template Path to a sample fit.models.consensus.json and switching
#' the $gene_symbol_$study identifier to $modelFile$.
#' @param rec_seq Consensus sequence to be recognized, use N for variable base positions
#' @param pos_index Position names of the consensus sequence
#' @param checkSymmetry To check symmetry like bHLH protein, should be changed for families other than bHLH
#' @param withTable If true, output with bHLH_index table with scoring and model number; if false, only output
#' @param useMode List of binding modes to use, if left empty, the highest scored mode according to rec_seq
#' will be used.
#' @return List of motifs in mono_motifs, also an info table if withTable == T
#' @export
loadMono_motifs <- function(index,
                            modelFile_Template = '~/$modelFile$/result/fit.models.consensus.json',
                            rec_seq = 'CANNTG',
                            pos_index = c('P-3','P-2','P-1','P1','P2','P3'),
                            checkSymmetry = FALSE,
                            withTable = TRUE,
                            useMode = c()){
  bHLH_index <- index
  mono_motifs <- list()
  motif_model <- c()
  symmetry <- c()
  model_score <- c()
  n <- 0
  for(i in 1:nrow(bHLH_index)){
    tryCatch({
      name <- paste0(bHLH_index$gene_symbol[i], '_',bHLH_index$study[i])
      modelFile <- gsub("\\$modelFile\\$", name, modelFile_Template)
      JSON_Lines <- readLines(modelFile)


      score <- c()
      JSON_matrix_motif <- list()
      loop <- TRUE
      m <- 1
      while(loop){
        JSON_matrix <- tryCatch({
          JSON2Matrix(JSON_Lines, mode = m)
        }, error = function(x){return(c())})
        if(length(JSON_matrix) == 0){
          loop <- FALSE
        }else{
          bind_pos <- find_binding_site(JSON_matrix, rec_seq)$max_pos[1]
          score <- c(score, max(find_binding_site(JSON_matrix, rec_seq)$scores))
          JSON_matrix_motif[[m]] <- JSON_matrix[,c(bind_pos:(bind_pos+nchar(rec_seq)-1))]
          m <- m+1
        }
      }
      if(length(useMode) == 0){
        bestIndex <- which(score == max(score))[1]
      }else{
        bestIndex <- useMode[i]
      }
      motif_model <- c(motif_model, bestIndex)
      model_score <- c(model_score, score[bestIndex])
      JSON_matrix_motif <- JSON_matrix_motif[[bestIndex]]

      colnames(JSON_matrix_motif) <- pos_index

      if(checkSymmetry){
        #check for symmetry, specific for bHLH can be deleted is not relevant
        if(abs(sum(JSON_matrix_motif[,1:(ncol(JSON_matrix_motif)/2)]) - sum(JSON_matrix_motif[,(ncol(JSON_matrix_motif)/2 + 1):(ncol(JSON_matrix_motif))])) < 1e-10){
          symmetry <- c(symmetry, 1)
        }else{
          symmetry <- c(symmetry, 0)
        }
        #end of check symmetry
      }

      n <- n + 1
      mono_motif <- list(name = name, matrix = JSON_matrix_motif)
      mono_motifs[[n]] <- mono_motif
    }, error = function(e){return(paste0('Error in ', i))})
  }
  if(checkSymmetry){
    bHLH_model_info <- data.frame(bHLH_index, motif_model, model_score, symmetry)
  }else{
    bHLH_model_info <- data.frame(bHLH_index, motif_model, model_score)
  }
  if(withTable){
    out <- list (motifs = mono_motifs, table = bHLH_model_info)
  }else{
    out <- mono_motifs
  }
  return(out)
}


#' Get base-line Accuracy
#'
#' Get the base-line accuracy for motif scoring according to different experiments performed on the same sample protein
#' @param mono_motifs List of motifs that are input as frequency values (0-1 with highest of 1)
#' @param pos Position in motif to extract
#' @param randomSample If true, randomly select from studies of the sample protein; if fasle, takes the first two studies
#' @return Data frame with two columns with pair-wise ddG values that can be plotted
#' @export
getBaseLineAccuracy <- function(mono_motifs, pos, randomSample = FALSE){
  posMotif <- pos
  func <- function(str){return(strsplit(str$name, '_')[[1]][1])}
  geneList <- unlist(lapply(mono_motifs, func))
  geneNames <- unique(geneList)

  P.1HT1 <- data.frame(placeHolder = c(0,0,0,0))
  P.1HT2 <- data.frame(placeHolder = c(0,0,0,0))
  for(i in 1:length(geneNames)){
    HTList <- grep(geneNames[i],geneList)
    if(length(HTList) < 2){
      next()
    }
    if(randomSample){
      sample <- sample(length(HTList),2)
      P.1HT1 <- cbind.data.frame(P.1HT1, mono_motifs[[HTList[sample[1]]]]$matrix[,posMotif])
      P.1HT2 <- cbind.data.frame(P.1HT2, mono_motifs[[HTList[sample[2]]]]$matrix[,posMotif])
    }else{
      P.1HT1 <- cbind.data.frame(P.1HT1, mono_motifs[[HTList[1]]]$matrix[,posMotif])
      P.1HT2 <- cbind.data.frame(P.1HT2, mono_motifs[[HTList[2]]]$matrix[,posMotif])
    }


  }

  P.1HT1 <- P.1HT1[,-1]
  P.1HT2 <- P.1HT2[,-1]

  P.1HT1 <- log(P.1HT1)
  P.1HT2 <- log(P.1HT2)
  HT1 <- unlist(data.frame(apply(P.1HT1, 2, function(column) column - mean(column))))
  HT2 <- unlist(data.frame(apply(P.1HT2, 2, function(column) column - mean(column))))

  HTdata <- data.frame(HT1, HT2)
  rownames(HTdata) <- 1:nrow(HTdata)
  return(HTdata)
}

#' Get p-value table
#'
#' Generate the p-value table of MANOVA between each amino acid residue-motif position combination
#' @param mono_motifs List of motifs that are input as frequency values (0-1 with highest of 1)
#' @param Alignment The alignment resulted from concatAli
#' @param pos_index Position names of the consensus sequence
#' @return Data frame of MANOVA p-values
#' @export
getPvalTable <- function(mono_motifs, Alignment, pos_index = c('P-3','P-2','P-1','P1','P2','P3')){
  pVal.Pos <- data.frame(NULL)
  for(AApos in 1:nchar(Alignment$alignment[1])){
    ProteinPos <- data.frame(name = Alignment$name, pos = substr(Alignment$alignment,AApos,AApos))
    PosAA <- c()
    for(i in 1:length(mono_motifs)){
      PosAA <- c(PosAA, ProteinPos[ProteinPos$name == mono_motifs[[i]]$name, 2])
    }
    if(length(unique(PosAA)) == 1){
      addLine <- NA
      pVal.Pos <- rbind.data.frame(pVal.Pos, addLine)
      next()
    }
    addLine <- c()
    for(n.pos in 1:length(pos_index)){
      pos_matrix <- gene2pos(mono_motifs, pos = pos_index[n.pos])
      df <- matrix2tetrahedron(pos_matrix)

      df <- data.frame(AA = PosAA, A = df[,1], C = df[,2], G = df[,3])
      man.test <- stats::manova(cbind(A,C,G) ~ AA, data = df)

      summary <- summary(man.test, tol=0)
      p.val <- summary$stats[1,6]
      #p.val <- man.test$`Pr(>F)`[2]
      addLine <- c(addLine, p.val)
    }
    pVal.Pos <- rbind.data.frame(pVal.Pos, addLine)
  }
  rownames(pVal.Pos) <- paste('AA', c(1:nchar(Alignment$alignment[1])), sep = '')
  colnames(pVal.Pos) <- pos_index
  return(pVal.Pos)
}

#' Amino Acid Colors
#'
#' Get the color that represent the vector of amino acids.
#' @param AAs Amino Acids to take
#' @return vector of HEX string representing the colors
#' @export
AAcolor <- function(AAs){
  AA <- c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','-')
  colors <- c('#ff9966',
               '#009999',
               '#ff0000',
               '#cc0033',
               '#00ff00',
               '#f2f20c',
               '#660033',
               '#cc9933',
               '#663300',
               '#ff9933',
               '#cc99cc',
               '#336666',
               '#0099ff',
               '#6666cc',
               '#990000',
               '#0000ff',
               '#00ffff',
               '#ffcc33',
               '#66cc66',
               '#006600',
               '#FEFEFE'
               )
  names(colors) <- AA
  return(colors[AAs])
}

#' Match Alignment and Motif
#'
#' Match the order of mono_motifs list and Alignment data frame
#' @param mono_motifs List of motifs output from loadMono_motifs
#' @param Alignment Data frame of Alignment from concatAli
#' @return A filtered and reordered Alignment table
#' @export
matchAliMotif <- function(mono_motifs, Alignment){
  geneNames <- unlist(lapply(mono_motifs, function(x) x$name))
  rownames(Alignment) <- Alignment$name
  out <- Alignment[geneNames,]
  return(out)
}

#' Amino acid residue type-base pair combination
#'
#' Creating a matrix of position-specific motif matrix at a given motif position,
#' colname labeled with the residue type at a given alignment position.
#' @param mono_motifs List of motifs output from loadMono_motifs()
#' @param Alignment Data frame of Alignment from concatAli
#' @param AApos Position of residue along the protien alignment
#' @param motifPos Position of motif
#' @return Matrix of binding motifs at a specific motif position, similar to result of gene2pos(),
#' but with residue type as colnames instead of gene name.
#' @export
AAbpCombination <- function(mono_motifs, Alignment, AApos, motifPos){
  pos_matrix <- gene2pos(mono_motifs, pos = motifPos)
  colnames(pos_matrix) <- substr(Alignment$alignment,AApos,AApos)
  return(pos_matrix)
}

#' Plot tetrahedron
#'
#' Generate a tetrahedron representation system and plot all data samples at a given motif position
#' @param posMatrix Matrix of binding motifs at a specific motif position, result of gene2pos()
#' @param base_colors Color of bases at the vertexes
#' @param size Size of dots in the tetrahedron
#' @param axis True to show axis
#' @param label Show label of data points according to colnames of posMatrix
#' @param color Color sample points according to AA identity, only when colnames of posMatrix are 1 letter AA identifiers
#' @return 3D Plotly plot
#' @export
plot_tetrahedron <- function(posMatrix, base_colors = c('green','blue','orange','red'), size = 5, axis = FALSE, label = F, color = F){
  JSON_matrix <- posMatrix
  if(color){
    resis <- unique(colnames(JSON_matrix))
    resiColors <- AAcolor(resis)
  }else{
    resiColors <- '#000000'
  }

  #transform to 3d space
  tetra_trans_matrix <- matrix(data = c(1,1,1,1,-1,-1,-1,1,-1,-1,-1,1), nrow = 4, ncol = 3, byrow = TRUE)
  df <- normalize_sum1(t(JSON_matrix)) %*% tetra_trans_matrix
  AAtype <- rownames(df)
  df <- as.data.frame(df)
  df <- cbind.data.frame(df, name = AAtype)
  colLen <- length(resiColors)
  colorNum <- as.factor(df$name)
  resiCol <- c()
  for(i in 1:length(colorNum)){
    index <- as.numeric(colorNum[i])%%colLen
    if(index == 0){
      index <- colLen
    }
    resiCol <- c(resiCol, resiColors[index])
  }

  #create tetrahedron
  tetrahedron <-  matrix(data = c(1,1,1,1,-1,-1,-1,1,-1,-1,-1,1,1,1,1,-1,1,-1,1,-1,-1,-1,-1,1), nrow = 8, ncol = 3, byrow = TRUE)
  tetrahedron <- as.data.frame(tetrahedron)
  bases <- matrix(data = c(1,1,1, 1,-1,-1, -1,1,-1, -1,-1,1), nrow = 4, ncol = 3, byrow = TRUE)
  bases <- as.data.frame(bases)
  #set parameters
  base_names <- rownames(JSON_matrix)
  len <- nrow(df)
  #set axis
  axx <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE
  )
  #plot 3d graph
  plot3D <- plot_ly()%>%
    add_trace(tetrahedron, x = tetrahedron[,1],
              y=tetrahedron[,2],z=tetrahedron[,3],color = I('blue'),
              type = 'scatter3d',
              mode = 'lines',
              line = list(width = 1),
              opacity = 0.5,
              name = 'tetrahedron')%>%
    add_trace(bases, x = as.numeric(bases[1,1]), y = as.numeric(bases[1,2]), z = as.numeric(bases[1,3]), color = I(base_colors[1]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[1])%>%
    add_trace(bases, x = as.numeric(bases[2,1]), y = as.numeric(bases[2,2]), z = as.numeric(bases[2,3]), color = I(base_colors[2]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[2])%>%
    add_trace(bases, x = as.numeric(bases[3,1]), y = as.numeric(bases[3,2]), z = as.numeric(bases[3,3]), color = I(base_colors[3]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[3])%>%
    add_trace(bases, x = as.numeric(bases[4,1]), y = as.numeric(bases[4,2]), z = as.numeric(bases[4,3]), color = I(base_colors[4]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[4])%>%
    add_trace(df, x = df[1:len,1], y=df[1:len,2],z=df[1:len,3], color = I(resiCol),
              type = 'scatter3d',
              mode = 'markers',
              size = size,
              opacity = 0.8,
              name = 'Entries',
              showlegend = F)

  if(label){
    plot3D <- plot3D%>%add_trace(df, x = df[1:len,1], y=df[1:len,2],z=df[1:len,3],
                                 type = 'scatter3d',
                                 mode = 'text',
                                 text = df[1:len, 4],
                                 textposition = 'top',
                                 name = 'label',
                                 textfont = list(color = 'grey', size = size*2),
                                 opacity = 0.5)
  }
  if(!axis){
    plot3D <- plot3D%>%layout(scene = list(xaxis = axx, yaxis = axx, zaxis = axx))
  }

  return(plot3D)
}

#' Find binding site
#'
#' Find the start position of a recognition sequence from a motif matrix and get the affinity score
#' @param JSON_matrix Matrix of binding motifs at a specific motif position, result of gene2pos()
#' @param rec_seq Consensus sequence to be recognized, use N for variable base positions
#' @param threshold Threshold for score to consider as a matching sequence
#' @return A list of two numbers. max_pos: the start position of the recognized
#' sequence with the highest score; score: the affinity score of the recognized
#' sequence at the resulting position
#' @export
find_binding_site <- function(JSON_matrix,rec_seq,threshold = 0.5){
  #set output data structure
  out <- list(start_pos = c(), high_scores = c(), scores = c(), max_pos = 0)
  #normalize JSON_matrix
  norm_matrix <- normalize_sum1(JSON_matrix, rows = FALSE)
  len_motif <- ncol(norm_matrix)
  len_rec_seq <- nchar(rec_seq)
  for(i in 1:(len_motif - len_rec_seq + 1)){
    score <- 0
    for(j in 1:len_rec_seq){
      if(substr(rec_seq,j,j) == 'N'){
        score <- score + threshold
      }else{
        score <- score + norm_matrix[substr(rec_seq,j,j),(i+j-1)]
      }
    }
    score <- score/len_rec_seq
    if(score >= threshold){
      out$start_pos <- c(out$start_pos, i)
      out$high_scores <- c(out$high_scores, score)
    }
    out$scores <- c(out$scores, score)
    out$max_pos <- grep(max(out$scores), out$scores)
  }
  return(out)
}

#' Gene to position-specific matrix
#'
#' Creating a matrix of position-specific motif matrix at a given motif position
#' @param motifs List of motifs with $matrix containing a motif matrix, result form loadMono_motif()
#' @param pos Binding matrix position to extract
#' @param nrow Number of rows in motif matrices
#' @return A position-specific motif matrix at the given motif position.
#' @export
gene2pos <- function(motifs, pos = 'P1', nrow = 4){
  ngene <- length(motifs)
  out <- matrix(nrow = nrow, ncol = ngene, data = 0)
  rownames(out) <- rownames(motifs[[1]]$matrix)
  names <- c()
  for(i in 1:ngene){
    motif <- motifs[[i]]
    names <- c(names, motif$name)
    out[,i] <- motif$matrix[,pos]
  }
  colnames(out) <- names
  return(out)
}

#' Plot 4-way tetrahedron
#'
#' Looking at a tetrahedron representation system from the 4 vertexes and plot all data samples at a given motif position
#' @param posMatrix Matrix of binding motifs at a specific motif position, result of gene2pos()
#' @param base_colors Color of bases at the vertexes
#' @param size Size of dots in the tetrahedron
#' @param label Show label of data points according to colnames of posMatrix
#' @param color Color sample points according to AA identity, only when colnames of posMatrix are 1 letter AA identifiers
#' @return 3D Plotly plot
#' @export
plot_4Graph <- function(posMatrix, base_colors = c('green','blue','orange','red'), size = 5, label = F, color = F){
  JSON_matrix <- posMatrix
  if(color){
    resis <- unique(colnames(JSON_matrix))
    resiColors <- AAcolor(resis)
  }else{
    resiColors <- '#000000'
  }

  #transform to 3d space
  tetra_trans_matrix <- matrix(data = c(1,1,1,1,-1,-1,-1,1,-1,-1,-1,1), nrow = 4, ncol = 3, byrow = TRUE)
  df <- normalize_sum1(t(JSON_matrix)) %*% tetra_trans_matrix
  AAtype <- rownames(df)
  df <- as.data.frame(df)
  df <- cbind.data.frame(df, name = AAtype)
  colLen <- length(resiColors)
  colorNum <- as.factor(df$name)
  resiCol <- c()
  for(i in 1:length(colorNum)){
    index <- as.numeric(colorNum[i])%%colLen
    if(index == 0){
      index <- colLen
    }
    resiCol <- c(resiCol, resiColors[index])
  }
  #create tetrahedron
  tetrahedron <-  matrix(data = c(1,1,1,1,-1,-1,-1,1,-1,-1,-1,1,1,1,1,-1,1,-1,1,-1,-1,-1,-1,1), nrow = 8, ncol = 3, byrow = TRUE)
  tetrahedron <- as.data.frame(tetrahedron)
  bases <- matrix(data = c(1,1,1, 1,-1,-1, -1,1,-1, -1,-1,1), nrow = 4, ncol = 3, byrow = TRUE)
  bases <- as.data.frame(bases)
  #set parameters
  base_names <- rownames(JSON_matrix)
  len <- nrow(df)
  #set axis
  axx <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE
  )
  #set cemera

  #plot 3d graph A
  plotA <- plot_ly(scene = 'scene1')%>%
    #layout(scene = list(xaxis = axx, yaxis = axx, zaxis = axx, camera = eyepos))%>%
    add_trace(tetrahedron[c(2,3,4,2),], x = tetrahedron[c(2,3,4,2),1],
              y=tetrahedron[c(2,3,4,2),2],z=tetrahedron[c(2,3,4,2),3],color = I('blue'),
              type = 'scatter3d',
              mode = 'lines',
              line = list(width = 1),
              opacity = 0.5,
              name = 'tetrahedron')%>%
    add_trace(df, x = df[1:len,1], y=df[1:len,2],z=df[1:len,3], color = I(resiCol),
              type = 'scatter3d',
              mode = 'markers',
              size = size,
              opacity = 0.8,
              name = 'Entries')%>%
    add_trace(bases, x = as.numeric(bases[2,1]), y = as.numeric(bases[2,2]), z = as.numeric(bases[2,3]), color = I(base_colors[2]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[2])%>%
    add_trace(bases, x = as.numeric(bases[3,1]), y = as.numeric(bases[3,2]), z = as.numeric(bases[3,3]), color = I(base_colors[3]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[3])%>%
    add_trace(bases, x = as.numeric(bases[4,1]), y = as.numeric(bases[4,2]), z = as.numeric(bases[4,3]), color = I(base_colors[4]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[4])
  if(label){
    plotA <- plotA%>%add_trace(df, x = df[1:len,1], y=df[1:len,2],z=df[1:len,3],
                               type = 'scatter3d',
                               mode = 'text',
                               text = df[1:len, 4],
                               textposition = 'top',
                               name = 'label',
                               textfont = list(color = 'grey', size = size*2),
                               opacity = 0.5)
  }


  #plot 3d graph C
  plotC <- plot_ly(scene = 'scene2')%>%
    #layout(scene = list(xaxis = axx, yaxis = axx, zaxis = axx, camera = eyepos))%>%
    add_trace(tetrahedron[c(1,3,4,1),], x = tetrahedron[c(1,3,4,1),1],
              y=tetrahedron[c(1,3,4,1),2],z=tetrahedron[c(1,3,4,1),3],color = I('blue'),
              type = 'scatter3d',
              mode = 'lines',
              line = list(width = 1),
              opacity = 0.5,
              name = 'tetrahedron')%>%
    add_trace(df, x = df[1:len,1], y=df[1:len,2],z=df[1:len,3], color = I(resiCol),
              type = 'scatter3d',
              mode = 'markers',
              size = size,
              opacity = 0.8,
              name = 'Entries')%>%
    add_trace(bases, x = as.numeric(bases[1,1]), y = as.numeric(bases[1,2]), z = as.numeric(bases[1,3]), color = I(base_colors[1]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[1])%>%
    add_trace(bases, x = as.numeric(bases[3,1]), y = as.numeric(bases[3,2]), z = as.numeric(bases[3,3]), color = I(base_colors[3]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[3])%>%
    add_trace(bases, x = as.numeric(bases[4,1]), y = as.numeric(bases[4,2]), z = as.numeric(bases[4,3]), color = I(base_colors[4]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[4])
  if(label){
    plotC <- plotC%>%add_trace(df, x = df[1:len,1], y=df[1:len,2],z=df[1:len,3],
                               type = 'scatter3d',
                               mode = 'text',
                               text = df[1:len, 4],
                               textposition = 'top',
                               name = 'label',
                               textfont = list(color = 'grey', size = size*2),
                               opacity = 0.5)
  }


  #plot 3d graph T
  plotT <- plot_ly(scene = 'scene4')%>%
    #layout(scene = list(xaxis = axx, yaxis = axx, zaxis = axx, camera = eyepos))%>%
    add_trace(tetrahedron[c(1,2,3,1),], x = tetrahedron[c(1,2,3,1),1],
              y=tetrahedron[c(1,2,3,1),2],z=tetrahedron[c(1,2,3,1),3],color = I('blue'),
              type = 'scatter3d',
              mode = 'lines',
              line = list(width = 1),
              opacity = 0.5,
              name = 'tetrahedron')%>%
    add_trace(df, x = df[1:len,1], y=df[1:len,2],z=df[1:len,3], color = I(resiCol),
              type = 'scatter3d',
              mode = 'markers',
              size = size,
              opacity = 0.8,
              name = 'Entries')%>%
    add_trace(bases, x = as.numeric(bases[1,1]), y = as.numeric(bases[1,2]), z = as.numeric(bases[1,3]), color = I(base_colors[1]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[1])%>%
    add_trace(bases, x = as.numeric(bases[2,1]), y = as.numeric(bases[2,2]), z = as.numeric(bases[2,3]), color = I(base_colors[2]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[2])%>%
    add_trace(bases, x = as.numeric(bases[3,1]), y = as.numeric(bases[3,2]), z = as.numeric(bases[3,3]), color = I(base_colors[3]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[3])
  if(label){
    plotT <- plotT%>%add_trace(df, x = df[1:len,1], y=df[1:len,2],z=df[1:len,3],
                               type = 'scatter3d',
                               mode = 'text',
                               text = df[1:len, 4],
                               textposition = 'top',
                               name = 'label',
                               textfont = list(color = 'grey', size = size*2),
                               opacity = 0.5)
  }


  #reverse G plot
  tetra_trans_matrix <- matrix(data = c(1,-1,-1,1,1,1,-1,-1,1,-1,1,-1), nrow = 4, ncol = 3, byrow = TRUE)
  df2 <- normalize_sum1(t(JSON_matrix)) %*% tetra_trans_matrix
  df2 <- as.data.frame(df2)
  df <- cbind.data.frame(df2, name = rownames(df2))

  plotG2 <- plot_ly(scene = 'scene3')%>%
    #layout(scene = list(xaxis = axx, yaxis = axx, zaxis = axx, camera = eyepos))%>%
    add_trace(tetrahedron[c(1,2,3,1),], x = tetrahedron[c(1,2,3,1),1],
              y=tetrahedron[c(1,2,3,1),2],z=tetrahedron[c(1,2,3,1),3],color = I('blue'),
              type = 'scatter3d',
              mode = 'lines',
              line = list(width = 1),
              opacity = 0.5,
              name = 'tetrahedron')%>%
    add_trace(df, x = df[1:len,1], y=df[1:len,2],z=df[1:len,3], color = I(resiCol),
              type = 'scatter3d',
              mode = 'markers',
              size = size,
              opacity = 0.8,
              name = 'Entries')%>%
    add_trace(bases, x = as.numeric(bases[1,1]), y = as.numeric(bases[1,2]), z = as.numeric(bases[1,3]), color = I(base_colors[2]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[1])%>%
    add_trace(bases, x = as.numeric(bases[2,1]), y = as.numeric(bases[2,2]), z = as.numeric(bases[2,3]), color = I(base_colors[1]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[2])%>%
    add_trace(bases, x = as.numeric(bases[3,1]), y = as.numeric(bases[3,2]), z = as.numeric(bases[3,3]), color = I(base_colors[4]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[3])
  if(label){
    plotG2 <- plotG2%>%add_trace(df, x = df[1:len,1], y=df[1:len,2],z=df[1:len,3],
                               type = 'scatter3d',
                               mode = 'text',
                               text = df[1:len, 4],
                               textposition = 'top',
                               name = 'label',
                               textfont = list(color = 'grey', size = size*2),
                               opacity = 0.5)
  }


  plot4 <- subplot(plotA, plotC, plotG2, plotT)%>%
    layout(scene = list(domain=list(x=c(0,0.5),y=c(0,0.5)),
                        xaxis=axx, yaxis=axx, zaxis=axx,
                        camera = list(eye = list(x = 1, y = 1, z = 1)),
                        aspectmode='cube'),
           scene2 = list(domain=list(x=c(0.25,0.75),y=c(0.12,0.62)),
                         xaxis=axx, yaxis=axx, zaxis=axx,
                         camera = list(eye = list(x = 1, y = -1, z = -1)),
                         aspectmode='cube'),
           scene3 = list(domain=list(x=c(0.25,0.75),y=c(0.43,0.93)),
                         xaxis=axx, yaxis=axx, zaxis=axx,
                         camera = list(eye = list(x = -1, y = -1, z = 1)),
                         aspectmode='cube'),
           scene4 = list(domain=list(x=c(0.5,1),y=c(0,0.5)),
                         xaxis=axx, yaxis=axx, zaxis=axx,
                         camera = list(eye = list(x = -1, y = -1, z = 1)),
                         aspectmode='cube'),
           showlegend = FALSE)
  return(plot4)
}

#' Matrix to Tetrahedron
#'
#' Change a frequency matrix into tetrahedron matrix
#' @param matrix A frequency matrix with 4 rows corresponding to the 4 bases ACGT
#' @return A tetrahedron matrix with 3 columns representing the xyz coordinates
#' @export
matrix2tetrahedron <- function(matrix){
  tetra_trans_matrix <- matrix(data = c(1,1,1,1,-1,-1,-1,1,-1,-1,-1,1), nrow = 4, ncol = 3, byrow = TRUE)
  colnames(tetra_trans_matrix) <- c('x','y','z')
  return(normalize_sum1(t(matrix)) %*% tetra_trans_matrix)
}

#' Mono-nucleotide logo
#'
#' Create the logo of a mono-nucleotide binding matrix
#' @param matrix Binding matrix in either energy, frequency, or information bit.
#' @param type One of 'energy', 'info', prob'
#' @param axes To show axes or not
#' @param reverse To show reverse strand or not
#' @param labels To show labels or not
#' @return A binding motif logo
#' @export
mononucleotide_logo <- function(matrix, type="energy", axes=TRUE, reverse=FALSE,
                                labels=TRUE) {
  suppressWarnings({
    pwm <- matrix
    letter_complement <- c('T','G','C','A')
    pwm <- apply(pwm, 2, function(column) column - mean(column))


    if (type == "info") {
      pwm <- exp(pwm)
    }

    if (type == "prob") {
      pwm <- apply(exp(pwm), 2, function(column) column / sum(column))
    }

    if (reverse) {
      pwm <- pwm[, rev(seq_len(ncol(pwm)))]

      if (is.null(dim(pwm))) {
        pwm <- matrix(pwm, length(pwm), 1)
      }

      rownames(pwm) <- letter_complement
    }
    col_scheme <- ggseqlogo::make_col_scheme(chars = c("A", "C", "G", "T"),
                                             cols = c("#5CC93B", "#0D00C4", "#F4B63F", "#BB261A"))
    alltheme   <- ggplot2::theme_bw() + ggplot2::theme(text = ggplot2::element_text(size = 23),
                                                       axis.line       = ggplot2::element_line(color = "black", size = 1),
                                                       axis.ticks      = ggplot2::element_line(color = "black", size = 1),
                                                       panel.border    = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(),
                                                       legend.position = "none")
    no_axes    <- list(ggplot2::theme(line            = ggplot2::element_blank(),
                             rect            = ggplot2::element_blank(),
                             text            = ggplot2::element_blank(),
                             axis.line       = ggplot2::element_blank(),
                             axis.text.x     = ggplot2::element_blank(),
                             axis.ticks      = ggplot2::element_blank(),
                             legend.position = "none",
                             panel.spacing   = ggplot2::unit(0, "lines"),
                             plot.margin     = ggplot2::margin(0, 0, 0, 0)),
                       ggplot2::scale_x_continuous(expand = c(0, 0)),
                       ggplot2::scale_y_continuous(expand = c(0, 0)))
    font <- "helvetica_regular"
    if (type == "energy") {
      plot <- ggseqlogo::ggseqlogo(pwm, method = "c", font = font, col_scheme = col_scheme) +
        list(ggplot2::labs(x = NULL, y = NULL),
             ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
                               alpha = 0.5, fill = "white"),
             ggplot2::scale_y_continuous(breaks = scales::pretty_breaks()),
             ggplot2::geom_hline(yintercept = 0), alltheme)

      if (labels) {
        plot <- plot + ggplot2::ylab(expression(paste(Delta, Delta, "G/RT")))
      }
    }

    if (type == "info") {
      plot <- ggseqlogo::ggseqlogo(pwm, method = "b", font = font, col_scheme = col_scheme) +
        list(ggplot2::labs(x = NULL, y = NULL),
             ggplot2::scale_y_continuous(breaks = scales::pretty_breaks()), alltheme)

      if (labels) {
        plot <- plot + ggplot2::ylab("Bits")
      }
    }

    if (type == "prob") {
      plot <- ggseqlogo::ggseqlogo(pwm, method = "c", font = font, col_scheme = col_scheme) +
        list(ggplot2::labs(x = NULL, y = NULL),
             ggplot2::scale_y_continuous(breaks = c(0.0, 0.5, 1.0), limits = c(0.0, 1.0)),
             alltheme)

      if (labels) {
        plot <- plot + ggplot2::ylab("Probability")
      }
    }

    suppressMessages(
      if (!axes) {
        plot <- plot + no_axes
      }
      else {
        plot <- plot + ggplot2::scale_x_continuous(expand = c(0.01, 0),
                                                   breaks = 1:dim(pwm)[2])
      }
    )
    return(plot)
  })

}

#' Matrix reverse
#'
#' Get the binding matrix of the reverse of a binding
#' @param matrix Binding matrix in either ddG or frequency.
#' @return The reversed binding motif matrix
#' @export
matrixReverse <- function(matrix){
  baseOrder <- rownames(matrix)
  posOrder <- colnames(matrix)
  matrix <- matrix[rev(rownames(matrix)), rev(colnames(matrix))]
  rownames(matrix) <- baseOrder
  colnames(matrix) <- posOrder
  return(matrix)
}

#get sequence for mesh plotting in tetrahedron
meshSeq <- function(n){
  combSeq <- expand.grid(c(1:n),c(1:n))
  combSeq <- combSeq[combSeq$Var1 > combSeq$Var2,]
  nrow(combSeq)
  seq <- c(n,n-1,n-2)
  continue <- TRUE
  while(nrow(combSeq) > 1){
    V1 <- which(combSeq$Var1 == seq[length(seq)-2])
    V2 <- which(combSeq$Var2 == seq[length(seq)-1])
    dele <- intersect(V1, V2)
    if(length(dele) > 0){
      combSeq <- combSeq[-dele,]
    }
    V1 <- which(combSeq$Var2 == seq[length(seq)-2])
    V2 <- which(combSeq$Var1 == seq[length(seq)-1])
    dele <- intersect(V1, V2)
    if(length(dele) > 0){
      combSeq <- combSeq[-dele,]
    }
    seq <- c(seq, as.numeric(combSeq[1,]))
    combSeq <- combSeq[-1,]
  }
  return(seq)

}

#' Tetrahedron transformation matrix
#'
#' Get the tetrahedron transformation matrix
#' @export
tetra_trans_matrix <- function(){
  return(matrix(data = c(1,1,1,1,-1,-1,-1,1,-1,-1,-1,1), nrow = 4, ncol = 3, byrow = TRUE))
}


#' Form dddG list
#'
#' Create dddG matrix of selected key positions
#' @param Alignment Data frame of Alignment from concatAli
#' @param mono_motifs List of motifs output from loadMono_motifs()
#' @param keyPos Key positions along the alignment to screen
#' @param posMotif Position in the matrix to screen
#' @param leaveOut Samples from the alignment to be left out
#' @return A list containing dddG matrices of key positions
#' @export
form.dddGList <- function(Alignment, mono_motifs, keyPos = c(13,14,5,26,8), posMotif = 'P1', leaveOut = c()){
  bHLH_pbAlignment <- Alignment
  bHLH_motifs <- mono_motifs
  bHLH_pbAlignment_backup <- bHLH_pbAlignment
  #start here
  bHLH_pbAlignment <- bHLH_pbAlignment_backup
  ittm <- MASS::ginv(tetra_trans_matrix())
  if(length(leaveOut) > 0){
    for(i in 1:length(leaveOut)){
      bHLH_pbAlignment <- bHLH_pbAlignment[-grep(leaveOut[i], bHLH_pbAlignment$name),]
    }
  }
  alignmentCheck <- nrow(bHLH_pbAlignment)
  if(alignmentCheck == 0){
    print('leaveOut not in Alignment')
    return()
  }

  bHLH_motifs_backup <- bHLH_motifs
  if(length(leaveOut > 0)){
    bHLH_motifs <- list()
    for(i in 1:length(bHLH_motifs_backup)){
      isLeave <- 0
      for(j in 1:length(leaveOut)){
        if(length(grep(leaveOut[j], bHLH_motifs_backup[[i]]$name)) > 0){
          isLeave <- 1
        }
      }
      if(isLeave == 0){
        add <- bHLH_motifs_backup[[i]]
        bHLH_motifs <- rlist::list.append(bHLH_motifs, add)
      }
    }
    if(length(bHLH_motifs) == length(bHLH_pbAlignment_backup)){
      print('nothing left out')
    }
  }

  #create dddG
  dddGList <- list()
  for(p in 1:length(keyPos)){
    AApos <- keyPos[p]
    for(p in posMotif){
      ProteinPos <- data.frame(name = bHLH_pbAlignment$name, pos = substr(bHLH_pbAlignment$alignment,AApos,AApos))
      PosAA <- c()
      for(i in 1:length(bHLH_motifs)){
        PosAA <- c(PosAA, ProteinPos[ProteinPos$name == bHLH_motifs[[i]]$name, 2])
      }
      pos_matrix <- gene2pos(bHLH_motifs, pos = p, nrow = 4)
      colnames(pos_matrix) <- PosAA
      AAs <- unique(colnames(pos_matrix))
      dddGdf <- data.frame(placeHolder = c(0,0,0,0))
      for(aa1 in AAs){
        for(aa2 in AAs){
          if(aa2 != aa1){

            mean1 <- data.frame(ddG = apply(data.frame(matrix2tetrahedron(pos_matrix[,colnames(pos_matrix) == aa1])), 2, function(column) mean(column)))
            mean2 <- data.frame(ddG = apply(data.frame(matrix2tetrahedron(pos_matrix[,colnames(pos_matrix) == aa2])), 2, function(column) mean(column)))
            mean1 <- unlist(mean1)%*%ittm
            mean2 <- unlist(mean2)%*%ittm
            mean1 <- data.frame(ddG = as.numeric(mean1))
            mean2 <- data.frame(ddG = as.numeric(mean2))

            mean1 <- mean1 - max(mean1) + 1
            mean2 <- mean2 - max(mean2) + 1
            mean1 <- log(mean1)
            mean2 <- log(mean2)
            dddG <- data.frame(mean2 - mean1)


            colnames(dddG) <- paste0(aa1, '>', aa2)
            dddGdf <- cbind.data.frame(dddGdf, dddG)

          }
        }
      }
      dddGdf <- dddGdf[,-1]
      dfName <- paste0('AA', AApos, '+', posMotif)
      dddGList[[dfName]] <- dddGdf
    }
  }

  out<- dddGList
  class(out) <- 'dddGList'

  return(out)
}

#' Predict dddGList
#'
#' Predict energy given sequence and reference energy matrix and dddGList
#' @param object dddGList model output from form.dddGList()
#' @param referenceMatrix Binding matrix of the reference sample in ddG measurements
#' @param referenceSequence Protein sequence along the alignment of the reference sample
#' @param targetSequence Protien sequence along the alignmetn of the mutant sample
#' @param keyPos Key positions along the alignment to screen, must be included the dddGList
#' @param posMotif Position in the matrix to screen, must be included in the dddGList
#' @param ... Place holder for generic function
#' @return Binding matrix of the mutant sample at the selected matrix position
#' @export
predict.dddGList <- function(object, referenceMatrix, referenceSequence, targetSequence,  keyPos = c(13,14,5,8),posMotif = 'P-1', ...){
  dddGList <- object
  AAfrom <- c()
  AAto <- c()
  for(AApos in keyPos){
    AAfrom <- c(AAfrom,substr(referenceSequence,AApos, AApos))
    AAto <- c(AAto,substr(targetSequence,AApos, AApos))
  }
  mutDt <- data.frame(keyPos, AAfrom, AAto)

  deleRows <- c()
  for(i in 1:nrow(mutDt)){
    if(mutDt$AAfrom[i] == mutDt$AAto[i]){
      deleRows <- c(deleRows, i)
    }
  }
  if(length(deleRows) == nrow(mutDt)){
    print('No difference in Key Positions')
    return(referenceMatrix)
  }else if(length(deleRows) > 0){
    mutDt <- mutDt[-deleRows,]
  }

  PredMatrix <- referenceMatrix

  for(i in 1:nrow(mutDt)){

    modifier <- 1
    if(length(grep(paste0(mutDt$AAfrom[i], '>', mutDt$AAto[i]), colnames(dddGList[[paste0('AA', mutDt$keyPos[i], '+', posMotif)]]))) != 0){
      PredMatrix <- PredMatrix + modifier * data.frame(dddGList[[paste0('AA', mutDt$keyPos[i], '+', posMotif)]][,paste0(mutDt$AAfrom[i], '>', mutDt$AAto[i])])
    }
  }
  outMatrix <- data.frame(apply(as.data.frame(PredMatrix), 2, function(column) column - mean(column)))
  return(outMatrix)
}

#' Get matrix
#'
#' Get matrix according to name from a list of motifs
#' @param motifs List of motifs, result of loadMono_motifs()
#' @param name Name of sample to retrieve
#' @param exactMatch If ture, the name has to be exactly matching the name in the motifs list;
#' if false, grep is used to search.
#' @return Motif Matrices with samples with the given name
#' @export
getMatrix <- function(motifs, name, exactMatch = T){
  if(exactMatch){
    return(motifs[[which(unlist(lapply(motifs, function(x) x$name == name)))]])
  }else{
    return(motifs[which(unlist(lapply(motifs, function(x) length(grep(name, x$name))  > 0 )))])
  }
}

#' Frequency matrix to ddG matrix
#'
#' Change motif matrix in frequency measurements to ddG measurements
#' @param matrix Motif matrix to be changed
#' @return Motif matrix in ddG measurements
#' @export
frequency2ddG <- function(matrix){
  return(apply(matrix, 2, function(x) log(x) - mean(log(x))))
}

#' ddG matrix to Frequency matrix
#'
#' Change motif matrix in ddG measurements to frequency measurements
#' @param matrix Motif matrix to be changed
#' @return Motif matrix in frequency measurements
#' @export
ddG2frequency <- function(matrix){
  return(apply(matrix, 2, function(x) exp(x - max(x))))
}

#' Weighted mean
#'
#' Compute weighted mean of a vector of numbers according to a vector of weight
#' @param vector A vector of numbers
#' @param weight weight of each number, should have the same length as vector
#' @return Weighted mean
#' @export
weightedMean <- function(vector, weight){
  w <- weight/sum(weight)
  return(sum(vector*w))
}

#' Root Mean Squared Difference
#'
#' Compute the RMSD between two vectors
#' @param a Vector 1
#' @param b Vector 2
#' @return RMSD of the 2 vectors
#' @export
RMSD <- function(a,b){
  return(mean((a-b)^2)^(1/2))
}

#' Amino Acid features
#'
#' Get Hydropathy, Volume, and Isoelectric Point features for Amino Acids
#' @return AA feature data frame
#' @export
AAfeatures <- function(){
  AA <- c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','-')
  Hydropathy <- c(1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8, 1.9, -3.5, -1.6, -3.5, -4.5, -0.8, -0.7, 4.2, -0.9, -1.3, 0)
  Volume <- c(88.6, 108.5,111.1, 138.4, 189.9, 60.1, 153.2, 166.7, 168.6, 166.7, 162.9, 114.1, 112.7, 143.8,
              173.4, 89.0, 116.1, 140.0, 227.8, 193.9, 141.275)
  PI <- c(6.11, 5.15, 2.98, 3.08, 5.76, 6.06, 7.64, 6.04, 9.47, 6.04, 5.71, 5.43, 6.30, 5.65, 10.76, 5.70,
          5.60, 6.02, 5.88, 5.63, 6.0505)

  AAfeature <- data.frame(Hydropathy, Volume, PI)
  rownames(AAfeature) <- AA
  return(AAfeature)
}

#' Matrix Singular Value Decomposition
#'
#' Compute the tetrahedron SVD matrices of a position specific motif matrix of multiple samples
#' @param posMatrix Matrix of binding motifs at a specific motif position, result of gene2pos()
#' @return u matrix, v matrix, d vector from SVD, and the mean of each column of the transformed
#' tetrahedron coordinate matrix
#' @export
matrixSVD <- function(posMatrix){
  tetra_matrix <- matrix2tetrahedron(posMatrix)
  tetra_means <- apply(tetra_matrix, 2,function(x) mean(x))
  tetra_matrix <- apply(tetra_matrix, 2,function(x) x-mean(x))
  svd <- svd(tetra_matrix)
  svd$tetra_mean <- tetra_means
  return(svd)
}

#' SVD ANOVA on residue types
#'
#' Perform ANOVA test between the u-matrix and the residue types along the alignment
#' @param svd Result from matrixSVD()
#' @param Alignment Alignment table with name and aligned sequences in the same order as the matrixSVD input
#' @return A data frame of -log10 p-values with cols = number of principal components in SVD, rows = number of positions in the alignment
#' @export
svdANOVA <- function(svd, Alignment){
  svd.aov <- matrix(nrow = ncol(svd$v), ncol = nchar(Alignment$alignment[1]), data = 0)
  for(i in 1:ncol(svd$v)){
    for(j in 1:nchar(Alignment$alignment[1])){
      u1 <- svd$u[,i]
      aa <- substr(Alignment$alignment,j,j)
      dt <- cbind.data.frame(aa,u1)
      if(length(unique(aa)) == 1){
        next()
      }
      aov <- aov(u1~aa,dt)
      svd.aov[i,j] <- summary(aov)[[1]][[5]][1]

    }
  }
  svd.aov <- -log(svd.aov,10)
  svd.aov[svd.aov == Inf] <- 0
  svd.pvalTable <- data.frame(t(svd.aov))
  svd.pvalTable$name <- c(1:nchar(Alignment$alignment[1]))
  return(svd.pvalTable)
}

#' Train SVD-regression model
#'
#' Train a SVD-regression model with matrixSVD and matched Alignment
#' @param svd Result from matrixSVD()
#' @param Alignment Alignment table with name and aligned sequences in the same order as the matrixSVD input
#' @param no.keyPos Number of key positions to use for each PC. If not provided,
#' the d vector of the svd result, rounded up, will be the number of key positions used for each PC.
#' If provided, the length of this vector should be the same as the length of the svd$d vector
#' @param keyPos A list of key positions to use for prediction, using this argument stops key position identification.
#' no.keyPos should still be entered for the number of keyPos used for each PC. If using 0 keyPos for a PC,
#' enter a valid keypos number and put 0 in the no.keyPos vector for this PC.
#' @return SVD-regression model
#' @export
trainSVD <- function(svd, Alignment, no.keyPos = c(), keyPos = c()){
  if(length(no.keyPos) == 0){
    no.keyPos <- (floor(svd$d) + 1)
  }
  if(length(svd$d) != length(no.keyPos)){
    print('no.keyPos has different length from svd, please fill unused PCs with 0')
    return()
  }
  skip <- which(no.keyPos == 0)
  alignment <- Alignment
  if(length(keyPos) == 0){
    svd.pvalTable <- svdANOVA(svd, alignment)
    keyPos <- list()
    for(i in 1:length(no.keyPos)){
      X1feature <- dplyr::arrange(svd.pvalTable, dplyr::desc(svd.pvalTable[,i]))[1:no.keyPos[i], 'name']
      sele <- rep(X1feature[1], max(no.keyPos))
      sele[1:no.keyPos[i]] <- X1feature
      keyPos[[i]] <- sele
    }
  }else{
    for(kp in 1:length(keyPos)){
      if(length(keyPos[[kp]]) < max(unlist(lapply(keyPos, function(x) length(x))))){
        keyPos[[kp]][(length(keyPos[[kp]])+1):max(unlist(lapply(keyPos, function(x) length(x))))] <- keyPos[[kp]][1]
      }
    }
  }



  #average list
  uList <- list()
  for(i in 1:length(no.keyPos)){
    rowList <- list()
    for(j in keyPos[[i]]){
      u1 <- svd$u[,i]
      aa <- substr(alignment$alignment,j,j)
      dt <- cbind.data.frame(aa,u1)
      aas <- unique(dt$aa)
      means <- c()
      for(a in aas){
        subDt <- dt[dt$aa == a,]
        means <- c(means, mean(subDt$u1))
      }
      meanDt <- data.frame(aa = aas, mean = means)
      rowList[[j]] <- meanDt
    }
    uList[[i]] <- rowList
  }

  #synthetic U matrix
  synUList <- list()
  for(ali in 1:nrow(alignment)){
    synU <- list()
    for(i in 1:length(no.keyPos)){
      addList <- c()
      for(j in keyPos[[i]]){
        aa <- substr(alignment$alignment[ali],j,j)
        uset <- uList[[i]][[j]]
        add <- uset[uset$aa == aa, 2]
        addList <- c(addList, add)
      }
      synU[[i]] <- addList
    }
    synUList[[ali]] <- synU
  }

  #get linear regression coefficients
  synUpred <- list()
  coefList <- list()
  for(uindex in 1:length(no.keyPos)){
    synthesizedU <- matrix(nrow = nrow(alignment), ncol = length(keyPos[[1]]))
    for(i in 1:nrow(alignment)){
      synthesizedU[i,] <- synUList[[i]][[uindex]]
    }
    synthesizedU <- data.frame(synthesizedU)
    synthesizedU$label <- svd$u[,uindex]
    lm <- lm(label~., synthesizedU)

    coef <- lm$coefficients
    coef[is.na(coef)] <- 0
    if(length(which(uindex == skip)) > 0){
      coef <- rep(0, length(coef))
    }
    coefList[[uindex]] <- coef

    newU <- c()
    for(i in 1:nrow(alignment)){
      add <- sum(synthesizedU[i,1:length(keyPos[[uindex]])]*coef[1:length(keyPos[[uindex]])+1]) + coef[1]
      newU <- c(newU,add)
    }
    synUpred[[uindex]] <- newU
  }


  Upred <- matrix(nrow = nrow(alignment), ncol = length(no.keyPos))
  for(i in 1:length(no.keyPos)){
    Upred[,i] <- synUpred[[i]]
  }

  reSvd <- Upred %*% diag(svd$d) %*% t(svd$v)
  true <- svd$u %*% diag(svd$d) %*% t(svd$v)
  trainingRMSD <- RMSD(reSvd, true)
  if(length(skip) != 0){
    keyPos[[skip]] <- rep(0, length(keyPos[[skip]]))
  }
  out <- list()
  out$model <- coefList
  out$svd <- svd
  out$keyPos <- keyPos
  out$trainRMSD <- trainingRMSD
  out$trainPred <- reSvd
  out$uList <- uList
  class(out) <- 'SVD'
  return(out)
}

#' Predict with SVD-regression model
#'
#' Predict binding matrix with SVD-regression model
#' @param object SVD-regression model result from trainSVD()
#' @param Alignment Alignment table with name and aligned sequences to be predicted
#' @param zero Lower limit to be assigned to predicted values
#' @param ... Place holder for generic function
#' @return Matrix of binding motif models in feaquency measurements.
#' @export
predict.SVD <- function(object, Alignment, zero = 0.001, ...){
  svdModel <- object
  if(length(svdModel$keyPos) != 3){
    print('Cannot perform tetrahedron transformation, number of PCs should be 3')
    return()
  }
  predList <- c()
  for(ali in 1:nrow(Alignment)){
    synUtest <- list()
    for(i in 1:length(svdModel$keyPos)){
      addList <- c()
      for(j in svdModel$keyPos[[i]]){
        aa <- substr(Alignment$alignment[ali],j,j)
        uset <- svdModel$uList[[i]][[j]]
        add <- uset[uset$aa == aa, 2]
        if(length(add) == 0){
          addList <- c(addList, 0)
        }else{
          addList <- c(addList, add)
        }
      }
      synUtest[[i]] <- addList
    }

    predTest <- c()
    for(i in 1:length(svdModel$keyPos)){
      predTest <- c(predTest, sum(svdModel$model[[i]][1:length(svdModel$keyPos[[i]])+1] * synUtest[[i]]) + svdModel$model[[i]][1])
    }
    pred <- predTest %*% diag(svdModel$svd$d) %*% t(svdModel$svd$v)
    predList <- c(predList, pred)
  }
  tetra_trans_matrix <- matrix(data = c(1,1,1,1,-1,-1,-1,1,-1,-1,-1,1), nrow = 4, ncol = 3, byrow = TRUE)
  pssm_trans_matrix <- MASS::ginv(tetra_trans_matrix)

  mean_matrix <- matrix(nrow = length(predList)/3, ncol = 3, data = svdModel$svd$tetra_mean, byrow = T)
  predMatrix <- matrix(nrow = length(predList)/3, ncol = 3, data = predList, byrow = T)
  predMatrix1 <- apply(t((predMatrix+mean_matrix)%*%pssm_trans_matrix + 0.25), 2, function(x) x/max(x))

  predMatrix1[predMatrix1 < zero] <- zero
  colnames(predMatrix1) <- Alignment$name
  rownames(predMatrix1) <- c('A','C','G','T')
  return(predMatrix1)
}

#' Grouped R-square
#'
#' Find average R-square value for two vectors separated in to groups with multiple members
#' @param x Vector 1
#' @param y Vector 2
#' @param member number of member of each group
#' @return Average R-square value
#' @export
groupedR2 <- function(x,y,member = 4){
  group <- member
  if(length(x)!=length(y)){
    print('x and y have different length')
    return()
  }
  if(length(x)%%group != 0){
    print('Cannot be devided evenly by group number')
    return()
  }
  xM <- matrix(nrow = length(x)/group, ncol = group, data = x, byrow = T)
  yM <- matrix(nrow = length(y)/group, ncol = group, data = y, byrow = T)
  R2s <- c()
  for(i in 1:nrow(xM)){
    dt <- data.frame(xx = xM[i,], yy = yM[i,])
    lm <- lm(yy~xx+0, dt)
    R2s <- c(R2s, suppressWarnings(summary(lm)$r.squared))
  }
  return(mean(R2s))
}

#' Closest sequence prediction
#'
#' Use the Closest sequence prediction method (Weirauch, 2014) to predict the binding matrix
#' @param mono_motifs List of motifs that are input as frequency values (0-1 with highest of 1)
#' @param Alignment The alignment resulted from concatAli
#' @param pos Binding matrix position to extract
#' @return Data frame with two columns containing the true and predicted values
#' @export
closestSeqPred <- function(mono_motifs, Alignment, pos = 'P-1'){
  pos_matrix <- gene2pos(mono_motifs, pos = pos)
  alignment <- Alignment
  trueList <- c()
  predList <- c()
  for(i in 1:ncol(pos_matrix)){
    seq <- alignment$alignment[i]
    maxStr <- Inf
    predIndex <- 0
    for(j in 1:ncol(pos_matrix)){
      if(j != i){
        Str <- stringdist::stringdist(seq, alignment$alignment[j])
        if(Str < maxStr){
          maxStr <- Str
          predIndex <- j
        }
      }
    }
    trueList <- c(trueList, unlist(data.frame(apply(log(as.data.frame(pos_matrix[,i])), 2, function(column) column - mean(column)))))
    predList <- c(predList, unlist(data.frame(apply(log(as.data.frame(pos_matrix[,predIndex])), 2, function(column) column - mean(column)))))
  }
  out <- data.frame(true = trueList, pred = predList)
  return(out)
}

#' Amino acid mapping against principal components
#'
#' Perform linear regression between the PC value of amino acid type at a given alignment position and
#' motif position and the selected amino acid property
#' @param mono_motifs List of motifs that are input as frequency values (0-1 with highest of 1)
#' @param Alignment The alignment resulted from concatAli
#' @param pos Binding matrix position to extract
#' @param AApos Residue position along the alignment
#' @param PC Principal componant to screen
#' @param property Biophysical property to screen. 1 = Hydropathy, 2 = Volume, 3 = Isoelectric Point
#' @return Plot PC value against amino acid prooperty, returns summary of linear regression
#' @export
aaPCmap <- function(mono_motifs, Alignment, pos = 'P-1', AApos = 13, PC = 1, property = 1){
  svd <- matrixSVD(gene2pos(mono_motifs, pos))
  i <- PC
  j <- AApos
  u1 <- svd$u[,i]
  p <- property
  aa <- substr(Alignment$alignment,j,j)
  AAf <- AAfeatures()
  AAf$AA <- rownames(AAf)
  dt <- cbind.data.frame(AA = aa,u1)
  dt <- merge(dt, AAf, by = 'AA')
  f <- stats::as.formula(paste0('u1~', colnames(dt)[p+2]))
  plot(dt[,p+2], dt$u1, pch = 19, xlab = colnames(AAf)[p], ylab = paste0('PC:', i), main = paste0('AA Position ', j))
  lm <- lm(f, dt)
  if(is.na(lm[1]$coefficients[2])){
    return()
  }
  graphics::abline(lm, lwd = 2, col = 'red')
  return(summary(lm))
}

#deprecated
pred.SVD.ridge <- function(alignment, tetra_matrix, testSeq = '-KSLRPLLEKRRRARINQSLSQLKGLI-L------PLLGRENS--NCSKLEKADVL', keyPos = c(1:55)){
  tetra_means <- apply(tetra_matrix, 2,function(x) mean(x))
  tetra_matrix <- apply(tetra_matrix, 2,function(x) x-mean(x))
  svd <- svd(tetra_matrix)
  trueList <- c()
  predList <- c()
  #keyPos <- unique(c(X1feature,X2feature,X3feature))
  coefTable <- data.frame(NULL)

  #average list
  uList <- list()
  for(i in 1:3){
    rowList <- list()
    for(j in keyPos){
      u1 <- svd$u[,i]
      aa <- substr(alignment$alignment,j,j)
      dt <- cbind.data.frame(aa,u1)
      aas <- unique(dt$aa)
      means <- c()
      for(a in aas){
        subDt <- dt[dt$aa == a,]
        means <- c(means, mean(subDt$u1))
      }
      meanDt <- data.frame(aa = aas, mean = means)
      rowList[[j]] <- meanDt
    }
    uList[[i]] <- rowList
  }

  #synthetic U matrix
  synUList <- list()
  for(ali in 1:nrow(alignment)){
    synU <- list()
    for(i in 1:3){
      addList <- c()
      for(j in keyPos){
        aa <- substr(alignment$alignment[ali],j,j)
        uset <- uList[[i]][[j]]
        add <- uset[uset$aa == aa, 2]
        addList <- c(addList, add)
      }
      synU[[i]] <- addList
    }
    synUList[[ali]] <- synU
  }

  #get linear regression coefficients
  synUpred <- list()
  coefList <- list()
  for(uindex in 1:3){
    synthesizedU <- matrix(nrow = nrow(alignment), ncol = length(keyPos))
    for(i in 1:nrow(alignment)){
      synthesizedU[i,] <- synUList[[i]][[uindex]]
    }
    synthesizedU <- data.frame(synthesizedU)
    synthesizedU$label <- svd$u[,uindex]
    cv.glm <- glmnet::cv.glmnet(as.matrix(synthesizedU[,-ncol(synthesizedU)]), synthesizedU$label, alpha = 0)
    best_lambda <- cv.glm$lambda.min
    best_model <- glmnet::glmnet(as.matrix(synthesizedU[,-ncol(synthesizedU)]), synthesizedU$label, alpha = 0, lambda = best_lambda)
    #print(best_lambda)
    summary(best_model)
    coef <- as.vector(coef(best_model))
    coef[coef == '.'] <- 0
    coefList[[uindex]] <- coef

    newU <- c()
    for(i in 1:nrow(alignment)){
      add <- sum(synthesizedU[i,1:length(keyPos)]*coef[1:length(keyPos)+1]) + coef[1]
      newU <- c(newU,add)
    }
    synUpred[[uindex]] <- newU
  }


  Upred <- matrix(nrow = nrow(alignment), ncol = 3)
  for(i in 1:3){
    Upred[,i] <- synUpred[[i]]
  }

  reSvd <- Upred %*% diag(svd$d) %*% t(svd$v)
  #print(RMSD(reSvd, tetra_matrix[-interation,]))

  #test
  synUtest <- list()
  for(i in 1:3){
    addList <- c()
    for(j in keyPos){
      aa <- substr(testSeq,j,j)
      uset <- uList[[i]][[j]]
      add <- uset[uset$aa == aa, 2]
      if(length(add) == 0){
        addList <- c(addList, 0)
      }else{
        addList <- c(addList, add)
      }
    }
    synUtest[[i]] <- addList
  }

  predTest <- c()
  for(i in 1:3){
    predTest <- c(predTest, sum(coefList[[i]][1:length(keyPos)+1] * synUtest[[i]]) + coefList[[i]][1])
  }
  pred <- predTest %*% diag(svd$d) %*% t(svd$v)
  #print(RMSD(pred,tetra_matrix[interation,]))


  tetra_trans_matrix <- matrix(data = c(1,1,1,1,-1,-1,-1,1,-1,-1,-1,1), nrow = 4, ncol = 3, byrow = TRUE)
  pssm_trans_matrix <- MASS::ginv(tetra_trans_matrix)

  mean_matrix <- matrix(nrow = 1, ncol = 3, data = tetra_means, byrow = T)

  #apply(t((tetra_matrix+mean_matrix)%*%pssm_trans_matrix + 0.25), 2, function(x) x/max(x)) - pos_matrix

  predMatrix <- matrix(nrow = 1, ncol = 3, data = pred, byrow = T)
  predMatrix1 <- apply(t((predMatrix+mean_matrix)%*%pssm_trans_matrix + 0.25), 2, function(x) x/max(x))

  predMatrix1[predMatrix1 < 0] <- 0.001

  pred <- unlist(data.frame(apply(log(predMatrix1), 2, function(column) column - mean(column))))
  pred <- as.vector(pred)
  return(pred)
}
