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
#' @param matrix motif matrix with row names labeled as base types
#' @return A string that contains the most preferred sequence of the scoring matrix
#' @export
matrix2seq <- function(matrix){
  scoringMatrix <- matrix
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
#' @param wait If true would wait for the pymol window to close
#' @return prompt: pymol started
#' @export
run.pymol <- function(pymol.dir = '~/pyMOLWin.exe',
                      pml = 'pymolScript.pml',
                      script,
                      wait = FALSE){
  cat('#run pymol \n', file = pml, append = F)
  for(i in 1:length(script)){
    cat(script[i],' \n', sep = '', file = pml, append = T)
  }
  if(file.exists(pymol.dir)){
    cmd <- paste(pymol.dir, ' ', getwd(), '/',pml, sep = '')
    system(cmd, wait = wait)

    return('Pymol started')
  }else{
    file.remove(pml)
    warning('Pymol executable not found. Please download pymol from https://pymol.org/2/ and assign the full path to the pymol executable to the pymol.dir argument.')
  }

}

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
#' @param matrix Matrix to be normalized
#' @param rows Logical that if true normalizes rows to 1, if false normalizes columns to 1
#' @return Normalized matrix
#' @export
normalize_sum1 <- function(matrix, rows = TRUE){
  mat <- matrix
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
#' @param Alignment Data frame with 2 columns, column 1 has names as identifiers and column 2 has a string of aligned sequence
#' @param start Starting position of the DNA-binding domain
#' @param length length of the DNA-binding domain
#' @return Data frame with 2 columns: name and aligned sequence of DNA-binding domain
#' @export
concatAli <- function(Alignment, start = 1, length = nchar(Ali[1,2])){
  Ali <- Alignment
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
#' Or can be a vector of protein names.
#' @param modelFile_Template Path to a sample fit.models.consensus.json and switching
#' the $gene_symbol_$study identifier to $modelFile$.
#' @param rec_seq Consensus sequence to be recognized, use N for variable base positions. Can be a vector of character strings
#' of the same length
#' @param pos_index Position names of the consensus sequence
#' @param weight Weight parameters to feed into find_binding_site().
#' @param threshold Threshold to evaluate recognized motifs.
#' @param rev If true, also screens the reverse matrix. Usually used for non-symmetrical motifs
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
                            weight = c(),
                            threshold = 0,
                            rev = F,
                            checkSymmetry = FALSE,
                            withTable = TRUE,
                            useMode = c()
                            ){
  if(is.null(ncol(index))){
    len <- length(index)
  }else{
    len <- nrow(index)
  }
  bHLH_index <- index
  mono_motifs <- list()
  motif_model <- c()
  symmetry <- c()
  model_score <- c()
  n <- 0
  for(i in 1:len){
    tryCatch({
      if(is.null(ncol(index))){
        name <- index[i]
      }else{
        name <- paste0(bHLH_index$gene_symbol[i], '_',bHLH_index$study[i])
      }
      modelFile <- gsub("\\$modelFile\\$", name, modelFile_Template)
      JSON_Lines <- readLines(modelFile)


      score <- c()
      JSON_matrix_motif <- list()
      loop <- TRUE

      m <- 1

      while(loop){
        JSON_matrix <- tryCatch({
          JSON2Matrix(JSON_Lines, mode = m)
        }, error = function(x){return(c())}, warning = function(w){})
        if(length(JSON_matrix) == 0){
          loop <- FALSE
        }else{
          for(rec in 1:length(rec_seq)){
            bestScore <- 0
            rec_seq1 <- rec_seq[rec]
            bindingSite <- find_binding_site(JSON_matrix, rec_seq1, threshold = threshold, weight = weight)
            bind_pos <- bindingSite$max_pos[1]
            maxScore <- max(bindingSite$scores)
            maxMatrix <- JSON_matrix
            if(rev){
              revMatrix <- matrixReverse(JSON_matrix)
              revBindingSite <- find_binding_site(revMatrix, rec_seq1, threshold = threshold, weight = weight)
              revbind_pos <- revBindingSite$max_pos[1]
              revmaxScore <- max(revBindingSite$scores)
              if(revmaxScore > maxScore){
                maxMatrix <- revMatrix
                maxScore <- revmaxScore
                bind_pos <- revbind_pos
              }
            }

            if((bind_pos+nchar(rec_seq1)-1) > ncol(maxMatrix) || bind_pos < 1){
              next()
            }
            if(maxScore > bestScore ){
              bestScore <- maxScore
              bestPos <- bind_pos
              bestMatrix <- maxMatrix[,c(bestPos:(bestPos+nchar(rec_seq1)-1))]
            }
          }
          score <- c(score, bestScore)
          JSON_matrix_motif[[m]] <- bestMatrix
          m <- m+1
        }
      }
      if(length(score) == 0){
        motif_model <- c(motif_model, NA)
        model_score <- c(model_score, NA)
        symmetry <- c(symmetry,NA)
        next()
      }
      if(length(useMode) == 0){
        bestIndex <- which(score == max(score))[1]
      }else{
        bestIndex <- useMode[i]
      }
      motif_model <- c(motif_model, bestIndex)
      model_score <- c(model_score, score[bestIndex])
      JSON_matrix_motif <- JSON_matrix_motif[[bestIndex]]

      colnames(JSON_matrix_motif) <- pos_index[1:ncol(JSON_matrix_motif)]

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
  if(length(which(is.na(bHLH_model_info$model_score))) != 0){
    mono_motifs <- mono_motifs[-which(is.na(bHLH_model_info$model_score))]
    bHLH_model_info <- bHLH_model_info[-which(is.na(bHLH_model_info$model_score)),]
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
#' @param mono_motifs List of motifs that are input as frequency values (0-1 with highest of 1), each element of the
#' list has two slots, $name containing protienName_study, $matrix containing the scoring matrix
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
      PosAA <- c(PosAA, ProteinPos[ProteinPos$name == mono_motifs[[i]]$name, 2][1])
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
      df[is.na(df)] <- 0
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
  AA <- c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','-')
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
               '#000000',
               '#006600',
               '#EEEEEE'
               )
  names(colors) <- AA
  return(colors[AAs])
}

#' Match Alignment and Motif
#'
#' Match the order of mono_motifs list and Alignment data frame
#' @param mono_motifs List of motifs output from loadMono_motifs
#' @param Alignment Data frame of Alignment from concatAli
#' @param both If True both mono_motifs list and Alignment data frame is ouput; if False
#' only Alignment data frame is output
#' @return A filtered and reordered Alignment table (and mono_motif list)
#' @export
matchAliMotif <- function(mono_motifs, Alignment, both = F){
  geneNames <- unlist(lapply(mono_motifs, function(x) x$name))
  dele <- c()
  for(i in 1:length(geneNames)){
    if(length(which(geneNames[i] == Alignment$name)) == 0){
      dele <- c(dele, i)
    }
  }
  if(length(dele) > 0){
    mono_motifs <- mono_motifs[-dele]
  }
  geneNames <- unlist(lapply(mono_motifs, function(x) x$name))
  out <- list()
  ali <- data.frame()
  for(j in 1:length(geneNames)){
    ali <- rbind.data.frame(ali, unlist(Alignment[Alignment$name == geneNames[j], ][1,]))
  }
  colnames(ali) <- c('name', 'alignment')
  out$alignment <- ali
  out$motifs <- mono_motifs
  if(!both){
    return(out$alignment)
  }else{
    return(out)
  }
}

#' Amino acid residue type-base pair combination
#'
#' Creating a matrix of position-specific motif matrix at a given motif position,
#' colname labeled with the residue type at a given alignment position.
#' @param mono_motifs List of motifs output from loadMono_motifs()
#' @param Alignment Data frame of Alignment from concatAli
#' @param AApos Position of residue along the protein alignment
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
#' @importFrom plotly add_text
#' @param posMatrix Matrix of binding motifs in frequency measurements at a specific motif position, result of gene2pos()
#' @param base_colors Color of bases at the vertexes
#' @param size Size of dots in the tetrahedron
#' @param axis True to show axis
#' @param label Show label of data points according to colnames of posMatrix
#' @param color Color sample points according to AA identity, only when colnames of posMatrix are 1 letter AA identifiers
#' @param vertex Show label of vertex
#' @return 3D Plotly plot
#' @export
plot_tetrahedron <- function(posMatrix, base_colors = c('green','blue','orange','red'), size = 5, axis = FALSE, label = F, color = F, vertex = F){
  oldw <- getOption("warn")
  options(warn = -1)
  JSON_matrix <- posMatrix
  if(color){
    resis <- unique(as.factor(colnames(JSON_matrix)))
    resis <- as.character(levels(resis))
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
              name = 'tetrahedron',
              showlegend = F)%>%
    add_trace(bases, x = as.numeric(bases[1,1]), y = as.numeric(bases[1,2]), z = as.numeric(bases[1,3]), color = I(base_colors[1]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[1],
              showlegend = F)%>%
    add_trace(bases, x = as.numeric(bases[2,1]), y = as.numeric(bases[2,2]), z = as.numeric(bases[2,3]), color = I(base_colors[2]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[2],
              showlegend = F)%>%
    add_trace(bases, x = as.numeric(bases[3,1]), y = as.numeric(bases[3,2]), z = as.numeric(bases[3,3]), color = I(base_colors[3]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[3],
              showlegend = F)%>%
    add_trace(bases, x = as.numeric(bases[4,1]), y = as.numeric(bases[4,2]), z = as.numeric(bases[4,3]), color = I(base_colors[4]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[4],
              showlegend = F)%>%
    add_trace(df, x = df[1:len,1], y=df[1:len,2],z=df[1:len,3], color = I(resiCol),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size/2, opacity = 0.8),
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
  if(vertex){
    plot3D <- plot3D%>%add_text(x = bases[,1], y = bases[,2], z = bases[,3], text = c('A','C', 'G', 'T'),color = I(base_colors[1:4]),
             textposition = 'top right', showlegend = F, textfont = list(family="sans serif",size = size*3)
    )
  }
  if(!axis){
    plot3D <- plot3D%>%layout(scene = list(xaxis = axx, yaxis = axx, zaxis = axx))
  }
  options(warn = oldw)
  return(plot3D)
}

#' Find binding site
#'
#' Find the start position of a recognition sequence from a motif matrix and get the affinity score
#' @param JSON_matrix Matrix of binding motifs in frequency measurements at a specific motif position, result of gene2pos()
#' @param rec_seq Consensus sequence to be recognized, use N for variable base positions
#' @param threshold Threshold for score to consider as a matching sequence
#' @param weight Weight of each position in rec_seq. Larger number resemble higher importance.
#' @return A list of two numbers. max_pos: the start position of the recognized
#' sequence with the highest score; score: the affinity score of the recognized
#' sequence at the resulting position
#' @export
find_binding_site <- function(JSON_matrix,rec_seq,threshold = 0.5, weight = c()){
  if(length(weight) == 0){
    weight <- rep(1,nchar(rec_seq))
  }
  if(length(weight) != nchar(rec_seq)){
    print('Weight length different from rec_seq')
    return(NA)
  }
  #set output data structure
  out <- list(start_pos = c(), high_scores = c(), scores = c(), max_pos = 0)
  #normalize JSON_matrix
  norm_matrix <- normalize_sum1(JSON_matrix, rows = F)
  len_motif <- ncol(norm_matrix)
  len_rec_seq <- nchar(rec_seq)
  for(i in 1:(len_motif - len_rec_seq + 1)){
    score <- 0
    for(j in 1:len_rec_seq){
      if(substr(rec_seq,j,j) == 'N'){
        score <- score + threshold*weight[j]
      }else{
        score <- score + norm_matrix[substr(rec_seq,j,j),(i+j-1)]*weight[j]
      }
    }
    score <- score/sum(weight)
    if(score >= threshold){
      out$start_pos <- c(out$start_pos, i)
      out$high_scores <- c(out$high_scores, score)
    }
    out$scores <- c(out$scores, score)
  }
  out$max_pos <- out$start_pos[which(out$high_scores == max(out$high_scores))]
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
#' @param posMatrix Matrix of binding motifs in frequency measurements at a specific motif position, result of gene2pos()
#' @param base_colors Color of bases at the vertexes
#' @param size Size of dots in the tetrahedron
#' @param label Show label of data points according to colnames of posMatrix
#' @param color Color sample points according to AA identity, only when colnames of posMatrix are 1 letter AA identifiers
#' @return 3D Plotly plot
#' @export
plot_4Graph <- function(posMatrix, base_colors = c('green','blue','orange','red'), size = 5, label = F, color = F){
  oldw <- getOption("warn")
  options(warn = -1)
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
  options(warn = oldw)
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
        plot <- plot + ggplot2::ylab(expression(paste('-',Delta, Delta, "G/RT")))
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
#' Get the binding matrix of the reverse complement of a binding motif
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
#' @param mono_motifs List of motifs in frequency measurements output from loadMono_motifs()
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
#' @param posMatrix Matrix of binding motifs in frequency measurements at a specific motif position, result of gene2pos()
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
#' @importFrom stats anova lm predict
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
  out <- list()
  out$model <- coefList
  out$svd <- svd
  out$keyPos <- keyPos
  out$trainRMSD <- trainingRMSD
  out$trainPred <- reSvd
  out$uList <- uList
  out$alignment <- Alignment
  class(out) <- 'SVD'
  return(out)
}

#' BLOSUM62
#'
#' protein substitution matrix BLOSUM62 used in svd prediction
#'
#' @name BLOSUM62
#' @docType data
#' @author H. Pagès and P. Aboyoun
#' @references \url{https://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/Biostrings/html/substitution_matrices.html}
#' @keywords imported BLOSUM62 matrix from K. Malde, The effect of sequence quality on sequence alignment, Bioinformatics, Feb 23, 2008.
'BLOSUM62'

#' Predict with SVD-regression model
#'
#' Predict binding matrix with SVD-regression model
#' @param object SVD-regression model result from trainSVD()
#' @param Alignment Alignment table with name and aligned sequences to be predicted
#' @param zero Lower limit to be assigned to predicted values
#' @param useSimilarAA Use the weighted sum of similar AA when no matching AA is found in training data
#' @param ... Place holder for generic function
#' @return Matrix of binding motif models in frequency measurements. attr('confidence') give confidence estimates of prediction
#' @export
predict.SVD <- function(object, Alignment, zero = 0.001, useSimilarAA = F, ...){
  svdModel <- object
  if(length(svdModel$keyPos) != 3){
    print('Cannot perform tetrahedron transformation, number of PCs should be 3')
    return()
  }
  predList <- c()
  confList <- c()
  for(ali in 1:nrow(Alignment)){
    synUtest <- list()
    confidences <- list()
    for(i in 1:length(svdModel$keyPos)){
      addList <- c()
      confidenceList <- c()
      for(j in svdModel$keyPos[[i]]){
        aa <- substr(Alignment$alignment[ali],j,j)
        uset <- svdModel$uList[[i]][[j]]
        add <- uset[uset$aa == aa, 2]
        if(length(add) == 0){
          if(!useSimilarAA){
            addList <- c(addList, 0)
            confidenceList <- c(confidenceList,0)
          }else{
            #weighted mean of uset, weighted by 1/(non-replacement blosum62 score - replacement blosum62 score)
            sim <- BLOSUM62[aa,uset$aa]
            max <- BLOSUM62[aa,aa]
            sim <- 1/(max - sim)
            add <- sum(uset$mean * sim) / sum(sim)
            addList <- c(addList, add)
            confidenceList <- c(confidenceList, max(sim)/max(table(substr(svdModel$alignment$alignment,j,j))))
          }

        }else{
          addList <- c(addList, add)
          confidenceList <- c(confidenceList, sum(substr(svdModel$alignment$alignment,j,j)==aa)/max(table(substr(svdModel$alignment$alignment,j,j))))
        }
      }
      synUtest[[i]] <- addList
      confidences[[i]] <- confidenceList
    }

    predTest <- c()
    for(i in 1:length(svdModel$keyPos)){
      predTest <- c(predTest, sum(svdModel$model[[i]][1:length(svdModel$keyPos[[i]])+1] * synUtest[[i]]) + svdModel$model[[i]][1])
    }
    pred <- predTest %*% diag(svdModel$svd$d) %*% t(svdModel$svd$v)
    predList <- c(predList, pred)
    #weighted sum for confidence sum of confidence * beta * d-vector
    confidence <- 0
    dev <- 0
    for(i in 1:length(confidences)){
      add <- sum(confidences[[i]] * abs(svdModel$model[[i]][-1]))# * svdModel$svd$d[i] / max(svdModel$svd$d))
      dev <- dev + sum(svdModel$model[[i]][-1])
      confidence <- sum(confidence, add)
    }
    confList <- c(confList, confidence/dev)
  }
  tetra_trans_matrix <- matrix(data = c(1,1,1,1,-1,-1,-1,1,-1,-1,-1,1), nrow = 4, ncol = 3, byrow = TRUE)
  pssm_trans_matrix <- MASS::ginv(tetra_trans_matrix)

  mean_matrix <- matrix(nrow = length(predList)/3, ncol = 3, data = svdModel$svd$tetra_mean, byrow = T)
  predMatrix <- matrix(nrow = length(predList)/3, ncol = 3, data = predList, byrow = T)
  predMatrix1 <- apply(t((predMatrix+mean_matrix)%*%pssm_trans_matrix + 0.25), 2, function(x) x/max(x))

  predMatrix1[predMatrix1 < zero] <- zero
  colnames(predMatrix1) <- Alignment$name
  rownames(predMatrix1) <- c('A','C','G','T')
  attr(predMatrix1,'confidence') <- confList
  return(predMatrix1)
}

#' Grouped R-square
#'
#' Find average R-square value for two vectors separated in to groups with multiple members
#' @param x Vector 1
#' @param y Vector 2
#' @param member number of member of each group
#' @param mean If true, returns the mean of all groups; if false, returns a vector of each R-square value
#' @param throughZero If true the regression line passes through zero
#' @return Average R-square value
#' @export
groupedR2 <- function(x,y,member = length(x), mean = T, throughZero = T){
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
    if(throughZero){
      lm <- lm(yy~xx+0, dt)
    }else{
      lm <- lm(yy~xx, dt)
    }

    if(is.na(suppressWarnings(summary(lm)$r.squared))){
      add <- 0
    }
    add <- suppressWarnings(summary(lm)$r.squared)
    if(length(lm$coefficients[1]) == 0){
      add <- 0
    }else if(lm$coefficients[1] <= 0){
      add <- 0
    }
    R2s <- c(R2s, add)
  }
  if(mean){
    return(mean(R2s))
  }else{
    return(R2s)
  }

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
  similarity <- c()
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
    similarity <- c(similarity, rep(maxStr,4))
    trueList <- c(trueList, unlist(data.frame(apply(log(as.data.frame(pos_matrix[,i])), 2, function(column) column - mean(column)))))
    predList <- c(predList, unlist(data.frame(apply(log(as.data.frame(pos_matrix[,predIndex])), 2, function(column) column - mean(column)))))
  }
  out <- data.frame(true = trueList, pred = predList, similarity = similarity)
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
  dt <- dt[dt$AA != '-',]
  f <- stats::as.formula(paste0('u1~', colnames(dt)[p+2]))
  plot(dt[,p+2], dt$u1, pch = 19, xlab = colnames(AAf)[p], ylab = paste0('PC:', i), main = paste0('AA Position ', j))
  lm <- lm(f, dt)
  if(is.na(lm[1]$coefficients[2])){
    return()
  }
  graphics::abline(lm, lwd = 2, col = 'red')
  return(summary(lm))
}

#' Motif list to ddG matrix list
#'
#' Change motif matrices in a list to ddG measurements
#' @param mono_motifs List of motif matrices
#' @return List of mono_motifs. A warning will be print if the mono_motifs list does not
#' have the $name tag.
#' @export
motifList2ddG <- function(mono_motifs){
  if(length(mono_motifs[[1]]$name) == 0){
    print('No name tag in found List')
  }
  temp <- function(x){
    if(sum(x$matrix[,1]) > 1){
      x$matrix <- frequency2ddG(x$matrix)
    }
    return(x)
  }
  out <- lapply(mono_motifs, temp)
  return(out)
}

#' Motif list to frequency matrix list
#'
#' Change motif matrices in a list to frequency measurements
#' @param mono_motifs List of motif matrices
#' @return List of mono_motifs. A warning will be print if the mono_motifs list does not
#' have the $name tag.
#' @export
motifList2frequency <- function(mono_motifs){
  if(length(mono_motifs[[1]]$name) == 0){
    print('No name tag in found List')
  }
  temp <- function(x){
    if(sum(x$matrix[,1]) < 1){
      x$matrix <- ddG2frequency(x$matrix)
    }
    return(x)
  }
  out <- lapply(mono_motifs, temp)
  return(out)
}

#' Inspect feature SVD
#'
#' Inspect the correlation between amino acid properties and SVD PC values
#' @param svd Output from matrixSVD()
#' @param Alignment The alignment resulted from concatAli, ideally after matchAliMotif().
#' @return List of 3 matrices, each matrix records the p-value of linear regression
#'  between SVD PCs and one of the three AAfeatures.
#' @export
inspectFeatureSVD <- function(svd, Alignment){
  AAfeature <- AAfeatures()
  properties <- c('Hydropathy', 'Volumn', 'Isoelectic Point')
  PCs <- c('PC1','PC2','PC3')
  AAfeature$AA <- rownames(AAfeature)
  upRsqMatrixList <- list()
  for (p in 1:3){
    upRsqMatrix <- matrix(nrow = nchar(Alignment$alignment[1]), ncol = 3, data = 0)
    for(i in 1:3){
      for(j in 1:nchar(Alignment$alignment[1])){
        u1 <- svd$u[,i]
        aa <- substr(Alignment$alignment,j,j)
        dt <- cbind.data.frame(AA = aa,u1)
        dt <- merge(dt, AAfeature, by = 'AA')
        f <- stats::as.formula(paste0('u1~', colnames(dt)[p+2]))
        lm <- lm(f, dt)
        if(is.na(lm[1]$coefficients[2])){
          next()
        }
        summary <- summary(lm(f, dt))
        if(nrow(summary$coefficients) == 1){
          next()
        }
        rsq <- -log(summary$coefficients[2,4], 10)
        upRsqMatrix[j,i] <- rsq
      }
    }
    upRsqMatrixList[[properties[p]]] <- upRsqMatrix
  }
  return(upRsqMatrixList)
}


#' SVD-regression Cross-validation
#'
#' Leave-one-out cross-validation for SVD-regression model
#' @param mono_motifs Motif samples to incldue in the CV test, result from loadMono_motifs() or filterMotifList().
#' @param Alignment The alignment resulted from concatAli, ideally after matchAliMotif().
#' @param pos Position of the binding matrix to predict.
#' @param useSimilarAA Use AA similarity assesment to fill unseen AAs
#' @return A data frame with two columns of true and predicted values
#' @export
SVDregression.CV <- function(mono_motifs, Alignment, pos = 'P-1', useSimilarAA = F){
  trainSVDModel <- trainSVD
  HDmotifs <- mono_motifs
  HDAlignment <- Alignment
  posM <- pos
  no.keyPos <- matrix(nrow = length(HDmotifs), ncol = 3, data = 0)
  predTrue <- data.frame(NULL)
  for(t in 1:length(HDmotifs)){
    train_motifs <- HDmotifs[-t]
    test_motifs <- HDmotifs[t]
    train_alignment <- matchAliMotif(train_motifs, HDAlignment, both = T)$alignment
    test_alignment <- matchAliMotif(test_motifs, HDAlignment, both = T)$alignment
    #Train svd-regression model
    svd <- matrixSVD(gene2pos(train_motifs, pos = pos))
    svdModel <- trainSVDModel(svd, train_alignment)
    no.keyPos[t,] <- unlist(lapply(svdModel$keyPos,function(x) length(x)))
    #Predict binding motifs for test set
    pred_motifs <- predict(svdModel, test_alignment, zero = 0.01, useSimilarAA = useSimilarAA)
    #Comparing between true and predicted testing set motifs
    true <- frequency2ddG(gene2pos(test_motifs, pos = pos))
    pred <- frequency2ddG(pred_motifs)
    add <- data.frame(true = unlist(as.numeric(true)), pred = unlist(as.numeric(pred)))
    predTrue <- rbind.data.frame(predTrue, add)
  }
  attr(predTrue, 'no.keyPos') <- no.keyPos
  return(predTrue)
}

#' Amino acid type ploted against principal components in box plot
#'
#' Make box plot of amino acid types and PC values
#' @param mono_motifs List of motifs that are input as frequency values (0-1 with highest of 1)
#' @param Alignment The alignment resulted from concatAli
#' @param pos Binding matrix position to extract
#' @param AApos Residue position along the alignment
#' @param PC Principal componant to screen
#' @return box plot
#' @export
aaPCboxPlot <- function(mono_motifs, Alignment, pos = 'P-1', AApos = 13, PC = 1){
  svd <- matrixSVD(gene2pos(mono_motifs, pos))
  i <- PC
  j <- AApos
  u1 <- svd$u[,i]
  aa <- substr(Alignment$alignment,j,j)
  AAf <- AAfeatures()
  AAf$AA <- rownames(AAf)
  dt <- cbind.data.frame(AA = aa,u1)
  dt <- dt[dt$AA != '-',]
  cols <- AAcolor(as.character(levels(as.factor(dt$AA))))
  plot <- graphics::boxplot(u1~AA,
          data=dt,
          main=paste0('AA Position ', j),
          xlab="Amino Acid resiue type",
          ylab=paste0('PC:', i),
          col=cols,
          border="black"
  )
  return(plot)
}

#' Train SVD-regression model with Iterative method
#'
#' Train a SVD-regression model with matrixSVD and matched Alignment with Iterative method
#' @param svd Result from matrixSVD()
#' @param Alignment Alignment table with name and aligned sequences in the same order as the matrixSVD input
#' @param Ftest_pVal Threshold for F-test
#' @return SVD-regression model
#' @export
trainSVD.Iterative <- function(svd, Alignment, Ftest_pVal = 0.1){
  nf <- nchar(Alignment$alignment[1])
  no.keyPos <- c(nf,nf,nf)
  alignment <- Alignment
  keypos <- c()
  svd.pvalTable <- svdANOVA(svd, alignment)
  keyPos <- list()
  for(i in 1:length(no.keyPos)){
    X1feature <- dplyr::arrange(svd.pvalTable, dplyr::desc(svd.pvalTable[,i]))[1:no.keyPos[i], 'name']
    sele <- rep(X1feature[1], max(no.keyPos))
    sele[1:no.keyPos[i]] <- X1feature
    keyPos[[i]] <- sele
  }
  #average list
  uList <- list()
  for(i in 1:length(no.keyPos)){
    rowList <- list()
    for(j in 1:nchar(alignment$alignment[1])){
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
      for(j in 1:nchar(alignment$alignment[1])){
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
  keyPosAdd <- list()
  for(uindex in 1:length(no.keyPos)){
    synthesizedU <- matrix(nrow = nrow(alignment), ncol = nchar(alignment$alignment[1]))
    for(i in 1:nrow(alignment)){
      synthesizedU[i,] <- synUList[[i]][[uindex]]
    }
    synthesizedU <- data.frame(synthesizedU)
    synthesizedU$label <- svd$u[,uindex]

    addedPos <- keyPos[[uindex]][1]
    for(ff in 1:(ncol(synthesizedU)-2)){
      trainEpoc <- synthesizedU[,c(addedPos,ncol(synthesizedU))]
      lm1 <- lm(label~., trainEpoc)
      nextPos <- keyPos[[uindex]][1+ff]
      nextEpoc <- synthesizedU[,c(addedPos,nextPos,ncol(synthesizedU))]
      lm2 <- lm(label~., nextEpoc)
      Ftest <- anova(lm1,lm2)
      pFtest <- Ftest$`Pr(>F)`[2]
      if(is.na(pFtest)){
        break()
      }
      if(pFtest < Ftest_pVal){
        addedPos <- c(addedPos, nextPos)
      }else{
        break()
      }
    }

    coef <- lm1$coefficients
    coef[is.na(coef)] <- 0
    coefList[[uindex]] <- coef

    newU <- c()
    for(i in 1:nrow(alignment)){
      add <- sum(synthesizedU[i,addedPos]*coef[1:length(addedPos)+1]) + coef[1]
      newU <- c(newU,add)
    }
    synUpred[[uindex]] <- newU
    keyPosAdd[[uindex]] <- addedPos
  }


  Upred <- matrix(nrow = nrow(alignment), ncol = length(no.keyPos))
  for(i in 1:length(no.keyPos)){
    Upred[,i] <- synUpred[[i]]
  }

  reSvd <- Upred %*% diag(svd$d) %*% t(svd$v)
  true <- svd$u %*% diag(svd$d) %*% t(svd$v)
  trainingRMSD <- RMSD(reSvd, true)
  out <- list()
  out$model <- coefList
  out$svd <- svd
  out$keyPos <- keyPosAdd
  out$trainRMSD <- trainingRMSD
  out$trainPred <- reSvd
  out$uList <- uList
  out$alignment <- Alignment
  class(out) <- 'SVD'
  return(out)
}

#' SVD-regression Cross-validation for Iterative method
#'
#' Leave-one-out cross-validation for SVD-regression model with Iterative feature selection method
#' @param mono_motifs Motif samples to incldue in the CV test, result from loadMono_motifs().
#' @param Alignment The alignment resulted from concatAli, ideally after matchAliMotif().
#' @param pos Position of the binding matrix to predict.
#' @param Ftest_pVal Threshold for F-test, only used when iterative is True.
#' @param zero Lower limit to be assigned to predicted values
#' @param useSimilarAA Use AA similarity assesment to fill unseen AAs
#' @return A data frame with two columns of true and predicted values
#' @export
SVDregression.Iterative.CV <- function(mono_motifs, Alignment, pos = 'P-1', Ftest_pVal = 0.001, zero = 0.01, useSimilarAA = F){

  trainSVDModel <- trainSVD.Iterative

  HDmotifs <- mono_motifs
  HDAlignment <- Alignment
  posM <- pos
  no.keyPos <- matrix(nrow = length(HDmotifs), ncol = 3, data = 0)
  predTrue <- data.frame(NULL)
  confidences <- c()
  for(t in 1:length(HDmotifs)){
    train_motifs <- HDmotifs[-t]
    test_motifs <- HDmotifs[t]
    train_alignment <- matchAliMotif(train_motifs, HDAlignment, both = T)$alignment
    test_alignment <- matchAliMotif(test_motifs, HDAlignment, both = T)$alignment
    #Train svd-regression model
    svd <- matrixSVD(gene2pos(train_motifs, pos = pos))
    svdModel <- trainSVDModel(svd, train_alignment, Ftest_pVal = Ftest_pVal)
    no.keyPos[t,] <- unlist(lapply(svdModel$keyPos,function(x) length(x)))
    #Predict binding motifs for test set
    pred_motifs <- predict(svdModel, test_alignment, zero = zero, useSimilarAA = useSimilarAA)
    #Comparing between true and predicted testing set motifs
    true <- frequency2ddG(gene2pos(test_motifs, pos = pos))
    pred <- frequency2ddG(pred_motifs)
    add <- data.frame(true = unlist(as.numeric(true)), pred = unlist(as.numeric(pred)))
    predTrue <- rbind.data.frame(predTrue, add)
    confidences <- c(confidences, attr(pred_motifs,'confidence'))
  }
  attr(predTrue, 'no.keyPos') <- no.keyPos
  attr(predTrue,'confidence') <- confidences
  return(predTrue)
}

#' Tetrahedron to Matrix
#'
#' Transform tetrahedron coordinate representation of PSSM to PSSM
#' @param tetraMatrix Tetrahedron coordinates, result from matrix2tetrahedron
#' @return A PSSM
#' @export
tetrahedron2matrix <- function(tetraMatrix){
  tetra_trans_matrix <- matrix(data = c(1,1,1,1,-1,-1,-1,1,-1,-1,-1,1), nrow = 4, ncol = 3, byrow = TRUE)
  pssm_trans_matrix <- MASS::ginv(tetra_trans_matrix)
  predMatrix1 <- apply(t((tetraMatrix)%*%pssm_trans_matrix + 0.25), 2, function(x) x/max(x))
  return(predMatrix1)
}

#' Predict PSAM
#'
#' predict for multiple positions and output PSAM
#' @param train_motifs motifs for training, result from loadMono_motifs()
#' @param train_alignment alignment for training, result for matchAliMotif()
#' @param predSeq sequences to predict, matching alignment of train_alignment
#' @param pos positions to show in the PSAM
#' @param zero Lower limit to be assigned to predicted values
#' @param useSimilarAA Use the weighted sum of similar AA when no matching AA is found in training data
#' @return A PSAM in frequency measruements
#' @export
SVDpredPSAM <- function(train_motifs, train_alignment, predSeq, pos = c('P-3','P-2','P-1','P1','P2','P3'), zero = 0.001, useSimilarAA = T){
  test_alignment <- data.frame(name = 'Sample', alignment = predSeq)
  PSAM <- matrix(nrow = 4, ncol = length(pos), data = 1)
  rownames(PSAM) <- c('A', 'C', 'G','T')
  colnames(PSAM) <- pos
  for(p in 1:length(pos)){
    thisPos <- pos[p]
    svd <- matrixSVD(gene2pos(train_motifs, pos = thisPos))
    svdModel <- trainSVD(svd, train_alignment)
    #Predict binding motifs for test set
    pred_motifs <- predict(svdModel, test_alignment, zero = zero, useSimilarAA = useSimilarAA)
    PSAM[,p] <- pred_motifs
  }
  return(PSAM)
}

#' DNA
#'
#' DNA rownames
#' @return DNA seq
#' @export
DNA <- function(){
  return(c('A','C','G','T'))
}

#' Score Matrix Seed
#'
#' Score a motif matrix with a seed, the matrix and seed needs to be in the same length
#' @param matrix motif matrix in frequency measurements
#' @param seed seed matrix describing binding motif.
#' @param weight numeric vector of weights of each motif position of the seed
#' @return Minimum distance score
#' @export
scoreMatrixSeed <- function(matrix, seed, weight = c()){
  matrix <- normalize_sum1(matrix, rows = F)
  if(length(weight) == 0){
    weight <- rep(1,ncol(seed))
  }else if(length(weight) != ncol(seed)){
    print('Weight length not matched')
    return()
  }
  if(ncol(matrix) != ncol(seed)){
    print('Matrix and Seed column number not matched')
    return()
  }
  score <- 0
  c <- 0
  for(i in 1:ncol(seed)){
    if(sum(seed[,i]) == 0){
      next()
    }
    posScores <- c()
    for(j in 1:nrow(seed)){
      if(seed[j,i] != 0){
        add <- (seed[j,i] - matrix[j,i])^2
        posScores <- c(posScores,add)
      }
    }
    score <- score + min(posScores)* weight[i]
    c <- c+1
  }
  score <- (score/c)^(1/2)
  return(score)
}

#' Score Motif Seed
#'
#' Score a motif matrix with a seed to find the position of a motif
#' @param matrix motif matrix in frequency measurements
#' @param seed seed matrix describing binding motif.
#' @param weight numeric vector of weights of each motif position of the seed
#' @return Minimum distance score and position
#' @export
scoreMotifSeed <- function(matrix, seed, weight = c()){
  result <- c()
  if(ncol(seed) > ncol(matrix)){
    #print('seed longer than matrix')
  }
  scores <- c()
  for(i in 1:(ncol(matrix)-ncol(seed)+1)){
    subMatrix <- normalize_sum1(matrix[,i:(i+ncol(seed)-1)], rows = F)
    scores <- c(scores, scoreMatrixSeed(subMatrix, seed, weight))
  }
  return(list(pos = which.min(scores), score = min(scores)))
}

#' Load From Index
#'
#' Load motif models in json format from a index form
#' @param index data frame of model information, contain to columns: gene_symbol and study
#' @param modelFile_Template dir of all model files, with specific identifier replaced as: $modelFile$
#' @return List of all motif models
#' @export
loadFromIndex <- function(index, modelFile_Template){
  len <- nrow(index)
  addList <- c()
  revList <- c()
  n <- 0
  for(i in 1:len){
    tryCatch({
      if(is.null(ncol(index))){
        name <- index[i]
      }else{
        name <- paste0(index$gene_symbol[i], '_',index$study[i])
      }
      modelFile <- gsub("\\$modelFile\\$", name, modelFile_Template)
      JSON_Lines <- readLines(modelFile)
      loop <- TRUE
      m <- 1
      while(loop){
        JSON_matrix <- tryCatch({
          JSON2Matrix(JSON_Lines, mode = m)
        }, error = function(x){return(c())}, warning = function(w){})
        if(length(JSON_matrix) == 0){
          loop <- FALSE
        }else{
          newEle <- list(gene_symbol = index$gene_symbol[i], study = index$study[i],
                         mode = m, matrix = JSON_matrix)
          addList[[length(addList)+1]] <- newEle
          revMotif <- JSON_matrix[nrow(JSON_matrix):1,ncol(JSON_matrix):1]
          rownames(revMotif) <- DNA()
          revEle <- list(gene_symbol = index$gene_symbol[i], study = index$study[i],
                         mode = -m, matrix = revMotif)
          addList[[length(addList)+1]] <- revEle
          m <- m+1
        }
      }
    },error = function(e){return(paste0('Error in ', i))})
  }
  return(addList)
}


#' Score Motif List
#'
#' Score a List of motifs
#' @param motifList List of motifs, results from loadFromIndex(), list of lists, each elements needs to have 'gene_symbol', 'study', 'mode', 'matrix'
#' @param seed seed matrix describing binding motif.
#' @param weight numeric vector of weights of each motif position of the seed
#' @return List of all motif models
#' @export
scoreMotifList <- function(motifList, seed, weight = c()){
  addList <- motifList
  motifsInfo <- data.frame(NULL)
  for(i in 1:length(addList)){
    tryCatch({
      score <- scoreMotifSeed(addList[[i]]$matrix, seed)
      addLine <- c(gene_symbol = addList[[i]]$gene_symbol, study = addList[[i]]$study,
                   mode = addList[[i]]$mode, score = score$score, pos = score$pos, id = i)
      motifsInfo <- rbind.data.frame(motifsInfo, addLine)
    },error = function(e){return(paste0('Error in ', i))})
  }
  colnames(motifsInfo) <- c('gene_symbol', 'study', 'mode', 'score', 'pos', 'id')
  return(motifsInfo)
}

#' Filter Motif List
#'
#' Filter from all motifs list for motifs according to a index data frame.
#' @param index filtered result for scoreMotifList()
#' @param motifList List of motifs, results from loadFromIndex(), used in scoreMotifList()
#' @param length length of motif alignment seed
#' @param posNames name of each motif position
#' @return List of filtered motif models (frequency measurements) according to index
#' @export
filterMotifList <- function(index, motifList, length = 8, posNames=paste0('P',1:length)){
  out <- list()
  for(i in 1:nrow(index)){
    mono <- motifList[[as.numeric(index$id[i])]]$matrix[,index$pos[i]:(as.numeric(index$pos[i])+length-1)]
    mono <- frequency2ddG(mono)
    mono <- ddG2frequency(mono)
    colnames(mono) <- posNames
    add <- list()
    add$name <- index$gene_symbol[i]
    add$matrix <- mono
    out[[i]] <- add
  }
  return(out)
}

#' Plot Prediction vs. True
#'
#' Make scatter plot between predicted and true value of motif prediction at a position
#' @param predTrue data frame with 2 columns true and pred, result from getBaseLineAccuracy(), closestSeqPred(), or SVDregression.Iterative.CV()
#' @param xlab x axis label
#' @param ylab y axis label
#' @param main plot title
#' @param xlim x axis limits
#' @param ylim y axis limits
#' @param throughZero If true the regression line passes through zero
#' @return Scatter plot showing prediction accuracy
#' @export
plotPredTrue <- function(predTrue, xlab = 'Experimental -ΔΔG/RT', ylab = 'Predicted -ΔΔG/RT',
                         main = 'prediction',xlim = c(1.05*min(predTrue$true),1.05*max(predTrue$true)),
                         ylim = c(1.05*min(predTrue$pred),1.05*max(predTrue$pred)),
                         throughZero = T){
  plot(predTrue$true,predTrue$pred, pch = 19, xlab = xlab, ylab = ylab, main = main, col = c('green','blue','orange','red'), xlim = xlim, ylim = ylim, cex = 1, cex.lab = 1.2, cex.axis = 1.2)
  if(throughZero){
    abline(lm(pred~true+0, predTrue), col = '#666666', lty = 2)
  }else{
    abline(lm(pred~true, predTrue), col = '#666666', lty = 2)
  }
  legend( x = "bottomright", legend = c('A', 'C', 'G', 'T'), col =c('green','blue','orange','red'), pch = 19, bty='n')
  legend('topleft', paste0('R^2 = ', round(groupedR2(predTrue$pred, predTrue$true, nrow(predTrue), throughZero=throughZero),4), '\n',
                           'RMSD = ', round(RMSD(predTrue$pred,predTrue$true),4)), bty = 'n')
}

#' make SVD model
#'
#' Make prediction of entire motif
#' @param mono_motifs Motif samples to incldue in the model, result from loadMono_motifs() or filterMotifList().
#' @param Alignment The alignment resulted from concatAli, ideally after matchAliMotif().
#' @param Positions positions to train the model for
#' @param Ftest_pVal threshold used for F-test
#' @return SVD regression model for entire motif
#' @export
makeSVDModel <- function(mono_motifs, Alignment,
                         Positions,
                         Ftest_pVal = 0.001){
  train_motifs <- mono_motifs
  train_alignment <- Alignment
  svdModel <- list()
  for(i in 1:length(Positions)){
    svd <- matrixSVD(gene2pos(train_motifs, pos = Positions[i]))
    svdModel[[Positions[i]]] <- trainSVD.Iterative(svd, train_alignment, Ftest_pVal = Ftest_pVal)
  }
  return(svdModel)
}

#' Predict PSAM with SVD model
#'
#' predict for multiple positions and output PSAM with pre-trained SVDmodel
#' @param SVDmodel pre-trained SVDmodel from makeSVDModel()
#' @param predSeq sequences to predict, matching alignment of train_alignment
#' @param pos positions to show in the PSAM
#' @param zero Lower limit to be assigned to predicted values
#' @return A PSAM in frequency measurements
#' @export
SVDmodelPredPSAM <- function(SVDmodel, predSeq, zero = 0.001){
  pos <- names(SVDmodel)
  test_alignment <- data.frame(name = 'Sample', alignment = predSeq)
  PSAM <- matrix(nrow = 4, ncol = length(pos), data = 1)
  rownames(PSAM) <- c('A', 'C', 'G','T')
  colnames(PSAM) <- pos
  for(p in 1:length(pos)){
    thisPos <- pos[p]
    #Predict binding motifs for test set
    pred_motifs <- predict(SVDmodel[[thisPos]], test_alignment, zero = zero)
    PSAM[,p] <- pred_motifs
  }
  return(PSAM)
}

#' predTrue ddG to frequency
#'
#' Change output of SVDregression.Iterative.CV() from ddG measurements to frequency measurements
#' @param predTrue predTrue table, output of SVDregression.Iterative.CV()
#' @param PFM If true, will output PFM (column sum = 1) instead of PSAM(column max = 1)
#' @return predTrue table
#' @export
predTrue.ddG2frequency <- function(predTrue, PFM = F){
  BaseLines <- predTrue
  colnames(BaseLines) <- c('true','pred')
  TrueMatrix <- matrix(ncol = nrow(BaseLines)/4, nrow = 4, data = BaseLines$true, byrow = F)
  PredMatrix <- matrix(ncol = nrow(BaseLines)/4, nrow = 4, data = BaseLines$pred, byrow = F)
  TrueMatrix <- ddG2frequency(TrueMatrix)
  PredMatrix <- ddG2frequency(PredMatrix)
  if(PFM){
    TrueMatrix <- apply(TrueMatrix, 2, function(x) x = x/sum(x))
    PredMatrix <- apply(PredMatrix, 2, function(x) x = x/sum(x))
  }
  out <- data.frame(true = as.numeric(TrueMatrix),pred = as.numeric(PredMatrix))
  return(out)
}

#' Similarity Regression Distance
#'
#' Calculate distance of similarity regression
#' @param seq1 Sequence 1
#' @param seq2 Sequence 2
#' @param weight pretrained weights for SR
#' @return sequence distance
#' @export
SRdist <- function(seq1, seq2, weight){
  ID <- c()
  for(i in 1:nchar(seq1)){
    if(substr(seq1,i,i) == substr(seq2,i,i)){
      ID <- c(ID, 1)
    }else{
      ID <- c(ID, 0)
    }
  }
  return(-sum(ID * weight))
}

#' Similarity Regression Prediction
#'
#' Predict motif model with Similarity regression
#' @param mono_motifs Motif samples to incldue in the model, result from loadMono_motifs() or filterMotifList().
#' @param Alignment The alignment resulted from concatAli, ideally after matchAliMotif().
#' @param pos positions to show in the PSAM
#' @param weightfile json file that contains pretrained weights from Lamber et. al.,
#' @return predTrue data frame
#' @export
SRpred <- function(mono_motifs, Alignment, pos, weightfile){
  weights <- RJSONIO::fromJSON(paste(readLines(weightfile), collapse = ''))
  weights <- weights$SR.Weights
  if(length(weights) != nchar(Alignment$alignment[1])){
    print('model and data not the same length')
  }

  pos_matrix <- gene2pos(mono_motifs, pos = pos)
  alignment <- Alignment
  trueList <- c()
  predList <- c()
  similarity <- c()
  for(i in 1:ncol(pos_matrix)){
    seq <- alignment$alignment[i]
    maxStr <- Inf
    predIndex <- 0
    for(j in 1:ncol(pos_matrix)){
      if(j != i){
        Str <- SRdist(seq, alignment$alignment[j], weights)
        if(Str < maxStr){
          maxStr <- Str
          predIndex <- j
        }
      }
    }
    similarity <- c(similarity, rep(maxStr,4))
    trueList <- c(trueList, unlist(data.frame(apply(log(as.data.frame(pos_matrix[,i])), 2, function(column) column - mean(column)))))
    predList <- c(predList, unlist(data.frame(apply(log(as.data.frame(pos_matrix[,predIndex])), 2, function(column) column - mean(column)))))
  }
  out <- data.frame(true = trueList, pred = predList, similarity = similarity)
  return(out)
}

#' make predictions with Similarity Regression
#'
#' Predict motif model with Similarity regression
#' @param predSeq sequences to predict, matching alignment of train_alignment
#' @param mono_motifs Motif samples to incldue in the model, result from loadMono_motifs() or filterMotifList().
#' @param Alignment The alignment resulted from concatAli, ideally after matchAliMotif().
#' @param weightfile json file that contains pretrained weights from Lamber et. al.,
#' @return predTrue data frame
#' @export
SRpredPSAM <- function(predSeq, mono_motifs, Alignment, weightfile){
  weights <- RJSONIO::fromJSON(paste(readLines(weightfile), collapse = ''))
  weights <- weights$SR.Weights
  if(length(weights) != nchar(Alignment$alignment[1])){
    print('model and data not the same length')
  }
  predList <- list()
  similarity <- c()
  source <- c()
  for(i in 1:nrow(predSeq)){
    seq <- predSeq$alignment[i]
    maxStr <- Inf
    predIndex <- 0
    for(j in 1:nrow(Alignment)){
      Str <- SRdist(seq, Alignment$alignment[j], weights)
      if(Str < maxStr){
        maxStr <- Str
        predIndex <- j
      }
    }
    similarity <- c(similarity, maxStr)
    source <- c(source, mono_motifs[[predIndex]]$name)
    add <- list()
    add$name <- predSeq$name[i]
    add$matrix <- mono_motifs[[predIndex]]$matrix
    predList[[i]] <- add
  }
  attr(predList, 'similarity') <- similarity
  attr(predList, 'source') <- source
  return(predList)
}

#' make predictions with Closest sequence
#'
#' Predict motif model with Closest sequence
#' @param predSeq sequences to predict, matching alignment of train_alignment
#' @param mono_motifs Motif samples to include in the model, result from loadMono_motifs() or filterMotifList().
#' @param Alignment The alignment resulted from concatAli, ideally after matchAliMotif().
#' @return predTrue data frame
#' @export
CSpredPSAM <- function(predSeq, mono_motifs, Alignment){
  predList <- list()
  similarity <- c()
  source <- c()
  for(i in 1:nrow(predSeq)){
    seq <- predSeq$alignment[i]
    maxStr <- Inf
    predIndex <- 0
    for(j in 1:nrow(Alignment)){
      Str <- stringdist::stringdist(seq, Alignment$alignment[j])
      if(Str < maxStr){
        maxStr <- Str
        predIndex <- j
      }
    }
    similarity <- c(similarity, maxStr)
    source <- c(source, mono_motifs[[predIndex]]$name)
    add <- list()
    add$name <- predSeq$name[i]
    add$matrix <- mono_motifs[[predIndex]]$matrix
    predList[[i]] <- add
  }
  attr(predList, 'similarity') <- similarity
  attr(predList, 'source') <- source
  return(predList)
}

#' make Amino Acid classification
#'
#' Classify AA according to biophysical properties
#' @param AAs A vector of amino acid in 1 letter code.
#' @return Vector of classification with attribute 'encode' encoding proprieties
#' @export
AAclassify <- function(AAs){
  props <- c()
  encodes <- c()
  for(i in 1:length(AAs)){
    AA <- AAs[i]
    H <- AAfeatures()[AA,'Hydropathy']
    E <- AAfeatures()[AA,'PI']
    if(AA == 'F'){
      prop <- 'aromatic'
      encode <- -100
    }else if(AA == 'Y'){
      prop <- 'aromatic'
      encode <- 100
    }else if(H > -0.6){
      prop <- 'non-polar'
      encode <- -10
    }else if(E > 7.6){
      prop <- 'positive'
      encode <- 1
    }else if(E < 4) {
      prop <- 'negative'
      encode <- -1
    }else{
      prop <- 'polar'
      encode <- 0
    }
    props <- c(props, prop)
    encodes <- c(encodes, encode)
  }
  attr(props, 'encode') <- encodes
  return(props)
}

#' Classify AA-AA interaction
#'
#' Classify AA-AA interaction according to AA propaties.
#' @param prop1s A vector of amino acid proprieties, atrribute 'encode' of AAclassify() output.
#' @param prop2s A vector of amino acid proprieties, atrribute 'encode' of AAclassify() output.
#' @return Vector of interaction classification
#' @export
AAinteraction <- function(prop1s, prop2s){
  outList <- c()
  for(i in 1:length(prop1s)){
    prop1 <- prop1s[i]
    prop2 <- prop2s[i]
    if(prop1 == prop2 && prop1 == 0){
      interaction <- 'h-bond'
    }else if (abs(prop1) + abs(prop2) == 200){
      interaction <- 'Pi-stacking'
    }else if ((prop1 + prop2) == -20 || (prop1 + prop2) == -110){
      interaction <- 'hydrophobic'
    }else if ((prop1 + prop2) == 0){
      interaction <- 'electro-static'
    }else if (abs(prop1 + prop2) < 2){
      interaction <- 'h-bond'
    }else if (abs(prop1 + prop2) < 3){
      interaction <- 'repel'
    }else if (prop1 + prop2 > 90){
      interaction <- 'h-bond'
    }else{
      interaction <- 'repel'
    }
    outList <- c(outList, interaction)
  }
  return(outList)
}

#' load motif from json
#'
#' Load mono motif from JSON
#' @param modelFile json file with consensus motif
#' @param gene gene symbol of the motif
#' @param study study of the motif
#' @return list of motifs in the json file forward and backward
#' @export
loadJSONmotif <- function(modelFile, gene = 'sampleGene', study = 'NA'){
  JSON_Lines <- readLines(modelFile)
  loop <- TRUE
  m <- 1
  addList <- list()
  while (loop) {
    JSON_matrix <- tryCatch({
      JSON2Matrix(JSON_Lines, mode = m)
    }, error = function(x) {
      return(c())
    }, warning = function(w) {
    })
    if (length(JSON_matrix) == 0) {
      loop <- FALSE
    }
    else {
      newEle <- list(gene_symbol = gene,
                     study = study, mode = m, matrix = JSON_matrix)
      addList[[length(addList) + 1]] <- newEle
      revMotif <- JSON_matrix[nrow(JSON_matrix):1,
                              ncol(JSON_matrix):1]
      rownames(revMotif) <- DNA()
      revEle <- list(gene_symbol = gene,
                     study = study, mode = -m, matrix = revMotif)
      addList[[length(addList) + 1]] <- revEle
      m <- m + 1
    }
  }
  return(addList)
}
