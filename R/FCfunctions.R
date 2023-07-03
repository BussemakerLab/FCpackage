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

#SNAP_gene related
bpAA <- function(SNAPList){
  bp <- substr(SNAPList$bp.aa,1,2)
  AA <- substr(SNAPList$bp.aa,4,6)
  bpAA <- data.frame(bp, AA)
  return(bpAA)
}

#SNAP_gene related
extractPair.OriginAA <- function(
  SNAPList.AA,
  pml,
  SNAPdir,
  inspect.file,
  dssr = 'E:/x3dna-dssr.exe',
  cutoff = 4.5
){
  if(!file.exists(inspect.file)){
    dir.create(inspect.file)
  }
  Files <- c()
  for(i in 1:nrow(SNAPList.AA)){
    pairFile <- extractPair(testSNAP = SNAPList.AA[i,], pml = pml, SNAPFile = SNAPdir,
                            pdbFile = inspect.file, tag = paste('_',SNAPList.AA$aa[i], sep = ''),
                            cutoff = cutoff)
    AAcode <- paste(substr(SNAPList.AA$aa[i],1,2), substr(SNAPList.AA$aa[i], 6, nchar(SNAPList.AA$aa[i])), sep = '')
    cmd <- paste(dssr, ' -i=', pairFile,' --frame-aa=', AAcode,' -o=',pairFile, sep = '')
    system(cmd)
    Files <- c(Files, pairFile)
  }
  return(Files)
}

#SNAP_gene related
screenSNAP.AA <- function(
  PDBdir,
  SNAPdir,
  pdbid ,
  hmmIndex,
  AApos,
  SNAP.update = TRUE,
  block = "pair/amino-acid interactions",
  cutoff = 4.5
){
  out <- data.frame(NULL)
  for(j in 1:length(pdbid)){
    snapName <- list.files(SNAPdir)[grep(pdbid[j], list.files(SNAPdir))]
    if(length(snapName) != 1 || SNAP.update){
      snap <- formSNAP(PDBdir, SNAPdir, pdbid[j], dssr, cutoff = cutoff)
    }else{
      snap <- paste(SNAPdir, '/', snapName, sep = '')
    }
    snap_df <- readSNAP(snap, block)
    indexPDB <- hmmIndex[hmmIndex$pdbid == pdbid[j],]
    for(i in 1:nrow(indexPDB)){
      chain <- indexPDB$chain[i]
      pos <- indexPDB[i,colnames(indexPDB) == paste('AA', AApos, sep = '')]
      add <- snap_df[substr(snap_df$aa,1,1)==chain,]
      add <- add[substr(add$aa,6,nchar(add$aa))==pos,]
      out <- rbind.data.frame(out,add)
    }
  }
  return(out)
}

#SNAP_gene related
#' Read HMM file
#'
#' Read in HMM file of alignment data and build data table
#' @param hmmFile Path to the HMM file to be parsed
#' @return Data frame containing PDBid, chain, start, end, and alignment
#' @export
readHmm <- function(hmmFile){
  hmmLines <- read.csv(hmmFile, header=FALSE, sep=";")
  hmmAlign <- read.csv(hmmFile,header = FALSE, quote="\"", comment = '#', sep = '')
  AlLen <- nchar(hmmAlign[1,2])
  #organize hmmLines
  GS <- substr(hmmLines[,1],1,4) == '#=GS'
  hmmLines <- hmmLines[GS,]

  hmmIndex <- data.frame(NULL)
  for(i in 1:length(hmmLines)){
    linetest <- hmmLines[i]
    split1 <- strsplit(linetest, '\\[')
    split2 <- strsplit(split1[[1]], '\\]')
    id <- split2[[1]][1]
    id <- strsplit(id, ' ')[[1]]
    id <- id[grep('/', id)]
    for(j in 4:length(split2)){
      addline <- c()
      entry <- split2[[j]][1]
      entry <- strsplit(entry, ', ')[[1]]
      pos <- strsplit(entry[3], '-')[[1]]
      addline <- c(addline, entry[1:2], pos, id)
      ind <- as.numeric(pos[1])
      seq <- hmmAlign[hmmAlign$V1==id,2]
      for(k in 1:AlLen){
        AA <- substr(seq, k,k)
        if(AA == '.' || AA == '-'){
          addline <- c(addline, NA)
        }else{
          addline <- c(addline, ind)
          ind <- ind+1
        }
      }
      hmmIndex <- rbind.data.frame(hmmIndex, addline)
    }
  }
  colnames(hmmIndex) <- c('pdbid', 'chain', 'start', 'end', 'alignment', paste('AA', c(1:AlLen), sep = ''))
  return(hmmIndex)
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

#SNAP_gene related
AnchorOriginSNAP <- function(PDBdir,
                             SNAPdir,
                             anchorID,
                             line,
                             block = 'pair/amino-acid interactions',
                             pml,
                             inspect.file,
                             cutoff = 4.5
                             ){
  Anchor <- selectAnchor(PDBdir, SNAPdir, pdbid = anchorID, block = block, line = line, cutoff = cutoff)
  Adt <- rbind.data.frame(Anchor,Anchor)
  colnames(Adt) <- names(Anchor)
  Adt <- Adt[1,]
  ApdbFile <- extractPair(Adt, pml, SNAPFile = SNAPdir, pdbFile = inspect.file, block = block, cutoff = cutoff)
  Anchor.pdb <- bio3d::read.pdb(ApdbFile)
  Anchor.origin <- Anchor.pdb$atom[Anchor.pdb$atom$elety == 'CA',c('x','y','z')]
  Anchor.new <- c(Anchor, Anchor.origin)
  attr(Anchor.new, 'block') <- attr(Anchor,'block')
  attr(Anchor.new, 'cutoff') <- attr(Anchor, 'cutoff')
  Anchor <- Anchor.new
  return(Anchor)
}

#SNAP_gene related
ScreenOriginSNAP <- function(PDBdir, SNAPdir, pdbid, pml, Anchor, inspect.file, structure = F,
                             cutoff = 4.5,
                             threshold = c(5,50),
                             tag = '',
                             block = 'pair/amino-acid interactions',
                             pymol.dir = 'E:/pymol/pymol_app/pyMOLWin.exe',
                             wait.time = 3,
                             SNAP.update = T
                             ){
    if(SNAP.update){
      snap.file <- formSNAP(PDBdir = PDBdir, SNAPdir = SNAPdir, pdbid = pdbid, auxfile = T, cutoff = cutoff)
    }else{
      snap.file <- paste()
    }

    snap <- readSNAP(snap.file, block = attr(Anchor,'block'))
    pairs <- unique(snap$bp.aa)
    an.ch <- unlist(Anchor[c('Cx','Cy','Cz','Rx','Ry','Rz')])
    if(length(unique(pairs)) < 1){
      return()
    }
    snap.xyz <- data.frame(NULL)
    for(i in 1:length(unique(pairs))){
      snap.pair <- snap[snap$bp.aa == pairs[i],]
      for(j in 1:nrow(snap.pair)){
        sa.ch <- snap.pair[j,c('Cx','Cy','Cz','Rx','Ry','Rz')]
        distOR <- CAdist(as.numeric(an.ch),as.numeric(sa.ch))
        if(distOR[1] <= threshold[1] && distOR[2] <= threshold[2]){
          add <- c(unlist(snap.pair[j,]),distO = distOR[1], distR = distOR[2])
          snap.xyz <- rbind.data.frame(snap.xyz, add)
          pdb.file <- paste(PDBdir, '/',pairs[i], '.pdb', sep = '')
          colnames(snap.xyz) <- names(add)
          if(structure){
            cat('#save state \n', file = pml, append = F)
            cat('load ', pdb.file,' \n', sep = '', file = pml, append = T)
            outfile <- paste(inspect.file, '/', snap.pair$id[1], '_', pairs[i], '_', j,'_', tag,'_snap', '.pdb', sep = '')
            cat('save ', outfile,', state=', j, ' \n', sep = '', file = pml, append = T)
            cat('quit \n', file = pml, append = T)
            cmd <- paste(pymol.dir, ' -qc ', pml, sep = '')
            system(cmd, wait = TRUE)
            Sys.sleep(wait.time)
          }
        }
      }
    }
    options(warn=-1)
    files <- list.files(PDBdir)
    deleFiles <- files[!complete.cases(as.numeric(substr(files,1,1)))]
    file.remove(paste(PDBdir, '/', deleFiles, sep = ''))
    options(warn=0)
  return(snap.xyz)
}

#SNAP_gene related
CAdist <- function(anchor, sample){
  distO <- dist(rbind(anchor[1:3], sample[1:3]))
  distR <- (rotationDiff(anchor[4], sample[4])^2 + rotationDiff(anchor[5], sample[5])^2 +
              rotationDiff(anchor[6], sample[6])^2)^(1/2)
  return(c(distO,distR))
}

#SNAP_gene related
rotationDiff <- function(R1, R2){
  P1 <- R1-R2
  if(P1 > 180){
    P1 <- -(360-P1)
  }else if(P1 < -180){
    P1 <- 360+P1
  }
  return(P1)
}


#SNAP_gene related
TR2matrix <- function(TR, R.first = T){ #colnames have to be Tx, Ty, Tz, Rx, Ry, Rz
  out <- matrix(nrow = 4, ncol = 3, data = c(0,0,0,1,0,0,0,1,0,0,0,1), byrow = T)
  colnames(out) <- c('x', 'y', 'z')
  rownames(out) <- c('origin', 'Xaxis', 'Yaxis', 'Zaxis')
  cols <- names(TR)
  xTh <- as.numeric(TR[4])/180*pi
  yTh <- as.numeric(TR[5])/180*pi
  zTh <- as.numeric(TR[6])/180*pi
  Tx <- as.numeric(TR[1])
  Ty <- as.numeric(TR[2])
  Tz <- as.numeric(TR[3])
  RxM <- matrix(nrow = 3, ncol = 3, data = c(1,0,0,0,cos(xTh), -sin(xTh), 0, sin(xTh), cos(xTh)), byrow = T)
  RyM <- matrix(nrow = 3, ncol = 3, data = c(cos(yTh),0,sin(yTh),0,1, 0, -sin(yTh), 0, cos(yTh)), byrow = T)
  RzM <- matrix(nrow = 3, ncol = 3, data = c(cos(zTh), -sin(zTh), 0,sin(zTh), cos(zTh), 0,0,0,1 ), byrow = T)
  if(R.first){
    for(i in 2:nrow(out)){
      out[i,] <- RzM%*%RyM%*%RxM%*%out[i,]
    }
    out[1,] <- out[1,] + Tx*out[2,] + Ty*out[3,] + Tz*out[4,]
  }else{
    out[1,] <- out[1,] + Tx*out[2,] + Ty*out[3,] + Tz*out[4,]
    for(i in 2:nrow(out)){
      out[i,] <- RzM%*%RyM%*%RxM%*%out[i,]
    }
  }
  return(out)
}

#SNAP_gene related
angle <- function(x,y){
  dot.prod <- x%*%y
  norm.x <- norm(x,type="2")
  norm.y <- norm(y,type="2")
  theta <- acos(dot.prod / (norm.x * norm.y))
  as.numeric(theta)
}

#SNAP_gene related
TRmatrixDist <- function(TRM1, TRM2){
  OriginDist <- dist(rbind.data.frame(TRM1[1,], TRM2[1,]))
  A1 <- angle(TRM1[2,], TRM2[2,])
  A2 <- angle(TRM1[3,], TRM2[3,])
  A3 <- angle(TRM1[4,], TRM2[4,])
  AngleDist <- (A1^2 + A2^2 + A3^2)^(1/2)
  return(c(OriginDist,AngleDist))
}

#SNAP_gene related
extractPair <- function(testSNAP,
                        pml,
                        SNAPFile,
                        pdbFile,
                        block = 'pair/amino-acid interactions',
                        pymol.dir = 'E:/pymol/pymol_app/pyMOLWin.exe',
                        wait.time = 3,
                        cutoff = 4.5,
                        tag = ''){
  out <- c()
  for(i in 1:nrow(testSNAP)){
    snap <- formSNAP(PDBdir = pdbFile,
                     SNAPdir = SNAPFile,
                     pdbid = testSNAP$id[i], auxfile = T, cutoff = cutoff)
    snap_df <- as.data.frame(readSNAP(snap,block))
    pair <- testSNAP$bp.aa[i]
    snap_df <- snap_df[snap_df$bp.aa == pair,]
    state <- intersect(which(as.numeric(snap_df$Tdst) == as.numeric(testSNAP$Tdst[i])), which(as.numeric(snap_df$Rdst) == as.numeric(testSNAP$Rdst[i])))
    pairFile <- paste(pair, '.pdb', sep = '')
    cat('#save state \n', file = pml, append = F)
    cat('load ', pdbFile, '/',pairFile,' \n', sep = '', file = pml, append = T)
    outfile <- paste(pdbFile, '/', testSNAP$id[i], '_', pair, '_', state, tag, '.pdb', sep = '')
    cat('save ', pdbFile, '/', testSNAP$id[i], '_', pair, '_', state, tag, '.pdb, state=', state, ' \n', sep = '', file = pml, append = T)
    cat('quit \n', file = pml, append = T)
    cmd <- paste(pymol.dir, ' -qc ', pml, sep = '')
    system(cmd, wait = TRUE)
    Sys.sleep(wait.time)
    snap <- formSNAP(PDBdir = pdbFile,
                     SNAPdir = SNAPFile,
                     pdbid = testSNAP$id[i], auxfile = F, cutoff = cutoff)
    out <- c(out, outfile)
  }
  options(warn=-1)
  files <- list.files(pdbFile)
  deleFiles <- files[!complete.cases(as.numeric(substr(files,1,1)))]
  file.remove(paste(pdbFile, '/', deleFiles, sep = ''))
  options(warn=0)
  return(out)
}


#SNAP_gene related
addBaseLink <- function(testScreen, PDBdir){
  baseLink <- c()
  for(i in 1:nrow(testScreen)){
    x <- tryCatch({
      a1 <- readPosition(bio3d::read.pdb(download.pdb(testScreen$pdb[i], PDBdir,F)), baseChain = testScreen$baseChain[i], baseResi = testScreen$baseNO[i],
                       AAChain = testScreen$AAChain[i], AAResi = testScreen$AANo[i],
                       baseAnchor = c('N1'), AAAnchor = c('CA'))

      a2 <- readPosition(bio3d::read.pdb(download.pdb(testScreen$pdb[i], PDBdir,F)), baseChain = testScreen$base2Chain[i], baseResi = testScreen$base2NO[i],
                         AAChain = testScreen$AAChain[i], AAResi = testScreen$AANo[i],
                         baseAnchor = c('N1'), AAAnchor = c('CA'))
    }, error = function(e) e)
    if(length(x) == 2){
      baseLink <- c(baseLink,-1)
    }else if(a1 < a2){
      baseLink <- c(baseLink, testScreen$baseType[i])
    }else{
      baseLink <- c(baseLink, testScreen$base2Type[i])
    }
  }
  return(cbind.data.frame(testScreen, baseLink))
}

#SNAP_gene related
getDNAPos <- function(string){
  posList <- c()
  posList <- c(posList, gregexpr('A', string)[[1]])
  posList <- c(posList, gregexpr('T', string)[[1]])
  posList <- c(posList, gregexpr('C', string)[[1]])
  posList <- c(posList, gregexpr('G', string)[[1]])
  if(length(posList) > 0){
    return(max(posList))
  }else{
    return(-1)
  }

}

#SNAP_gene related
SNAP2screenList <- function(SNAPList){
  out <- data.frame(NULL)
  if(length(grep('nt2', colnames(SNAPList))) > 0){
    for(i in 1:nrow(SNAPList)){
      pdb <- SNAPList$id[i]
      base1 <- strsplit(SNAPList$nt1[i],'\\.')[[1]]
      baseChain <- base1[1]
      pos <- getDNAPos(base1[2])
      baseType <- substr(base1[2], 1, pos)
      baseNO <- substr(base1[2], pos+1, nchar(base1[2]))
      base1 <- strsplit(SNAPList$nt2[i],'\\.')[[1]]
      base2Chain <- base1[1]
      pos <- getDNAPos(base1[2])
      base2Type <- substr(base1[2], 1, pos)
      base2NO <- substr(base1[2], pos+1, nchar(base1[2]))
      base1 <- strsplit(SNAPList$aa[i],'\\.')[[1]]
      AAChain <- base1[1]
      AAType <- gsub('[[:digit:]]', '', base1[2])
      AANo <- gsub(AAType, '', base1[2])

      add <- data.frame(pdb, baseChain, baseNO, baseType, base2Chain, base2NO, base2Type,
                        AAChain, AANo, AAType, maxDiff = as.numeric(SNAPList$dist[i]))
      out <- rbind.data.frame(out, add)
    }
  }else{
    for(i in 1:nrow(SNAPList)){
      pdb <- SNAPList$id[i]
      base1 <- strsplit(SNAPList$nt[i],'\\.')[[1]]
      baseChain <- base1[1]
      baseType <- gsub('[[:digit:]]', '', base1[2])
      baseNO <- gsub(baseType, '', base1[2])
      base1 <- strsplit(SNAPList$aa[i],'\\.')[[1]]
      AAChain <- base1[1]
      AAType <- gsub('[[:digit:]]', '', base1[2])
      AANo <- gsub(AAType, '', base1[2])

      add <- data.frame(pdb, baseChain, baseNO, baseType,
                        AAChain, AANo, AAType, maxDiff = as.numeric(SNAPList$dist[i]))
      out <- rbind.data.frame(out, add)
    }
  }
  return(out)
}


#SNAP_gene related
pymolInspect <- function(screenList, projectFile, outFile, scriptFile, outTag = '', base.no = 1, adj = F){
  cat('#screenList for pyMol \n', file = scriptFile, append = F)
  for(i in 1:nrow(screenList)){
    pdbname <- screenList$pdb[i]
    cat('load ', projectFile, '/',pdbname,'.pdb \n', sep = '', file = scriptFile, append = T)
    if(base.no == 1){
      cat('select (chain ',screenList$baseChain[i], ' and resi ', screenList$baseNO[i], ') or (chain ',
          screenList$AAChain[i], ' and resi ', screenList$AANo[i], ') \n', sep = '', file = scriptFile, append = T)
    }else if(base.no == 2){
      if(adj){
        cat('select (chain ',screenList$baseChain[i], ' and resi ', as.numeric(screenList$baseNO[i]), '-',
            as.numeric(screenList$baseNO[i])+1, ') or (chain ',
            screenList$AAChain[i], ' and resi ', screenList$AANo[i], ') or (chain ', screenList$base2Chain[i], ' and resi ',
            as.numeric(screenList$base2NO[i]), '-', as.numeric(screenList$base2NO[i])+1,
             ') \n', sep = '', file = scriptFile, append = T)
      }else{
        cat('select (chain ',screenList$baseChain[i], ' and resi ', screenList$baseNO[i], ') or (chain ',
            screenList$AAChain[i], ' and resi ', screenList$AANo[i], ') or (chain ', screenList$base2Chain[i], ' and resi ',
            screenList$base2NO[i], ') \n', sep = '', file = scriptFile, append = T)
      }
    }else{
      print('error in base.no')
      return()
    }

    cat('create obj01, sele \n', sep = '', file = scriptFile, append = T)
    cat('delete ', pdbname, ' \n', sep = '', file = scriptFile, append = T)
    cat('show sticks, obj01 \n', file = scriptFile, append = T)
    cat('save ', outFile, '/', pdbname, '_', outTag, '.pdb \n', sep = '', file = scriptFile, append = T)
    cat('delete all \n', file = scriptFile, append = T)
  }
  cat('quit \n', file = scriptFile, append = T)
}


#SNAP_gene related
screenPosition <- function(pdbList,
                           projectFile,
                           threshold = 0.5,
                           anchorMatrix){
  out <- data.frame(NULL)
  for(i in 1:length(pdbList)){
    tryCatch({
      pdbname <- pdbList[i]
      pdb <- download.pdb(pdbname, projectFile = projectFile, as.file = FALSE)
      pdb <- read.pdb(pdb)
      atom <- pdb$atom
      baseAtoms <- data.frame(NULL)
      AAAtoms <- data.frame(NULL)
      for(j in 1:nrow(atom)){
        Line <- atom[j,]
        if(Line$resid == 'DA' || Line$resid == 'DC' || Line$resid == 'DG'|| Line$resid == 'DT'){
          baseAtoms <- rbind.data.frame(baseAtoms, Line)
        }else if(Line$resid == 'HOH'){
          next()
        }else{
          AAAtoms <- rbind.data.frame(AAAtoms, Line)
        }
      }
      baseUse <- data.frame(NULL)
      for(k in 1:nrow(anchorMatrix)){
        baseAdd <- baseAtoms[grep(rownames(anchorMatrix)[k], baseAtoms$elety),]
        baseAdd <- baseAdd[nchar(baseAdd$elety) == nchar(rownames(anchorMatrix)[k]),]
        baseUse <- rbind.data.frame(baseUse,baseAdd)
      }
      basePairs <- unique(baseUse[,c('chain','resno')])
      AAUse <- data.frame(NULL)
      for(k in 1:ncol(anchorMatrix)){
        AAAdd <- AAAtoms[grep(colnames(anchorMatrix)[k], AAAtoms$elety),]
        AAAdd <- AAAdd[nchar(AAAdd$elety) == nchar(colnames(anchorMatrix)[k]),]
        AAUse <- rbind.data.frame(AAUse,AAAdd)
      }
      AAPairs <- unique(AAUse[,c('chain','resno')])
      for(k in 1:nrow(basePairs)){
        baseNow <- baseUse[baseUse[,'chain'] == basePairs[k,'chain'],]
        baseNow <- baseNow[baseNow[,'resno'] == basePairs[k,'resno'],]
        if(nrow(baseNow) < nrow(anchorMatrix)){
          next()
        }
        for(l in 1:nrow(AAPairs)){
          AANow <- AAUse[AAUse[,'chain'] == AAPairs[l,'chain'],]
          AANow <- AANow[AANow[,'resno'] == AAPairs[l,'resno'],]
          if(nrow(AANow) < ncol(anchorMatrix)){
            next()
          }
          fit <- 0
          maxDiff <- -Inf
          for(A1 in 1:nrow(anchorMatrix)){
            for(A2 in 1:ncol(anchorMatrix)){
              true <- anchorMatrix[A1,A2]
              XYZbase <- baseNow[grep(rownames(anchorMatrix)[A1],baseNow$elety),c('x','y','z')]
              XYZAA <- AANow[grep(colnames(anchorMatrix)[A2],AANow$elety),c('x','y','z')]
              if(nrow(XYZbase) != 1 || nrow(XYZAA) != 1){
                next
              }
              test <- dist(rbind(XYZbase,XYZAA))
              if(abs(true - test) < threshold){
                fit <- fit + 1
                if(abs(true-test) > maxDiff){
                  maxDiff <- abs(true-test)
                }
              }
            }
          }
          if(fit == nrow(anchorMatrix)*ncol(anchorMatrix)){
            add <- data.frame(pdb = pdbList[i], baseChain = baseNow$chain[1], baseNO = baseNow$resno[1], baseType = baseNow$resid[1],
                              AAChain = AANow$chain[1], AANo = AANow$resno[1], AAType = AANow$resid[1], maxDiff = as.numeric(maxDiff))
            out <- rbind.data.frame(out, add)
            print('Base-AA pair found')
          }
        }
      }
      #print(paste('screened over pdb: ', pdbname, '; pdb number:', i, sep = ''))
    }, error = function(e) e)
  }
  out <- unique(out)
  return(out)
}


#SNAP_gene related
readPosition <- function(pdb,
                         baseChain = 'B',
                         baseResi = 9,
                         AAChain = 'A',
                         AAResi = 217,
                         baseAnchor = c('P','N9','C1','O4'),
                         AAAnchor = c('CA','CZ')){
  atom <- pdb$atom
  base <- atom[atom$chain == baseChain,]
  base <- base[base$resno == baseResi,]
  AA <- atom[atom$chain == AAChain,]
  AA <- AA[AA$resno == AAResi,]
  outMatrix <- matrix(nrow = length(baseAnchor), ncol = length(AAAnchor), data = NA)
  rownames(outMatrix) <- baseAnchor
  colnames(outMatrix) <- AAAnchor
  for(i in 1:length(baseAnchor)){
    baseAtom <- base[grep(baseAnchor[i], base$elety),]
    baseAtom <- baseAtom[nchar(baseAtom$elety) == min(nchar(baseAtom$elety)),]
    baseXYZ <- c(baseAtom$x,baseAtom$y,baseAtom$z)
    for(j in 1:length(AAAnchor)){
      AAAtom <- AA[grep(AAAnchor[j], AA$elety),]
      AAAtom <- AAAtom[nchar(AAAtom$elety) == min(nchar(AAAtom$elety)),]
      AAXYZ <- c(AAAtom$x,AAAtom$y,AAAtom$z)
      addDist <- dist(rbind(baseXYZ,AAXYZ))
      outMatrix[i,j] <- addDist
    }
  }
  return(outMatrix)
}

#SNAP_gene related
formSNAP <- function(PDBdir, SNAPdir = PDBdir, pdbid, dssr = 'E:/x3dna-dssr.exe', auxfile = F, cutoff = 4.5){
  wd <- getwd()
  setwd(PDBdir)
  proteinStructureBoost::download.pdb(pdbid, projectFile = PDBdir,as.file = FALSE)
  pdbfile <- paste(pdbid,'.pdb',sep = '')
  SNAPfile <- paste(SNAPdir, '/', pdbid, '_SNAP.txt', sep = '')
  if(auxfile){
    cmd <- paste(dssr, ' snap -i=', pdbfile, ' --c-alpha --pair-order=as-is --cutoff=',cutoff,' --auxfile > ', SNAPfile, sep = '')
  }else{
    cmd <- paste(dssr, ' snap -i=', pdbfile, ' --c-alpha --pair-order=as-is --cutoff=',cutoff,' > ', SNAPfile, sep = '')
  }
  shell(cmd)
  setwd(wd)
  return(SNAPfile)
}

#SNAP_gene related
readSNAP <- function(SNAPfile, block = 'nucleotide/amino-acid interaction'){
  options(warn=-1)
  data <- read.table(file = SNAPfile,sep='\n',fill=TRUE,header = F)
  start <- grep(block, data$V1) + 2
  if(length(start) != 1){
    print('error in block value')
    return()
  }
  end <- min(grep('\\*', data$V1)[grep('\\*', data$V1) > start]) - 1
  if(length(end) == 0){
    end <- nrow(data)
  }
  trim <- strsplit(gsub('\\s+', '\t', data[start:end,]), split = '\t')
  if(length(trim) > 999){
    for(tt in 1000:length(trim)){
      trim[[tt]] <- c('',trim[[tt]])
    }
  }
  trim <- as.data.frame(trim)
  trim <- t(trim)
  trim <- trim[,3:ncol(trim)]
  colnames(trim) <- strsplit(gsub('\\s+', '\t', data[start-1,]), split = '\t')[[1]][-1]
  colnames(trim)[2] <- 'bp.aa'
  rownames(trim) <- 1:nrow(trim)
  options(warn=0)
  return(as.data.frame(trim))
}

#SNAP_gene related
selectAnchor <- function(PDBdir, SNAPdir = PDBdir, pdbid, dssr = 'E:/x3dna-dssr.exe', line = NA, block = 'nucleotide/amino-acid interaction', cutoff = 4.5){
  snap <- formSNAP(PDBdir, SNAPdir, pdbid, dssr, cutoff = cutoff)
  snap_df <- readSNAP(snap, block)
  if(is.na(line)){
    print(snap_df)
    anchorLine <- my.name <- readline(prompt="Enter anchor interaction line number: ")
  }else{
    anchorLine <- line
  }
  anchorLine <- snap_df[anchorLine,]
  attr(anchorLine, 'block') <- block
  attr(anchorLine, 'cutoff') <- cutoff
  return(anchorLine)
}

#SNAP_gene related
screenSNAP <- function(PDBdir, SNAPdir = PDBdir, anchorLine, pdbid,  dssr = 'E:/x3dna-dssr.exe',
                       threshold = c(5,1), method = 'matrix'){ #avaliable methods: raw (8 thresh), abs (8 thresh), matrix (2 thresh)
  block <- attr(anchorLine, 'block')
  cutoff <- attr(anchorLine, 'cutoff')
  out <- data.frame(NULL)
  snapName <- list.files(SNAPdir)[grep(pdbid, list.files(SNAPdir))]
  if(length(snapName) != 1){
    snap <- formSNAP(PDBdir, SNAPdir, pdbid, dssr, cutoff = cutoff)
  }else{
    snap <- paste(SNAPdir, '/', snapName, sep = '')
  }
  snap_df <- readSNAP(snap, block)
  if(is.null(snap_df)){
    attr(out, 'is.dep') <- 1
    return(out)
  }
  anchorCoord <- anchorLine[c('Tdst', 'Rdst', 'Tx', 'Ty', 'Tz', 'Rx', 'Ry', 'Rz')]
  snapCoords <- snap_df[,c('Tdst', 'Rdst', 'Tx', 'Ty', 'Tz', 'Rx', 'Ry', 'Rz')]
  for(i in 1:nrow(snap_df)){
    if(method == 'raw'){
      dist <- SNAP8dDist.abs(anchorCoord, snapCoord = snapCoords[i,], threshold = threshold)
    }else if(method == 'abs'){
      dist <- SNAP8dDist(anchorCoord, snapCoord = snapCoords[i,], threshold = threshold)
    }else if(method == 'matrix'){
      dist <- SNAPMatrixDist(anchorCoord, snapCoord = snapCoords[i,], threshold = threshold)
      if(as.numeric(anchorCoord[2])*as.numeric(snapCoords[i,2]) < 0){
        dist <- Inf
      }
    }
    if(dist <= 1){
      dist <- as.numeric(dist)
      line <- c(snap_df[i,], dist = dist)
      line <- t(as.data.frame(line))
      out <- rbind.data.frame(out, line)
    }
  }
  return(out)
}

#SNAP_gene related
SNAP8dDist <- function(anchorCoord, snapCoord, threshold = c(Inf, Inf, 1,1,1,Inf,Inf,Inf)){
  return(sum(((as.numeric(anchorCoord) - as.numeric(snapCoord))/threshold)^2)^(1/2))
}

#SNAP_gene related
SNAP8dDist.abs <- function(anchorCoord, snapCoord, threshold = c(Inf, Inf, 1,1,1,Inf,Inf,Inf)){
  return(sum(((abs(as.numeric(anchorCoord)) - abs(as.numeric(snapCoord)))/threshold)^2)^(1/2))
}

#SNAP_gene related
SNAPMatrixDist <- function(anchorCoord, snapCoord, threshold = c(5,0.5)){
  TRM1 <- TR2matrix(anchorCoord[3:8])
  TRM2 <- TR2matrix(snapCoord[3:8])
  dist <- TRmatrixDist(TRM1, TRM2)
  if(dist[1] <= threshold[1] && dist[2] <= threshold[2]){
    return(dist[1]/threshold[1])
  }else{
    return(Inf)
  }
}

#SNAP_gene related
SNAPsummary <- function(data){
  SNAP.summary <- matrix(nrow = 3, ncol = 8, data = 0)
  rownames(SNAP.summary) <- c('mean', 'median', 'sd')
  colnames(SNAP.summary) <- c('Tdst', 'Rdst', 'Tx', 'Ty', 'Tz', 'Rx', 'Ry', 'Rz')
  for(i in 1:length(colnames(SNAP.summary))){
    SNAP.summary[1,i] <- mean(as.numeric(data[,colnames(SNAP.summary)[i]]))
    SNAP.summary[2,i] <- median(as.numeric(data[,colnames(SNAP.summary)[i]]))
    SNAP.summary[3,i] <- sd(as.numeric(data[,colnames(SNAP.summary)[i]]))
  }
  return(SNAP.summary)
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
#' @param bHLH_index data frame containing the identifiers of each experiment,
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
#' list mono_motifs.
#' @return List of motifs in mono_motifs, also an info table if withTable == T
#' @export
loadMono_motifs <- function(bHLH_index,
                            modelFile_Template = '~/$modelFile$/result/fit.models.consensus.json',
                            rec_seq = 'CANNTG',
                            pos_index = c('P-3','P-2','P-1','P1','P2','P3'),
                            checkSymmetry = FALSE,
                            withTable = TRUE,
                            useMode = c()){
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
#' @param posMotif Position in motif to extract
#' @param randomSample If true, randomly select from studies of the sample protein; if fasle, takes the first two studies
#' @return Data frame with two columns with pair-wise ddG values that can be plotted
#' @export
getBaseLineAccuracy <- function(mono_motifs, posMotif, randomSample = FALSE){
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
#' @param motifPos To show reverse strand or not
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
    no_axes    <- list(theme(line            = ggplot2::element_blank(),
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
#' @param matrix Binding matrix in either energy, frequency, or information bit.
#' @param type One of 'energy', 'info', prob'
#' @param axes To show axes or not
#' @param motifPos To show reverse strand or not
#' @param labels To show labels or not
#' @return A binding motif logo
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

#how much is AA1>AA2 mutation in List1 explained by List2
AAChiSquare <- function(List1, List2, AA1, AA2, table = F){
  defaultW <- getOption("warn")
  options(warn = -1)
  pos1 <- which(List1 == AA1)
  pos2 <- which(List1 == AA2)
  chiSeq1 <- c()
  ChiS1 <- List1
  ChiS1[-pos1] <- '+'

  for(i2 in unique(List2[pos1])){
    ChiS2 <- List2
    ChiS2[ChiS2 != i2] <- '+'
    if(table){
      print(table(ChiS1, ChiS2))
    }

    chi <- chisq.test(ChiS1, ChiS2, correct = F)
    chiSeq1 <- c(chiSeq1, chi$p.value)
  }
  associatedAA1 <- unique(List2[pos1])[which(chiSeq1 == min(chiSeq1))]
  chiSeq2 <- c()
  ChiS1 <- List1
  ChiS1[-pos2] <- '+'
  for(i2 in unique(List2[pos2])){
    ChiS2 <- List2
    ChiS2[ChiS2 != i2] <- '+'
    if(table){
      print(table(ChiS1, ChiS2))
    }
    chi <- chisq.test(ChiS1, ChiS2, correct = F)
    chiSeq2 <- c(chiSeq2, chi$p.value)
  }
  associatedAA2 <- unique(List2[pos2])[which(chiSeq2 == min(chiSeq2))]
  options(warn = defaultW)
  if(associatedAA1 == associatedAA2){
    return(1 - max(c(min(chiSeq1),min(chiSeq2) )))
  }
  return(max(chiSeq))
}

#dependency between two mutations
AAmutChi <- function(List1, List2, AAfrom1, AAfrom2, AAto1, AAto2, table = F){
  defaultW <- getOption("warn")
  options(warn = -1)
  pos1 <- which(List1 == AAfrom1)
  pos2 <- which(List1 == AAto1)

  ChiS1 <- List1
  ChiS1[-pos1] <- '+'
  ChiS2 <- List2
  ChiS2[ChiS2 != AAfrom2] <- '+'
  if(table){
    print(table(ChiS1, ChiS2))
  }
  if(length(unique(ChiS1)) + length(unique(ChiS2)) != 4){
    chifrom <- 0.5
  }else{
    chifrom <- chisq.test(ChiS1, ChiS2, correct = F)$p.value
  }


  ChiS1 <- List1
  ChiS1[-pos2] <- '+'
  ChiS2 <- List2
  ChiS2[ChiS2 != AAto2] <- '+'
  if(table){
    print(table(ChiS1, ChiS2))
  }
  if(length(unique(ChiS1)) + length(unique(ChiS2)) != 4){
    chito <- 0.5
  }else{
    chito <- chisq.test(ChiS1, ChiS2, correct = F)$p.value
  }



  options(warn = defaultW)
  return(max(c(chifrom, chito)))
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
#' @param dddGList dddGList model output from form.dddGList()
#' @param referenceMatrix Binding matrix of the reference sample in ddG measurements
#' @param referenceSequence Protein sequence along the alignment of the reference sample
#' @param targetSequence Protien sequence along the alignmetn of the mutant sample
#' @param keyPos Key positions along the alignment to screen, must be included the dddGList
#' @param posMotif Position in the matrix to screen, must be included in the dddGList
#' @return Binding matrix of the mutant sample at the selected matrix position
#' @export
predict.dddGList <- function(dddGList, referenceMatrix, referenceSequence, targetSequence,  keyPos = c(13,14,5,8),posMotif = 'P-1'){
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


#ddG matrix of single AA
form.ddGFeatureList <- function(bHLH_pbAlignment, bHLH_motifs, keyPos = c(13,14,5,26,8), leaveOut = c()){
  bHLH_pbAlignment_backup <- bHLH_pbAlignment
  #start here
  bHLH_pbAlignment <- bHLH_pbAlignment_backup

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

  #create AAtypesList
  AAtypes <- data.frame(placeHolder = c(1:nrow(bHLH_pbAlignment)))
  for(i in 1:length(keyPos)){
    AApos <- keyPos[i]
    addPos <- substr(bHLH_pbAlignment$alignment,AApos,AApos)
    AAtypes <- cbind.data.frame(AAtypes, addPos)
  }
  AAtypes <- AAtypes[,-1]
  colnames(AAtypes) <- paste0('AA', keyPos)

  #create dddG
  dddGList <- list()
  svdList <- list()
  for(p in 1:length(keyPos)){
    AApos <- keyPos[p]
    for(posMotif in c('P1', 'P-1')){
      ProteinPos <- data.frame(name = bHLH_pbAlignment$name, pos = substr(bHLH_pbAlignment$alignment,AApos,AApos))
      PosAA <- c()
      for(i in 1:length(bHLH_motifs)){
        PosAA <- c(PosAA, ProteinPos[ProteinPos$name == bHLH_motifs[[i]]$name, 2])
      }
      pos_matrix <- gene2pos(bHLH_motifs, pos = posMotif, nrow = 4)
      colnames(pos_matrix) <- PosAA
      AAs <- unique(colnames(pos_matrix))
      dddGdf <- data.frame(placeHolder = c(0,0,0,0))
      for(aa1 in AAs){
        filteredPosMatrix <- pos_matrix[,colnames(pos_matrix) == aa1]
        mean1 <- data.frame(ddG = apply(data.frame(matrix2tetrahedron(filteredPosMatrix)), 2, function(column) mean(column)))
        mean1 <- unlist(mean1)%*%ittm
        mean1 <- data.frame(ddG = as.numeric(mean1))

        mean1 <- mean1 - max(mean1) + 1
        mean1 <- log(mean1)
        dddG <- data.frame(mean1)


        colnames(dddG) <- paste0(aa1)
        dddGdf <- cbind.data.frame(dddGdf, dddG)
        PMsvd <- svd(filteredPosMatrix)
        svdList[[paste0(AApos, '+', posMotif, '+', aa1)]] <- PMsvd
      }
      dddGdf <- dddGdf[,-1]
      dfName <- paste0('AA', AApos, '+', posMotif)
      dddGList[[dfName]] <- dddGdf
    }
  }
  out<- list()
  out$dddGList <- dddGList
  out$AAtypes <- AAtypes
  out$svd <- svdList
  return(out)
}

#unused
predEigen <- function(targetSequence, ddGFeature, keyPos = c(13,14,5,26,8), motifPos = 'P-1', print.dist = F, distType = 'Euclidean', singularity = 'one'){
  predKeyAA <- c()
  for(i in keyPos){
    predKeyAA <- c(predKeyAA, substr(targetSequence, i, i))
  }

  Energys <- data.frame(NULL)
  Eigens <- data.frame(NULL)
  for(i in 1:length(keyPos)){
    aa1 <- predKeyAA[i]
    pos <- keyPos[i]
    if(length(grep(aa1, colnames(ddGFeature$dddGList[[paste0('AA',pos,'+',motifPos)]]))) == 0){
      next()
    }
    Energy1 <- ddGFeature$dddGList[[paste0('AA',pos,'+',motifPos)]][,aa1]

    EigenM <- ddGFeature$svd[[paste0(pos, '+', motifPos, '+', aa1)]]$u
    EigenV <- ddGFeature$svd[[paste0(pos, '+', motifPos, '+', aa1)]]$d
    if(singularity == 'one'){
      Eigen1 <- EigenM[,1]
    }else if(singularity == 'all'){
      Eigen1 <- c(0,0,0,0)
      for(e in 1:ncol(EigenM)){
        Eigen1 <- Eigen1 + EigenM[,e] * EigenV[e]
      }
    }
    Energys <- rbind.data.frame(Energys, Energy1)
    Eigens <- rbind.data.frame(Eigens, Eigen1)
  }

  if(print.dist){
    print(dist(Eigens))
  }
  tempSum <- Energys[1,]
  for(i in 2:nrow(Energys)){
    if(distType == 'Euclidean'){
      dists <- dist(Eigens[c(1:i),])
    }else if (distType == 'Cosine'){
      dists <- c()
      cosM <- lsa::cosine(t(as.matrix(Eigens)))
      for(r in 1:nrow(cosM)){
        for(c in 1:ncol(cosM)){
          if(r > c){
            dists <- c(dists, cosM[r,c])
          }
        }
      }
      dists <- 1-dists
    }else if (distType == 'Angular'){
      dists <- c()
      cosM <- acos(lsa::cosine(t(as.matrix(Eigens))))*2/pi
      for(r in 1:nrow(cosM)){
        for(c in 1:ncol(cosM)){
          if(r > c){
            dists <- c(dists, min(cosM[r,c],1))
          }
        }
      }
    }

    minDist <- min(dists[(length(dists)-i+2):length(dists)])
    tempSum <- tempSum + minDist * Energys[i,]
  }
  pred <- data.frame(ddG = as.numeric(tempSum))
  PredM <- data.frame(apply(pred, 2, function(column) column - mean(column)))
  return(PredM)
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
#' @return A matrix of p-values with nrow = number of principal components in SVD, ncol = number of positions in the alignment
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
  return(svd.aov)
}


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
