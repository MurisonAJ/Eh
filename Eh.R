JSD<-function(input, base) {
  # Compute the Jensen-Shannon Divergence for a matrix of discrete probability distributions
  # Each row in the matrix should be a different probability distribution
  # All probability distributions must be continuous with respect to one another
  # Equal weighting is assumed for all probability vectors
  N<-nrow(input)
  weighted_input<-input
  for(i in 1:N) {
    
    weighted_input[i,]<-input[i,]/N
    
  }

  
  left_term<-H(apply(weighted_input, 2, sum), base)
  right_term<-0
  for (i in 1:N) {
    right_term<-right_term+H(input[i,], base)
  }
  return(left_term-(right_term/N))
}

H<-function(input, base) {  
  ## Compute Shannon's Entryopy to a given base for a given input vector
  current<-0
  for (i in 1:length(input)) {
    if(input[i]!=0) {
      current<-current+(input[i]*log(input[i],base))
    }
  }
  return(-1 * current)   
  }

create_track<-function(chr, start, end, length, sample) {
# Given a chromosome, start, end, region length and sample coverage table
# Output a track for that sample and region.
  
  track<-vector(length=length)
  region<-sample[which(sample[,1]==chr),]
  region<-region[which(as.numeric(region[,3])<=end),]
  region<-region[which(as.numeric(region[,2])>=start),]

  for(i in 1:nrow(region)) {
    a<-max(start, as.numeric(region[i,2]))
    a<-a-start+1
    b<-min(end, as.numeric(region[i,3]))
    b<-min(b-start+1, length)
    if(is.na(a)) { return(rep(0,length)) }
    else if (is.na(b)) { return(rep(0,length))}
    else{
      track[a:b]<-as.numeric(region[i,4])
    }
    }
  
  return(track[1:length])
  
}

### Read Region List

read_regions<-function(filename) {
  ### Read in a .bed format file 
  ### Chr, Start, End (1 based), 4 is an identifier and any number of arbitrary decorative columns are allowed
  dat<-read.table(file=filename, sep="\t", header=FALSE, stringsAsFactors=FALSE) 
  return(dat)
}

get_unique_regions<-function(regions) {
    unr<-rownames(unique(regions[,1:3]))
    return(regions[unr,])
  }
  
continuity_check<-function(coverage_matrix) {
  prune<-c()
  for (i in 1:ncol(coverage_matrix)) {
    if(min(coverage_matrix[,i])==0) {
      prune<-c(prune,i)
    }
  }
  if (length(prune)>0) {coverage_matrix<-coverage_matrix[,-prune]}
  return(coverage_matrix)
}

average_bins<-function(track, binsize) {
    new_track<-vector(length=(length(track)/binsize))
    for(i in 1:(length(track)/binsize)) {
      subtrack<-track[(1+((i-1)*binsize)):(i*binsize)]
      sum<-sum(subtrack)
      sum<-sum/binsize
      new_track[i]<-sum
    }
    return(new_track)
}

BatchJSD<-function(N, casecontrol, regions, regionsize, base=2, average=FALSE, rolling_bin_size=0, correct=FALSE, delta=1*(10^-10), 
                   prune=TRUE, minbins=(regionsize/max(rolling_bin_size, 1)/2)) {

  results<-vector(length=nrow(regions))
  support<-vector(length=nrow(regions))
  ## Work out if we are doing averages
  ## If so adjust bin size
  full_regionsize<-regionsize
  if(average) {
    ##Check we have a valid bin size
    if (regionsize%%rolling_bin_size>0) {print("Invalid Bin Size: Must exactly divide region size"); return(-1)}
    else {regionsize<-regionsize/rolling_bin_size}
  }
  
  for(i in 1:nrow(regions)) {
    ##Build Tracks for each specified region                        
    chr<-regions[i,1]
    start<-regions[i,2]
    end<-regions[i,3]
    current<-matrix(nrow=N, ncol=regionsize)
    
    for (j in 1:N) {
      
      track<-create_track(chr=chr, start=start, end=end, length=full_regionsize, sample=get(paste0(casecontrol,j)))
      track[which(is.na(track))]=0
      if (average) {
        current[j,]<-average_bins(track,rolling_bin_size)
      } else { current[j,]<-track}
    }
   
    ## Perform Sanity Checks on the region
    ErrorYN<-FALSE  
    if(is.na(sum(apply(current, 1,sum)))) { print(paste0(chr, ":",start,"-",end, "  ", "Region has NA values")); ErrorYN<-TRUE}
    
    # Apply Continuity Correction if requested
    if(!(ErrorYN)) { 
      if(correct) {
        current<-current+delta
      }
    
      # Prune discontinuous bins if requested
      if(prune) {
        current<-continuity_check(current)
      }
    }
    ## Check we still have sufficient bins
    z<-ncol(current)
    if (class(z)!="integer") {
      print(paste0(chr, ":",start,"-",end, "  ","Error Calculating Number of Columns"))
      ErrorYN<-TRUE
    } else if (z < minbins) { print(paste0(chr, ":",start,"-",end, "  ", regionsize-z, " bins have 0 values")); support[i]<-regionsize-z; ErrorYN<-TRUE}
    
    if(ErrorYN) {
      results[i]<-"NA"
    } else {
      
      currentP<-current
      for(j in 1:N) {
        currentP[j,]<-current[j,]/sum(current[j,])
      }
      
      results[i]<-JSD(currentP,base)
      print(paste0(chr, ":", start, "-", end, " : ", results[i]))
      support[i]<-regionsize-z;
    }
    
  }
  return(cbind(results, support))
}


contrast_matrix<-read.csv("contrasts.csv")
contrast_matrix[,2]<-contrast_matrix[,2]+1
name_vector<-c("TLR","THR")

  contrasts<-unique(contrast_matrix[,2])
  for (contrast in contrasts) {
    file_list<-contrast_matrix[which(contrast_matrix[,2]==contrast),1]
    name<-name_vector[contrast]
    for (i in 1:length(file_list)) {
      assign(paste0(name, i), read.table(as.character(file_list[i]), sep="\t", stringsAsFactors=FALSE))
    } 
  }

regions<-read.table("", sep="\t", stringsAsFactors=FALSE)
regions<-get_unique_regions(regions)
z<-unique(regions[,1])
for (i in 25:51) {
  cull<-which(regions[,1]==z[i])
  regions<-regions[-cull,]  
}
