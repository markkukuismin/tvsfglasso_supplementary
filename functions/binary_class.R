binary_class = function(A, A_est){
  
  A_est = A_est[upper.tri(A_est)]
  A = A[upper.tri(A)]
  
  #TP = sum(A_est*A)/2
  TP = sum(A_est*A)
  
  #FP = sum(A_est*(!A))/2
  FP = sum(A_est*(!A))
  
  #FN = sum((!A_est)*A)/2
  FN = sum((!A_est)*A)
  
  #TN = sum((!A_est)*(!A))/2
  TN = sum((!A_est)*(!A))
  
  TP = as.numeric(TP)
  FP = as.numeric(FP)
  FN = as.numeric(FN)
  TN = as.numeric(TN)
  
  FDR = FP/(FP + TP)
  
  if(FP + TP == 0) FDR = 0
  
  TPR = TP/(TP + FN) 
  
  Precision = TP/(TP + FP) 
  
  FPR = FP/(FP + TN)
  
  B = sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))
  
  MCC = (TP*TN - FP*FN)/B
  
  if(B == 0) MCC = 0
  
  F1 = (2*TP)/(2*TP + FP + FN)
  
  JI = TP/(TP + FN + FP) #Jaccard index
  
  ED = FP + FN # Edge distance
  
  results = list(TP = TP, 
                 FP = FP, 
                 FN = FN, 
                 TN = TN,
                 FDR = FDR,
                 TPR = TPR,
                 Precision = Precision,
                 FPR = FPR,
                 MCC = MCC,
                 F1 = F1,
                 JI = JI,
                 ED = ED)
  
  return(results)
  
}