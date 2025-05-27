cv.med.rd.int.mensaje <- function(x,y,x.eval,val,xl,kr){
  # cv bandwidth treating the cutoof point as an interior point
  cut   <- quantile(abs(x-x.eval),probs=xl)
  index <- which(abs(x-x.eval) <= cut)
  cri   <- array(0,c(length(val),1))
  for (i in 1:length(val)){
    h.tau <- val[i]
    cri.h <- NULL
    for (j in index){
      xx=x[-j]
      yy=y[-j]
      Qy  <- lprq.pro.rdd(xx,yy,tau=0.5,x[j],h.tau,mon=0,kr)$b0
      cri.h <- c(cri.h,abs(y[j]-Qy))
      
      mensaje <- paste("Iteración", j, "completada con éxito.")
      cat(mensaje, "\n")
    }
    cri[i] <- mean(cri.h)
  }
  return(list(h.cv = val[which.min(cri)], cand = cbind(val,cri)))
}

cv.med.rd.int.mensaje(x = running2rmt1, y = outcomepre1_d2,
                      x.eval = 0, val = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08),
                      xl = 0.5, kr = 3) # se queda atrapado en un caso

# Probando con método que fuerza monotonicidad
bwd2preo1_afj <- rdd.bandwidth.AFJ(running2rmt1, outcomepre1_d2, d2rmt1,
                                   x.eval = 0, tt = c(0.25,0.75), m=3, method = c(3,4,5), 
                                   val=c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08), kr = 3)
