################################## embo-popgen day 2: genetic drift ##################################
N <- 50
fA <- 0.5
Ngen <- 100
FA_Ngen <- rep(0, Ngen)

################# random sampling, 100 replicates, from gen 0

for (i in 1:Ngen) {
				 FA_Ngen[i] <- rbinom(1, 2*N, fA) / (2*N)
}
mean(FA_Ngen)
gigi_100rep_plot <- plot(FA_Ngen, type="l", ylim=c(0,1), lwd=2)
gigi_100rep_hist <- hist(FA_Ngen)

################# random sampling from previous generation starting from gen 0 to Ngen

for (i in 1:Ngen) {
                 if (i==1) {
				 FA_Ngen[i] <- rbinom(1, 2*N, fA) / (2*N)
				 }
				 else {
				 FA_Ngen[i] <- rbinom(1, 2*N, FA_Ngen[i-1]) / (2*N)
				 }
}

mean(FA_Ngen)
gigi_100gen_plot <- plot(FA_Ngen, type="l", ylim=c(0,1), lwd=2)
gigi_100gen_hist <- hist(FA_Ngen)

## 20 different replicates of the experiment
Xrep <- 20
matrix_drift_20_rep <- matrix(nrow= Xrep, ncol = Ngen)
for (k in 1:Xrep) {
 for (i in 1:Ngen) {
                 if (i==1) {
				 matrix_drift_20_rep[k,i] <- rbinom(1, 2*N, fA) / (2*N)
				 }
				 else {
				 matrix_drift_20_rep[k,i] <- rbinom(1, 2*N, matrix_drift_20_rep[k,i-1]) / (2*N)
				 }
 }
}

matplot(t(matrix_drift_20_rep), type = "l", lty = 1, col = 1:5,
        xlab = "generation", ylab = "fA", main = "Replicates of random sampling")

Nfixations <- 0
for (k in 1:Xrep) {
 for (i in Ngen) {
                 if (matrix_drift_20_rep[k,i]==1) {
				 Nfixations <- Nfixations + 1
				 }
				 else if (matrix_drift_20_rep[k,i]==0){
				 Nfixations <- Nfixations + 1
				 }
 }
}

print(paste("Number of fixations of one allele:", Nfixations))


################# now we change Ne!

Xrep <- 20
matrix_drift_20_rep_variableNe <- matrix(nrow= Xrep, ncol = Ngen)
N <- 2

## trial 1 with incremental Ne
for (k in 1:Xrep) {
       if (k!=1) {N <- N + 10}
 for (i in 1:Ngen) {
                 if (i==1) {
				 matrix_drift_20_rep_variableNe[k,i] <- rbinom(1, 2*N, fA) / (2*N)
				 }
				 else {
				 matrix_drift_20_rep_variableNe[k,i] <- rbinom(1, 2*N, matrix_drift_20_rep_variableNe[k,i-1]) / (2*N)
				 }
 }
}

matplot(t(matrix_drift_20_rep_variableNe), type = "l", lty = 1, col = 1:5,
        xlab = "generation", ylab = "fA", main = "Replicates of random sampling")

## trial 2 with 10 small and 10 big Ne
N <- 4
for (k in 1:Xrep) {
       if (k>=10) {N <- 1000}
 for (i in 1:Ngen) {
                 if (i==1) {
				 matrix_drift_20_rep_variableNe[k,i] <- rbinom(1, 2*N, fA) / (2*N)
				 }
				 else {
				 matrix_drift_20_rep_variableNe[k,i] <- rbinom(1, 2*N, matrix_drift_20_rep_variableNe[k,i-1]) / (2*N)
				 }
 }
}

matplot(t(matrix_drift_20_rep_variableNe), type = "l", lty = 1, col = 1:5,
        xlab = "generation", ylab = "fA", main = "Replicates of random sampling")

################################## embo-popgen day 2: genetic drift and natural selection ##################################

simulateTrajectory <- function(s, N, t=500, nrepl=100, initFreq=1/(2*N)) {

        cat("2Ns =",2*N*s,"\n")

        # initialise frequencies
        fA <- matrix(NA, nrow=nrepl, ncol=t)
        # fA[,1] <- 1/(2*N)
        fA[,1] <- initFreq

        # viability
        vAA <- 1
        vAa <- 1 - s
        vaa <- 1 - (2*s)

        for (r in 1:nrepl) {

                for (i in 2:t) {

                        # selection
                        fpA <- fA[r,i-1] * (vAA*fA[r,i-1] + (vAa*(1-fA[r,i-1]))) / (vAA*fA[r,i-1]^2 + 2*vAa*fA[r,i-1]*(1-fA[r,i-1]) + vaa*(1-fA[r,i-1])^2)

                        if (fpA <= 0) { fA[r,i:t] <- 0; break} # lost
                        if (fpA >= 1) { fA[r,i:t] <- 1; break} # fixed

                        # drift
                        # fA[r,i] <- sum(sample(x=c(0,1), size=(2*N), replace=T, prob=c((1-fpA),fpA))) / (2*N)
                        fA[r,i] <- rbinom (n=1, size=2*N, prob=fpA) / (2*N)

                }

        }

        u <- 0
        if ((2*N*s) > -1) u <- 1/(2*N)
        if ((2*N*s) > 1) u <- 2*s
    
        cat("Lost = ", length(which(fA[,t]==0)), "\n")
        cat("Fixed = ", length(which(fA[,t]==1)), "\t (expected = ", (u*nrepl), ")\n")
    
        return(invisible(fA));

}

plotTrajectory <- function(fA, ylim=c(0,1), tlim=c(1,NA)) {
        cols <- colors()
        if (is.na(tlim[2])) tlim <- c(1,ncol(fA))
        plot(fA[1,],ylim=ylim,ty="l",xlim=tlim,col=cols[2],xlab="generations",ylab="frequency",lwd=2)
        for (i in 2:nrow(fA)) lines(fA[i,],type="l",col=cols[i+1],lwd=2)
}

plotTrajectory(simulateTrajectory(s=0.001, N=100, t=100, nrepl=100))