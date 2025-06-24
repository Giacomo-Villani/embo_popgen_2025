#### let create neutral evolution , by Oscar Lao####


## forward method, brownian motion

next_gen_neutral <- function(Ne, F){
    coin <- rep(0,Ne)
    for (i in 1:Ne) {
	coin[i] <- rbinom(1, 1, F)
	}
	return(coin)
}

recording_of_frequency <- function(Ne, F, Ngen){
    record = rep(0,Ngen)
	for (i in 1:Ne) {
	pop= next_gen_neutral(Ne, F)
    f <- mean(pop)
	record[i] <- f
	}
	return(record)
}

Ne=1000
F= 1/Ne
Ngen = 1000
record <- recording_of_frequency(Ne, F, Ngen)

library(data.table)

record_to_plot <- data.table(gen = 1:Ngen, frequency= record)

plot(record_to_plot, type="b")

## backward method, coalescence