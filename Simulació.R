library(ggplot2)
library(dplyr)
library(tidyr)

#----Simulació-----

#GENERAL----

#Funció per comptar estats
estats_possibles <- function(N) {
  estats_possibles <- expand.grid(replicate(N, c(0, 1), simplify = FALSE)) # Totes les combinacions binàries possibles de mida N
  return(estats_possibles)
}


#Funció per comptar fronteres
comptar_fronteres <- function(estat) {
  sum(estat != c(tail(estat, -1), estat[1])) # Compta canvis de valor entre veïns (amb frontera circular)
}

#Funció per generar ck
genera_ck <- function(L, epsilon = NULL, alpha = NULL) {
  if (L == 2) {
    if (is.null(epsilon)) stop("Per L = 2 cal especificar 'epsilon'")
    ck <- c(epsilon, 0.5, 1 - epsilon)
  } else if (L == 3) {
    if (is.null(epsilon) || is.null(alpha)) stop("Per L = 3 calen 'epsilon' i 'alpha'")
    ck <- c(epsilon, alpha, 1 - alpha, 1 - epsilon)
  } else {
    dim_ck_mig <- floor((L + 1) / 2)
    ck_mig <- runif(dim_ck_mig, min = 1e-10, max = 1 - 1e-10) #Creem la meitat del vector c de manera aleatòria
    ck_mig2 <- 1 - rev(ck_mig) #Segona meitat del vector
    
    if ((L + 1) %% 2 == 1) {
      ck <- c(ck_mig, 0.5, ck_mig2) #Si L+1 és senar el valor del mig és 0.5
    } else {
      ck <- c(ck_mig, ck_mig2)
    }
  }
  return(ck)
}


#--------------------------------------------------------------------------
#APARTAT 3. Mida del barri igual a la mida de l'autòmat (L=N)
#--------------------------------------------------------------------------
#funció agrupada per nombre d'1's----
midaLigualN <- function(N, iteracions) {
  set.seed(123)
  
  ck <- genera_ck(N)
  
  # Estat inicial
  estat <- sample(c(0, 1), N, replace = TRUE)
  
  # Simulació
  freq_k <- rep(0, N + 1) #inicialitzem el comptador de freqüències
  for (iter in 1:iteracions) {
    cela_triada <- sample(1:N, 1) #triem una cel·la aleatòriament
    k <- sum(estat) #sumem el veinatge de la cel·la triada, en aquest cas L=N
    estat[cela_triada] <- rbinom(1, 1, ck[k + 1]) #actualització segons ck 
    freq_k[sum(estat) + 1] <- freq_k[sum(estat) + 1] + 1 #actualitzem freqüència
  }
  pi_simulada <- freq_k / sum(freq_k) #normalitzem
  
  # Càlcul de la distribució teòrica
  if (N %% 2 == 0) {
    k_mig <- N / 2 #valor central
    pi_mig <- rep(0, k_mig) #inicialitzem el vector
    for (k in 0:(k_mig - 1)) {
      binomi <- choose(N, k)
      prod1 <- if (k == 0) 1 else prod(ck[1:k])
      prod2 <- if (k == k_mig) 1 else prod(1 - ck[(k + 2):(k_mig + 1)])
      pi_mig[k + 1] <- binomi * prod1 * prod2
    }
    pi_mig2 <- rev(pi_mig) #part simètrica
    pi_central <- choose(N, k_mig) * prod(ck[1:k_mig]) #definim el valor central ja que L+1 és senar
    pi_teorica <- c(pi_mig, pi_central, pi_mig2)
  } else { #cas quan L+1 és parell
    k_mig <- floor( N / 2)
    pi_mig <- rep(0, k_mig + 1)
    for (k in 0:k_mig) {
      binomi <- choose(N, k)
      prod1 <- if (k == 0) 1 else prod(ck[1:k])
      prod2 <- if (k == k_mig) 1 else prod(1 - ck[(k + 2):(k_mig + 1)])
      pi_mig[k + 1] <- binomi * prod1 * prod2
    }
    pi_mig2 <- rev(pi_mig)
    pi_teorica <- c(pi_mig, pi_mig2)
  }
  
  pi_teorica <- pi_teorica / sum(pi_teorica) #normalitzem
  
  return(list(pi_simulada = pi_simulada, pi_teorica = pi_teorica))
}

#exemple----
N <- 8
iteracions <- 100000
resultats <- midaLigualN(N, iteracions)
resultats
#GRÀFICA:----
k_valors <- 0:N #eix x
pi_simulada <- resultats$pi_simulada
pi_teorica <- resultats$pi_teorica


plot(k_valors, pi_simulada, type = "b", col = "orange", pch = 16, lwd = 2,
     xlab = "k (nombre de 1s)", ylab = expression(widehat(pi)[k]),
     main = paste("Distribució simulada vs teòrica (N = L =", N, ",iter =", iteracions, ")"),
     ylim = range(c(pi_simulada, pi_teorica))) #gràfica simulada

lines(k_valors, pi_teorica, type = "b", col = "steelblue", lty = 2, pch = 1, lwd = 2) #gràfica teòrica

legend("topright", legend = c("Simulada", "Teòrica"),
       col = c("orange", "steelblue"), pch = c(16, 1), lty = c(1, 2),
       bty = "n") #Llegenda


#--------------------------------------------------------------------------
#APARTAT 4. Cas especial: només dos veïns per cel·la (L=2)
#--------------------------------------------------------------------------
#funció agrupada per fronteres----
midaLigual2 <- function(N, epsilon, iteracions = 100000){
  set.seed(124)
  
  ck <- genera_ck(2, epsilon = epsilon)
  ck
  
  # Estat inicial
  estat <- sample(c(0, 1), N, replace = TRUE)
  
  # Simulació
  fronteres_possibles <- seq(0, ifelse(N %% 2 == 0, N, N - 1), by = 2) #Nombres de fronteres possibles (números parells)
  freq_fronteres <- setNames(rep(0, length(fronteres_possibles)), fronteres_possibles) #inicialitzacio frequencia
  
  for (iter in 1:iteracions) {
    cela_triada <- sample(1:N, 1) #triem una cel·la aleatoriament
    seguent <- ifelse(cela_triada == N, 1, cela_triada + 1) #cel·la veïna
    k <- estat[cela_triada] + estat[seguent] #suma d'1's en el veinatge
    estat[cela_triada] <- rbinom(1, 1, ck[k + 1]) #actualitzacio segons ck
    f <- comptar_fronteres(estat) #comptar el nombre de fronteres del nou estat
    freq_fronteres[as.character(f)] <- freq_fronteres[as.character(f)] + 1 #actualitzem frequencia de fronteres
  }
  pi_simulada_fronteres <- freq_fronteres / sum(freq_fronteres) #normalitzar
  
  # Càlcul de la distribució teòrica
  estats <- estats_possibles(N) #generacio de tots els estats possibles
  fronteres <- apply(estats, 1, comptar_fronteres) #comptar el nombre de fronteres per cadacun dels estats
  taula_fronteres <- table(fronteres) #taula de frequencies
  nfronteres_possibles <- as.numeric(names(taula_fronteres)) #nombre de fronteres
  freq_fronteres <- as.numeric(taula_fronteres) #frequencia
  
  epsilon <- ck[1]
  pi_i <- freq_fronteres * (2*epsilon)^(nfronteres_possibles/2)
  
  pi_teorica_fronteres <- pi_i / sum(pi_i) #normalitzar
  names(pi_teorica_fronteres) <- names(pi_simulada_fronteres) #assignem noms
  
  return(list(pi_simulada = pi_simulada_fronteres, pi_teorica = pi_teorica_fronteres))
}

#exemple----
N <- 8
epsilon <- runif(1, min = 1e-10, max = 1 - 1e-10)
resultats <- midaLigual2(N, epsilon)
resultats

#GRÀFICA:----
fronteres_possibles <- as.numeric(names(resultats$pi_simulada))
pi_simulada <- resultats$pi_simulada
pi_teorica <- resultats$pi_teorica

# Gràfic
plot(fronteres_possibles, pi_simulada, type = "b", col = "orange", pch = 16, lwd = 2,
     xlab = "Nombre de fronteres", ylab = "Probabilitat",
     main = bquote("Distribució estacionària" ~ 
                     (N == .(N) ~ "," ~ L == .(2) ~ "i" ~ epsilon == .(round(epsilon, 3)))),
     xaxt = "n",
     ylim = range(c(pi_simulada, pi_teorica))) #grafic distribucio simulada agrupada segons el nombre de frequencies

axis(1, at = seq(0, N, by = 2))
lines(fronteres_possibles, pi_teorica, type = "b", col = "steelblue", lty = 2, pch = 1, lwd = 2) #grafic distribucio teorica agrupada segons el nombre de frequencies

legend("topright", legend = c("Simulada", "Teòrica"),
       col = c("orange", "steelblue"), pch = c(16, 1), lty = c(1, 2),
       bty = "n") #llegenda

#Exemple 4.9----
#1. Comparació del nombre de fronteres----

epsilon <- seq(0, 1, by = 0.01) #valors d'epsilon

# Càlcul comparació (4 fronteres vs 2 fronteres, 4 fronteres vs 0 fronteres)
fronteres_4_vs_2 <- 2 * epsilon #quocient 4/2
fronteres_4_vs_0 <- 4 * epsilon^2 #quocient 4/0

plot(epsilon, fronteres_4_vs_2, type = "l", col = "purple", lwd = 2,
     xlab = expression(epsilon), ylab = "Ràtios de probabilitats estacionàries",
     main = "Comparació de probabilitats estacionàries",
     ylim = c(0, 4))
lines(epsilon, fronteres_4_vs_0, col = "orange", lwd = 2)
legend("topleft", legend = c("4 fronteres vs 2 fronteres", "4 fronteres vs 0 fronteres"),
       col = c("purple", "orange"), lwd = 2)
abline(h = 1, v = 0.5, lty = 2, col = "gray")


#2. Distribució estacionària segons el nombre de fronteres----
N <- 4
estats <- estats_possibles(N) #possibles estats
fronteres <- apply(estats, 1, comptar_fronteres) #nombre de fronteres per cada estat
taula_fronteres <- table(fronteres) #taula de frequencies
nfronteres_possibles <- as.numeric(names(taula_fronteres)) #nombre de fronteres
freq_fronteres <- as.numeric(taula_fronteres) #frequència

epsilon_valors <- seq(0.1, 0.9, by = 0.1) #Valors d’epsilon
pi_list <- list() #creem llista per per guardar el valors
for (i in seq_along(nfronteres_possibles)) {
  f <- nfronteres_possibles[i]
  pi_list[[paste0("f", f)]] <- (2 * epsilon_valors)^(f / 2) #distribucio estacionaria no normalitzada
}
print(pi_list)

# Calcular Z per despres normalizar
Z <- rep(0, length(epsilon_valors))
for (i in seq_along(nfronteres_possibles)) {
  f <- nfronteres_possibles[i]
  freq <- freq_fronteres[i]
  pi_f <- pi_list[[paste0("f", f)]]
  Z <- Z + freq * pi_f
}
print(Z)

# Normalitzar per cada valor de pi_f
pi_norm_list <- lapply(pi_list, function(pi_f) pi_f / Z) 

df_plot <- as.data.frame(pi_norm_list) #probabilitat estacionaria per un nombre concret de fronteres i epsilon
df_plot$epsilon <- epsilon_valors  # afegim columna amb els valors d'epsilon
df_plot <- df_plot[, c("epsilon", names(pi_norm_list))]  # posem epsilon al davant

# Gràfic
matplot(df_plot$epsilon, df_plot[, -1], type = "o", lty = 1, pch = 1,
        col = rainbow(ncol(df_plot) - 1),
        xlab = expression(epsilon), ylab = "Probabilitat estacionària",
        main = paste("Distribució estacionària (N =", N, ")"))

# Llegenda
legend("topright", legend = gsub("f", "", names(pi_norm_list)), 
       col = rainbow(ncol(df_plot) - 1), lty = 1, pch = 1,
       title = "Fronteres")


#3. Simulació vs teòrica----
set.seed(124)
N <- 4
L <- 2
iteracions <- 100000
ck_mig_valors <- c(0.025, 0.075, 0.1, 0.5, 0.9) #triem els epsilons
colors <- c("darkorange", "forestgreen", "purple", "steelblue", "darkgrey")
fronteres_possibles <- seq(0, ifelse(N %% 2 == 0, N, N - 1), by = 2) #nombre de fronteres possibles

#Gràfic
plot(NA, xlim = range(fronteres_possibles), ylim = c(0, 1),
     xlab = "Nombre de fronteres", ylab = "Probabilitat",
     main = bquote("Simulació vs teoria per diferents valors d'" * epsilon * ", amb N =" ~ .(N)),
     xaxt = "n") #base del gràfic
axis(1, at = fronteres_possibles[fronteres_possibles %% 2 == 0])


for (i in seq_along(ck_mig_valors)) {
  epsilon <- ck_mig_valors[i]
  color <- colors[i]
  
  resultat <- midaLigual2(N, epsilon, iteracions) #apliquem funcio
  
  x_valors <- as.numeric(names(resultat$pi_simulada)) #convertim en nombre l'etiqueta de les fronteres pel gràfic
  
  lines(x_valors, resultat$pi_simulada, type = "b", col = color, pch = 16, lwd = 2) #grafic simulada
  lines(x_valors, resultat$pi_teorica, type = "b", col = "black", lty = 2, pch = 1, lwd = 1) #grafic teorica
}


legend_labels <- paste0("Simulada: ε = ", ck_mig_valors)
legend("topright",
       inset = c(-0.001, 0),
       xpd = TRUE,
       legend = c(legend_labels, "Teòriques"),
       col = c(colors, "black"),
       lty = c(rep(1, length(ck_mig_valors)), 2),
       pch = c(rep(16, length(ck_mig_valors)), 1),
       lwd = c(rep(2, length(ck_mig_valors)), 1),
       bty = "n",
       cex = 0.75,
       y.intersp = 0.8) #llegenda


#--------------------------------------------------------------------------
#APARTAT 5. Mida del barri igual a tres (L=3)
#--------------------------------------------------------------------------
#funció agrupada per fronteres----
midaLigual3 <- function(N, epsilon, alpha, iteracions = 100000){
  set.seed(127)
  
  ck <- genera_ck(3, epsilon = epsilon, alpha = alpha)
  ck
  # Estat inicial
  estat <- sample(c(0, 1), N, replace = TRUE)
  
  # Simulació
  fronteres_possibles <- seq(0, ifelse(N %% 2 == 0, N, N - 1), by = 2) 
  freq_fronteres <- setNames(rep(0, length(fronteres_possibles)), fronteres_possibles) 
  
  for (iter in 1:iteracions) {
    cela_triada <- sample(1:N, 1)
    anterior <- ifelse(cela_triada == 1, N, cela_triada - 1)
    seguent  <- ifelse(cela_triada == N, 1, cela_triada + 1)
    k <- estat[anterior] + estat[cela_triada] + estat[seguent]
    estat[cela_triada] <- rbinom(1, 1, ck[k + 1])
    f <- comptar_fronteres(estat)
    freq_fronteres[as.character(f)] <- freq_fronteres[as.character(f)] + 1
  }
  
  pi_simulada_fronteres <- freq_fronteres / sum(freq_fronteres)
  
  # Càlcul de la distribució teòrica
  estats <- estats_possibles(N)
  fronteres <- apply(estats, 1, comptar_fronteres)
  taula_fronteres <- table(fronteres)
  nfronteres_possibles <- as.numeric(names(taula_fronteres))
  freq_fronteres <- as.numeric(taula_fronteres)
  
  epsilon <- ck[1]
  alpha <- ck[2]
  pi_i <- freq_fronteres * (epsilon/(1-alpha))^(nfronteres_possibles/2)
  
  pi_teorica_fronteres <- pi_i / sum(pi_i)
  names(pi_teorica_fronteres) <- names(pi_simulada_fronteres)
  
  return(list(pi_simulada = pi_simulada_fronteres, pi_teorica = pi_teorica_fronteres))
}
#exemple----
N <- 7
epsilon <- runif(1, min = 1e-10, max = 1 - 1e-10)
alpha <- runif(1, min = 1e-10, max = 1 - 1e-10)
resultats <- midaLigual3(N, epsilon, alpha)
resultats


#GRÀFICA:----
fronteres_possibles <- as.numeric(names(resultats$pi_simulada))
pi_simulada <- resultats$pi_simulada
pi_teorica <- resultats$pi_teorica

# Gràfic
plot(fronteres_possibles, pi_simulada, type = "b", col = "orange", pch = 16, lwd = 2,
     xlab = "Nombre de fronteres", ylab = "Probabilitat",
     main = bquote("Distribució estacionària" ~ 
                     (N == .(N) ~ "," ~ L == .(3) ~ "," ~ epsilon == .(round(epsilon, 3)) ~ "i" ~ alpha == .(round(alpha, 3)))),
     xaxt = "n",
     ylim = range(c(pi_simulada, pi_teorica)))

axis(1, at = seq(0, N, by = 2))
lines(fronteres_possibles, pi_teorica, type = "b", col = "steelblue", lty = 2, pch = 1, lwd = 2)

legend("topleft", legend = c("Simulada", "Teòrica"),
       col = c("orange", "steelblue"), pch = c(16, 1), lty = c(1, 2),
       bty = "n")




#--------------------------------------------------------------------------
#APARTAT 6. Simulació d'altres escenaris
#--------------------------------------------------------------------------
#1. Canvi de la ck (sense simetria) ----
#Vector prob aleatòri
genera_ck_aleatori <- function(L){
  if(L<1)stop("L ha de ser un enter més gran que 1")
  ck <- runif(L+1)
  return(ck)
}

#Vector prob ascendent
genera_ck_ascendent <- function(L){
  if(L<1)stop("L ha de ser un enter més gran que 1")
  valors <- runif(L+1)
  ck <- sort(valors)#ordenar
  return(ck)
}

#Vector prob descendent
genera_ck_descendent <- function(L){
  if(L<1)stop("L ha de ser un enter més gran que 1")
  valors <- runif(L+1)
  ck <- rev(sort(valors))#ordenar i posar del reves
  return(ck)
}

#Aplicació 1: Mida del barri igual a la mida de l'autòmat (L=N)----
#funció agrupada per nombre d'1's----
midaLigualN_ck <- function(N, iteracions) {
  
  ck <- genera_ck_aleatori(N) #ho podem canviar per alguna de les funcions anteriors (genera_ck_aleatori, genera_ck_ascendent, genera_ck_descendent)
  
  # Estat inicial
  estat <- sample(c(0, 1), N, replace = TRUE)
  
  # Simulació
  freq_k <- rep(0, N + 1)
  for (iter in 1:iteracions) {
    cela_triada <- sample(1:N, 1)
    k <- sum(estat)
    estat[cela_triada] <- rbinom(1, 1, ck[k + 1])
    freq_k[sum(estat) + 1] <- freq_k[sum(estat) + 1] + 1
  }
  pi_simulada <- freq_k / sum(freq_k)
  
  # Càlcul de la distribució teòrica
  pi_teorica <- rep(0, N + 1)
  
  for (k in 0:N) { #igual que la secció tres per el valor k arriba fins a N, ja que no tenim simetria
    binomi <- choose(N, k)
    prod1 <- if (k == 0) 1 else prod(ck[1:k])
    prod2 <- if (k == N) 1 else prod(1 - ck[(k + 2):(N + 1)])
    pi_teorica[k + 1] <- binomi * prod1 * prod2
  }
  pi_teorica <- pi_teorica / sum(pi_teorica)
  
  return(list(pi_simulada = pi_simulada, pi_teorica = pi_teorica, ck=ck))
}

#exemple----
N <- 7
iteracions <- 100000
resultats <- midaLigualN_ck(N, iteracions)
resultats
#GRÀFICA:----
k_valors <- 0:N #eix x
pi_simulada <- resultats$pi_simulada
pi_teorica <- resultats$pi_teorica


plot(k_valors, pi_simulada, type = "b", col = "orange", pch = 16, lwd = 2,
     xlab = "k (nombre de 1s)", ylab = expression(widehat(pi)[k]),
     main = paste("Distribució simulada vs teòrica (N = L =", N, ")"),
     ylim = range(c(pi_simulada, pi_teorica)))

lines(k_valors, pi_teorica, type = "b", col = "steelblue", lty = 2, pch = 1, lwd = 2)

legend("topleft", legend = c("Simulada", "Teòrica"),
       col = c("orange", "steelblue"), pch = c(16, 1), lty = c(1, 2),
       bty = "n")










#Aplicació 2: Cas especial: només dos veïns per cel·la (L=2)----
#funció agrupada per fronteres----
midaLigual2_ck <- function(N, iteracions = 100000){
  
  ck <- genera_ck_descendent(2) #ho podem canviar per alguna de les funcions anteriors (genera_ck_aleatori, genera_ck_ascendent, genera_ck_descendent)
  
  
  # Estat inicial
  estat <- sample(c(0, 1), N, replace = TRUE)
  
  # Simulació
  fronteres_possibles <- seq(0, ifelse(N %% 2 == 0, N, N - 1), by = 2)
  freq_fronteres <- setNames(rep(0, length(fronteres_possibles)), fronteres_possibles)
  
  for (iter in 1:iteracions) {
    cela_triada <- sample(1:N, 1)
    seguent <- ifelse(cela_triada == N, 1, cela_triada + 1)
    k <- estat[cela_triada] + estat[seguent]
    estat[cela_triada] <- rbinom(1, 1, ck[k + 1])
    f <- comptar_fronteres(estat)
    freq_fronteres[as.character(f)] <- freq_fronteres[as.character(f)] + 1
  }
  pi_simulada_fronteres <- freq_fronteres / sum(freq_fronteres)
  
  # Càlcul de la distribució teòrica
  estats <- estats_possibles(N)
  fronteres <- apply(estats, 1, comptar_fronteres)
  taula_fronteres <- table(fronteres)
  nfronteres_possibles <- as.numeric(names(taula_fronteres))
  freq_fronteres <- as.numeric(taula_fronteres)
  
  epsilon <- ck[1]
  pi_i <- freq_fronteres * (2*epsilon)^(nfronteres_possibles/2)
  
  pi_teorica_fronteres <- pi_i / sum(pi_i)
  names(pi_teorica_fronteres) <- names(pi_simulada_fronteres)
  
  return(list(pi_simulada = pi_simulada_fronteres, pi_teorica = pi_teorica_fronteres, ck=ck))
}

#exemple----
N <- 8
resultats <- midaLigual2_ck(N)
resultats

#GRÀFICA:----
fronteres_possibles <- as.numeric(names(resultats$pi_simulada))
pi_simulada <- resultats$pi_simulada
pi_teorica <- resultats$pi_teorica
epsilon <- resultats$ck[1]


# Gràfic
plot(fronteres_possibles, pi_simulada, type = "b", col = "orange", pch = 16, lwd = 2,
     xlab = "Nombre de fronteres", ylab = "Probabilitat",
     main = bquote("Distribució estacionària" ~ 
                     (N == .(N) ~ "," ~ L == .(2) ~ "i" ~ epsilon == .(round(epsilon, 3)))),
     xaxt = "n",
     ylim = range(c(pi_simulada, pi_teorica)))

axis(1, at = seq(0, N, by = 2))
lines(fronteres_possibles, pi_teorica, type = "b", col = "steelblue", lty = 2, pch = 1, lwd = 2)

legend("topright", legend = c("Simulada", "Teòrica"),
       col = c("orange", "steelblue"), pch = c(16, 1), lty = c(1, 2),
       bty = "n")




#Aplicació 3: Mida del barri igual a tres (L=3)----
#funció agrupada per fronteres----
midaLigual3_ck <- function(N, iteracions = 100000){
  
  ck <- genera_ck_ascendent(3) #ho podem canviar per alguna de les funcions anteriors (genera_ck_aleatori, genera_ck_ascendent, genera_ck_descendent)
  ck
  # Estat inicial
  estat <- sample(c(0, 1), N, replace = TRUE)
  
  # Simulació
  fronteres_possibles <- seq(0, ifelse(N %% 2 == 0, N, N - 1), by = 2)
  freq_fronteres <- setNames(rep(0, length(fronteres_possibles)), fronteres_possibles)
  
  for (iter in 1:iteracions) {
    cela_triada <- sample(1:N, 1)
    anterior <- ifelse(cela_triada == 1, N, cela_triada - 1)
    seguent  <- ifelse(cela_triada == N, 1, cela_triada + 1)
    k <- estat[anterior] + estat[cela_triada] + estat[seguent]
    estat[cela_triada] <- rbinom(1, 1, ck[k + 1])
    f <- comptar_fronteres(estat)
    freq_fronteres[as.character(f)] <- freq_fronteres[as.character(f)] + 1
  }
  
  pi_simulada_fronteres <- freq_fronteres / sum(freq_fronteres)
  
  # Càlcul de la distribució teòrica
  estats <- estats_possibles(N)
  fronteres <- apply(estats, 1, comptar_fronteres)
  taula_fronteres <- table(fronteres)
  nfronteres_possibles <- as.numeric(names(taula_fronteres))
  freq_fronteres <- as.numeric(taula_fronteres)
  
  epsilon <- ck[1]
  alpha <- ck[2]
  pi_i <- freq_fronteres * (epsilon/(1-alpha))^(nfronteres_possibles/2)
  
  pi_teorica_fronteres <- pi_i / sum(pi_i)
  names(pi_teorica_fronteres) <- names(pi_simulada_fronteres)
  
  return(list(pi_simulada = pi_simulada_fronteres, pi_teorica = pi_teorica_fronteres, ck=ck))
}
#exemple----
N <- 8
resultats <- midaLigual3_ck(N)
resultats

#GRÀFICA:----
fronteres_possibles <- as.numeric(names(resultats$pi_simulada))
pi_simulada <- resultats$pi_simulada
pi_teorica <- resultats$pi_teorica
epsilon <- resultats$ck[1]
alpha <- resultats$ck[2]


# Gràfic
plot(fronteres_possibles, pi_simulada, type = "b", col = "orange", pch = 16, lwd = 2,
     xlab = "Nombre de fronteres", ylab = "Probabilitat",
     main = bquote("Distribució estacionària" ~ 
                     (N == .(N) ~ "," ~ L == .(3) ~ "," ~ epsilon == .(round(epsilon, 3)) ~ "i" ~ alpha == .(round(alpha, 3)))),
     xaxt = "n",
     ylim = range(c(pi_simulada, pi_teorica)))

axis(1, at = seq(0, N, by = 2))
lines(fronteres_possibles, pi_teorica, type = "b", col = "steelblue", lty = 2, pch = 1, lwd = 2)

legend("topright", legend = c("Simulada", "Teòrica"),
       col = c("orange", "steelblue"), pch = c(16, 1), lty = c(1, 2),
       bty = "n")



#--------------------------------------------------------------------------
#2. Augmentem el veinatge (L=5)----
#funció per comptar les 2fronteres
comptar_2fronteres <- function(estat) {
  N <- length(estat)
  comptador <- 0
  for (i in 1:N) {
    j <- (i + 2 - 1) %% N + 1  # posició a distància 2, circular
    if (estat[i] != estat[j]) {
      comptador <- comptador + 1
    }
  }
  return(comptador)
}

midaLigual5 <- function(N, epsilon, alpha, gamma, iteracions = 100000){
  if (N < 5) stop("El valor de N ha de ser un enter major que 4")
  ck <- c(epsilon, alpha, gamma, 1-gamma, 1-alpha, 1-epsilon) #creem el vector de probabilitats
  
  # Estat inicial
  estat <- sample(c(0, 1), N, replace = TRUE)
  
  # Simulació
  fronteres_possibles <- seq(0, ifelse(N %% 2 == 0, N, N - 1), by = 2)
  freq_fronteres <- setNames(rep(0, 2 * N + 1), 0:(2 * N)) #com que no sabem quins possibles valors pot tenir b1+b2, fem un vector general de 2N+1 components
  
  for (iter in 1:iteracions) {
    cela_triada <- sample(1:N, 1)
    i_v1 <- ifelse(cela_triada - 2 < 1, N + cela_triada - 2, cela_triada - 2)#indexem tots els veïns
    i_v2 <- ifelse(cela_triada - 1 < 1, N + cela_triada - 1, cela_triada - 1)
    i_v3 <- ifelse(cela_triada + 1 > N, cela_triada + 1 - N, cela_triada + 1)
    i_v4 <- ifelse(cela_triada + 2 > N, cela_triada + 2 - N, cela_triada + 2)
    k <- estat[i_v1] + estat[i_v2] + estat[cela_triada] + estat[i_v3] + estat[i_v4] #modifiquem el valor k per tenir en compte els 5 veïns
    estat[cela_triada] <- rbinom(1, 1, ck[k + 1])
    b1 <- comptar_fronteres(estat) #comptar quantes 1fronteres té el nou estat
    b2 <- comptar_2fronteres(estat) #comptar quantes 2fronteres té el nou estat
    frontera_total <- b1 + b2 #agrupem segons la suma de 1fronteres + 2fronteres
    
    freq_fronteres[as.character(frontera_total)] <- freq_fronteres[as.character(frontera_total)] + 1 #sumem un al comptador
  }
  
  pi_simulada_fronteres <- freq_fronteres / sum(freq_fronteres) #normalitzem
  pi_simulada_fronteres <- pi_simulada_fronteres[pi_simulada_fronteres > 0] #treiem les els valors de freq_fronteres que siguin 0
  
  # Càlcul de la distribució teòrica
  estats <- estats_possibles(N)
  # Comptem les 1-fronteres i 2-fronteres per cada estat
  b1 <- apply(estats, 1, comptar_fronteres)       # comptem els possibles valors de 1-fronteres
  b2 <- apply(estats, 1, comptar_2fronteres)      # comptem els possibles valors de 2-fronteres
  
  # Càlcul no normalitzat segons el teorema
  alpha <- ck[2]
  gamma <- ck[3]
  pi_no_normalitzada <- (alpha / (1 - gamma))^((b1 + b2)/2)
  
  # Agrupem pel nombre total de fronteres b1+b2
  fronteres_totals <- b1 + b2 
  taula <- tapply(pi_no_normalitzada, fronteres_totals, sum) #agrupem segons el nombre total de fronteres
  
  # Normalitzem
  pi_teorica_fronteres <- taula / sum(taula)
  names(pi_teorica_fronteres) <- names(pi_simulada_fronteres) #assignem noms
  
  return(list(pi_simulada = pi_simulada_fronteres, pi_teorica = pi_teorica_fronteres))
}

#exemple amb relacio----
N <- 7
alpha <- runif(1, min = 1e-10, max = 1 - 1e-10)
gamma <- runif(1, min = 1e-10, max = 1 - 0.1) # no pot ser molt gran perquè sinó epsilon dona error
epsilon <- (alpha^2 * (1 - alpha)) / (1 - gamma)^2
resultats <- midaLigual5(N, epsilon, alpha, gamma)
resultats

#exemple sense relacio----
N <- 8
alpha <- runif(1, min = 1e-10, max = 1 - 1e-10)
gamma <- runif(1, min = 1e-10, max = 1 - 1e-10)
epsilon <- runif(1, min = 1e-10, max = 1 - 1e-10) #epsilon s'escull aleatoriament
#si volem que es compleixi la relacio hem de posar: 
#epsilon <- (alpha^2 * (1 - alpha)) / ((1 - gamma)^2) 
resultats <- midaLigual5(N, epsilon, alpha, gamma)
resultats

#GRÀFICA: ----
fronteres_possibles <- as.numeric(names(resultats$pi_simulada))
pi_simulada <- resultats$pi_simulada
pi_teorica <- resultats$pi_teorica

# Gràfic
plot(fronteres_possibles, pi_simulada, type = "b", col = "orange", pch = 16, lwd = 2,
     xlab = "1-fronteres + 2-fronteres", ylab = "Probabilitat",
     main = bquote(
       atop("Distribució estacionària",
            N == .(N) ~ "," ~ L == .(5) ~ "," ~ 
              epsilon == .(round(epsilon, 3)) ~ "," ~ 
              alpha == .(round(alpha, 3)) ~ "," ~ 
              gamma == .(round(gamma, 3)))
     ),
     xaxt = "n",
     ylim = range(c(pi_simulada, pi_teorica))
)

axis(1, at = as.numeric(names(resultats$pi_simulada)))
lines(fronteres_possibles, pi_teorica, type = "b", col = "steelblue", lty = 2, pch = 1, lwd = 2)

legend("topleft", legend = c("Simulada", "Teòrica"),
       col = c("orange", "steelblue"), pch = c(16, 1), lty = c(1, 2),
       bty = "n")







#--------------------------------------------------------------------------
#3. Actualització síncrona----
#Aplicació 1: Mida del barri igual a la mida de l'autòmat (L=N)----
midaLigualN_sincron <- function(N, iteracions) {
  
  ck <- genera_ck(N)
  
  # Estat inicial
  estat <- sample(c(0, 1), N, replace = TRUE)
  
  # Simulació
  freq_k <- rep(0, N + 1)
  
  for (iter in 1:iteracions) {
    estat_nou <- estat
    
    for (i in 1:N) {#a cada iteració recorrem totes les cel·les i mirem si canviarien segons ck
      k <- sum(estat)
      
      p_canvi <- ck[k + 1]
      
      estat_nou[i] <- rbinom(1, 1, p_canvi)#ho guardem a un vector nou per tal de que la decisió no afecti als altres estats
    }
    
    estat <- estat_nou #un cop acabat actualitzem de cop
    
    freq_k[sum(estat) + 1] <- freq_k[sum(estat) + 1] + 1 #sumem el nombre d'1's per cada estat resultant
  }
  
  pi_simulada <- freq_k / sum(freq_k) #normalitzem
  
  # Distribució teòrica igual que la asíncrona (no tenim resultat per la síncrona)
  if (N %% 2 == 0) {
    k_mig <- N / 2
    pi_mig <- rep(0, k_mig)
    for (k in 0:(k_mig - 1)) {
      binomi <- choose(N, k)
      prod1 <- if (k == 0) 1 else prod(ck[1:k])
      prod2 <- if (k == k_mig) 1 else prod(1 - ck[(k + 2):(k_mig + 1)])
      pi_mig[k + 1] <- binomi * prod1 * prod2
    }
    pi_mig2 <- rev(pi_mig)
    pi_central <- choose(N, k_mig) * prod(ck[1:k_mig])
    pi_teorica <- c(pi_mig, pi_central, pi_mig2)
  } else {
    k_mig <- floor(N / 2)
    pi_mig <- rep(0, k_mig + 1)
    for (k in 0:k_mig) {
      binomi <- choose(N, k)
      prod1 <- if (k == 0) 1 else prod(ck[1:k])
      prod2 <- if (k == k_mig) 1 else prod(1 - ck[(k + 2):(k_mig + 1)])
      pi_mig[k + 1] <- binomi * prod1 * prod2
    }
    pi_mig2 <- rev(pi_mig)
    pi_teorica <- c(pi_mig, pi_mig2)
  }
  
  pi_teorica <- pi_teorica / sum(pi_teorica)
  
  return(list(pi_simulada = pi_simulada, pi_teorica = pi_teorica))
}

#exemple----
N <- 7
iteracions <- 100000
resultats <- midaLigualN_sincron(N, iteracions)
resultats
#GRÀFICA:----
k_valors <- 0:N
pi_simulada <- resultats$pi_simulada
pi_teorica <- resultats$pi_teorica


plot(k_valors, pi_simulada, type = "b", col = "orange", pch = 16, lwd = 2,
     xlab = "k (nombre de 1s)", ylab = expression(widehat(pi)[k]),
     main = paste("Distribució estacionària (N = L =", N, ")"),
     ylim = range(c(pi_simulada, pi_teorica)))

lines(k_valors, pi_teorica, type = "b", col = "steelblue", lty = 2, pch = 1, lwd = 2)

legend("topright", legend = c("Simulada Síncrona", "Teòrica Asincrona"),
       col = c("orange", "steelblue"), pch = c(16, 1), lty = c(1, 2),
       bty = "n")


#Aplicació 2: Cas especial: només dos veïns per cel·la (L=2)----
midaLigual2_sincron <- function(N, epsilon, iteracions = 100000){
  
  ck <- genera_ck(2, epsilon = epsilon)
  ck
  
  # Estat inicial
  estat <- sample(c(0, 1), N, replace = TRUE)
  
  # Simulació
  fronteres_possibles <- seq(0, ifelse(N %% 2 == 0, N, N - 1), by = 2) #nombre de fronteres possibles
  freq_fronteres <- setNames(rep(0, length(fronteres_possibles)), fronteres_possibles) #inicialitzem frequencia
  
  for (iter in 1:iteracions) {
    estat_nou <- estat
    
    for (i in 1:N) { #per cada iteració recorrem totes les cel·les mirant si canviarien d'estat i ho guardem a un vector nou per no influir a la resta de decisions
      seguent <- ifelse(i == N, 1, i + 1)
      k <- estat[i] + estat[seguent]
      p_canvi <- ck[k + 1]
      
      estat_nou[i] <- rbinom(1, 1, p_canvi)
    }
    estat <- estat_nou #actualitzem de cop
    f <- comptar_fronteres(estat) #comptem el nombre de fronteres l'estat resultant
    freq_fronteres[as.character(f)] <- freq_fronteres[as.character(f)] + 1 #comptabilitzem a la frequencia
  }
  pi_simulada_fronteres <- freq_fronteres / sum(freq_fronteres)
  
  # Càlcul de la distribució teòrica igual que la asíncrona (no tenim resultat teòric per la síncrona)
  estats <- estats_possibles(N)
  fronteres <- apply(estats, 1, comptar_fronteres)
  taula_fronteres <- table(fronteres)
  nfronteres_possibles <- as.numeric(names(taula_fronteres))
  freq_fronteres <- as.numeric(taula_fronteres)
  
  epsilon <- ck[1]
  pi_i <- freq_fronteres * (2*epsilon)^(nfronteres_possibles/2)
  
  pi_teorica_fronteres <- pi_i / sum(pi_i)
  names(pi_teorica_fronteres) <- names(pi_simulada_fronteres)
  
  return(list(pi_simulada = pi_simulada_fronteres, pi_teorica = pi_teorica_fronteres))
}

#exemple----
N <- 7
epsilon <- runif(1, min = 1e-10, max = 1 - 1e-10)
resultats <- midaLigual2_sincron(N, epsilon)
resultats

#GRÀFICA:----
fronteres_possibles <- as.numeric(names(resultats$pi_simulada))
pi_simulada <- resultats$pi_simulada
pi_teorica <- resultats$pi_teorica

# Gràfic
plot(fronteres_possibles, pi_simulada, type = "b", col = "orange", pch = 16, lwd = 2,
     xlab = "Nombre de fronteres", ylab = "Probabilitat",
     main = bquote("Distribució estacionària" ~ 
                     (N == .(N) ~ "," ~ L == .(2) ~ "i" ~ epsilon == .(round(epsilon, 3)))),
     xaxt = "n",
     ylim = range(c(pi_simulada, pi_teorica)))

axis(1, at = seq(0, N, by = 2))
lines(fronteres_possibles, pi_teorica, type = "b", col = "steelblue", lty = 2, pch = 1, lwd = 2)

legend("topright", legend = c("Simulada Síncrona", "Teòrica Asíncrona"),
       col = c("orange", "steelblue"), pch = c(16, 1), lty = c(1, 2),
       bty = "n")





#Línia comparativa - (x =epsilon) ----
N <- 4
epsilon_valors <- seq(0.1, 0.9, by = 0.02)#prenem molts possibles valors que pot prendre epsilon

diferenciesvt <- numeric(length(epsilon_valors)) #inicialitzem variació total
diferenciesl2 <- numeric(length(epsilon_valors)) #inicialitzem distancia l_2
names(diferenciesvt) <- paste0("ε=", round(epsilon_valors, 2))#creem etiquetes pels valors d'epsilon
names(diferenciesl2) <- paste0("ε=", round(epsilon_valors, 2))

for (j in seq_along(epsilon_valors)) {
  epsilon <- epsilon_valors[j] #per cada valor d'epsilon pres calculem la distribució estacionaria simulada sincrona i la teòrica asíncrona
  
  resultat <- midaLigual2_sincron(N, epsilon)
  pi_s <- resultat$pi_simulada
  pi_t <- resultat$pi_teorica
  
  # Distància de variació total
  dist_vt <- 0.5 * sum(abs(pi_s - pi_t)) #calculem la variacio total
  diferenciesvt[j] <- dist_vt #guardem al vector diferencies
  dist_l2 <- sqrt(sum((pi_s - pi_t)^2)) #calculem la distancia l_2
  diferenciesl2[j] <- dist_l2 #guardem al vector diferencies
}
print(diferenciesvt)
print(diferenciesl2)

# Convertim a dataframe per visualitzar
diferencies_df <- data.frame(
  epsilon = epsilon_valors,
  Diferenciavt = diferenciesvt,
  Diferencial2 = diferenciesl2
)


# Gràfic variació total
ggplot(diferencies_df, aes(x = epsilon, y = Diferenciavt)) +
  geom_col(fill = "orchid4", width = 0.015) +
  labs(
    title = "Variació total entre síncron i asíncron segons ε (N=4, L=2)",
    x = expression(epsilon),
    y = "Distància de variació total"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold")
  )

#Gràfic distància L2
ggplot(diferencies_df, aes(x = epsilon, y = Diferencial2)) +
  geom_col(fill = "orchid4", width = 0.015) +
  labs(
    title = "Distància L2 entre síncron i asíncron segons ε (N=4, L=2)",
    x = expression(epsilon),
    y = "Distància de variació total"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold")
  )

#Gràfic comparació de distàncies
ggplot(diferencies_df, aes(x = epsilon)) +
  geom_line(aes(y = Diferenciavt, color = "Variació total"), size = 1.2) +
  geom_line(aes(y = Diferencial2, color = "ℓ_2"), size = 1.2) +
  labs(
    title = "Comparació de distàncies",
    x = expression(epsilon),
    y = "Valor de la distància",
    color = "Distància \nestadística"
  ) +
  theme_minimal()

#Aplicació 3: Mida del barri igual a tres (L=3)----
midaLigual3_sincron <- function(N, epsilon, alpha, iteracions = 100000){
  
  ck <- genera_ck(3, epsilon = epsilon, alpha = alpha)
  ck
  # Estat inicial
  estat <- sample(c(0, 1), N, replace = TRUE)
  
  # Simulació
  fronteres_possibles <- seq(0, ifelse(N %% 2 == 0, N, N - 1), by = 2) #nombres possibles de fronteres
  freq_fronteres <- setNames(rep(0, length(fronteres_possibles)), fronteres_possibles) #inicialitzem freq fronteres
  
  for (iter in 1:iteracions) {
    estat_nou <- estat
    
    for (i in 1:N) { #per cada iteració mirem totes les cel·les i guardem el canvi segons ck en un vector nou per no afectar a la resta de resultats
      anterior <- ifelse(i == 1, N, i - 1)
      seguent <- ifelse(i == N, 1, i + 1)
      k <- estat[anterior] + estat[i] + estat[seguent]
      p_canvi <- ck[k + 1]
      
      estat_nou[i] <- rbinom(1, 1, p_canvi)
    }
    estat <- estat_nou #actualitzem tot l'estat de cop
    f <- comptar_fronteres(estat) #comptem quantes fronteres te l'estat nou
    freq_fronteres[as.character(f)] <- freq_fronteres[as.character(f)] + 1 #comptabilitzem el valor
  }
  pi_simulada_fronteres <- freq_fronteres / sum(freq_fronteres) #normalitzem
  
  # Càlcul de la distribució teòrica igual que la asíncrona (no tenim resultat teoric per la síncrona)
  estats <- estats_possibles(N)
  fronteres <- apply(estats, 1, comptar_fronteres)
  taula_fronteres <- table(fronteres)
  nfronteres_possibles <- as.numeric(names(taula_fronteres))
  freq_fronteres <- as.numeric(taula_fronteres)
  
  epsilon <- ck[1]
  alpha <- ck[2]
  pi_i <- freq_fronteres * (epsilon/(1-alpha))^(nfronteres_possibles/2)
  
  pi_teorica_fronteres <- pi_i / sum(pi_i)
  names(pi_teorica_fronteres) <- names(pi_simulada_fronteres)
  
  return(list(pi_simulada = pi_simulada_fronteres, pi_teorica = pi_teorica_fronteres))
}
#exemple----
N <- 7
epsilon <- 0.485
alpha <- 0.539
resultats <- midaLigual3_sincron(N, epsilon, alpha)
resultats


#GRÀFICA:----
fronteres_possibles <- as.numeric(names(resultats$pi_simulada))
pi_simulada <- resultats$pi_simulada
pi_teorica <- resultats$pi_teorica

# Gràfic
plot(fronteres_possibles, pi_simulada, type = "b", col = "orange", pch = 16, lwd = 2,
     xlab = "Nombre de fronteres", ylab = "Probabilitat",
     main = bquote("Distribució estacionària" ~ 
                     (N == .(N) ~ "," ~ L == .(3) ~ "," ~ epsilon == .(round(epsilon, 3)) ~ "i" ~ alpha == .(round(alpha, 3)))),
     xaxt = "n",
     ylim = range(c(pi_simulada, pi_teorica)))

axis(1, at = seq(0, N, by = 2))
lines(fronteres_possibles, pi_teorica, type = "b", col = "steelblue", lty = 2, pch = 1, lwd = 2)

legend("topleft", legend = c("Simulada Síncrona", "Teòrica Asíncrona"),
       col = c("orange", "steelblue"), pch = c(16, 1), lty = c(1, 2),
       bty = "n")

#Taula comparativa - (x = epsilon, y = alpha) - heatmap----
N <- 4
epsilon_valors <- seq(0.1, 0.9, by = 0.1) #possibles valors que pot prendre epsilon
alpha_valors <- seq(0.1, 0.9, by = 0.1) #possibles valors que pot prendre alpha

fronteres_possibles <- seq(0, ifelse(N %% 2 == 0, N, N - 1), by = 2) #nombre de fronteres possibles

diferenciesvt <- matrix(NA, nrow = length(alpha_valors), ncol = length(epsilon_valors),
                        dimnames = list(paste0("α=", alpha_valors), paste0("ε=", epsilon_valors)))#inicialitzem matriu per guardar variacio total
diferenciesl2 <- matrix(NA, nrow = length(alpha_valors), ncol = length(epsilon_valors),
                        dimnames = list(paste0("α=", alpha_valors), paste0("ε=", epsilon_valors)))#inicialitzem matriu per guarda l_2

for (i in seq_along(alpha_valors)) {
  for (j in seq_along(epsilon_valors)) {#calculem per cada valor epsilon i alpha la distribució estacionaria simulada síncrona i la teòrica asíncrona
    epsilon <- epsilon_valors[j]
    alpha <- alpha_valors[i]
    
    resultat <- midaLigual3_sincron(N, epsilon, alpha)
    pi_s <- resultat$pi_simulada
    pi_t <- resultat$pi_teorica
    
    # Calcular distància de variació total
    dist_vt <- 0.5 * sum(abs(pi_s - pi_t))#calculem la variació total
    diferenciesvt[i, j] <- dist_vt
    dist_l2 <- sqrt(sum((pi_s - pi_t)^2)) #calculem l2
    diferenciesl2[i, j] <- dist_l2
  }
}
print(diferenciesvt)
print(diferenciesl2)

# Mapes de calor

#passem primer les dades a data.frame
diferenciesvt_df <- as.data.frame(diferenciesvt)
diferenciesl2_df <- as.data.frame(diferenciesl2)
diferenciesvt_df$alpha <- rownames(diferenciesvt_df)
diferenciesl2_df$alpha <- rownames(diferenciesl2_df)

#les passem a pivot_longer per graficar
taulavt <- diferenciesvt_df %>% 
  pivot_longer(cols = starts_with("ε="), names_to = "epsilon", values_to = "Diferencia")

taulal2 <- diferenciesl2_df %>%
  pivot_longer(cols = starts_with("ε="), names_to = "epsilon", values_to = "Diferencia")

#variació total
ggplot(taulavt, aes(x = epsilon, y = factor(alpha, levels = rev(paste0("α=", alpha_valors))), fill = Diferencia)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "white", mid = "darkmagenta", high = "darkslateblue",
    midpoint = 0.3  
  ) +
  labs(
    title = "Variació total entre síncron i asíncron segons ε i α (N=4, L=3)",
    x = expression(epsilon),
    y = expression(alpha),
    fill = "Distància"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 12, face = "bold")
  )

#distància L2
ggplot(taulal2, aes(x = epsilon, y = factor(alpha, levels = rev(paste0("α=", alpha_valors))), fill = Diferencia)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "white", mid = "darkmagenta", high = "darkslateblue",
    midpoint = 0.3  
  ) +
  labs(
    title = "Distància L2 entre síncron i asíncron segons ε i α (N=4, L=3)",
    x = expression(epsilon),
    y = expression(alpha),
    fill = "Distància"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 12, face = "bold")
  )

#comparació entre distàncies:
diff_mat <- abs(diferenciesvt - diferenciesl2) #calculem les diferencies element per element per cada parell (epsilon, alpha)
print(diff_mat)

# Passar a dataframe per ggplot
diff_df <- as.data.frame(diff_mat)
diff_df$alpha <- rownames(diff_df)
taula_long <- diff_df %>%
  pivot_longer(cols = starts_with("ε="), names_to = "epsilon", values_to = "Diferencia")

ggplot(taula_long, aes(x = epsilon, y = factor(alpha, levels = rev(paste0("α=", alpha_valors))), fill = Diferencia)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "white", mid = "darkolivegreen", high = "black",
    midpoint = 0.115  
  ) +
  labs(
    title = "Diferència entre distàncies",
    x = expression(epsilon),
    y = expression(alpha),
    fill = bquote("|VT - ℓ" * symbol("")[2] * "|")
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 12, face = "bold")
  )

#--------------------------------------------------------------------------
#APARTAT 7. Aplicació al model exponencial de votant
#--------------------------------------------------------------------------
#funcions generals ----
estats_possibles2 <- function(N) { #creem la funció estats possibles per el cas on els valors són 1 o -1
  estats_possibles <- expand.grid(replicate(N, c(-1, 1), simplify = FALSE))
  return(estats_possibles)
}

energia <- function(s) {#funció energia tenint en compte circularitat
  N <- length(s)
  suma <- 0
  for (i in 1:N) {
    si <- s[i]
    si1 <- if (i == N) s[1] else s[i + 1]
    suma <- suma + si * (si - si1)
  }
  return(suma)
}

prob_canvi <- function(beta, av){ #funció de probabilitat de canvi en funció de la mitjana
  prob <- (1/(1+exp(-log(exp(beta)-1)*av)))
  return(prob)
}

#funció agrupada per fronteres (L=2)----
mitjana2 <- function(s, i){ #funció mitjana tenint en compte circularitat
  N <- length(s)
  av <- if(i==N) mean(c(s[N], s[1])) else mean(c(s[i], s[i+1]))
  return(av)
}

votantL2 <- function(N, beta, iteracions = 100000){
  set.seed(132)
  
  # Estat inicial
  estat <- sample(c(-1, 1), N, replace = TRUE)
  
  # Simulació
  fronteres_possibles <- seq(0, ifelse(N %% 2 == 0, N, N - 1), by = 2) #nombres possibles de fronteres 
  freq_fronteres <- setNames(rep(0, length(fronteres_possibles)), fronteres_possibles) #inicialitzem freq fronteres
  
  for (iter in 1:iteracions) {
    cela_triada <- sample(1:N, 1)
    k <- mitjana2(estat, cela_triada) #k és la mitjana de la cel·la triada
    estat[cela_triada] <- sample(c(-1, 1), 1, prob = c(1 - prob_canvi(beta, k), prob_canvi(beta, k))) #actualitzem cel·la triada,-1 amb probabilitat 1 - prob_canvi i 1 amb probabilitat prob_canvi 
    f <- comptar_fronteres(estat) #compta nombre fronteres que té el nou estat
    freq_fronteres[as.character(f)] <- freq_fronteres[as.character(f)] + 1 #comptabilitza el nombre de frontera
  }
  
  pi_simulada_fronteres <- freq_fronteres / sum(freq_fronteres) #normalizem
  
  # Càlcul de la distribució teòrica
  estats <- estats_possibles2(N) #generació de tots els possibles estats
  energies <- apply(estats, 1, energia) #calcul de l'energia per cada estat
  taula_energies <- table(energies) #taula de freq de l'energia
  nenergies_possibles <- as.numeric(names(taula_energies)) #valor possible d'energia
  freq_energies <- as.numeric(taula_energies) #freq de l'energia
  
  pi_i <- freq_energies * exp(-(beta - log(2))*(nenergies_possibles/4))
  pi_teorica_fronteres <- pi_i / sum(pi_i)
  names(pi_teorica_fronteres) <- names(pi_simulada_fronteres)
  
  return(list(pi_simulada = pi_simulada_fronteres, pi_teorica = pi_teorica_fronteres))
}

#Exemple----
N <- 8
beta <- 1.6
resultats <- votantL2(N, beta)

#GRÀFICA:----
fronteres_possibles <- as.numeric(names(resultats$pi_simulada))
pi_simulada <- resultats$pi_simulada
pi_teorica <- resultats$pi_teorica

plot(fronteres_possibles, pi_simulada, type = "b", col = "orange", pch = 16, lwd = 2,
     xlab = "Nombre de fronteres", ylab = "Probabilitat",
     main = bquote("Distribució estacionària:" ~ 
                     (N == .(N) ~ "," ~ L == .(2) ~ "i" ~ beta == .(round(beta, 3)))),
     xaxt = "n",
     ylim = range(c(pi_simulada, pi_teorica)))

axis(1, at = seq(0, N, by = 2))
lines(fronteres_possibles, pi_teorica, type = "b", col = "steelblue", lty = 2, pch = 1, lwd = 2)

legend("topright", legend = c("Simulada", "Teòrica"),
       col = c("orange", "steelblue"), pch = c(16, 1), lty = c(1, 2),
       bty = "n")


#funció agrupada per fronteres (L=3)----
mitjana3 <- function(s, i) {
  N <- length(s)
  anterior <- if (i == 1) s[N] else s[i - 1]
  actual <- s[i]
  seguent <- if (i == N) s[1] else s[i + 1]
  return(mean(c(anterior, actual, seguent)))
}


votantL3 <- function(N, beta, iteracions = 100000){
  set.seed(122)
  
  # Estat inicial
  estat <- sample(c(-1, 1), N, replace = TRUE)
  
  # Simulació
  fronteres_possibles <- seq(0, ifelse(N %% 2 == 0, N, N - 1), by = 2)
  freq_fronteres <- setNames(rep(0, length(fronteres_possibles)), fronteres_possibles)
  
  for (iter in 1:iteracions) {
    cela_triada <- sample(1:N, 1)
    k <- mitjana3(estat, cela_triada)
    estat[cela_triada] <- sample(c(-1, 1), 1, prob = c(1 - prob_canvi(beta, k), prob_canvi(beta, k)))
    f <- comptar_fronteres(estat)
    freq_fronteres[as.character(f)] <- freq_fronteres[as.character(f)] + 1
  }
  
  pi_simulada_fronteres <- freq_fronteres / sum(freq_fronteres)
  
  # Càlcul de la distribució teòrica
  estats <- estats_possibles2(N)
  energies <- apply(estats, 1, energia)
  taula_energies <- table(energies)
  nenergies_possibles <- as.numeric(names(taula_energies))
  freq_energies <- as.numeric(taula_energies)
  
  pi_i <- freq_energies * exp(-(beta - log(1 + (exp(beta)-1)^(-1/3)))*(nenergies_possibles/4))
  pi_teorica_fronteres <- pi_i / sum(pi_i)
  names(pi_teorica_fronteres) <- names(pi_simulada_fronteres)
  
  return(list(pi_simulada = pi_simulada_fronteres, pi_teorica = pi_teorica_fronteres))
}

#Exemple----
N <- 8
beta <- 1.2
resultats <- votantL3(N, beta)
resultats

#GRÀFICA:----
fronteres_possibles <- as.numeric(names(resultats$pi_simulada))
pi_simulada <- resultats$pi_simulada
pi_teorica <- resultats$pi_teorica


plot(fronteres_possibles, pi_simulada, type = "b", col = "orange", pch = 16, lwd = 2,
     xlab = "Nombre de fronteres", ylab = "Probabilitat",
     main = bquote("Distribució estacionària" ~ 
                     (N == .(N) ~ "," ~ L == .(3) ~ "i" ~ beta == .(round(beta, 3)))),
     xaxt = "n",
     ylim = range(c(pi_simulada, pi_teorica)))

axis(1, at = seq(0, N, by = 2))
lines(fronteres_possibles, pi_teorica, type = "b", col = "steelblue", lty = 2, pch = 1, lwd = 2)

legend("topright", legend = c("Simulada", "Teòrica"),
       col = c("orange", "steelblue"), pch = c(16, 1), lty = c(1, 2),
       bty = "n")

#COMPARACIÓ BETA CANVIANT----
betas <- c(0.1, 1, 4) #valors de beta
N <- 8

par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))  #per mostrar 6 gràfics alhora

# L = 2
for (beta in betas) {
  resultats <- votantL2(N, beta)
  fronteres_possibles <- as.numeric(names(resultats$pi_simulada))
  pi_simulada <- resultats$pi_simulada
  pi_teorica <- resultats$pi_teorica
  epsilon <- 1/exp(beta)
  
  plot(fronteres_possibles, pi_simulada, type = "b", col = "orange", pch = 16, lwd = 2,
       xlab = "Fronteres", ylab = "Probabilitat",
       main = bquote(L == 2 ~ "i" ~ beta == .(beta)),
       xaxt = "n",
       ylim = range(c(pi_simulada, pi_teorica)))
  axis(1, at = seq(0, N, by = 2))
  lines(fronteres_possibles, pi_teorica, type = "b", col = "steelblue", lty = 2, pch = 1, lwd = 2)
}

# L = 3
for (beta in betas) {
  resultats <- votantL3(N, beta)
  fronteres_possibles <- as.numeric(names(resultats$pi_simulada))
  pi_simulada <- resultats$pi_simulada
  pi_teorica <- resultats$pi_teorica
  ck <- genera_ck(3, 1/exp(beta), 1 / (1 + (exp(beta) - 1)^(1/3)))
  epsilon <- 1/exp(beta)
  alpha <- 1/(1+(exp(beta)-1)^(1/3))
  
  plot(fronteres_possibles, pi_simulada, type = "b", col = "orange", pch = 16, lwd = 2,
       xlab = "Fronteres", ylab = "Probabilitat",
       main = bquote(L == 3 ~ "i" ~ beta == .(beta)),
       xaxt = "n",
       ylim = range(c(pi_simulada, pi_teorica)))
  axis(1, at = seq(0, N, by = 2))
  lines(fronteres_possibles, pi_teorica, type = "b", col = "steelblue", lty = 2, pch = 1, lwd = 2)
}



