library(urca)
library(FitAR)
library(forecast)

#Projet de Séries temporelles
#Indice brut de la production d'ordinateurs et d'équipements périphériques en France de 1990 à 2017

#Partie 1 : Les données
#Chargement des donnees
donnees <-read.csv(file="/Users/johannalalou/Desktop/Valeurs.csv",encoding="latin1", sep=";",dec=",", header = TRUE)
#On crée un nouveau format de date
donnees$Date<-as.Date(paste0("01/",paste0(paste0(donnees$Mois,"/"),donnees$Annee)),format="%d/%m/%Y")
print(donnees)
#On change l'ordre des données pour que la 1ere composante du vecteur soit bien la 1ère donnée de 1990
donnees <- donnees[rev(rownames(donnees)) ,]
donnees

#On trace la série brute
plot(donnees$Date,donnees$Indice ,xlab="Date",ylab="Indice brut", type="l", main="Série des indices du nombre d'ordinateurs et d'équipements périphériques produits")
#Il semble qu'il y ait une saisonnalité (plutôt annuelle) et une tendance non linéaire (quadratique)
#La tendance à la hausse du début des années 90 jusqu'au début des années 2000 correspond à la décennie du développement d'Internet et au début de l'expansion des "technologies de l'information et de la communication". Les sociétés informatiques ont été très sollicitées durant cette période, et la production d'ordinateurs s'est développée à grande échelle.   

#On trace la transformation logarithmique de la série.
#En effet, les modeles ARMA sont assez restrictifs dans la mesure où ils supposent que les observations sont stationnaires. Il faut donc transformer ces données.
#La transformation log permet ici de stabiliser la variance de la série car celle-ci croit avec les valeurs la série. 
plot(donnees$Date,log(donnees$Indice),xlab="Date",ylab="Log de l'indice brut",type="l", main="Log de la série des indices du nombre d'ordinateurs et d'équipements périphériques produits")

#On travaille sur cette dernière série, et on regarde de plus près quelle peut être la période de la saisonnalité, en "zoomant" sur le comportement de la série.
plot(donnees$Date[190:230],log(donnees$Indice[190:230]),xlab="Date",ylab="Log de l'indice",type="l")
plot(donnees$Date[140:180],log(donnees$Indice[140:180]),xlab="Date",ylab="Log de l'indice",type="l")
#On observe bien une saisonnalité annuelle mais il semble qu'il y ait également une saisonnalité de période plus petite (trimestrielle), autrement dit deux saisonnalités.

#On trace les autocorrélogrammes (ACF et PACF de la série log-linéarisée)
#On observe des autocorrélations qui corroborent plutôt bien notre intuition sur la saisonnalité trimestrielle, mais moins sur la saisonnalité annuelle.
acf(log(donnees$Indice))
pacf(log(donnees$Indice))

#On cherche maintenant à stationnariser la série. On la différencie à l'ordre 3 pour enlever la saisonnalité trimestrielle et de surcroit la tendance, qu'on soupçonne polynomiale d'ordre 2.
plot(donnees$Date[1:322],diff(log(donnees$Indice),lag=3),xlab="Date",ylab="Log de l'indice brut en différences premières",type="l", main="Différences premières du log de la série de l'indice du nombre d'ordinateurs et d'équipements périphériques produits")

#On trace les autocorrélogrammes de la série différenciée à l'ordre 3
acf(diff(log(donnees$Indice),lag=3)) 
pacf(diff(log(donnees$Indice),lag=3))
#On remarque sur l'ACF de la série différenciée à l'ordre 3 qu'il subsiste des autocorrélations d'ordre 12 (provenant de la saisonnalité annuelle)

#On décide donc de différencier la série brute (d'origine) à l'odre 12 pour enlever la saisonnalité annuelle, ce qui permettra d'enlever également la saisonnalité trimestrielle.
plot(donnees$Date[1:313],diff(log(donnees$Indice),lag=12),xlab="Date",ylab="Log de l'indice brut en différences premières",type="l", main="Différences premières du log de la série de l'indice du nombre d'ordinateurs et d'équipements périphériques produits")

#On trace l'acf de la série différenciée à l'ordre 12. Les autocorrélations décroissant linéairement, on suppose qu'il y a une racine unitaire.
acf(diff(log(donnees$Indice),lag=12))

#La série désaisonnalisée semble avoir une faible tendance linéaire, nous redifférencions donc la série à l'ordre 1 (nous allons donc étudier un modèle ARMA plutôt qu'un modèle ARIMA avec d=1)
plot(donnees$Date[1:312],diff(diff(log(donnees$Indice),lag=12),1),xlab="Date",ylab="Log de l'indice brut en différences premières",type="l", main="Différences premières du log de la série de l'indice du nombre d'ordinateurs et d'équipements périphériques produits")
#Cette série semble plutôt stationnaire.

#On trace les autocorrélogrammes de la série stationnarisée ( qui nous confirment bien que la série semble bien stationnaire ) afin de déterminer les bornes suppérieures des ordres p et q, que nous noterons pmax et qmax, et qui permettront d'estimer les ordres p et q de notre ARMA(p,q). 
acf(diff(diff(log(donnees$Indice),lag=12),1)) #qmax= 12 car après, les autocorrélations ne sortent pas significativement de l'intervalle de confiance
pacf(diff(diff(log(donnees$Indice),lag=12),1)) #pmax= 12, même si le 24ème élément sort de l'intervalle de confiance, nous ne le considérons par soucis de simplicité et dans la mesure où celui-ci n'est pas extrêmement significatif.
diff <- diff(diff(log(donnees$Indice),lag=12),1)
plot(diff,type='l')


#Vérifier que la série est stationnaire (Dickey-Fuller simple, test de racine unitaire)
df1=ur.df(diff,type="none") #La série ne présentant pas de tendance a priori, on considère le cas où la série n’a pas de tendance déterministe. Nous n'incluons pas non plus de constante dans le test.
summary(df1) 
#On rejette très largement l'hypothèse nulle de non stationnarité, la série est donc probablement stationnaire.
df2=ur.df(diff,type="drift") #Nous procédons à test où nous incluons une constante
summary(df2) 
#On rejette toujours l'hypothèse nulle de non stationnarité.
#On fait un test de Dickey Fuller plus compliqué, en intégrant les retards et en cherchant le bon nombre de retards dans la modélisation 
df3=ur.df(diff,type="none",lags=10)
summary(df3)
#Les 5 derniers retards sont significatifs, il faudrait regarder les critères d'information pour sélectionner le meilleur modèle de test.
#On rejette toujours l'hypothèse nulle de non stationnarité.

#On regarde qu'on a pas sur-différencié à l'aide de l'autocorrélogramme inverse 
#Création de la fonction
ACFinverse <- function(z,p=15){
    g<-TacvfMA(GetFitARpLS(z-mean(z),1:p)$phiHat, lag.max=p)
    g/g[1]
  }

AcfPlot(ACFinverse(diff))
#On a tracé l'autocorrélogramme inverse pour notre série 
#Il est décroissant exponentiellement, signe que nous n'avons pas sur-différencié

#Partie 2 : Modèles ARMA (ARIMA dans notre cas)
#Estimation et sélection des modèles par le critère d'information AIC
#On regarde les résultats du critère AIC pour tous les ordres et les modèles retenus sont ceux qui minimisent ce critère.
pmax <-12
qmax <-12
#On crée une matrice de 0 avec en colonnes les ordres 0<=q<=12 et en lignes 0<=p<=12
matriceAIC<- array(0, dim=c(pmax+1,qmax+1), dimnames =list(c(0:pmax), c(0:qmax)))
for (p in 0:pmax){
  for (q in 0:qmax) {
    if (p==11 & q==2 | p==11 & q==10| p==11 & q==11 | p==11 & q==12 ){
      matriceAIC[p+1,q+1] <- NA #car il y a des erreurs dans l'optimisation
    }
    else{
      arma <- arima( diff, order=c(p,0,q),method="CSS-ML", optim.control = list(maxit=2000))
      matriceAIC[p+1,q+1] <- arma$aic
    }
  }
}
matriceAIC
#On sélectionne les modèles ARMA(5,12), ARMA(6,12), ARMA(7,12) , ARMA(10,12) ,ARMA(6,11) et ARMA(9,11), ARMA(6,10) et ARMA(7,10),

#On regarde également le critère BIC, plus pénalisant en terme de complexité du modèle.
pmax <-12
qmax <-12
#On crée une matrice de 0 avec en colonnes les ordres 0<=q<=12 et en lignes 0<=p<=12
matriceBIC<- array(0, dim=c(pmax+1,qmax+1), dimnames =list(c(0:pmax), c(0:qmax)))
for (p in 0:pmax){
  for (q in 0:qmax) {
    if (p==11 & q==2 | p==11 & q==10| p==11 & q==11 | p==11 & q==12 ){
      matriceBIC[p+1,q+1] <- NA #car il y a des erreurs dans l'optimisation
    }
    else{
     bic <- BIC(arima( diff, order=c(p,0,q),method="ML"))
      matriceBIC[p+1,q+1] <- bic
    }
  }
}
matriceBIC
#Les modèles sélectionnés sont parmi les BIC les plus faibles.


#On regarde les coefficients significatifs dans les ARIMA sélectionnées (autrement dit les statistiques de test), pour affiner notre modélisation.
#On commence par le modèle Arma(10,12)
arma1012 = arima(diff, order = c(10,0,12)) 
print(arma1012$coef/sqrt(diag(arma1012$var.coef)))
#On divise la valeur du coefficient par l'écart type du coefficient pour obtenir la statistique de Student
#Le deuxième coefficient de la partie AR est significatif (en valeur absolue plus grand que 1,96) et les autres non, donc on peut réduire l'odre de la partie AR.
#Le douxième élément de la partie MA est très significatif, le 3ème est significatif et les autres non.

#On se penche maintenant sur le modèle Arma(7,12).
arma712 = arima(diff, order = c(7,0,12)) 
print(arma712$coef/sqrt(diag(arma712$var.coef))) 
#Le douzième coefficient de la partie MA est très significatif mais pas le 7ème de la partie AR. Nous pouvons diminuer l'ordre de la partie AR.

#On regarde maintenant le modèle Arma(6,12).
arma612 = arima(diff, order = c(6,0,12)) 
print(arma612$coef/sqrt(diag(arma612$var.coef))) 
#Le 6ème coefficient de la partie AR n'est pas significatif.

#On regarde maintenant le modèle Arma(5,12).
arma512 = arima(diff, order = c(5,0,12)) 
print(arma512$coef/sqrt(diag(arma512$var.coef))) 
#Le 5ème coefficient de la partie AR est significatif, ainsi que le 12ème coefficient de la partie MA.
#On sélectionne ce modèle, surtout dans la mesure où c'est celui qui minimise le critère AIC

#Validité du modèle Arma(5,12) (Tests Portmanteau)
#Test du portemanteau : test de blancheur des résidus (type Ljung-Box)
residual<-residuals(arma512)
plot(residual)
Box.test(residual,type ="Ljung-Box",lag=30) 
#Donc on ne rejette pas H0 (qui est "les résidus sont un bruit blanc"), on passe donc le test Portmanteau.
#On sélectionne ce modèle.

#Test de blancheur des résidus (type Box-Pierce, très proche du test de Ljung-Box)
Box.test(residual,type ="Box-Pierce",lag=30)                            
#Donc on ne rejette pas H0 également (qui est "les résidus sont un bruit blanc").


#Prévision
# On dit à R que l'objet traité est bien une série temporelle ( nécessaire pour utiliser les fonctions forecast et predict)
donnees <- ts(donnees["Indice"] , start=c(1990,1) , end=c(2017,1) , frequency=12)


#5. Equation de l'I.C.

#6. HP: modèle stationnaire + ARIMA vrai moodèle + le bruit est gaussien et matrice de covariance résidus inversible 

#7.
#1ère méthode avec la fonction predict
predict(arima(diff, order = c(5,0,12)),2, se.fit = TRUE)
#On récupère les 2 valeurs de la prédiction: 0.47714514 0.04333077
#On récupère les standards-errors pour les 2 périodes: 0.1446143 0.1726248

#2ème méthode avec la fonction forecast, on vérifie que les 2 valeurs prédites sont bien les mêmes
modele = arima(as.matrix(diff), order = c(5,0,12))
forecast.Arima(modele, h=2)
#On ajoute les deux valeurs de prédiction à la suite de la série
diff2 <- c(diff, 0.47714514,  0.04333077) 


#On trace les 12 dernières valeurs + les deux valeurs obtenues par prédiction
ans = c(2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)
plot(ans,diff2[300:314], type = "l", xlab = "Années", ylim=c(-0.5,0.7), main="Valeur de la série estimée et prévisions")
abline(v=2017, col="blue", lty = 2) #pour marquer le début de la prévision 

#On trace à présent l'intervalle de confiance 
#En segments: (moche)
segments(2018,diff2[313]-0.1473518,2019,diff2[314]- 0.1726248, col="red")
segments(2018,diff2[313]+0.1473518,2019,diff2[314]+ 0.1726248, col="red")
#Par points:
points(2018,diff2[313]-0.1473518, col= "red")
points(2018,diff2[313]+0.1473518, col= "red")
points(2019,diff2[314]- 0.1726248, col= "blue")
points(2019,diff2[314]+ 0.1726248, col="blue")

#8. Causalité au sens de Granger






