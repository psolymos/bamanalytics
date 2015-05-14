WetCoverA <- list( # 1
        .~. + ,
        .~. + ,
        .~. + 
    )
WetCoverB <- list( # 1
        .~. + ,
        .~. + ,
        .~. + 
    )
WetCoverC <- list( # 1
        .~. + ,
        .~. + 
    )
modsX <- list(
    "Struct" = list( # 2
        .~. + ,
        .~. + ,
        .~. + ,
        .~. + ,
        .~. + ,
        .~. + 
    ),
    "Complex" = list( # 3
        .~. + 
    ),
    "Disturb" = list( # 4
        .~. + ,
        .~. + ,
        .~. + 
    ),
    "Road" = list( # 5
        .~. + 
    ),
    "Connect" = list( # 6
        .~. + 
    ),
    "Protect" = list( # 7
        .~. + 
    )
)
modsA <- list("Wet"=WetCoverA, modsX)
modsB <- list("Wet"=WetCoverB, modsX)
modsC <- list("Wet"=WetCoverC, modsX)
