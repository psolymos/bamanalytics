#habitat
#HAB = habitat without open/sparse/dense modifier
mod_Hab <- list(
    . ~ . + HAB,
    . ~ . + HGT,
    . ~ . + HGT + HGT2,
    . ~ . + TR3,
    . ~ . + HGT + TR3,
    . ~ . + HGT + HGT2 + TR3,
    . ~ . + HAB + HGT,
    . ~ . + HAB + HGT + HGT2,
    . ~ . + HAB + HGT + HGT:isDM + HGT:isNF,
    . ~ . + HAB + HGT + HGT2 + HGT:isDM + HGT:isNF + HGT2:isDM + HGT2:isNF,
    . ~ . + HAB + TR3,
    . ~ . + HAB + HGT + TR3,
    . ~ . + HAB + HGT + HGT2 + TR3,
    . ~ . + HAB + HGT + TR3 + HGT:isDM + HGT:isNF,
    . ~ . + HAB + HGT + HGT2 + TR3 + HGT:isDM + HGT:isNF + HGT2:isDM + HGT2:isNF)
mod_Road <- list(
    . ~ . + ROAD,
    . ~ . + ROAD + ROAD:HAB,
    . ~ . + ROAD + ROAD:TR3)
mod_ARU <- list(
    . ~ . + ARU)
mod_Wet <- list(
    . ~ . + CTI,
    . ~ . + CTI + CTI2,
    . ~ . + SLP,
    . ~ . + SLP + SLP2)
mod_Climate <- list(
    . ~ . + CMIJJA + DD0 + DD5 + EMT + TD,
    . ~ . + CMI + DD0 + DD5 + EMT + TD,
    . ~ . + CMI + CMIJJA + DD5 + MSP + TD,
    . ~ . + CMI + CMIJJA + DD0 + DD5 + MSP)
## all years
mod_DisturbFire <- list(
    ## 10 yrs post disturbance
    . ~ . + LIN,
    . ~ . + POL,
    . ~ . + BRN,
    . ~ . + LIN + POL,
    . ~ . + LIN + BRN,
    . ~ . + BRN + POL,
    . ~ . + LIN + BRN + POL,
    ## years since disturbance
    . ~ . + LIN,
    . ~ . + POL,
    . ~ . + YSF,
    . ~ . + LIN + POL,
    . ~ . + LIN + YSF,
    . ~ . + YSF + POL,
    . ~ . + LIN + YSF + POL)
## >1999 years
mod_DisturbGFWFire <- list(
    ## 10 yrs post
    . ~ . + DTB,
    . ~ . + BRN + LSS,
    . ~ . + LIN + POL,
    . ~ . + LIN + POL + DTB,
    . ~ . + LIN + POL + LSS + BRN,
    ## year since
    . ~ . + YSD,
    . ~ . + YSF + YSL,
    . ~ . + LIN + POL,
    . ~ . + LIN + POL + YSD,
    . ~ . + LIN + POL + YSL + YSF)
mod_HS <- list(
    . ~ . + HSH,
    . ~ . + HSH + HSH2)
mod_Year <- list(
    . ~ . + YR,
    . ~ . + YR + YR:REG)

mods_fire <- list(Hab=mod_Hab, 
    Road=mod_Road, ARU=mod_ARU, Wet=mod_Wet, Clim=mod_Climate,
    Dist=mod_DisturbFire,
    HS=mod_HS,
    Year=mod_Year)
mods_gfw <- list(Hab=mod_Hab, 
    Road=mod_Road, ARU=mod_ARU, Wet=mod_Wet, Clim=mod_Climate,
    Dist=mod_DisturbGFWFire,
    HS=mod_HS,
    Year=mod_Year)
