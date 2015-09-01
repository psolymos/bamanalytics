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
    . ~ . + CMIJJA + DD0 + DD5 + EMT + TD + DD02 + DD52,
    . ~ . + CMI + DD0 + DD5 + EMT + TD + DD02 + DD52,
    . ~ . + CMI + CMIJJA + DD0 + MSP + TD + DD02,
    . ~ . + CMI + CMIJJA + DD5 + MSP + TD + DD52,
    . ~ . + CMIJJA + DD0 + DD5 + EMT + TD + MSP + DD02 + DD52,
    . ~ . + CMI + DD0 + DD5 + EMT + TD + MSP + DD02 + DD52,
    . ~ . + CMI + CMIJJA + DD0 + MSP + TD + EMT + DD02,
    . ~ . + CMI + CMIJJA + DD5 + MSP + TD + EMT + DD52)
## all years
mod_DisturbFire <- list(
    . ~ . + BRN,
    . ~ . + YSF,
    . ~ . + LIN + POL,
    . ~ . + LIN + POL + BRN,
    . ~ . + LIN + POL + YSF)
## 2001-2013
mod_DisturbGFWFire <- list(
    . ~ . + DTB,
    . ~ . + YSD,
    . ~ . + BRN + LSS,
    . ~ . + YSF + YSL,
    . ~ . + LIN + POL,
    . ~ . + LIN + POL + DTB,
    . ~ . + LIN + POL + LSS + BRN,
    . ~ . + LIN + POL + YSD,
    . ~ . + LIN + POL + YSL + YSF)
mod_HS <- list(
    . ~ . + HSH,
    . ~ . + HSH + HSH2)
mod_ND <- list(
    . ~ . + ND2)
mod_Year <- list(
    . ~ . + YR,
    . ~ . + YR + YR:EW)

mods_fire <- list(
    Hab=mod_Hab, 
    Road=mod_Road, 
    ARU=mod_ARU, 
    Wet=mod_Wet, 
    Dist=mod_DisturbFire,
    Clim=mod_Climate,
    Year=mod_Year)
mods_gfw <- list(
    Hab=mod_Hab, 
    Road=mod_Road, 
    ARU=mod_ARU, 
    Wet=mod_Wet, 
    Dist=mod_DisturbGFWFire,
    Clim=mod_Climate,
    Year=mod_Year)
