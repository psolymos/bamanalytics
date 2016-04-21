mod_Hab <- list(
    . ~ . + HAB,
    . ~ . + HABTR)
mod_Road <- list(
    . ~ . + ROAD,
    . ~ . + ROAD + ROAD:isNF,
    . ~ . + ROAD + ROAD:isDev + ROAD:isWet + ROAD:isOpn)
mod_Hgt <- list(
    . ~ . + HGT,
    . ~ . + HGT + HGT2,
    . ~ . + HGT + HGT:isDM + HGT:isWet,
    . ~ . + HGT + HGT:isDec + HGT:isMix + HGT:isWet,
    . ~ . + HGT + HGT2 + HGT:isDM + HGT:isWet + HGT2:isDM + HGT2:isWet,
    . ~ . + HGT + HGT2 + HGT:isDec + HGT:isMix + HGT:isWet + HGT2:isDec + HGT2:isMix + HGT2:isWet,
    . ~ . + HGT05,
    . ~ . + HGT05 + HGT,
    . ~ . + HGT05 + HGT05:isDM + HGT05:isWet,
    . ~ . + HGT05 + HGT05:isDec + HGT05:isMix + HGT05:isWet,
    . ~ . + HGT05 + HGT + HGT05:isDM + HGT05:isWet + HGT:isDM + HGT:isWet,
    . ~ . + HGT05 + HGT + HGT05:isDec + HGT05:isMix + HGT05:isWet + 
            HGT:isDec + HGT:isMix + HGT:isWet)
## 2001-2013
mod_Dist <- list(
    . ~ . + DTB,
    . ~ . + BRN + LSS,
    . ~ . + YSD,
    . ~ . + YSF + YSL,
    . ~ . + LIN + POL,
    . ~ . + LIN + POL + DTB,
    . ~ . + LIN + POL + BRN + LSS,
    . ~ . + LIN + POL + YSD,
    . ~ . + LIN + POL + YSF + YSL)
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
mod_Year <- list(
    . ~ . + YR)

mods <- list(
    Hab=mod_Hab, 
    Road=mod_Road, 
    Hgt=mod_Hgt, 
    Dist=mod_Dist,
    Wet=mod_Wet, 
    Clim=mod_Climate,
    Year=mod_Year)
