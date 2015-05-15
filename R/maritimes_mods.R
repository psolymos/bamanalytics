WetA <- list(list( # 1
        .~. + LOC_WET_VEG, # Dominant wetland type within local buffer
        .~. + WET_LENGTH, # Perimeter of wetland + stream length within territory buffer
        .~. + WET_TOTAL)) # Proportion of territory buffer that is wetland

WetB <- list(list( # 1
        .~. + LOC_DTW_PROP, # Prop w/i local buffer with depth to water table classified as wet
        .~. + DTW_PROP, # Prop w/i territory buffer with depth to water table classified as wet
        .~. + DTW_STD)) # Std dev of depth to water table within territory buffer

CoverAB <- list(list( # 2
        .~. + LOC_ltree, # Tree species coverage (CASFRI) within local buffer
        .~. + ltree)) # Tree species coverage (CASFRI) within territory buffer


WetCoverC <- list(list( # 1 & 2
        .~. + LOC_ltree + LOC_DTW_PROP + LOC_ltree:LOC_DTW_PROP, # ltree x wetness w/i local buffer
        .~. + ltree + DTW_PROP + ltree:DTW_PROP)) # ltree x wetness w/i territory buffer

modsX <- list(
    "Struct" = list( # 3
        .~. + LOC_CROWNCL_AV, # Average canopy closure within local buffer
        .~. + CRCL_AV, # Average canopy closure within territory buffer
        .~. + LOC_CROWNCL_STD, # Standard deviation of canopy closure within local buffer
        .~. + CRCL_STD, # Standard deviation of canopy closure within territory buffer
        .~. + HT_AV, # Average canopy height within within local buffer
        .~. + HT_STD), # Standard deviation of canopy height within territory buffer
    "Complex" = list( # 4
        .~. + COMPLEXITY), # Average landscape complexity within territory buffer
    "Disturb" = list( # 5
        .~. + ldist, # Leading disturbances within territory buffer  (CASFRI)
        #.~. + , # Average human footprint index within local buffer -- not extracted !!!
        .~. + HUMAN_FOOTPRINT), # Average human footprint index within territory buffer
    "Road" = list( # 6
        .~. + LOC_ROAD01), # Distance from road
    "Connect" = list( # 7
        .~. + CONNECTEDNESS), # Average connectivity index within territory buffer
    "Protect" = list( # 8
        .~. + PROTECT) # Protected/Unprotected
)
modsA <- c("Wet"=WetA, "Cover"=CoverAB, modsX)
modsB <- c("Wet"=WetB, "Cover"=CoverAB, modsX)
modsC <- c("WetCover"=WetCoverC, modsX)
