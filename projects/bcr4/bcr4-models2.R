mods <- list(
    "hab"=list(
        . ~ . + hab,
        . ~ . + hab + ROAD,
        . ~ . + hab + ROAD + ROAD:isForest
    ),
    "ysd"=list(
        . ~ . + ysd:isForest,
        . ~ . + ysd:isForest + ysd2:isForest,
        . ~ . + ysd:isForest + ysd2:isForest + ysd3:isForest,
        . ~ . + ysd05:isForest,
        . ~ . + ysd05:isForest + ysd:isForest,
        . ~ . + ysd05:isForest + ysd:isForest + ysd2:isForest
    ),
    "terrain"=list(
        . ~ . + slp,
        . ~ . + slp + slp2,
        . ~ . + slp + slp2 + slp3
    ),
    "clim"=list(
        . ~ . + CMIJJA + DD0 + DD5 + EMT + MSP + DD02 + DD52 + CMIJJA2 +
            CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP,
        . ~ . + CMI + DD0 + DD5 + EMT + MSP + DD02 + DD52 + CMI2 +
            CMI:DD0 + CMI:DD5 + EMT:MSP,
        . ~ . + CMI + CMIJJA + DD0 + MSP + TD + DD02 + CMI2 + CMIJJA2 +
            CMI:DD0 + CMIJJA:DD0 + MSP:TD,
        . ~ . + CMI + CMIJJA + DD5 + MSP + TD + DD52 + CMI2 + CMIJJA2 +
            CMI:DD5 + CMIJJA:DD5 + MSP:TD,
        . ~ . + CMIJJA + DD0 + DD5 + EMT + TD + MSP + DD02 + DD52 + CMIJJA2 +
            CMIJJA:DD0 + CMIJJA:DD5 + MSP:TD + MSP:EMT,
        . ~ . + CMI + DD0 + DD5 + EMT + TD + MSP + DD02 + DD52 + CMI2 +
            CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT,
        . ~ . + CMI + CMIJJA + DD0 + MSP + TD + EMT + DD02 + CMI2 + CMIJJA2 +
            CMI:DD0 + CMIJJA:DD0 + MSP:TD + MSP:EMT,
        . ~ . + CMI + CMIJJA + DD5 + MSP + TD + EMT + DD52 + CMI2 + CMIJJA2 +
            CMI:DD5 + CMIJJA:DD5 + MSP:TD + MSP:EMT
        )
)
