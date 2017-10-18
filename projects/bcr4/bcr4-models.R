mods <- list(
    "hab"=list(
        . ~ . + hab,
        . ~ . + hab + ROAD,
        . ~ . + hab + ROAD + ROAD:isForest
    ),
    "ysd"=list(
        . ~ . + ysd,
        . ~ . + ysd + ysd2,
        . ~ . + ysd05,
        . ~ . + ysd05 + ysd,
        . ~ . + ysd + ysd:isForest,
        . ~ . + ysd + ysd2 + ysd:isForest + ysd2:isForest,
        . ~ . + ysd05 + ysd05:isForest,
        . ~ . + ysd05 + ysd + ysd05:isForest + ysd:isForest,
#        . ~ . + ysd + ysd:isConif + ysd:isDecid,
#        . ~ . + ysd + ysd2 + ysd:isConif + ysd2:isConif + ysd:isDecid + ysd2:isDecid,
#        . ~ . + ysd05 + ysd05:isConif + ysd05:isDecid,
#        . ~ . + ysd05 + ysd + ysd05:isConif + ysd:isConif + ysd05:isDecid + ysd:isDecid,
        . ~ . + ysd10,
        . ~ . + ysd60,
        . ~ . + ysd100,
        . ~ . + ysd10 + ysdmid,
        . ~ . + ysd10 + ysd10:isForest,
        . ~ . + ysd60 + ysd60:isForest,
        . ~ . + ysd100 + ysd100:isForest,
        . ~ . + ysd10 + ysdmid + ysd10:isForest + ysdmid:isForest
#        . ~ . + ysd10 + ysd10:isConif + ysd10:isDecid,
#        . ~ . + ysd60 + ysd60:isConif + ysd60:isDecid,
#        . ~ . + ysd100 + ysd100:isConif + ysd100:isDecid,
#        . ~ . + ysd10 + ysdmid + ysd10:isConif + ysd10:isDecid + ysdmid:isConif + ysdmid:isDecid
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
