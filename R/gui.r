# 
# if (!require("gWidgets")) stop("gWidgets required for GUI")
# if (!require("cairoDevice")) stop("cairoDevice required for GUI")
# 
# win <- gwindow("localmds")
# group <- ggroup(container=win)
# 
# table <-  glayout(container=group)
# table[1,1] <- glabel("k")
# sl1 <- gslider(1, 10)
# table[1,2] <- sl1
# table[2,1] <- glabel("mu")
# table[2,2] <- gslider(0, 10)
# table[3,1] <- glabel("tau")
# table[3,2] <- gslider(0, 10)
# table[4,1] <- glabel("lambda")
# table[4,2] <- gslider(0, 10)
# table[5,1] <- glabel("nu")
# table[5,2] <- gslider(0, 10)
# visible(table) <- TRUE
# 
# ggraphics(container=group)