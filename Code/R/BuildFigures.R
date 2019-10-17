library(Calahanlab)

fns <- paste0("Visual Elements/Images/",
              c("Building for Mounting Gas Tank.png",
                "Gas Tank Support.png",
                "Dump bucket bearing and side support.png",
                "Dump Bucket Front View.png")
)

AssemblePanels("Visual Elements/Figures/1.png",
               c(2,1,1),
               fns,
               6.5,
               1/16,
               300,
               LETTERS[1:4],
               c("white", "black", "black", "black"),
               xoff = 35,
               yoff = -35,
               cex = 1)

fns <- paste0("Visual Elements/Images/",
              c("Module Frame (0.25 m).png",
                "Module Frame (0.5 m).png",
                "Headworks Tray (0.25 m).png",
                "Headworks Tray (0.5 m).png",
                "Floway Tray (0.25 m).png",
                "Floway Tray (0.5 m).png")
)

AssemblePanels("Visual Elements/Figures/2.png",
               c(2,2,2),
               fns,
               6.5,
               1/16,
               300,
               LETTERS[1:4],
               c("black", "black", "black", "black", "black", "black"),
               xoff = 35,
               yoff = -35,
               cex = 1)