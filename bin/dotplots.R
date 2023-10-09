# Generate dotplots to compare Starships
#===========================================

#colour palette
gene_pal1 <- c("firebrick", "darkgoldenrod2", "dark grey")
gene_pal2 <- c("darkgoldenrod2", "dark grey")

hephaestus_promer <- read.table(file = "/home/aaronvogan/Documents/Grants/ERC2023/hephaestus_protein.out", header = F) # read mummer file
names(hephaestus_promer) <- c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "PER_IDY", "PER_SIM", "PER_STP", "TOT_LEN1", "TOT_LEN2", "FRM1", "FRM2", "TAG1", "TAG2")

hephaestus_promer <- hephaestus_promer %>% mutate(position = seq(1, by = 1, length.out = n()))
hpros <- hephaestus_promer %>% select(S1, S2, position) # Split table into reference coordinates
names(hpros) <- c("Seq1", "Seq2", "Position")
hproe <- hephaestus_promer %>% select(E1, E2, position) # Split table into reference coordinates
names(hproe) <- c("Seq1", "Seq2", "Position")
hpro <- rbind(hpros, hproe)

Voyager_MvB <- read.table(file = "/home/aaronvogan/Documents/Grants/ERC2023/Voyager_MvB.out", header = F) # read mummer file
names(Voyager_MvB) <- c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "PER_IDY", "TOT_LEN1", "TOT_LEN2", "TAG1", "TAG2")

Voyager_MvB <- Voyager_MvB %>% mutate(position = seq(1, by = 1, length.out = n()))
VMBnucs <- Voyager_MvB %>% select(S1, S2, position) # Split table into reference coordinates
names(VMBnucs) <- c("Seq1", "Seq2", "Position")
VMBnuce <- Voyager_MvB %>% select(E1, E2, position) # Split table into reference coordinates
names(VMBnuce) <- c("Seq1", "Seq2", "Position")
VMBnuc <- rbind(VMBnucs, VMBnuce)

Voyager_MvL <- read.table(file = "/home/aaronvogan/Documents/Grants/ERC2023/Voyager_MvL.out", header = F) # read mummer file
names(Voyager_MvL) <- c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "PER_IDY", "TOT_LEN1", "TOT_LEN2", "TAG1", "TAG2")

Voyager_MvL <- Voyager_MvL %>% mutate(position = seq(1, by = 1, length.out = n()))
VMLnucs <- Voyager_MvL %>% select(S1, S2, position) # Split table into reference coordinates
names(VMLnucs) <- c("Seq1", "Seq2", "Position")
VMLnuce <- Voyager_MvL %>% select(E1, E2, position) # Split table into reference coordinates
names(VMLnuce) <- c("Seq1", "Seq2", "Position")
VMLnuc <- rbind(VMLnucs, VMLnuce)

Voyager_MvB_pro <- read.table(file = "/home/aaronvogan/Documents/Grants/ERC2023/Voyager_MvB_pro.out", header = F) # read mummer file
names(Voyager_MvB_pro) <- c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "PER_IDY", "PER_SIM", "PER_STP", "TOT_LEN1", "TOT_LEN2", "FRM1", "FRM2", "TAG1", "TAG2")

Voyager_MvB_pro <- Voyager_MvB_pro %>% mutate(position = seq(1, by = 1, length.out = n()))
VMBpros <- Voyager_MvB_pro %>% select(S1, S2, position) # Split table into reference coordinates
names(VMBpros) <- c("Seq1", "Seq2", "Position")
VMBproe <- Voyager_MvB_pro %>% select(E1, E2, position) # Split table into reference coordinates
names(VMBproe) <- c("Seq1", "Seq2", "Position")
VMBpro <- rbind(VMBpros, VMBproe)

Voyager_MvL_pro <- read.table(file = "/home/aaronvogan/Documents/Grants/ERC2023/Voyager_MvL_pro.out", header = F) # read mummer file
names(Voyager_MvL_pro) <- c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "PER_IDY", "PER_SIM", "PER_STP", "TOT_LEN1", "TOT_LEN2", "FRM1", "FRM2", "TAG1", "TAG2")

Voyager_MvL_pro <- Voyager_MvL_pro %>% mutate(position = seq(1, by = 1, length.out = n()))
VMLpros <- Voyager_MvL_pro %>% select(S1, S2, position) # Split table into reference coordinates
names(VMLpros) <- c("Seq1", "Seq2", "Position")
VMLproe <- Voyager_MvL_pro %>% select(E1, E2, position) # Split table into reference coordinates
names(VMLproe) <- c("Seq1", "Seq2", "Position")
VMLpro <- rbind(VMLpros, VMLproe)

xlim_VoyagerM <- Voyager_MvB[1,8]
ylim_VoyagerL <- Voyager_MvL[1,9]
ylim_VoyagerB <- Voyager_MvB[1,9]

#Hephaestus vs. Aristaeus

#Mummer data
hephaestus_nucmer <- read.table(file = "/home/aaronvogan/Documents/Grants/ERC2023/hephaestus.out", header = F) # read mummer file
names(hephaestus_nucmer) <- c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "PER_IDY", "TOT_LEN1", "TOT_LEN2", "TAG1", "TAG2")

hephaestus_nucmer <- hephaestus_nucmer %>% mutate(position = seq(1, by = 1, length.out = n()))
hnucs <- hephaestus_nucmer %>% select(S1, S2, position) # Split table into reference coordinates
names(hnucs) <- c("Seq1", "Seq2", "Position")
hnuce <- hephaestus_nucmer %>% select(E1, E2, position) # Split table into reference coordinates
names(hnuce) <- c("Seq1", "Seq2", "Position")
hnuc <- rbind(hnucs, hnuce)

xlim_hephaestus <-  hephaestus_nucmer[1,8]
ylim_hephaestus <- hephaestus_nucmer[1,9]

# Segments
# genes
hephaestus_CDS <- read.table(file = "/home/aaronvogan/Starships/Hephaestus.gff", header = T) # read gff
aristaeus_CDS <- read.table(file = "/home/aaronvogan/Starships/Aristaeus.gff", header = T) # read gff


#dotplot


dotplot <- ggplot() +
  geom_segment(data = hephaestus_CDS, aes(x = CDSStart, y = -1, xend = CDSEnd, yend = -1, colour= Gene), size = 12) +
  geom_segment(data = aristaeus_CDS, aes(x = -1, y = CDS_start, xend = -1, yend = CDS_stop, colour= Gene), size = 12) +
  geom_line(data=hnuc, aes(x=Seq1, y=Seq2, group=Position), size = 1.5) +
  geom_point(data=hnuc, aes(x=Seq1, y=Seq2, group=Position), size=3) +
  scale_colour_manual(values = gene_pal1) +
  scale_x_continuous(limits = c(-2,xlim_hephaestus), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-2,ylim_hephaestus), expand = c(0, 0)) +
  xlab("HEPHAESTUS") +
  ylab("ARISTAEUS") +
  theme(text = element_text(size = 20))

#Figure 4 Starfish paper

##Pvar

###mummerdata

Pvar_nucmer <- read.table(file = "/home/aaronvogan/Starships/Dotplots/Ships/Pvar.out", header = F) # read mummer file
names(Pvar_nucmer) <- c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "PER_IDY", "TOT_LEN1", "TOT_LEN2", "TAG1", "TAG2")

Pvar_nucmer <- Pvar_nucmer %>% mutate(position = seq(1, by = 1, length.out = n()))
Pvarnucs <- Pvar_nucmer %>% select(S1, S2, position) # Split table into reference coordinates
names(Pvarnucs) <- c("Seq1", "Seq2", "Position")
Pvarnuce <- Pvar_nucmer %>% select(E1, E2, position) # Split table into reference coordinates
names(Pvarnuce) <- c("Seq1", "Seq2", "Position")
Pvarnuc <- rbind(Pvarnucs, Pvarnuce)

xlim_Pvar <-  Pvar_nucmer[1,8]
ylim_Pvar <- Pvar_nucmer[1,9]

# Segments
# genes
Pvar_genes <- read.table(file = "/home/aaronvogan/Starships/Dotplots/Paevar.gff", header = T) # read gff
Apro_genes <- read.table(file = "/home/aaronvogan/Starships/Dotplots/Asppro.gff", header = T) # read gff

###dotplot
PvarDotplot <- ggplot() +
  geom_segment(data = Pvar_genes, aes(x = CDSStart, y = -700, xend = CDSEnd, yend = -700, colour= Gene), size = 8) +
  geom_segment(data = Apro_genes, aes(x = -1000, y = CDSStart, xend = -1000, yend = CDSEnd, colour= Gene), size = 8) +
  geom_line(data=Pvarnuc, aes(x=Seq1, y=Seq2, group=Position), size = 1.5) +
  geom_point(data=Pvarnuc, aes(x=Seq1, y=Seq2, group=Position), size=3) +
  scale_colour_manual(values = gene_pal2) +
  scale_x_continuous(limits = c(-2000,xlim_Pvar), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1500,ylim_Pvar), expand = c(0, 0)) +
  xlab("paevar1_s08201") +
  ylab("asppro1_s02300") +
  theme(text = element_text(size = 20))  +
  theme(legend.position = "none")

ggsave(file = "PvarDotplot.svg", width = 12, height = 12, scale = 2.5, units = "cm")

##Fusvir

###mummerdata

Fvir_nucmer <- read.table(file = "/home/aaronvogan/Starships/Dotplots/Ships/Fusvir.out", header = F) # read mummer file
names(Fvir_nucmer) <- c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "PER_IDY", "TOT_LEN1", "TOT_LEN2", "TAG1", "TAG2")

Fvir_nucmer <- Fvir_nucmer %>% mutate(position = seq(1, by = 1, length.out = n()))
Fvirnucs <- Fvir_nucmer %>% select(S1, S2, position) # Split table into reference coordinates
names(Fvirnucs) <- c("Seq1", "Seq2", "Position")
Fvirnuce <- Fvir_nucmer %>% select(E1, E2, position) # Split table into reference coordinates
names(Fvirnuce) <- c("Seq1", "Seq2", "Position")
Fvirnuc <- rbind(Fvirnucs, Fvirnuce)

xlim_Fvir <-  Fvir_nucmer[1,8]
ylim_Fvir <- Fvir_nucmer[1,9]

###dotplot
ggplot() +
  geom_line(data=Fvirnuc, aes(x=Seq1, y=Seq2, group=Position), size = 1.5) +
  geom_point(data=Fvirnuc, aes(x=Seq1, y=Seq2, group=Position), size=3)

##Aspfum

###mummerdata

Afum_nucmer <- read.table(file = "/home/aaronvogan/Starships/Dotplots/Ships/Aspfum.out", header = F) # read mummer file
names(Afum_nucmer) <- c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "PER_IDY", "TOT_LEN1", "TOT_LEN2", "TAG1", "TAG2")

Afum_nucmer <- Afum_nucmer %>% mutate(position = seq(1, by = 1, length.out = n()))
Afumnucs <- Afum_nucmer %>% select(S1, S2, position) # Split table into reference coordinates
names(Afumnucs) <- c("Seq1", "Seq2", "Position")
Afumnuce <- Afum_nucmer %>% select(E1, E2, position) # Split table into reference coordinates
names(Afumnuce) <- c("Seq1", "Seq2", "Position")
Afumnuc <- rbind(Afumnucs, Afumnuce)

xlim_Afum <-  Afum_nucmer[1,8]
ylim_Afum <- Afum_nucmer[1,9]

# Segments
# genes
Afum1174_genes <- read.table(file = "/home/aaronvogan/Starships/Dotplots/Aspfum3_s1174.gff", header = T) # read gff
Afum1168_genes <- read.table(file = "/home/aaronvogan/Starships/Dotplots/Aspfum3_s1168.gff", header = T) # read gff

###dotplot
AfumDotplot <- ggplot() +
  geom_segment(data = Afum1174_genes, aes(x = CDSStart, y = -700, xend = CDSEnd, yend = -700, colour= Gene), size = 8) +
  geom_segment(data = Afum1168_genes, aes(x = -700, y = CDSStart, xend = -700, yend = CDSEnd, colour= Gene), size = 8) +
  geom_line(data=Afumnuc, aes(x=Seq1, y=Seq2, group=Position), size = 1.5) +
  geom_point(data=Afumnuc, aes(x=Seq1, y=Seq2, group=Position), size=3) +
  scale_colour_manual(values = gene_pal1) +
  scale_x_continuous(limits = c(-1500,xlim_Afum), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1500,ylim_Afum), expand = c(0, 0)) +
  xlab("aspfum3_s01174") +
  ylab("aspfum3_s01168") +
  theme(text = element_text(size = 20))  +
  theme(
    legend.position = "none")

ggsave(file = "AspfumDotplot.svg", width = 12, height = 12, scale = 2.5, units = "cm")


##Aspasp

###mummerdata

Aasp_nucmer <- read.table(file = "/home/aaronvogan/Starships/Dotplots/Ships/Aspasp.out", header = F) # read mummer file
names(Aasp_nucmer) <- c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "PER_IDY", "TOT_LEN1", "TOT_LEN2", "TAG1", "TAG2")

Aasp_nucmer <- Aasp_nucmer %>% mutate(position = seq(1, by = 1, length.out = n()))
Aaspnucs <- Aasp_nucmer %>% select(S1, S2, position) # Split table into reference coordinates
names(Aaspnucs) <- c("Seq1", "Seq2", "Position")
Aaspnuce <- Aasp_nucmer %>% select(E1, E2, position) # Split table into reference coordinates
names(Aaspnuce) <- c("Seq1", "Seq2", "Position")
Aaspnuc <- rbind(Aaspnucs, Aaspnuce)

xlim_Afum <-  Aasp_nucmer[1,8]
ylim_Afum <- Aasp_nucmer[1,9]

###dotplot
ggplot() +
  geom_line(data=Aaspnuc, aes(x=Seq1, y=Seq2, group=Position), size = 1.5) +
  geom_point(data=Aaspnuc, aes(x=Seq1, y=Seq2, group=Position), size=3)

##Fusoxy

###mummerdata

Foxy_nucmer <- read.table(file = "/home/aaronvogan/Starships/Dotplots/Ships/Fusoxy.out", header = F) # read mummer file
names(Foxy_nucmer) <- c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "PER_IDY", "TOT_LEN1", "TOT_LEN2", "TAG1", "TAG2")

Foxy_nucmer <- Foxy_nucmer %>% mutate(position = seq(1, by = 1, length.out = n()))
Foxynucs <- Foxy_nucmer %>% select(S1, S2, position) # Split table into reference coordinates
names(Foxynucs) <- c("Seq1", "Seq2", "Position")
Foxynuce <- Foxy_nucmer %>% select(E1, E2, position) # Split table into reference coordinates
names(Foxynuce) <- c("Seq1", "Seq2", "Position")
Foxynuc <- rbind(Foxynucs, Foxynuce)

xlim_Foxy<-  Foxy_nucmer[1,8]
ylim_Foxy <- Foxy_nucmer[1,9]

# Segments
# genes
Foxy5919_genes <- read.table(file = "/home/aaronvogan/Starships/Dotplots/Fusoxy21_5919.gff", header = T) # read gff
Foxy5925_genes <- read.table(file = "/home/aaronvogan/Starships/Dotplots/Fusoxy22_5925.gff", header = T) # read gff

###dotplot
FoxyDotplot <- ggplot() +
  geom_segment(data = Foxy5919_genes, aes(x = CDSStart, y = -1500, xend = CDSEnd, yend = -1500, colour= Gene), size = 12) +
  geom_segment(data = Foxy5925_genes, aes(x = -1500, y = CDSStart, xend = -1500, yend = CDSEnd, colour= Gene), size = 12) +
  geom_line(data=Foxynuc, aes(x=Seq1, y=Seq2, group=Position), size = 1.5) +
  geom_point(data=Foxynuc, aes(x=Seq1, y=Seq2, group=Position), size=3) +
  scale_colour_manual(values = gene_pal1) +
  scale_x_continuous(limits = c(-3000,xlim_Foxy), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-3000,ylim_Foxy), expand = c(0, 0)) +
  xlab("Fusoxy21_s05919") +
  ylab("Fusoxy22_s05925") +
  theme(text = element_text(size = 20)) +
  theme(
    legend.position = c(.95, .45),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )

ggsave(file = "FoxyDotplot.svg", width = 12, height = 12, scale = 2.5, units = "cm")

##Aspfla

###mummerdata

Afla_nucmer <- read.table(file = "/home/aaronvogan/Starships/Dotplots/Ships/Afla.out", header = F) # read mummer file
names(Afla_nucmer) <- c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "PER_IDY", "TOT_LEN1", "TOT_LEN2", "TAG1", "TAG2")

Afla_nucmer <- Afla_nucmer %>% mutate(position = seq(1, by = 1, length.out = n()))
Aflanucs <- Afla_nucmer %>% select(S1, S2, position) # Split table into reference coordinates
names(Aflanucs) <- c("Seq1", "Seq2", "Position")
Aflanuce <- Afla_nucmer %>% select(E1, E2, position) # Split table into reference coordinates
names(Aflanuce) <- c("Seq1", "Seq2", "Position")
Aflanuc <- rbind(Aflanucs, Aflanuce)

xlim_Afla<-  Afla_nucmer[1,8]
ylim_Afla <- Afla_nucmer[1,9]

# Segments
# genes
Afla1131_genes <- read.table(file = "/home/aaronvogan/Starships/Dotplots/Aspfla7_s01131.gff", header = T) # read gff
Afla1155_genes <- read.table(file = "/home/aaronvogan/Starships/Dotplots/Aspfla12_s01155.gff", header = T) # read gff

###dotplot
AflaDotplot <- ggplot() +
  geom_segment(data = Afla1131_genes, aes(x = CDSStart, y = -2000, xend = CDSEnd, yend = -2000, colour= Gene), size = 12) +
  geom_segment(data = Afla1155_genes, aes(x = -2000, y = CDSStart, xend = -2000, yend = CDSEnd, colour= Gene), size = 12) +
  geom_line(data=Aflanuc, aes(x=Seq1, y=Seq2, group=Position), size = 1.5) +
  geom_point(data=Aflanuc, aes(x=Seq1, y=Seq2, group=Position), size=3) +
  scale_colour_manual(values = gene_pal1) +
  scale_x_continuous(limits = c(-4000,xlim_Afla), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-4000,ylim_Afla), expand = c(0, 0)) +
  xlab("Aspfla7_s01131") +
  ylab("Aspfla12_s01155") +
  theme(text = element_text(size = 20)) +
  theme(
    legend.position = c(.95, .45),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )

ggsave(file = "AflaDotplot.png", width = 12, height = 12, scale = 2.5, units = "cm")


ggplot(data=hpro, aes(x=Seq1, y=Seq2, group=Position)) +
  geom_line()+
  geom_point() +
  xlim(0,xlim_hephaestus) +
  ylim(0,ylim_hephaestus)

ggplot(data=VMBnuc, aes(x=Seq1, y=Seq2, group=Position)) +
  geom_line()+
  geom_point() +
  xlim(0,xlim_VoyagerM) +
  ylim(0,ylim_VoyagerB)

ggplot(data=VMBpro, aes(x=Seq1, y=Seq2, group=Position)) +
  geom_line()+
  geom_point() +
  xlim(0,xlim_VoyagerM) +
  ylim(0,ylim_VoyagerB)

ggplot(data=VMLnuc, aes(x=Seq1, y=Seq2, group=Position)) +
  geom_line()+
  geom_point() +
  xlim(0,xlim_VoyagerM) +
  ylim(0,ylim_VoyagerL)

ggplot(data=VMLpro, aes(x=Seq1, y=Seq2, group=Position)) +
  geom_line()+
  geom_point() +
  xlim(0,xlim_VoyagerM) +
  ylim(0,ylim_VoyagerL)