Data5$P1=as.numeric(Data5$P1)
Data6$P1=as.numeric(Data6$P1)

tiff('Predator Densities Phase No Evolution.tiff', units="in", width=8, height=8.5, res=1000)
ggplot(subset(Data6, Organism=="1P"), aes(B, P1, width=1)) +
  geom_tile(aes(fill=DiffeP), color="grey90") + coord_cartesian(ylim=c(0.12,0.98), xlim=c(1.10,10.90)) +
  geom_tile_pattern(aes(pattern=State), fill=NA, pattern_fill="black", pattern_color="black", pattern_density=0.1, pattern_spacing=0.02) +
  ylab(expression('A'['1']~'ingestion probability'~'('*p[1]*')')) + xlab(expression('Microplastic density'~'('*10^4~particles~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.y.right=element_text(angle=90)) +
  scale_y_continuous(breaks=seq(0.2,1.0,length.out=5), labels=sprintf(seq(0.2,1.0,length.out=5), fmt="%.1f"),
  sec.axis=sec_axis(~rev(.)*1/1, expression('A'['2,1']~'ingestion probability difference'~'('*p[Delta]*')'), breaks=seq(0.9,0.1,length.out=5), labels=sprintf(seq(0.9,0.1,length.out=5), fmt="%.1f"))) +
  scale_x_discrete(breaks=seq(0,0.05,length.out=6), labels=sprintf(seq(0,5.0,length.out=6), fmt="%.1f")) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.line=element_line(colour="black", linetype="solid", size=0.5)) +
  scale_fill_gradientn(name=expression('Predator density'~'('*rotifers~mL^-1*')'), 
  colours=c("firebrick3","white","royalblue2"), values=rescale(x=c(0,6.4,8),from=c(0,8)), na.value="white", breaks=c(-30,8), limits=c(-30,8)) +
  scale_pattern_manual(values=c("EX"="circle","CRC"="none","EC"="none","SS"="none")) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(plot.title=element_text(face="plain", colour="black", size=22, vjust=2, hjust=0.5)) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1.2) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1.2) +
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf, size=1.4) +
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, size=1.4) +
  guides(fill=guide_colorbar(title.position="top", title.hjust=0.5, ticks.colour=NA, barwidth=15), pattern="none") +
  theme(legend.title=element_text(size=14), legend.text=element_text(size=14)) +
  theme(legend.key.size=unit(1.5,"line"), legend.position="top") +
  facet_wrap(~Organism, ncol=1, nrow=1)
dev.off()

tiff('Predator Densities Phase Evolution.tiff', units="in", width=8, height=8.5, res=1000)
ggplot(subset(Data6, Organism=="2P"), aes(B, P1, width=1)) +
  geom_tile(aes(fill=DiffeP), color="grey90") + coord_cartesian(ylim=c(0.12,0.98), xlim=c(1.10,10.90)) +
  geom_tile_pattern(aes(pattern=State), fill=NA, pattern_fill="black", pattern_color="black", pattern_density=0.1, pattern_spacing=0.02) +
  ylab(expression('A'['1']~'ingestion probability'~'('*p[1]*')')) + xlab(expression('Microplastic density'~'('*10^4~particles~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(breaks=seq(0.2,1.0,length.out=5), labels=sprintf(seq(0.2,1.0,length.out=5), fmt="%.1f"),
  sec.axis=sec_axis(~rev(.)*1/1, expression('A'['2,1']~'ingestion probability difference'~'('*p[Delta]*')'), breaks=seq(0.9,0.1,length.out=5), labels=sprintf(seq(0.9,0.1,length.out=5), fmt="%.1f"))) +
  scale_x_discrete(breaks=seq(0,0.05,length.out=6), labels=sprintf(seq(0,5.0,length.out=6), fmt="%.1f")) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.line=element_line(colour="black", linetype="solid", size=0.5)) +
  scale_fill_gradientn(name=expression('Predator density'~'('*rotifers~mL^-1*')'), 
  colours=c("firebrick3","white","royalblue2"), values=rescale(x=c(0,6.4,8),from=c(0,8)), na.value="white", breaks=c(-30,8), limits=c(-30,8)) +
  scale_pattern_manual(values=c("EX"="circle","CRC"="none","EC"="none","SS"="none")) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(plot.title=element_text(face="plain", colour="black", size=22, vjust=2, hjust=0.5)) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1.2) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1.2) +
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf, size=1.4) +
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, size=1.4) +
  guides(fill=guide_colorbar(title.position="top", title.hjust=0.5, ticks.colour=NA, barwidth=15), pattern="none") +
  theme(legend.title=element_text(size=14), legend.text=element_text(size=14)) +
  theme(legend.key.size=unit(1.5,"line"), legend.position="top") +
  facet_wrap(~Organism, ncol=1, nrow=1)
dev.off()

tiff('Prey Frequencies Phase Evolution.tiff', units="in", width=8, height=8.5, res=1000)
ggplot(subset(Data5, Phase=="Post"), aes(B, P1, width=1)) +
  geom_tile(aes(fill=FA2End), color="grey90") + coord_cartesian(ylim=c(0.12,0.98), xlim=c(1.10,10.90)) +
  geom_tile_pattern(aes(pattern=State), fill=NA, pattern_fill="black", pattern_color="black", pattern_density=0.1, pattern_spacing=0.02) +
  ylab(expression('A'['1']~'ingestion probability'~'('*p[1]*')')) + xlab(expression('Microplastic density'~'('*10^4~particles~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(breaks=seq(0.2,1.0,length.out=5), labels=sprintf(seq(0.2,1.0,length.out=5), fmt="%.1f"),
  sec.axis=sec_axis(~rev(.)*1/1, expression('A'['2,1']~'ingestion probability difference'~'('*p[Delta]*')'), breaks=seq(0.9,0.1,length.out=5), labels=sprintf(seq(0.9,0.1,length.out=5), fmt="%.1f"))) +
  scale_x_discrete(breaks=seq(0,0.05,length.out=6), labels=sprintf(seq(0,5.0,length.out=6), fmt="%.1f")) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.line=element_line(colour="black", linetype="solid", size=0.5)) +
  scale_fill_gradient2(name=expression('Prey frequency'~'(%)'), 
  low="white", mid=alpha("royalblue2",0.1), high=alpha("royalblue2",0.8), na.value="white", breaks=c(0,100), limits=c(0,100)) +
  scale_pattern_manual(values=c("EX"="circle","CRC"="none","EC"="none","SS"="none")) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(plot.title=element_text(face="plain", colour="black", size=22, vjust=2, hjust=0.5)) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1.2) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1.2) +
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf, size=1.4) +
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, size=1.4) +
  guides(fill=guide_colorbar(title.position="top", title.hjust=0.5, ticks.colour=NA, barwidth=15), pattern="none") +
  theme(legend.title=element_text(size=14), legend.text=element_text(size=14)) +
  theme(legend.key.size=unit(1.5,"line"), legend.position="top") +
  facet_wrap(~Phase, ncol=1, nrow=1)
dev.off()
