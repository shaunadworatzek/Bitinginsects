
library(maps)
library(ggmap)

CBAY_mapdata <- CBAY2024_sites %>%
  select(Lon, Lat, ExactSite) %>%
  group_by(ExactSite) %>%
  summarize(
    avg_Lon = mean(Lon, na.rm = TRUE),
    avg_Lat = mean(Lat, na.rm = TRUE)
  ) %>%
  filter(!is.na(avg_Lat), !is.na(avg_Lon)) 
  
  


Coord_images<- CBAY_mapdata %>%
  group_by(avg_Lon, avg_Lat)%>%
  summarise(coord= n())%>%
  filter(!is.na(avg_Lat))%>%
  filter(!is.na(avg_Lon))

world_map <- map_data("world")

map_data_cr <- map_data('world')[map_data('world')$region == "Nunavut",]
p <- ggplot() + coord_fixed(1.3, xlim = c(-105, -103), ylim = c(68, 70))



base_world_messy <- p + 
  geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
               colour="Gainsboro", fill="Gainsboro")+
  geom_polygon(data = map_data_cr, aes(x=long, y=lat, group = group),
               colour = 'red', fill = 'pink')
cleanup <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.background = element_rect(fill = 'white', colour = 'white'), 
                 axis.line = element_line(colour = "white"), legend.position="none",
                 axis.ticks=element_blank(), axis.text.x=element_blank(),
                 axis.text.y=element_blank())
base_world <- base_world_messy + cleanup
base_world+
  geom_point(data=Coord_images, aes(x= avg_Lon, y= avg_Lat), colour="Dim Gray", fill="Dim Gray", size=2)


ll_means <- sapply(googlemapdata_all[1:2], mean)
sq_map2 <- get_map(location = ll_means,  maptype = "satellite", source = "google", zoom = 9)


map <- ggmap(sq_map2) + 
  geom_point(data = CBAY_mapdata, aes(x = avg_Lon, y = avg_Lat, colour = ExactSite), size = 3) +
  scale_color_manual(values = c("#FF7F0E",  # Orange
                                "#2CA02C",  # Green
                                "#D62728",  # Red
                                "#9467BD",  # Purple
                                "#8C564B",  # Brown
                                 "#BCBD22",  # Olive
                                "#17BECF",  # Cyan
                                "#FFD700",  # Gold
                                "#FF69B4",  # Hot Pink
                                "red",  # Turquoise
                                "darkblue",  # Sienna
                                "black",  # Medium Spring Green
                                "#6495ED",  # Cornflower Blue
                                "orange", 
                                "green","maroon","grey", "pink", "purple3" # Crimson
  )) +  
  xlab("Longitude") + 
  ylab("Latitude") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # Remove x-axis label
        axis.text.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = "bold"))   # Remove y-axis label

map

ggsave("mapall.png", map, width = 5, height = 5, dpi = 300)


#from abundence 

CBAY_mapdataabun <- cbay2024_abundence %>%
  select(Long, Lat, Site) %>%
  group_by(Site) %>%
  summarize(
    avg_Lon = mean(Long, na.rm = TRUE),
    avg_Lat = mean(Lat, na.rm = TRUE)
  ) %>%
  filter(!is.na(avg_Lat), !is.na(avg_Lon)) 




Coord_images<- CBAY_mapdataabun %>%
  group_by(avg_Lon, avg_Lat)%>%
  summarise(coord= n())%>%
  filter(!is.na(avg_Lat))%>%
  filter(!is.na(avg_Lon))

world_map <- map_data("world")

map_data_cr <- map_data('world')[map_data('world')$region == "Nunavut",]
p <- ggplot() + coord_fixed(1.3, xlim = c(-105, -103), ylim = c(68, 70))



base_world_messy <- p + 
  geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
               colour="Gainsboro", fill="Gainsboro")+
  geom_polygon(data = map_data_cr, aes(x=long, y=lat, group = group),
               colour = 'red', fill = 'pink')
cleanup <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.background = element_rect(fill = 'white', colour = 'white'), 
                 axis.line = element_line(colour = "white"), legend.position="none",
                 axis.ticks=element_blank(), axis.text.x=element_blank(),
                 axis.text.y=element_blank())
base_world <- base_world_messy + cleanup
base_world+
  geom_point(data=Coord_images, aes(x= avg_Lon, y= avg_Lat), colour="Dim Gray", fill="Dim Gray", size=2)


ll_means <- sapply(googlemapdata_all[1:2], mean)
sq_map2 <- get_map(location = ll_means,  maptype = "satellite", source = "google", zoom = 9)


mapall <- ggmap(sq_map2) + 
  geom_point(data = CBAY_mapdataabun, aes(x = avg_Lon, y = avg_Lat, colour = Site), size = 3) +
  scale_color_manual(values = c("#FF7F0E",  # Orange
                                "#2CA02C",  # Green
                                "#D62728",  # Red
                                "#9467BD",  # Purple
                                "#8C564B",  # Brown
                                "#BCBD22",  # Olive
                                "#17BECF",  # Cyan
                                "#FFD700",  # Gold
                                "#FF69B4",  # Hot Pink
                                "red",  # Turquoise
                                "darkblue",  # Sienna
                                "black",  # Medium Spring Green
                                "#6495ED",  # Cornflower Blue
                                "orange",
                                "pink"# Crimson
  )) +  
  xlab("Longitude") + 
  ylab("Latitude") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # Remove x-axis label
        axis.text.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = "bold"))   # Remove y-axis label

mapall

ggsave("mapall.png", mapall, width = 5, height = 5, dpi = 300)
