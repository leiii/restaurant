# Lei Dong
# Email: arch.dongl@gmail.com
library(glmnet)
library(tidyr)
library(dplyr)    
library(ggplot2) 
library(ggpubr)
library(stringr)

setwd("~/downloads/restaurant/replicate")

########### Figure 1 ##############
food <- read.table("data_dianping/food_hengyang_v2.csv", header = T, sep = ",")
category <- data.frame(table(food$mainCategoryNameEn))
category <- category[order(category$Freq), ]

ggbarplot(tail(category,20), x = "Var1", y = "Freq",
          fill = "#4292c6",           # change fill color by cyl
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "desc",          # Sort the value in dscending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          xlab = "Hengyang"
)
#ggsave("hengyang_cnt_v2.pdf", width = 5, height = 3)



########### Figure 2 ##############
### Figure 2ABCD group size
city <- c("beijing", "shenzhen", "chengdu", 
          "shenyang", "zhengzhou", "kunming",
          "baoding", "yueyang", "hengyang")

# read files
r2 <- data.frame()
for (i in city){
  dt <- read.table(paste0("rst/", i, "_lasso_r2_rolling_v2.csv"), header = TRUE, sep = ",")
  r2 <- rbind(r2, dt)
}

# filter, summarise, and plot
for (i in city){
  d <- r2 %>% 
    filter(city == i) %>%
    group_by(cellsize, y) %>% 
    summarise(r2_test_mean = mean(r2_test),
              r2_test_sd = sd(r2_test))
  
  p <- ggplot(d, aes(cellsize/1000, r2_test_mean)) + 
    geom_pointrange(aes(color = factor(y),
                        ymin = r2_test_mean - r2_test_sd, 
                        ymax = r2_test_mean + r2_test_sd), 
                    alpha = 0.8, size = 0.8) + 
    geom_line(aes(color = factor(y)), size = 1.5) + 
    scale_color_manual(name = "legend",
                       breaks = c("employment", "population", "firm", "consumption"),
                       values = c("employment" = "red", "population" = "#E69F00", 
                                  "firm" = "#56B4E9", "consumption" = "#999999")) + 
    geom_vline(xintercept = 4.5, linetype = 2) + 
    annotate("text", label = i, x = 4, y = 0.5) + 
    ylim(c(0.3, 1.1)) + 
    xlab("Grid size (km)") + 
    ylab("Accuracy (R2)") +
    theme_classic(base_size = 18)
  p
  
  ggsave(paste0(i, "_r2_v2.pdf"), width = 6, height = 4)
}



### Figure 2E group size
big_city <- r2 %>% 
  filter(city %in% "beijing" | city %in% "shenzhen" | city %in% "chengdu") %>%
  filter(cellsize == 5000 | cellsize == 1000) %>%
  group_by(cellsize, y) %>% 
  summarise(r2_test_mean = mean(r2_test),
            r2_test_sd = sd(r2_test))

mid_city <- r2 %>%
  filter(city %in% "shenyang" | city %in% "zhengzhou" | city %in% "kunming") %>%
  filter(cellsize == 5000 | cellsize == 1000) %>%
  group_by(cellsize, y) %>% 
  summarise(r2_test_mean = mean(r2_test),
            r2_test_sd = sd(r2_test))

sma_city <- r2 %>% 
  filter(city %in% "baoding" | city %in% "yueyang" | city %in% "hengyang") %>%
  filter(cellsize == 5000 | cellsize == 1000) %>%
  group_by(cellsize, y) %>% 
  summarise(r2_test_mean = mean(r2_test),
            r2_test_sd = sd(r2_test))

# set x axis
sma_city$x <- 0
mid_city$x <- 1
big_city$x <- 2
rst <- as.data.frame(rbind(sma_city, mid_city, big_city))
rst[rst$y == "employment", ]$x <- rst[rst$y == "employment", ]$x - 0.3
rst[rst$y == "population", ]$x <- rst[rst$y == "population", ]$x - 0.15
rst[rst$y == "consumption", ]$x <- rst[rst$y == "consumption", ]$x + 0.15

# plot
ggplot(rst, aes(x, r2_test_mean)) + 
  geom_pointrange(aes(color = factor(y), 
                      ymin = r2_test_mean - r2_test_sd, 
                      ymax = r2_test_mean + r2_test_sd),
                  size = 0.8) + 
  scale_color_manual(name = "legend",
                     breaks = c("employment", "population", "firm", "consumption"),
                     values = c("employment" = "red", "population" = "#E69F00", 
                                "firm" = "#56B4E9", "consumption" = "#999999")) + 
  #ylim(c(0.4, 1.2)) + 
  xlab("City size group") + 
  ylab("Accuracy (R2)") +
  theme_classic(base_size = 18)
#ggsave("citygroup_r2_v2.pdf", width = 6, height = 4, useDingbats=FALSE)



### Figure 2F subsample
bj <- read.table("rst/beijing_subsample_lasso_r2_rolling_v2.csv", header = TRUE, sep = ",")
bj1 <- bj %>% 
  group_by(smp_per, y) %>% 
  summarise(r2_test_mean = mean(r2_test),
            r2_test_sd = sd(r2_test))

ggplot(bj1, aes(smp_per*100, r2_test_mean)) + 
  geom_pointrange(aes(color = factor(y), 
                      ymin = r2_test_mean - r2_test_sd, 
                      ymax = r2_test_mean + r2_test_sd),
                  size = 0.8, alpha = 0.8) + 
  scale_color_manual(name = "legend",
                     breaks = c("employment", "population", "firm", "consumption"),
                     values = c("employment" = "red", "population" = "#E69F00", 
                                "firm" = "#56B4E9", "consumption" = "#999999")) +
  ylim(c(0.6, 1.0)) + 
  xlab("Percentage (%)") + 
  ylab("Accuracy (R2)") +
  theme_classic(base_size = 18)
#ggsave("beijing_subsample_r2_v2.pdf", width = 6, height = 4, useDingbats=FALSE)



########### Figure 3 ##############
### Figure 3A spatial distribution of erros
dt0 <- read.table("feature/beijing_employment_3000_0_0.csv", header = TRUE, sep = ",")
dt1 <- read.table("feature/beijing_employment_3000_0.5_0.csv", header = TRUE, sep = ",")
dt2 <- read.table("feature/beijing_employment_3000_0_0.5.csv", header = TRUE, sep = ",")
dt3 <- read.table("feature/beijing_employment_3000_0.5_0.5.csv", header = TRUE, sep = ",")
dt <- rbind(dt1, dt2, dt3)
dt <- separate(data = dt, col = X, into = c("x", "y"), sep = "_")
dt_ <- separate(data = dt0, col = X, into = c("x", "y"), sep = "_")
dt[is.na(dt)] <- 0
dt_[is.na(dt_)] <- 0

# data format
train_x <- as.matrix(dt[,c(4:ncol(dt))])
train_y <- as.matrix(dt[,c(3)])
test_x <- as.matrix(dt_[,c(4:ncol(dt_))])
test_y <- as.matrix(dt_[,c(3)])

# LASSO fit (alpha=1: LASSO)
cvfit <- cv.glmnet(train_x, train_y, alpha = 1, nfolds = 5)
plot(cvfit)

# r2 of training data
mse <- cvfit$cvm[cvfit$lambda == cvfit$lambda.min]
r2_train <- 1- mse/var(train_y)
print(r2_train)

# r2 of testing data
pred_lasso <- predict(cvfit, newx = test_x, s = "lambda.min")
mse <- sum((pred_lasso - test_y)^2) / length(test_y)
r2_test <- 1- mse/var(test_y)
print(r2_test)

# calculate errors
rst <- data.frame(cbind(dt_$x, dt_$y, dt_$employment, pred_lasso))
colnames(rst) <- c("X", "Y", "observed", "predicted")
rst <- mutate_all(rst, function(x) as.numeric(as.character(x)))
rst$error <- (rst$predicted - rst$observed)/rst$observed
rst$error[rst$error > 1] <- 1
rst$error[rst$error < -1] <- -1

# plot ground truth vs. prediction
ggplot(aes(log10(observed), log10(predicted)), 
       data = rst[rst$predicted > 100 & rst$observed > 100,]) +
  geom_point(color = "blue", alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0) + 
  ggtitle("Beijing daytime population (3km)") +
  xlab("Observed  (log10)") +
  ylab("Predicted (log10)") 
#ggsave("beijing_employment_3000.pdf", width = 4, height = 3)



# plot spatial distribution of errors
library(rgdal) #coord. transform
library(OpenStreetMap) #base map

# convert projection to wgs84
cord_dec <- SpatialPoints(cbind(rst$X, rst$Y), proj4string = CRS("+init=epsg:2436"))
cord_wgs <- spTransform(cord_dec, CRS("+init=epsg:4326"))
lonlat <- as.data.frame(cord_wgs@coords)
colnames(lonlat) <- c("lon", "lat")
rst <- cbind(rst, lonlat)

# boundary of Beijing
LAT1 =  39.45 ; LAT2 = 40.5
LON1 = 115.8 ; LON2 = 117

# boundary of Chengdu
#LAT1 =  30.3; LAT2 = 30.96
#LON1 = 103.73; LON2 = 104.47

# boundary of Zhengzhou
#LAT1 =  34.55; LAT2 = 34.97
#LON1 = 113.46; LON2 = 113.95

# boundary of Baoding
#LAT1 =  38.6; LAT2 = 39.1
#LON1 = 115.17; LON2 = 115.69

map <- openmap(c(LAT2, LON1), c(LAT1, LON2), zoom = NULL,
               type = "stamen-toner",
               mergeTiles = TRUE)
map_latlon <- openproj(map, projection = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# with basemap
osmmap <- autoplot(map_latlon)+
  #labs(x = "Lon", y = "Lat")+
  geom_point(data = rst[rst$observed > 2000,], 
             aes(x = lon, y = lat, color = error*100, size = observed)) +
  scale_size(range = c(1, 5), name = "Daytime popu.") + 
  scale_color_gradient2(low = "steelblue", high = "#f46d43", name = "Error(%)") +
  xlim(LON1, LON2) +
  ylim(LAT1, LAT2)
osmmap
#ggsave("beijing_employment_3000_error_spatial_v2.pdf", width = 5, height = 4, useDingbats=FALSE)



### Figure 3B and Table S3 jiedao demographic
demo <- read.table("rst/beijing_demographic_r2_jiedao.csv", header = TRUE, sep = ",")
demo_feature <- read.table("rst/beijing_demographic_feature_jiedao.csv", header = TRUE, sep = ",")

demo1 <- demo %>% 
  group_by(y) %>% 
  summarise(r2_mean = mean(r2),
            r2_sd = sd(r2))

ggplot(demo1, aes(y, r2_mean)) + 
  geom_pointrange(aes(color = factor(y), 
                      ymin = r2_mean - r2_sd, 
                      ymax = r2_mean + r2_sd),
                  size = 0.8, alpha = 0.8) + 
  xlab("Group") + 
  ylab("Accuracy (R2)") +
  ylim(0.4, 0.8) +
  theme_classic(base_size = 18)
#ggsave("beijing_demographic_r2_jiedao.pdf", width = 4, height = 4, useDingbats=FALSE)



# Table S3
young <- demo_feature[demo_feature$y == "14",]
middle <- demo_feature[demo_feature$y == "15-64",]
old <- demo_feature[demo_feature$y == "65",]
migrant <- demo_feature[demo_feature$y == "migrant",]
wealth <- demo_feature[demo_feature$y == "wealth",]

migrant$split <- as.data.frame(str_split_fixed(migrant$name, "_", 2))[[1]]
wealth$split <- as.data.frame(str_split_fixed(wealth$name, "_", 2))[[1]]
young$split <- as.data.frame(str_split_fixed(young$name, "_", 2))[[1]]
middle$split <- as.data.frame(str_split_fixed(middle$name, "_", 2))[[1]]
old$split <- as.data.frame(str_split_fixed(old$name, "_", 2))[[1]]

migrant_sum <- as.data.frame(table(migrant$split))
wealth_sum <- as.data.frame(table(wealth$split))
young_sum <- as.data.frame(table(young$split))
middle_sum <- as.data.frame(table(middle$split))
old_sum <- as.data.frame(table(old$split))



### Figure 3C Firm
firm <- read.table("rst/beijing_lasso_firm_r2_jiedao_v2.csv", header = TRUE, sep = ",")
firm_feature <- read.table("rst/beijing_lasso_firm_feature_jiedao_v2.csv", header = TRUE, sep = ",")

firm1 <- firm %>% 
  group_by(y) %>% 
  summarise(r2_train_mean = mean(r2_train),
            r2_train_sd = sd(r2_train))

ggplot(firm1, aes(y, r2_train_mean)) + 
  geom_pointrange(aes(color = factor(y), 
                      ymin = r2_train_mean - r2_train_sd, 
                      ymax  =r2_train_mean + r2_train_sd)) + 
  xlab("Percentage (%)") + 
  ylab("Accuracy (R2)") +
  ylim(0.4, 0.6) +
  theme_classic(base_size = 18)
#ggsave("beijing_firm_r2_jiedao.pdf", width = 6, height = 4, useDingbats=FALSE)



# Table S3
business_service <- firm_feature[firm_feature$y == "retail",]
business_service$split <- as.data.frame(str_split_fixed(business_service$name, "_", 2))[[1]]
business_service_sum <- as.data.frame(table(business_service$split))



########### Transfer ##############
### Figure 4
city <- c("beijing", "shenzhen", "chengdu", 
          "shenyang", "zhengzhou", "kunming",
          "baoding", "yueyang", "hengyang")

r2_matrix <- data.frame()
feature <- data.frame()
for (i in city){
  print(i)
  dt <- read.table(paste0("rst/", i, "_lasso_r2_transfer_v2.csv"), header = TRUE, sep = ",")
  fe <- read.table(paste0("rst/" , i, "_lasso_feature_transfer_v2.csv"), header = TRUE, sep = ",")
  feature <- rbind(feature, fe)
  dt1 <- dt %>%
    group_by(city_b, y) %>%
    summarise(r2_test_mean = mean(r2_test),
              r2_test_sd = sd(r2_test))
  dt1 <- as.data.frame(dt1)
  dt1$city_a <- i
  r2_matrix <- rbind(r2_matrix, dt1)
}

# add diag of the matrix
r2_matrix_diag <- data.frame()
for (i in city){
  print(i)
  dt <- read.table(paste0("rst/", i, "_lasso_r2_rolling_v2.csv"), header = TRUE, sep = ",")
  dt1 <- dt %>%
    filter(cellsize == 3000) %>%
    group_by(y) %>%
    summarise(r2_test_mean = mean(r2_test),
              r2_test_sd = sd(r2_test))
  dt1 <- as.data.frame(dt1)
  dt1$city_a <- i
  dt1$city_b <- i
  r2_matrix_diag <- rbind(r2_matrix_diag, dt1)
}

r2_matrix_sum <- rbind(r2_matrix, r2_matrix_diag)
r2_matrix_sum$city_a_order <- factor(r2_matrix_sum$city_a, levels = city)
r2_matrix_sum$city_b_order <- factor(r2_matrix_sum$city_b, levels = city)

# plot
ggplot(data = r2_matrix_sum[r2_matrix_sum$y == "population", ], 
       aes(city_a_order, city_b_order, fill = r2_test_mean)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(r2_test_mean, 2)), alpha = 0.8) +
  scale_fill_gradient(high = "steelblue", low = "white", 
                      limits = c(0.5, 1), name = "R2") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1)) +
  coord_fixed()
#ggsave("transfer_population_r2_v2.pdf", width = 6, height = 5, useDingbats=FALSE)


transfer_feature <- feature[feature$y == "consumption",]
transfer_feature$split <- as.data.frame(str_split_fixed(transfer_feature$name, "_", 2))[[1]]
transfer_feature_sum <- as.data.frame(table(transfer_feature$split))




########### Supplementary Information ##############
### Figure S3-5
city <- c("baoding", "yueyang", "hengyang")
empl <- data.frame()
popu <- data.frame()
firm <- data.frame()
cons <- data.frame()
poi.cnt <- data.frame()
price <- data.frame()

for (i in city){
  empl_ <-  read.table(paste0("feature/", i ,"_employment_3000_0_0.csv"), header = TRUE, sep = ",")
  popu_ <-  read.table(paste0("feature/", i ,"_population_3000_0_0.csv"), header = TRUE, sep = ",")
  firm_ <-  read.table(paste0("feature/", i ,"_firm_3000_0_0.csv"), header = TRUE, sep = ",")
  cons_ <-  read.table(paste0("feature/", i ,"_consumption_3000_0_0.csv"), header = TRUE, sep = ",")
  empl <- rbind(empl, cbind(empl_$employment, rep(i, nrow(empl_))))
  popu <- rbind(popu, cbind(popu_$population, rep(i, nrow(popu_))))
  firm <- rbind(firm, cbind(firm_$firm, rep(i, nrow(firm_))))
  cons <- rbind(cons, cbind(cons_$consumption, rep(i, nrow(cons_))))
  poi.cnt <- rbind(poi.cnt, cbind(empl_$poi_cnt, rep(i, nrow(empl_))))
  price <- rbind(price, cbind(empl_$price_mean, rep(i, nrow(empl_))))
}

ggplot(empl, aes(x=log(as.numeric(as.character(V1))), color=as.factor(V2))) +
  geom_density(size = 1) + 
  scale_color_brewer(palette="Paired") + 
  xlab("Daytime population (log)") + 
  theme_classic(base_size = 16)
#ggsave("smallcity_employment.pdf", width = 5, height = 4, useDingbats=FALSE)

ggplot(popu, aes(x=log(as.numeric(as.character(V1))), color=as.factor(V2))) +
  geom_density(size = 1) + 
  scale_color_brewer(palette="Paired") + 
  xlab("Nightime population (log)") + 
  theme_classic(base_size = 16)
#ggsave("smallcity_population.pdf", width = 5, height = 4, useDingbats=FALSE)

ggplot(firm, aes(x=log(as.numeric(as.character(V1))), color=as.factor(V2))) +
  geom_density(size = 1) + 
  scale_color_brewer(palette="Paired") + 
  xlab("Firm (log)") + 
  theme_classic(base_size = 16)
#ggsave("smallcity_firm.pdf", width = 5, height = 4, useDingbats=FALSE)

ggplot(cons, aes(x=log(as.numeric(as.character(V1))), color=as.factor(V2))) +
  geom_density(size = 1) + 
  scale_color_brewer(palette="Paired") + 
  xlab("Consumption (log)") + 
  theme_classic(base_size = 16)
#ggsave("smallcity_consumption.pdf", width = 5, height = 4, useDingbats=FALSE)

ggplot(poi.cnt, aes(x=as.numeric(as.character(V1)), color=as.factor(V2))) +
  geom_density(size = 1) + 
  scale_color_brewer(palette="Paired") + 
  xlab("Restaurant number (log)") + 
  theme_classic(base_size = 16)
#ggsave("smallcity_restaurant.pdf", width = 5, height = 4, useDingbats=FALSE)

ggplot(price, aes(x=as.numeric(as.character(V1)), color=as.factor(V2))) +
  geom_density(size = 1) + 
  scale_color_brewer(palette="Paired") + 
  xlab("Average price (log)") + 
  theme_classic(base_size = 16)
#ggsave("smallcity_price.pdf", width = 5, height = 4, useDingbats=FALSE)



### Figure S6
city <- c("beijing")
r2 <- data.frame()
for (i in city){
  dt <- read.table(paste0("rst/", i, "_lasso_r2_norolling_v2.csv"), header = TRUE, sep = ",")
  r2 <- rbind(r2, dt)
}

d <- r2 %>% 
  group_by(city, cellsize, y) %>% 
  summarise(r2_test_mean = mean(r2_test),
            r2_test_sd = sd(r2_test))

ggplot(d, aes(cellsize/1000, r2_test_mean)) + 
  geom_pointrange(aes(colour = factor(y),
                      ymin=r2_test_mean-r2_test_sd, 
                      ymax=r2_test_mean+r2_test_sd), 
                  alpha = 0.8, size = 0.8) + 
  geom_line(aes(colour = factor(y)), size = 1.5) + 
  scale_color_manual(name = "legend",
                     breaks = c("employment", "population", "firm", "consumption"),
                     values = c("employment" = "red", "population" = "#E69F00", 
                                "firm" = "#56B4E9", "consumption" = "#999999")) + 
  geom_vline(xintercept = 4.5, linetype = 2) + 
  annotate("text", label = i, x = 4, y = 0.5) + 
  ylim(c(0.3, 1.1)) + 
  xlab("Grid size (km)") + 
  ylab("Accuracy (R2)") +
  theme_classic(base_size = 18)
#ggsave("beijing_r2_norolling_v2.pdf", width = 6, height = 4)



### Table S2
r2_rst <- expand.grid(
  city     = c("beijing"),
  cellsize = seq(1000, 5000, by = 500),
  y        = c("population", "employment", "firm", "consumption"),
  r2       = 0
)

for (i in 1:nrow(r2_rst)){
  print (i)
  # readfile
  basefilename <- paste(r2_rst[i, "city"], r2_rst[i, "y"], r2_rst[i, "cellsize"], sep = "_")
  dt <- read.table(paste0("feature/", basefilename, "_0_0.csv"), header = T, sep = ",")
  dt1 <- dt[,c(1,2)]
  
  light <- read.table(paste0("feature/beijing_light_", r2_rst[i, "cellsize"], "_0_0.csv"), header = TRUE, sep = ",")
  light1 <- light %>%
    group_by(idx) %>%
    summarise(light_sum = sum(light))
  m <- merge(dt1, light1, by.x = "X", by.y = "idx")
  model <- summary(lm(log(m[, 2]) ~ log(m$light_sum)))
  r2_rst[i, "r2"] <- model$r.squared
}


