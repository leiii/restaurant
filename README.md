# Predicting neighborhoods' socioeconomic attributes using restaurant data

<div align="center">
  <img src="https://github.com/leiii/restaurant/blob/master/food-small.jpg" width = "1000"><br><br>
</div>

Proceedings of the National Academy of Sciences (PNAS)

by Lei Dong (MIT), Carlo Ratti (MIT), and Siqi Zheng (MIT)


### Abstract

Accessing high-resolution, timely socioeconomic data such as data on population, employment, and enterprise activity at the neighborhood level is critical for social scientists and policy makers to design and implement location-based policies. However, in many developing countries or cities, reliable local-scale socioeconomic data remain scarce. Here, we show an easily accessible and timely updated location attribute—restaurant—can be used to accurately predict a range of socioeconomic attributes of urban neighborhoods. We merge restaurant data from an online platform with 3 microdatasets for 9 Chinese cities. Using features extracted from restaurants, we train machine-learning models to estimate daytime and nighttime population, number of firms, and consumption level at various spatial resolutions. The trained model can explain 90 to 95% of the variation of those attributes across neighborhoods in the test dataset. We analyze the tradeoff between accuracy, spatial resolution, and number of training samples, as well as the heterogeneity of the predicted results across different spatial locations, demographics, and firm industries. Finally, we demonstrate the cross-city generality of this method by training the model in one city and then applying it directly to other cities. The transferability of this restaurant model can help bridge data gaps between cities, allowing all cities to enjoy big data and algorithm dividends.

[Web](http://senseable.mit.edu/tasty-data/) | [Paper](https://www.pnas.org/content/116/31/15447) | [Appendix](https://www.pnas.org/content/suppl/2019/07/09/1903064116.DCSupplemental)

### Replication data and code

- restaurant_replication.R
    * R code to replicate the results of the paper
    
- data_dianping
    * Dianping restaurant data of nine cities
    * Baoding, Beijing, Chengdu, Hengyang, Kunming, Shenyang, Shenzhen, Yueyang, and Zhengzhou
    
- rst 
    * Model training results
    * [Download](https://drive.google.com/open?id=1O8rIy6CDWjapFu1YOYmmOqte8WfHe4Q-)
    
- feature 
    * Feature for training
    * [Download](https://drive.google.com/open?id=1VbWKqrkNU6MIZb17xH8B1y-k0PWdHUVw)
    
    
Contact: arch.dongl@gmail.com
