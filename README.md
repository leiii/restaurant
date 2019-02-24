# Predicting neighborhoods' socioeconomic attributes using restaurant data

<div align="center">
  <img src="https://github.com/leiii/restaurant/blob/master/food-small.jpg" width = "1000"><br><br>
</div>

### Abstract

Accessing high resolution, timely socioeconomic data such as data on population, employment, and enterprise activity at neighborhood level is critical for social scientists and policy makers to design and implement location-based policies. Yet, in many developing countries or cities, reliable local scale socioeconomic data remain scarce. Here we show an easily accessible and timely updated location attribute - restaurant - can be used to accurately predict a range of socioeconomic attributes of urban neighborhoods. We merge restaurant data from an online platform with three novel micro-datasets for nine Chinese cities. Using features extracted from restaurants, we train machine learning models to estimate daytime and nighttime population, number of firms, and consumption level at various spatial resolutions. The trained model can explain 90%-95% of the variation of those attributes across neighborhoods in the test dataset. We analyze the trade-off between accuracy, spatial resolution, and number of training samples, as well as the heterogeneity of the predicted results across different spatial locations, demographics, and firm industries. Finally, we demonstrate the cross-city generality of this method by training the model in one city and then applying it directly to other cities. The transferability of this restaurant model can help bridge data-gaps between cities, allowing all cities to enjoy big data and algorithm dividends.

### Replicate data and code

- data_dianping
    * Dianping restaurant data of nine cities
    * Baoding, Beijing, Chengdu, Hengyang, Kunming, Shenyang, Shenzhen, Yueyang, and Zhengzhou
    
- rst
    * Model training results
    
- feature
    * Feature for training
    
    
Contact: arch.dongl@gmail.com
