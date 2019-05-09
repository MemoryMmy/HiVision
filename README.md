# HiVision

(DEMO) Rapid Visualization of Large-Scale Spatial Vector Data

## Setting

***Tab1. Datasets: Roads, POI and Farmland of Mainland China (10-million-scale)***

| Name           | Type       | Records    | Size                 |
| -------------- | ---------- | ---------- | -------------------- |
| China_Road     | LineString | 21,898,508 | 163,171,928 segments |
| China_POI      | Point      | 20,258,450 | 20,258,450 points    |
| China_Farmland | Polygon    | 10,520,644 | 133,830,561 edges    |

***Tab2.  Demo Environment***

| Item             | Description                                    |
| ---------------- | ---------------------------------------------- |
| CPU              | 4core, Intel(R) Xeon(R) CPU E5-2680 v3@2.50GHz |
| Memory           | 32 GB                                          |
| Operating System | Centos7                                        |



## Application Scenarios

### [Demo](http://www.higis.org.cn:8080/hivision/) (Rapid Visualization of large-scale spatial vector data)

The datasets (see Tab 1) used in the demonstration are provided by map service providers. As the datasets are not open published, the raw datasets are encrypted by adding offsets. The interface of the demonstration is simple to use, choose a dataset and click the Enter button, then the visualizing results will be added to the map in real time. Fig 1-3 show the visualizing results.

![fig1](./figures/f1.jpg)

*Fig1. China_POI*



![fig2](./figures/f2.jpg)

*Fig2. China_Road*



![fig3](./figures/f3.jpg)

*Fig3. China_Farmland*



