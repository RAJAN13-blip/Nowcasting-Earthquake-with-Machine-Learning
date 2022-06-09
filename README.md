# Nowcasting-Earthquake-with-Machine-Learning

This was the Design Project(MATH F376) which was done under the supervision of **Dr Sumanta Pasari**, *Department of Mathematics*, BITS Pilani.

The project considers the method of nowcasting to determine the current state of fault system and its current progress through the earthquake cycle in seismically active regions. We have implemented the nowcasting of earthquake on 4 different regions (*Taiwan, New Zealand, Himalayas and Sumatra*) using an unsupervised clustering approach divided into five major steps which are:
- Automated definition and classification of seismicity into candidates of seismic bursts.
- Rejection of outliers
- Selects the members of the group of accepted bursts, which will then be displayed as a time series.
- Apply Exponential Moving Average to bursts to construct burst time series.
- Optimization of the group of possible bursts with a cost function.

Detailed implementation and analysis of results are present in the [**Final Report**](https://github.com/thelords1007/Nowcasting-Earthquake-with-Machine-Learning/blob/main/FINAL%20REPORT.pdf). 

To run the code, run the notebook [Burst CSV Input](https://github.com/thelords1007/Nowcasting-Earthquake-with-Machine-Learning/blob/main/Burst_CSV_Input.ipynb) after downloading [dop_notebooks](https://github.com/thelords1007/Nowcasting-Earthquake-with-Machine-Learning/tree/main/dop_notebooks) folder.
