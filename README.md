# Nowcasting-Earthquake-with-Machine-Learning

This was the Design Project which was done under the supervision of **Dr Sumanta Pasari**, *Department of Mathematics*, BITS Pilani.

The project considers the method of nowcasting to determine the current state of fault system and its current progress through the earthquake cycle in seismically active regions. We have implemented the nowcasting of earthquake on 4 different regions (*Taiwan, New Zealand, Himalayas and Sumatra*) using an unsupervised clustering approach divided into five major steps which are:
- Automated definition and classification of seismicity into candidates of seismic bursts.
- Rejection of outliers
- Selects the members of the group of accepted bursts, which will then be displayed as a time series.
- Apply Exponential Moving Average to bursts to construct burst time series.
- Optimization of the group of possible bursts with a cost function
Detailed implementation and results are present in the file **Final Report**. 
