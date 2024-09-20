# Mobile-Human-Ad-Hoc-Networks (MoHANETs)

## **Overview**

This repository contains MATLAB codes developed for the paper "Mobile human ad hoc networks: A communication engineering viewpoint on interhuman airborne pathogen transmission" published in Nano Communication Networks. The research presents a novel communication engineering approach to model the transmission of airborne pathogens, combining various disciplines such as epidemiology, biology, medicine, and fluid dynamics. It introduces the concept of Mobile Human Ad Hoc Networks (MoHANETs), where pathogen-laden droplets are modeled using molecular communication principles. The codes in this repository correspond to the proof-of-concept results, validated with empirical COVID-19 data from Italy. The MATLAB code file names correspond to the figures in the paper (Figures 6-9).

## **Background**

Researchers from various disciplines, such as fluid dynamics, biology, medicine, and epidemiology, work independently to model airborne pathogen transmission and the behavior of epidemics. Fluid dynamics focuses on the propagation of pathogen-laden droplets and their interaction with the air, while biologists study the survival of airborne pathogens and their interaction with human cells. Medical research primarily focuses on finding drugs to treat infectious diseases at the cellular level. Epidemiologists rely on empirical data to develop mathematical models for the spread of epidemics.

However, these epidemiological models often do not integrate detailed information from fluid dynamics, biology, and medicine, resulting in estimations based purely on statistical data. A more comprehensive model would consider factors like droplet dynamics, the air distribution in indoor environments, human-pathogen interactions, drug efficacy, and the locations and mobility of humans.

To address this, we propose a unified framework called MoHANET (Mobile Human Ad Hoc Network) that combines insights from multiple fields with a communication engineering perspective. In this model, mobile human groups are viewed similarly to mobile ad hoc networks (MANETs), where the infectious human is treated as the transmitter (TX), the susceptible human as the receiver (RX), and pathogen-laden droplets act as information carriers.

Molecular communication (MC) enables this framework, as it uses chemical signals instead of electrical ones, providing biocompatibility with the human body and multiscale applicability. This communication engineering approach partitions MoHANET into layers, with each layer representing a research area such as fluid dynamics, biology, medicine, or epidemiology. Similar to conventional telecommunications networks, each layer feeds information to the next, resulting in a more accurate model of disease transmission.

The MoHANET framework allows the use of communication theory to model the complexity of airborne pathogen transmission, offering a new perspective for researchers across different disciplines. This repository provides the codes used to develop and test the framework, along with a proof-of-concept study using the Omnidirectional Multicast Transmission (OMT) algorithm. This algorithm models human mobility with a truncated Lévy walk and estimates the number of effective contacts and infection states based on the relative distance between humans. The infection probability is then used to compute the effective contact rate, which is incorporated into an SIR (susceptible-infectious-recovered) epidemiological model to estimate the time course of an epidemic.

Numerical results, validated by empirical COVID-19 data, demonstrate that MoHANET can accurately predict the number of infected individuals over time by accounting for pathogen-laden droplet transmission and human mobility. The results also show that strengthening population immunity or reducing the number of received droplets can lead to a milder outbreak. The empirical COVID-19 data used for validation is downloaded from:

Dong, E., Du, H., & Gardner, L. (2020). "An interactive web-based dashboard to track COVID-19 in real-time." The Lancet Infectious Diseases, 20(5), 533-534. https://doi.org/10.1016/S1473-3099(20)30120-1

## **Code Overview**

The MATLAB code file names correspond to the figures in the paper (Figures 6-9). The file infection2.m is used for Omnidirectional Multicast Transmission algorithm. Also, data files ara given in .mat format.

## **Citation Requirement**

If you use or build upon this code in your research, or if this code is used for any academic work, publication, or research, proper attribution and citation of the paper is **required**. Please cite the paper below  in any related publications, presentations, or derived works.

Gulec, F., Atakan, B., & Dressler, F. (2022). "Mobile human ad hoc networks: A communication engineering viewpoint on interhuman airborne pathogen transmission." Nano Communication Networks, 32, 100410. https://doi.org/10.1016/j.nancom.2022.100410

For the empirical COVID-19 data, please cite:

Dong, E., Du, H., & Gardner, L. (2020). "An interactive web-based dashboard to track COVID-19 in real time." The Lancet Infectious Diseases, 20(5), 533-534. https://doi.org/10.1016/S1473-3099(20)30120-1

