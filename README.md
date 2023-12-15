Code for paper "A three-state coupled Markov switching model for COVID-19 outbreaks across Quebec based on hospital admissions" by Dirk Douwes-Schultz, Alexandra M. Schmidt, Yannan Shen and David Buckeridge. See paper for more details.

R version 3.6.3

Nimble version 0.9.1
some code may not work on newer versions of Nimble

"Quebec_hospital_study_data.Rdata" - cleaned data

"Full_Coupled_Model.R" - code for the Full Coupled Model from Section 5 of the main text

"No_Absence_Clone_State_Model.R" - code for the No Absence/Clone State Model from Section 5 of the main text 

"Non_Coupled_Model.R" - code for the Non-coupled Model from Section 5 of the main text

"Quebec_hospital_study_data_EE.Rdata" - cleaned data for the Endemic-epidemic
model, with a maximum temporal lag of p=2, from Section 5 of the main text. It is the same as the other data file but moved up a week as the model is autoregressive of order 2.

"Endemic_epidemic.R" - code for the Endemic-epidemic model from Section 5 of the main text with a maximum temporal lag of p=2

"simulation_study_SM.R" - will run a single replication of the simulation study from Section 2 of the supplementary materials

"simulation_study_main_create.R" - code for creating the simulated data set used in the simulation study in Section 4 of the main text.

"simulations1_spatial.Rdata" - the simulated data set used in the simulation study in Section 4 of the main text. Created by "simulation_study_main_create.R".

"simulation_study_main_run.R" - will run the spatial model on the entire simulated data set from Section 4 of the main text (simulations1_spatial.Rdata). Produces and evaluates the retrospective state estimates. To evaluate the real-time state estimates, need to run this code up to times 100,101,...,120 and save the state estimates for the time the model is run up to.

