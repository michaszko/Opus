# User guide

Only MINERvA_nu is supported.

1. Make two NuWro simulations:
	- of MEC eventswith paramters given in 'MINERvA_nu/dataVanilla/onlyMEC_{SF or LFG}/params.txt' and name it 'only_mec.root'
	- and the rest kinematics with paramters given in 'MINERvA_nu/dataVanilla/allEvents_{SF or LFG}/params.txt' ans name it 'allEvents.root'
2. Run two scripts:
	- 'MINERvA_nu/dataVanilla/onlyMEC_{SF or LFG}/makeMatrix.c' which will produce matrix of MEC events used in later analysis and save it to 'MINERvA_nu/data'
	- 'MINERvA_nu/dataVanilla/allEvents_{SF or LFG}/dataProcessing.c' which will produce result of NuWro in ($p_T,p_L$) used to calulate $\chi^2$ and save it to 'MINERvA_nu/data'
3. Use Mathematica or Python script (C++ does not work)
