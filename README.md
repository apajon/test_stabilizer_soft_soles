# test_stabilizer_soft_soles

WARNING : work in progress

Developped for Matlab 2012a and more recent

Need Simulink toolbox.

Soft soles stabilizer using a WPG based on deboor bspline and a deformation estimator to take into account soft sole deformation detailed in https://www.researchgate.net/publication/318412580_Optimized_Humanoid_Walking_with_Soft_Soles


In main.m is the WPG built-in as a programming object in matlab. Just RUN => open any WPG_results_OK_* to test stabilizer


matlab compiler option for the deformation estimator:

mex '-IC:\Users\Giovanni\Documents\eigen-eigen-1306d75b4a21\eigen-eigen-1306d75b4a21' Simulation_3D\splineBasis.cpp

mex '-IC:\Users\Giovanni\Documents\eigen-eigen-1306d75b4a21\eigen-eigen-1306d75b4a21' GaussFtotZMP.cpp Gausstot.cpp

In \fsqp : 
mex mexfsqp.c cfsqp.c qld.c