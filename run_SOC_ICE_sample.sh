#!/bin/tcsh

#To run use:
#nohup run_test &

	(echo .run SOC-ICE_sample.pro;\
	echo SOC_ICE;\
	echo exit)|\
	idl>&Sample_SOC_ICE_Water_Fairbanks.log&

