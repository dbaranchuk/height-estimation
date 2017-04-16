global lab_calibration = {
    'intrinsic' :  [[ 1.4141663290265992e+03, 0., 9.4543373843712743e+02],
                    [ 0., 1.4186903309255270e+03, 4.9882814117191981e+02],
                    [ 0., 0., 1.]],
    'skew' : 0,
    'rotation' : [[ 8.4205874916034151e-01, -5.3904494260258895e-01, 1.9173231785290568e-02],
                 [ 2.4392531329014799e-01, 4.1226573035971098e-01, 8.7780260258629728e-01],
                 [ -4.8107951993285275e-01, -7.3448452497349836e-01, 4.7863867173017927e-01]],
    'translation' : [ 1.0079020121414087e+03, 7.4787263478618445e+02, 3.7828233199012411e+03],
    'distortion' : [ 5.9773280266265704e-03, -6.8812400681584937e-02, -3.2746095516290839e-03,
                    -3.4638808119244616e-03, 3.4946717700923085e-03 ]
}

def read_tc_calibration(filename):
	tc_calibration = {}
	f = open(filename, 'r')
	for line in f:
		splitted_line = line.split(' ')
		param_name = splitted_line[0]
		param_value = float(splitted_line[2])
		tc_calibration[param_name] = param_value 
	f.close()
	return tc_calibration


def get_calibration_json(tc_calibration):
	return {
		'focal_length' : [tc_calibration['FocalLengthX'],
						  tc_calibration['FocalLengthY']],
		'principal_point' : [tc_calibration['PrincipalPointX'],
							 tc_calibration['PrincipalPointY']],
		'skew' : tc_calibration['Skew'],
		'translation' : [tc_calibration['TranslationX'],
						 tc_calibration['TranslationY'],
						 tc_calibration['TranslationZ']],
		'rotation' : [tc_calibration['RotationX'],
					  tc_calibration['RotationY'],
					  tc_calibration['RotationZ'],
					  tc_calibration['RotationW']],
		'distortion' : [tc_calibration['DistortionK1'],
						tc_calibration['DistortionK2'],
						tc_calibration['DistortionP1'],
						tc_calibration['DistortionP2']]
	}
