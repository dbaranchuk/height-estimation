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
		'focal_length' : (tc_calibration['FocalLengthX'],
						  tc_calibration['FocalLengthY']),
		'principal_point' : (tc_calibration['PrincipalPointX'],
							 tc_calibration['PrincipalPointY']),
		'skew' : tc_calibration['Skew'],
		'translation' : (tc_calibration['TranslationX'],
						 tc_calibration['TranslationY'],
						 tc_calibration['TranslationZ']),
		'rotation' : (tc_calibration['RotationX'],
					  tc_calibration['RotationY'],
					  tc_calibration['RotationZ'],
					  tc_calibration['RotationW']),
		'distortion' : (tc_calibration['DistortionK1'],
						tc_calibration['DistortionK2'],
						tc_calibration['DistortionP1'],
						tc_calibration['DistortionP2'])
	}
