import numpy as np
import math
from tc_calibration import read_tc_calibration, get_calibration_json

def quaternion2matrix(quaternion):
    # http://stackoverflow.com/questions/1556260/convert-quaternion-rotation-to-rotation-matrix
    # Normalization
    n = 1./np.sqrt(qx**2 + qy**2 + qz**2 + qw**2);
    qx, qy, qz, qw = quaternion * n
    return np.matrix([
        [1. - 2.*qy*qy - 2.*qz*qz, 2.*qx*qy - 2.*qz*qw, 2.*qx*qz + 2.*qy*qw],
        [2.*qx*qy + 2.*qz*qw, 1. - 2.*qx*qx - 2.*qz*qz, 2.*qy*qz - 2.*qx*qw],
        [2.*qx*qz - 2.*qy*qw, 2.*qy*qz + 2.*qx*qw, 1. - 2.*qx*qx - 2.*qy*qy]
    ])



class Calibration(object):
    __slots__ = ['__intrinsic', '__translation', '__rotation', 
				 '__projection', '__distortion']

    def __init__(self, calibration_data):

        focal_length_x, focal_length_y = calibration_data['focal_length']
        principal_point_x, principal_point_y = calibration_data['principal_point']
        skew = float(calibration_data['skew'])
        self.__intrinsic = np.matrix([
                [focal_length_x, math.tan(skew) * focal_length_y, principal_point_x],
                [0, focal_length_y, principal_point_y],
                [0, 0, 1]
        ])

        self.__translation = np.matrix(calibration_data['translation'])
        self.__rotation = quaternion2matrix(calibration_data['rotation'])
        self.__projection = self.__intrinsic * np.hstack((
            self.__rotation,
            -self.__rotation * self.__translation.T
        ))
        self.__distortion = calibration_json['distortion']
        if len(self.__distortion) == 4:
            # Set k3 = 0.
            self.__distortion.append(0.)

    def undistort_points(self, pts):
        k1,k2,p1,p2, k3 = self.__distortion
        cx, cy = (K[0,2], K[1,2])
        fx, fy = (K[0,0], K[1,1])
        _fx = 1./fx
        _fy = 1./fy
        undistort_pts = []
        for pt in pts:
            u, v = pt
            x = (u - cy)*_fx
            y = (v - cx)*_fy
            r = np.sqrt(x**2 + y**2)
            u_undistort = (x * (1+ (k1*r**2) + (k2*r**4) + (k3*r**6))) + 2*p1*x*y + p2*(r**2 + 2*x**2)
            v_undistort = (y * (1+ (k1*r**2) + (k2*r**4) + (k3*r**6))) + 2*p2*y*x + p1*(r**2 + 2*y**2)
            x_undistort = fx*u_undistort+ cx
            y_undistort = fy*v_undistort+ cy
            undistort_pts.append(x_undistort[0], y_undistort[0])
        return np.array(undistort_pts)

    def get_groundplane_proj(self, point):
        # noinspection PyPep8Naming
        A = np.vstack((self.__projection, np.array([0, 0, 1, 0])))
        b = np.array([point[0], point[1], 1, 0])
        buf = np.linalg.solve(A, b)
        return buf[:3] / buf[3]

    def convert_world_to_image(self, point):
        point = np.append(point, [1])
        q = np.squeeze(
            np.array(self.__projection * point[:, np.newaxis])
        )
        # noinspection PyUnresolvedReferences
        return q[:2] / q[2]
    
    def get_world_part(self, feet_im, part_im):
        feet_w = [self.get_groundplane_proj(foot_im) for foot_im in feet_im]
        foot_vector = feet_w[0] - feet_w[1]
        bias = feet_w[1][1]*feet_w[0][0] - feet_w[0][1]*feet_w[1][0]
        A = np.vstack((self.__projection, [[-foot_vector[1],
                                            foot_vector[0], 0, bias]]))
        b = np.append(part_im, [1, 0])
        
        res = np.linalg.solve(A, b)
        return res[:3] / res[3]
    
    def compute_height(self, feet_im, parts_im):
        if len(feet_im) != 2 or len(parts_im) != 4:
			raise 'Wrong points collection'
        parts_w = [self.get_world_part(feet_im, part_im) for part_im in parts_im]
        dist = lambda x, y: np.sqrt(np.sum(np.power(x - y, 2)))
        head_size = dist(parts_w[0], parts_w[1])
        neck_size = dist(parts_w[1], parts_w[2])
        torso_size = dist(parts_w[2], parts_w[3])
        # Floor point - the projection of the last part on Z=0
        floor_pt = parts_w[-1][:]
        floor_pt[2] = 0.
        leg_size = dist(parts_w[-1], floor_pt)
        # Considering ankle points as foot points,
		# we have to add ankle size to total height
        ankle_size = 0.1
        # noinspection PyUnresolvedReferences
        return head_size + neck_size + torso_size + leg_size + ankle_size
    

tc_calibration = read_tc_calibration('TownCentre-calibration.txt')
calibration_json = get_calibration_json(tc_calibration)
calibration = Calibration(calibration_json)
