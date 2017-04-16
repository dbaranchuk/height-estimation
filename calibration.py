import numpy as np
import math
from tc_calibration import read_tc_calibration, get_calibration_json
import os
import json
from PIL import Image, ImageDraw, ImageFont, ImageOps
import time

DATA_DIR = 'towncentre'
OUTPUT_DIR = 'result'


class Calibration(object):
    __slots__ = ['__intrinsic', '__translation', '__C',
                 '__rotation', '__distortion', '__projection']

    def __init__(self, calibration_data):
        fx, fy = calibration_data['focal_length']
        cx, cy = calibration_data['principal_point']
        skew = float(calibration_data['skew'])
        self.__intrinsic = np.matrix([
                [fx, math.tan(skew)*fy, cx],
                [0, fy, cy],
                [0, 0, 1]
        ])
        self.__translation = np.matrix(calibration_data['translation'])
        # Rotation Matrix
        rotation = calibration_data['rotation']
        if len(rotation) == 4:
            self.__rotation = self.quaternion2matrix(rotation)
        else:
            self.__rotation = np.matrix(rotation)
        # Projection Matrix
        self.__projection = self.__intrinsic * np.hstack((
            self.__rotation,
            self.__translation.T
        ))
        # Distortion K1, K2, P1, P2, K3
        self.__distortion = calibration_data['distortion']
        if len(self.__distortion) == 4:
            # Set k3 = 0.
            self.__distortion.append(0.)

    def undistort_points(self, pts):
        k1,k2,p1,p2, k3 = self.__distortion
        K = self.__intrinsic
        cx, cy = (K[0,2], K[1,2])
        fx, fy = (K[0,0], K[1,1])
        _fx = 1./fx
        _fy = 1./fy
        undistort_pts = []
        for pt in pts:
            u, v = pt
            x = (u - cx) * _fx
            y = (v - cy) * _fy
            r = np.sqrt(x**2 + y**2)
            u_undistort = (x * (1+ (k1*r**2) + (k2*r**4) + (k3*r**6))) + 2*p1*x*y + p2*(r**2 + 2*x**2)
            v_undistort = (y * (1+ (k1*r**2) + (k2*r**4) + (k3*r**6))) + 2*p2*y*x + p1*(r**2 + 2*y**2)
            x_undistort = fx*u_undistort + cx
            y_undistort = fy*v_undistort + cy
            undistort_pts.append([x_undistort, y_undistort])
        return np.array(undistort_pts)

    def get_groundplane_proj(self, point):
        # noinspection PyPep8Naming
        A = np.vstack((self.__projection,
                       [0, 0, 1, 0]))
        b = np.array([point[0], point[1], 1, 0])
        buf = np.linalg.solve(A, b)
        self.__C = buf[3]
        return buf[:3] / buf[3]

    def get_world_part(self, part_im):
        A = np.vstack((self.__projection,
                       [0, 0, 0, 1]))
        b = np.append(part_im, [1, self.__C])
        res = np.linalg.solve(A, b)
        return res[:3] / res[3]

    def compute_height(self, feet_im, parts_im):
        if len(feet_im) != 2 or len(parts_im) != 4:
			raise 'Wrong points collection'
        foot_im = ((feet_im[0] + feet_im[1])/2).tolist()
        feet_w = [self.get_groundplane_proj(foot_im)
                  for foot_im in feet_im]
        foot_w = self.get_groundplane_proj(foot_im)
        parts_w = [self.get_world_part(part_im)
                   for part_im in parts_im]

        dist = lambda x, y: np.sqrt(np.sum(np.power(x - y, 2)))
        head_size = dist(parts_w[0], parts_w[1])
        neck_size = dist(parts_w[1], parts_w[2])
        torso_size = dist(parts_w[2], parts_w[3])

        #if np.mean(output[pairs[i,0]]) > activThresh and np.mean(output[pairs[i,1]]) > activThresh:
        MAX_LEG_SIZE = 1.0
        MIN_LEG_SIZE = 0.5
        leg_sizes = np.zeros(3)
        leg_sizes[0] = dist(parts_w[3], feet_w[0])
        leg_sizes[1] = dist(parts_w[3], feet_w[1])
        leg_sizes[2] = dist(parts_w[3], foot_w)
        leg_sizes = leg_sizes[leg_sizes < MAX_LEG_SIZE]
        leg_sizes = leg_sizes[leg_sizes > MIN_LEG_SIZE]
        leg_size = np.mean(leg_sizes) if len(leg_sizes)>0 else (MAX_LEG_SIZE+MIN_LEG_SIZE)/2
        # Considering ankle points as foot points,
		# we have to add ankle size to total height
        ankle_size = 0.2
        # noinspection PyUnresolvedReferences
        return head_size + neck_size + torso_size + leg_size + ankle_size

    def _world2camera(self, X_w):
        R = self.__rotation
        T = self.__translation
        X_w = np.matrix(X_w)
        X_c = R*X_w.T - R*T.T
        return X_c

    def _camera2world(self, X_c):
        R = self.__rotation
        T = self.__translation
        X_c = np.matrix(X_c)
        X_w = R.T*X_c.T + T
        return X_w[0]

    @staticmethod
    def quaternion2matrix(quaternion):
        # http://stackoverflow.com/questions/1556260/convert-quaternion-rotation-to-rotation-matrix
        # Normalization
        quaternion = np.array(quaternion)
        n = 1./np.sqrt(np.sum(quaternion**2))
        qx, qy, qz, qw = quaternion * n
        return np.matrix([
            [1. - 2.*qy*qy - 2.*qz*qz, 2.*qx*qy - 2.*qz*qw, 2.*qx*qz + 2.*qy*qw],
            [2.*qx*qy + 2.*qz*qw, 1. - 2.*qx*qx - 2.*qz*qz, 2.*qy*qz - 2.*qx*qw],
            [2.*qx*qz - 2.*qy*qw, 2.*qy*qz + 2.*qx*qw, 1. - 2.*qx*qx - 2.*qy*qy]
        ])

    def convert_world_to_image(self, point):
        point = np.append(point, [1])
        q = np.squeeze(
            np.array(self.__projection * point[:, np.newaxis])
        )
        # noinspection PyUnresolvedReferences
        return q[:2] / q[2]


if __name__ == '__main__':
    tc_calibration = read_tc_calibration(os.path.join(DATA_DIR, 'TownCentre-calibration.txt'))
    calibration_data = get_calibration_json(tc_calibration)
    calibration = Calibration(calibration_data)

    with open(os.path.join(DATA_DIR, 'annotation.json'), 'r') as f:
        annotation = json.load(f)

    pairs = np.array([[1,2], [2,3], [3,7], [4,5], [4,7], [5,6], [7,9], [9,10],
                      [14,9], [11,12], [12,13], [13,9], [14,15], [15,16]])-1
    start_time = time.time()
    for image_name in annotation.keys():
        pil_image = Image.open(os.path.join(DATA_DIR, image_name))
        draw = ImageDraw.Draw(pil_image)
        for person in annotation[image_name]:
            pts = np.array(person['points'])
            pts = calibration.undistort_points(pts)
            feet_im = [pts[0], pts[5]]
            parts_im = pts[6:10][::-1].tolist()
            height = calibration.compute_height(feet_im, parts_im)

            bbox = person['bbox']
            bbox['p1'][1] -= 25
            text_color = (255, 255, 255)
            font = ImageFont.truetype("/Library/Fonts/arial.ttf", 25)
            draw.text(bbox['p1'], '{0:.2f}'.format(height), text_color, font)

        if not os.path.exists(OUTPUT_DIR):
            os.makedirs(OUTPUT_DIR)
        pil_image.save(os.path.join(OUTPUT_DIR, image_name))

    end_time = time.time()
    print("Running time: %.2f s (%.2f minutes)" %
          (round(end_time - start_time, 2), round((end_time - start_time) / 60, 2)))

