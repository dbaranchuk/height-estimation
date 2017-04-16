import numpy as np

web_calibration = {
    'intrinsic' : [[ 1.4141663290265992e+03, 0., 9.4543373843712743e+02],
                  [ 0., 1.4186903309255270e+03, 4.9882814117191981e+02],
                  [ 0., 0., 1.]],
    'rotation' : [[ 8.4205874916034151e-01, -5.3904494260258895e-01, 1.9173231785290568e-02],
                 [ 2.4392531329014799e-01, 4.1226573035971098e-01, 8.7780260258629728e-01],
                 [ -4.8107951993285275e-01, -7.3448452497349836e-01, 4.7863867173017927e-01]],
    'translation' : [ 1.0079020121414087e+03, 7.4787263478618445e+02, 3.7828233199012411e+03],
    'distortion' : [ 5.9773280266265704e-03, -6.8812400681584937e-02, -3.2746095516290839e-03,
                    -3.4638808119244616e-03, 3.4946717700923085e-03 ]
}

class Calibration(object):
    __slots__ = ['__intrinsic', '__translation',
                 '__rotation', '__distortion', '__projection']

    def __init__(self, calibration_data):
        self.__intrinsic = np.matrix(calibration_data['intrinsic'])
        self.__translation = np.matrix(calibration_data['translation'])
        # Rotation Matrix
        self.__rotation = np.matrix(calibration_data['rotation'])
        # Projection Matrix
        self.__projection = self.__intrinsic * np.hstack((
            self.__rotation,
            self.__translation.T
        ))
        # Distortion K1, K2, P1, P2, K3
        self.__distortion = calibration_data['distortion']


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
        A = np.vstack((self.__projection,
                       [0, 0, 1, 0]))
        b = np.append(point, [1, 0])
        buf = np.linalg.solve(A, b)
        return buf[:3] / buf[3]

    def get_world_part(self, foot_w, part_im):
        A = np.vstack((self.__projection,
                       [-foot_w[1], foot_w[0], 0, 0]))
        b = np.append(part_im, [1, 0])
        res = np.linalg.solve(A, b)
        return res[:3] / res[3]

    def compute_height(self, foot_im, part_im):
        foot_w = self.get_groundplane_proj(foot_im)
        part_w = self.get_world_part(foot_w, part_im)
        height = np.linalg.norm(part_w-foot_w)
        return height


if __name__ == '__main__':
    calibration = Calibration(web_calibration)

    top_pt = (776, 217)
    bottom_pt = (776, 1046)

    top_pt, bottom_pt = calibration.undistort_points([top_pt, bottom_pt])
    height = calibration.compute_height(bottom_pt, top_pt)
    print height
