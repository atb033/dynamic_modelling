from sympy import symbols, Matrix, eye, cos, sin, simplify, pi


class SerialRobotKinematics:
    def __init__(self, dh_table):
        self.dh_table = dh_table
        self.transformation_matrices = []
        self.T0end = eye(4)

        self.compute_transformation_matrices()

    def compute_transformation_matrices(self):
        """
        This function computes a whole bunch of matrices from the DH-table
        """
        for dh_parameters in dh_table:
            T = self.get_trans_matrix(dh_parameters)
            self.transformation_matrices.append(T)
            self.T0end *= T
            self.T0end = simplify(self.T0end)

    def get_trans_matrix(self, dh_parameters):
        """
        This function returns a simplified symbolic transformation matrix based on 
        DH-parametrs (alpha, a, d, theta)
        """

        alpha_in, a_in, d_in, theta_in = dh_parameters

        alpha, a, d, theta = symbols('alpha a d theta')
        trans_matrix = eye(4)
        trans_matrix[0] = cos(theta)
        trans_matrix[1] = -sin(theta)
        trans_matrix[2] = 0
        trans_matrix[3] = a
        trans_matrix[4] = sin(theta) * cos(alpha)
        trans_matrix[5] = cos(theta) * cos(alpha)
        trans_matrix[6] = -sin(alpha)
        trans_matrix[7] = -sin(alpha) * d
        trans_matrix[8] = sin(theta) * sin(alpha)
        trans_matrix[9] = cos(theta) * sin(alpha)
        trans_matrix[10] = cos(alpha)
        trans_matrix[11] = cos(alpha) * d
        trans_matrix[12] = 0
        trans_matrix[13] = 0
        trans_matrix[14] = 0
        trans_matrix[15] = 1

        trans_matrix = trans_matrix.subs([
            (alpha, alpha_in),
            (a, a_in),
            (d, d_in),
            (theta, theta_in),
        ])

        trans_matrix = simplify(trans_matrix)

        return trans_matrix


if __name__ == "__main__":

    print('3RRR kinematics')
    l1, l2, theta1, theta2, theta3 = symbols('l1 l2 theta1 theta2 theta3')
    dh_table = [[0, 0, 0, theta1],
                [0, l1, 0, theta2],
                [0, l2, 0, theta3]]
    srk = SerialRobotKinematics(dh_table)
    print(srk.T0end)

    print('PUMA kinematics')
    theta1, theta2, theta3, a2, d3 = symbols('theta1 theta2 theta3 a2 d3')

    dh_table = [[0, 0, 0, theta1],
                [-pi/2, 0, 0, theta2],
                [0, a2, d3, theta3], ]
    srk = SerialRobotKinematics(dh_table)

    print(srk.T0end)
