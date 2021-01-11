from sympy import symbols, Matrix, eye, cos, sin, simplify, pi, zeros, diff


class SerialRobotKinematics:
    def __init__(self, dh_table):
        self.t = symbols('t')
        self.dh_table = dh_table
        self.T0end = eye(4)

        self.theta1dot, self.theta2dot, self.theta3dot = symbols(
            'theta1dot theta2dot theta3dot')

        self.transformation_matrices = []
        self.compute_transformation_matrices()

        self.rotation_matrices = []
        self.compute_rotation_matrices()

        self.origins = []
        self.compute_origins()

        self.omegas = []
        self.angular_vel_propagation()

        self.vs = []
        self.v_tool = None # Velocity of the tool wrt ground

        self.linear_vel_propagation()

        # self.compute_jacobian()

    def compute_transformation_matrices(self):
        """
        This function computes a whole bunch of matrices from the DH-table
        """
        for dh_parameters in self.dh_table:
            T = self.get_trans_matrix(dh_parameters)
            self.transformation_matrices.append(T)
            self.T0end *= T
            self.T0end = simplify(self.T0end)

    def compute_rotation_matrices(self):
        for T in self.transformation_matrices:
            self.rotation_matrices.append(T[:3, :3])

    def compute_origins(self):
        for T in self.transformation_matrices:
            self.origins.append(T[:3, 3])

    def angular_vel_propagation(self):
        z_vec = Matrix([0, 0, 1])
        omega1 = z_vec * self.theta1dot
        omega2 = self.rotation_matrices[1].transpose() * omega1 + z_vec * self.theta2dot
        omega3 = self.rotation_matrices[2].transpose() * omega2 + z_vec * self.theta3dot

        self.omegas.append(omega1)
        self.omegas.append(omega2)
        self.omegas.append(omega3)

        self.omega_tool = simplify( self.T0end[:3,:3] *  omega3)

    def linear_vel_propagation(self):
        # Velocity at the joints
        v1 = zeros(3,1)
        v2 = simplify(self.rotation_matrices[1].transpose() * (v1  + self.omegas[0].cross(self.origins[1])))
        v3 = simplify(self.rotation_matrices[2].transpose() * (v2  + self.omegas[1].cross(self.origins[2])))

        # velocity of the tool
        v4 = simplify((v3  + self.omegas[2].cross(self.origins[3])))
        
        self.vs.append(v1)
        self.vs.append(v2)
        self.vs.append(v3)
        self.vs.append(v4)

        # Velocity of the tool wrt ground
        self.v_tool = simplify(self.T0end[:3,:3] * v4)

    def compute_jacobian(self):
        """This function computes the Jacobian of the tool wrt to the ground"""

        twist = Matrix([self.v_tool, self.omega_tool])

        col1 = diff(twist,self.theta1dot)
        col2 = diff(twist,self.theta2dot)
        col3 = diff(twist,self.theta3dot)

        self.jacobian = col1
        self.jacobian = self.jacobian.col_insert(1,col2)
        self.jacobian = self.jacobian.col_insert(2,col3)

        print(self.jacobian)

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
    l1, l2, l3, theta1, theta2, theta3 = symbols(
        'l1 l2 l3 theta1 theta2 theta3')
    dh_table = [[0, 0, 0, theta1],
                [0, l1, 0, theta2],
                [0, l2, 0, theta3],
                [0, l3, 0, 0], ]
    srk = SerialRobotKinematics(dh_table)
    # print(srk.T0end)

    # print('PUMA kinematics')
    # theta1, theta2, theta3, a2, d3 = symbols('theta1 theta2 theta3 a2 d3')

    # dh_table = [[0, 0, 0, theta1],
    #             [-pi/2, 0, 0, theta2],
    #             [0, a2, d3, theta3], ]
    # srk = SerialRobotKinematics(dh_table)

    # print(srk.T0end)
