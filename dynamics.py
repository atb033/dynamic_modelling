from kinematics import SerialRobotKinematics
from sympy import symbols, diag, diff, eye, simplify, Matrix, pi


class SerialRobotDynamics:
    def __init__(self, kinematics):

        self.kinematics = kinematics

        # define coordinates of COM (here we assume that the COM lies on the link, and xCOM1 is zero)
        self.xCOM2, self.xCOM3 = symbols('xCOM2 xCOM3')

        # define masses
        self.m1, self.I1, self.m2, self.I2, self.m3, self.I3 = symbols(
            'm1 I1 m2 I2 m3 I3')

        self.compute_total_KE()
        self.compute_mass_matrix()

        self.compute_PE()
        self.compute_Q()

    def compute_total_KE(self):

        # KE1
        I1 = diag(self.m1, self.I1, self.I1)
        KE1 = self.compute_KE(
            0, I1, self.kinematics.vs[0], self.kinematics.omegas[0])

        # KE2
        I2 = diag(0, self.I2, self.I2)
        KE2 = self.compute_KE(
            self.m2, I2, self.kinematics.vs[1], self.kinematics.omegas[1])

        # KE3
        I3 = diag(0, self.I3, self.I3)
        KE3 = self.compute_KE(
            self.m3, I3, self.kinematics.vs[2], self.kinematics.omegas[2])

        self.KE = simplify(KE1 + KE2 + KE3)

    def compute_mass_matrix(self):

        self.mass_mat = eye(3)
        self.mass_mat[0, 0] = simplify(diff(
            self.KE, self.kinematics.theta1dot, self.kinematics.theta1dot))
        self.mass_mat[1, 1] = simplify(diff(
            self.KE, self.kinematics.theta2dot, self.kinematics.theta2dot))
        self.mass_mat[2, 2] = simplify(diff(
            self.KE, self.kinematics.theta3dot, self.kinematics.theta3dot))
        self.mass_mat[0, 1] = simplify(diff(
            self.KE, self.kinematics.theta1dot, self.kinematics.theta2dot))
        self.mass_mat[0, 2] = simplify(diff(
            self.KE, self.kinematics.theta1dot, self.kinematics.theta3dot))
        self.mass_mat[1, 2] = simplify(diff(
            self.KE, self.kinematics.theta2dot, self.kinematics.theta3dot))
        self.mass_mat[1, 0] = self.mass_mat[0, 1]
        self.mass_mat[1, 2] = self.mass_mat[2, 1]

    def compute_PE(self):
        g = Matrix([0, 0, 9.8])

        T01 = self.kinematics.transformation_matrices[0]
        p1 = T01[3, :3].transpose() + T01[:3, :3] * Matrix([0, 0, 0])
        PE1 = -self.m1 * g.dot(p1)

        T02 = T01 * self.kinematics.transformation_matrices[1]
        p2 = T02[3, :3].transpose() + T02[:3, :3] * Matrix([self.xCOM2, 0, 0])
        PE2 = -self.m2 * g.dot(p2)

        T03 = T02 * self.kinematics.transformation_matrices[2]
        p3 = T03[3, :3].transpose() + T03[:3, :3] * Matrix([self.xCOM3, 0, 0])
        PE3 = -self.m2 * g.dot(p3)

        self.PE = PE1 + PE2 + PE3

    def compute_Q(self):
        theta1, theta2, theta3 = symbols('theta1 theta2 theta3')
        joint_variables = Matrix([theta1, theta2, theta3])

        print(simplify(diff(self.PE, joint_variables)).transpose())

        self.Qmat = simplify(diff(self.PE, joint_variables)).transpose()

    def compute_KE(self, m, I, v, omega):
        return 0.5 * (m * v.transpose() * v + omega.transpose() * I * omega)[0]

    def do_substitutions(self, substitutions):
        for key, value in substitutions.items():
            self.mass_mat = self.mass_mat.subs(key, value)
            self.Qmat = self.Qmat.subs(key, value)
        print(self.mass_mat)
        print(self.Qmat)

    def compute_torque(self):
        # print(self.mass_mat)
        import numpy as np

        thrange = np.linspace(-np.pi,np.pi,10)

        th_combinations = np.array(np.meshgrid(thrange,thrange,thrange)).T.reshape(-1,3)

        torque = []
        max_ang_acc = 15 / 0.55
        max_ang_acc_vec = Matrix([max_ang_acc, max_ang_acc,max_ang_acc])
        torque_eq = (self.mass_mat * max_ang_acc_vec).transpose() + self.Qmat

        tmax = 3 * [-np.inf]

        for th in th_combinations:
            t = torque_eq.subs('theta1', th[0])
            t = t.subs('theta2', th[1])
            t = t.subs('theta3', th[2])
            
            torque.append(t)
            

            for i in range(3):
                if tmax[i] < abs(t[i]):
                    tmax[i] = abs(t[i])

        print(tmax)

if __name__ == "__main__":

    # print('3RRR kinematics')
    # l1, l2, l3, theta1, theta2, theta3 = symbols(
    #     'l1 l2 l3 theta1 theta2 theta3')
    # dh_table = [[0, 0, 0, theta1],
    #             [0, l1, 0, theta2],
    #             [0, l2, 0, theta3],
    #             [0, l3, 0, 0], ]
    # print("PUMA")
    # theta1, theta2, theta3, a2, a3, d3 = symbols(
    #     'theta1 theta2 theta3 a2 a3 d3')

    # dh_table = [[0, 0, 0, theta1],
    #             [-pi/2, 0, 0, theta2],
    #             [0, a2, d3, theta3],
    #             [0, a3, 0, 0]]

    print("Our robot")
    l2, l3, theta1, theta2, theta3 = symbols(
        'l2 l3 theta1 theta2 theta3')
    dh_table = [[0, 0, 0, theta1],
                [-pi/2, 0, 0, theta2],
                [pi/2, l2, 0, theta3],
                [0, l3, 0, 0], ]

    srk = SerialRobotKinematics(dh_table)
    srd = SerialRobotDynamics(srk)

    substitutions = {}
    # length substitutions
    substitutions[l2] = 0.3
    substitutions[l3] = 0.25

    # dynamics substitutions
    m1, I1, m2, I2, m3, I3, xCOM2, xCOM3= symbols('m1 I1 m2 I2 m3 I3 xCOM2 xCOM3')
    substitutions[m1] = 0.5
    substitutions[m2] = 0.5
    substitutions[m3] = 2
    
    substitutions[I1] = 0
    substitutions[I2] = 0
    substitutions[I3] = 2 * (0.25 **2)

    substitutions[xCOM2] = 0.15
    substitutions[xCOM3] = 0.125


    srd.do_substitutions(substitutions)
    srd.compute_torque()
