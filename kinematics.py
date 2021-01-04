from sympy import symbols, Matrix, eye, cos, sin, simplify, pi


def get_trans_matrix(alpha_in, a_in, d_in, theta_in):
    """
    This function returns a simplified symbolic transformation matrix based on DH
    parametrs (alpha, a, d, theta)
    """

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
    T01 = get_trans_matrix(0, 0, 0, theta1)
    T12 = get_trans_matrix(0, l1, 0, theta2)
    T23 = get_trans_matrix(0, l2, 0, theta3)

    T03 = T01 * T12 * T23
    T03 = simplify(T03)
    
    print(T03)
    print('3RRR kinematics')

    print('PUMA kinematics')
    theta1, theta2, theta3, a2, d3 = symbols('theta1 theta2 theta3 a2 d3')
    T01 = get_trans_matrix(0, 0, 0, theta1)
    T12 = get_trans_matrix(-pi/2, 0, 0, theta2)
    T23 = get_trans_matrix(0, a2, d3, theta3)

    T03 = simplify(T01 * T12 * T23)
    print(T03)