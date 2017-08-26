import matplotlib.pyplot as plt
import numpy as np


def wing_grapher(wing, value):
    x = [wing.length * i / wing.num_sec for i in range(wing.num_sec + 1)]
    if value == 'distributed loads':
        y = wing.internal_wing_loads[:, 0]
        units = '[N/m]'
    elif value == 'loads':
        y = wing.internal_wing_loads[:, 1]
        units = '[N]'
    elif value == 'moments':
        y = wing.internal_wing_loads[:, 2]
        units = '[Nm]'
    elif value == 'inertias':
        y = list(reversed(wing.moment_areas[:, 1]))
        units = '[m^4]'
    elif value == 'deflections':
        y = wing.bending_deflections[:,1]
        units = '[m]'
    elif value == 'torques':
        y = wing.torques
        units = '[Nm]'
    elif value == 'twists':
        y = wing.twists
        units = '[rad]]'
    elif value == 'polar areas':
        y = list(reversed(wing.polar_areas))
        units = '[m^4]'
    else:
        raise ValueError("Invalid value")

    plt.plot(x, y)
    plt.xlabel("distance from tip [m]")
    plt.ylabel(value + units)
    plt.title(value + " spanwise")
    plt.show()


def run_sim(wing):
    print('NEW SIM ----------------------')
    # print('new sim thick {0}' .format(wing.sections[1].thickness))
    wing.interal_loading()
    wing.load_moment_nd_area()
    # wing.deflections()
    # wing.twist()
    # print(wing.torques)
    # print("twist {0}" .format(wing.twists))
    # wing_grapher(wing, 'loads')
    # wing_grapher(wing, 'moments')
    # wing_grapher(wing, 'deflections')
    # wing_grapher(wing, 'torques')
    # wing_grapher(wing, 'twists')


#
# def optimal_thickness(num_sec, init_thick):
#     wing = main.Wing(num_sec, 10, 4.815)
#     twist = 1
#     t = init_thick + np.zeros(wing.num_sec + 1)
#     while twist > 0.175:
#         t += 0.00001
#         wing = main.Wing(num_sec, 10, 4.815)
#         for i, section in enumerate(wing.sections):
#             section.thickness = t[i]
#             print("enumerate thick {0}" .format(section.thickness))
#         run_sim(wing)
#         twist = wing.twists[-2]
#     t_optimal = t[0]
#     return t_optimal, wing.polar_areas

#
# def optimal_polar_requirements(wing):
#     x = [wing.length * i / wing.num_sec for i in range(wing.num_sec + 1)]
#     x = np.array(x)
#     exp = [0,1,1.3,1.6,2,2.5]
#     good_distributions = len(exp)*[[0]]
#     req_met = False
#     for i, exponent in enumerate(exp):
#         twist = 1
#         wing.exp_twisting = exponent
#         for j in range(len(wing.critical_cases_circulation)):
#             wing.critical_case = j
#             print('exponent {0}' .format(i))
#             while twist > 0.175:
#                 if exponent < 2:
#                     wing.c_twisting += 10**-8
#                 elif 2 <= exponent <= 2.3:
#                     wing.c_twisting += 10**-7
#                 elif 2.3 < exponent <= 2.5:
#                     wing.c_twisting += 10**-6
#                 run_sim(wing)
#                 twist = wing.twists[1]
#                 print('twist: {0}' .format(twist))
#
#         good_distributions[i].append(wing.polar_areas)

# return good_distributions


def optimal_polar_requirements(wing):
    x = [wing.length * i / wing.num_sec for i in range(wing.num_sec + 1)]
    x = np.array(x)
    # exp = [0.78, 0, 1, 1.3, 1.6, 2, 2.5]
    exp = [i / 10 for i in range(30)]
    exp = [1.5]

    good_distributions = len(exp) * [[0]]
    for i, exponent in enumerate(exp):
        print('EXP : {0} -----------'.format(exponent))
        twist = 1
        wing.exp_twisting = exponent
        wing.reset_coefficients()
        for j in range(len(wing.critical_cases_circulation)):
            wing.critical_case = j
            wing.load_sections()
            if j > 0:
                wing.reset_deformations()
                run_sim(wing)
                twist = wing.twists[1]
                print('twist new case {0:.4f}'.format(twist))

            while abs(twist) > 0.175:
                wing.reset_deformations()
                wing.c_twisting += 10 ** -6
                run_sim(wing)
                twist = wing.twists[1]
                print('twist : {0:.3f} \n'
                      'case : {1}'.format(twist, j))
            print('c_twisting: {0} \n'
                  'c_exponent: {1}'.format(wing.c_twisting, wing.exp_twisting))
            wing_grapher(wing, 'twists')
            wing_grapher(wing, 'polar areas')
        good_distributions[i] = wing.polar_areas
        # wing_grapher(wing,'twists')
        # wing_grapher(wing, 'polar areas')
        save_distribution(exponent, 'polar', wing.polar_areas, x)
    return good_distributions


def optimal_bending_requirements(wing):
    x = [wing.length * i / wing.num_sec for i in range(wing.num_sec + 1)]
    x = np.array(x)
    # exp = [1.5, 2, 2.5, 3]
    exp = [2]
    # exp = [i/10 for i in range(30)]
    good_distributions = len(exp) * [[0]]
    for i, exponent in enumerate(exp):
        print('EXP : {0} -----------'.format(exponent))
        deflection = 10
        wing.exp_bending = exponent
        wing.reset_coefficients()
        if 1 < exponent < 2:
            wing.c_bending = 10 ** -6
        elif 0 < exponent < 1:
            wing.c_bending = 10**-5
        elif exponent > 2:
            wing.c_bending = 5 * 10 ** -7

        for j in range(len(wing.critical_cases_circulation)):
            wing.critical_case = j
            wing.load_sections()
            if j > 0:
                wing.reset_deformations()
                run_sim(wing)
                deflection = wing.bending_deflections[1, 1]
                print('deflection new case {0:.4f}'.format(deflection))

            while abs(deflection) > wing.length * 0.15 * 2:
                wing.reset_deformations()
                if exponent == 0:
                    wing.c_bending += 10 ** -6
                else:
                    wing.c_bending += 10 ** -7 / (5 * exponent)
                run_sim(wing)
                deflection = wing.bending_deflections[1, 1]
                print('deflection : {0:.3f} \n'
                      'case : {1}'.format(deflection, j))
            print('c_bending: {0} \n'
                  'c_exponent: {1}'.format(wing.c_bending, wing.exp_bending))
            # wing_grapher(wing, 'deflections')
            wing_grapher(wing, 'inertias')
        good_distributions[i] = wing.moment_areas[:, 0]
        # wing_grapher(wing,'deflectons')
        # wing_grapher(wing, 'inertiasi')
        save_distribution(exponent, 'bending', wing.moment_areas[:, 0], x)
    return good_distributions


def save_distribution(exp, type, polar_distribution, x):
    array = np.vstack((polar_distribution, x))
    path = type + '_distributions/'
    file_name = type + '_distribution_exp' + str(exp) + '.txt'
    np.savetxt(path + file_name, array, delimiter=',')
