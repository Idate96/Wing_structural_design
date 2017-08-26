import numpy as np
import sections
import helpers
import matplotlib.pyplot as plt
# import matplotlib.pyplot as plt
# Moments of inertia should not be loaded but calculated from the geometry of the section
# unfortunately it works like that now: to change for the next work package


class Wing():
    def __init__(self, num_sections, length, root_chord, fuel_case, crit_case):
        # wing geometry
        self.num_sec = num_sections
        self.length = length
        self.root_chord = root_chord
        self.sections = []

        # fuel case
        self.fuel_case = fuel_case

        # critical cases circulations
        self.critical_cases_circulation = [243.44*115.6, 179*140, -81.5*123, 71*141]

        # critical case
        self.critical_case = crit_case

        # init loads
        self.distributed_load = []
        self.forces = []
        self.moments = []
        self.force()
        self.moment()
        self.sections = []
        self.polar_areas = np.zeros((num_sections + 1))
        self.load_sections()
        self.torques = np.zeros((num_sections + 1))

        # structural parameters
        self.internal_wing_loads = np.zeros((num_sections + 1, 3))
        self.moment_areas = np.zeros((num_sections + 1, 2))
        self.E = 72 * 10 ** 9
        self.G = 26.9 * 10 ** 9
        self.bending_deflections = np.zeros((num_sections + 1, 2))
        self.twists = np.zeros((num_sections + 1))

        # flight conditions
        self.speed = [94.37,115,101]
        self.density = 1.25

        # stiffness parameters
        self.exp_twisting = 1
        self.c_bending = 10**-8
        self.c_twisting = 5 * 10**-6
        self.exp_bending = 1

    def reset_coefficients(self):
        self.c_bending = 0.0001
        self.c_twisting = 5 * 10**-6

    def reset_deformations(self):
        self.bending_deflections = np.zeros((num_sec + 1, 2))
        self.twists = np.zeros((self.num_sec + 1))

    def local_chord(self, pos):
        chord = -0.19624 * pos + 4.815
        return chord

    def local_cross_area(self, pos):
        cross_area = 2.177 - (2.177 - 0.348) * (pos/self.length)
        return cross_area

    def local_perimeter(self,pos):
        perimeter = 9.880 - (9.880-3.950) * (pos/self.length)
        return perimeter

    def load_sections(self):
        '''
        Load the subsection of the wing starting from the tip
        and assign to each subsecition the associated moments, forces,
        and distributed loads imputted by the users.
        '''
        self.sections = []
        for i in range(self.num_sec):
            loads = [[], [], []]
            sec_lenght = self.length / self.num_sec
            self.distributed_loads(sec_lenght * i)
            for element in self.distributed_load:
                if element[1] <= sec_lenght * i <= element[2]:
                    loads[0].append(element[0])
            for force in self.forces:
                if float(sec_lenght * i) <= force[1] <= float(sec_lenght * (i + 1)):
                    loads[1].append(force[0])
            for moment in self.moments:
                if sec_lenght * i <= moment[1] <= sec_lenght * (i + 1):
                    loads[2].append(moment[0])
            section = sections.Section(i, sec_lenght * i, sec_lenght * (i + 1), loads)
            section.chord = self.local_chord(sec_lenght * i)
            section.cross_area = self.local_cross_area(sec_lenght*i)
            section.perimeter = self.local_perimeter(sec_lenght*i)
            # section.polar_area_init()
            # self.polar_areas[i] = section.polar_area_nd
            self.sections.append(section)

    def distributed_loads(self, pos):
        '''
        :param pos: position
        list here the dist. loads as follows:
        dist loads = [value, start, end]
        calculate the distributed load value at specified position
        '''

        self.distributed_load = []
        # check correction factor for lift
        lift = np.array(([self.critical_cases_circulation[self.critical_case] *
                          (1-(2*((pos-self.length)/self.length/2))**2)**0.5, 0, self.length]))
        # lift = np.array(([1, 0, self.length]))
        # lift = np.array(([243*(1-(2))]))
        self.distributed_load.append(lift)
        fuel_1 = np.array(([self.fuel_case*(-15386), self.length-0.33, self.length]))
        fuel_2 = np.array(([self.fuel_case*(-10483), self.length - 1.7, self.length - 0.33]))
        fuel_3 = np.array(([self.fuel_case*(-13859), self.length - 6.23, self.length - 1.7]))
        self.distributed_load.append(fuel_1)
        self.distributed_load.append(fuel_2)
        self.distributed_load.append(fuel_3)
        # self.distributed_load = []
        wing_weight = np.array(([-1080*pos, 0,self.length]))
        wing_weight = np.array(([-8161/(self.length - pos + 1)**2, 0, self.length]))
        self.distributed_load.append(wing_weight)



    def moment_area_nd_dist(self, pos):
        '''
        :param pos: position
        :return:
        '''
        value_I_skin = 0.036 * (-0.19624 * pos + 4.815) ** 4 * 0.1303/1000 * 0
        # values for the skin are too high
        moment_area = np.array([value_I_skin + self.c_bending * (self.length - pos) ** self.exp_bending, 0, self.length])
        # moment_area = np.array([value_I_skin + self.c_bending, 0, self.length])
        # print('second moment {0}'.format(moment_area))
        return moment_area

    def polar_area(self, pos):
        # print(self.exp_twisting)
        value_polar = self.c_twisting * (self.length - pos)**self.exp_twisting
        polar_area = np.array([value_polar, 0, self.length])
        # print("polar area {0}" .format(self.polar_areas))
        return polar_area

    def load_moment_nd_area(self):
        # change for the future: each section should have its own moment area not the one with opposite index
        # remove therefore all the reversed statements
        sec_lenght = self.length / self.num_sec
        for i, section in enumerate(self.sections):
                # print("moment area")
                # print(self.moment_area_nd_dist(i * self.length / self.num_sec))
            # i = self.num_sec - i
            if self.moment_area_nd_dist(i * self.length / self.num_sec)[1] <= sec_lenght * i <= \
                    self.moment_area_nd_dist(i * self.length / self.num_sec)[2]:
                section.moment_area_nd = self.moment_area_nd_dist(i * self.length / self.num_sec)[0]
            self.moment_areas[i, :] = section.moment_area_nd
            if self.polar_area(i * self.length / self.num_sec)[1] <= sec_lenght * i <= \
                    self.polar_area(i * self.length / self.num_sec)[2]:
                section.polar_area_nd = self.polar_area(i * self.length / self.num_sec)[0]
            self.polar_areas[i] = section.polar_area_nd
        self.moment_areas[-1, :] = self.moment_areas[-2, :]
        self.polar_areas[-1] = self.polar_areas[-2]

    def force(self):
        '''
        list here forces as follow:
        force = [value, application pt.]
        '''
        force1 = np.array(([10, 0]))
        force2 = np.array(([15, 0.2]))
        force3 = np.array(([20, 0.7]))
        # landing_gear = np.array(([2, 2]))
        # self.forces.append(force1)
        # self.forces.append(force2)
        # self.forces.append(force3)

    def moment(self):
        '''
        list here moments as follows:
        moment = [value, application pt.]
        '''
        # moment1 = np.array(([11, 3]))
        # self.moments.append(moment1)

    def interal_loading(self):
        '''
        Calculate internal loading for each section
        '''
        for i, section in enumerate(self.sections):
            for value in self.internal_wing_loads[i,:]:
                if abs(value) > 10**11:
                    raise ValueError('Internal loads diverged {0}' .format(self.internal_wing_loads))
            section.update_internal_loads(self.internal_wing_loads[i, :])
            if i < self.num_sec+1:
                self.internal_wing_loads[i + 1] = section.internal_loads
                # aero_moment = 0.5 * self.density * section.cm * self.speed[self.critical_case] ** 2 * section.chord ** 2\
                #               * (section.pos_e - section.pos_i)
                # aero_moment = 0
                # print('loads {0}' .format(section.loads[0][0]))
                # print('shear and aero {0}, {1}' .format(section.shear_center, section.aero_center))
                self.torques[i] = section.loads[0][0] * (section.pos_e * section.pos_i) * (section.shear_center - section.aero_center)
                self.torques[i] += (section.loads[0][1]) * (section.pos_e - section.pos_i) * \
                                    (section.cg - section.shear_center)
                if i*(section.pos_e-section.pos_i) > (self.length - 6.2):
                    self.torques[i] += (section.loads[0][2]) * (section.pos_e - section.pos_i) * \
                                       (section.cg - section.shear_center)
                if abs(self.torques[i+1]) > 10**11:
                    raise ValueError("torque diverged {0}\n" +
                                     "section chord {1} \n" +
                                     "section pos e {2} \n" +
                                     "section pos i {3}" .format(self.torques, section.chord,
                                                                 section.pos_e, section.pos_i))

            self.internal_wing_loads[-1] = self.internal_wing_loads[-2]
            self.torques[-1] = self.torques[-2]
        # print("wing loads {0}".format(self.internal_wing_loads))
        # print("torques {0}" .format(self.torques))

    def deflections(self):
        '''
        Calculate orientation of the beam (bending_deflection[0]) and deflection
        of the beam in the vertical direction (bending_deflection[1])
        '''
        # bending deflections should be integrated from root to tip
        # while moments are calculated from tip to root
        # print(self.moment_areas)
        rev_moments = list(reversed(self.internal_wing_loads[:, 2]))
        for i, section in enumerate(self.sections):
            self.bending_deflections[i + 1, 0] = self.bending_deflections[i, 0] + (
                rev_moments[i] / (self.E * self.moment_areas[i, 0])) * (section.pos_e - section.pos_i)
            if i == self.num_sec:
                self.bending_deflections[i,0] = self.bending_deflections[i-1,0]
            if abs(self.bending_deflections[i,0]) > 10000:
                raise ValueError("angle bending simulation diverged")
        for i, section in enumerate(self.sections):
            self.bending_deflections[i + 1, 1] = self.bending_deflections[i, 1] + self.bending_deflections[i, 0] * (
                section.pos_e - section.pos_i)
            if i == self.num_sec:
                self.bending_deflections[i,1] = self.bending_deflections[i-1,1]
            if abs(self.bending_deflections[i,1]) > 100 and i == 0:
                self.bending_deflections[i,1] = 0
            if abs(self.bending_deflections[i,1]) > 10000:
                raise ValueError("deflection simulaltion diverged {0}".format(self.bending_deflections))
        # print('bending deflecitons {0}' .format(self.bending_deflections))
        # to make the consistent the x axis let's flip the deflections
        self.bending_deflections[:,1] = list(reversed(self.bending_deflections[:,1]))
        self.bending_deflections[:,0] = list(reversed(self.bending_deflections[:,0]))

    def twist(self):
        # twist calculation
        for i, section in enumerate(self.sections):
            # if i == 0:
            #     pass
            self.twists[i + 1] = self.twists[i] + self.torques[self.num_sec - i] * \
                                                      (section.pos_e - section.pos_i) / (section.polar_area_nd * self.G)
            # print('twist and torque {0}, {1}' .format(self.twists[i], self.torques[self.num_sec - i]))
            if abs(self.twists[i]) > 1000:
                raise ValueError("twist simulation diverged {0}" .format(self.twists))
        self.twists = list(reversed(self.twists))
        # print('twists {0}' .format(self.twists))


def loads_plotter(value, crit_case):
    wing = Wing(200, 14.35, 4.815, 0, crit_case)
    fuel_cases = [0.1, 0.5, 0.9]
    distributions = 3*[[0]]
    index = None
    units = None
    if value == 'distributed loads':
        index = 0
        units = '[N/m]'
    elif value == 'loads':
        index = 1
        units = '[N]'
    elif value == 'moments':
        index = 2
        units = '[Nm]'
    elif value == 'torques':
        index = 3
        units = '[Nm]'
    else:
        raise ValueError('Invalid input value')
    lines = []
    x = [wing.length * i / wing.num_sec for i in range(wing.num_sec + 1)]
    for i, fuel_case in enumerate(fuel_cases):
        wing_case = Wing(200, 14.35, 4.815, fuel_case, crit_case)
        helpers.run_sim(wing_case)
        if index < 3:
            distributions[i] = wing_case.internal_wing_loads[:,index]
            array = np.vstack((distributions[i], x))
            np.savetxt('loads/internal_load' + value + str(fuel_case) + 'crit' + str(crit_case) + '.txt', array, delimiter=',')
        if index == 3:
            distributions[i] = wing_case.torques
            array = np.vstack((distributions[i], x))
            np.savetxt('loads/internal_load' + value + str(fuel_case) + 'crit' + str(crit_case) + '.txt', array,
                       delimiter=',')
        legend = str(fuel_case*100) + "% of fuel"
        line1, = plt.plot(x,distributions[i], label=legend, linewidth=2)
        lines.append(line1)
        # plt.plot(x, distributions[i], label=legend)
    # first_legend = plt.legend(handles=lines, loc=1)
    # ax = plt.gca().add_artist(first_legend)
    plt.legend(handles=lines, loc=4)

    plt.xlabel("distance from tip [m]")
    plt.ylabel(value + units)
    plt.title(value + " spanwise")
    path = 'plots/'
    # plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
    #            ncol=2, mode="expand", borderaxespad=0.)
    plt.savefig(path + value + ' ' + 'case ' + str(crit_case) +'.png', bbox_inches='tight')
    # close figure
    # plt.show()
    plt.close()
    # plt.show()


def test(wing):
    # num_sec = 2000
    # wing = Wing(num_sec, 14.35, 4.815)
    wing.interal_loading()
    # print("Section num, dist. load, load, moments")
    # for i in range(num_sec):
    #     wing.sections[i].printer()

    # helpers.wing_grapher(wing, 'loads')
    wing.load_moment_nd_area()
    wing.deflections()
    wing.twist()
    print('wing polar area {0}' .format(wing.polar_areas))
    # print("torques {0}" .format(wing.torques))
    helpers.wing_grapher(wing, 'distributed loads')
    helpers.wing_grapher(wing, 'loads')
    helpers.wing_grapher(wing, 'moments')
    # helpers.wing_grapher(wing, 'deflections')
    helpers.wing_grapher(wing, 'torques')
    helpers.wing_grapher(wing, 'twists')
    helpers.wing_grapher(wing, 'inertias')


if __name__ == '__main__':
    num_sec = 200
    wing = Wing(num_sec, 14.35, 4.815, 0.2, 0)
    # loads_plotter('torques', 0)
    loads = ['distributed loads', 'loads', 'moments', 'torques']
    # for load in loads:
    #     loads_plotter(load, 3)
    for i in range(4):
        for load in loads:
            loads_plotter(load, i)

    # distributions = helpers.optimal_polar_requirements(wing)
    # distributions_bending = helpers.optimal_bending_requirements(wing)

    # print("tip section loads: {0}, moment_area : {1}"
    #       .format(wing.sections[1].internal_loads, wing.sections[1].moment_area_nd))
    # print("root section loads: {0}, moment_area : {1}"
    #       .format(wing.sections[-1].internal_loads, wing.sections[-1].moment_area_nd))

    # x = [wing.length*i/wing.num_sec for i in range(wing.num_sec+1)]

    # for distribution in distributions:
    #     plt.plot(x, np.array(distribution)[0])
    #     plt.show()
    # thickness, polar_distribution = helpers.optimal_thickness(num_sec, 0.0000005)
    # print('Thickness {0} \n'
    #       'polar_distribution {1}' .format(thickness, polar_distribution))