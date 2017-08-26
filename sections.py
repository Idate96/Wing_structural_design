import numpy as np


class Section:
    def __init__(self, id, pos_i, pos_e, loads):

        self.id = id
        self.pos_i = pos_i
        self.pos_e = pos_e
        self.loads = loads
        self.thickness = 0.0001
        self.cm = - 0.
        self.geo_inertia = None
        self.cross_area = None
        self.perimeter = None
        self._polar_area_nd = None
        self._moment_area_nd = None
        self._distributed_load = None
        self._forces = None
        self._moments = None
        self.parse_loads()
        self.internal_loads = None
        self._chord = None
        self.shear_center = None
        self.aero_center = None
        self.lift = None
        self.cg = None

    def shear_areo_centers(self):
        self.shear_center = 0.4 * self.chord
        self.aero_center = 0.25 * self.chord
        self.cg = 0.4 * self.chord

    @property
    def chord(self):
        return self._chord

    @chord.setter
    def chord(self, chord):
        self._chord = chord
        self.shear_areo_centers()

    @property
    def moment_area_nd(self):
        return self._moment_area_nd

    @moment_area_nd.setter
    def moment_area_nd(self, moment_area_nd):
        self._moment_area_nd = moment_area_nd

    @property
    def polar_area_nd(self):
        return self._polar_area_nd

    @polar_area_nd.setter
    def polar_area_nd(self, polar_area):
        self._polar_area_nd = polar_area

    @property
    def distributed_load(self):
        return self._distributed_load

    @distributed_load.setter
    def distributed_load(self, distributed_load):
        self.distributed_load = distributed_load

    @property
    def forces(self):
        return self._forces

    @forces.setter
    def forces(self, forces):
        self._forces = forces

    @property
    def moments(self):
        return self._moments

    @moments.setter
    def moments(self, moments):
        self._moments = moments

    def polar_area_init(self):
        print("thickness {0}" .format(self.thickness))
        self.polar_area_nd = 4 * self.cross_area**2 * self.thickness/self.perimeter

    def parse_loads(self):
        # print("loads:{0}".format(self.loads))
        self._distributed_load = np.sum(self.loads[0])

        self._forces = np.sum(self.loads[1])
        # print(self._forces)
        self._moments = np.sum(self.loads[2])
        if not self._moments:
            self._moments = 0
        if not self._forces:
            self._forces = 0
        if not self._distributed_load:
            self._distributed_load = 0
            # print(self._distributed_load, self._forces, self._moments)

    def update_internal_loads(self, loads):
        """
        :type loads: np.array
        """
        loads[0] = self._distributed_load
        loads[1] += self._forces + self.distributed_load * (self.pos_e - self.pos_i)
        loads[2] += self.moments + loads[1] * (self.pos_e - self.pos_i)
        self.internal_loads = loads
        # print('loads section' .format(print(loads)))
        return loads

    def printer(self):
        pass
        # print(self.internal_loads)
        # print("Section {0}: {1:.2f} N/m, {2:.2f} N, {3:.2f} Nm \n"
        #       .format(self.id,self.internal_loads[0], list(self.internal_loads[1]), list(self.internal_loads[2])))

#
# def test():
#     section = Section(0,3, np.array(([0.],[1.],[2.])))
#     loads = np.array(([1.],[1.],[1.]))
#     print('loads')
#     print(section.distributed_load, section.forces, section.moments)
#     print("updated loads")
#     print(section.update_internal_loads(loads))
#
# test()
#
