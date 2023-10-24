from fiction.pyfiction import *
import unittest
import os

dir_path = os.path.dirname(os.path.realpath(__file__))


class TestCriticalTemperature(unittest.TestCase):

    def test_perturber_and_DB_pair(self):
        layout = sidb_layout((10, 10))
        layout.assign_cell_type((0, 1), sidb_technology.cell_type.NORMAL)
        layout.assign_cell_type((4, 1), sidb_technology.cell_type.NORMAL)
        layout.assign_cell_type((6, 1), sidb_technology.cell_type.NORMAL)

        # params = critical_temperature_params()
        # params.engine = simulation_engine.EXACT
        # params.temperature_mode = critical_temperature_mode.NON_GATE_BASED_SIMULATION
        #
        # stats = critical_temperature_stats()
        #
        # cds = charge_distribution_surface(layout)
        #
        # critical_temperature(cds, params, stats)
        #
        # self.assertEqual(stats.algorithm_name, "ExGS")
        # self.assertEqual(stats.critical_temperature, 400)
        # self.assertEqual(stats.num_valid_lyt, 1)

    def test_bestagon_inv(self):
        layout = read_sqd_layout(dir_path + "/../../../resources/hex_inv_diag_0.sqd", "inverter_input_0")

        # params = critical_temperature_params()
        # params.engine = simulation_engine.APPROXIMATE
        # params.temperature_mode = critical_temperature_mode.GATE_BASED_SIMULATION
        #
        # stats = critical_temperature_stats()
        #
        # cds = charge_distribution_surface(layout)
        #
        # critical_temperature(cds, params, stats)
        #
        # self.assertEqual(stats.algorithm_name, "QuickSim")
        # self.assertLessEqual(stats.critical_temperature, 400)
        # self.assertGreater(stats.num_valid_lyt, 1)


if __name__ == '__main__':
    unittest.main()
