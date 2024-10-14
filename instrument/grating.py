import math
import numpy as np


class Grating:
    def __init__(
        self,
        line_density: float,            # [l/um] #rho
        blaze_angle: float,             # [deg]
        eps_angle: float,               # [deg]
        # Efficiency
        ilg: float,
        n_p: float,
    ) -> None:
        self.line_density = line_density
        self.blaze_angle = blaze_angle
        self.eps_angle = eps_angle
        self.n_p = n_p * 3

    def set_SIT(self, SIT):
        self.SIT = SIT

    def set_wavelength_orders(self):
        try:
            self.order_table = np.loadtxt(
                self.orders_table_file, delimiter=" ", skiprows=1)

            self.sx_wavelegnth_per_order = []
            self.dx_wavelegnth_per_order = []

            # Create from order table two list containing the wavelength min and max for each order
            for n_ord in self.n_orders:
                self.sx_wavelegnth_per_order.append(
                    self.order_table[self.order_table[:, 0] == n_ord, :][0][2+self.SIT])
                self.dx_wavelegnth_per_order.append(
                    self.order_table[self.order_table[:, 0] == n_ord, :][-1][2+self.SIT])
        except Exception as e:
            print(e)
            exit()


class Echelle(Grating):
    def __init__(
        self,
        type: str,

        line_density: float,            # [l/um] #rho
        blaze_angle: float,             # [deg]
        eps_angle: float,               # [deg]
        # Efficiency
        ilg: float,
        n_p: float,

        n_order_start,
        n_order_end,
        orders_table_file: str,
        orders_table_spectral_shift: float = None,
        orders_table_spatial_shift: float = None
    ) -> None:
        super().__init__(line_density, blaze_angle, eps_angle, ilg, n_p)
        self.type_name = type

        self.n_order_start = n_order_start
        self.n_order_end = n_order_end

        # Orders Information
        if n_order_start > n_order_end:
            self.n_orders = [*range(n_order_start, n_order_end-1, -1)]
        elif n_order_start < n_order_end:
            self.n_orders = [*range(n_order_start, n_order_end+1, 1)]
        else:
            self.n_orders = [n_order_start]
        self.len_n_orders = len(self.n_orders)

        # Load Orders Table file
        self.orders_table_file = orders_table_file

        self.orders_table_spectral_shift = orders_table_spectral_shift if orders_table_spectral_shift else 0
        self.orders_table_spatial_shift = orders_table_spatial_shift if orders_table_spatial_shift else 0

    # Different for each grating

    def get_efficiency(self):
        wave_matrix = np.zeros((self.len_n_orders, self.n_p))
        b_phase = np.zeros((self.len_n_orders, self.n_p))
        b = np.zeros((self.len_n_orders, self.n_p))

        # d=grating const and N=number of lines
        d = 1/self.line_density  # [um/l]

        for i in range(0, self.len_n_orders):
            wave_matrix[i] = np.linspace(
                self.sx_wavelegnth_per_order[i], self.dx_wavelegnth_per_order[i], self.n_p)
            for index, lmbd in enumerate(wave_matrix[i]):
                beta = math.asin(((self.n_orders[i] * self.line_density * lmbd) / math.cos(
                    math.radians(self.eps_angle))) - math.sin(math.radians(self.blaze_angle)))

                b_phase[i, index] = (math.pi*d*math.cos(math.radians(self.blaze_angle))/lmbd) * (math.sin(
                    math.radians(self.blaze_angle-self.blaze_angle)) + math.sin(beta-math.radians(self.blaze_angle)))

                b[i, index] = np.sinc(b_phase[i][index]/math.pi)**2

        b = cut_spurius_efficiency(b, self.len_n_orders)

        # B is now deprecated, was used to adjust the grating efficiency
        # NOW grating must be included in instrument efficiency file already adjusted!
        return wave_matrix, self.len_n_orders, self.n_p


class Binary(Grating):
    def __init__(
        self,
        type: str,

        line_density: float,            # [l/um] #rho
        blaze_angle: float,             # [deg]
        eps_angle: float,               # [deg]
        # Efficiency
        ilg: float,
        n_p: float,

        orders_table_file: str,
        orders_table_spectral_shift: float = None,
        orders_table_spatial_shift: float = None
    ) -> None:
        super().__init__(line_density, blaze_angle, eps_angle, ilg, n_p)
        self.type_name = type

        # Orders Information
        self.n_orders = [1]
        self.len_n_orders = len(self.n_orders)

        # Load Orders Table file
        self.orders_table_file = orders_table_file

        self.orders_table_spectral_shift = orders_table_spectral_shift if orders_table_spectral_shift else 0
        self.orders_table_spatial_shift = orders_table_spatial_shift if orders_table_spatial_shift else 0

    # Different for each grating

    # TDB: Change name in get wavematrix, after CUBES RUN smooth
    def get_efficiency(self):
        wave_matrix = np.zeros((self.len_n_orders, self.n_p))

        for i in range(0, self.len_n_orders):
            wave_matrix[i] = np.linspace(
                self.sx_wavelegnth_per_order[i], self.dx_wavelegnth_per_order[i], self.n_p)

        # B is now deprecated, was used to adjust the grating efficiency
        # NOW grating must be included in instrument efficiency file already adjusted!
        return wave_matrix, self.len_n_orders, self.n_p


def cut_spurius_efficiency(b, len_n_orders):
    for i in range(0, len_n_orders):
        max_idx = np.where(b[i] == np.amax(b[i]))[0]
        for k in range(0, max_idx[0]):
            if b[i][k] > b[i][k+1]:
                b[i][k] = 0
            else:
                break
        for k in reversed(range(max_idx[0], len(b[i]))):
            if b[i][k] > b[i][k-1]:
                b[i][k] = 0
            else:
                break
    return b
