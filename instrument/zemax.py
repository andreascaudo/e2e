import math
import numpy as np
from ..tool import tools
from ..zos.ZOS_connection import ZOS_connection

DEBUG = False


class Zemax:
    def __init__(
        self,
        file_path: str,
        order_table_output = None,
        PSF_map_output = None,
        order_table_flag: bool = True,
        PSF_map_flag: bool = True,
        CFG_path: str = None
    ) -> None:
        self.N_SRE = 50
        self.file_path = file_path
        if order_table_flag != None:
            self.order_table_flag = order_table_flag
        if order_table_output != None and order_table_flag:
            self.order_table_output = order_table_output
        if PSF_map_flag != None:
            self.PSF_map_flag = PSF_map_flag
            if CFG_path == None:
                raise Exception("CFG_path not specified")
            else:
                self.CFG_path = CFG_path
        if PSF_map_output != None and PSF_map_flag:
            self.PSF_map_output = PSF_map_output

        self.order_table = []
        self.psf_map_final = np.empty(0)
        self.zos = ZOS_connection.ZOS_connection()
        self.TheSystem = self.load_zemax()
        
    def load_zemax(self):
        # load local variables
        result = self.zos.TheSystem.LoadFile(self.file_path, True)
        if result:
            print("Zemax loaded successfully")
        else:
            raise Exception("Zemax failed to load")
        return self.zos.TheSystem
    
    def close_connection(self):
        del self.zos
        self.zos = None

    def get_orders_wavelegth_range(self, orders, detector_size, pixel_size):
        print("Getting orders wavelength range from zemax...")
        self.orders_wavelegth_range = {}
        orders_pixel_range = {} # For debugging purposes

        # make orders descending and cast it to float
        orders = [float(order) for order in orders]
        orders.sort(reverse=True)

        wave_start = 0

        delta_target = 0.000001 # Enough precision

        if DEBUG:   
                                                                                                                                # to dict
            self.orders_wavelegth_range = np.load("C:\\Users\\andreascaudo\\Desktop\\orders_wavelegth_range.npy", allow_pickle=True).item()
            return self.orders_wavelegth_range

        for idx, order in enumerate(orders):
            
            delta = 0.1
            lambda_start = 0
            px_start = 0
            lambda_end = 0
            px_end = 0
            lambda_temp = -1
            
            while lambda_start == 0:

                self.TheSystem.MCE.GetOperandAt(3).GetOperandCell(1).DoubleValue = order
                self.TheSystem.MCE.GetOperandAt(4).GetOperandCell(1).DoubleValue = wave_start
                self.TheSystem.MFE.CalculateMeritFunction()

                CenX = self.TheSystem.MFE.GetOperandAt(2).Value
                CenY = self.TheSystem.MFE.GetOperandAt(3).Value

                if CenX == 0:
                    lambda_temp = wave_start
                    wave_start = wave_start + delta
                    continue
                
                R_mat = [[1,0,0], [0,-1,0], [0, 0, -1]]
                T_vect = np.dot((detector_size/2)*pixel_size, [1, 1, 0])

                v = np.array([CenX, CenY, 0])
                v1 = np.dot(R_mat, v)
                v1 = v1 * 1000 #mm to micron
                v2 = v1 + T_vect
                v_tab_x_pix = v2[0] / pixel_size

                if v_tab_x_pix > 0:
                    if delta > delta_target:
                        delta = delta / 10
                        wave_start = lambda_temp
                    else:
                        lambda_start = wave_start
                        px_start = v_tab_x_pix
                else:
                    lambda_temp = wave_start
                    wave_start = wave_start + delta

            
            self.orders_wavelegth_range[order] = [lambda_start]
            orders_pixel_range[order] = [px_start]
            print("Order: " + str(order) + " - lambda_start: " + str(self.orders_wavelegth_range[order][0]) + " - Pixel: " + str(orders_pixel_range[order][0]))

            lambda_temp = -1
            v_tab_x_pix_temp = -1
            delta = 0.01
            while lambda_end == 0:
                self.TheSystem.MCE.GetOperandAt(3).GetOperandCell(1).DoubleValue = order
                self.TheSystem.MCE.GetOperandAt(4).GetOperandCell(1).DoubleValue = wave_start
                self.TheSystem.MFE.CalculateMeritFunction()

                CenX = self.TheSystem.MFE.GetOperandAt(2).Value
                CenY = self.TheSystem.MFE.GetOperandAt(3).Value

                if CenX == 0:
                    lambda_temp = wave_start
                    wave_start = wave_start + delta
                    continue
                
                R_mat = [[1,0,0], [0,-1,0], [0, 0, -1]]
                T_vect = np.dot((detector_size/2)*pixel_size, [1, 1, 0])

                v = np.array([CenX, CenY, 0])
                v1 = np.dot(R_mat, v)
                v1 = v1 * 1000
                v2 = v1 + T_vect
                v_tab_x_pix = v2[0] / pixel_size

                if v_tab_x_pix > detector_size:
                    if delta > delta_target:
                        delta = delta / 10
                        wave_start = lambda_temp
                        v_tab_x_pix = v_tab_x_pix_temp
                    else:
                        lambda_end = lambda_temp
                        px_end = v_tab_x_pix_temp
                        wave_start = lambda_start
                else:
                    lambda_temp = wave_start
                    v_tab_x_pix_temp = v_tab_x_pix
                    wave_start = wave_start + delta

            self.orders_wavelegth_range[order].append(lambda_end)
            orders_pixel_range[order].append(px_end)
            print("Order: " + str(order) + " - lambda_end: " + str(self.orders_wavelegth_range[order][1]) + " - Pixel: " + str(orders_pixel_range[order][1]))
            print("Done!")
        
        return self.orders_wavelegth_range
        

    def get_order_table(self, orders, detector_size, pixel_size, psf_map_pixel_number):
        #N pixels for box convolution: THE EDGE-FRAME TO BE ADDED IN EACH SINGLE BOX
        N_pix_box_x = psf_map_pixel_number / 2
        N_pix_box_y = psf_map_pixel_number / 2

        order_index_MCE = 3
        wavelength_1_index_MCE = 4

        n_points = 5

        for order_index, order_number in enumerate(orders):
            print("Order: " + str(order_number))
            for i in range(0, n_points):
                print("Point: " + str(i+1))
                self.TheSystem.MCE.GetOperandAt(order_index_MCE).GetOperandCell(i+1).DoubleValue = order_number
                self.TheSystem.MFE.GetOperandAt(i*3+1).GetCellAt(2).IntegerValues = i+1

            lam_MCE_vect = np.linspace(self.orders_wavelegth_range[order_number][0], self.orders_wavelegth_range[order_number][1], self.N_SRE)

            for idx_wl in range(self.N_SRE):
                lambda_MCE = lam_MCE_vect[idx_wl]
                for i in range(0, n_points):
                    self.TheSystem.MCE.GetOperandAt(wavelength_1_index_MCE).GetOperandCell(i+1).DoubleValue = lambda_MCE
                    
                self.TheSystem.MFE.CalculateMeritFunction()

                cenX = []
                cenY = []

                for i in range(0, n_points):
                    cenX.append(self.TheSystem.MFE.GetOperandAt(i*3+2).Value)
                    cenY.append(self.TheSystem.MFE.GetOperandAt(i*3+3).Value)

                sX_lam = 1000 * abs(cenX[2] - cenX[3])/pixel_size # Pixel
                sY_lam = 1000 * abs(cenY[2] - cenY[1])/pixel_size # Pixel
                slit_tilt = math.degrees(math.atan((abs(cenX[2] - cenX[1]))/abs(cenY[2] - cenY[1])))
                matX = cenX[1:-1]
                matY = cenY[1:-1]

                R_mat = [[1,0,0], [0,-1,0], [0, 0, -1]]
                T_vect = np.dot((detector_size/2)*pixel_size, [1, 1, 0])

                v = np.array([cenX[0], cenY[0], 0])
                v1 = np.dot(R_mat, v)
                v1 = v1 * 1000 #mm to micron
                v2 = v1 + T_vect
                v_tab_x_pix = v2[0] / pixel_size
                v_tab_y_pix = v2[1] / pixel_size

                #vector 0
                v0 = np.array([matX[0], matY[0], 0])
                v0 = np.dot(R_mat, v0)
                v0 = v0 * 1000 #mm to micron
                v0_2 = v0 + T_vect
                v_tab_x_pix_0 = v0_2[0] / pixel_size
                v_tab_y_pix_0 = v0_2[1] / pixel_size
                v_tab_x_pix_box_0 = np.ceil(v_tab_x_pix_0) + np.floor(N_pix_box_x) # pix
                v_tab_y_pix_box_0 = np.ceil(v_tab_y_pix_0) - np.floor(N_pix_box_y) # pix

                #Vector 2
                v2 = np.array([matX[2], matY[2], 0])
                v2 = np.dot(R_mat, v2)
                v2 = v2 * 1000 #mm to micron
                v2_2 = v2 + T_vect
                v_tab_x_pix_2 = v2_2[0] / pixel_size
                v_tab_y_pix_2 = v2_2[1] / pixel_size
                v_tab_x_pix_box_2 = np.ceil(v_tab_x_pix_2) - np.floor(N_pix_box_x) # pix
                v_tab_y_pix_box_2 = np.ceil(v_tab_y_pix_2) + np.floor(N_pix_box_y) # pix

                # X and Y coord initial pixel of the box, and Nx Ny
                v_tab_x_pix_box = v_tab_x_pix_box_2
                v_tab_y_pix_box = v_tab_y_pix_box_0
                v_tab_Nx_pix_box = np.abs(v_tab_x_pix_box_2 - v_tab_x_pix_box_0)
                v_tab_Ny_pix_box = np.abs(v_tab_y_pix_box_2 - v_tab_y_pix_box_0)

                row = [int(order_number), int(idx_wl+1), lambda_MCE, v_tab_x_pix, v_tab_y_pix, sX_lam, sY_lam, v_tab_x_pix_box, v_tab_y_pix_box, v_tab_Nx_pix_box, v_tab_Ny_pix_box, slit_tilt]
                self.order_table.append(row)
    
        if self.order_table_output != None:
            self.save_order_table()

        #return self.order_table as an array
        return np.array(self.order_table)


    def save_order_table(self):
        if self.order_table_flag and self.order_table != []:
            np.savetxt(self.order_table_output, self.order_table,  
                    fmt="%d %d %f %f %f %f %f %f %f %f %f %f", 
                    header="m_ordd SRE-ID Lam-ID X_coor Y_coor Samp_x Samp_y X_box1 Y_box1 X_BoxN Y_BoxN Slit-T", 
                    comments='', delimiter=' ')
            print("Order table saved successfully")
        else:
            print("Order table not saved")

    def get_PSF_map(self, orders, psf_field_sampling):
        #Init PSF file
        slice = np.dtype([('slice_n', np.uint8, (1,)), #Like this doesn't work with CUBES
                        ('slice', np.ndarray, (psf_field_sampling, psf_field_sampling, self.N_SRE))])
        
        order = np.dtype([('order_n', np.uint8, (1,)), ('order', np.ndarray)])

        self.psf_map_final = np.empty(16, dtype=order)

        order_index_MCE = 3
        wavelength_1_index_MCE = 4

        n_points = 5

        #PSF ACQUISITION
        analysis = self.TheSystem.Analyses.New_HuygensPsf()
        # Settings
        analysis_Settings = analysis.GetSettings()
        analysis_Settings.LoadFrom(self.CFG_path)

        psf_per_order = np.empty((psf_field_sampling, psf_field_sampling, self.N_SRE), dtype=np.ndarray)

        self.TheSystem.MCE.SetCurrentConfiguration(1)

        for order_index, order_number in enumerate(orders):
            print("Order: " + str(order_number))
            self.TheSystem.MCE.GetOperandAt(order_index_MCE).GetOperandCell(self.TheSystem.MCE.CurrentConfiguration).DoubleValue = order_number
                #self.TheSystem.MFE.GetOperandAt(i*3+1).GetCellAt(2).IntegerValues = i+1 CHECK IN ZEMAX IF NEEDED

            print("Current configuration: " + str(self.TheSystem.MCE.CurrentConfiguration))

            lam_MCE_vect = np.linspace(self.orders_wavelegth_range[order_number][0], self.orders_wavelegth_range[order_number][1], self.N_SRE)

            for idx_wl in range(self.N_SRE):
                lambda_MCE = lam_MCE_vect[idx_wl]
                self.TheSystem.MCE.GetOperandAt(wavelength_1_index_MCE).GetOperandCell(self.TheSystem.MCE.CurrentConfiguration).DoubleValue = lambda_MCE
        
                # Run the analysis with the current settings and pull the results
                analysis.ApplyAndWaitForCompletion()
                huygensResults = analysis.GetResults()

                # The results will be split into multple data structures: Header, Intensity data, Metadata
                # We want the intensity data. This is the same data you will see in the Text tab of the analysis
                matrixData = huygensResults.GetDataGrid(0).Values
                xLength = matrixData.GetLength(0)
                yLength = matrixData.GetLength(1)

                # Reformat the data so that it matches the output we see in OS
                # We must also convert the matrix into a Python array for post-processing
                huygensData = self.zos.reshape(matrixData, xLength, yLength, False)

                # From list to numpy array
                np_psf = np.array(huygensData)

                psf_per_order[:,:,idx_wl] = np_psf
        
            slice_final = np.empty(1, dtype=slice)
            print(psf_per_order.shape)

            psf = np.array([(1, psf_per_order)], dtype=slice)
            slice_final[0] = psf
            self.psf_map_final[order_index] = np.array([(order_number, slice_final)], dtype=order)

        if self.PSF_map_flag and self.psf_map_final.size != 0:
            np.save(self.PSF_map_output, self.psf_map_final, allow_pickle=True)
            print("PSF map saved successfully")
        else:
            print("PSF map not saved")

        return self.psf_map_final









