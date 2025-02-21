# E2E Simulator

A comprehensive end-to-end (E2E) simulator for spectrographs that models the expected outcome of spectroscopic observations (or calibration).

## Features

- Full spectroscopic simulation pipeline including:
  - Spectral Energy Distribution (SED) generation
  - Sky background and atmospheric transmission
  - Telescope and instrument efficiency
  - Optical effects (PSF, ..)
  - Detector characteristics
- Support for multiple spectrograph configurations
- Integration with Zemax optical design software
- Parallel processing capabilities for faster simulations
- Configurable through YAML files
- Output in FITS format

## Requirements

- Python 3.x
- Required Python packages (install via pip):
  ```
  numpy
  scipy
  astropy
  pyyaml
  tqdm
  ```
- Zemax OpticStudio (if using Zemax integration)
- skycalc_cli in PATH (for sky background calculations)

## Installation

1. Clone the repository
2. Install the required Python packages:
   ```bash
   pip install -r requirements.txt
   ```
3. Ensure skycalc_cli is properly installed and accessible in your PATH

## Usage

The simulator is configured using YAML files that specify all simulation parameters. Here's how to run a simulation:

```python
from e2e.simulation import load_config
from e2e.main import run

# Load configuration from YAML file
config = load_config('path/to/your/config.yaml')

# Run the simulation
run(config)
```

## Configuration File Structure

The YAML configuration file should contain the following main sections:

```yaml
output:
  output_folder: "path/to/output"  # Where to save the results
  paraxial_model: false           # Enable/disable paraxial model
  drs: false                      # Enable/disable DRS

telescope:
  diameter: 8.2                   # Telescope diameter in meters
  # Additional telescope parameters...

spectrograph:
  name: "INSTRUMENT_NAME"         # Name of the instrument
  arm: 
  
  grating:
    # Choose one of the following grating types:
    
    # For Echelle grating: SOXS example
    type: "Echelle"              # Grating type
    n_order_start: 24            # Starting order number
    n_order_end: 9               # Ending order number
    line_density: 0.072          # Line density in lines/μm (rho)
    blaze_angle: 44              # Blaze angle in degrees
    eps_angle: 4                 # Epsilon angle in degrees
    ilg: 64.713383393226500      # ILG parameter
    n_p: 585                     # Number of pixels
    orders_table_file: "path/to/order/table.dat"  # Path to orders table file

    # For Binary grating: CUBES example
    type: "Binary"               # Grating type
    line_density: 3.60           # Line density in lines/μm (rho)
    blaze_angle: 36.07           # Blaze angle in degrees
    eps_angle: 0                 # Epsilon angle in degrees
    ilg: 197.4723512055243       # ILG parameter
    n_p: 2900                    # Number of pixels
    orders_table_file: "path/to/table.csv"  # Path to orders table file

  slit_size_x: [0.5, 1, 1.5, 5] # Slit size X in arcsec
  slit_size_y: [12]             # Slit size Y in arcsec
  slit_size_x_calibration: [0.0055, 0.011, 0.0165, 0.055]  # Slit size X calibration in cm
  slit_size_y_calibration: [0.132]  # Slit size Y calibration in cm

  # Efficiency and PSF configuration files
  telescope_fdr_file: "path/to/TELESCOPEfdr.dat"  # Telescope efficiency data file
  instrument_fdr_file: "path/to/INSTRUMENT_fdr.dat"  # Instrument efficiency data file
                                                    # Includes: common path * grating * collimator * 
                                                    # field mirror * cross-disperser * fold mirror * 
                                                    # camera * detector QE
  psf_map_file: "path/to/PSFmap.npy"  # PSF map data file

  n_pixels: 4096                  # Number of detector pixels
  dimension_pixel: 15             # Pixel size in microns
  # Additional spectrograph parameters...

acquisition:
  sed:
    sed_type: "Blackbody"        # Available types: "Blackbody", "Powerlaw", "Flat", "Lamp", "ThermalRadiation", "Spectrum"
    # Parameters based on sed_type:
    # For "Blackbody":
    band: "V"                    # Photometric band
    magnitude: 0.0               # Magnitude in the specified band
    magnitude_system: "Vega"     # "Vega" or "AB"
    temperature: 5800            # Temperature in Kelvin
    bandpass_normalization: true # Optional: normalize using bandpass

    # For "Powerlaw":
    band: "V"                    # Photometric band
    magnitude: 0.0               # Magnitude in the specified band
    magnitude_system: "Vega"     # "Vega" or "AB"
    index: -2.0                  # Power law index
    bandpass_normalization: true # Optional: normalize using bandpass

    # For "Flat":
    energy: 1.0                  # Energy level

    # For "Lamp":
    spectrum_file: "path/to/spectrum.dat"  # Path to lamp spectrum file

    # For "ThermalRadiation":
    temperature: 273.15          # Temperature in Kelvin
    Fn: 8.0                      # F-number
    filter: "path/to/filter.dat" # Path to filter file

    # For "Spectrum":
    band: "V"                    # Photometric band
    magnitude: 0.0               # Magnitude in the specified band
    magnitude_system: "Vega"     # "Vega" or "AB"
    spectrum_file: "path/to/spectrum.dat"  # Path to spectrum file
    bandpass_normalization: true # Optional: normalize using bandpass
    z: 0.0                      # Optional: redshift

  characteristics:
    slit_size_x_index: 0         # Index for slit size X
    slit_size_y_index: 0         # Index for slit size Y
    detector_integration_time: 1200  # Exposure time in seconds
    number_of_dit: 1             # Number of exposures
    number_of_integration: 1      # Number of integrations
    bin_x: 1                     # Binning in X
    bin_y: 1                     # Binning in Y

  sky:
    airmass: 1.16                # Atmospheric airmass
    pwv: 10                      # Precipitable water vapor
    mfli: 0                      # Moon fractional lunar illumination
    seeing: 0.8                  # Seeing in arcseconds

simulation:
  pixel_oversampling: 5          # Pixel oversampling factor
  wavelength_overscanning: 40   # Wavelength oversampling
  mask_oversampling: 5           # Mask oversampling
  psf_pupil_sampling: 128       # PSF pupil sampling
  psf_field_sampling: 64         # PSF field sampling
  psf_map_pixel_number: 10       # PSF map pixel number
  orders_index: [3, 4, 5]        # Optional: specific orders to simulate

zemax:                           # Optional: Zemax integration
  file_path: "path/to/zemax/file"
  order_table_flag: true
  PSF_map_flag: true
  CFG_path: "path/to/cfg/file"
```

## Output

The simulator generates:
- FITS files containing the simulated detector image

## Notes

- Processing time depends on the complexity of the simulation and the chosen parameters
- Memory usage can be significant with high oversampling values
- For large simulations, parallel processing can be enabled by setting `PARALLEL = True` in the code

## Contributing

Feel free to submit issues and enhancement requests.
