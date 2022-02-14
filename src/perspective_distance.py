# -*- coding: utf-8 -*-

"""Raster Distance Calculator

Requirements:
    - numpy
    - gdal

This script allows the user to take an elevation raster and a reference 
point and create a new raster with distances from each cell in the 
elevation raster to the reference point.

This is done using the PerspectiveDistance class, which has the following
functions

    * PerspectiveDistance.qgis_3d_view_crs_to_crs()
        [static method]
        Convert camera coordinates from a QGIS 3D view to real world 
        coordinates.

    * PerspectiveDistance.dist_from_point()
        [static method]
        Calculate the distances between each point in an array to a
        reference point using numpy.

    * PerspectiveDistance().set_point()
        Set the x, y, z coordinates to use for calculating distances
        using this instance of the PerspectiveDistance class. 
        Coordinates must be in the same CRS as the input raster(s).

    * PerspectiveDistance().calc_perspective_difference()
        Process an input raster against set reference point.

Example:

    >>> # Import PerspectiveDistance class
    >>> from perspective_distance import PerspectiveDistance

    >>> # Input camera coordinates
    >>> Cxyz = (122848, 7704.62, -291798)

    >>> # Real world camera coordinates
    >>> xyz = PerspectiveDistance.qgis_3d_view_crs_to_crs(*Cxyz)

    >>> # Initialize PerspectiveDistance class
    >>> PD = PerspectiveDistance()

    >>> # Set real world coordinates
    >>> PD.set_point(*xyz)

    >>> # Definte input/output raster filenames
    >>> input_raster = "~/DEM.vrt"

    >>> output_raster = input_raster.replace(".vrt", "_distance-from-camera.tif")

    >>> # Process raster
    >>> PD.calc_perspective_difference(
    >>>     input_raster,
    >>>     output_raster,
    >>>     scale=True,
    >>> )

To do:
    - Fix bug with setting NoDataValue on output file.
    - Add size check and raise an error if input file is too large.
    - Create base class and add:
        - light gdal.Translate wrapper
        - method for processing raster line by line/tile by tile.
    - Create class for calulating distances from an orthographic view.

"""

import os

import numpy as np

from os import PathLike
from osgeo import gdal

gdal.UseExceptions()

################################################################################

__author__ = "James Ford"
__copyright__ = "Copyright (c) 2022 James Ford"
__credits__ = ["James Ford"]
__license__ = """MIT License
Copyright (c) 2022 James Ford
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
__version__ = "0.1.0"
__maintainer__ = "James Ford"
__email__ = "irvine.ford@gmail.com"
__status__ = "Dev"

################################################################################


class PerspectiveDistance:
    """Class for calulating the distance from a point to each cell of a raster."""

    ############################################################################

    def __init__(self):
        """Initialize an instance of PerspecitveDistance class."""

        pass

    ############################################################################

    @staticmethod
    def qgis_3d_view_crs_to_crs(
        cx: float, cy: float, cz: float, x_origin: float, y_origin: float,
    ):
        """Convert camera coordinates from QGIS 3D view to real world coordinates.

        Args:
            cx (float): Camera x coordinate.
            cy (float): Camera y coordinate.
            cz (float): Camera z coordinate (is returned unchanged).
            x_origin (float, optional): Origin of X coordinate in destination CRS.
            y_origin (float, optional): Origin of Y coordinate in destination CRS.


        Returns:
            tuple(float, float, float): A tuple of x, y, z float coordinates.
        """

        return (x_origin + cx, y_origin - cz, cy)

    ############################################################################

    @staticmethod
    def dist_from_point(a: np.array, x: float, y: float, z: float) -> np.array:
        """Calculate the euclidean distance from a point (x,y,z) for each point in an array

        Args:
            a (np.array): Array of values to calulate distance to.
            x (float): X coordinate of point to calculate distance from.
            y (float): Y coordinate of point to calculate distance from.
            z (float): Z coordinate of point to calculate distance from.

        Returns:
            np.array: Array of distances from xyz point to each value in the input array.
        """

        # Array of reference point coordinates
        b = np.array([x, y, z])

        # Choose which axis of input array to calculate distances on
        _axis = len(a.shape) - 1

        # Calculate difference between arrays
        AB = a - b

        # Calculate distance from each point in array a to the reference point
        distance = np.linalg.norm(AB, axis=_axis)

        return distance

    ############################################################################

    @staticmethod
    def xyz_array_from_raster(raster_src: gdal.Dataset, band_no: int = 1) -> np.array:
        """[summary]

        Args:
            raster_src (gdal.Dataset): The raster dataset to use. Output of gdal.Open(dataset_filename).
            band_no (int, optional): The raster band to read from. Defaults to 1.

        Returns:
            np.array: A 3d numpy array of 3 stacked 2d arrays, x coordinates, y coordinates, as band values.
        """

        # Default to band number 1 if band_no not with acceptable range
        if band_no not in range(1, raster_src.RasterCount + 1):
            band_no = 1

        # Get band object from gdal dataset
        band = raster_src.GetRasterBand(band_no)

        # Get band as a numpy array
        band_array = band.ReadAsArray()

        # Get and unpack raster geotransform to upper left coordinate, skew, and resolution
        lx, rx, sx, uy, sy, ry = raster_src.GetGeoTransform()

        # Get raster dimensions
        w = raster_src.RasterXSize
        h = raster_src.RasterYSize

        # Calculate upper x and lower y coordinates.
        ux = lx + w * rx
        ly = uy + h * ry

        # Create 1d arrays of x & y coordinates. Y coordinates need to be inverted.
        x_range = np.linspace(lx, ux, w)
        y_range = np.linspace(ly, uy, h)[::-1]

        # Create 2d arrays of x & y coordinates
        x_array, y_array = np.meshgrid(x_range, y_range)

        # Stack x, y, and raster band 2d arrays
        stacked_arrays = np.array([x_array, y_array, band_array])

        # Zip into a 2d array of (x, y, value) tuples
        zipped_array = np.dstack(stacked_arrays)

        return zipped_array

    ############################################################################

    def set_point(self, x: float, y: float, z: float):
        """Set the reference point for distance functions.
        
        
        Coordinates should be in the same CRS as the raster.

        Args:
            x (float): X coordinate of point.
            y (float): Y coordinate of point.
            z (float): Z coordinate of point.
        """

        setattr(self, "x", x)
        setattr(self, "y", y)
        setattr(self, "z", z)

    ############################################################################

    def calc_perspective_difference(
        self,
        input_raster: PathLike,
        output_raster: PathLike,
        band_no=1,
        scale=False,
        output_format="GTIFF",
        creation_options=["COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES"],
    ):
        """Calculate a perspective distance raster from an input raster.


        Perspective point must be set first using PerspectiveDistance.set_point(x, y, z)
        or PerspectiveDistance(x, y, z)

        Writes out to a geotiff by default.

        Args:
            input_raster (PathLike): Input raster filename to calculate distance to from perspective point.
            output_rater (PathLike): Output raster filename
            band_no (int, optional): Band number of input raster to process. Defaults to 1.
            scale (bool, optional): Whether or not to scale distance values between 0 - 1. Defaults to False.
            output_format (str, optional): Raster format of output file. Defaults to GTIFF.
            creation_options (list, optional): A list of valid creation options for raster format.
                Defaults to creation_options=[
                    "COMPRESS=DEFLATE",
                    "PREDICTOR=3",
                    "TILED=YES"
                ]

        """

        # Verify input raster exist
        if not os.path.exists(input_raster):
            raise ValueError(f"Input raster '{input_raster}' does not exists.")

        # Verify point exists
        if not all([hasattr(self, v) for v in ("x", "y", "z")]):
            raise AttributeError(
                f"PerspectiveDistance object must have a point set before '.calc_perspective_difference()' can be called."
            )

        # Create gdal.Dataset object in read-only mode from input raster filename
        src = gdal.Open(input_raster, 0)

        # Get (x, y, value) array from raster
        value_array = PerspectiveDistance.xyz_array_from_raster(src, band_no=band_no)

        # Calculate distance array
        distance_array = PerspectiveDistance.dist_from_point(
            value_array, self.x, self.y, self.z
        )

        # Scale distances
        if scale:
            d_min = np.min(distance_array)
            d_max = np.max(distance_array)
            distance_array -= d_min
            distance_array /= d_max - d_min

        # Get raster format driver
        driver = gdal.GetDriverByName(output_format)

        # Write out to output file
        dst = driver.Create(
            output_raster,  # output filename
            src.RasterXSize,  # width
            src.RasterYSize,  # heigh
            1,  # number of bands
            gdal.GDT_Float32,  # data type
            options=creation_options,  # creation options
        )

        # Set geotransform and projection
        dst.SetGeoTransform(src.GetGeoTransform())
        dst.SetProjection(src.GetProjection())

        # Write data to output band
        output_band = dst.GetRasterBand(1)
        _ = output_band.WriteArray(distance_array)

        # Set output band no data value
        # output_band.SetNoDataValue(src.GetRasterBand(band_no).GetNoDataValue())

        # Write output raster to file
        dst = None


################################################################################


def main():
    pass


if __name__ == "__main__":
    main()
