from typing import Tuple

import numpy as np
import scipy

def great_circle_pairwise(
    longitude_a: np.ndarray,
    latitude_a: np.ndarray,
    longitude_b: np.ndarray,
    latitude_b: np.ndarray,
    earth_radius: float = 6378.137,
    mod_bearing: bool = True
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Computes the great circle distance (km) and true fore bearing (deg) between
    pairs of observations in input arrays `longitude_a` and `longitude_b` and
    `latitude_a` and `latitude_b`.

    For two longitude and latitude pairs, the great circle distance is the
    shortest distance between the two points along the Earth's surface. This
    distance is calculated using the Haversine formula. The instances in
    `longitude_a` and `latitude_a` are designated as point `a`; the instances
    in `longitude_b` and `latitude_b` then form point `b`. The true fore
    bearing is the bearing, measured from true north, of `b` as seen from `a`.

    Note:
        When given `latitude_a/b` and `longitude_a/b` of shape (n,), n > 1,
        the great circle distance and fore bearing will be calculated between
        `a` and `b` entries such that the returned arrays will be of shape
        (n,). To compute the great circle distance and bearings between
        adjacent coordinates of single longitude and latitude arrays (i.e.,
        along a trajectory), use `great_circle_pathwise`.

    Args:
        longitude_a (np.array): of shape (n,) in units of decimal degrees
        latitude (np.array): of shape (n,) in units of decimal degrees
        earth_radius (float, optional): earth's radius in units of km. Defaults to 6378.137 km (WGS-84)
        mod_bearing (bool, optional): return bearings modulo 360 deg. Defaults to True.

    Returns:
        Tuple[np.array, np.array]: great circle distances (in km) and true fore
        bearings between adjacent longitude and latitude pairs; shape (n,)
"""
    # Convert decimal degrees to radians
    longitude_a_rad, latitude_a_rad = map(np.radians, [longitude_a, latitude_a])
    longitude_b_rad, latitude_b_rad = map(np.radians, [longitude_b, latitude_b])

    # Difference longitude and latitude
    longitude_difference = longitude_b_rad - longitude_a_rad
    latitude_difference = latitude_b_rad - latitude_a_rad

    # Haversine formula
    a_1 = np.sin(latitude_difference / 2) ** 2
    a_2 = np.cos(latitude_a_rad)
    a_3 = np.cos(latitude_b_rad)
    a_4 = np.sin(longitude_difference / 2) ** 2
    c = 2 * np.arcsin(np.sqrt(a_1 + a_2 * a_3 * a_4))
    distance_km = earth_radius * c

    # True bearing
    bearing_num = np.cos(latitude_b_rad) * np.sin(-longitude_difference)
    bearing_den_1 = np.cos(latitude_a_rad) * np.sin(latitude_b_rad)
    bearing_den_2 = - np.sin(latitude_a_rad) * np.cos(latitude_b_rad) * np.cos(longitude_difference)
    bearing_deg = -np.degrees(np.arctan2(bearing_num, bearing_den_1 + bearing_den_2))

    if mod_bearing:
        bearing_deg = bearing_deg % 360

    return distance_km, bearing_deg


def match_model_and_buoy_by_interpolation(
    buoy_time,
    buoy_longitude,
    buoy_latitude,
    model_time,
    model_field,
    model_longitude,
    model_latitude,
    temporal_tolerance: np.timedelta64 = np.timedelta64(30, 'm'),
    **interpn_kwargs,
):
    """
    Match model and buoy observations using linear interpolation in time
    and bilinear interpolation in space.

    Note: the `buoy_time` and `model_time` arrays must be sorted.

    Args:
        buoy_time (np.array[datetime64]): buoy observation times.
        buoy_longitude (np.array[float]): buoy longitudes.
        buoy_latitude (np.array[float]): buoy latitudes.
        model_time (np.array[datetime64]): model observation times.
        model_field (np.array): model field variable to be matched onto
            the buoy coordinates.
        model_longitude (np.array[float]): model longitudes.
        model_latitude (np.array[float]): model latitudes.
        temporal_tolerance (np.timedelta64, optional): maximum allowable
        time difference between a model and observation point. Defaults
        to np.timedelta64(30, 'm').
        **interpn_kwargs: Remaining keyword arguments passed to scipy.interpn.

    Returns:
        np.ndarray: interpolated field values for each time in `buoy_time`.
    """
    if 'method' not in interpn_kwargs:
        interpn_kwargs['method'] = 'linear'
    if 'bounds_error' not in interpn_kwargs:
        interpn_kwargs['bounds_error'] = True

    t_sort_indices = np.searchsorted(model_time, buoy_time)

    # Adjust the sort indices so that the final index is not greater
    # than the length of the array.  If so, replace with the last index.
    n = model_time.size
    t_sort_indices[t_sort_indices >= n] = n - 1

    field_matches = []

    points = (model_latitude, model_longitude)
    for i, j in enumerate(t_sort_indices):

        time_difference = np.abs(buoy_time[i] - model_time[j])

        if time_difference > temporal_tolerance:
            value = np.nan
        else:
            x_i = (buoy_latitude[i], buoy_longitude[i])

            field_values_jm1 = model_field[j-1]  # left
            field_values_j = model_field[j]  # right

            bilinear_value_jm1 = scipy.interpolate.interpn(points,
                                                           field_values_jm1,
                                                           x_i,
                                                           **interpn_kwargs)

            bilinear_value_j = scipy.interpolate.interpn(points,
                                                         field_values_j,
                                                         x_i,
                                                         **interpn_kwargs)

            value = np.interp(buoy_time[i].astype("float"),
                              np.array([model_time[j-1],
                                        model_time[j]]).astype("float"),
                              np.concatenate([bilinear_value_jm1,
                                              bilinear_value_j]))

        field_matches.append(value)

    return np.array(field_matches)
