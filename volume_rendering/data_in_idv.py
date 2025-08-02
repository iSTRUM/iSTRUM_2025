from yt_xarray import transformations as tf
import yt_xarray
import yt
import xarray as xr
import os
from scipy.io import loadmat
import yt_idv
import numpy as np

# data dict
load_dir = '/Users/chavlin/data/vbr_related/christian_data'


# get the right limits
lon_x_range = (-170, -130)
lat_y_range = (50, 75)
depth_range = (0, 500)

bbox = np.array([lat_y_range, lon_x_range, depth_range])


data = {}

for file, field in zip(['S_vel_to_share.mat', 'P_vel_to_share.mat', 'S_atten_to_share.mat'],
                       ['grd_s_vel', 'grd_p_vel', 'grd_atten_s']):
    fname = os.path.join(load_dir, file)

    vels = loadmat(fname,
                struct_as_record=False,
                squeeze_me=True,)
    data[field] = vels[field].val


shp = data[field].shape


coords = []
axis_order=('latitude', 'longitude', 'depth')
for iax, ax in enumerate(axis_order):
    rangevals = bbox[iax]
    coords.append(np.linspace(rangevals[0], rangevals[1], shp[iax]))

das = {}
for field, vals in data.items():
    das[field] = xr.DataArray(vals, coords=coords, dims=axis_order)

ds_x = xr.Dataset(data_vars=das)

grid_resolution = (32, 32, 32)
gc = tf.GeocentricCartesian(radial_type='depth', r_o=6371., use_neg_lons=True)
ds_yt = tf.build_interpolated_cartesian_ds(
    ds_x,
    gc,
    fields = ['grd_s_vel', 'grd_p_vel', 'grd_atten_s'] ,
    grid_resolution = grid_resolution,
    refine_grid=True,
    refine_max_iters=2000,
    refine_min_grid_size=4,
    refine_by=4,
    interp_method='interpolate',
)


def _fast_vels_s(field, data):
    vals = data['stream', 'grd_s_vel'].copy()
    vals[vals<0] = 0
    vals[np.isnan(vals)] = 0.0
    return vals

def _slow_vels_s(field, data):
    vals = data['stream', 'grd_s_vel'].copy()
    vals[vals>0] = 0
    vals[np.isnan(vals)] = 0.0
    return np.abs(vals)

ds_yt.add_field(
    name=("stream", "fast_s"),
    function=_fast_vels_s,
    sampling_type="local",
    units="",
    take_log=False,
)

ds_yt.add_field(
    name=("stream", "slow_s"),
    function=_slow_vels_s,
    sampling_type="local",
    units="",
    take_log=False,
)

rc = yt_idv.render_context(height=800, width=800, gui=True)
sg = rc.add_scene(ds_yt, ("stream", "slow_s"), no_ghost=True)
rc.run()
