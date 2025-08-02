
from os import PathLike
import os
from scipy.io import loadmat
import numpy as np
import numpy.typing as npt


def load_christians_alaska_model(load_dir: PathLike) -> tuple(dict[npt.NDArray, tuple(int, int,int)]):
    data = {}

    for file, field in zip(['S_vel_to_share.mat', 'P_vel_to_share.mat', 'S_atten_to_share.mat'],
                        ['grd_s_vel', 'grd_p_vel', 'grd_atten_s']):
        fname = os.path.join(load_dir, file)

        vels = loadmat(fname,
                    struct_as_record=False,
                    squeeze_me=True,)
        data[field] = vels[field].val


    shp = data[field].shape
    return data, shp


