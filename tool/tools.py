from scipy import interpolate


def interp(x, y, new_x, fill_value, kind="linear"):
    function = interpolate.interp1d(x, y, kind=kind, fill_value=fill_value)
    return function(new_x)
