from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats
from scipy import interpolate

SMALL_SIZE = 9
MEDIUM_SIZE = 12
BIGGER_SIZE = 20


plt.rc("font", size=MEDIUM_SIZE)  # controls default text sizes
plt.rc("axes", titlesize=MEDIUM_SIZE)  # fontsize of the axes title
plt.rc("axes", labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc("ytick", labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc("legend", fontsize=SMALL_SIZE)  # legend fontsize
plt.rc("figure", titlesize=BIGGER_SIZE)  # fontsize of the figure title

def plot_up(station:str, t:np.array, u:np.array, noise_func, N:int=1):

    silent = u - noise_func(t)

    u_moving_avg = np.convolve(u, np.ones((N,)) / N, mode="same")
    silent_moving_avg = np.convolve(silent, np.ones((N,)) / N, mode="same")

    u_gradient, u_intercept, u_r, u_p, u_stderr = stats.linregress(t, u)
    silent_gradient, silent_intercept, silent_r, silent_p, silent_stderr = stats.linregress(
        t, silent
    )
    u_fit = u_gradient * t + u_intercept
    silent_fit = silent_gradient * t + silent_intercept

    plt.figure(figsize=(14,8))
    plt.title(station)
    plt.plot(
        t,
        u,
        ".",
        markersize=2,
        linewidth=0.25,
        color="lightsalmon",
        label="Data points",
    )
    plt.plot(
        t,
        silent,
        "+",
        markersize=2,
        linewidth=0.25,
        color="lightskyblue",
        label="Data points, de-noised",
    )
    plt.plot(
        t,
        u_moving_avg,
        "-",
        markersize=2,
        linewidth=1,
        color="indianred",
        label=f"Moving average (N={N})",
    )
    plt.plot(
        t,
        silent_moving_avg,
        "-",
        markersize=2,
        linewidth=1,
        color="royalblue",
        label=f"Moving average, de-noised (N={N})",
    )
    plt.plot(
        t,
        u_fit,
        linewidth=1.5,
        color="red",
        label=f"Linear fit, {u_gradient:.2f} mm/yr, r={u_r:.2f}, std_dev={u_stderr:.2f}",
    )
    plt.plot(
        t,
        silent_fit,
        "-.",
        linewidth=1.5,
        color="navy",
        label=f"Linear fit, de-noised, {silent_gradient:.2f} mm/yr, r={silent_r:.2f}, std_dev={silent_stderr:.2f}",
    )

    plt.grid()
    plt.xlabel("Time [year]", fontsize=12)
    plt.ylabel("Up [mm]", fontsize=12)
    plt.xlim([2002, 2020])
    plt.ylim([-20, 40])
    plt.legend(loc="lower right")
    print(f"{station}  {silent_gradient:.2f} mm/yr +/- {silent_stderr:.2f} mm/yr")


def load_data(filename: Path):
    """
    Load the data into a numpy array consisting of tuples.
    Each tuple contains: (decimal year,up,easting) of a reading from a station.
    Return the array as well as a 4 letter identifier of the station.
    """

    with open(filename, "r") as f:
        timeseries = f.readlines()

    try: # Is the data in NEU format?
        ts = np.loadtxt(
            timeseries,
            comments="%",
            dtype={
                "names": ("year", "N", "N_e", "E", "E_e", "U", "U_e"),
                "formats": ("f4", "f4", "f4", "f4", "f4", "f4", "f4"),
            },
        )
    except (IndexError, ValueError): # ... or in year/up format?
        ts = np.loadtxt(
            timeseries,
            comments="%",
            dtype={"names": ("year", "U", "U_e"), "formats": ("f4", "f4", "f4")},
        )

    gradient, intercept, _, _, _ = stats.linregress(
        ts["year"], ts["U"]
    )
    offset = gradient * ts["year"][0] + intercept
    ts["U"] = ts["U"] - offset
    return ts, filename.stem[0:4]


def detrend(t:np.array, u:np.array):
    gradient, intercept, r_value, p_value, std_err = stats.linregress(t, u)
    u_trend = gradient * t + intercept

    return u - u_trend


def find_avg_signal(data:dict, stations:tuple):
    """
    1. Interpoler rå uplift signal, så datapunkter falder sammen i alle tidsserier (samme tidspunkter)
    2. Detrend uplift-signaler
    3. Find gennemsnit af signaler - det giver et udtryk for den generelle støj
    4. Return interpolationsfunktion til støj-signalet
    """
    t_min = 999999.9
    t_max = 0.0
    for station in stations:
        t_min = min(t_min, np.min(data[station]["year"]))
        t_max = max(t_max, np.max(data[station]["year"]))
    t = np.arange(t_min, t_max, 1 / 365.25)

    ts = []
    for station in stations:
        f = interpolate.interp1d(
            data[station]["year"],
            data[station]["U"],
            bounds_error=False,
            fill_value=(data[station]["U"][0], data[station]["U"][-1]),
        )
        u = f(t)
        u = detrend(t, u)

        ts.append(u)

    tss = np.vstack(ts)
    ts_avg = np.average(tss, axis=0)

    return interpolate.interp1d(
        t, ts_avg, bounds_error=False, fill_value=(ts_avg[0], ts_avg[-1])
    )


if __name__ == "__main__":

    files = Path(r"data/GPS").glob("*.txt")
    Path("out/gnssts").mkdir(exist_ok=True)
    data = {}
    for filename in files:
        ts, station = load_data(filename)
        data[station] = ts

    noise_func = find_avg_signal(data, ("BUDP", "SMID", "SULD"))

    # save timeseries
    for station in data:
        with open(Path("out/gnssts") / Path(f"{station}.txt"), 'w') as f:
            f.write("% site name BUDP    component: U\n")
            f.write("% Time [Year],    Uplift [mm]\n")
            for d in data[station]:
                f.write(f"{d['year']:.8f}    {d['U']-noise_func(d['year']): .3f}\n")
    
    # plot timeseries
    mpl.rc('figure', max_open_warning = 0)
    for station in data:
        u = data[station]["U"]
        t = data[station]["year"]
        plot_up(station, t, u, noise_func, N=30)
        plt.savefig(Path("out/gnssts") / Path(f"{station}.png"), bbox_inches="tight")


