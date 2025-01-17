from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats
from scipy import interpolate

SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12


plt.rc("font", size=SMALL_SIZE)  # controls default text sizes
plt.rc("axes", titlesize=SMALL_SIZE)  # fontsize of the axes title
plt.rc("axes", labelsize=SMALL_SIZE)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc("ytick", labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc("legend", fontsize=SMALL_SIZE)  # legend fontsize
plt.rc("figure", titlesize=BIGGER_SIZE)  # fontsize of the figure title

def plot_sample(station, t, u, intervals=4, samples=100, criterion=1):
    
    u_gradient, u_intercept, _, _, _ = stats.linregress(t, u)
    sample_fit,sample_gradient,sample_intercept=[],[],[]
    for _ in range(samples):
        sample = sample_data(data[station],intervals)
        gradient, intercept, _, _, _ = stats.linregress(
            sample["year"], sample["U"]
        )
        sample_fit.append(t*gradient + intercept)
        sample_gradient.append(gradient)
        sample_intercept.append(intercept)
    
    alpha=0.99
    mean,std,ci_lower,ci_upper = base_stats(
        sample_fit,
        samples,
        alpha=alpha,
        distribution=0
    )

    dif=np.abs(sample_gradient-u_gradient)
    max = np.argmax(dif)
    N_above = np.sum(dif>criterion)

    plt.figure(figsize=(14,8))
    plt.title(f"{station} {N_above/len(sample_gradient):.2f}% samples above gradient difference criterion of {criterion:.1f} mm/yr")
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
        mean,
        "--",
        linewidth=1.5,
        color="indigo",
        label=f"Sample mean ({samples} samples)"
    )
    plt.plot(
        t,
        ci_lower,
        "--",
        linewidth=1,
        color="darkorchid",
        label=f"Gaussian {100*alpha:.2f}% confidence interval"
    )
    plt.plot(
        t,
        ci_upper,
        "--",
        linewidth=1,
        color="darkorchid",
    )
    plt.plot(
        t,
        mean-2*std,
        "--",
        linewidth=1,
        color="orchid",
        label="2*sample std"
    )
    plt.plot(
        t,
        mean+2*std,
        "--",
        linewidth=1,
        color="orchid",
    )
    plt.plot(
        t,
        sample_fit[max],
        color="lightskyblue",
        label="Biggest outlier sample"
    )
    

    plt.grid()
    plt.xlabel("Time [year]", fontsize=12)
    plt.ylabel("Up [mm]", fontsize=12)
    plt.xlim([2002, 2020])
    plt.ylim([-20, 40])
    plt.legend(loc="lower right")
    print(f"{station} {N_above/len(sample_gradient):.2f}% above criterion of {criterion:.1f} mm/yr")


def plot_up(station, t, u, noise_func, N=1):

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

    gradient, intercept, r_value, p_value, std_err = stats.linregress(
        ts["year"], ts["U"]
    )
    offset = gradient * ts["year"][0] + intercept
    ts["U"] = ts["U"] - offset
    return ts, filename.stem[0:4]


def detrend(t, u):
    gradient, intercept, r_value, p_value, std_err = stats.linregress(t, u)
    u_trend = gradient * t + intercept

    return u - u_trend


def find_avg_signal(data, stations):
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

def sample_data(data, N):
    """
    Opdel data i N tidsintervaller og sample et punkt fra hvert interval
    """
    sample = []
    for subdata in np.array_split(data,N):
        sample.append(np.random.choice(subdata))
    return np.array(sample)

def base_stats(fit, N_samples, alpha=0.95, distribution=0):
    fit = np.array(fit)
    mean = fit.mean(axis=0)
    std = fit.std(axis=0)

    if distribution == 0: # Gaussian
        t = stats.t.ppf(1-alpha/2, N_samples)
        ci_lower = mean - t*std/np.sqrt(N_samples)
        ci_upper = mean + t*std/np.sqrt(N_samples)
    elif distribution == 1: # Laplacian
        median = np.median(fit,axis=0)
        decay = np.mean(np.abs(fit-median))
        exp = stats.expon.ppf(1-alpha,N_samples)
        ci_lower = median-((N_samples-1)*decay / exp)
        ci_upper = median+((N_samples-1)*decay / exp)
    return mean,std,ci_lower,ci_upper

if __name__ == "__main__":

    files = Path(r"data/GPS").glob("*.txt")
    Path("out").mkdir(exist_ok=True)
    data = {}
    for filename in files:
        ts, station = load_data(filename)
        data[station] = ts

    noise_func = find_avg_signal(data, ("BUDP", "SMID", "SULD"))

    # save timeseries
    for station in data:
        with open(Path("out") / Path(f"{station}.txt"), 'w') as f:
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
        plt.savefig(Path("out") / Path(f"{station}.png"), bbox_inches="tight")

        plot_sample(station, t, u, intervals=4, samples=1000, criterion=1)
        plt.savefig(Path("out") / Path(f"{station}_sampling.png"), bbox_inches="tight")



