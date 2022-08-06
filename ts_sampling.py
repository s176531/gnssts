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

def plot_sample(station, t, u, intervals=4, samples=100, criterion=1, alpha=0.05):
    
    u_gradient, _, _, _, _ = stats.linregress(t, u)
    sample_fit,sample_gradient=[],[]
    for _ in range(samples):
        sample = sample_data(data[station],intervals)
        gradient, intercept, _, _, _ = stats.linregress(
            sample["year"], sample["U"]
        )
        sample_fit.append(t*gradient + intercept)
        sample_gradient.append(gradient)
    
    mean,std,ci_lower,ci_upper = base_stats(
        sample_fit,
        intervals,
        alpha=alpha,
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
        "-",
        linewidth=.75,
        color="indigo",
        label=f"Sample mean ({samples} samples)"
    )
    plt.plot(
        t,
        ci_lower,
        "--",
        linewidth=1,
        color="darkorchid",
        label=f"Gaussian {100*(1-alpha):.2f}% confidence interval"
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


def plot_intervals(station, t, u, intervals=np.arange(4,11,dtype=int), samples=100, alpha=0.95):
    
    plt.figure(figsize=(14,8))
    plt.plot(
            t,
            u,
            ".",
            markersize=2,
            linewidth=0.25,
            color="lightsalmon",
            label="Data points",
        )
    for intv in intervals:
        u_gradient, _, _, _, _ = stats.linregress(t, u)
        sample_fit,sample_gradient=[],[]
        for _ in range(samples):
            sample = sample_data(data[station],intv)
            gradient, intercept, _, _, _ = stats.linregress(
                sample["year"], sample["U"]
            )
            sample_fit.append(t*gradient + intercept)
            sample_gradient.append(gradient)
        
        mean,std,ci_lower,ci_upper = base_stats(
            sample_fit,
            intv,
            alpha=alpha,
        )
        denom = np.max(intervals)*1.1
        cicol = (1-intv/denom,1-intv/denom,intv/denom)
        
        plt.plot(
            t,
            ci_lower,
            "--",
            linewidth=1,
            color=cicol,
            label=f"2*standard deviation at {intv} intervals"
        )
        plt.plot(
            t,
            ci_upper,
            "--",
            linewidth=1,
            color=cicol,
        )

    plt.plot(
        t,
        mean,
        "-",
        linewidth=0.75,
        color="crimson",
        label=f"Sample mean"
    )

    plt.grid()
    plt.title(f"{(1-alpha)*100:.2f}% confidence intervals at {intervals[0:]} intervals in sampling")
    plt.xlabel("Time [year]", fontsize=12)
    plt.ylabel("Up [mm]", fontsize=12)
    plt.xlim([2002, 2020])
    plt.ylim([-20, 40])
    plt.legend(loc="upper left")
    print(f"{station} done with multiple intervals")


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

    gradient, intercept, _, _, _ = stats.linregress(
        ts["year"], ts["U"]
    )
    offset = gradient * ts["year"][0] + intercept
    ts["U"] = ts["U"] - offset
    return ts, filename.stem[0:4]

def sample_data(data, N):
    """
    Opdel data i N tidsintervaller og sample et punkt fra hvert interval
    """
    sample = []
    for subdata in np.array_split(data,N):
        sample.append(np.random.choice(subdata))
    return np.array(sample)

def base_stats(fit, N_samples, alpha=0.05):
    fit = np.array(fit)
    mean = fit.mean(axis=0)
    std = fit.std(axis=0)

    t = stats.t.ppf(1-alpha/2, N_samples)
    ci_lower = mean - t*std/np.sqrt(N_samples)
    ci_upper = mean + t*std/np.sqrt(N_samples)

    return mean,std,ci_lower,ci_upper

def weeklify(t,u):
    """
    return weekly average of time series
    """
    excess = len(u) % 7 # days beyond last week
    if excess > 0:
        u = u[:-excess]
        t = t[:-excess]
    u = np.mean(u.reshape(-1,7),axis=1)
    t = np.mean(t.reshape(-1,7),axis=1)
    return t,u



if __name__ == "__main__":

    files = Path(r"data/GPS").glob("*.txt")
    Path("out/sampling/single").mkdir(exist_ok=True)
    Path("out/sampling/multiple").mkdir(exist_ok=True)
    data = {}
    for filename in files:
        ts, station = load_data(filename)
        data[station] = ts

    
    # plot timeseries
    mpl.rc('figure', max_open_warning = 0)
    for station in data:
        u = data[station]["U"]
        t = data[station]["year"]

        t,u = weeklify(t,u)

        plot_sample(station, t, u, intervals=4, samples=1000, criterion=1, alpha=0.01)
        plt.savefig(Path("out/sampling/single") / Path(f"{station}_sampling.png"), bbox_inches="tight")

        plot_intervals(station, t, u, intervals=np.arange(3,16,dtype=int), samples=1000, alpha=0.01)
        plt.savefig(Path("out/sampling/multiple") / Path(f"{station}_variable_interval_sampling.png"),bbox_inches="tight")

