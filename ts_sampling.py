from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats
from gnssts import SMALL_SIZE,BIGGER_SIZE,load_data

plt.rc("font", size=SMALL_SIZE)  # controls default text sizes
plt.rc("axes", titlesize=SMALL_SIZE)  # fontsize of the axes title
plt.rc("axes", labelsize=SMALL_SIZE)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc("ytick", labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc("legend", fontsize=SMALL_SIZE)  # legend fontsize
plt.rc("figure", titlesize=BIGGER_SIZE)  # fontsize of the figure title

def plot_sample(
        station:str,
        t:np.array,
        u:np.array,
        intervals:int=4,
        samples:int=100,
        criterion:float=1.0,
        alpha:float=0.05,
        state:str="u"
    ):
    """
    Sample the data and plot the mean of the samples with confidence interval
    and the largest outlier of the samples.
    """
    u_gradient, _, _, _, _ = stats.linregress(t, u)
    sample_fit,sample_gradient=[],[]
    for _ in range(samples):
        sample = sample_data(data[station],intervals)
        if state == "e":
            gradient, intercept, _, _, _ = stats.linregress(
                sample["year"], sample["U_e"]
            )
        else:
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
    plt.title(f"{station} {N_above/len(sample_gradient):.2f}% samples above gradient difference criterion of {criterion:.2f} mm/yr with sample size {intervals}")
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
        sample_fit[max],
        color="lightskyblue",
        label="Biggest outlier sample"
    )
    

    plt.grid()
    plt.xlabel("Time [year]", fontsize=12)
    plt.ylabel("Up [mm]", fontsize=12)
    if state == "e": plt.ylabel("Easting [mm]", fontsize=12)
    plt.xlim([2002, 2020])
    plt.ylim([-20, 40])
    if state == "e": plt.ylim([2,6])
    plt.legend(loc="lower right")
    print(f"{station} {N_above/len(sample_gradient):.2f}% above criterion of {criterion:.2f} mm/yr")


def plot_intervals(
        station:str,
        t:np.array,
        u:np.array,
        intervals:np.array=np.arange(4,11,dtype=int),
        samples:int=100,
        alpha:float=0.05
    ):
    """
    Sample different sized samples and plot (1-alpha)% confidence interval
    for each sample size around the mean
    """
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
            label=f"Gaussian {100*(1-alpha):.2f}% confidence interval at sample size {intv}"
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


def sample_data(data:np.array, N:int):
    """
    Separate data in N time intervals and sample a point from each interval
    """
    sample = []
    for subdata in np.array_split(data,N):
        sample.append(np.random.choice(subdata))
    return np.array(sample)

def base_stats(fit:list, N_samples:int, alpha:float=0.05):
    """
    Compute the mean, std and confidence interval for the sample fit
    """
    fit = np.array(fit)
    mean = fit.mean(axis=0)
    std = fit.std(axis=0)

    t = stats.t.ppf(1-alpha/2, N_samples)
    ci_lower = mean - t*std/np.sqrt(N_samples)
    ci_upper = mean + t*std/np.sqrt(N_samples)

    return mean,std,ci_lower,ci_upper

def weeklify(t:np.array,u:np.array,ue:np.array):
    """
    return weekly average of time series
    """
    excess = len(u) % 7 # days beyond last week
    if excess > 0:
        u = u[:-excess]
        t = t[:-excess]
        ue = ue[:-excess]

    u = np.mean(u.reshape(-1,7),axis=1)
    ue = np.mean(ue.reshape(-1,7),axis=1)
    t = np.mean(t.reshape(-1,7),axis=1)
    return t,u,ue



if __name__ == "__main__":

    files = Path(r"data/GPS").glob("*.txt")
    Path("out/sampling/single").mkdir(exist_ok=True,parents=True)
    Path("out/sampling/multiple").mkdir(exist_ok=True,parents=True)
    Path("out/sampling/easting").mkdir(exist_ok=True,parents=True)
    data = {}
    for filename in files:
        ts, station = load_data(filename)
        data[station] = ts

    
    # plot timeseries
    mpl.rc('figure', max_open_warning = 0)
    for station in data:
        u = data[station]["U"]
        t = data[station]["year"]
        ue = data[station]["U_e"]

        t,u,ue = weeklify(t,u,ue)

        # Plot sample with "intervals" sample size sampled "samples" times
        plot_sample(station, t, u, intervals=10, samples=1000, criterion=1, alpha=0.01)
        plt.savefig(Path("out/sampling/single") / Path(f"{station}_sampling.png"), bbox_inches="tight")

        # Plot samples with multiple sample sizes
        plot_intervals(station, t, u, intervals=np.arange(3,13,dtype=int), samples=1000, alpha=0.05)
        plt.savefig(Path("out/sampling/multiple") / Path(f"{station}_variable_interval_sampling.png"),bbox_inches="tight")

        # Same as first plot for easting instead of uplift
        plot_sample(station, t, ue, intervals=4, samples=1000, criterion=.01, alpha=0.01,state="e")
        plt.savefig(Path("out/sampling/easting") / Path(f"{station}_sampling_Ue.png"), bbox_inches="tight")
