# gnssts

Filtrer og plot tidsserier fra GNSS danske permanente stationer.

## Installationsvejledning

Klon git repository og hop i repositoriet:

```
$ git clone https://github.com/Kortforsyningen/gnssts.git
$ cd gnssts
```

Opret miljø med conda:

```
$ conda env create --file environment.yml
```

Aktivér conda miljø:

```
$ conda activate gnssts
```

## Kør programmet

Programmet er todelt. Med kommandoen

```
$ python gnssts.py
```

afvikles gnssts programmet, som viser tidsserier for hver station med moving average og et linear fit. Samtidigt udføres en de-noising af data og de samme parametre der blev vist for den rå tidsserie bliver også vist for disse.

herefter outputtes de filtrerede tidsserier i mappen `gnssts/out/gnssts`.

Køres herimod kommandoen

```
$ python ts_sampling.py
```

afvikles ts_sampling programmet som først tager 1-uges gennemsnit af tidsserierne. Herefter samples disse for at undersøge hvor mange samples der skal til fra en station, før man kan antage at data er pålidelige. Programmet outputter 3 plots for hver station:

1. Plot med gennemsnitlig hældning af samples samt konfidensinterval for denne for en sample størrelse. I titlen vises hvor stor en andel af samples der er over en selvvalgt grænseværdi. Samtidig vises det største outlier sample. Disse plots outputtes til mappen `gnssts/out/sampling/single`
2. Plot med gennemsnitlig hældning af samples med konfidensintervaller for denne for flere sample størrelser. Disse plots outputtes til mappen `gnssts/out/sampling/multiple`
3. Samme plot som det første, men for easting i stedet for up. Disse plots outputtes til mappen `gnssts/out/sampling/easting`

Repositoriet kommer med data fra 2019 som standard. Data fra 2019 befinder sig i
mappen `data/GPS/` og kan blot overskrives med filer i samme format hvis opdaterede
tidsserier skal filtreres og plottes.

Som set i plot 2, ser det ud til at ved en sample size større 7-8 er der ikke længere større ændringer i 99% konfidensintervallet, mens det for et 95% konfidensinterval er ved sample size 5-7. Derfor anbefales det som absolut minimum at anvende 5 målinger til fremtidige beregninger, mens en større sikkerhed vil opnås ved mere end 7 målinger pr. station.


