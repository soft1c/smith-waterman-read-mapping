# Apple CPU profiling: v1.0.0–v1.7.0

## Mjerenje cijelog workload-a

Verzije v1.0.0–v1.6.0 imaju po jedan Release run; v1.7.0 ima pet runova i prikazan je medijan. Svi koriste isti deterministički ulaz na Apple M5 Pro, 18 CPU jezgri, 24 GB RAM-a, macOS 26.7.2, Apple clang 21, `-O3 -DNDEBUG -g`. Ulaz je referenca od 100.000 bp i 800 readova od 160 bp, pet zamjena po readu, seed `21021`; to je ukupno 12,8 milijardi DP ćelija. Bodovanje je `match=2`, `mismatch=-1`, `gap=-2`.

`wall` je ukupno proteklo vrijeme. `DP` i `GCUPS` su iz programa. CPU vrijeme je zbir `user + system` iz `/usr/bin/time -l`; prosječne jezgre računam kao CPU vrijeme / wall. Procenat je udio od 18 dostupnih jezgri. Paralelna verzija može zato imati više CPU sekundi nego wall sekundi.

| Verzija | Režim | Wall (s) | DP (s) | GCUPS | CPU vrijeme (s) | Prosjek jezgri (% od 18) |
|---|---|---:|---:|---:|---:|---:|
| v1.0.0 | `full` | 35,11 | 33,77 | 0,365 | 34,97 | 1,00 (5,5%) |
| v1.1.0 | `rolling` | 39,69 | 39,67 | 0,323 | 39,51 | 1,00 (5,5%) |
| v1.2.0 | `trace` | 20,40 | 20,37 | 0,628 | 20,33 | 1,00 (5,5%) |
| v1.3.0 | `trace` | 21,33 | 21,32 | 0,600 | 21,28 | 1,00 (5,5%) |
| v1.4.0 | `trace` | 15,97 | 15,96 | 0,802 | 15,94 | 1,00 (5,5%) |
| v1.5.0 | `wavefront` | 5,13 | 5,01 | 2,498 | 31,01 | 6,05 (33,6%) |
| v1.6.0 | `wavefront` | 2,46 | 2,40 | 5,211 | 30,14 | 12,25 (68,2%) |
| v1.7.0 | `batch` | 0,95 | 0,94 | 13,553 | 16,09 | 16,93 (94,1%) |

v1.7.0 raspoređuje nezavisne readove na jedan zajednički worker pool, a svaki worker računa svoj read `trace` kernelom. To izbjegava pravljenje wavefront worker grupe za svaki kratki read. Batch medijan je 2,59× brži od prethodnog v1.6.0 runa na ovom ulazu. Peak RSS raste sa 15,1 MiB na 99,2 MiB jer do 18 readova istovremeno drži traceback podatke. Raniji v1.0.0–v1.6.0 brojevi su pojedinačni runovi, pa je ovo praktično poređenje, ne strogo statističko poređenje.

## `auto` dispatcher follow-up

U trenutnom radnom stablu je `auto` podrazumijevani režim. Za jedan read ili batch do jedne četvrtine broja radnika bira wavefront; veće grupe idu kroz jedan batch worker pool. Na 18-thread mašini to znači wavefront za najviše četiri read-a. Wavefront obrađuje te readove jedan za drugim, a svaki koristi više threadova. Sedam single-read i pet batch ponavljanja su pokrenuti na istoj Apple M5 Pro mašini i istim Release zastavicama navedenim iznad. CSV izlaz je bio identičan između `auto` i eksplicitnog `wavefront` ili `batch`, te između `auto` i sekvencijalnog `trace` na single-read ulazu.

| Ulaz | Zatražen režim → izabrani | Program / DP (s) | GCUPS | CPU (s) / prosjek jezgri | Peak RSS / footprint (MiB) | Page faults / reclaims | Instrukcije / ciklusi (milijarde) | IPC |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 1 × 10 kbp read × 100 kbp reference | `auto` → `wavefront` | 0,1042 / 0,1041 | 9,599 | 1,71 / 16,4 (91,2%) | 303,7 / 303,4 | 1 / 19.620 | 38,9 / 7,3 | 5,31 |
| Isti single-read ulaz, kontrola | `trace` | 1,0142 / 1,0142 | 0,986 | 1,01 / 1,0 (5,5%) | 241,2 / 240,9 | 1 / 15.614 | 28,7 / 4,4 | 6,49 |
| 800 × 160 bp read × 100 kbp reference | `auto` → `batch` | 0,9120 / 0,9120 | 14,035 | 15,75 / 17,3 (96,0%) | 99,2 / 98,8 | 1 / 6.539 | 365,1 / 66,0 | 5,53 |

CPU prosjek je CPU vrijeme podijeljeno programskim vremenom; procenat je udio od 18 jezgri. `auto` smanjio je single-read vrijeme 9,7× u odnosu na `trace`, uz 62,5 MiB više peak RSS-a. U batch slučaju automatski odabrani pool već koristi oko 96% dostupnog CPU kapaciteta. Paralelizacija wavefront-a za svaki kratki read ne bi bila korisna za ovaj slučaj: ranije v1.6.0 mjerenje na istom 800-read workload-u trajalo je 2,46 s, dok je batch v1.7.0 trajao 0,94 s.

Prag malog batcha provjeren je na dodatnim ulazima; svi režimi su dali identičan CSV. Medijani su iz pet runova po režimu.

| Ulaz | `auto` izbor | `auto` vrijeme / GCUPS | Batch vrijeme / GCUPS | Ubrzanje | Peak RSS `auto` / batch |
|---|---|---:|---:|---:|---:|
| 4 × 160 bp, reference 100 kbp | `wavefront` | 0,0116 s / 5,508 | 0,0204 s / 3,136 | 1,76× | 14,4 / 20,5 MiB |
| 4 × 512 bp, reference 200 kbp | `wavefront` | 0,0471 s / 8,694 | 0,1097 s / 3,735 | 2,33× | 52,4 / 106,1 MiB |
| 4 × 2 kbp, reference 100 kbp | `wavefront` | 0,0929 s / 8,616 | 0,2079 s / 3,849 | 2,24× | 63,4 / 195,9 MiB |
| 4 × 10 kbp, reference 100 kbp | `wavefront` | 0,3927 s / 10,186 | 1,0268 s / 3,896 | 2,61× | 303,8 / 959,0 MiB |
| 8 × 160 bp, reference 100 kbp | `batch` | 0,0202 s / 6,328 | 0,0200 s / 6,397 | 0,99× | 39,0 / 39,0 MiB |

Kod četiri read-a, wavefront koristi oko 12–17 CPU jezgri prosječno, naspram oko 3–4 za batch. Kod osam kratkih readova batch je oko 1,13× brži od wavefront-a, pa ostaje odabran. Puni CPU, memorijski i VM brojači za 75 runova nalaze se u [auto-dispatch-small-batches.csv](runtimes/auto-dispatch-small-batches.csv).

Swap i block input/output bili su nula u ovim mjerenjima. IPC i page reclaims nisu memorijski stall brojači. `sample` i `/usr/bin/time -l` ne pokazuju čekanje na pristup RAM-u/cache-u; za to su potrebni Instruments ili hardverski performance counteri. Svi pojedinačni podaci su u [auto-dispatch-followup.csv](runtimes/auto-dispatch-followup.csv).

Isti v1.7.0 batch median pod `--threads` sweepom: 4 niti 3,504 s; 8 niti 1,869 s; 12 niti 1,312 s; 16 niti 1,021 s; 18 niti 0,937 s. Zato `--threads 0` bira preporučeni broj od 18 niti na ovoj mašini. [Sirovi podaci sweepa](runtimes/v1.7.0-thread-scaling.csv).

Postojeći `simd` mod je takođe odmjeren na ovom dugom-reference/kratki-read ulazu: 16,80 s i 0,762 GCUPS ([raw `/usr/bin/time -l`](runtimes/simd-control.time.txt)). Ovdje batch trace pristup bolje koristi jezgre.

### Zašto ne pokretati wavefront pool za svaki read istovremeno?

`batch` koristi jedan thread po aktivnom read-u, do broja traženih radnika. Ako su četiri read-a, obično se koriste četiri thread-a; ostali CPU threadovi mogu mirovati. `wavefront` koristi više radnika za pločice jednog read-a. Za četiri duga read-a zato `auto` pokreće jedan wavefront po read-u redom, koristeći skoro sve jezgre i ograničavajući memoriju na jedan traceback/border set.

Pet runova na 4 × 10 kbp ulazu potvrđuju da više wavefront poola istovremeno ne pomaže: jedan proces sa 18 radnika i rednom obradom trajao je 0,410 s, a četiri odvojena procesa sa po četiri radnika trajala su 0,451 s. Zbir maksimalnih RSS-a za paralelne procese bio je 1.212 MiB, naspram 304 MiB za rednu obradu. Batch trace je trajao 1,050 s uz 959 MiB. Jedan wavefront već koristi oko 17 jezgri, pa dodatni pooli povećaju memoriju i ne povećaju brzinu. Zajednički pool koji raspoređuje spremne pločice više readova bio bi zasebna implementacija; ova mjerenja ne pokazuju razlog da bude brži od trenutnog `auto`. [Svi runovi](runtimes/wavefront-concurrency-experiment.csv).

### Dostupni Apple brojači

Sistemski `/usr/bin/powermetrics` je dostupan, ali zahtijeva root. Dodani [`profile_apple_counters.py`](../../scripts/profile_apple_counters.py) snima njegove `tasks` i `cpu_power` uzorke, per-process IPC, frekvenciju i rezultate `/usr/bin/time -l`. Komanda za pokretanje je u README-u. `powermetrics` ne daje L1/L2 cache miss ili memory-stall brojače. `xctrace` postoji kao sistemski stub, ali ne radi s aktivnim `/Library/Developer/CommandLineTools`; za Apple `CPU Counters` instrument treba punu Xcode instalaciju. Nju nisam instalirao jer je ranije odabrano da se Xcode/Instruments ne instalira.

## Memorija i VM brojači

Peak RSS je maksimalni resident set iz kompletnog run-a, pretvoren u MiB. `sample` footprint je trenutni i maksimalni physical footprint iz odvojenog, dvosekundnog stack capture-a. `page faults / reclaims` su `/usr/bin/time -l` brojači. Swap i block input/output bili su nula za sve runove.

| Verzija | Peak RSS (MiB) | `sample` footprint / peak | Page faults / reclaims | Swap / block I/O |
|---|---:|---:|---:|---:|
| v1.0.0 | 370,8 | 370,7M / 370,7M | 3 / 23.978 | 0 / 0 |
| v1.1.0 | 6,9 | 6752K / 6752K | 2 / 618 | 0 / 0 |
| v1.2.0 | 8,5 | 6768K / 6768K | 3 / 749 | 0 / 0 |
| v1.3.0 | 7,7 | 6752K / 6752K | 3 / 664 | 0 / 0 |
| v1.4.0 | 6,9 | 6752K / 6752K | 3 / 613 | 0 / 0 |
| v1.5.0 | 15,0 | 14,3M / 14,5M | 4 / 28.305 | 0 / 0 |
| v1.6.0 | 15,1 | 14,3M / 14,6M | 4 / 28.307 | 0 / 0 |
| v1.7.0 | 99,2 | nema (sample izvještaj nije dao taj podatak) | 1 / 6.537 | 0 / 0 |

## CPU instructions and wait observations

Apple `/usr/bin/time -l` also recorded retired instructions and elapsed cycles. The ratio is instructions/cycle for this process; it is not a memory-stall counter.

| Verzija | Instrukcije (milijarde) | Ciklusi (milijarde) | Instrukcija/ciklus |
|---|---:|---:|---:|
| v1.0.0 | 533,7 | 157,1 | 3,40 |
| v1.1.0 | 651,3 | 171,9 | 3,79 |
| v1.2.0 | 640,5 | 88,2 | 7,26 |
| v1.3.0 | 648,8 | 93,4 | 6,95 |
| v1.4.0 | 432,0 | 71,3 | 6,06 |
| v1.5.0 | 517,9 | 130,6 | 3,96 |
| v1.6.0 | 509,0 | 128,0 | 3,98 |
| v1.7.0 | 365,1 | 69,0 | 5,29 |

Apple `sample` profili su zasebni dvosekundni stack capture-i na 10 ms intervalu. Profilni buildovi koriste `-O3 -DNDEBUG -g -fno-inline` da bi se u stackovima vidjele funkcije algoritma; nemoj porediti njihovu brzinu s Release runovima iz prve tabele.

| Verzija | Glavni rad u stackovima | Thread stackovi u sync waitu |
|---|---|---:|
| v1.0.0 | `choose_cell`, pristup punoj DP matrici | 0 / 160 |
| v1.1.0 | `choose_cell`, pristup rolling redovima | 0 / 164 |
| v1.2.0 | `std::max`/`max_element`, `smith_waterman_trace` | 0 / 180 |
| v1.3.0 | `smith_waterman_trace`, pristup bazama i redovima | 0 / 180 |
| v1.4.0 | trace cell lambda i recurrence | 0 / 179 |
| v1.5.0 | wavefront tile funkcija; 1009 condition-variable i 19 mutex wait uzoraka | 1028 / 2221 (46,3%) |
| v1.6.0 | wavefront ready-tile funkcija; 675 condition-variable i 16 mutex wait uzoraka | 691 / 2205 (31,3%) |
| v1.7.0 | `smith_waterman_trace` cell loop; nema sync wait stackova u sažetku | 0 / 2023 |

Ovi brojevi su stack uzorci svih threadova tokom capture-a, ne precizno sabrane sekunde čekanja. `sample` pokazuje sinhronizacijsko čekanje na condition variable/mutex; ne mjeri koliko se čeka na RAM/cache pristup. Tabela zato ne izmišlja memorijsko čekanje: **cache misses, memorijski bandwidth i memory-stall ciklusi nisu izmjereni**. Broj page faults, nizak RSS ili IPC ne mogu zamijeniti te brojače. Za njih bi trebao Instruments/hardverski counter profil, koji nije instaliran na ovoj mašini.

## Izvorni podaci

| Verzija | Izvor i profilni izlaz |
|---|---|
| v1.0.0 | Početni commit `bc38b86`; [`sample`](reports/v1.0.0.sample.txt), [`time -l`](runtimes/v1.0.0.time.txt) |
| v1.1.0 | Rolling mod početnog commita `bc38b86`; [`sample`](reports/v1.1.0.sample.txt), [`time -l`](runtimes/v1.1.0.time.txt) |
| v1.2.0 | Tag `034fe05`; [`sample`](reports/v1.2.0.sample.txt), [`time -l`](runtimes/v1.2.0.time.txt) |
| v1.3.0 | Tag `f7a85fa`; [`sample`](reports/v1.3.0.sample.txt), [`time -l`](runtimes/v1.3.0.time.txt) |
| v1.4.0 | Tag `97310d3`; [`sample`](reports/v1.4.0.sample.txt), [`time -l`](runtimes/v1.4.0.time.txt) |
| v1.5.0 | Tag `6f56859`; [`sample`](reports/v1.5.0.sample.txt), [`time -l`](runtimes/v1.5.0.time.txt) |
| v1.6.0 | Tag `07e2e4e`; [`sample`](reports/v1.6.0.sample.txt), [`time -l`](runtimes/v1.6.0.time.txt) |
| v1.7.0 | Batch worker-pool varijanta u trenutnom radnom stablu, još bez Gita taga; [`sample`](reports/v1.7.0.sample.txt), [pet runova i medijani](runtimes/v1.7.0.time.txt), [pojedinačni podaci](runtimes/v1.7.0-runs.csv) |

FASTA ulaz je sačuvan u [`data/`](data/). Release runovi v1.0.0–v1.6.0 imaju kompletan `/usr/bin/time -l` izlaz u svojim logovima. v1.7.0 ima pet pojedinačnih mjerenja u CSV-u i izračunate medijane u svom sažetku.
