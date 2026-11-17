# Smith–Waterman: read mapping

Početni C++ projekat za temu **P21: Sequence Alignment and Read Mapping**. Program računa lokalno poravnanje svakog read-a prema jednoj referentnoj sekvenci. Za svaki read ispisuje najbolji score i interval poravnanja.

## Implementacije

- `auto` (podrazumijevani režim): bira `wavefront` za jedan read ili mali batch do četvrtine broja radnika; veće grupe idu kroz `batch`. Time koristi paralelizam unutar poravnanja kad bi batch pool imao premalo nezavisnih readova.
- `full`: direktna Smith–Waterman DP matrica. Jednostavna referentna verzija, sa `O(nm)` memorije po read-u.
- `rolling`: ista rekurencija, ali čuva samo prethodni i trenutni red DP matrice. Ovo je sekvencijalna memorijska optimizacija i koristi `O(m)` memorije po read-u.
- `trace`: rolling score redovi uz 2-bitni smjer za svaku DP ćeliju; početne koordinate se dobiju naknadnim traceback prolazom. Time se izbjegava propagiranje koordinata kroz svaku ćeliju, uz `O(nm)` bitova traceback memorije.
- `wavefront`: tiled wavefront Smith–Waterman za paralelno računanje jednog read-a na više CPU threadova. Pločica ulazi u bounded ready queue čim završe njene sjeverna i zapadna zavisnost; rezultat i traceback smjerovi su isti kao u `trace` modu.
- `batch`: persistent worker pool dijeli nezavisne readove između CPU threadova, a svaki worker računa svoj read `trace` kernelom. Ovo izbjegava pokretanje wavefront worker grupe za svaki kratki read.
- `simd`: ARM NEON inter-sequence SIMD. Jedna 128-bitna vektorska instrukcija obrađuje ćelije iz četiri nezavisna poravnanja istovremeno. Read-ovi se grupišu po četiri; read-ovi različitih dužina koriste neaktivne lane-ove nakon završetka kraćih sekvenci.

## Verzije algoritma

- **v1.0.0 — `full`**: sekvencijalna puna DP matrica i početna referentna implementacija.
- **v1.1.0 — `rolling`**: sekvencijalno računanje sa rolling redovima i manjom memorijskom potrošnjom.
- **v1.2.0 — `trace`**: single-read optimizacija sa kompaktno zapakovanim traceback smjerovima.
- **v1.3.0 — `trace`**: smjerovi se pakuju po redovima, bez početnog nuliranja traceback buffera.
- **v1.4.0 — `trace`**: raw pointer pristupi u kernelu i traceback pakovanje po blokovima od četiri ćelije.
- **v1.5.0 — `wavefront`**: paralelno računanje pločica unutar jednog read-a, uz worker pool i traceback kompatibilan s `trace`.
- **v1.6.0 — `wavefront`**: uklonjene globalne barijere između anti-dijagonala; pločice se pokreću čim se završe njihovi preduslovi.
- **v1.7.0 — `batch`**: nezavisni readovi se raspoređuju na jedan pool threadova, uz `trace` računanje po readu. Kandidat je u radnom stablu i još nema Git tag.

Ove verzije označavaju algoritamske varijante; `full`, `rolling` i `simd` su već prisutni u početnom Git commitu. `simd` obrađuje više read-ova zajedno i zato nije dio single-read serije.

Detaljan govorni tekst, objašnjenje bottleneckova i pitanja za odbranu su u [prezentacijskom vodiču](docs/prezentacija-rad-bs.md).

Apple CPU profili i rezultati single-read iteracija v1.1.0–v1.6.0 su u [profiling izvještaju](docs/profiling-v1.1.0.md). Na testiranom read-u od 10 kbp, uklanjanje wavefront barijera u v1.6.0 je oko 1.25× brže od v1.5.0. Izvještaj sadrži i sweep tile dimenzija te memorijski svjesne automatske postavke. `sample` profil ne daje hardverske brojače za memorijske stallove; za Instruments trace potreban je puni Xcode, dok je na ovoj mašini dostupan samo Command Line Tools.

Podrazumijevani `auto` režim koristi wavefront paralelizam za jedan read ili mali batch. Na Apple M5 Pro, 4 × 10 kbp read-a i 100 kbp reference, trajalo je 0,393 s naspram 1,027 s za `batch` (2,6× brže) uz peak RSS od 303,8 MiB naspram 959,0 MiB. Za 8 × 160 bp i 800 × 160 bp bira batch pool, koji bolje iskorištava nezavisne readove. Ponavljanja i ostale metrike su u [Apple sample sažetku](docs/apple-profiler-samples/README.md).

Tabela wall/DP vremena, GCUPS-a, CPU zauzeća, peak RSS-a, VM brojača i profiler wait uzoraka za v1.0.0–v1.7.0 je u [Apple sample sažetku](docs/apple-profiler-samples/README.md).

`full`, `rolling` i `trace` računaju jedan DP po read-u sekvencijalno. `wavefront` paralelizuje jedan DP preko threadova; `batch` paralelizuje više nezavisnih readova. `trace` koristi manje memorije od `full`; u odnosu na `rolling`, potrošnja zavisi od dužine read-a i reference. `simd` paralelizuje ćelije između četiri read-a na ARM64 procesorima koji podržavaju NEON. Sve implementacije računaju `O(nm)` DP ćelija po read-u. `simd` nije dostupan na x86-64 build-ovima.

Koristi se lokalno poravnanje s linearnom gap kaznom. Podrazumijevane vrijednosti su `match=2`, `mismatch=-1`, `gap=-2`. Koordinate su 0-based i half-open. Program vraća score i intervale, ali ne rekonstruiše traceback putanju niti CIGAR string. Ulazni format je FASTA; referentni fajl mora imati tačno jednu sekvencu.

## Izgradnja

Potrebni su CMake 3.16+ i C++17 kompajler.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Na macOS-u s Homebrewom, instaliraj CMake pomoću `brew install cmake`.

Na ARM64 macOS-u može se kompajlirati i direktno, bez CMake-a:

```sh
clang++ -O3 -DNDEBUG -std=c++17 -Wall -Wextra -Wpedantic -pthread src/main.cpp -o build/sw_mapper
```

## Pokretanje svih iteracija za prezentaciju

`scripts/demo_versions.py` kompajlira i pokreće v1.0–v1.7 i trenutni `auto` nad istim FASTA fajlovima i scoring parametrima. Historijske verzije se izvlače iz Git commita bez promjene trenutne grane. V1.7 batch i `final-auto` su lokalne netagovane iteracije i grade se iz trenutnog `src/main.cpp`.

Za uporedivi benchmark iz izvještaja prvo napravi fiksni ulaz:

```sh
python3 scripts/generate_benchmark.py --output-dir build/presentation-data \
  --reference-length 100000 --read-count 800 --read-length 160 \
  --substitutions 5 --seed 21021
```

Zatim kompajliraj sve verzije i pokreni ih jednu po jednu:

```sh
python3 scripts/demo_versions.py \
  --reference build/presentation-data/reference.fa \
  --reads build/presentation-data/reads.fa \
  --output-dir build/presentation-versions --repeats 1
```

Skripta na kraju ispiše prvo poravnanje svake verzije i tabelu vremena, DP vremena, GCUPS-a, CPU vremena/zauzeća, peak RSS-a, IPC-a i provjere jednakosti izlaza. Puni CSV izlaz svakog read-a i `/usr/bin/time -l` logovi ostaju u `build/presentation-versions/results/`; sažetak je `summary.csv`. Za brzu probu možeš koristiti `examples/reference.fa` i `examples/reads.fa`, ali su toliko mali da vrijeme nije korisno za poređenje.

Za jednu određenu verziju, dodaj `--versions`:

```sh
python3 scripts/demo_versions.py --reference build/presentation-data/reference.fa \
  --reads build/presentation-data/reads.fa --versions v1.4.0
```

Ili koristi Bash helper koji samo gradi i pokreće izabranu verziju:

```sh
scripts/run_version.sh v1.4.0 \
  build/presentation-data/reference.fa build/presentation-data/reads.fa
```

Zamijeni `v1.4.0` s oznakom željene verzije. Za posljednju lokalnu varijantu koristi `final-auto`; za batch iteraciju `v1.7.0-batch`. Bash helper svaki put kompajlira taj izvor i ispisuje CSV rezultate i `/usr/bin/time -l` log. Python runner je za izgradnju i poređenje više verzija odjednom: lakše zapisuje per-version CSV fajlove, parsira CPU/RSS podatke i pravi zbirni CSV. Obje skripte kompajliraju i pokreću C++ mapper; Python nije implementacija algoritma.

Samo build bez pokretanja:

```sh
python3 scripts/demo_versions.py --build-only --output-dir build/presentation-versions
```

Kasnije pokreni već izgrađene binarne fajlove sa drugim ulazima uz `--run-only`, `--reference` i `--reads`. Opcije `--match`, `--mismatch`, `--gap`, `--threads` i `--repeats` mogu se podesiti; podrazumijevani scoring je `2, -1, -2`, a `--threads 0` znači automatski broj threadova.

## Pokretanje

```sh
./build/sw_mapper --reference examples/reference.fa --reads examples/reads.fa --mode auto
./build/sw_mapper --reference examples/reference.fa --reads examples/reads.fa --mode rolling
./build/sw_mapper --reference examples/reference.fa --reads examples/reads.fa --mode full
./build/sw_mapper --reference examples/reference.fa --reads examples/reads.fa --mode trace
./build/sw_mapper --reference examples/reference.fa --reads examples/reads.fa --mode wavefront --threads 0
./build/sw_mapper --reference examples/reference.fa --reads examples/reads.fa --mode batch --threads 0
./build/sw_mapper --reference examples/reference.fa --reads examples/reads.fa --mode simd
```

Za `wavefront` i `batch`, `--threads 0` bira broj radnika automatski. `auto` koristi isti broj: `wavefront` bira za jedan read ili do jedne četvrtine radnika u broju readova, a `batch` za veće grupe. Wavefront tile postavke važe samo za `wavefront`: `--tile-rows 0` bira 8 redova za readove do 160 bp, 16 do 512 bp i 64 za duže readove. `--tile-columns 0` automatski bira 4096 kolona za readove duže od 512 bp uz reference od najmanje 200 kbp, a 2048 u ostalim slučajevima. Ručno zadana širina mora biti višekratnik broja 4.

Opcije bodovanja:

```sh
./build/sw_mapper --reference examples/reference.fa --reads examples/reads.fa \
  --mode rolling --match 2 --mismatch -1 --gap -2
```

Rezultati su CSV na stdout. Vrijeme računanja i procijenjeni GCUPS ispisuju se na stderr; mjerenje ne uključuje parsiranje FASTA fajlova ni CSV ispis. Koristi Release build i isti ulaz/scoring za poređenje režima.

## Apple performance sampling

Na macOS-u `scripts/profile_apple_counters.py` snima `powermetrics` CPU/frekvenciju i per-process IPC uz `/usr/bin/time -l` instrukcije, cikluse, RSS i page brojače. Pokreni ga iz terminala preko `sudo`; mjeri se ciljna komanda kao trenutni korisnik:

```sh
sudo python3 scripts/profile_apple_counters.py \
  --output-dir build/performance-counters/auto \
  -- ./build/sw_mapper --reference build/single-profile-data/reference.fa \
  --reads build/single-profile-data/reads.fa
```

Izlazi se čuvaju u navedenom direktoriju. Uzorkovanje je sistemsko, pa je najbolje mjeriti workload-e koji traju barem nekoliko intervala od 50 ms. `powermetrics` ne izlaže cache misses ili RAM-stall cikluse. Apple Instruments `CPU Counters` instrument zahtijeva punu Xcode instalaciju; aktivan Command Line Tools direktorij trenutno nema `xctrace` podršku.

## Poređenje i mjerenje

`batch` znači da svaki radnik uzme drugi read i sam ga računa sekvencijalnim `trace` kernelom. `wavefront` dijeli pločice jednog read-a između radnika. `auto` koristi wavefront kad batch ima premalo readova da zaposli procesor; za veće grupe obrađuje readove kroz zajednički batch pool. Na testiranom 18-thread M5 Pro, četiri read-a su bila brža kroz wavefront, a osam kratkih readova kroz batch. To nije isto što i pokretanje posebnog wavefront poola za svaki read.

`auto`, `full`, `rolling`, `trace`, `wavefront`, `simd` i `batch` trebaju dati identične score-ove i koordinate. `full`, `rolling` i `trace` računaju jedan read; `wavefront` koristi više threadova za taj read; `simd` obrađuje read-ove u grupi; `batch` raspoređuje readove na threadove. `full` koristi mnogo memorije, pa ga pokreći samo nad malim ulazima. Za batch poređenje pokreni:

```sh
python3 scripts/generate_benchmark.py --output-dir build/benchmark-data
python3 scripts/benchmark.py --binary build/sw_mapper \
  --reference build/benchmark-data/reference.fa \
  --reads build/benchmark-data/reads.fa --repeats 7
```

`scripts/benchmark.py` provjeri da su CSV izlazi rolling, trace, SIMD i batch verzija identični, zagrije režime i rotira redoslijed mjerenja. Za fokus na jedan read koristi `scripts/benchmark_single.py`, koji upoređuje full, rolling i trace; dodaj `--include-wavefront` da uključi paralelni režim unutar jednog poravnanja. Vrijeme u programu obuhvata DP računanje i alokacije, ali ne uključuje učitavanje FASTA ni ispis. Zabilježi CPU, kompajler, build zastavice, veličine sekvenci i broj ponavljanja.

```sh
python3 scripts/generate_benchmark.py --output-dir build/single-read-data \
  --reference-length 100000 --read-count 1 --read-length 160 \
  --substitutions 5 --seed 21021
python3 scripts/benchmark_single.py --binary build/sw_mapper \
  --reference build/single-read-data/reference.fa \
  --reads build/single-read-data/reads.fa --repeats 7
```

Poređenje `wavefront` moda sa sekvencijalnim `trace` izlazom može se pokrenuti ovako:

```sh
python3 scripts/benchmark_single.py --binary build/sw_mapper \
  --baseline-binary build/releases/sw_mapper-v1.4.0 --trace-only \
  --include-wavefront --reference build/single-read-data/reference.fa \
  --reads build/single-read-data/reads.fa --repeats 7
```

Za poređenje `wavefront` verzija kroz benchmark skriptu, dodaj `--baseline-mode wavefront` uz prethodnu verziju:

```sh
python3 scripts/benchmark_single.py --binary build/profile-release/sw_mapper \
  --baseline-binary build/releases/sw_mapper-v1.5.0 \
  --baseline-mode wavefront --trace-only --include-wavefront \
  --reference build/single-profile-data/reference.fa \
  --reads build/single-profile-data/reads.fa --repeats 7
```

Na ovoj mašini (Apple M5 Pro, ARM64; Apple clang 21.0.0; CMake Release, `-O3 -DNDEBUG`) lokalni rolling/SIMD mjerni ciklus dao je:

| Režim | Medijana vremena | GCUPS |
|---|---:|---:|---|
| `rolling` | 0.358093 s | 0.343 |
| `simd` | 0.157451 s | 0.780 |

SIMD je na ovom konkretnom batchu bio **2.27× brži** od rolling verzije. Batch je sadržavao jednu referencu dužine 12,000 bp i 64 read-a dužine 160 bp, svaki s pet zamjena baza; sjeme generatora je `21021`. To je rezultat jedne mašine i jednog ulaza, ne garancija ubrzanja za druge podatke.

## Poređenje s javnom sekvencijalnom implementacijom

Za nezavisno poređenje koristili smo sekvencijalnu klasu iz [giorgapost/smith-waterman-algorithm](https://github.com/giorgapost/smith-waterman-algorithm), commit `98660d57f818a1d18a6129c6e66947d598acf878` (GPL-3.0). Repo se kompajlira lokalno samo iz njegovih `Framework.cpp` i `SequentialImplementation.cpp` fajlova uz mali benchmark driver; paralelne implementacije se ne buildaju. Kod vanjskog projekta nije kopiran u ovaj repo.

Na istom skupu i scoring parametrima, medijane DP-fill vremena iz sedam mjerenja bile su:

| Implementacija | DP vrijeme | GCUPS | Poređenje s vanjskom sekvencijalnom |
|---|---:|---:|---:|
| Vanjska sekvencijalna, puna matrica | 0.392593 s | 0.313 | 1.00× |
| Naš `full` | 0.315261 s | 0.390 | 1.25× |
| Naš `rolling` | 0.364843 s | 0.337 | 1.08× |
| Naš `simd` (ARM NEON) | 0.155273 s | 0.791 | 2.53× |

Vanjski kod je koristio matricu `vector<vector<int>>`, bilježio sve ćelije s maksimalnim score-om i rekonstruisao traceback. Njegov objavljeni `calcCellsTime` mjeri punjenje DP matrice bez matrice alokacije; naši `dp_seconds` također zaustavljaju sat oko DP petlji, nakon što su scratch bufferi alocirani. Naše ćelije dodatno propagiraju početne koordinate. Zato ovo poređenje je bliže kernels-only poređenju, ali se implementacije ipak razlikuju u količini rada po ćeliji.

Provjera rezultata: svi 64 score-a su identični za `full`, `rolling` i `simd`, a svaki izabrani interval bio je među optimalnim intervalima vanjske implementacije za 64/64 read-a. Vanjski `totalTime` uključuje traceback i iznosio je 0.400892 s medijanom; ne poredi se direktno s našim vremenom jer naš alat ne konstruiše CIGAR/poravnate sekvence.

Ponovi online poređenje ovako nakon kloniranja izvora na zabilježeni commit:

```sh
git clone https://github.com/giorgapost/smith-waterman-algorithm build/online-source
git -C build/online-source checkout 98660d57f818a1d18a6129c6e66947d598acf878
python3 scripts/benchmark_online.py --binary build/sw_mapper \
  --online-source build/online-source \
  --reference build/benchmark-data/reference.fa \
  --reads build/benchmark-data/reads.fa --repeats 7
```

### Veći benchmark suite

Za tri dataset profila i mnogo ponavljanja pokreni:

```sh
git clone https://github.com/giorgapost/smith-waterman-algorithm build/online-source
git -C build/online-source checkout 98660d57f818a1d18a6129c6e66947d598acf878
python3 scripts/benchmark_suite.py --online-source build/online-source --repeats 7
```

Suite automatski izgradi CMake Release verziju, napravi determinističke podatke, provjeri rezultate svih četiriju implementacija i mjeri ih redom izmiješanim ponavljanjima. Profili rastu od 4 kbp reference i 64 read-a do 32 kbp i 192 read-a; srednji i veliki profil imaju promjenjive dužine read-ova radi mjerenja SIMD lane iskorištenja. Izvršena su 84 vremenska mjerenja, ukupno oko 44,45 milijardi DP ćelija. [Detaljni izvještaj](docs/benchmark-results.md), [sirovi podaci za svaki run](docs/benchmark-raw-results.csv), [sažeti CSV](docs/benchmark-summary.csv) i [peak RSS mjerenja](docs/benchmark-memory.csv) su verzionisani u repozitoriju.

U ovim mjerenjima NEON SIMD je bio najbrži u sva tri profila. Njegovo ubrzanje u odnosu na vanjski sekvencijalni kod palo je s 2.64× na malom profilu na 1.87× na velikom. `rolling` troši mnogo manje memorije od našeg `full` pristupa, ali nije bio brži u DP petlji na ovoj mašini. Veliki profil je imao peak RSS oko 1.28 GiB za `full`, naspram 3.44 MiB za `rolling` i 4.97 MiB za SIMD.

| Režim | Vrijeme | Memorija po read-u |
|---|---:|---:|
| `full` | `O(nm)` | `O(nm)` |
| `rolling` | `O(nm)` | `O(m)` |
| `simd` | `O(nm)` | `O(m)` po SIMD batchu |

Ovdje je `n` dužina read-a, a `m` dužina reference.

## Teorijski GPU dio

GPU implementacija nije dio projekta. Teorijski wavefront/tiling pristup i detaljniji CPU SIMD kontekst nalaze se u odvojenom `docs.md` dokumentu izvan git repozitorija.

## Izvor algoritma

Smith, T. F. & Waterman, M. S. (1981), “Identification of Common Molecular Subsequences”, *Journal of Molecular Biology*, 147(1), 195–197. [DOI: 10.1016/0022-2836(81)90087-5](https://doi.org/10.1016/0022-2836(81)90087-5).

CPU inter-sequence SIMD pristup je inspirisan Rognesom, “Faster Smith-Waterman database searches with inter-sequence SIMD parallelisation”, *BMC Bioinformatics* 12, 221 (2011). [DOI: 10.1186/1471-2105-12-221](https://doi.org/10.1186/1471-2105-12-221).
