# Smith–Waterman: read mapping

Početni C++ projekat za temu **P21: Sequence Alignment and Read Mapping**. Program računa lokalno poravnanje svakog read-a prema jednoj referentnoj sekvenci. Za svaki read ispisuje najbolji score i interval poravnanja.

## Implementacije

- `full`: direktna Smith–Waterman DP matrica. Jednostavna referentna verzija, sa `O(nm)` memorije po read-u.
- `rolling`: ista rekurencija, ali čuva samo prethodni i trenutni red DP matrice. Ovo je sekvencijalna memorijska optimizacija i koristi `O(m)` memorije po read-u.
- `trace`: rolling score redovi uz 2-bitni smjer za svaku DP ćeliju; početne koordinate se dobiju naknadnim traceback prolazom. Time se izbjegava propagiranje koordinata kroz svaku ćeliju, uz `O(nm)` bitova traceback memorije.
- `simd`: ARM NEON inter-sequence SIMD. Jedna 128-bitna vektorska instrukcija obrađuje ćelije iz četiri nezavisna poravnanja istovremeno. Read-ovi se grupišu po četiri; read-ovi različitih dužina koriste neaktivne lane-ove nakon završetka kraćih sekvenci.

## Verzije algoritma

- **v1.0.0 — `full`**: sekvencijalna puna DP matrica i početna referentna implementacija.
- **v1.1.0 — `rolling`**: sekvencijalno računanje sa rolling redovima i manjom memorijskom potrošnjom.
- **v1.2.0 — `trace`**: single-read optimizacija sa kompaktno zapakovanim traceback smjerovima.
- **v1.3.0 — `trace`**: smjerovi se pakuju po redovima, bez početnog nuliranja traceback buffera.

Ove verzije označavaju algoritamske varijante; `full`, `rolling` i `simd` su već prisutni u početnom Git commitu. `simd` obrađuje više read-ova zajedno i zato nije dio single-read serije.

Apple CPU profili i rezultati single-read iteracija v1.1.0–v1.3.0 su u [profiling izvještaju](docs/profiling-v1.1.0.md). `trace` v1.3.0 je 1.04× brži od v1.2.0 na read-u od 160 bp i 1.07× na read-u od 10 kbp. `sample` profil ne daje hardverske brojače za memorijske stallove; za Instruments trace potreban je puni Xcode, dok je na ovoj mašini dostupan samo Command Line Tools.

`full`, `rolling` i `trace` računaju jedan DP po read-u sekvencijalno. `trace` koristi manje memorije od `full`; u odnosu na `rolling`, potrošnja zavisi od dužine read-a i reference. `simd` paralelizuje ćelije između četiri read-a na ARM64 procesorima koji podržavaju NEON. Sve implementacije računaju `O(nm)` DP ćelija po read-u. `simd` nije dostupan na x86-64 build-ovima.

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
clang++ -O3 -DNDEBUG -std=c++17 -Wall -Wextra -Wpedantic src/main.cpp -o build/sw_mapper
```

## Pokretanje

```sh
./build/sw_mapper --reference examples/reference.fa --reads examples/reads.fa --mode rolling
./build/sw_mapper --reference examples/reference.fa --reads examples/reads.fa --mode full
./build/sw_mapper --reference examples/reference.fa --reads examples/reads.fa --mode trace
./build/sw_mapper --reference examples/reference.fa --reads examples/reads.fa --mode simd
```

Opcije bodovanja:

```sh
./build/sw_mapper --reference examples/reference.fa --reads examples/reads.fa \
  --mode rolling --match 2 --mismatch -1 --gap -2
```

Rezultati su CSV na stdout. Vrijeme računanja i procijenjeni GCUPS ispisuju se na stderr; mjerenje ne uključuje parsiranje FASTA fajlova ni CSV ispis. Koristi Release build i isti ulaz/scoring za poređenje režima.

## Poređenje i mjerenje

`full`, `rolling`, `trace` i `simd` trebaju dati identične score-ove i koordinate. `full`, `rolling` i `trace` računaju jedan read; `simd` obrađuje read-ove u grupi. `full` koristi mnogo memorije, pa ga pokreći samo nad malim ulazima. Za batch poređenje pokreni:

```sh
python3 scripts/generate_benchmark.py --output-dir build/benchmark-data
python3 scripts/benchmark.py --binary build/sw_mapper \
  --reference build/benchmark-data/reference.fa \
  --reads build/benchmark-data/reads.fa --repeats 7
```

Batch benchmark provjeri da su CSV izlazi rolling, trace i SIMD verzija identični. Za fokus na jedan read koristi `scripts/benchmark_single.py`, koji upoređuje full, rolling i trace, provjeri jednakost izlaza, zagrije režime i rotira redoslijed mjerenja. Vrijeme u programu obuhvata DP računanje i alokacije, ali ne uključuje učitavanje FASTA ni ispis. Zabilježi CPU, kompajler, build zastavice, veličine sekvenci i broj ponavljanja.

```sh
python3 scripts/generate_benchmark.py --output-dir build/single-read-data \
  --reference-length 100000 --read-count 1 --read-length 160 \
  --substitutions 5 --seed 21021
python3 scripts/benchmark_single.py --binary build/sw_mapper \
  --reference build/single-read-data/reference.fa \
  --reads build/single-read-data/reads.fa --repeats 7
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
