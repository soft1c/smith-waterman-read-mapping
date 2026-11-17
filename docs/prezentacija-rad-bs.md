# Tekst za prezentaciju: optimizacija Smith–Waterman algoritma

## Govorni tekst

U ovom radu sam implementirao lokalno poravnanje sekvenci Smith–Waterman algoritmom i zatim ga optimizirao kroz više iteracija. Cilj nije bio samo da program daje ispravan score, nego da razumijem gdje troši vrijeme i memoriju, pa da svaka promjena odgovara izmjerenom problemu.

Smith–Waterman popunjava matricu dinamičkog programiranja. Svaka ćelija zavisi od dijagonalne ćelije, ćelije iznad i ćelije lijevo. Score se nikad ne spušta ispod nule, jer se traži lokalno, a ne globalno poravnanje. Kada se matrica izračuna, najveći score daje kraj najboljeg poravnanja, a traceback smjerovi vode do početka. Program vraća score i početne i završne koordinate za read i referencu; ne pravi CIGAR string.

Prvo sam napravio `full` verziju. Ona čuva cijelu matricu, a ćelija pored score-a nosi i koordinate početka poravnanja. To je jednostavna referentna implementacija i pogodna je za poređenje, ali memorija raste s proizvodom dužine read-a i reference. Na referenci od 100.000 baza i 800 readova od 160 baza, izmjereni peak RSS bio je oko 370,8 MiB.

U `rolling` verziji sam zadržao samo prethodni i trenutni red matrice. Memorija je pala na oko 6,9 MiB. Međutim, manja memorija sama po sebi nije značila i kraće vrijeme: rolling verzija je na velikom ulazu trajala 39,69 sekundi, naspram 35,11 za `full`. Razlog je bio u radu po ćeliji: rolling je i dalje prenosio koordinate početka kroz svaku ćeliju, a smanjen pristup memoriji nije nadoknadio taj dodatni posao.

Tu sam promijenio pristup u `trace` verziji. U glavnoj petlji čuvam score samo za dva reda, a za svaku ćeliju zapakujem smjer tracebacka u dva bita. Poslije prolaza se početne koordinate pronađu jednim tracebackom. Tako se izbjegava prenošenje dvije koordinate kroz sve ćelije. Na jednom read-u od 160 baza i referenci od 100.000 baza, trace je bio 1,8 puta brži od rollinga. Na cijelom 800-read benchmarku v1.2 je trajao 20,40 sekundi.

Da bih našao sljedeća uska grla, koristio sam Apple `sample` za uzorkovanje call stackova i `/usr/bin/time -l` za ukupno vrijeme, CPU vrijeme, RSS, page faults, instrukcije i cikluse. U početnom profilu zabilježeno je 6.846 stack uzoraka: 79,5% bilo je u izboru ćelije i prenošenju koordinata, 10,4% u ažuriranju najboljeg poravnanja, a 6,4% u računanju score-a zamjene baza. To je pokazalo da prva meta treba biti rad po ćeliji, ne parsiranje FASTA fajlova ili ispis. U trace kernelu sampler je pokazao dodatni rad pri inicijalizaciji smjerova i pozive `std::max` preko initializer liste. Za čitanje izvornog call stacka koristio sam i profilni build sa `-fno-inline`; taj build je služio za pronalaženje hot spotova, ne za poređenje brzine.

U v1.3 sam smjerove pakovao po četiri ćelije u jedan bajt, upisivao bajt jednom i izostavio nuliranje buffera koji se u cijelosti prepisuje. Maksimum sam računao direktnim poređenjima. Na pojedinačnom read-u poboljšanje je bilo malo, oko četiri do sedam procenata. To je važan rezultat: optimizacija koja izgleda dobro u profileru nije automatski brža na svakom workloadu.

U v1.4 sam iz unutrašnje petlje podigao pristup bazama i score redovima, tako da se pointeri i baza read-a ne računaju ponovo za svaku ćeliju. Smjerovi su se i dalje pakovali u blokovima od četiri. U pojedinačnim mjerenjima v1.4 je bio oko 1,28 puta brži od v1.3 za kratki read i oko 1,26 puta za read od 10.000 baza.

Kada je jedna sekvenca duga, ostaje veliki serijski DP račun. Zavisnosti gore i lijevo znače da ne mogu proizvoljno podijeliti ćelije između threadova. U v1.5 sam zato podijelio matricu na pločice. Pločica može početi kada završe njena sjeverna i zapadna susjedna pločica; nezavisne pločice se izvršavaju paralelno. Na read-u od 10.000 baza, wavefront je bio oko 8,3 puta brži od v1.4 trace verzije. Cijena je bila dodatna memorija za granice pločica i sinhronizacija između radnika.

Profil v1.5 je pokazao da coordinator dio vremena čeka da se završi cijela anti-dijagonala. U v1.6 sam zamijenio globalne barijere brojačima zavisnosti po pločici i spremnim redom. Pločica se stavlja u red čim su završeni njeni preduslovi, bez čekanja na nepovezane pločice. Vrijeme čekanja na condition variable u uzorcima glavnog threada palo je sa oko 19 procenata na oko 4,6 procenata. Na read-u od 10.000 baza v1.6 je bio 1,25 puta brži od v1.5.

Zatim sam pogledao oblik glavnog benchmarka: 800 kratkih readova prema istoj referenci. Za taj ulaz je prirodnije paralelizovati različite readove nego svako kratko poravnanje posebno. V1.7 uvodi jedan zajednički batch worker pool. Svaki radnik uzima sljedeći read i računa ga trace kernelom. Pet mjerenja dalo je medijanu od 0,944 sekunde i 13,55 milijardi ćelija u sekundi, odnosno GCUPS. To je oko 2,6 puta brže od ranijeg v1.6 rezultata na tom workloadu. Peak RSS je porastao na oko 99 MiB jer više radnika istovremeno drži traceback buffere.

Trenutni `auto` režim bira paralelizam prema broju readova. Ako ima jedan read ili mali batch, do četvrtine dostupnog broja radnika, koristi wavefront unutar jednog poravnanja. Na M5 Pro sa 18 CPU radnika to znači najviše četiri read-a. Za veće grupe koristi batch pool. Na četiri read-a od po 10.000 baza wavefront je trajao 0,393 sekunde, a batch trace 1,027 sekundi, uz peak RSS od oko 304 MiB naspram 959 MiB. Za osam kratkih readova batch je ponovo bio malo brži, pa ga auto bira. Za 800 kratkih readova auto je odabrao batch i izmjerio 0,912 sekundi.

Ovo objašnjava i zašto ne pokrećem wavefront pool za svaki read istovremeno. Jedan wavefront već koristi oko 17 jezgri. Četiri nezavisna wavefront poola sa po četiri radnika bila su sporija od redne obrade wavefrontom i zajedno imala oko 1,2 GiB zbirnog maksimalnog RSS-a. Batch i wavefront rješavaju dva različita oblika paralelizma: batch raspoređuje readove, a wavefront pločice unutar read-a. Auto bira onaj koji je na izmjerenom ulazu bolje iskoristio procesor bez nepotrebnog memorijskog rasta.

Glavni uporedivi benchmark koristio je istu determinističku referencu od 100.000 baza i 800 readova od 160 baza, pet zamjena po readu i seed 21021. To je ukupno 12,8 milijardi DP ćelija; bodovanje je match 2, mismatch −1 i gap −2. Build je bio Release sa Apple clangom i `-O3 -DNDEBUG -g`, na Apple M5 Pro sa 18 CPU jezgri i 24 GB RAM-a. Sve verzije su dale iste score-ove i koordinate.

Rezultat za cijeli 800-read workload išao je od 35,11 sekundi za `full` do 0,94 sekunde za batch v1.7. Trenutni auto kandidat je na istom batch ulazu izmjerio 0,912 sekundi. Memorija se u tom poređenju smanjila sa oko 371 MiB za `full` na oko 99 MiB za batch. Za jedan dugi read, automatski wavefront je izmjeren 9,7 puta brži od trace verzije, uz veći peak RSS: približno 304 MiB naspram 241 MiB.

Važno ograničenje mjerenja je da `sample` prikazuje gdje su threadovi izvršavali kod i koliko su uzoraka proveli u sinhronizaciji, ali ne mjeri direktno cache misses ni čekanje na RAM. `/usr/bin/time -l` daje instrukcije, cikluse i IPC, ali ni IPC nije memory-stall brojač. Zato ne tvrdim da smo izmjerili memorijsko čekanje. Za takve hardverske brojače potreban je Apple Instruments CPU Counters iz pune Xcode instalacije. Umjesto toga sam zaključke vezao za ono što smo stvarno izmjerili: runtime, CPU zauzeće, memoriju i stack uzorke.

Zaključak je da ne postoji jedan najbolji paralelni režim za svaki ulaz. Za mnogo kratkih readova najbolja jedinica posla je read i koristi se batch; za jedan ili nekoliko dugih readova bolji je wavefront. Najvažniji rezultat iteracija je kombinacija kompaktnog tracebacka, uklonjenih wavefront barijera i automatskog izbora paralelizma prema obliku ulaza.

## Sažetak iteracija za brzo ponavljanje

| Iteracija | Glavna izmjena | Šta je pokazalo mjerenje |
|---|---|---|
| v1.0.0 `full` | Cijela DP matrica i koordinate u ćelijama | 35,11 s; 370,8 MiB na 800 × 160 bp |
| v1.1.0 `rolling` | Samo dva score reda | Memorija pada na 6,9 MiB; vrijeme 39,69 s na velikom workloadu |
| v1.2.0 `trace` | 2-bitni traceback smjerovi; koordinate se vraćaju naknadno | 20,40 s; 0,628 GCUPS |
| v1.3.0 `trace` | Pakovanje četiri smjera po bajtu; izostavljeno inicijalizovanje | Malo ubrzanje na single-read slučajevima |
| v1.4.0 `trace` | Pointeri i baze izvan unutrašnje petlje | 15,97 s; 0,802 GCUPS |
| v1.5.0 `wavefront` | Paralelne pločice po anti-dijagonalama | 5,13 s na 800-read ulazu; 8,29× za jedan 10 kbp read naspram trace |
| v1.6.0 `wavefront` | Per-tile zavisnosti i spreman red | 2,46 s; manje čekanja na barijere |
| v1.7.0 `batch` | Jedan worker pool preko nezavisnih readova | 0,944 s medijana; 13,55 GCUPS |
| Lokalni `auto` kandidat | Wavefront za male, batch za veće grupe | 0,912 s za 800 kratkih readova; 0,104 s za jedan 10 kbp read |

Brojevi v1.0.0–v1.6.0 u velikom benchmarku su pojedinačni runovi, v1.7.0 je medijana pet runova, a lokalni auto rezultati su zasebna serija. Koristi ih kao rezultat na toj mašini i tom ulazu, ne kao univerzalni speedup.

## Pitanja koja mogu doći na odbrani

**Zašto je rolling sporiji iako koristi manje memorije?**  
Zato što smanjenje memorije ne uklanja rad po ćeliji. Rolling verzija i dalje propagira koordinate početka kroz sve ćelije; trace to radi tek u kratkom traceback prolazu.

**Šta znači GCUPS?**  
Milijarde DP ćelija izračunate u sekundi. Računa se kao broj readova puta dužina read-a puta dužina reference, podijeljeno vremenom i milijardom.

**Zašto wavefront ne radi uvijek?**  
Za mnogo kratkih readova pokretanje i koordinacija wavefronta svakog read-a troše više od korisnog paralelnog rada. Batch dijeli nezavisne readove i u tom slučaju efikasnije puni jezgre.

**Da li su izmjereni cache misses ili čekanje RAM-a?**  
Ne. Imamo RSS, page brojače, instrukcije, cikluse, IPC i stack uzorke. Nijedan od tih brojeva sam po sebi nije memorijski stall counter.

**Da li su v1.7 i auto zvanični Git tagovi?**  
Ne. Git tagovi postoje do v1.6. V1.7 batch i trenutni auto dispatcher su lokalne, netagovane iteracije. Demo skripta to posebno označava i gradi ih iz radnog stabla.
