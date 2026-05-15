# BLAST Visualiser

Program w Pythonie, który uruchamia lokalne porównanie **blastn** dwóch
sekwencji DNA (query i subject) i tworzy graficzną wizualizację uzyskanych
dopasowań, naniesionych na obie sekwencje wraz z ich opisanymi elementami
(geny, CDS-y, eksony, introny, tRNA, rRNA, misc_feature).

Program został zaprojektowany do analiz porównawczych blisko spokrewnionych
genomów, ze szczególnym uwzględnieniem:

- wykrywania horyzontalnego transferu genów (HGT) między pasożytniczymi
  roślinami a ich gospodarzami (np. *Cuscuta*, *Orobanche*, *Phelipanche*),
- porównywania genomów organellarnych (plastydowych, mitochondrialnych),
- wskazywania regionów pochodzenia plastydowego w genomach mitochondrialnych
  (intracellular gene transfer),
- wizualizacji zachowanych i rozbieżnych regionów między blisko spokrewnionymi
  gatunkami.

---

## Możliwości programu

- Lokalne uruchamianie **blastn** — bez konieczności użycia BLAST online.
- Wizualizacja dopasowań na obu sekwencjach wraz z liniami łączącymi.
- Elementy sekwencji (features) rysowane jako kolorowe prostokąty nad/pod linią
  każdej sekwencji.
- Kolor każdego dopasowania koduje procent identyczności (`pident`)
  w gradiencie „rainbow" z osadzoną legendą skali.
- Specjalne wyróżnianie zdefiniowanych przez użytkownika klas features
  (domyślnie misc_features oznaczone jako `plastid-derived` — rysowane dalej
  od linii sekwencji i w dedykowanym kolorze).
- Trzy formaty wyjściowe: **SVG** (skalowalny), **PDF** (do druku), **PNG**
  (rastrowy).
- Pełna konfiguracja przez plik JSON (geometria, czcionki, kolory, filtry
  features).
- Sekwencje można podawać na trzy sposoby — różne dla query i subject.

---

## Wymagania

- Python 3.8 lub nowszy
- **NCBI BLAST+** (`makeblastdb`, `blastn`) zainstalowane i dostępne w `PATH`
- Pakiety Pythona: `biopython`, `pandas`, `pycairo`, `matplotlib`

### Instalacja na Fedorze

```bash
sudo dnf install ncbi-blast+ python3-cairo
pip install biopython pandas matplotlib
```

### Instalacja na Debianie / Ubuntu

```bash
sudo apt install ncbi-blast+ libcairo2-dev pkg-config python3-dev
pip install biopython pandas pycairo matplotlib
```

### Instalacja przez conda / mamba (dowolny system)

```bash
mamba install -c bioconda blast
mamba install -c conda-forge biopython pandas pycairo matplotlib
```

---

## Pliki w projekcie

| Plik | Rola |
|---|---|
| `blast_and_show.py` | Główny skrypt — interfejs linii poleceń. |
| `blasting.py` | Tworzenie bazy BLAST i uruchamianie `blastn`. |
| `dataimport.py` | Pobieranie sekwencji z GenBank, parsowanie plików `.gb`, eksport features. |
| `draw.py` | Generowanie grafiki (SVG / PDF / PNG) przy użyciu Cairo. |
| `settings.json` | Domyślna konfiguracja rysunku. |
| `settings_short.json` | Konfiguracja dopasowana do krótkich sekwencji (1–15 kb). |
| `settings_long.json` | Konfiguracja dopasowana do długich sekwencji (200–800 kb). |

---

## Sposób użycia

```
python blast_and_show.py [argumenty]
```

### Sposoby podawania sekwencji

Dla każdej z dwóch sekwencji (query i subject) trzeba podać **dokładnie jeden**
z następujących wariantów:

1. **Numer dostępu (accession) GenBank** (`-q` / `-s`) — sekwencja zostanie
   pobrana automatycznie; NCBI wymaga w tym przypadku adresu e-mail (`-e`).
2. **Lokalny plik GenBank** (`.gb`) (`-a` / `-c`) — daje pełną adnotację
   features na rysunku.
3. **Lokalny plik FASTA** (`-b` / `-d`) — features nie będą rysowane dla tej
   sekwencji.

### Argumenty linii poleceń

| Skrót | Pełna nazwa | Opis | Domyślnie |
|:---:|---|---|---|
| `-q` | `--query` | Numer dostępu query w GenBank. | |
| `-Q` | `--query_gb_file` | Lokalny plik GenBank dla query (z adnotacjami). | |
| `-a` | `--query_fasta_file` | Lokalny plik FASTA dla query (bez adnotacji). | |
| `-s` | `--subject` | Numer dostępu subject w GenBank. | |
| `-S` | `--subject_gb_file` | Lokalny plik GenBank dla subject (z adnotacjami). | |
| `-b` | `--subject_fasta_file` | Lokalny plik FASTA dla subject (bez adnotacji). | |
| `-e` | `--email` | Prawdziwy adres e-mail. Wymagany przez NCBI przy `-q`/`-s`. | |
| `-o` | `--output_dir` | Katalog wyjściowy (zostanie utworzony, jeśli nie istnieje). | `output` |
| `-c` | `--config` | Plik JSON z konfiguracją rysunku. | `settings.json` |
| `-t` | `--title` | Tytuł wyświetlany na rysunku. | `"Results of blastn"` |
| `-m` | `--minlen` | Minimalna długość dopasowania (bp) do narysowania. | `100` |
| `-v` | `--e_value` | Próg E-value dla BLAST. | `1e-10` |

> **Uwaga o konwencji skrótów:** symetria między query a subject jest
> wbudowana w skróty — małe litery dla query (`-q`, `-Q`, `-a`),
> odpowiadający im przełącznik w tej samej kolumnie dla subject
> (`-s`, `-S`, `-b`). Numery dostępu i pliki GenBank używają tej samej
> litery (`-q`/`-Q`, `-s`/`-S`), pliki FASTA — sąsiednich liter (`-a`/`-b`).

Pełna pomoc dostępna po uruchomieniu `python blast_and_show.py --help`.

---

## Przykłady

```bash
# Obie sekwencje jako lokalne pliki GenBank (gospodarz vs. plastyd pasożyta)
python blast_and_show.py \
    -Q Laburnum_anagyroides-OZ176120.gb \
    -S Cuscuta_pedicellata-NC_052872.gb \
    -o results

# Query z lokalnego pliku GB, subject pobierany po numerze dostępu
python blast_and_show.py \
    -Q Laburnum_anagyroides-OZ176120.gb \
    -s NC_052872 \
    -o results \
    -e twój@email.adres

# Obie sekwencje pobierane po numerach dostępu
python blast_and_show.py \
    -q OZ176120 -s NC_052872 \
    -o results_gb \
    -e twój@email.adres

# Porównanie dwóch genomów mitochondrialnych (Cuscuta) — po accession
python blast_and_show.py \
    -q OZ176120 -s BK059237 \
    -o results_gb \
    -e twój@email.adres

# Dwa genomy plastydowe Cuscuta z lokalnych plików GB
python blast_and_show.py \
    -Q Cuscuta_pedicellata-NC_052872.gb \
    -S Cuscuta_epilinum-BK059237.gb \
    -o results

# Własny FASTA jako query, GB jako subject
python blast_and_show.py \
    -a 223-LaF4R4-223-LaF1R1.fasta \
    -S Cuscuta_epilinum-BK059237.gb \
    -o results

# Fragment genomu mitochondrialnego w FASTA przeciwko sekwencji z GenBank,
# z podwyższonym progiem e-value (100) — dopuszczamy słabsze dopasowania
python blast_and_show.py \
    -a Cuscuta_epithymum_mt_328700-340700.fasta \
    -s NC_070194 \
    -o results \
    -v 100 \
    -e twój@email.adres

# Użycie profilu konfiguracyjnego dla długich sekwencji
python blast_and_show.py \
    -q NC_052872 -s BK059237 \
    -o results \
    -e twój@email.adres \
    -c settings_long.json

# Użycie profilu dla krótkich sekwencji i własnego tytułu
python blast_and_show.py \
    -a region_5kb.fasta -b region_5kb_inny.fasta \
    -o results \
    -c settings_short.json \
    -t "Region Pytheas: Cuscuta vs. Orobanche"
```

---

## Pliki wyjściowe

Wszystkie pliki trafiają do `OUTPUT_DIR` (domyślnie: `output/`):

| Wzorzec pliku | Zawartość |
|---|---|
| `*.gb`, `*.fasta`, `*.json` | Pliki sekwencji query i subject (pobrane lub skopiowane) oraz dane features w JSON. |
| `*.tsv` | Eksport features wyekstrahowanych z każdego pliku GenBank, w formacie tabelarycznym. |
| `blast_results.tsv` | Wyniki BLAST w formacie tabelarycznym (`-outfmt 6`) z wierszem nagłówków. |
| `blast_results.xml` | Wyniki BLAST w formacie XML (`-outfmt 5`) do dalszego parsowania. |
| `aligned_query-*.fasta` | Dopasowane regiony sekwencji query w formacie FASTA. |
| `aligned_subject-*.fasta` | Dopasowane regiony sekwencji subject w formacie FASTA. |
| `out_draw.svg`, `out_draw.pdf`, `out_draw.png` | Wizualizacja w trzech formatach. |
| `settings.json` | Kopia konfiguracji użytej w danym uruchomieniu. |

Tymczasowe pliki bazy BLAST (`blast_db*`) są tworzone i usuwane automatycznie.

---

## Konfiguracja

Układ rysunku, czcionki, kolory oraz to, które typy features mają być
wyświetlane, kontrolowane są przez plik konfiguracyjny w formacie JSON.
W projekcie znajdują się trzy gotowe pliki konfiguracyjne:

- **`settings.json`** — profil domyślny, odpowiedni dla sekwencji rzędu kilkudziesięciu kb.
- **`settings_short.json`** — dopasowany do krótkich sekwencji (1–15 kb):
  większe czcionki, wyższe features, gęstsza podziałka w bp.
- **`settings_long.json`** — dopasowany do długich sekwencji (200–800 kb),
  takich jak genomy mitochondrialne roślin: mniejsze features, rzadsza
  podziałka, większa wysokość rysunku i mniej zaszumione etykiety.

Plik konfiguracyjny wskazujesz argumentem `-c PLIK`. Domyślnie `settings.json`.

### Struktura pliku

Plik konfiguracyjny składa się z dwóch głównych sekcji:

- **`settings`** — geometria, czcionki, marginesy, dane do wyświetlenia.
- **`feature_colors`** — kolory poszczególnych typów features (CDS, intron, tRNA itd.).

### Konwencje ogólne

- **Kolory** — listy trzech liczb `[R, G, B]` w zakresie `0.0–1.0`
  (format Cairo, nie `0–255`). Parametr `background` może mieć dodatkowy
  czwarty element — kanał alfa (przezroczystość), np. `[1, 1, 1, 1]` =
  nieprzezroczysta biel.
- **Przezroczystość (alpha)** — wartości od `0.0` (całkowicie przezroczyste)
  do `1.0` (pełna krycia).
- **Jednostki długości** — wszystkie wymiary geometryczne (marginesy,
  wysokości, długości) są podawane w **pikselach Cairo**. Rzeczywista
  szerokość rysunku w pikselach zależy dodatkowo od parametru `ratio`.
- **Czcionki** — `"Sans-Serif"` oznacza domyślną czcionkę bezszeryfową
  systemu. Można też podać konkretną nazwę, np. `"DejaVu Sans"` czy `"Arial"`.
  Jeśli czcionka nie jest dostępna, Cairo wybierze podstawienie.

### Parametry — sekcja `settings`

**Skala i przelicznik długości**

| Parametr | Wart. domyślna | Znaczenie |
|---|---|---|
| `ratio` | `50` | **Liczba par zasad na 1 piksel.** Większe = krótszy rysunek. Dla genomów mitochondrialnych roślin (~500 kb) typowo 50–200; dla krótkich sekwencji (~10 kb) zejdź do 5–10. |
| `tick` | `100` | **Krok podziałki osi w parach zasad.** Wartości typowe: 1 000, 10 000, 100 000, zależnie od długości sekwencji. |
| `tick_length` | `20` | Długość kreski podziałki (w pikselach), prostopadłej do linii sekwencji. |
| `tick_line_width` | `0.5` | Grubość kreski podziałki. |
| `tick_font` | `"Sans-Serif"` | Czcionka liczb przy podziałce. |
| `tick_font_size` | `2` | Rozmiar czcionki podziałki. **Uwaga:** domyślnie bardzo mała — etykiety są celowo dyskretne. Zwiększ, jeśli nieczytelne. |

**Marginesy i ogólne wymiary**

Wymiary rysunku liczone są jako:

```
draw_width  = left_margin + (długość_sekwencji / ratio) + right_margin
draw_height = top_margin  + main_height                  + bottom_margin
```

| Parametr | Wart. domyślna | Znaczenie |
|---|---|---|
| `left_margin` | `100` | Lewy margines — miejsce na etykiety `"query"`/`"subject"`. |
| `right_margin` | `300` | Prawy margines — mieści też legendę gradientu kolorów. |
| `top_margin` | `100` | Górny margines — nad tytułem i pierwszą linią sekwencji. |
| `bottom_margin` | `20` | Dolny margines. Zwiększ, jeśli features pod linią subject są ucinane. |
| `main_height` | `600` | Wysokość głównego obszaru rysunku. Określa miejsce na obie linie sekwencji i wszystkie features. |

**Pozycje linii query i subject**

| Parametr | Wart. domyślna | Znaczenie |
|---|---|---|
| `query_y_position` | `300` | Współrzędna Y (od góry) linii sekwencji **query**. Features rysowane nad linią. |
| `subject_y_position` | `500` | Współrzędna Y linii **subject**. Features rysowane pod linią. |

> Odstęp między liniami (`subject_y_position − query_y_position`) określa
> miejsce na linie łączące dopasowania. Zbyt mały → linie nakładają się
> i są nieczytelne.

**Features (geny, CDS, eksony, tRNA, rRNA, misc_feature)**

| Parametr | Wart. domyślna | Znaczenie |
|---|---|---|
| `feature_height` | `20` | Wysokość prostokąta pojedynczego features. |
| `feature_margin` | `65` | Pionowy odstęp od linii sekwencji do pasów z features. Features z listy `special_features` są rysowane w trzykrotnie większej odległości (`feature_margin × 3`) dla wyróżnienia. |
| `features_rendered` | `["CDS", "exon", "intron", "tRNA", "rRNA", "misc_feature"]` | **Lista typów features do narysowania.** Usuń typ z listy, żeby go ukryć. Przykłady: `["CDS", "tRNA", "rRNA"]` ukrywa eksony, introny i misc_feature; `["CDS"]` pokazuje tylko sekwencje kodujące; `["misc_feature"]` przydatne do podejrzenia samych regionów plastid-derived. |
| `special_features` | `["plastid-derived"]` | Lista features traktowanych jako „specjalne" — rysowane dalej od linii i w dedykowanym kolorze z `feature_colors`. |
| `feature_data` | `["name", "type", "start", "end"]` | **Pola wyświetlane w opisie features** nad każdym prostokątem. Kolejność ma znaczenie. Pole `end` poprzedzone myślnikiem (`start-end`), pozostałe spacją. |

**Alignmenty (dopasowania BLAST)**

| Parametr | Wart. domyślna | Znaczenie |
|---|---|---|
| `alignment_height` | `10` | Wysokość prostokąta pojedynczego dopasowania. |
| `alignment_alpha` | `0.9` | Krycie prostokątów dopasowań (`0.0`–`1.0`). |
| `alignment_line_alpha` | `0.3` | Krycie linii łączących dopasowania query/subject. Niższe = bardziej dyskretne. |
| `alignment_line_width` | `0.5` | Grubość linii łączącej. |
| `alignment_font` | `"Sans-Serif"` | Czcionka opisów dopasowań. |
| `alignment_font_size` | `5` | Rozmiar czcionki opisów dopasowań. |
| `query_alignment_data` | `["pident", "length", "gaps", "qstart", "qend"]` | **Pola wyświetlane w opisie dopasowania po stronie query.** `pident` jest automatycznie poprzedzane `%`, `qend`/`send` myślnikiem. |
| `subject_alignment_data` | `["pident", "length", "gaps", "sstart", "send"]` | Analogicznie po stronie subject. |

> **Dostępne pola** w `query_alignment_data` / `subject_alignment_data`
> odpowiadają kolumnom BLAST: `length`, `qlen`, `slen`, `qcovs`, `qcovhsp`,
> `pident`, `mismatch`, `gapopen`, `evalue`, `bitscore`, `gaps`, `qseqid`,
> `sseqid`, `qstart`, `qend`, `sstart`, `send`.

**Legenda skali kolorów (gradient procentowej identyczności)**

Dopasowania są kolorowane wg `pident` (% identyczności) gradientem „rainbow"
od `min_scale` do 100 %.

| Parametr | Wart. domyślna | Znaczenie |
|---|---|---|
| `min_scale` | `70` | **Dolna granica skali kolorów** (procent pident). Dopasowania o `pident < min_scale` nie mają przypisanego koloru — wartość musi być ≤ minimalnego pident w wynikach BLAST, inaczej skrypt rzuci `KeyError`. |
| `scale_move_left` | `100` | Pozioma pozycja legendy — odległość od **prawej krawędzi**. |
| `scale_move_down` | `50` | Pionowa pozycja legendy — odległość od górnej krawędzi. |
| `scale_step_height` | `10` | Wysokość jednego stopnia gradientu w legendzie. |
| `scale_step_width` | `10` | Szerokość jednego stopnia gradientu w legendzie. |
| `scale_font` | `"Sans-Serif"` | Czcionka liczb przy legendzie. |
| `scale_font_size` | `8` | Rozmiar czcionki liczb przy legendzie. |
| `scale_font_color` | `[0.1, 0.1, 0.1]` | Kolor liczb legendy (prawie czarny). |

**Tytuł i tło**

| Parametr | Wart. domyślna | Znaczenie |
|---|---|---|
| `title_font` | `"Sans-Serif"` | Czcionka tytułu. |
| `title_font_size` | `12` | Rozmiar czcionki tytułu. |
| `background` | `[1, 1, 1, 1]` | Kolor tła w formacie RGBA `[R, G, B, alpha]`. Domyślnie biel. Dla przezroczystego PNG: `[1, 1, 1, 0]`. |

### Parametry — sekcja `feature_colors`

Kolory features w formacie RGB (`[R, G, B]`, każdy w zakresie `0.0–1.0`).

| Klucz | Wart. domyślna | Znaczenie |
|---|---|---|
| `alpha` | `0.5` | **Wspólna przezroczystość** dla wszystkich features. Stosowana zarówno do wypełnienia prostokąta, jak i opisu. |
| `CDS` | `[0, 0.8, 0]` | Zielony — sekwencje kodujące. |
| `exon` | `[0, 0.8, 0]` | Zielony — eksony (jak CDS). |
| `intron` | `[0.6, 0.8, 0.6]` | Bladozielony — introny. |
| `tRNA` | `[0, 0, 0.9]` | Niebieski — geny tRNA. |
| `rRNA` | `[0.4, 0.4, 0]` | Oliwkowy — geny rRNA. |
| `misc_feature` | `[0.8, 0.8, 0.8]` | Jasnoszary — pozostałe oznaczenia. |
| `plastid-derived` | `[1, 0, 0]` | Czerwony — **misc_features pochodzące z plastydu** (HGT). Rysowane dalej od linii sekwencji. |

> **Dodawanie własnych typów:** żeby pokolorować nowy typ features (np. `repeat_region`),
> dodaj klucz w `feature_colors` oraz dopisz typ do listy `features_rendered`.
> Dodatkowo trzeba rozszerzyć listę akceptowanych typów w funkcji
> `export_features_from_gb_to_tsv()` w `dataimport.py`.

### Typowe scenariusze edycji

**Rysunek wychodzi za szeroki**
Zwiększ `ratio` (np. `50` → `100` — 2× kompresja w poziomie). Następnie
proporcjonalnie zwiększ `tick`, żeby etykiet podziałki nie było za dużo.

**Rysunek wychodzi za krótki / mało szczegółowy**
Zmniejsz `ratio` (np. `50` → `20`). Zmniejsz `tick`, żeby gęstsza podziałka
była czytelna.

**Linie query i subject zachodzą na features**
Zwiększ `feature_margin` (np. `65` → `100`), lub rozsuń linie zmieniając
`query_y_position` i `subject_y_position`.

**Features pod linią subject są ucinane**
Zwiększ `bottom_margin` (np. `20` → `100`), lub zwiększ `main_height`, jeśli
chcesz zostawić więcej miejsca także na features powyżej.

**Etykiety dopasowań są nieczytelne**
Zwiększ `alignment_font_size` (np. `5` → `8`), lub uprość zawartość usuwając
pola z `query_alignment_data` / `subject_alignment_data`.

**Zbyt wiele dopasowań → bałagan kolorystyczny**
Filtruj wyniki już na etapie wczytywania (argument `-m` / `--minlen`
w `blast_and_show.py`). Zmniejsz `alignment_line_alpha` (np. do `0.1`),
żeby linie łączące były dyskretniejsze.

**Legenda gradientu wychodzi poza obraz**
Zwiększ `right_margin` (np. `300` → `400`), lub zmniejsz `scale_move_left`
(np. `100` → `60`), żeby legenda była bliżej krawędzi.

**Chcę pokazać tylko wybrane typy elementów**
Edytuj `features_rendered`. Na przykład żeby zobaczyć same regiony
plastid-derived bez tła z genów: `"features_rendered": ["misc_feature"]`.

**Błąd `KeyError: <liczba>` podczas rysowania dopasowań**
Jakieś dopasowanie ma `pident < min_scale`. Zmniejsz `min_scale`
(np. `70` → `60`).

### Edytowanie pliku

Plik jest standardowym JSON-em — można go edytować dowolnym edytorem
tekstu. Po edycji warto sprawdzić poprawność składni:

```bash
python -m json.tool settings.json
```

Polecenie wypisze plik z formatowaniem i zwróci błąd, jeśli składnia jest
nieprawidłowa (brakujący przecinek, niedomknięty nawias itp.).

---

## Rozwiązywanie problemów

**`sh: makeblastdb: nie znaleziono polecenia`**
NCBI BLAST+ nie jest zainstalowany lub nie jest w `PATH`. Zainstaluj poprzez
`sudo dnf install ncbi-blast+` (Fedora) lub
`sudo apt install ncbi-blast+` (Debian/Ubuntu).

**`ModuleNotFoundError: No module named 'cairo'`**
Brakuje PyCairo. Zainstaluj poprzez `sudo dnf install python3-cairo` (Fedora)
albo `pip install pycairo` po wcześniejszym zainstalowaniu nagłówków Cairo
(`libcairo2-dev` na Debianie/Ubuntu, `cairo-devel` na Fedorze).

**`KeyError: <liczba>` podczas rysowania dopasowań**
Jakieś dopasowanie ma `pident` mniejszy niż `min_scale` (domyślnie 70).
Zmniejsz `min_scale` w pliku konfiguracyjnym.

**Features pod linią subject są ucinane**
Zwiększ `bottom_margin` (lub `main_height`) w pliku konfiguracyjnym.
