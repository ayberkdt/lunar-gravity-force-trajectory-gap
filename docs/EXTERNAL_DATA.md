# External data products

None of these are redistributed here. `fetch_data.py` downloads all of them
into the expected paths and verifies each gravity file against the SHA-256 the
campaign manifest recorded, so a superseded or corrupted product is caught on
arrival rather than becoming a wrong number later.

```bash
python fetch_data.py --list     # status and source of every product
python fetch_data.py            # fetch everything (~430 MB)
```

Every URL below was checked against the live archive. Where a product is
served from more than one place, the source listed is the one whose bytes match
the recorded digest — for several of these the PDS copy and the Tudat resource
copy are **not** byte-identical, and only the listed one verifies.

## Gravity fields

Expected path: `data/gravity_models/`

| Product | Filename | Bytes | Source | SHA-256 |
|---|---|---|---|---|
| GRAIL JGGRX_1800F (primary field) | `jggrx_1800f_sha.tab` | 197,969,644 | [PDS Geosciences](https://pds-geosciences.wustl.edu/grail/grail-l-lgrs-5-rdr-v1/grail_1001/shadr/jggrx_1800f_sha.tab) | `d2a552067a78bf1d2755807ae14ee1d6843a8f6a4228e01ce59a665516738fec` |
| GSFC GRGM1200A | `gggrx_1200a_sha.tab` | 88,059,844 | [PDS Geosciences](https://pds-geosciences.wustl.edu/grail/grail-l-lgrs-5-rdr-v1/grail_1001/shadr/gggrx_1200a_sha.tab) | recorded in `metrics/r16_final_experiment_manifest.json` |
| GSFC GGGRX_1200L | `gggrx_1200l_sha.tab` | 87,913,291 | [Tudat resources](https://raw.githubusercontent.com/tudat-team/tudat-resources/master/resource/gravity_models/Moon/gggrx_1200l_sha.tab) | recorded in `metrics/r16_final_experiment_manifest.json` |
| Earth GOCO05c | `GOCO05c.txt` | 23,100,684 | [Tudat resources](https://raw.githubusercontent.com/tudat-team/tudat-resources/master/resource/gravity_models/Earth/GOCO05c.txt) | recorded in `metrics/r16_final_experiment_manifest.json` |
| Earth EGM96 | `egm96.txt` | 3,201,589 | [Tudat resources](https://raw.githubusercontent.com/tudat-team/tudat-resources/master/resource/gravity_models/Earth/egm96.txt) | recorded in `metrics/r16_final_experiment_manifest.json` |
| Mars JGMRO120D | `jgmro120d.txt` | 893,027 | [Tudat resources](https://raw.githubusercontent.com/tudat-team/tudat-resources/master/resource/gravity_models/Mars/jgmro120d.txt) | recorded in `metrics/r16_final_experiment_manifest.json` |
| Mercury JGMESS_160A | `jgmess_160a_sha.tab` | 1,591,122 | [Tudat resources](https://raw.githubusercontent.com/tudat-team/tudat-resources/master/resource/gravity_models/Mercury/jgmess_160a_sha.tab) | recorded in `metrics/r16_final_experiment_manifest.json` |
| Venus SHGJ180U | `shgj180u.a01` | 2,009,581 | [Tudat resources](https://raw.githubusercontent.com/tudat-team/tudat-resources/master/resource/gravity_models/Venus/shgj180u.a01) | recorded in `metrics/r16_final_experiment_manifest.json` |

The digests are not duplicated into this table on purpose. `fetch_data.py`
reads them from `metrics/r16_final_experiment_manifest.json` at run time, so
there is one authority for each value and no chance of a transcription drift
between the manifest and the documentation. The primary field's digest is
quoted above because the manuscript quotes it too.

Lunar Prospector LP150Q (`data/jgl150q1.sha`, 1,400,194 bytes) is small enough
to ship and is already in the repository; it is what `verify_snapshot.py` uses,
so the quick check needs no download. Its archive copy is at
[PDS Geosciences](https://pds-geosciences.wustl.edu/lunar/lp-l-rss-5-gravity-v1/lp_1001/sha/jgl150q1.sha).

## SPICE kernels

Expected path: `data/spice_kernels/`

| Kernel | Filename | Bytes | Source |
|---|---|---|---|
| Planetary ephemeris | `de440s.bsp` | 32,726,016 | [NAIF](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp) |
| Leap seconds | `naif0012.tls` | 5,257 | [NAIF](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls) |
| Gravitational parameters | `gm_de440.tpc` | 12,406 | [NAIF](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de440.tpc) |
| Lunar orientation | `moon_pa_de440_200625.bpc` | 12,863,488 | [NAIF](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de440_200625.bpc) |
| Lunar reference frames | `moon_de440_250416.tf` | 19,478 | [NAIF](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/moon_de440_250416.tf) |

The campaign manifests record no digests for the SPICE kernels, so
`fetch_data.py` cannot verify them. It prints the SHA-256 it computed for each
so the values can be pinned downstream. As fetched on 2026-07-29:

| Filename | SHA-256 |
|---|---|
| `de440s.bsp` | `c1c7feeab882263fc493a9d5a5b2ddd71b54826cdf65d8d17a76126b260a49f2` |
| `naif0012.tls` | `678e32bdb5a744117a467cd9601cd6b373f0e9bc9bbde1371d5eee39600a039b` |
| `gm_de440.tpc` | `924ddf4fb9ead9fe8a1aa55780bcabde40b09d00065d58226e24b68d8092f140` |
| `moon_pa_de440_200625.bpc` | `60cd55aa401ea2ea97360636f567554bfe4e37bb829f901b4460a455dfaf783f` |
| `moon_de440_250416.tf` | `a47c71e9c9f33796bdafb2c9d69a7ee447b6016ecad80f71cd6f3e479f9cf768` |

These are a record of what the archives served on that date, not a claim that
the campaigns ran against exactly these bytes.

## Licensing

The PDS and NAIF products are public domain NASA data. The Tudat resource
archive is distributed under its own terms; consult that project before
redistributing its copies.
