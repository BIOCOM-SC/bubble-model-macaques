# data/ - experimental data and model targets

## Raw data (`raw/`)

- `NogueiraSuppl_File_2_RawData.xlsx` - CT measurements of every lesion in the rhesus and
  cynomolgus macaques of Nogueira, I. *et al.* Pathogens **12**, 236 (2023),
  [doi:10.3390/pathogens12020236](https://doi.org/10.3390/pathogens12020236), Supplementary File S2
  (CC BY 4.0), unchanged. One sheet per animal with the tracked lesions (lobe, segment, and per CT:
  number of daughter lesions, volume, main and short axis), and one sheet per animal with the
  isolated micronodules per segment ("Isol Mic").
- `lung_volumes.csv` - lung volume of each animal at the first CT (week 3), in cm^3, used for the
  volume fraction of consolidations (VFC). Provided by the authors of Nogueira et al. (2023).

## Processing

```sh
python data/process.py                 # unvaccinated rhesus macaques (RUV) -> processed/
python data/process.py --cohort RIDV   # intradermal-BCG-vaccinated macaques -> processed/ridv/
python data/process.py --help
```

Animals and time points: RUV01-RUV11 without RUV05 (miliary TB, outside the scope of the model),
n = 10, CTs at days 21, 53, 77; RIDV01-RIDV05, CTs at days 21, 49, 77 (Suppl. Section S5).

Per animal and CT:

- **NC**: tracked lesions with main axis >= 4.5 mm (consolidations).
- **NMD**: daughter micronodules of the tracked lesions (counts written as `10+` are read as 10;
  `--plus-value` changes this).
- **NMI**: isolated micronodules, plus tracked lesions that are present but smaller than 4.5 mm.
- **NM** = NMD + NMI.
- **VFC**: total volume of the consolidations / lung volume.
- **MAC**: mean main axis of the consolidations (undefined without consolidations; averaged over
  the animals that have at least one).

The targets are the mean and the sample SD over the animals.

## Outputs (`processed/`, `processed/ridv/`)

| file | content |
|---|---|
| `experimental_results.txt` | model input: mean and SD lines of NC, NM, NMD, NMI, VFC, MAC |
| `summary.txt` | readable version (Table 2 of the paper), lobe frequencies, SA/MA, lung volumes |
| `per_animal.csv` | one line per animal and CT |
| `consolidations.csv` | one line per tracked lesion and CT |
| `lobe_frequency_week3.csv` | share of the week-3 consolidations in each lobe |
