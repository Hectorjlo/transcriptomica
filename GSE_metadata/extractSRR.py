import pandas as pd

META = "GSE132040_MACA_Bulk_metadata.csv"

TISSUE = "Bone"
AGES = {3, 24}
SEX = "m"  # machos

df = pd.read_csv(META)

# Normaliza nombres de columnas (porque tienen espacios/dos puntos)
col_source = "source name"
col_age = "characteristics: age"
col_sex = "characteristics: sex"
col_srr = "raw file"

# Validación rápida
missing = [c for c in [col_source, col_age, col_sex, col_srr] if c not in df.columns]
if missing:
    raise SystemExit(f"Faltan columnas en el CSV: {missing}\nColumnas disponibles: {df.columns.tolist()}")

# Limpieza / casting
df[col_sex] = df[col_sex].astype(str).str.strip().str.lower()
df[col_age] = pd.to_numeric(df[col_age], errors="coerce")

# Filtrado:
# 1) source name empieza con "Bone_" (evita cosas raras)
# 2) age en {3,24}
# 3) sex == "m"
mask = (
    df[col_source].astype(str).str.startswith(f"{TISSUE}_") &
    df[col_age].isin(AGES) &
    (df[col_sex] == SEX)
)

sub = df.loc[mask, [col_source, col_age, col_sex, col_srr]].copy()

# Ordena para que quede bonito
sub = sub.sort_values([col_age, col_source, col_srr])

# Extrae SRR
srrs = sub[col_srr].astype(str).str.strip()

# Reporte
print("=== Resumen (Hector: Bone, 3 vs 24, machos) ===")
print(sub.to_string(index=False))
print("\n=== SRR (uno por línea) ===")
for s in srrs:
    print(s)

# También guarda a archivo por conveniencia
out_txt = "hector_bone_3_vs_24_srr.txt"
srrs.to_csv(out_txt, index=False, header=False)
print(f"\nGuardado: {out_txt}  (lista de SRR)")