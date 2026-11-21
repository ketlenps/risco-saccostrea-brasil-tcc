# Teste de sensibilidade da escala de decaimento da pressão (DECAIMENTO_KM do R).
# Refaz só a pressão e o risco do Módulo 11 com 25, 50 e 100 km, usando a
# adequabilidade e as infraestruturas da rodada final (não refaz a modelagem):
#
# Conta igual à do R:
#   - cada infraestrutura cai na célula de 0,05° que a contém;
#   - d = distância (km) do centro da célula da faixa costeira até a célula com
#     infraestrutura mais próxima, separada por tipo (portos, poços, plataformas);
#   - p = exp(−d / L) para cada tipo; P = média dos três; R = A × P;
#   - classes: Baixo A < 0,542; Médio A ≥ 0,542 e P < e^−1; Alto A ≥ 0,542 e P ≥ e^−1.
# Primeiro confere se L = 50 km reproduz o pressao_antropica.tif do R.
#
# Saída: <RESULTADOS>\sensibilidade_decaimento\ (tabelas .csv e figura .png)
#   Rode com o Python que vem com o QGIS (3.44), na pasta do repositório, por exemplo:
#   "C:\Program Files\QGIS 3.44.9\bin\python-qgis-ltr.bat" scripts_python\sensibilidade_decaimento.py
import csv, os
import numpy as np
from osgeo import gdal, ogr
from scipy import ndimage
from scipy.spatial import cKDTree
gdal.UseExceptions(); ogr.UseExceptions()

AQUI = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(AQUI)
RES = os.path.join(BASE, os.environ.get("RESULTADOS", "resultados"))
SAIDA = os.path.join(RES, "sensibilidade_decaimento")
GPKG_ECO = os.path.join(BASE, "vetores", "ecorregioes_marinhas.gpkg")
LIMIAR = 0.542
ESCALAS = (25, 50, 100)
TIPOS = ("portos", "pocos", "plataformas")
R_TERRA = 6371.0  # km
os.makedirs(SAIDA, exist_ok=True)


def ler(nome):
    ds = gdal.Open(os.path.join(RES, "rasters", nome + ".tif"))
    b = ds.GetRasterBand(1)
    a = b.ReadAsArray().astype(float)
    nd = b.GetNoDataValue()
    if nd is not None:
        a[a == nd] = np.nan
    return ds, a


def xyz(lon, lat):
    lo, la = np.radians(lon), np.radians(lat)
    return np.column_stack((np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)))


ref, adeq = ler("adequabilidade_brasil")
_, pressao_r = ler("pressao_antropica")
_, risco_r = ler("risco_continuo")
gt = ref.GetGeoTransform(); nx, ny = ref.RasterXSize, ref.RasterYSize
faixa = ~np.isnan(adeq)
ii, jj = np.where(faixa)
lon_c = gt[0] + gt[1] * (jj + 0.5); lat_c = gt[3] + gt[5] * (ii + 0.5)
arvore_cel = xyz(lon_c, lat_c)

# infraestruturas usadas no R (as mesmas da tabela 12), colocadas no centro da sua célula
pts = {t: [] for t in TIPOS}
with open(os.path.join(RES, "tabelas", "12_risco_por_infraestrutura.csv"), encoding="utf-8") as fh:
    for r in csv.DictReader(fh):
        pts[r["tipo"]].append((float(r["lon"]), float(r["lat"])))
dist = {}
for t in TIPOS:
    p = np.array(pts[t])
    cj = np.floor((p[:, 0] - gt[0]) / gt[1]); ci = np.floor((p[:, 1] - gt[3]) / gt[5])
    cel = np.unique(np.column_stack((ci, cj)), axis=0)
    lon_p = gt[0] + gt[1] * (cel[:, 1] + 0.5); lat_p = gt[3] + gt[5] * (cel[:, 0] + 0.5)
    corda, _ = cKDTree(xyz(lon_p, lat_p)).query(arvore_cel)
    dist[t] = 2 * R_TERRA * np.arcsin(np.clip(corda / 2, 0, 1))  # km pelo arco
    print(f"{t:12s} {len(p):5d} pontos em {len(cel):5d} células")

# ecorregião de cada célula (mesma regra da tabela 14)
eco_ds = ogr.Open(GPKG_ECO); camada = eco_ds.GetLayerByName("ecorregioes_marinhas")
nomes = {f.GetField("ordem"): f.GetField("ecorregiao") for f in camada}
mem = gdal.GetDriverByName("MEM").Create("", nx, ny, 1, gdal.GDT_Int16)
mem.SetGeoTransform(gt); mem.SetProjection(ref.GetProjection())
gdal.RasterizeLayer(mem, [1], camada, options=["ATTRIBUTE=ordem"])
ordem = mem.ReadAsArray()
_, (ai, aj) = ndimage.distance_transform_edt(ordem == 0, return_indices=True)
eco = ordem[ai, aj][faixa]

A = adeq[faixa]
adequado = A >= LIMIAR - 1e-6


def calcular(L):
    P = np.mean([np.exp(-dist[t] / L) for t in TIPOS], axis=0)
    R = A * P
    classe = np.where(~adequado, 1, np.where(P < np.exp(-1), 2, 3))
    return P, R, classe


# 1) conferência com o R
P50, R50, C50 = calcular(50)
dif_p = np.abs(P50 - pressao_r[faixa]); dif_r = np.abs(R50 - risco_r[faixa])
print(f"\nconferência L = 50 km x R: pressão |dif| média {np.nanmean(dif_p):.4f}, máx {np.nanmax(dif_p):.4f}; "
      f"risco |dif| média {np.nanmean(dif_r):.4f}; correlação da pressão {np.corrcoef(P50, pressao_r[faixa])[0, 1]:.4f}")

# 2) as três escalas
res = {L: calcular(L) for L in ESCALAS}
n = len(A)
q = np.quantile(R50, 1 - np.mean(C50 == 3))  # mesmo nº de células do "alto" de 50 km


def rank(x):
    return np.argsort(np.argsort(x))


linhas = []
for L in ESCALAS:
    P, R, C = res[L]
    top = R >= np.quantile(R, 1 - np.mean(C50 == 3))
    linhas.append({
        "decaimento_km": L,
        "pressao_media": round(float(P.mean()), 3),
        "risco_medio": round(float(R.mean()), 3),
        "pct_baixo": round(100 * float(np.mean(C == 1)), 1),
        "pct_medio": round(100 * float(np.mean(C == 2)), 1),
        "pct_alto": round(100 * float(np.mean(C == 3)), 1),
        "spearman_risco_vs_50km": round(float(np.corrcoef(rank(R), rank(R50))[0, 1]), 3),
        "pct_alto_50km_que_segue_alto": round(100 * float(np.mean(C[C50 == 3] == 3)), 1),
        "pct_maiores_riscos_em_comum_com_50km": round(100 * float(np.mean(top[R50 >= q])), 1),
    })
with open(os.path.join(SAIDA, "S1_sensibilidade_geral.csv"), "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(linhas[0])); w.writeheader(); w.writerows(linhas)

# por ecorregião: % de risco alto e posição no ranking do risco médio
eco_linhas = []
for o in sorted(nomes):
    m = eco == o
    if not m.any():
        continue
    l = {"ecorregiao": nomes[o], "celulas": int(m.sum())}
    for L in ESCALAS:
        P, R, C = res[L]
        l[f"risco_medio_{L}km"] = round(float(R[m].mean()), 3)
        l[f"pct_alto_{L}km"] = round(100 * float(np.mean(C[m] == 3)), 1)
    eco_linhas.append(l)
for L in ESCALAS:
    ordem_eco = sorted(eco_linhas, key=lambda l: -l[f"risco_medio_{L}km"])
    for k, l in enumerate(ordem_eco, 1):
        l[f"posicao_{L}km"] = k
with open(os.path.join(SAIDA, "S2_sensibilidade_por_ecorregiao.csv"), "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(eco_linhas[0])); w.writeheader(); w.writerows(eco_linhas)

# rasters de classes para o QGIS (1 Baixo, 2 Médio, 3 Alto)
for L in ESCALAS:
    a = np.zeros((ny, nx), dtype=np.uint8); a[faixa] = res[L][2]
    d = gdal.GetDriverByName("GTiff").Create(os.path.join(SAIDA, f"risco_classes_{L}km.tif"), nx, ny, 1,
                                             gdal.GDT_Byte, options=["COMPRESS=DEFLATE"])
    d.SetGeoTransform(gt); d.SetProjection(ref.GetProjection())
    d.GetRasterBand(1).WriteArray(a); d.GetRasterBand(1).SetNoDataValue(0); d = None

print("\nGERAL")
for l in linhas:
    print(l)
print("\nPOR ECORREGIÃO")
for l in eco_linhas:
    print(l["ecorregiao"][:34].ljust(34), "  ".join(
        f"{L}km: risco {l[f'risco_medio_{L}km']:.3f} alto {l[f'pct_alto_{L}km']:5.1f}% #{l[f'posicao_{L}km']}" for L in ESCALAS))
print("ok", SAIDA)
