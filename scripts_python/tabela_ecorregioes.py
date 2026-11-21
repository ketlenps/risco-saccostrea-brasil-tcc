# Tabela de risco por ECORREGIÃO marinha (Spalding et al., 2007), no lugar das
# regiões do IBGE: mesmas colunas do 12_risco_por_regiao.csv do R + infraestruturas
# por classe de risco (como a Tabela 3 do TCC de referência).
# Variável RESULTADOS = pasta da rodada (padrão: resultados).
# Saída: <RESULTADOS>\tabelas\14_risco_por_ecorregiao.csv e 14_infraestruturas_por_ecorregiao.csv
#
# Cada célula de 0,05° da faixa costeira entra na ecorregião que contém o centro
# dela; as poucas que caem fora dos polígonos (a faixa entra um pouco em terra)
# vão para a ecorregião mais próxima. Infraestruturas: o ponto entra na
# ecorregião que o contém ou, se estiver em terra (portos), na mais próxima.
#   Rode com o Python que vem com o QGIS (3.44), na pasta do repositório, por exemplo:
#   "C:\Program Files\QGIS 3.44.9\bin\python-qgis-ltr.bat" scripts_python\tabela_ecorregioes.py
import csv, os
import numpy as np
from osgeo import gdal, ogr
from scipy import ndimage
gdal.UseExceptions(); ogr.UseExceptions()

AQUI = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(AQUI)
RES = os.path.join(BASE, os.environ.get("RESULTADOS", "resultados"))
GPKG = os.path.join(BASE, "vetores", "ecorregioes_marinhas.gpkg")
LIMIAR = float(os.environ.get("LIMIAR", "0.542"))


def ler(nome):
    ds = gdal.Open(os.path.join(RES, "rasters", nome + ".tif"))
    b = ds.GetRasterBand(1)
    a = b.ReadAsArray().astype(float)
    nd = b.GetNoDataValue()
    if nd is not None:
        a[a == nd] = np.nan
    return ds, a


ref, adeq = ler("adequabilidade_brasil")
_, pressao = ler("pressao_antropica")
_, risco = ler("risco_continuo")
_, classe = ler("risco_classes")
valida = ~np.isnan(adeq)

# ecorregião de cada célula (campo "ordem")
eco_ds = ogr.Open(GPKG)
camada = eco_ds.GetLayerByName("ecorregioes_marinhas")
nomes = {f.GetField("ordem"): f.GetField("ecorregiao") for f in camada}
mem = gdal.GetDriverByName("MEM").Create("", ref.RasterXSize, ref.RasterYSize, 1, gdal.GDT_Int16)
mem.SetGeoTransform(ref.GetGeoTransform()); mem.SetProjection(ref.GetProjection())
gdal.RasterizeLayer(mem, [1], camada, options=["ATTRIBUTE=ordem"])
ordem = mem.ReadAsArray()
sem = ordem == 0
_, (ii, jj) = ndimage.distance_transform_edt(sem, return_indices=True)
ordem_completa = ordem[ii, jj]
fora = int((valida & sem).sum())

linhas = []
for o in sorted(nomes):
    m = valida & (ordem_completa == o)
    n = int(m.sum())
    if n == 0:
        continue
    linhas.append({
        "ecorregiao": nomes[o], "celulas": n,
        "adeq_media": round(float(np.nanmean(adeq[m])), 3),
        "pct_costa_adequada": round(100 * float(np.mean(adeq[m] >= LIMIAR - 1e-6)), 1),
        "pressao_media": round(float(np.nanmean(pressao[m])), 3),
        "risco_medio": round(float(np.nanmean(risco[m])), 3),
        "pct_risco_alto": round(100 * float(np.mean(classe[m] == 3)), 1),
    })

# infraestruturas
pol = {}
for f in camada:
    pol[f.GetField("ordem")] = f.GetGeometryRef().Clone()
cont = {}
with open(os.path.join(RES, "tabelas", "12_risco_por_infraestrutura.csv"), encoding="utf-8") as fh:
    for r in csv.DictReader(fh):
        if r["classe"] in ("", "NA"):
            continue
        p = ogr.Geometry(ogr.wkbPoint); p.AddPoint_2D(float(r["lon"]), float(r["lat"]))
        dentro = [o for o, g in pol.items() if g.Contains(p)]
        o = dentro[0] if dentro else min(pol, key=lambda k: pol[k].Distance(p))
        chave = (nomes[o], r["tipo"], r["classe"])
        cont[chave] = cont.get(chave, 0) + 1

tipos = ("portos", "pocos", "plataformas")
for l in linhas:
    for t in tipos:
        for c in ("Baixo", "Médio", "Alto"):
            l[f"{t}_{c.lower().replace('é', 'e')}"] = cont.get((l["ecorregiao"], t, c), 0)

saida = os.path.join(RES, "tabelas", "14_risco_por_ecorregiao.csv")
with open(saida, "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(linhas[0].keys()))
    w.writeheader(); w.writerows(linhas)

print(f"células fora dos polígonos (foram para a ecorregião mais próxima): {fora}")
print(f"{'ecorregião':38s} {'cél':>5s} {'adeq':>5s} {'%adeq':>6s} {'press':>6s} {'risco':>6s} {'%alto':>6s}  portos A/M/B  poços A/M/B  plat A/M/B")
for l in linhas:
    print(f"{l['ecorregiao']:38s} {l['celulas']:5d} {l['adeq_media']:5.3f} {l['pct_costa_adequada']:6.1f} "
          f"{l['pressao_media']:6.3f} {l['risco_medio']:6.3f} {l['pct_risco_alto']:6.1f}  "
          + "  ".join(f"{l[t+'_alto']}/{l[t+'_medio']}/{l[t+'_baixo']}" for t in tipos))
print("total de células:", sum(l["celulas"] for l in linhas), "| ok", saida)
