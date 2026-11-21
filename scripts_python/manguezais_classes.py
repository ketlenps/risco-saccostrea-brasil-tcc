# Manguezais (Atlas dos Manguezais do Brasil, ICMBio 2018) por classe de risco,
# com a MESMA conta do Módulo 12B do R, para o mapa do apêndice no QGIS:
#   cada célula de 0,05° com manguezal recebe o MAIOR valor de adequabilidade e de
#   pressão da faixa costeira num raio de 2 células (~11 km, janela 5 × 5);
#   classe: Sem valor (sem faixa por perto), Baixo (adeq < limiar), Médio
#   (adequado e pressão < e^-1), Alto (adequado e pressão >= e^-1).
# Saída: resultados\vetores\manguezais_classes.gpkg (camada
# manguezais_classes: um quadrado de 0,05° por célula, campos classe, codigo,
# area_km2) e a soma por classe, para conferir com 12B_manguezais_por_estado.csv.
#   Rode com o Python que vem com o QGIS (3.44), na pasta do repositório, por exemplo:
#   "C:\Program Files\QGIS 3.44.9\bin\python-qgis-ltr.bat" scripts_python\manguezais_classes.py
import os
import numpy as np
from osgeo import gdal, ogr, osr
from scipy import ndimage
gdal.UseExceptions(); ogr.UseExceptions()

AQUI = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(AQUI)
RES = os.path.join(BASE, os.environ.get("RESULTADOS", "resultados"), "rasters")
MANGUES = os.path.join(BASE, "dados", "habitats", "Mangues.shp")
SAIDA = os.path.join(BASE, "resultados", "vetores", "manguezais_classes.gpkg")
LIMIAR = 0.542 - 1e-6
RAIO = 2
FINO = 10  # subdivisões por célula para calcular a fração coberta


def ler(nome):
    ds = gdal.Open(os.path.join(RES, nome))
    a = ds.GetRasterBand(1).ReadAsArray().astype(float)
    nd = ds.GetRasterBand(1).GetNoDataValue()
    if nd is not None:
        a[a == nd] = np.nan
    return ds, a


def maximo_vizinho(a):
    b = np.where(np.isnan(a), -np.inf, a)
    m = ndimage.maximum_filter(b, size=2 * RAIO + 1, mode="constant", cval=-np.inf)
    return np.where(np.isinf(m), np.nan, m)


ref, adeq = ler("adequabilidade_brasil.tif")
_, pressao = ler("pressao_antropica.tif")
gt = ref.GetGeoTransform(); nx, ny = ref.RasterXSize, ref.RasterYSize
adeq_v, pres_v = maximo_vizinho(adeq), maximo_vizinho(pressao)

# fração de cada célula coberta por manguezal (rasteriza numa grade 10× mais fina)
fino = gdal.GetDriverByName("MEM").Create("", nx * FINO, ny * FINO, 1, gdal.GDT_Byte)
fino.SetGeoTransform((gt[0], gt[1] / FINO, 0, gt[3], 0, gt[5] / FINO))
srs = osr.SpatialReference(); srs.ImportFromEPSG(4326)
srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
fino.SetProjection(srs.ExportToWkt())
mg = ogr.Open(MANGUES)
gdal.RasterizeLayer(fino, [1], mg.GetLayer(), burn_values=[1])
fracao = fino.ReadAsArray().reshape(ny, FINO, nx, FINO).mean(axis=(1, 3))
fino = None

# área da célula em km² (esfera do WGS84, como o cellSize do terra)
lat = gt[3] + gt[5] * (np.arange(ny) + 0.5)
km_lat = abs(gt[5]) * 110.574
km_lon = gt[1] * 111.320 * np.cos(np.radians(lat))
area = fracao * (km_lat * km_lon)[:, None]

com = area > 0
classe = np.full(adeq.shape, "", dtype=object)
classe[com & np.isnan(adeq_v)] = "Sem valor"
ok = com & ~np.isnan(adeq_v)
classe[ok & (adeq_v < LIMIAR)] = "Baixo"
adequado = ok & (adeq_v >= LIMIAR)
classe[adequado & (pres_v < np.exp(-1))] = "Médio"
classe[adequado & (pres_v >= np.exp(-1))] = "Alto"
codigo = {"Baixo": 1, "Médio": 2, "Alto": 3, "Sem valor": 4}

if os.path.exists(SAIDA):
    os.remove(SAIDA)
os.makedirs(os.path.dirname(SAIDA), exist_ok=True)
ds = ogr.GetDriverByName("GPKG").CreateDataSource(SAIDA)
l = ds.CreateLayer("manguezais_classes", srs, ogr.wkbPolygon)
for campo, t in (("classe", ogr.OFTString), ("codigo", ogr.OFTInteger), ("area_km2", ogr.OFTReal),
                 ("adequabilidade", ogr.OFTReal), ("pressao", ogr.OFTReal)):
    l.CreateField(ogr.FieldDefn(campo, t))
soma = {}
l.StartTransaction()
for i, j in zip(*np.where(com)):
    x0 = gt[0] + gt[1] * j; y0 = gt[3] + gt[5] * i
    anel = ogr.Geometry(ogr.wkbLinearRing)
    for x, y in ((x0, y0), (x0 + gt[1], y0), (x0 + gt[1], y0 + gt[5]), (x0, y0 + gt[5]), (x0, y0)):
        anel.AddPoint_2D(x, y)
    p = ogr.Geometry(ogr.wkbPolygon); p.AddGeometry(anel)
    f = ogr.Feature(l.GetLayerDefn())
    c = classe[i, j]
    f.SetField("classe", c); f.SetField("codigo", codigo[c]); f.SetField("area_km2", round(float(area[i, j]), 3))
    if not np.isnan(adeq_v[i, j]):
        f.SetField("adequabilidade", round(float(adeq_v[i, j]), 3))
        f.SetField("pressao", round(float(pres_v[i, j]), 3))
    f.SetGeometry(p); l.CreateFeature(f)
    soma[c] = soma.get(c, 0) + float(area[i, j])
l.CommitTransaction()
ds = None
total = sum(soma.values())
print("células com manguezal:", int(com.sum()))
for c in ("Baixo", "Médio", "Alto", "Sem valor"):
    print(f"  {c:10s} {soma.get(c, 0):9.1f} km²")
print(f"  total      {total:9.1f} km²  | adequada (médio+alto): "
      f"{100 * (soma.get('Médio', 0) + soma.get('Alto', 0)) / total:.1f} %")
print("ok", SAIDA)
