# Material Suplementar

Este repositório é o material suplementar do Trabalho de Conclusão de Curso II em Oceanologia da Universidade Federal do Rio Grande (FURG), requisito para a obtenção do título de Bacharel em Oceanologia da aluna Ketlen Pirangi da Silva, sob orientação da Dra. Paula Spotorno de Oliveira e coorientação da Dra. Juliana Costi.

O trabalho é intitulado **"Análise do risco de introdução da ostra exótica invasora *Saccostrea cucullata* em áreas costeiras ligadas à indústria de óleo e gás na costa brasileira"**.

**O código explicado, módulo por módulo, está na página: https://ketlenps.github.io/risco-saccostrea-brasil-tcc/**

## Sobre o trabalho

A ostra *Saccostrea cucullata*, originária do Indo-Pacífico, foi introduzida na costa brasileira provavelmente pela água de lastro e pela bioincrustação em navios, e hoje já ocupa costões, manguezais e estruturas portuárias do Sudeste e do Sul. Neste trabalho, a adequabilidade ambiental da costa brasileira para a espécie foi estimada por um modelo de distribuição calibrado com registros de todo o mundo e quatro algoritmos combinados em um *ensemble* (biomod2). A adequabilidade foi combinada com a pressão exercida por portos, plataformas e poços de petróleo e gás para mapear o risco de introdução, que também foi avaliado nas ecorregiões marinhas, nos manguezais e em cenários de mudanças climáticas.

## Conteúdo

| Arquivo | O que é |
|---|---|
| `modelagem-saccostrea.Rmd` | Código do trabalho com a explicação de cada etapa (R Markdown) |
| `index.html` | A mesma coisa em página web, publicada no endereço acima |
| `Saccostrea_modelagem_TCC2.R` | Script completo em R, pronto para rodar |
| `scripts_python/` | Análises complementares feitas com o Python do QGIS (ecorregiões, manguezais e sensibilidade da escala de pressão) |
| `risco-saccostrea-brasil-tcc.Rproj` | Projeto do RStudio |

O trabalho usou o R 4.3.3 com o pacote biomod2 4.3-4-6. Os dados de entrada não estão no repositório. Para repetir a análise, coloque-os em `dados/`, nas pastas indicadas abaixo, abra o projeto no RStudio e rode o script. Os registros do GBIF e do OBIS, as camadas do Bio-ORACLE, os estados do IBGE e a linha de costa (GSHHG) são baixados pelo próprio script.

## Fontes dos dados

| Pasta | Dado | Fonte |
|---|---|---|
| `dados/ocorrencias/` | Registros da espécie | Literatura brasileira; GBIF (gbif.org); OBIS (obis.org) |
| `dados/ambientais/` | Variáveis ambientais, presente e cenários SSP | Bio-ORACLE v3 (Assis et al., 2024) |
| `dados/antropicas/` | Portos; plataformas em operação e poços marítimos | ANTAQ; ANP (dados abertos) |
| `dados/habitats/` | Manguezais | Atlas dos Manguezais do Brasil (ICMBio, 2018) |
| `dados/limites/` | Estados do Brasil e linha de costa | IBGE (pacote geobr); GSHHG 2.3.7 |

## Licença

O código está sob a licença MIT (veja [LICENSE](LICENSE)). Os dados de terceiros seguem as licenças das fontes originais.
