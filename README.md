# Reconstrução de Espectros de Absorção de Tecido Renal (Saudável vs. ChRCC)

Trabalho desenvolvido na UC de Métodos Óticos em Medicina, no 1ºano de Mestrado em Engenharia Biomédica, no ISEP. 

## Autores

Catarina Marques, Guilherme Pereira, Mariana Sá, Marta Seixas e Rita Portugal  
Departamento de Física, ISEP, Porto, Portugal

## Conceito

O espectro de absorção de um tecido biológico pode ser descrito como uma soma pesada dos espectros dos seus cromóforos (DNA, proteínas, HbO₂, Hb, água, lípidos, melanina e lipofuscina). Neste projeto, o coeficiente de absorção µa(λ) é calculado a partir de medições de transmitância e refletância total (200–1000 nm) e reconstruído por mínimos quadrados não negativos (NNLS, `lsqnonneg` do MATLAB), o que permite estimar a concentração de cada componente.

## Aplicação

Comparação de tecido renal saudável com carcinoma renal cromófobo (ChRCC), usando 10 amostras de cada tipo. O tecido tumoral apresenta mais DNA e proteínas, o que mostra o potencial da espectroscopia de absorção como técnica ótica, rápida e minimamente invasiva, de apoio ao diagnóstico de tumores renais.

## Diretórios

```
.
├── trabalho_5_n.m     # Script: tecido saudável
├── trabalho_5_c.m     # Script: tecido ChRCC
├── espectros_n/       # Dados do tecido saudável (Rt, Tt e espectros de referência)
├── espectros_c/       # Dados do tecido ChRCC (Rt, Tt e espectros de referência)
└── G5_MOTMBI.pdf      # Artigo
```

Requer MATLAB com Optimization Toolbox e a função `stdshade` (MATLAB File Exchange).

