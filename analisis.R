---
title: "Árbol Filogenético de Variantes de COVID-19 en 20 países"
author: "Ana Paula Figueroa - A00835792"
date: "2026-02-06"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Configuración Inicial y Librerías
```{r}
library(Biostrings) 
library(seqinr)  
library(adegenet)  
library(ape)  
library(ggtree)  
library(DECIPHER)  
library(viridis)  
library(ggplot2) 
library("seqinr") 
```

## Descarga de secuencias
```{r}
# Accesiones GenBank de variantes SARS-CoV-2
virus_ids <- c("OQ871050","OM918221","ON434752","OX448653","ON583472",
               "LC654482","OL966992","MW491232","OP777347","OQ059020",
               "MT478019","MW976780","OP748488","MZ410618","ON005388",
               "MW633894","OK356626","MZ314345","OK547173","MW273788")

# Descargar secuencias desde GenBank
virus_sequences <- read.GenBank(virus_ids)

# Asignar nombres más interpretables (países)
names(virus_sequences) <- c("EU","India","Francia","Alemania","Brasil",
                            "Japón","Corea_del_Sur","Italia","Reino_Unido",
                            "Rusia","Turquía","España","Vietnam",
                            "Australia","Taiwán","Argentina",
                            "Países_Bajos","Irán","México","Polonia")

# Guardar en formato FASTA
write.dna(virus_sequences, file = "virus_seqs.fasta",
          format = "fasta", nbcol = 6, colsep = " ", colw = 10) 


```


## Lectura y alineamiento
```{r}

# Leer FASTA como objeto Biostrings
virus_raw <- readDNAStringSet("virus_seqs.fasta")

# Orientar automáticamente las secuencias
virus_oriented <- OrientNucleotides(virus_raw)

# Alineamiento múltiple
virus_aligned <- AlignSeqs(virus_oriented)

# Guardar alineamiento
writeXStringSet(virus_aligned, file = "virus_aligned.fasta")

# Convertir a formato compatible con ape
virus_alignment <- read.alignment("virus_aligned.fasta", format = "fasta")

```


## Longitud de secuencias
```{r}

# Función para calcular longitud real de una secuencia
calc_length <- function(seq) {
  seq <- gsub("\n| ", "", seq)
  nchar(seq)
}

# Leer FASTA con seqinr para inspección directa
fasta_seq <- read.fasta("virus_seqs.fasta", as.string = TRUE)

cat("Longitud de cada secuencia:\n")
for (i in seq_along(fasta_seq)) {
  cat(names(fasta_seq)[i], "-", calc_length(fasta_seq[[i]]), "\n")
}

```

## Composicion de bases
```{r}

# Función que calcula porcentajes A, G, T, C
base_percentages <- function(seq) {
  seq <- tolower(gsub("\n| ", "", seq))
  bases <- strsplit(seq, NULL)[[1]]
  round(100 * table(factor(bases, levels = c("a","g","t","c"))) / length(bases), 2)
}

# Calcular para todas las secuencias
base_matrix <- sapply(fasta_seq, base_percentages)

# Convertir a data frame largo para ggplot
base_df <- as.data.frame(t(base_matrix))
base_df$Pais <- rownames(base_df)
base_df_long <- reshape2::melt(base_df, id.vars = "Pais",
                               variable.name = "Base",
                               value.name = "Porcentaje")

```


## Grafica de composicion de bases
```{r}

ggplot(base_df_long, aes(x = Pais, y = Porcentaje, fill = Base)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d() +
  theme_minimal() +
  labs(title = "Composición de Bases de Variantes SARS-CoV-2 por País",
       x = "País",
       y = "Porcentaje (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

```

## Matriz de distancias y arbol filogenetica
```{r}

# Calcular matriz de distancia (basada en similitud)
dist_matrix <- dist.alignment(virus_alignment, matrix = "identity")

# Construcción del árbol usando Neighbor Joining
virus_tree <- nj(dist_matrix)

# Ordenar visualmente el árbol
virus_tree <- ladderize(virus_tree)

```


## Visualizacion del arbol
```{r}

plot(virus_tree, cex = 0.8)
title("Árbol filogenético de SARS-CoV-2 en 20 países")

```

