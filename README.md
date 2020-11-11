# seqBio: An R package

El paquete seqBio contiene la funcion seq_alignment(), util para estudiar la significacion estadistica de un alineamiento por parejas.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisitos

El paquete necesita una instalación báscia del software R. Además, será necesaria también la instalación de varios paquetes: 

```
install.packages("devtools")
install.packages("BiocManager")
install.packages("seqinr")
install.packages("reliaR")
install.packages("gplots")
```

### Instalación

Para instalar este paquete, utiliza la siguiente linia de codigo:

```
devtools::install_github("laurajuliamelis/Bioinformatics_Task3", subdir="seqBio", build_vignettes=TRUE)
```

Luego, necesitarás cargar el paquete:

```
library(seqBio)
```


### Ayuda 

Para obtener más información acerca del paquete:
```
?seqBio
```

Para buscar en la ayuda de la función seq_alignment(), utilizar el siguiente código:
```
?seq_alignment
```
En la ayuda se podrá encontrar una breve descripción del paquete, los argumentos y sus definiciones, los valores de salida, las referencias, y varios ejemplos.  

Para obtener más información de la función generateSeqsWithMultinomialModel(), visitar la ayuda utilizando el siguiente código:

```
?generateSeqsWithMultinomialModel
```
Esta es una función necesaria para la ejecucion de seq_alignment().

### Vignette (informe dinámico)

Adicionalmente, se puede acceder al informe dinámico del paquete para obtener una explicación más detallada y extensiva de la función seq_align().

Para acceder a esta viñeta, se debe utilizar el siguiente código:

```
browseVignettes('seqBio')
```

A continuación, pulsar el boton "**HTML**".
s

## Licencia

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

