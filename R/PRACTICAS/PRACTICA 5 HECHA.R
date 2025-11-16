head#############################################################################
#
# PRACTICA R
#
# Expresión diferencial de genes de ratón
# Microarray de Affymetrix (Affymetrix Murine Genome U74A version 2 MG_U74Av2
# Origen de los datos: GEO GSE5583 (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5583)
# Publicación: Mol Cell Biol 2006 Nov;26(21):7913-28.  16940178 (http://www.ncbi.nlm.nih.gov/pubmed/16940178)
#
# Muestras: 3 Wild Type x 3 Histone deacetylase 1 (HDAC1)
#
# R código original (credits): Ahmed Moustafa
#
#
##############################################################################

# Instalar RCurl

if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("RCurl")

# Si esto falla, que seguro lo hace tratar de instalarlo usando el menú, Paquetes, Servidor Spain A Coruña, RCurl

# Cargamos el paquete y los datos
library(RCurl)
url = getURL ("http://bit.ly/GSE5583_data", followlocation = TRUE)
data = as.matrix(read.table (text = url, row.names = 1, header = T))

# Chequeamos las dimensiones de los datos, y vemos las primeras y las últimas filas
dim(data)
head(data)
tail(data)

# Hacemos un primer histograma para explorar los datos
hist(data, col = "red", main="GSE5583 - Histogram")

# Transformamos los datos con un logaritmo 
# ¿Qué pasa si hacemos una transformación logarítima de los datos? ¿Para qué sirve?
data2 = log2(data)

#la expreion genica suele tener una distribucion muy asimétrica. el logaritmo trasforma esos datos para que parezcan a una distribucion noramal, lo que mejora los analisis estadísticos y las comparaciones

hist(data2, col = "pink", main="GSE5583 (log2) - Histogram")
#cambia de forma el grafico, a una distribucion normal
#sirve para mejorar la visualizacion de los datos


# Hacemos un boxplot con los datos transformados. ¿Qué significan los parámetros que hemos empleado?
# ¿Qué es un boxplot? 
#es una representacion grafica que resume la distribucion estadistica de los datos
#la caja representa el rango intercuartilico
#la linea dentro es la mediana
#los bigotes muestran la variabilidad fuera de los cuartiles
#los puntos son valores atipicos

boxplot(data2, col=c("blue", "blue", "blue",
	"orange", "orange", "orange"),
	main="GSE5583 - boxplots", las=2)#las=2 sirve para poner los nombres de las columnas en vertical
#las=2 es para rotar las etiquetas en el eje de las x en vertical

#la caja coge el rango intercuarticulo, los bigotes son las desviaciones, los puntos qeu salen son los puntos con valores extremos

# Hacemos un hierarchical clustering de las muestras basándonos en un coeficiente de correlación ç
#clustering jerarquico (agrupa por similitud de expresion)
# de los valores de expresión. ¿Es correcta la separación? #si, es correcta

#si, la separacion es correcta, las muestras WT y los KO tienden a agruparse
#ademass, las muestras pueden quedar mal agrupadas si tienen valores atipicos

hc = hclust(as.dist(1-cor(data2)))
plot(hc, main="GSE5583 - Hierarchical Clustering")#nos agrupa por perfiles parecidos por ejemplo de genes parecidas)#por ejemplo aqui separa las de WT y las de knockout, habiendo una muestra que se entremezcla y eso es un perfil raro


#######################################
# Análisis de Expresión Diferencial 
#######################################
head(data)
# Primero separamos las dos condiciones. ¿Qué tipo de datos has generado?
wt <- data[,1:3]
ko <- data[,4:6]
class(wt)

# son matrices con las intesidades de expresion por gen y por replica

head(wt)
head(ko)

#es un listado de genes que tienen expresiones en replica de wt y knockout, la idea es detectar los genes que tienen una diferencia significativa entre uno y otro, indicando que tienen un papel y que cambia la expresion.
# Calcula las medias de las muestras para cada condición. Usa apply
wt.mean = apply(wt, 1, mean)
ko.mean = apply(ko, 1, mean)
head(wt.mean)#sale cada gen y su media de los sies primeros genes
head(ko.mean)
#al hacer resta y ver el simbolo que nos da sabremos si es mayor la wt y la knockout

# ¿Cuál es la media más alta? #37460.5
limit = max(wt.mean, ko.mean)#tambien se puede sacar el minimo
limit

# Ahora hacemos un scatter plot (gráfico de dispersión)
plot(ko.mean ~ wt.mean, xlab = "WT", ylab = "KO",#significa que asocie unas con otras, x lab es para poner el texto del eje de las x y el y lo mismo. y el lim es para poner el valor limite
	main = "GSE5583 - Scatter", xlim = c(0, limit), ylim = c(0, limit))
# Añadir una línea diagonal con abline
abline(0, 1, col = "red")

#es un grafico de dispersion en la que se comparan los datos de x con los de y
#en el eje x los del wt y en el y los de knockout

# ¿Eres capaz de añadirle un grid?#si, ademas al eliminar el comentario en las tres ultimas aparecen las lineas en el grafico de dispersion
grid()
#abline(a, b): línea de pendiente b y ordenada en el origen a
#abline(h=y): línea horizontal
#abline(v=x): línea vertical
abline(1, 2, col = "red")     # línea y = 2x + 1
abline(h = 2, col = "green")  # línea y = 2
abline(v = 3, col = "violet") # línea x = 3

# Calculamos la diferencia entre las medias de las condiciones
diff.mean = wt.mean - ko.mean

# Hacemos un histograma de las diferencias de medias
hist(diff.mean, col = "gray")

# Calculamos la significancia estadística con un t-test.#PARA CADA GEN, USAMOS LOS DATOS SIN TRASNFORMAR PORQUE T-TEST ASUME NOMRLAIDAD DE DATOS ORIGINALES
# Primero crea una lista vacía para guardar los p-values
# Segundo crea una lista vacía para guardar las estadísticas del test.
# OJO que aquí usamos los datos SIN TRANSFORMAR. ¿Por qué?
# ¿Cuántas valores tiene cada muestra?
#cada muestra tiene 3 valores

pvalue = NULL 
tstat = NULL 
for(i in 1 : nrow(data)) { #Para cada gen
	x = wt[i,] # gene wt número i
	y = ko[i,] # gene ko número i
	
	# Hacemos el test
	t = t.test(x, y)
	
	# Añadimos el p-value a la lista
	pvalue[i] = t$p.value
	# Añadimos las estadísticas a la lista
	tstat[i] = t$statistic
}

head(pvalue)

# Ahora comprobamos que hemos hecho TODOS los cálculos
length(pvalue)

# Hacemos un histograma de los p-values.
# ¿Qué pasa si le ponemos con una transformación de -log10?
#si usamos -log10 (pvalue), los valores mas peuqeños (significativos) se veran mas arriba
hist(pvalue,col="gray")
hist(-log10(pvalue), col = "gray")

# Hacemos un volcano plot. Aquí podemos meter la diferencia de medias y la significancia estadística
plot(diff.mean, -log10(pvalue), main = "GSE5583 - Volcano")

# Queremos establecer que el mínimo para considerar una diferencia significativa, es con una diferencia de 2 y un p-value de 0.01
# ¿Puedes representarlo en el gráfico?
diff.mean_cutoff = 2
pvalue_cutoff = 0.01
abline(v = diff.mean_cutoff, col = "blue", lwd = 3)
#abline(v = -diff.mean_cutoff, col = "red", lwd = 3)
abline(h = -log10(pvalue_cutoff), col = "green", lwd = 3)

# Ahora buscamos los genes que satisfagan estos criterios
# Primero hacemos el filtro para la diferencia de medias (fold)
filter_by_diff.mean = abs(diff.mean) >= diff.mean_cutoff
dim(data[filter_by_diff.mean, ])

# Ahora el filtro de p-value
filter_by_pvalue = pvalue <= pvalue_cutoff
dim(data[filter_by_pvalue, ])

# Ahora las combinamos. ¿Cuántos genes cumplen los dos criterios? lo cumplen 426 genes 
filter_combined = filter_by_diff.mean & filter_by_pvalue
filtered = data[filter_combined,]
dim(filtered)
head(filtered)

# Ahora generamos otro volcano plot con los genes seleccionados marcados en rojo
plot(diff.mean, -log10(pvalue), main = "GSE5583 - Volcano #2")
points (diff.mean[filter_combined], -log10(pvalue[filter_combined]),col = "red")

# Ahora vamos a marcar los que estarían sobreexpresados (rojo) y reprimidos (azul). ¿Por qué parece que están al revés?

#porque el signo de diff.mean depende del orden Wt-KO, si invertimos el orden, se invierte el color
plot(diff.mean, -log10(pvalue), main = "GSE5583 - Volcano #3")
points (diff.mean[filter_combined & diff.mean < 0],
	-log10(pvalue[filter_combined & diff.mean < 0]), col = "red")
points (diff.mean[filter_combined & diff.mean > 0],
	-log10(pvalue[filter_combined & diff.mean > 0]), col = "blue")

# el eje x es la diferencia de medias, es decir cuanto de mucho se asemejan el wt del knockout. por tanto al estar cerca de cero significa que se parecen. y eje de las y al expresar el p valor y este es significativo a partir de 0.05.
#el eje y significancia estadistica (-log10pvalue)
#cuanto mas a la derecha o izquiera mas diferencia de expresion
#cuanto mas arriba mas significativo

# Ahora vamos a generar un mapa. Para ello primero tenemos que hacer un cluster de las columnas y los genes 
# ¿Qué es cada parámetro que hemos usado dentro de la función heatmap? 

#Rowv=dendrograma de genes (filas)
#colv=dendrograma de muestras (columnas)
#cexcol=tamaño de texto de columnas
#labRow=FALSE= oculta nombres de genes

# ¿Eres capaz de cambiar los colores del heatmap? Pista: usar el argumento col y hcl.colors
rowv = as.dendrogram(hclust(as.dist(1-cor(t(filtered)))))
colv = as.dendrogram(hclust(as.dist(1-cor(filtered))))
heatmap(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,labRow=FALSE)

heatmap(filtered)


# Ahora vamos a crear un heatmap más chulo. Para ello necesitamos dos paquetes: gplots y RcolorBrewer
#if (!requireNamespace("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install(c("gplots","RColorBrewer"))
install.packages("gplots")		
install.packages("RColorBrewer")	

library(gplots)
library(RColorBrewer)

# Hacemos nuestro heatmap
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,
	col = rev(redblue(256)), scale = "row")

# Lo guardamos en un archivo PDF
pdf ("GSE5583_DE_Heatmap.pdf")
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,
	col = rev(redblue(256)), scale = "row",labRow=FALSE)
dev.off()
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,col = redgreen(75), scale = "row",labRow=FALSE)

# Guardamos los genes diferencialmente expresados y filtrados en un fichero
write.table (filtered, "GSE5583_DE.txt", sep = "\t",quote = FALSE)
