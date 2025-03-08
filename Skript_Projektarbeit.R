#Installationen
cwbtools::corpus_install(doi = "10.5281/zenodo.12794676")
install.packages("polmineR")
library(polmineR)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cwbtools)
library(topicmodels)
library(slam)
library(DT)
library(LDAvis)
library(tsne)
install.packages("quanteda.textmodels@0.9.1")
install.packages("topicmodels")
install.packages("ldatuning")
library("ldatuning")
library(quanteda)
use("polmineR")

# Eingrenzung der Daten

#Partitionierung in Sprecher*innen 19.Wahlperiode Bundestag
lp19 <- partition("GERMAPARL2", protocol_lp = 19, p_type = "speech")
lp19_speakers <- partition_bundle(lp19, s_attribute = "speaker_name", progress = TRUE)
#Anreicherung mit Wörtern und Festlegung auf Substantive
lp19_speakers <- enrich(lp19_speakers, p_attribute = c("word", "xpos"))
pb <- subset(lp19_speakers, xpos == "NN")
str(pb)
size(pb)
#nach Filterung, Ausschluss von xpos
pb@objects <- lapply(pb@objects, function(x){x@stat[, "xpos" := NULL]; x@p_attribute <- "word"; x})

#Erstellung einer Document-Term-Matrix

dtm <- as.DocumentTermMatrix(pb, col = "count")

#Bereinigung der DTM

#Ausschluss Dokumente unter 100 Wörter

short_docs <- which(slam::row_sums(dtm) < 100)
if (length(short_docs) > 0) dtm <- dtm[-short_docs,]

#Ausschluss seltener Wörter
rare_words <- which(slam::col_sums(dtm) < 5)
if (length(rare_words) > 0) dtm <- dtm[,-rare_words]

#Ausschluss leerer Dokumente
empty_docs <- which(slam::row_sums(dtm) == 0)
if (length(empty_docs) > 0) dtm <- dtm[-empty_docs,]

#Ausschluss von Stopwörtern
stopit <- tm::stopwords("de")
stopit_upper <- paste(toupper(substr(stopit, 1, 1)), substr(stopit, 2, nchar(stopit)), sep = "")
stopit_upper_where <- which(stopit_upper %in% colnames(dtm))
if (length(stopit_upper_where) > 0)dtm <- dtm[, -stopit_upper_where]

#Ausschließen 100 häufigste Wörter (da hier viele technische Begriffe Bundestag)
word_frequencies <- slam::col_sums(dtm)
sorted_words <- sort(word_frequencies, decreasing = TRUE)
top100_words <- names(sorted_words)[1:100]
print(top100_words)
dtm <- dtm[, !(colnames(dtm) %in% top100_words)]
#Erstellen und Anzeigen der Matrix
dtm_matrix <- as.matrix(dtm)
View(dtm_matrix)
#speicherung dtm
saveRDS(dtm, "\\Users\\LENOVO\\OneDrive\\Desktop\\GermaParl-Arbeit\\dtm.rds")
#Berechnung lda

#Berechnung welche Themensetzung ideal; leider hier unendliche Ladezeit; erst mal weggelassen
result <- FindTopicsNumber(
  dtm,
  topics = seq(from = 5, to = 6, by = 1),
  metrics = c("Griffiths2004", "CaoJuan2009", "Arun2010", "Deveaud2014"),
  method = "Gibbs",
  control = list(seed = 77),
  mc.cores = 2L,
  verbose = TRUE
)
#Berechnung 20 Themen mithilfe lda
lda <- topicmodels::LDA(
  dtm, k = 20, method = "Gibbs",
  control = list(burnin = 200, iter = 500, keep = 50, verbose = 20, alpha = NULL)
)
#Anzeige 20 Wörterzu den Themen
topics_terms <- terms(lda, 20)
print(topics_terms)

#erstellen Wahrscheinlichkeiten Themen aus lda Modell und Ausgabe

tmResult <- posterior(lda)
theta <- tmResult$topics  # Verteilung Topic und Dokumente
beta <- tmResult$terms
top5termsPerTopic <- terms(lda, 5)
topicNames <- apply(top5termsPerTopic, 2, paste, collapse = " ")
print(topicNames)

#Partei BSW erstellen, da damals noch nicht existent

# generieren Liste der späteren BSW-Abgeordneten
later_BSW_abgeordnete <- c("Sevim Dağdelen", "Klaus Ernst", "Andrej Hunko", "Fabio De Masi", "Amira Mohamed Ali", "Zaklin Nastic", 
                           "Heike Hänsel", "Jessica Tatti", "Alexander Ulrich", "Sahra Wagenknecht", "Sabine Zimmermann", "Friedrich Straetmanns")

# Durch alle Sprechernamen gehen und Partei zuweisen
pb@objects <- lapply(pb@objects, function(x) {
  speaker_names <- s_attributes(x, "speaker_name")  # Sprecher-Namen
  speaker_party <- s_attributes(x, "speaker_party") # Parteizugehörigkeit
  
  #Vergleich der Namen und Zuweisung neuer Partei
  speaker_party[speaker_names %in% later_BSW_abgeordnete] <- "later_BSW"
  
  # speichern neuer Zuweisung
  x@stat$speaker_party <- speaker_party
  return(x)
})
parties <- sapply(pb@objects, function(x) {
  unique(x@stat$speaker_party)
})
names(parties) <- rownames(dtm)
##Überprüfung ob alle zugeordnet

zugeordnete_abgeordnete <- unique(unlist(sapply(pb@objects, function(x) {
  if ("later_BSW" %in% x@stat$speaker_party) {
    return(s_attributes(x, "speaker_name"))
  }
})))

# Zeige die zugeordneten Abgeordneten an
print(zugeordnete_abgeordnete)

#Da Dokumente  gefiltert, haben sie weniger Spalten als speaker_party, darum nur die Dokumente, in dem beides vorhanden
common_docs <- intersect(rownames(dtm), names(parties))
parties_filtered <- parties[common_docs]
theta_filtered <- theta[common_docs, ]

#Ausschluss von Topics, die nicht aussagekräftig genug für einen Vergleich sind, Ausgabe und Überführung in Datatable
unwanted_topics <- c(6,7,10,12,14,16)
theta_cleaned <- theta_filtered[, -unwanted_topics]
beta_cleaned <- beta[-unwanted_topics, ]
top5termsPerTopic_cleaned <- terms(lda, 5)[, -unwanted_topics]
topicNames_cleaned <- apply(top5termsPerTopic_cleaned, 2, paste, collapse = " ")
print(topicNames_cleaned)
topicNames_cleaned_df <- as.data.frame(topicNames_cleaned)
datatable(topicNames_cleaned_df)
print(topicNames_cleaned)


# DataFrame der gefilterten Daten
party_topic_distributions <- data.frame(
  speaker_party = unlist(parties_filtered),
  as.data.frame(theta_cleaned),
  stringsAsFactors = FALSE
)
#Topic-Verteilung pro partei berechnen

party_topic_means <- aggregate(
  . ~ speaker_party,
  data = party_topic_distributions,
  FUN = mean
)
#Änderung Fomrat für Visualisierung
party_topic_means_long <- pivot_longer(
  party_topic_means,
  cols = -speaker_party,
  names_to = "topic",
  values_to = "mean_probability"
)
#Ausschluss Fraktionslose für bessere Übersichtlichkeit Balkendiagramm
party_topic_means_long <- party_topic_means_long[!(party_topic_means_long$speaker_party %in% c("Die PARTEI", "NA", "parteilos", "LKR")), ]
#Anzeige Partei-Topic-Beziehungen
print(party_topic_means)
#Darstellung als Tabelle
datatable(party_topic_means)


#Visualisierung als Balkendiagramm
ggplot(party_topic_means_long, aes(x = topic, y = mean_probability, fill = speaker_party)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ speaker_party) +
  theme_minimal() +
  labs(
    title = "Durchschnitt Topic  Parteien",
    x = "Topic",
    y = "Durchschnitt Wahrscheinlichkeit"
  )

#Berechnung der Ähnlichkeit von Topics
library("proxy")
aehnlichkeit_matrix <- as.matrix(simil(party_topic_means[, -1], method = "cosine"))
rownames(aehnlichkeit_matrix) <- party_topic_means$speaker_party
print(aehnlichkeit_matrix)
datatable(aehnlichkeit_matrix)
colnames(aehnlichkeit_matrix) <- party_topic_means$speaker_party
heatmap(as.matrix(aehnlichkeit_matrix), main = "Themenähnlichkeit der Parteien")
#Ähnlichkeit BSW
bsw_aehnlichkeit <- data.frame(
  party = colnames(aehnlichkeit_matrix),
  similarity = aehnlichkeit_matrix["later_BSW", ]
)

#Filterung Parteien
bsw_aehnlichkeit <- bsw_aehnlichkeit[!bsw_aehnlichkeit$party %in% c("Die PARTEI", "NA", "parteilos", "LKR"), ]
bsw_aehnlichkeit <- bsw_aehnlichkeit[order(-bsw_aehnlichkeit$similarity), ]
bsw_aehnlichkeit$party <- factor(bsw_aehnlichkeit$party, levels = bsw_aehnlichkeit$party)

#Visualisierung Balkendiagramm Ähnlichkeit BSW
ggplot(bsw_aehnlichkeit, aes(x = party, y = similarity)) +
  geom_col() +
  labs(title = "Ähnlichkeit von 'later_BSW' mit anderen Parteien",
       x = "Partei",
       y = "Kosinus-Ähnlichkeit")

