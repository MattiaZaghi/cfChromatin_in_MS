library(tabulapdf)                        # parses tables in PDFs
pdf   <- "cfChromatin_in_MS/single_cell_data_processing/41587_2022_1250_MOESM1_ESM.pdf"
tab   <- extract_tables(pdf, pages = 33:37)
df    <- do.call(rbind, lapply(tab, as.data.frame))
colnames(df) <- c("DNA_ID","Description","Clone","Barcode")

adt <- df[, c("DNA_ID","Barcode")]       # 173 antibody tags
write.table(adt,
            file = "resources/adt.tsv",
            sep  = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

hto <- data.frame(
  ID = paste0("Hashtag_", 1:8),
  BC = c("GAACGACTACAGTCC","TGTCAGGATACCTTC","CAGTCATAGAGCTCT",
         "ACTGGTGTAGCAAAC","AGTTGCCGTGGATTG","TCGAGTCCAAGGCTA",
         "TCCAATTTCAGTGAT","CTGCAAGTACGACTG"))
write.table(hto,
            file = "resources/hto.tsv",
            sep  = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)
