match_metabolite_names <- function(matchTo){
  
  to.match <- data.table(
    old = c("Aba", "acges", "Ala", "Arg", "Asn", "Asp", "C0", "C10", "C101", 
            "C12", "C14", "C141", "C14OH", "C16", "C161", "C161OH", "C16OH", 
            "C18", "C181", "C181OH", "C182", "C182OH", "C18OH", "C2", "C201", 
            "C202", "C203", "C3", "C3DC", "C4", "C4OH", "C5", "C51", "C5OHHMG", 
            "C6", "C6DC", "C8", "C81", "Carnosin", "Cit", "Gln", "Glu", "Glut", 
            "Gly", "His", "LeuIle", "Lys", "MeGlut", "MeHis", "Met", "MMA", 
            "OHProl", "Orn", "Phe", "PiPA", "Pro", "Sarc", "Ser", "Tau", 
            "Thr", "Trp", "Tyr", "Val"),
    new = c("Aba", "AC-total", "Ala", "Arg", "Asn", "Asp", "C0", "C10", 
            "C10:1", "C12", "C14", "C14:1", "C14OH", "C16", "C16:1", "C16:1OH", 
            "C16OH", "C18", "C18:1", "C18:1OH", "C18:2", "C18:2OH", "C18OH", 
            "C2", "C20:1", "C20:2", "C20:3", "C3", "C3DC", "C4", "C4OH", 
            "C5", "C5:1", "C5OH+HMG", "C6", "C6DC", "C8", "C8:1", "Carn", 
            "Cit", "Gln", "Glu", "Glut", "Gly", "His", "Leu|Ile", "Lys", 
            "MeGlut", "MeHis", "Met", "MMA", "OH-Prol", "Orn", "Phe", "PiPA", 
            "Pro", "Sarc", "Ser", "Tau", "Thr", "Trp", "Tyr", "Val"),
    long = c("Aminobutyric acid","Acylcarnitine total", "Alanine", "Arginine", 
             "L-Asparagine", "Aspartic acid", "Carnitine free", "Decanoylcarnitine", 
             "Decenoylcarnitine", "Dodecanoylcarnitine", "Myristoylcarnitine", 
             "Tetradecenoylcarnitine", "3-Hydroxy-tetradecanoylcarnitine", 
             "Palmitoylcarnitine", "Hexadecenoylcarnitine", "3-Hydroxy-hexadecenoylcarnitine", 
             "3-Hydroxy-hexadecanoylcarnitine", "Stearoylcarnitine", "Octadecenoylcarnitine", 
             "Hydroxy-octadec-1-enoylcarnitine", "Trans,trans-9,12-octadecadienoic acid (Linoelaidic)", 
             "3-Hydroxylinoleoylcarnitine", "3-Hydroxy-octadecanoylcarnitine", 
             "Acetylcarnitine", "Cis-11-eicosenoic acid", "Cis-11,14-eicosadienoic acid", 
             "Cis-11,14,17-eicosatrienoic acid", "Propionylcarnitine", "Malonylcarnitine", 
             "Butyrylcarnitine", "3-Hydroxy-butyryl-carnitine", "Isovalerylcarnitine", 
             "Tiglylcarnitine", "2-Hydroxyisovalerylcarnitine", "Hexanoylcarnitine", 
             "Adipoylcarnitine", "Octanoylcarnitine", "Octenoylcarnitine", 
             "Carnosine", "Citrulline", "Glutamine", "Glutamic acid", "Glutarylcarnitine", 
             "Glycine", "Histidine", "Leucine|Isoleucine", "Lysine", "Methylmalonylcarnitine", 
             "3-Methylglutarylcarnitine", "Methylhistidine", "Methionine", 
             "Hydroxyproline", "Ornithine", "Phenylalanine", "Pipecolic acid", 
             "Proline", "Sarcosine", "Serine", "Taurine", "Threonine", "Tryptophan", 
             "Tyrosine", "Valine")
  )
  
  matched.names <- match(matchTo, to.match$old)
  res <- to.match[matched.names, new]
  return(res)
}
