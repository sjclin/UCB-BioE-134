package org.ucb.c5.composition.checkers;

import org.ucb.c5.composition.CompositionToDNA;
import org.ucb.c5.composition.model.Composition;
import org.ucb.c5.composition.model.Construct;
import org.ucb.c5.composition.model.Host;
import org.ucb.c5.composition.model.Transcript;
import org.ucb.c5.utils.FileUtils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * A class to check the Codon Adaptation Index (CAI)
 * for DNA sequences. The CAI is defined as the geometric
 * mean of the individual relative adaptiveness parameters for
 * the amino acids in a peptide of a translated DNA sequence.
 * The codon usage table is computed from a subset of the data in
 * coli_genes.txt based on the GENE_LIMIT.
 * @author Stephen Lin
 */
public class CAIChecker {

    private static Map<String, Double> relAdaptParam;

    public void initiate() throws Exception {
        String coli_genes = FileUtils.readResourceFile("composition/data/coli_genes.txt");
        String[] lines = coli_genes.split("\\r|\\r?\\n");
        List<String> topEColiGenes = new ArrayList<>();
        //Parse the top GENE_LIMIT expressing genes in coli_genes.txt
        int GENE_LIMIT = 700;
        for (int i = 0; i < GENE_LIMIT; i++) {
            String line = lines[i];
            String[] values = line.split("\t");
            String seq = values[6];
            topEColiGenes.add(seq);
        }
        String[][] allCodons = new String[][]{{"GCA", "GCT", "GCC", "GCG"}, {"AGA", "AGG", "CGA", "CGT", "CGC", "CGG"},
                {"AAT", "AAC"}, {"GAT", "GAC"}, {"TGT", "TGC"}, {"CAA", "CAG"}, {"GAA", "GAG"}, {"GGA", "GGT", "GGC", "GGG"},
                {"CAT", "CAC"}, {"ATA", "ATT", "ATC"}, {"TTA", "TTG", "CTA", "CTT", "CTC", "CTG"}, {"AAA", "AAG"}, {"ATG"},
                {"TTT", "TTC"}, {"CCA", "CCT", "CCC", "CCG"}, {"AGT", "AGC", "TCA", "TCT", "TCC", "TCG"},
                {"ACA", "ACT", "ACC", "ACG"}, {"TGG"}, {"TAT", "TAC"}, {"GTA", "GTT", "GTC", "GTG"}};
        Map<String, Integer> codonCounts = new HashMap<>();
        //Count occurrences of codons in topEColiGenes
        for (String seq: topEColiGenes) {
            for (int i = 0; i < seq.length(); i += 3) {
                String codon = seq.substring(i, i + 3);
                if (codonCounts.containsKey(codon)) {
                    codonCounts.put(codon, codonCounts.get(codon) + 1);
                } else {
                    codonCounts.put(codon, 1);
                }
            }
        }
        //Compute the relative adaptive parameters for each codon for quick retrieval
        relAdaptParam = new HashMap<>();
        for (String[] codons: allCodons) {
            int maxCount = 0;
            for (String codon: codons) {
                int count = codonCounts.get(codon);
                if (count > maxCount) {
                    maxCount = count;
                }
            }
            for (String codon: codons) {
                relAdaptParam.put(codon, (double) codonCounts.get(codon) / maxCount);
            }
        }
    }

    /**
     * Given a DNA sequence, calculate the CAI
     * based on relAdaptParam (the relative
     * adaptiveness parameter) for each codon
     * in the sequence.
     *
     * @param dna A DNA sequence encoding a peptide
     * @return CAI of dna
     * @author Stephen Lin
     */

    public double run(String dna) throws IllegalArgumentException {
        if (dna.isEmpty()) {
            throw new IllegalArgumentException("CAI not defined for empty DNA sequence");
        }
        if (dna.length() % 3 != 0) {
            throw new IllegalArgumentException("DNA sequence length must be multiple of 3");
        }
        double geoMean = 1;
        double N = dna.length() / 3;
        for (int i = 0; i < dna.length(); i += 3) {
            String codon = dna.substring(i, i + 3);
            geoMean *= Math.pow(relAdaptParam.get(codon), 1/N);
        }
        return geoMean;
    }
}
