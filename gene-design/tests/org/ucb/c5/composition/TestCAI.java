package org.ucb.c5.composition;

import org.junit.BeforeClass;
import org.junit.Test;
import org.ucb.c5.composition.model.Composition;
import org.ucb.c5.composition.model.Construct;
import org.ucb.c5.composition.model.Host;
import org.ucb.c5.composition.model.Transcript;
import org.ucb.c5.utils.FileUtils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.*;

/**
 * A class to calculate and test the Codon Adaptation Index (CAI)
 * for DNA sequences. The CAI is defined as the geometric
 * mean of the individual relative adaptiveness parameters for
 * the amino acids in a peptide of a translated DNA sequence.
 * The codon usage table is computed from a subset of the data in
 * coli_genes.txt based on the GENE_LIMIT.
 * @author Stephen Lin
 */
public class TestCAI {

    private static CompositionToDNA c2d;
    private static Map<String, Double> relAdaptParam;

    @BeforeClass
    public static void setUpClass() throws Exception {
        c2d = new CompositionToDNA();
        c2d.initiate();

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
     * Tests if the CAI for a trivial
     * sequence is correct.
     *
     * @author Stephen Lin
     */
    @Test
    public void testTrivialCAI() {
        //ATG and TGG have relative adaptiveness parameters equal to 1
        String dna = "ATGTGG";
        assertEquals(1.0, CAI(dna), 0);
    }

    /**
     * Tests if the CAI for a non-trivial
     * sequence (with different relative
     * adaptiveness parameters) is correct.
     *
     * @author Stephen Lin
     */
    @Test
    public void testNontrivialCAI() {
        //Codons CTT, GAC have relative adaptiveness parameters of ~0.1612 and ~0.6892 respectively
        String dna = "ATGCTTGAC";
        assertEquals(0.4807, CAI(dna), 1e-4);
    }

    /**
     * Tests if the CAI is well-defined,
     * i.e. the DNA sequence provided can
     * be fully translated to a peptide,
     * from which the CAI can be calculated.
     *
     * @author Stephen Lin
     */
    @Test(expected = IllegalArgumentException.class)
    public void testIllegalCAI() {
        //Translation of complete sequence is impossible
        String dna = "ATGAAGTT";
        CAI(dna);
    }

    /**
     * Tests if the CAI is permissible
     * for a sequence designed via
     * TranscriptDesigner on the crtE protein.
     * Permissibility is based on a threshold of CAI = 0.3.
     *
     * @author Stephen Lin
     */
    @Test
    public void testPermissible1CAI() throws Exception {
        String promoter = "GTACCAGTTACGCAGTAGCATGCATGCTAGCTATGCGTAGTCGATCTAGCTGATCCTAGCTGCTG";
        String terminator = "AGTCGTAGTCGAAGCTTTTACGATCGCTAGCCCTGTGTGAAAACGTAGCTGTGCCAGTCG";
        List<String> proteins = new ArrayList<>();
        String crtE = "MYPFIRTARMTVCAKKHVHLTRDAAEQLLADIDRRLDQLLPVEGERDVVGAAMREGALAPGKRIRPMLLLLTARDLGCAVSHDGLLDLACAVEMVHAASLILDDMPCMDDAKLRRGRPTIHSHYGEHVAILAAVALLSKAFGVIADADGLTPLAKNRAVSELSNAIGMQGLVQGQFKDLSEGDKPRSAEAILMTNHFKTSTLFCASMQMASIVANASSEARDCLHRFSLDLGQAFQLLDDLTDGMTDTGKDSNQDAGKSTLVNLLGPRAVEERLRQHLQLASEHLSAACQHGHATQHFIQAWFDKKLAAVS";
        proteins.add(crtE);
        Composition comp = new Composition(Host.Ecoli, promoter, proteins, terminator);
        Construct dna = c2d.run(comp);
        List<Transcript> mRNAs = dna.getmRNAs();
        String[] codons = mRNAs.get(0).getCodons();
        StringBuilder cds = new StringBuilder();
        for (String codon: codons) {
            cds.append(codon);
        }
        assertTrue(CAI(cds.toString()) > 0.3);
    }

    /**
     * Tests if the CAI is permissible
     * for a sequence designed via
     * TranscriptDesigner on the PAIDPS protein.
     * Permissibility is based on a threshold of CAI = 0.8.
     *
     * @author Stephen Lin
     */
    @Test
    public void testPermissible2CAI() throws Exception {
        String promoter = "GTACCAGTTACGCAGTAGCATGCATGCTAGCTATGCGTAGTCGATCTAGCTGATCCTAGCTGCTG";
        String terminator = "AGTCGTAGTCGAAGCTTTTACGATCGCTAGCCCTGTGTGAAAACGTAGCTGTGCCAGTCG";
        List<String> proteins = new ArrayList<>();
        String PaIPDS = "MALLSSSLSSQIPTGSHPLTHTQCIPHFSTTINAGISAGKPRSFYLRWGKGSNKIIACVGEGTTSLPYQSAEKTDSLSAPTLVKREFPPGFWKDHVIDSLTSSHKVSAAEEKRMETLISEIKNIFRSMGYGETNPSAYDTAWVARIPAVDGSEHPEFPETLEWILQNQLKDGSWGEGFYFLAYDRILATLACIITLTLWRTGETQIRKGIEFFKTQAGKIEDEADSHRPSGFEIVFPAMLKEAKVLGLDLPYELPFIKQIIEKREAKLERLPTNILYALPTTLLYSLEGLQEIVDWEKIIKLQSKDGSFLTSPASTAAVFMRTGNKKCLEFLNFVLKKFGNHVPCHYPLDLFERLWAVDTVERLGIDHHFKEEIKDALDYVYSHWDERGIGWARENPIPDIDDTAMGLRILRLHGYNVSSDVLKTFRDENGEFFCFLGQTQRGVTDMLNVNRCSHVAFPGETIMQEAKLCTERYLRNALEDVGAFDKWALKKNIRGEVEYALKYPWHRSMPRLEARSYIEHYGPNDVWLGKTMYMMPYISNLKYLELAKLDFNHVQSLHQKELRDLRRWWKSSGLSELKFTRERVTEIYFSAASFIFEPEFATCRDVYTKISIFTVILDDLYDAHGTLDNLELFSEGVKRWDLSLVDRMPQDMKICFTVLYNTVNEIAVEGRKRQGRDVLGYIRNVLEILLAAHTKEAEWSAARYVPSFDEYIENASVSISLGTLVLISVLFTGEILTDDVLSKIGRGSRFLQLMGLTGRLVNDTKTYEAERGQGEVASAVQCYMKEHPEISEEEALKHVYTVMENALDELNREFVNNRDVPDSCRRLVFETARIMQLFYMEGDGLTLSHEMEIKEHVKNCLFQPVA";
        proteins.add(PaIPDS);
        Composition comp = new Composition(Host.Ecoli, promoter, proteins, terminator);
        Construct dna = c2d.run(comp);
        List<Transcript> mRNAs = dna.getmRNAs();
        String[] codons = mRNAs.get(0).getCodons();
        StringBuilder cds = new StringBuilder();
        for (String codon: codons) {
            cds.append(codon);
        }
        assertTrue(CAI(cds.toString()) > 0.3);
    }

    /** 
     * Tests if the CAI is impermissible
     * for a sequence with minimal CAI.
     *
     * @author Stephen Lin
     */
    @Test
    public void testImpermissibleCAI() {
        //Codon AGG has the lowest CAI (≈ 0.0032085561497326204); a naive algorithm for TranscriptDesigner could output this DNA sequence
        String dna = "ATGAGG";
        System.out.println(CAI(dna));
        assertFalse(CAI(dna) > 0.3);
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

    private double CAI(String dna) throws IllegalArgumentException {
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
