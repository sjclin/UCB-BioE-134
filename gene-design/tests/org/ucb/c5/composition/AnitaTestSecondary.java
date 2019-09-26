package org.ucb.c5.composition;

import org.junit.BeforeClass;
import org.junit.Test;
import org.ucb.c5.composition.model.Composition;
import org.ucb.c5.composition.model.Construct;
import org.ucb.c5.composition.model.Host;
import org.ucb.c5.composition.model.Transcript;
import org.ucb.c5.sequtils.HairpinCounter;
import org.ucb.c5.utils.FileUtils;

import java.util.*;

import static org.junit.Assert.*;

/**
 * Project 3b: tests, specifically targeting secondary structure.
 * Closely based on the previous secondary structure test. 
 * But, there are a couple more proteins and the threshold is calculated, NOT arbitrary.
 * @author Anita Silver aphinneys
 */
public class AnitaTestSecondary {

    private static CompositionToDNA c2d;
    private static HairpinCounter hc;
    private static double meanSecondary;
    private static double meanSecondary37;

    @BeforeClass
    public static void setUpClass() throws Exception {
        c2d = new CompositionToDNA();
        c2d.initiate();

        hc = new HairpinCounter();
        hc.initiate();

        String data = FileUtils.readResourceFile("composition/data/coli_genes.txt");
                String[] entries = data.split("\\r|\\r?\\n");
        //Create and populate a list of 100 highest expressed coli gene strings to be run in hairpin counter.
        List<String> hundred_coli_genes = new ArrayList<>();
        for (int i = 0; i < 100; i++) {
            String entry = entries[i];
            String[] entry_array = entry.split("\t");
            String gene_seq = entry_array[6];
            hundred_coli_genes.add(gene_seq);
        }

        double secondarySum = 0;
        double secondary37Sum = 0;
        for (String gene : hundred_coli_genes) {
            double hairpinScore = hc.run(gene);
            double hairpinScore37 = hc.run(gene.substring(0, 37));
            secondarySum += hairpinScore;
            secondary37Sum += hairpinScore37;
        }
    double num = 100;
    meanSecondary = secondarySum / num;
    meanSecondary37 = secondary37Sum / num;
}

    /** Tests whether the output has more than twice as much secondary structure as the mean amount for
     * E Coli's top 100 expressed genes.
     * @throws Exception
     */
    @Test
    public void testSecondary() throws Exception {
        String PaIPDS = "MALLSSSLSSQIPTGSHPLTHTQCIPHFSTTINAGISAGKPRSFYLRWGKGSNKIIACVGEGTTSLPYQSAEKTDSLSAPTLVKREFPPGFWKDHVIDSLTSSHKVSAAEEKRMETLISEIKNIFRSMGYGETNPSAYDTAWVARIPAVDGSEHPEFPETLEWILQNQLKDGSWGEGFYFLAYDRILATLACIITLTLWRTGETQIRKGIEFFKTQAGKIEDEADSHRPSGFEIVFPAMLKEAKVLGLDLPYELPFIKQIIEKREAKLERLPTNILYALPTTLLYSLEGLQEIVDWEKIIKLQSKDGSFLTSPASTAAVFMRTGNKKCLEFLNFVLKKFGNHVPCHYPLDLFERLWAVDTVERLGIDHHFKEEIKDALDYVYSHWDERGIGWARENPIPDIDDTAMGLRILRLHGYNVSSDVLKTFRDENGEFFCFLGQTQRGVTDMLNVNRCSHVAFPGETIMQEAKLCTERYLRNALEDVGAFDKWALKKNIRGEVEYALKYPWHRSMPRLEARSYIEHYGPNDVWLGKTMYMMPYISNLKYLELAKLDFNHVQSLHQKELRDLRRWWKSSGLSELKFTRERVTEIYFSAASFIFEPEFATCRDVYTKISIFTVILDDLYDAHGTLDNLELFSEGVKRWDLSLVDRMPQDMKICFTVLYNTVNEIAVEGRKRQGRDVLGYIRNVLEILLAAHTKEAEWSAARYVPSFDEYIENASVSISLGTLVLISVLFTGEILTDDVLSKIGRGSRFLQLMGLTGRLVNDTKTYEAERGQGEVASAVQCYMKEHPEISEEEALKHVYTVMENALDELNREFVNNRDVPDSCRRLVFETARIMQLFYMEGDGLTLSHEMEIKEHVKNCLFQPVA";
        String crtE = "MYPFIRTARMTVCAKKHVHLTRDAAEQLLADIDRRLDQLLPVEGERDVVGAAMREGALAPGKRIRPMLLLLTARDLGCAVSHDGLLDLACAVEMVHAASLILDDMPCMDDAKLRRGRPTIHSHYGEHVAILAAVALLSKAFGVIADADGLTPLAKNRAVSELSNAIGMQGLVQGQFKDLSEGDKPRSAEAILMTNHFKTSTLFCASMQMASIVANASSEARDCLHRFSLDLGQAFQLLDDLTDGMTDTGKDSNQDAGKSTLVNLLGPRAVEERLRQHLQLASEHLSAACQHGHATQHFIQAWFDKKLAAVS";
        String peptide = "MYPFIRTARMTVCAKKHVHLTRDAAEQLLADIDRRLDQLLPVEGERDVVGAAMREGALAPGKRIRPMLLLLTARDLGCAVSHDGLLDLACAVEMVHAASLILDDMPCMDDAKLRRGRPTIHSHYGEHVAILAAVALLSKAFGVIADADGLTPLAKNRAVSELSNAIGMQGLVQGQFKDLSEGDKPRSAEAILMTNHFKTSTLFCASMQMASIVANASSEARDCLHRFSLDLGQAFQLLDDLTDGMTDTGKDSNQDAGKSTLVNLLGPRAVEERLRQHLQLASEHLSAACQHGHATQHFIQAWFDKKLAAVS";
        String fromPdf = "GLSDGEWQQVLNVWGK";
        String promoter = "ttatgacaacttgacggctacatcattcactttttcttcacaa";
        String terminator = "TGCCTGGCGGCAGTAGCGCGGTGGTCCCACCTGACCCCATGCC";
        ArrayList<String> list = new ArrayList<>();
        list.add(peptide);
        list.add(PaIPDS);
        list.add(crtE);
        list.add(fromPdf);
        Composition cp = new Composition(Host.Ecoli, promoter, list, terminator);
        Construct dna = c2d.run(cp);

        //TODO: make changes to this code?? (taken directly from past test)
        List<Transcript> protein2dna = dna.getmRNAs();
        String[] proteinsDnas = protein2dna.get(0).getCodons();
        StringBuilder check = new StringBuilder();
        for (String codon : proteinsDnas) {
            check.append(codon);
        }

        String seq = check.toString();
        double threshold = meanSecondary * 2;
        double score = hc.run(seq);
        assertTrue(score < threshold);

    }

    /** Tests whether the secondary structure levels in the the first 37 base pairs of the result
     * are no more than twice that of the first 37 base pairs of the top 100 genes.
     * Focusing on this region could be a more useful test than focusing on the whole thing.
     * @throws Exception
     */
    @Test
    public void testSecondary37() throws Exception {
        String PaIPDS = "MALLSSSLSSQIPTGSHPLTHTQCIPHFSTTINAGISAGKPRSFYLRWGKGSNKIIACVGEGTTSLPYQSAEKTDSLSAPTLVKREFPPGFWKDHVIDSLTSSHKVSAAEEKRMETLISEIKNIFRSMGYGETNPSAYDTAWVARIPAVDGSEHPEFPETLEWILQNQLKDGSWGEGFYFLAYDRILATLACIITLTLWRTGETQIRKGIEFFKTQAGKIEDEADSHRPSGFEIVFPAMLKEAKVLGLDLPYELPFIKQIIEKREAKLERLPTNILYALPTTLLYSLEGLQEIVDWEKIIKLQSKDGSFLTSPASTAAVFMRTGNKKCLEFLNFVLKKFGNHVPCHYPLDLFERLWAVDTVERLGIDHHFKEEIKDALDYVYSHWDERGIGWARENPIPDIDDTAMGLRILRLHGYNVSSDVLKTFRDENGEFFCFLGQTQRGVTDMLNVNRCSHVAFPGETIMQEAKLCTERYLRNALEDVGAFDKWALKKNIRGEVEYALKYPWHRSMPRLEARSYIEHYGPNDVWLGKTMYMMPYISNLKYLELAKLDFNHVQSLHQKELRDLRRWWKSSGLSELKFTRERVTEIYFSAASFIFEPEFATCRDVYTKISIFTVILDDLYDAHGTLDNLELFSEGVKRWDLSLVDRMPQDMKICFTVLYNTVNEIAVEGRKRQGRDVLGYIRNVLEILLAAHTKEAEWSAARYVPSFDEYIENASVSISLGTLVLISVLFTGEILTDDVLSKIGRGSRFLQLMGLTGRLVNDTKTYEAERGQGEVASAVQCYMKEHPEISEEEALKHVYTVMENALDELNREFVNNRDVPDSCRRLVFETARIMQLFYMEGDGLTLSHEMEIKEHVKNCLFQPVA";
        String crtE = "MYPFIRTARMTVCAKKHVHLTRDAAEQLLADIDRRLDQLLPVEGERDVVGAAMREGALAPGKRIRPMLLLLTARDLGCAVSHDGLLDLACAVEMVHAASLILDDMPCMDDAKLRRGRPTIHSHYGEHVAILAAVALLSKAFGVIADADGLTPLAKNRAVSELSNAIGMQGLVQGQFKDLSEGDKPRSAEAILMTNHFKTSTLFCASMQMASIVANASSEARDCLHRFSLDLGQAFQLLDDLTDGMTDTGKDSNQDAGKSTLVNLLGPRAVEERLRQHLQLASEHLSAACQHGHATQHFIQAWFDKKLAAVS";
        String peptide = "MYPFIRTARMTVCAKKHVHLTRDAAEQLLADIDRRLDQLLPVEGERDVVGAAMREGALAPGKRIRPMLLLLTARDLGCAVSHDGLLDLACAVEMVHAASLILDDMPCMDDAKLRRGRPTIHSHYGEHVAILAAVALLSKAFGVIADADGLTPLAKNRAVSELSNAIGMQGLVQGQFKDLSEGDKPRSAEAILMTNHFKTSTLFCASMQMASIVANASSEARDCLHRFSLDLGQAFQLLDDLTDGMTDTGKDSNQDAGKSTLVNLLGPRAVEERLRQHLQLASEHLSAACQHGHATQHFIQAWFDKKLAAVS";
        String fromPdf = "GLSDGEWQQVLNVWGK";
        String promoter = "ttatgacaacttgacggctacatcattcactttttcttcacaa";
        String terminator = "TGCCTGGCGGCAGTAGCGCGGTGGTCCCACCTGACCCCATGCC";
        ArrayList<String> list = new ArrayList<>();
        list.add(peptide);
        list.add(PaIPDS);
        list.add(crtE);
        list.add(fromPdf);
        Composition cp = new Composition(Host.Ecoli, promoter, list, terminator);
        Construct dna = c2d.run(cp);
        List<Transcript> protein2dna = dna.getmRNAs();
        String[] proteinsDnas = protein2dna.get(0).getCodons();
        StringBuilder check = new StringBuilder();
        for (String codon : proteinsDnas) {
            check.append(codon);
        }

        String seq = check.toString();
        double threshold = meanSecondary37 * 2;
        double score = hc.run(seq.substring(0, 37));
        assertTrue(score < threshold);
    }

}
 

