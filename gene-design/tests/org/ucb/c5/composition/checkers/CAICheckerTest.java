package org.ucb.c5.composition.checkers;

import org.junit.BeforeClass;
import org.junit.Test;
import org.ucb.c5.composition.CompositionToDNA;
import org.ucb.c5.composition.model.Composition;
import org.ucb.c5.composition.model.Construct;
import org.ucb.c5.composition.model.Host;
import org.ucb.c5.composition.model.Transcript;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * A class to test the Codon Adaptation Index
 * (CAI) for DNA sequences.
 * @author Stephen Lin
 */
public class CAICheckerTest {

    private static CompositionToDNA c2d;
    private static CAIChecker checker;

    @BeforeClass
    public static void setUpClass() throws Exception {
        c2d = new CompositionToDNA();
        c2d.initiate();

        checker = new CAIChecker();
        checker.initiate();
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
        assertEquals(1.0, checker.run(dna), 0);
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
        assertEquals(0.4807, checker.run(dna), 1e-4);
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
        checker.run(dna);
    }

    /**
     * Tests if the CAI is permissible
     * for a sequence designed via
     * TranscriptDesigner on the crtE protein.
     * Permissibility is based on a threshold of CAI = 0.5.
     * That is, the best codon is used at least 50% of the time.
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
        for (String codon : codons) {
            cds.append(codon);
        }
        System.out.println(checker.run(cds.toString()));
        assertTrue(checker.run(cds.toString()) >= 0.5);
    }

    /**
     * Tests if the CAI is permissible
     * for a sequence designed via
     * TranscriptDesigner on the PAIPDS protein.
     * Permissibility is based on a threshold of CAI = 0.5.
     * That is, the best codon is used at least 50% of the time.
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
        for (String codon : codons) {
            cds.append(codon);
        }
        System.out.println(checker.run(cds.toString()));
        assertTrue(checker.run(cds.toString()) >= 0.5);
    }

    /**
     * Tests if the CAI is permissible
     * for a sequence containing a
     * high frequency of arginine, which
     * has a relatively low CAI in the
     * reference set used.
     *
     * @author Stephen Lin
     */
    @Test
    public void testImpermissibleCAI() throws Exception {
        //Codon AGG encoding R has the lowest CAI (≈ 0.0032085561497326204)
        String promoter = "GTACCAGTTACGCAGTAGCATGCATGCTAGCTATGCGTAGTCGATCTAGCTGATCCTAGCTGCTG";
        String terminator = "AGTCGTAGTCGAAGCTTTTACGATCGCTAGCCCTGTGTGAAAACGTAGCTGTGCCAGTCG";
        List<String> proteins = new ArrayList<>();
        String polyRs = "MRRRHPRRRFLYRRRSRRRR";
        proteins.add(polyRs);
        Composition comp = new Composition(Host.Ecoli, promoter, proteins, terminator);
        Construct dna = c2d.run(comp);
        List<Transcript> mRNAs = dna.getmRNAs();
        String[] codons = mRNAs.get(0).getCodons();
        StringBuilder cds = new StringBuilder();
        for (String codon : codons) {
            cds.append(codon);
        }
        System.out.println(checker.run(cds.toString()));
        assertTrue(checker.run(cds.toString()) > 0.1);
    }
}

