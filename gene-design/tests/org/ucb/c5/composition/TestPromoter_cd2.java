package org.ucb.c5.composition;

import org.junit.BeforeClass;
import org.junit.Test;
import org.ucb.c5.composition.model.Composition;
import org.ucb.c5.composition.model.Construct;
import org.ucb.c5.composition.model.Host;
import org.ucb.c5.composition.model.Transcript;

import java.util.ArrayList;
import java.util.List;
import org.ucb.c5.composition.CompositionToDNA;
import org.ucb.c5.composition.checkers.PromoterChecker;

import static org.junit.Assert.assertTrue;

/**
 * @author Joanne Chang student joanne91218
 */
public class TestPromoter_cd2 {

    private static CompositionToDNA c2d;

    @BeforeClass
    public static void setUpClass() throws Exception {
        c2d = new CompositionToDNA();
        c2d.initiate();

    }

    /**
     * Checks if avoids forbidden sequences during long polymeric peptide runs
     *
     * @throws Exception
     */
    @Test
    public void avoidPromoter() throws Exception {
        PromoterChecker p = new PromoterChecker();
        p.initiate();

        String promoter = "TAGCGGATCCAGTTTATCCCTCGAAACTATTAGACGTACAGGTCGAAATCTTAAGTCAAAT";
        String terminator = "AGTCGTAGTCGATGCATTTACGATCGCTAGCCCTGTGTGAAAACGTAGCTGTGCTAGTCG";
        String peptide = "MKSGKWACTWACTATFMAT";
        List<String> proteins = new ArrayList<>();
        proteins.add(peptide);

        Composition comp = new Composition(Host.Ecoli, promoter, proteins, terminator);
        Construct dna = c2d.run(comp);
        List<Transcript> mRNAs = dna.getmRNAs();
        //System.out.println(mRNAs);
        String[] dnaCodons = mRNAs.get(0).getCodons();

        String dnaSeq = "";
        for (int i = 0; i < dnaCodons.length; i++) {
            dnaSeq = dnaSeq + dnaCodons[i];
        }
        assertTrue(p.run(dnaSeq));
    }

    @Test       //check JUnit test running true statement
    public void testRun() throws Exception {
        assertTrue(true);
    }

}
