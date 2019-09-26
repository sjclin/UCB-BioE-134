package org.ucb.c5.composition.checkers;

import org.junit.BeforeClass;
import org.junit.Test;

public class InternalRBSCheckerTest {

    private static InternalRBSChecker checker;

    @BeforeClass
    public static void setUpClass() throws Exception {

        checker = new InternalRBSChecker();
        checker.initiate();
    }
    /**
     * Inputs known sequences and checks for Internal RBS sites. False if match is present.
     * Written by Naveen Kumaran
     * t's were used so no accidental sites would be created.
     * @throws Exception    /**
     * Inputs known sequences and checks for Internal RBS sites. False if match is present.
     * Written by Naveen Kumaran
     * t's were used so no accidental sites would be created.
     * @throws Exception
     */

    @Test
    public void NoInternalRBS() throws Exception {
        //Inputs sequence into checker with no internal RBS.
        // Should return positive
        String DNA = "tttttttttttttttttttttttttttATG";
        boolean result = checker.run(DNA);
        assert(result == true);
    }

    @Test
    public void YesInternalRBS() throws Exception {
        //Inputs sequence into checker with A KNOWN internal RBS.
        // Should return FALSE
        String DNA = "ttttttttttttttttttttttAGGAGGtttATGtttt";
        boolean result = checker.run(DNA);
        assert(result == false);
    }

    @Test
    public void InternalRBSAfterStartCodon() throws Exception {
        //Inputs sequence into checker with A KNOWN internal RBS after the start codon.
        // Should return TRUE
        String DNA = "tttttttttttttATGtttttAAGGAGGGttttttttttt";
        boolean result = checker.run(DNA);
        assert(result == true);
    }

    @Test
    public void InternalRBSBetweenStartCodons() throws Exception {
        //Inputs sequence into checker with A KNOWN internal RBS between 2 start codons.
        // Should return False
        String DNA = "tttttttttttttATGtttttAAGGAGGGtttATGtttttttt";
        boolean result = checker.run(DNA);
        assert(result == false);
    }

    @Test
    public void InternalRBSBeforeAndAfterStartCodon() throws Exception {
        //Inputs sequence into checker with A KNOWN internal RBS before and after a start codon.
        // Should return False
        // Checks to see if algorithm works in the presence of multiple RBS sites.
        String DNA = "ttttttAGGAGGtttttttATGtttttAAGGAGGGttttttttttt";
        boolean result = checker.run(DNA);
        assert(result == false);
    }

    @Test
    public void CloseByInternalRBSSite() throws Exception {
        //Inputs small sequence into checker with A KNOWN internal RBS immediately before start codon.
        // Should return TRUE b/c window starts 4 bp before start codon.
        String DNA = "ttttttttttttttttttAGGAATGttttttttttt";
        boolean result = checker.run(DNA);
        assert(result == true);
    }
}
