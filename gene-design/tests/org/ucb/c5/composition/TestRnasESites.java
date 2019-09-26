package org.ucb.c5.composition;

/**
 * Todo:

    * Move RNAaseEChecker into SeqUtils, invoke it for this test, use in algorithms 
 
    * More clearly define the two violation conditions
 */



import java.util.ArrayList;
import java.util.List;
import org.junit.Test;
import org.junit.BeforeClass;
import static org.junit.Assert.*;
import org.ucb.c5.composition.model.Composition;
import org.ucb.c5.composition.model.Host;
import org.ucb.c5.composition.model.Transcript;
import org.ucb.c5.composition.model.Construct;

public class TestRnasESites {

    private static CompositionToDNA c2d;

    @BeforeClass
    public static void setUpClass() throws Exception {
        c2d = new CompositionToDNA();
        c2d.initiate();
    }

    @Test
    /**
     * RNAse E is an enzyme responsible for a lot of the mRNA degradation in
     * Ecoli. RNAse E acts on regions with a high content of A and U in a short
     * window. While there is some variation in the sequence, sites detected by
     * different species often share a string of 4-5 U's and an overall 80%+ AU
     * content within a 10 BP window. In checking DNA, we check AT
     * content/strings to see if there exists a match for the above traits.
     * http://www.jbc.org/content/272/1/609.full.pdf
     * http://www.jbc.org/content/269/14/10790.long
     *
     * @author Sahil Hansalia
     */
    public void RNaseEsitechecker() throws Exception {

        String promoter = "ttatgacaacttgacggctacatcattcactttttcttcacaa";
        String terminator = "TGCCTGGCGGCAGTAGCGCGGTGGTCCCACCTGACCCCATGCC";
        String protein = "MARCGGYCTMPRAKNKPCPLQTENKCEKPEFSWSICLYFCLILRSHKQELDDSWKPTSDDVFKAQWSIPGNTTSVPRWPFNMFRMLRAVLGSMFQVSKCDDVNNMLEIQYQHFNYEVYPAPYLKYNKIKRFYPYSLNWMQYMFMDVVYYADVCFCCKHVVMCCDIVQHIGTRMNNGYPVCRGMSRQNDHHAYTLFCVVMDSRPAKGCKC";
        ArrayList<String> proteins = new ArrayList<>();
        proteins.add(protein);
        Composition comp = new Composition(Host.Ecoli, promoter, proteins, terminator);
        Construct dna = c2d.run(comp);

        List<Transcript> outputDNA = dna.getmRNAs();
        String[] proteinsDnas = outputDNA.get(0).getCodons();
        StringBuilder out = new StringBuilder();
        for (String codon : proteinsDnas) {
            out.append(codon);
        }
        String sequence = out.toString();

        Boolean violatesT = false;
        Boolean violates80 = false;
        int window_size = 10;
        for (int i = 0; i < sequence.length() - window_size; i++) {
            String window = sequence.substring(i, i + window_size);
            for (int j = 0; j < window_size - 4; j++) {
                if (window.substring(j, j + 4).equals("TTTT")) {
                    violatesT = true;
                }
            }
            int counter = 0;
            for (int k = 0; k < window_size; k++) {
                if (window.substring(k).equals("A") | window.substring(k).equals("T")) {
                    counter++;
                }
            }
            if (counter >= 8) {
                violates80 = true;
            }
            assertTrue("both of the conditions are violated", !(violates80 && violatesT));
        }

    }
}
