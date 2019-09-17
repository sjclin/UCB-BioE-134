package org.ucb.c5.composition;

import java.util.*;

import org.ucb.c5.composition.model.RBSOption;
import org.ucb.c5.sequtils.CalcEditDistance;
import org.ucb.c5.sequtils.HairpinCounter;
import org.ucb.c5.sequtils.Translate;
import org.ucb.c5.utils.FileUtils;

/**
 * Second generation RBSChooser algorithm
 *
 * Employs a list of genes and their associated ribosome binding sites for
 * highly-expressed proteins in E. coli.
 *
 * @author J. Christopher Anderson
 */
public class RBSChooser2 {

    private List<RBSOption> rbss;
    private Translate translator;

    public void initiate() throws Exception {
        rbss = new ArrayList<>();
        translator = new Translate();
        translator.initiate();

        String coli_genes = FileUtils.readFile("src/org/ucb/c5/composition/data/coli_genes.txt");
        String[] lines = coli_genes.split("\\r|\\r?\\n");
        Map<String, String> nameToCds = new HashMap<>();
        for (String line: lines) {
            String[] values = line.split("\t");
            String name = values[1];
            String cds = values[6];
            nameToCds.put(name, cds);
        }

        String rbs_options = FileUtils.readFile("src/org/ucb/c5/composition/data/rbs_options.txt");
        lines = rbs_options.split("\\r|\\r?\\n");
        for (String line: lines) {
            String[] values = line.split("\t");
            String name =  values[0];
            String rbs = values[1];
            String cds = nameToCds.get(name);
            String first6aas = getFirst6aas(cds, translator);
            rbss.add(new RBSOption(name, null, rbs, cds, first6aas));
        }
    }

    private String getFirst6aas(String cds, Translate translator) {
        StringBuilder aas = new StringBuilder();
        for (int i = 0; i < 6; i++) {
            String aa = translator.run(cds.substring(3 * i, 3 * i + 3));
            aas.append(aa);
        }
        return aas.toString();
    }

    /**
     * Provided an ORF of sequence 'cds', this computes the best ribosome
     * binding site to use from a list of options.
     * 
     * It also permits a list of options to exclude.
     *
     * @param cds The DNA sequence, ie ATGCATGAT...
     * @param ignores The list of RBS's to exclude
     * @return The optimal ribosome binding site
     * @throws Exception when cds is an invalid DNA sequence
     */
    public RBSOption run(String cds, Set<RBSOption> ignores) throws Exception {
        if (!cds.matches("([ATCG])+"))
            throw new Exception();
        Map<Integer, List<RBSOption>> RbsScores = new HashMap<>();
        CalcEditDistance calc = new CalcEditDistance();
        calc.initiate();
        for (RBSOption rbs: rbss) {
            if (ignores.contains(rbs)) {
                continue;
            }
            String cdsFirst6aas = getFirst6aas(cds, translator);
            String rbsFirst6aas = rbs.getFirst6aas();
            int editDist = calc.run(cdsFirst6aas, rbsFirst6aas);
            if (!RbsScores.containsKey(editDist)) {
                List<RBSOption> rbsList = new ArrayList<>();
                rbsList.add(rbs);
                RbsScores.put(editDist, rbsList);
            } else {
                RbsScores.get(editDist).add(rbs);
            }
        }
        List<RBSOption> bestMatchRbss = new ArrayList<>();
        for (int i = 0; i < 7; i++) {
            if (RbsScores.containsKey(i)) {
                bestMatchRbss = RbsScores.get(i);
                break;
            }
        }
        double leastHairpins = Double.POSITIVE_INFINITY;
        RBSOption bestRbs = null;
        HairpinCounter counter = new HairpinCounter();
        counter.initiate();
        for (RBSOption rbs: bestMatchRbss) {
            double hairpins = counter.run(rbs.getRbs());
            if (hairpins < leastHairpins) {
                leastHairpins = hairpins;
                bestRbs = rbs;
            }
        }
        return bestRbs;
    }


    public static void main(String[] args) throws Exception {
        //Create an example
        String cds = "ATGGTAAGAAAACAGTTGCAGAGAGTTGAATTATCACCATCGTTATATGACACAGCTTGGGTGGCTATGGTGCCGGAGCGTAGTTCTTCTCAA";

        //Initiate the chooser
        RBSChooser2 chooser = new RBSChooser2();
        chooser.initiate();

        //Make the first choice with an empty Set of ignores
        Set<RBSOption> ignores = new HashSet<>();
        RBSOption selected1 = chooser.run(cds, ignores);

        //Add the first selection to the list of things to ignore
        ignores.add(selected1);

        //Choose again with an ignore added
        RBSOption selected2 = chooser.run(cds, ignores);

        //Print out the two options, which should be different
        System.out.println("CDS starts with:");
        System.out.println(cds.substring(0, 18));
        System.out.println();
        System.out.println("Selected1:\n");
        System.out.println(selected1.toString());
        System.out.println();
        System.out.println("Selected2:\n");
        System.out.println(selected2.toString());
    }
}
