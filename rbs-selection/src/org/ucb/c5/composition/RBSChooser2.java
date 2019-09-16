package org.ucb.c5.composition;

import java.util.*;

import org.ucb.c5.composition.model.RBSOption;
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
    //TODO:  Fill in

    public void initiate() throws Exception {
        rbss = new ArrayList<>();
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
        Translate translator = new Translate();
        translator.initiate();
        for (String line: lines) {
            String[] values = line.split("\t");
            String name =  values[0];
            String rbs = values[1];
            String cds = nameToCds.get(name);
            StringBuilder first6aas = new StringBuilder();
            for (int j = 0; j < 6; j++) {
                String aa = translator.run(cds.substring(3 * j, 3 * j + 3));
                first6aas.append(aa);
            }
            rbss.add(new RBSOption(name, null, rbs, cds, first6aas.toString()));
        }
        System.out.println("hello world");
        //TODO:  Fill in
    }

    /**
     * Provided an ORF of sequence 'cds', this computes the best ribosome
     * binding site to use from a list of options.
     * 
     * It also permits a list of options to exclude.
     *
     * @param cds The DNA sequence, ie ATGCATGAT...
     * @param ignores The list of RBS's to exclude
     * @return
     * @throws Exception
     */
    public RBSOption run(String cds, Set<RBSOption> ignores) throws Exception {
        //TODO:  Fill in
        
        return null;
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
