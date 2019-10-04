package org.ucb.c5.composition;

import javafx.util.Pair;
import org.ucb.c5.composition.checkers.ForbiddenSequenceChecker;
import org.ucb.c5.composition.model.RBSOption;
import org.ucb.c5.composition.model.Transcript;
import org.ucb.c5.sequtils.Translate;
import org.ucb.c5.utils.FileUtils;

import java.util.*;

/**
 * This reverse translates a protein sequence to a DNA and chooses an RBS. It
 * uses the highest CAI codon for each amino acid in the specified Host.
 *
 * @author J. Christopher Anderson
 * @author Stephen Lin
 */
public class TranscriptDesigner {

    private RBSChooser2 rbsChooser;
    private Translate translator;
    private ForbiddenSequenceChecker forbiddenSeqChecker;
    private Map<Character, String[]> aaMap;
    private String aaCharacters;

    public void initiate() throws Exception {
        //Initialize the RBSChooser
        rbsChooser = new RBSChooser2();  //Use new algorithm to choose RBS
        rbsChooser.initiate();

        translator = new Translate();
        translator.initiate();

        forbiddenSeqChecker = new ForbiddenSequenceChecker();
        forbiddenSeqChecker.initiate();

        aaMap = new HashMap<>();
        Character[] aminoAcids = new Character[]{'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                                           'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};
        String[][] codons = new String[][]{{"GCA", "GCT", "GCC", "GCG"}, {"AGA", "AGG", "CGA", "CGT", "CGC", "CGG"},
            {"AAT", "AAC"}, {"GAT", "GAC"}, {"TGT", "TGC"}, {"CAA", "CAG"}, {"GAA", "GAG"}, {"GGA", "GGT", "GGC", "GGG"},
            {"CAT", "CAC"}, {"ATA", "ATT", "ATC"}, {"TTA", "TTG", "CTA", "CTT", "CTC", "CTG"}, {"AAA", "AAG"}, {"ATG"},
            {"TTT", "TTC"}, {"CCA", "CCT", "CCC", "CCG"}, {"AGT", "AGC", "TCA", "TCT", "TCC", "TCG"},
            {"ACA", "ACT", "ACC", "ACG"}, {"TGG"}, {"TAT", "TAC"}, {"GTA", "GTT", "GTC", "GTG"}};
        for (int i = 0; i < 20; i++) {
            aaMap.put(aminoAcids[i], codons[i]);
        }

        aaCharacters = "([ARNDCQEGHILKMFPSTWYV])+";
    }

    public Transcript run(String peptide, Set<RBSOption> ignores) throws Exception {
        //Choose codons for each amino acid to generate cds
        if (!peptide.matches(aaCharacters)) {
            throw new IllegalArgumentException("peptide contains non-amino acid character");
        }
        String[] codons = new String[peptide.length()];
        StringBuilder cds = new StringBuilder();
        for(int i = 0; i < peptide.length(); i++) {
            String aaWindow = peptide.substring(i, Math.min(i + 3, peptide.length()));
            String preamble;
            if (i >= 3) {
                preamble = cds.substring(i * 3 - 9, i * 3);
            } else {
                preamble = "";
            }
            List<String> synCodonsList = synonymousCodons(aaWindow);
            Comparator<String> compareByGC = Comparator.comparingDouble(this::gcContent);
            if (gcContent(cds.toString()) < 0.5) {
                synCodonsList.sort(compareByGC.reversed());
            } else {
                synCodonsList.sort(compareByGC);
            }
            for (String synCodons: synCodonsList) {
                String testSeq = preamble + synCodons;
                if (forbiddenSeqChecker.run(testSeq)) {
                    cds.append(synCodons, 0, 3);
                    codons[i] = synCodons.substring(0, 3);
                    break;
                }
            }
        }
        //Choose an rbs
        RBSOption selectedRBS = rbsChooser.run(cds.toString(), ignores);
        //Construct the Transcript and return it
        Transcript out = new Transcript(selectedRBS, peptide, codons);
        return out;
    }

    private double gcContent(String seq) {
        double gc = 0;
        for (int i = 0; i < seq.length(); i++) {
            char base = seq.charAt(i);
            if (base == 'G' || base == 'C') {
                gc++;
            }
        }
        return gc / seq.length();
    }

    private List<String> synonymousCodons(String aaWindow) {
        return synonymousCodons(aaWindow, new ArrayList<>());
    }

    /**
     * Given a String aaWindow of amino acids,
     * return all combinations of synonymous
     * codons that can possibly encode it in a List.
     * @param aaWindow Peptide sequence
     * @return synCodons
     */
    private List<String> synonymousCodons(String aaWindow, List<String> synCodons) {
        String[] currCodons = aaMap.get(aaWindow.charAt(0));
        if (aaWindow.length() == 1) {
            Collections.addAll(synCodons, currCodons);
            return synCodons;
        }
        synonymousCodons(aaWindow.substring(1), synCodons);
        String[] synCodonArr = synCodons.toArray(new String[synCodons.size()]);
        synCodons.clear();
        for (String synCodon: synCodonArr) {
            for (String validCodon: currCodons) {
                synCodons.add(validCodon + synCodon);
            }
        }
        return synCodons;
    }
}
