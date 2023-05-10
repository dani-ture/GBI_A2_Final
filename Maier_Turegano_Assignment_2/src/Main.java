/**
 * Assignment 02
 * Authors: Annika Maier & Daniel Turegano
 */

/* *******************************************************************************
 *
 * This program was designed to be executed in a console/terminal
 *
 * In order to do that, the user should follow the instructions on the pdf "Maier_Turegano_Assignment_2.pdf"
 *
 *********************************************************************************/

import java.io.FileWriter;
import java.io.IOException;

public class Main {

    public static void main(String[] args) {

        try {

            //*********************************** INPUT *******************************

            int matchScore = Integer.parseInt(args[0]);
            int mismatchScore = Integer.parseInt(args[1]);
            int gapPenalty = Integer.parseInt(args[2]);
            String querySeqPath = args[3];
            String targetSeqPath = args[4];
            String outputDirectoryPath = args[5];

            /* Testing values for IntelliJ
            int matchScore = 3;
            int mismatchScore = -3;
            int gapPenalty = 5;
            String querySeqPath = "./resources/MTB_query.fasta"; ///////////////CHANGE
            String targetSeqPath = "./resources/MTB_target.fasta";
            String outputDirectoryPath = "./resources/";

             */

            System.out.println("****************** GBI Assignment 2 ******************");
            System.out.println("****************** Maier, Annika & Turégano, Daniel ******************\n");
            System.out.println("****************** EXERCISE 3 Smith-Waterman Algorithm ******************\n");

            // FastaReader to read FASTA files
            FastaReader fr = new FastaReader();

            // The corresponding sequences
            String targetSequence = fr.readInFasta(targetSeqPath).get(0).getSequence();
            String querySequence = fr.readInFasta(querySeqPath).get(0).getSequence();

            SmithWatermanAligner sma = new SmithWatermanAligner(matchScore, mismatchScore, gapPenalty, querySequence, targetSequence);

            System.out.println("These are your parameters:\n");
            System.out.println("Match score: " + sma.getMatchScore());
            System.out.println("Mismatch score: " + sma.getMismatchScore());
            System.out.println("GapPenalty: " + sma.getGapPenalty());
            System.out.println("Query sequence obtained from: " + querySeqPath);
            System.out.println("Target sequence obtained from: " + targetSeqPath);

            System.out.println("\nOptimal local alignment score: " + sma.alignSequences());
            System.out.println("It was reached " + sma.getOptimalIndices().size() + " time(s), in the following location(s) in the DP matrix:");

            for (Integer[] i : sma.getOptimalIndices()) {
                System.out.printf("(%d,%d)\n", i[0], i[1]);
            }

            System.out.println();

            //These are disabled because they look awful in the terminal
                //sma.printF();
                //sma.printT();

            //Output to file
            writeOutput(sma, querySeqPath, targetSeqPath, outputDirectoryPath);

        }
        catch (Exception e) {
            e.printStackTrace();
        }

    }

    /**
     * writeOutput() writes the output to a file in a similar way to the on displayed in exercise 4
     * @param sma - SmithWatermanAligner object
     * @param queryPath Path were the query sequence (.fasta file) is located
     * @param targetPath Path were the target sequence (.fasta file) is located
     * @param outputDirectoryPath Path to the directory were the user wants to save the output file
     */

    static void writeOutput(SmithWatermanAligner sma, String queryPath, String targetPath, String outputDirectoryPath) {

        try {

            FileWriter fw = new FileWriter(outputDirectoryPath+"Maier_Turegano_Assignment_2_Output.txt");

            fw.write("For the calculation of the optimal local alignment between the two sequences read from\n" +
                    queryPath + "\n" + targetPath +  "\n" +
            "the Smith-Waterman Algorithm was implemented using the following parameters:\n");
            fw.write("------------------------------------------------------------------------------------------------------------------------\n");
            fw.write("Parameters:\n");
            fw.write("s(a,b) if a == b (match): " + sma.getMatchScore() + "\n");
            fw.write("s(a,b) if a != b (mismatch): " + sma.getMismatchScore() + "\n");
            fw.write("d (gap penalty for linear scoring): " + sma.getGapPenalty() + "\n");
            fw.write("------------------------------------------------------------------------------------------------------------------------\n");

            fw.write("Optimal local alignment score: " + sma.alignSequences() + "\n");
            fw.write("(symbol key: ’+’ = match, ’-’ in symbol line = mismatch, ’-’ in sequence line = gap in sequence, ’ ’ = gap)" + "\n");

            fw.write("------------------------------------------------------------------------------------------------------------------------\n");

            int numberOfAlignments = 1;

            //For each cell where the max score was reached dwe write one possible alignment
            for(String[] alignment : sma.traceback()) {

                String yAligned = alignment[0];
                String casesAligned = alignment[1];
                String xAligned = alignment[2];
                String jOfAlignmentInTarget = Integer.toString(Integer.parseInt(alignment[3])+1);

                fw.write("The aligned sequences of alignment number " + numberOfAlignments + " has a length of "
                                + yAligned.length() + " characters. \n");
                fw.write("The alignment starts at the position number " + jOfAlignmentInTarget + " in the target sequence.\nThe alignment also has:\n");
                fw.write("------------------------------------------------------------------------------------------------------------------------\n");

                int[] symbols = calculateSymbolFrequency(alignment[1]);
                fw.write(symbols[0] + " matches,\n" + symbols[1] + " mismatches and\n" + symbols[2] + " gaps.\n");

                fw.write("------------------------------------------------------------------------------------------------------------------------\n");

                //Current alignment
                int p = 0;
                int l = 26; //Number of bases per line written to the file

                //System to write "l" bases per line
                while (p+l < yAligned.length()) {
                    //Position in the alignment
                    fw.write("Position:\t");
                    for (int i = p+1; i<=(p+l); i++) {
                        fw.write(i+"\t");
                    }
                    fw.write("\n");
                    fw.write("Target\t:\t");
                    fw.write(writeSubstringWithFormat(yAligned.substring(p, p + l)));
                    fw.write("Symbol\t:\t");
                    fw.write(writeSubstringWithFormat(casesAligned.substring(p, p + l)));
                    fw.write("Query\t:\t");
                    fw.write(writeSubstringWithFormat(xAligned.substring(p, p + l)));

                    p += l;
                }

                //This is needed to write the last line of every alignment
                fw.write("Position:\t");
                for (int i = p+1; i<=xAligned.length(); i++) {
                    fw.write((i)+"\t");
                }
                fw.write("\n");
                fw.write("Target\t:\t");
                fw.write(writeSubstringWithFormat(yAligned.substring(p)));
                fw.write("Symbol\t:\t");
                fw.write(writeSubstringWithFormat(casesAligned.substring(p)));
                fw.write("Query\t:\t");
                fw.write(writeSubstringWithFormat(xAligned.substring(p)));

                numberOfAlignments++;
                fw.write("------------------------------------------------------------------------------------------------------------------------\n");
            }

            System.out.println("Output successfully written to file");

            fw.close();

        }

        catch(IOException e){
            e.printStackTrace();
        }

    }

    /**
     * writeSubstringWithFormat()
     * @param s (String that will be formatted)
     * @return (formatted string)
     */
    static String writeSubstringWithFormat(String s){
            StringBuilder sb= new StringBuilder();
            for(String c : s.split("")) {
                sb.append(c + "\t");
            }
            sb.append("\n");
            return sb.toString();
    }

    /**
     * calculateSymbolFrequency()
     * @param s (String with symbols after the alignment)
     * @return int[] (contains the frequency of each symbol)
     */
    static int[] calculateSymbolFrequency(String s) {
        int[] symbolFrequencies = new int[3];

        int matchCount=0;
        int mismatchCount=0;
        int gapCount=0;

        for (String c : s.split("")) {
            switch (c) {
                case "+" -> matchCount++;
                case "-" -> mismatchCount++;
                case " " -> gapCount++;
            }
        }

        symbolFrequencies[0] = matchCount;
        symbolFrequencies[1] = mismatchCount;
        symbolFrequencies[2] = gapCount;

        return symbolFrequencies;

    }
}