import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.TreeMap;

/**
 * Assignment 02
 * Authors: Annika Maier & Daniel Turegano
 */

public class SmithWatermanAligner {

    /**
     * Integer value used to score matches during DP recurrence.
     */
    private final int matchScore;

    /**
     * Integer value used to score mismatches during DP recurrence.
     */
    private final int mismatchScore;

    /**
     * Integer value used to penalize gaps during DP recurrence.
     */
    private final int gapPenalty;

    /**
     * First sequence to align.
     */
    private final String s1;

    /**
     * Second sequence to align.
     */
    private final String s2;

    // We will use a 2D-Array for the Distance Matrix (F)
    private int[][] F;

    // and a 3D-Array for the Traceback matrix (T)
    private int[][][] T;

    // List with the indices pairs (i,j) which indicate the position of the cells which contain the max score
    // (as it can be achieved multiple times)
    private ArrayList<Integer[]> optimalIndices;

    //Constructor
    public SmithWatermanAligner(int matchScore, int mismatchScore, int gapPenalty, String s1, String s2) {
        this.matchScore = matchScore;
        this.mismatchScore = mismatchScore;
        this.gapPenalty = gapPenalty;
        this.s1 = s1;
        this.s2 = s2;
    }

    /**
     * Simple getters
     */

    public int getMatchScore() {
        return matchScore;
    }

    public int getMismatchScore() {
        return mismatchScore;
    }

    public int getGapPenalty() {
        return gapPenalty;
    }

    public String getS1() {
        return s1;
    }

    public String getS2() {
        return s2;
    }

    public int[][] getF() {
        return F;
    }

    public int[][][] getT() {
        return T;
    }
    public ArrayList<Integer[]> getOptimalIndices() {
        return optimalIndices;
    }

    /**
     * alignSequences() builds the DP matrix (F) following the Smith-Waterman algorithm, the traceback matrix (T)
     * @return Integer (score of the alignment)
     */

    public Integer alignSequences() {
        int optimalAlignmentScore = 0;

        //X = first sequence ; Y = second sequence
        String X = s1;
        String Y = s2;

        //n = length of the X sequence ; m = length of the Y sequence
        int n = X.length();
        int m = Y.length();

        //Scoring (F) and Traceback (T) matrices initialization
        //They have size (n+1)x(m+1) because of the initialization of the first row and column
        this.F = new int[n + 1][m + 1];
        //The third dimension of the T matrix will store the i value in k=0 and the j value ein k=1
        this.T = new int[n + 1][m + 1][2];

        //Initialize first cell
        F[0][0] = 0;
        T[0][0][0] = 0;
        T[0][0][1] = 0;

        //Initialize top row:
        for (int j = 1; j <= m; j++) {
            F[0][j] = 0;

            T[0][j][0] = 0;
            T[0][j][1] = j-1;
        }

        //Initialize left-most column:
        for (int i = 1; i <= n; i++) {
            F[i][0] = 0;

            T[i][0][0] = i-1;
            T[i][0][1] = 0;
        }

        //4 cases:

        //Rows: 4 possible cases (Match(M)/Replacement(R), Insertion(I) , Deletion(D), Zero(Z))
        //Columns: in the first one (cases[i][0]) store the possible score
        //         in the second one (cases[i][1]) store the i index of the origin cell
        //          in the third one (cases[i][2]) store the j index of the origin cell
        int[][] cases = new int[4][3];

        //Here we'll place the highest score in matrix F
        int max = 0;

        //List of pairs of indices where the optimal score is located
        this.optimalIndices = new ArrayList<>();

        //Fill matrices
        for (int i = 1; i <= n; i++) {                  //For each row
            for (int j = 1; j <= m; j++) {             //Move across the columns

                //Case 1 (mis)match
                cases[0][0] = F[i - 1][j - 1];
                cases[0][1] = i - 1;
                cases[0][2] = j - 1;

                //s(a,b)
                //Compare the character (i-1) from sequence X
                //to the character (j-1) from sequence Y, because the first
                //chars of both Strings are located at X[0] and Y[0]
                if (X.charAt(i - 1) == Y.charAt(j - 1)) {   //if xi == yj; add match score
                    cases[0][0] += matchScore;
                } else {  //if xi != yj; add negative misMatch score
                    cases[0][0] += mismatchScore;
                }

                //Case 2 (deletion)
                cases[1][0] = F[i - 1][j] - gapPenalty;
                cases[1][1] = i - 1;
                cases[1][2] = j;

                //Case 3 (insertion)
                cases[2][0] = F[i][j - 1] - gapPenalty;
                cases[2][1] = i;
                cases[2][2] = j - 1;

                //Case 4 (0)
                cases[3][0] = 0;
                cases[3][1] = 0;
                cases[3][2] = 0;

                //Sort the four rows of cases by the value of the first cell in each row (the score),
                // the maximum will be placed at the first cell of the last row

                Arrays.sort(cases, new Comparator<int[]>() {
                    public int compare(int[] a, int[] b) {
                        return Integer.compare(a[0], b[0]);
                    }
                });

                //Place it in the scoring matrix
                F[i][j] = cases[3][0];

                //Traceback
                T[i][j][0] = cases[3][1];
                T[i][j][1] = cases[3][2];

                //Check if it's the highest score (it can be reached multiple times)
                if(F[i][j] > max){
                    max = F[i][j];
                    optimalIndices.clear(); //New max scored reached, delete previous indices
                    optimalIndices.add(new Integer[]{i, j});
                }
                else if(F[i][j] == max){
                    optimalIndices.add(new Integer[]{i, j});
                }

            }

        }
        optimalAlignmentScore = max;
        return optimalAlignmentScore;
    }

    /**
     * traceback() reconstructs one possible alignment from each cell where the max score was achieved
     * @return ArrayList<String[]> , which is a list with all the reconstructed alignments. Each reconstructed
     * alignment has the aligned target subsequence (Y'), a string of symbols with info about the aligned bases,
     * the aligned query subsequence (X') and the index j in the target where the alignment starts.
     */
    public ArrayList<String[]> traceback() {

        ArrayList<String[]> allAlignments = new ArrayList<>();

            //for each pair of indices (i,j) where the max score was achieved
            for (Integer[] indices : optimalIndices) {

                // Array of size 4 for: Target, Symbol, Query, Index in the target where the alignment starts
                String[] alignment = new String[4];

                StringBuilder xAlignedReversed = new StringBuilder();
                StringBuilder yAlignedReversed = new StringBuilder();
                StringBuilder casesStringReversed = new StringBuilder();

                // Position of the hightest score
                int i = indices[0];
                int j = indices[1];

                // temporary variable to exchange values
                int temp;

                // index of alignment start in target
                int jOfAlignmentInTarget = indices[1];

                // while F is not 0
                while (F[i][j] != 0) {

                    // case 1 (mis)match
                    // T[i][j][0] = Row Position of where the value was computed from
                    // T[i][j][1] = Column position of where the value was computed from
                    if (T[i][j][0] == i - 1 && T[i][j][1] == j - 1) {
                        xAlignedReversed.append(s1.charAt(i - 1));
                        yAlignedReversed.append(s2.charAt(j - 1));

                        jOfAlignmentInTarget--;
                        // Checking if we have a match or mismactch
                        if (s1.charAt(i - 1) == s2.charAt(j - 1)) {   //if xi == yj; add match score
                            casesStringReversed.append("+");
                        } else {
                            casesStringReversed.append("-");
                        }
                    // case 2 - deletion
                    } else if (T[i][j][0] == i - 1 && T[i][j][1] == j) {
                        xAlignedReversed.append(s1.charAt(i - 1));
                        yAlignedReversed.append("-");
                        casesStringReversed.append(" ");

                    // case 3 - insertion
                    } else {
                        xAlignedReversed.append("-");
                        yAlignedReversed.append(s2.charAt(j - 1));
                        casesStringReversed.append(" ");
                        jOfAlignmentInTarget--;
                    }

                    // Update indices
                    temp=i;
                    i = T[i][j][0];
                    j = T[temp][j][1];

                }

                //Reverse al StringBuilders so that they keep the correct order
                alignment[0] =  yAlignedReversed.reverse().toString();     //Target seq
                alignment[1] =  casesStringReversed.reverse().toString();  //Cases
                alignment[2] =  xAlignedReversed.reverse().toString();     //Query seq

                alignment[3] =  Integer.toString(jOfAlignmentInTarget);    //Track target position

                allAlignments.add(alignment);

            }
            return allAlignments;

        }


    /**
     * printF() prints the whole F matrix
     */
    public void printF(){

        System.out.println("This is your DP Matrix:");
        //Print first row (Y sequence)
        System.out.print("\t\t");
        for (int i = 0; i < s2.length(); i++){
            System.out.print(s2.charAt(i) + "\t");
        }
        System.out.println();

        //Print line of "_"
        System.out.print("\t");
        for (int i = 0; i <= s2.length(); i++){
            System.out.print("_" + "\t");
        }
        System.out.println();

        //Print rest of the matrix
        int index = 0;
        for(int[] row : F){
            System.out.print((" " + s1).charAt(index) + " | ");
            for(int cell : row){
                System.out.print(cell + "\t");
            }
            System.out.println();
            index++;
        }
        System.out.println();
    }

    /**
     * printT() prints the whole Traceback matrix
     */
    public void printT(){

        System.out.println("This is your Traceback Matrix");
        //Print first row (Y sequence)
        System.out.print("\t\t\t\t");
        for (int i = 0; i < s2.length(); i++){
            System.out.print(s2.charAt(i) + "\t\t");
        }
        System.out.println();

        //Print line of "_"
        System.out.print("\t");
        for (int i = 0; i <= s2.length(); i++){
            System.out.print(" ___" + "\t");
        }
        System.out.println();

        //Print rest of the matrix
        int index = 0;
        for(int i = 0; i<=s1.length(); i++) {
            System.out.print((" " + s1).charAt(index) + " | ");
            for(int j = 0; j<=s2.length(); j++) {

                System.out.printf( "(%d,%d)\t", T[i][j][0],T[i][j][1]);
            }
            index++;
            System.out.println();
        }

    }
}