/**
 * 
 */
package lib.structure.capr;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import gui.aptatrace.logo.Logo;

/**
 * @author Jan Hoinka
 * Standalone command line interface for CapR
 */
public class CapRCLI {

	private static CapR capr = new CapR();
	
	private static String dir = System.getProperty("user.dir");
	
	/**
	 * Main
	 * @param args list of RNA or DNA strings. Will produce 
	 * results named sequence1, sequence2, ... sequenceN in order
	 * of apperence in args
	 */
	public static void main(String[] args) {
		
		// We need to be able to create graphical instances in pure console mode
		System.setProperty("java.awt.headless", "true");

		int sequence_counter = 1;
		
		for (String sequence : args) {
			
			System.out.println("Processing Sequence " + sequence_counter);
			
			// Make sure we have a valid sequence
			String seq = validateAlphabet(sequence);
			
			// Compute profile
			capr.ComputeStructuralProfile(seq.getBytes(), seq.length());

			// Get profile
			double[][] profile =  getMatrix(capr.getStructuralProfile(), sequence.length());
			
			// Output profile
			saveTxt(sequence_counter, String.format("sequence%s_profile.txt", sequence_counter), sequence, profile);
			saveLogo(sequence_counter, String.format("sequence%s_logo.pdf", sequence_counter), sequence, profile);
			
			sequence_counter++;
		}
		
		System.out.println("Prediction completed. Exiting.");

	}

	
	/**
	 * Beaks if alphabet is invalid, converts Us to Ts 
	 */
	private static String validateAlphabet(String seq) {
		
		//convert to upper
		seq = seq.toUpperCase();
		
		//validate
		StringBuilder sb = new StringBuilder();
		for ( byte c : seq.getBytes()) {
			
			//valid?
			if (c != (byte)'A' && c != (byte)'C' && c != (byte)'G' && c != (byte)'T' && c != (byte)'U') {
				
				throw new RuntimeException(String.format("ERROR: Sequence %s contains invalid character %s", seq, (char)c));
				
			}
			
			//convert (if required) and append
			sb.append( c == (byte)'U' ? 'T' : (char)c);
			
		}
		
		return sb.toString();
		
	}
	
	private static void saveTxt(int number, String filename, String sequence, double[][] probabilities) {
		
		Path path = Paths.get(dir , filename);
		
		//Use try-with-resource to get auto-closeable writer instance
		try (BufferedWriter writer = Files.newBufferedWriter(path))
		{
			// Store sequence ID
		    writer.write(String.format(">Sequence%s\n", number));

		    // Store nucleotides, tab separated
		    writer.write("\t");
		    writer.write(String.join("\t", sequence.split("")));
		    writer.write("\n");
		    
		    // Write probability for each context 
		    String[] hairpin  = new String[sequence.length()];
		    String[] inner    = new String[sequence.length()];
		    String[] bulge    = new String[sequence.length()];
		    String[] multi    = new String[sequence.length()];
		    String[] dangling = new String[sequence.length()];
		    String[] paired   = new String[sequence.length()];
		    
		    for (int pos=0; pos<probabilities.length; pos++) {
		    	
		    	hairpin[pos] = String.format("%.5f", probabilities[pos][0]);
			    inner[pos] = String.format("%.5f", probabilities[pos][1]);
			    bulge[pos] = String.format("%.5f", probabilities[pos][2]);
			    multi[pos] = String.format("%.5f", probabilities[pos][3]);
			    dangling[pos] = String.format("%.5f", probabilities[pos][4]);
			    paired[pos] = String.format("%.5f", probabilities[pos][5]);
		    	
		    }
		    
		    writer.write("H\t");
		    writer.write(String.join("\t", hairpin));
		    writer.write("\n");
		    
		    writer.write("I\t");
		    writer.write(String.join("\t", inner));
		    writer.write("\n");
		    
		    writer.write("B\t");
		    writer.write(String.join("\t", bulge));
		    writer.write("\n");
		    
		    writer.write("M\t");
		    writer.write(String.join("\t", multi));
		    writer.write("\n");
		    
		    writer.write("D\t");
		    writer.write(String.join("\t", dangling));
		    writer.write("\n");
		    
		    writer.write("P\t");
		    writer.write(String.join("\t", paired));
		    
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		System.out.println("\tSaved profile to " + Paths.get(filename).toAbsolutePath().toString());
		
	}

	/**
	 * Generates a logo containing the probability distribution of the structural 
	 * contexts for each nucleotide position.
	 * @param number
	 * @param filename
	 * @param sequence
	 * @param probabilities
	 */
	private static void saveLogo(int number, String filename, String sequence, double[][] probabilities) {
		
		String[] sequence_split = sequence.split("");
		String[] nucleotide_positions = new String[sequence.length()];
		for (int x=0; x<sequence.length(); x++){ nucleotide_positions[x] = (x+1)+":"+sequence_split[x]; }
		
		Logo logo = new Logo(probabilities, nucleotide_positions);
		
		logo.setAlphabetContexts();
		logo.setBit(false);
		logo.saveAsPDF(25*sequence.length(), 150, filename);
		
		System.out.println("\tSaved profile to " + Paths.get(filename).toAbsolutePath().toString());
	}
	
	
	/**
	 * Converts linear representation of matrix to rectangular
	 * @return [[h1,h2,...,hn],
	 *          [i1,i2,...,in],
	 *          [b1,b2,...,bn],
	 *          [m1,m2,...,mn],
	 *          [d1,d2,...,dn],
	 *          [p1,p2,...,pn]]
	 */ 
	private static double[][] getMatrix(double[] profile, int length){
		
		double[][] matrix = new double[length][6];
		
		for(int x = 0; x<length; x++) {
			
			matrix[x][0] = profile[x+(0*length)]; //H
			matrix[x][1] = profile[x+(1*length)]; //I
			matrix[x][2] = profile[x+(2*length)]; //B
			matrix[x][3] = profile[x+(3*length)]; //M
			matrix[x][4] = profile[x+(4*length)]; //D
			matrix[x][5] = 1.0 - (profile[x+(0*length)]+profile[x+(1*length)]+profile[x+(2*length)]+profile[x+(3*length)]+profile[x+(4*length)]); //P
			
		}
	
		return matrix;
		
	}
	
}
