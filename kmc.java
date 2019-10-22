import java.io.IOException;
import java.io.PrintWriter;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.lang.Exception;
import java.util.*;

public class kmc {

	/*
	 * The main function is where the filenames are added to 
	 * a list and then the runClustering function is called
	 * for two configurations for each data set. 1 means
	 * to run it as described in 2-2, and 2 means to run it
	 * as described in 2-3.
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		List<String> filepaths = new ArrayList<String>();
		filepaths.add("./data/artdata0.5.arff");
		filepaths.add("./data/artdata1.arff");
		filepaths.add("./data/artdata2.arff");
		filepaths.add("./data/artdata3.arff");
		filepaths.add("./data/artdata4.arff");
		filepaths.add("./data/ionosphere.arff");
		filepaths.add("./data/iris.arff");
		filepaths.add("./data/soybean-processed.arff");
		for(int i = 0; i < filepaths.size(); i++) {
			runClustering(1, filepaths.get(i));
			runClustering(2, filepaths.get(i));
			System.out.println("DONE WITH " + filepaths.get(i));
		}
	}
	
	/*
	 * This function takes in a mode (1 or 2) and a filename.
	 * This is where the different cluster runs are configured.
	 * Mode 1 represents section 2-2 of the assignment, and mode
	 * 2 represents section 2-3 of the assignment. All of the 
	 * functions used in this code are called from this function.
	 * Once the arff file is parsed and the data is processed, based
	 * on the mode the clustering is performed differently, and the
	 * cluster scatter and NMI is returned in a list for mode 1, and
	 * the cluster scatter is returned in a list for mode 2 (-1 is 
	 * used for the NMI so it is ignored in this part). Then the 
	 * CS and NMI (for 2-2) is written to a text file named after
	 * the data set.
	 * 
	 */
	public static void runClustering(int mode, String filename) throws IOException {
		List<String> arffContents = arff_parser(filename); //returns a list of each line from the arff file
		List<String> attributes = getAttributes(arffContents); //returns a list of the attributes
		List<ArrayList<String>> data = cleanData(getData(arffContents)); //returns a list of lists of the data
		String[] classes = getClasses(arffContents); //returns the classes in a data set
		int numClasses = classes.length; //number of classes in a data set
		List<ArrayList<String>> dataWithoutClass = data; //dataWithoutClass is used to strip class info from data
		List<ArrayList<ArrayList<String>>> trueClusters = generateTrueClusters(data, classes); //generate the true clusters based on class data
		for(int i = 0; i < dataWithoutClass.size(); i++) {
			dataWithoutClass.get(i).remove(dataWithoutClass.get(i).size()-1); //remove class from each data point
		}
		if(mode == 1) { //2-2 random restart + smart initialization
			List<ArrayList<Double>> clusterStats = new ArrayList<ArrayList<Double>>(); //stores the CS and NMI for each run
			String outputFileName = "./output/" +filename.substring(7, filename.length() - 5) + "output22.txt"; //naming convention for output file
			PrintWriter writer = new PrintWriter(outputFileName, "UTF-8"); //writer for output file
			for(int i = 0; i < 10; i++) { //10 runs of random restart clustering
				List<ArrayList<Double>> rMeans = generateRandomMeans(dataWithoutClass, numClasses); //generate random initial means based on number of classes
				List<ArrayList<ArrayList<String>>> clusters = createClusters(rMeans, dataWithoutClass); //create clusters
				List<Double> cs_nmi = clustering(attributes, data, numClasses, clusters, trueClusters, mode); //output CS and NMI
				clusterStats.add((ArrayList<Double>) cs_nmi); //add CS and NMI to list
			}
			List<ArrayList<Double>> sMeans = generateSmartMeans(dataWithoutClass, numClasses); //generate the smart initialization means
			List<ArrayList<ArrayList<String>>> sClusters = createClusters(sMeans, dataWithoutClass); //create initial clusters
			List<Double> sCS_nmi = clustering(attributes, data, numClasses, sClusters, trueClusters, mode); //run clustering
			clusterStats.add((ArrayList<Double>) sCS_nmi); //add smart CS and NMI to list
			for(int i = 0; i < clusterStats.size(); i++) { //run through clusterStats list and output CS and NMI for each run
				writer.println(i + " " + clusterStats.get(i).get(0) + " " + clusterStats.get(i).get(1));
			}
			writer.close();
		}
		if(mode == 2) { //2-3 picking k
			List<Double> allCS = new ArrayList<Double>(); //list of all of the cluster scatters for a data set
			String outputFileName = "./output/" +filename.substring(7, filename.length() - 5) + "output23.txt"; //naming convention for output file
			PrintWriter writer = new PrintWriter(outputFileName, "UTF-8"); //writer for output file
			for(int k = 2; k <= 22; k++) {
				List<ArrayList<Double>> clusterStats = new ArrayList<ArrayList<Double>>(); //stores the CS for each run
				for(int i = 0; i < 10; i++) {
					List<ArrayList<Double>> rMeans = generateRandomMeans(dataWithoutClass, k); //generate random initial means based on k value
					List<ArrayList<ArrayList<String>>> clusters = createClusters(rMeans, dataWithoutClass); //create initial clusters
					List<Double> cs_nmi = clustering(attributes, data, k, clusters, trueClusters, mode); //output CS 
					clusterStats.add((ArrayList<Double>) cs_nmi);  //add cs to clusterStats
				}
				List<Double> clusterScatters = new ArrayList<Double>(); //create list to store only cluster scatters
				for(int i = 0; i < clusterStats.size(); i++) {
					clusterScatters.add(clusterStats.get(i).get(0)); //add cluster scatters to list
				}
				Double minCS = Collections.min(clusterScatters); //find minimum cluster scatter from the 10 runs
				allCS.add(minCS); //add minimum to list of cluster scatters for each k
			}
			for(int i = 0; i < allCS.size(); i++) {
				writer.println((i+2) + " " + allCS.get(i)); //print k and cluster scatters to text file
			}
			writer.close();
		}
	}
	
	/*
	 * This function reads in a file path and then stores each line in a list
	 * and returns that list.
	 */
	public static List<String> arff_parser(String filename) throws IOException{
		BufferedReader reader = new BufferedReader(new FileReader(filename)); //reads in arff file
		String temp; //represents a line in arff file
		List<String> fileList = new ArrayList<String>(); //list of each line in arff file
		while((temp = reader.readLine()) != null){
		    fileList.add(temp); //adds each line of arff file to list
		}
		reader.close();
		return fileList;
	}

	/*
	 * This function takes in a list of each line from the arff file and outputs
	 * a list of attributes for the data set.
	 */
	public static List<String> getAttributes(List<String> fileContents) {
		List<String> attributes = new ArrayList<String>(); //list of attributes
		for(int i = 0; i < fileContents.size(); i++) {
			if(fileContents.get(i).startsWith("@ATTRIBUTE") || fileContents.get(i).startsWith("@attribute")) { //keeps track of # attributes
				attributes.add(fileContents.get(i)); //add attribute to list of attributes
		    }
		}
		return attributes;
	}
	
	/*
	 * This function takes in a list of each line from the arff file and outputs
	 * a list of each line of data.
	 */
	public static List<String> getData(List<String> fileContents) {
		int dataIndex = 0; //will keep track of number of data in a set
		for(int i = 0; i < fileContents.size(); i++) {
			if(fileContents.get(i).startsWith("@DATA") || fileContents.get(i).startsWith("@data")) { //keeps track of # data
		        dataIndex = i+1; //index where the data starts
		    }
		}
		List<String> data = new ArrayList<String>(); //list of only the data in an arff file
		for(int j = 0; j < fileContents.size() - dataIndex; j++) {
			data.add(fileContents.get(dataIndex+j)); //adds data point to list of data
		}
		return data;
	}
	
	/*
	 * This function takes in a list of each line from the arff file and outputs
	 * the classes for the data set.
	 */
	public static String[] getClasses(List<String> fileContents) {
		int classIndex = 0;
		for(int i = 0; i < fileContents.size(); i++) {
			if(fileContents.get(i).startsWith("@ATTRIBUTE class") || 
					fileContents.get(i).startsWith("@attribute class")) { //find line where class info is
		        classIndex = i; //line number for class info
		    }
		}
		String classLine = fileContents.get(classIndex).substring(18); //line of the class info
		classLine = classLine.replaceAll(",", " "); //remove commas from line
		classLine = classLine.replaceAll("}", " "); //remove brackets from line
		String[] classes = classLine.split("\\s+"); //split line by spaces to get each class
		return classes;
	}
	
	/*
	 * This function takes in a list of the attributes, data, the number of clusters, the clusters
	 * created, the true clusters, and the mode (which section of the assignment), and outputs a 
	 * list of the CS (both sections and modes) and NMI (2-2, mode 1 only). 
	 */
	public static List<Double> clustering(List<String> attributes, List<ArrayList<String>> tempData, int numClusters,
			List<ArrayList<ArrayList<String>>> clusters, List<ArrayList<ArrayList<String>>> trueClusters, int mode) {
		int numAttributes = attributes.size()-1; //removes the class attribute from the counter
		Double oldCS = 0.0; //initialize CS value from previous run
		Double diffCS = 1.0; //value to calculate difference in new CS and old CS (when this is 0, then there
		//is no difference in the CS so that means clustering is over
		Double CS = 1.0; //current CS value
		int numRuns = 0; //number of iterations of the clustering loop
		while(diffCS > 0) { //loop where it continues until CS stops decreasing
			List<ArrayList<Double>> tempMeans = new ArrayList<ArrayList<Double>>(); //means for a particular run
			for(int q = 0; q < clusters.size(); q++) { //loop through each cluster
				Double[][] values = new Double[clusters.get(q).size()][numAttributes]; //2d array of data values
				for(int r = 0; r < clusters.get(q).size(); r++) { //convert from 2D list to 2D array
					for(int s = 0; s < numAttributes; s++) { //
						values[r][s] = Double.parseDouble(clusters.get(q).get(r).get(s));
					}
				}
				ArrayList<Double> means = new ArrayList<Double>(); //calculate mean of each cluster
				for(int a = 0; a < numAttributes; a++) { //sum up values for each attribute to get mean for each dimension
					Double tempSum = 0.0; //sum of values of a particular attribute
					for(int b = 0; b < values.length; b++) {
						tempSum += values[b][a];
					}
					Double numValues = (double) values.length; //number of values in a cluster
					Double mean = tempSum / numValues; //calculate mean for each attribute in a cluster
					means.add(mean); //list of means for a cluster
				}
				tempMeans.add(means); //list of means for all clusters
			}
			CS = calculateCS(clusters, tempMeans, numAttributes); //calculate the CS of a run
			clusters = createClusters(tempMeans, tempData); //create initial clusters for the next run
			diffCS = Math.abs(oldCS - CS); //calculate difference in cluster scatter values
			oldCS = CS; //make the CS equal to the oldCS variable
			numRuns++; //increment the number of runs
		}
		Double NMI = -1.0; //initial value for the NMI
		if(mode == 1) { //calculate NMI
			NMI = calculateNMI(clusters, trueClusters);
		}
		List<Double> cs_nmi = new ArrayList<Double>(); //list of the CS and NMI
		cs_nmi.add(CS); //add CS
		cs_nmi.add(NMI); //add NMI
		//System.out.println("CS: " + CS + " numRuns: " + numRuns + " numClusters: " + numClusters);
		return cs_nmi;
	}
	
	/*
	 * This function takes in a list of the data points and outputs a cleaned
	 * version of the list of data points. Commas are removed and spaces are 
	 * split so each value in a line can be stored separately. 
	 */
	public static List<ArrayList<String>> cleanData(List<String> data){
		List<ArrayList<String>> dData = new ArrayList<ArrayList<String>>(); //list of list of each data point
		for(int i = 0; i < data.size(); i++) {
			String dataPoint = data.get(i); //get a line of data
			dataPoint = dataPoint.replaceAll(",", " "); //remove commas
			dataPoint = dataPoint.substring(0, dataPoint.length()); 
			String[] dataPoints = dataPoint.split("\\s+"); //split based on spaces to store individual values
			ArrayList<String> attributeValues = new ArrayList<String>();
			for(int j = 0; j < dataPoints.length; j++) {
				attributeValues.add(dataPoints[j]); //store each attribute value for a data set in a list
			}
			dData.add(attributeValues); //add list of attribute values to list (data)
		}
		return dData;
	}
	
	/*
	 * This function takes in a list of initial means and a list of the data
	 * and outputs a list of clusters. 
	 */
	public static List<ArrayList<ArrayList<String>>> createClusters(
			List<ArrayList<Double>> means, List<ArrayList<String>> data){
		List<ArrayList<ArrayList<String>>> clusters = new ArrayList<ArrayList<ArrayList<String>>>(); //list of clusters
		for(int i = 0; i < means.size(); i++) {
			ArrayList<ArrayList<String>> initCluster = new ArrayList<ArrayList<String>>(); //initial cluster
			clusters.add(initCluster); //add initialized cluster to list of clusters
		}
		for(int j = 0; j < data.size(); j++) {
			List<Double> distances = new ArrayList<Double>(); //list of distances of a data point to each mean
			for(int k = 0; k < means.size(); k++) { //calculates distance for each mean
				Double distance = 0.0;
				for(int m = 0; m < data.get(j).size(); m++) { //adds distances for each dimension
					Double difference = Double.parseDouble(data.get(j).get(m)) - means.get(k).get(m);
					distance += Math.pow(difference, 2); //calculate sum of distances 
				}
				distance = Math.sqrt(distance); //square root sum of distances to get true distance of a data point from a mean
				distances.add(distance);
			}
			int minIndex = distances.indexOf(Collections.min(distances)); //find minimum distance
			clusters.get(minIndex).add(data.get(j)); //add the data point to the cluster where it had the minimum distance
		}
		return clusters;
	}
	
	/*
	 * This function takes in a clustering, a list of means, and the number of attributes
	 * for a data set, and outputs the cluster scatter of a clustering. 
	 */
	public static Double calculateCS(List<ArrayList<ArrayList<String>>> clusters, 
			List<ArrayList<Double>> means, int numAttributes) {
		Double CS = 0.0; //initialize cluster scatter
		for(int i = 0; i < clusters.size(); i++) { //loop through each cluster
			Double tempSumDiff = 0.0; //temporary sum of differences of each data point to mean
			for(int r = 0; r < clusters.get(i).size(); r++) { //loop through each data point in cluster
				for(int s = 0; s < numAttributes; s++) { //loop through each feature value in a data point
					Double tempValue = Double.parseDouble(clusters.get(i).get(r).get(s)); //convert value to double
					tempSumDiff += Math.pow(tempValue - means.get(i).get(s), 2); //add squared difference to sum
				}
			}
			CS += tempSumDiff; //cluster scatter equals sum of squared differences
		}
		CS = (double)Math.round(CS * 1000d) / 1000d; //round cluster scatter to three decimal places
		return CS;
	}
	
	/*
	 * This function generates random examples to be set as the means
	 * for the first run of the random restart clustering. This function
	 * takes in the data and the number of clusters to be made. First,
	 * it removes duplicates from the picking of the means, so there can't
	 * be two clusters set to the same mean. Then it shuffles the list of
	 * unique examples and selects the first numClusters of them. Then it
	 * adds those examples to the list of means, and returns the list of
	 * means rMeans.
	 */
	public static List<ArrayList<Double>> generateRandomMeans(List<ArrayList<String>> data, int numClusters) {
		Set<ArrayList<String>> uniqueDataPoints = new HashSet<ArrayList<String>>(data); //set of unique data points (no duplicates)
		List<ArrayList<String>> randomMeans = new ArrayList<ArrayList<String>>(); //list of random initialized means
		List<ArrayList<String>> dataPoints = new ArrayList<ArrayList<String>>(); //list of unique data points
		dataPoints.addAll(uniqueDataPoints); //add data points to list
		Collections.shuffle(dataPoints); //shuffle the list of data points
		randomMeans = dataPoints.subList(0, numClusters); //let the random means equal the first numClusters of the shuffled list
		List<ArrayList<Double>> rMeans = new ArrayList<ArrayList<Double>>(); //list of random means with string values converted to doubles
		for(int i = 0; i < randomMeans.size(); i++) {
			ArrayList<Double> tempMean = new ArrayList<Double>();
			for(int j = 0; j < randomMeans.get(i).size(); j++) {
				Double temp = Double.parseDouble(randomMeans.get(i).get(j)); //convert string to double
				tempMean.add(temp);
			}
			rMeans.add(tempMean);
		}
		return rMeans;
	}
	
	/*
	 * This function generates the smart initialization means described in 2-2. The
	 * function takes in the data, and the number of clusters. It then randomly picks
	 * an example to be the first mean, then repeats the following process until all
	 * means are picked. First, pick 10 random examples, then calculate the distance
	 * of each example to the means. Pick the example with the maximum value for 
	 * distance as the next mean, and repeat the process until all means are selected.
	 */
	public static List<ArrayList<Double>> generateSmartMeans(List<ArrayList<String>> data, int numClusters) {
		Collections.shuffle(data); //shuffle the data
		List<ArrayList<Double>> sMeans = new ArrayList<ArrayList<Double>>(); //list of the smart means
		ArrayList<String> firstMean = data.get(0); //shuffled data point, first data point selected
		ArrayList<Double> fMean = new ArrayList<Double>(); //first mean stored with doubles
		for(int i = 0; i < firstMean.size(); i++) { //convert first mean from string values to doubles
			Double temp = Double.parseDouble(firstMean.get(i));
			fMean.add(temp);
		}
		sMeans.add(fMean); //add first mean to list of smart means
		for(int j = 0; j < numClusters-1; j++) {
			List<ArrayList<String>> tenExamples = data.subList((10*j)+1, (10*j)+11); //select ten random examples
			List<Double> distances = new ArrayList<Double>();
			for(int k = 0; k < tenExamples.size(); k++) { //calculate distances for each example
				Double tempDistance = 0.0;
				for(int m = 0; m < sMeans.size(); m++) {
					for(int n = 0; n < sMeans.get(0).size(); n++) {
						tempDistance += Math.pow((Double.parseDouble(tenExamples.get(k).get(n)) //calculate distance
								- sMeans.get(m).get(n)),2);
					}
				}
				distances.add(Math.sqrt(tempDistance));
			}
			//select biggest distance, and add that example to mean
			int meanIndex = distances.indexOf(Collections.max(distances));
			ArrayList<String> nextMean = tenExamples.get(meanIndex); //choose example with greatest distance
			ArrayList<Double> nMean = new ArrayList<Double>();
			for(int i = 0; i < nextMean.size(); i++) {
				Double temp = Double.parseDouble(nextMean.get(i));
				nMean.add(temp);
			}
			sMeans.add(nMean); //add to list of means
		}
		return sMeans;
	}
	
	/*
	 * This function takes in a data set and a list of the classes for that data set,
	 * and generates what the clustering should be based on the class label, and
	 * returns the "true" clusters of a data set.
	 */
	public static List<ArrayList<ArrayList<String>>> generateTrueClusters(List<ArrayList<String>> data, String[] classes){
		List<ArrayList<ArrayList<String>>> trueClusters = new ArrayList<ArrayList<ArrayList<String>>>();
		//loop through each data point
		//see what last value in data point is 
		//based on that assign data to a cluster
		for(int j = 0; j < classes.length; j++) { //loop through each class to find values that are in that class
			ArrayList<ArrayList<String>> aCluster = new ArrayList<ArrayList<String>>();
			for(int k = 0; k < data.size(); k++) {
				String tempClass = data.get(k).get(data.get(k).size()-1); //get class value from a data point
				if(classes[j].equals(tempClass)) { //if class value equals a particular class, add to cluster for that class
					aCluster.add(data.get(k));
				}
			}
			trueClusters.add(aCluster); //add cluster to list of clusters
		}
		return trueClusters;
	}
	
	/*
	 * This function takes in two clusters and outputs the NMI of the
	 * two clusters. 
	 */
	public static Double calculateNMI(List<ArrayList<ArrayList<String>>> cluster1, 
			List<ArrayList<ArrayList<String>>> cluster2) {
		int c1Size = cluster1.size(); //number of clusters in the first clustering
		int c2Size = cluster2.size(); //number of clusters in the second clusting
		Double N = 0.0; //initialize N value for NMI calculation (number of data points)
		Double aValues[] = new Double[c1Size]; //stores sizes of each cluster in first clustering
		Double bValues[] = new Double[c2Size]; //stores sizes of each cluster in second clustering
		for(int i = 0; i < c1Size; i++) {
			N += cluster1.get(i).size(); //use to calculate N
			aValues[i] = (double) cluster1.get(i).size(); //store size of each cluster in first clustering
			bValues[i] = (double) cluster2.get(i).size(); //store size of each cluster in second clustering
		}
		Double valuesTable[][] = new Double[c1Size][c2Size]; //calculate counts for intersection of values in two clusters
		Double I_C1_C2 = 0.0; //mutual information for first cluster and second cluster
		for(int i = 0; i < c1Size; i++) { //loop through first clustering
			for(int j = 0; j < c2Size; j++) { //loop through second clustering
				ArrayList<ArrayList<String>> tempCluster1 = cluster1.get(i); //temporary cluster
				ArrayList<ArrayList<String>> tempCluster2 = cluster2.get(j); //temporary cluster
				Double tempCount = 0.0; //count number of times the two clusters share a data point
				for(int m = 0; m < tempCluster1.size(); m++) {
					for(int n = 0; n < tempCluster2.size(); n++) {
						if(tempCluster1.get(m).equals(tempCluster2.get(n))) {
							tempCount++; 
						}
					}
				}
				valuesTable[i][j] = tempCount; //store number of times two clusters share a data point
				if(valuesTable[i][j] == 0) { //don't include zero counts in NMI calculation since 0log0 is not real
					//ignore term
				}
				else {
					Double num = valuesTable[i][j]/N; //numerator of I calculation
					Double dem = aValues[i]*bValues[j]/Math.pow(N, 2); //denominator of I calculation
					I_C1_C2 += (valuesTable[i][j]/N) * Math.log(num/dem); //calculate mutual information
				}
			}
		}
		Double c1H = 0.0; //entropy for first cluster
		Double c2H = 0.0; //entropy for second cluster
		for(int i = 0; i < c1Size; i++) {
			c1H -= (aValues[i]/N)*Math.log(aValues[i]/N); //calculate entropy for first cluster
		}
		for(int i = 0; i < c2Size; i++) {
			c2H -= (bValues[i]/N)*Math.log(bValues[i]/N); //calculate entropy for second cluster
		}
		Double NMI = (2*I_C1_C2)/(c1H + c2H); //calculate NMI
		NMI = (double)Math.round(NMI * 1000d) / 1000d; //make there be three decimal places
		return NMI;
	}
}
