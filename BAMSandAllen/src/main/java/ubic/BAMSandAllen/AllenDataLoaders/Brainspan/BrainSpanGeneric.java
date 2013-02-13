/*******************************************************************************
 * The ABAMS project
 * 
 * Copyright (c) 2012 University of British Columbia
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *       http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ******************************************************************************/
package ubic.BAMSandAllen.AllenDataLoaders.Brainspan;

import java.io.File;
import java.io.FileReader;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.ABAMSDataMatrix;
import ubic.BAMSandAllen.RankedGeneListLoader;
import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.AllenDataLoaders.ImageSeriesInfoLoader;
import ubic.BAMSandAllen.Homologene.HomoloGeneLoader;
import ubic.BAMSandAllen.adjacency.CorrelationAdjacency;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.math.CorrelationStats;
import ubic.basecode.math.DescriptiveWithMissing;
import ubic.basecode.math.RandomChooser;
import ubic.basecode.util.FileTools;
import au.com.bytecode.opencsv.CSVReader;
import cern.colt.list.DoubleArrayList;

public class BrainSpanGeneric {

	protected String AGE_KEY = "age";

	protected String GENE_SYMBOL_KEY = "gene_symbol";

	private static Log log = LogFactory
			.getLog(BrainSpanGeneric.class.getName());

	protected String baseFolder = "/home/abams/Desktop/human/BrainSpan/";

	protected DoubleMatrix<Map<String, String>, Map<String, String>> data;
	protected String type;

	public Set<Map<String, String>> getRowBySymbol(String geneSymbol) {
		Set<Map<String, String>> result = new HashSet<Map<String, String>>();
		for (Map<String, String> row : data.getRowNames()) {
			if (row.get(getGeneSymbolKey()).equals(geneSymbol)) {
				result.add(row);
			}
		}
		return result;
	}

	public BrainSpanGeneric() throws Exception {

	}

	public BrainSpanGeneric(String type) throws Exception {
		this.type = type;
		baseFolder += type + "/";

		String prefix = type + "_";
		if (type.equals("array") || type.equals("genes")) {
			prefix = "";
		}

		// read in columns
		List<Map<String, String>> colMaps = getLines(baseFolder + prefix
				+ "columns_metadata.csv");

		// read in row
		List<Map<String, String>> rowMaps = getLines(baseFolder + prefix
				+ "rows_metadata.csv");

		log.info("rowMaps.size(), colMaps.size() " + rowMaps.size() + "x"
				+ colMaps.size());

		data = new DenseDoubleMatrix<Map<String, String>, Map<String, String>>(
				rowMaps.size(), colMaps.size());

		data.setColumnNames(colMaps);

		log.info("setting row maps");
		data.setRowNames(rowMaps);
		log.info("done rownames");

		if (type.equals("array") || type.equals("genes")) {
			prefix = "expression";
		}
		String dataFile = baseFolder + prefix + "_matrix.csv";
		CSVReader reader = new CSVReader(new FileReader(dataFile));
		String line[];
		int row = 0;
		while ((line = reader.readNext()) != null) {
			// log.info( line.length );
			if (row % 10000 == 0) {
				log.info("Reading row:" + row);
				// break;
			}
			// shift over by one because the first column is the gene ID
			// (probably should check that against the row
			// data
			// for blueprint the first col is row number
			for (int i = 1; i < line.length; i++) {
				data.set(row, i - 1, Double.parseDouble(line[i]));
			}
			row++;
		}

	}

	public String getType() {
		return type;
	}

	/**
	 * Gets the variables for each row or column
	 * 
	 * @param filename
	 * @return
	 * @throws Exception
	 */

	public List<Map<String, String>> getLines(String filename) throws Exception {
		log.info("Reading:" + filename);
		CSVReader reader = new CSVReader(new FileReader(filename));
		List<Map<String, String>> result = new LinkedList<Map<String, String>>();
		List<String[]> lines = (List<String[]>) reader.readAll();
		// grab the header line
		String[] keys = lines.get(0);
		lines.remove(0);
		for (String key : keys)
			log.info(key);
		for (String[] line : lines) {
			Map<String, String> lineMap = new HashMap<String, String>();
			for (int i = 0; i < line.length; i++) {
				lineMap.put(keys[i], line[i]);
			}
			result.add(lineMap);
		}
		reader.close();

		return result;
	}

	public Set<String> getUniqueGenes() {
		return getUniqueValues(data.getRowNames(), getGeneSymbolKey());
	}

	public Set<String> getUniqueRegions() {
		return getUniqueValues(data.getColNames(), "structure_name");
	}

	public Set<String> getUniqueValues(List<Map<String, String>> maps,
			String key) {
		Set<String> result = new HashSet<String>();
		for (Map<String, String> map : maps) {
			result.add(map.get(key));
		}
		return result;
	}

	public void writeAgeMatrices(Collection<String> genes, String name)
			throws Exception {
		writeAgeMatrices(genes, name, null, false);
	}

	public void writeAgeMatrices(Collection<String> genes, String name,
			boolean rank) throws Exception {
		writeAgeMatrices(genes, name, null, rank);
	}

	public void writeRegionMatrix(Collection<String> genes, String name,
			String age) throws Exception {
		Set<String> uniqueRegions = getUniqueRegions(); // not specific to age
		DoubleMatrix<String, String> matrix = new DenseDoubleMatrix<String, String>(
				genes.size(), uniqueRegions.size());
		matrix.setRowNames(new LinkedList<String>(genes));
		matrix.setColumnNames(new LinkedList<String>(uniqueRegions));

		Set<String> uniqueGenes = getUniqueGenes();
		log.info("Genes not in set:" + Util.subtract(genes, uniqueGenes).size());
		genes.retainAll(uniqueGenes);

		// for a given age and a given gene, we have a list - search all?
		for (String gene : genes) {
			for (String region : uniqueRegions) {
				List<Double> values = findData(gene, age, region);
				double[] valueArray = ArrayUtils.toPrimitive(values
						.toArray(new Double[0]));
				DoubleArrayList expValuesDAL = new DoubleArrayList(valueArray);

				double mean = DescriptiveWithMissing.mean(expValuesDAL);

				matrix.setByKeys(gene, region, mean);
			}
		}

		ABAMSDataMatrix result = new ABAMSDataMatrix(matrix, "result",
				new CorrelationAdjacency(matrix));
		result = result.removeZeroRows();

		String filename = baseFolder + name + ".byRegion.age." + age + "."
				+ type + ".size." + result.rows() + ".txt";
		Util.writeRTable(filename, result);
		log.info("File written:" + filename);
	}

	public void writeAgeMatrices(Collection<String> genes, String name,
			String region, boolean rank) throws Exception {
		Set<String> trimmedGenes = new HashSet<String>();
		for (String gene : genes)
			trimmedGenes.add(gene.trim());
		genes = trimmedGenes;

		// how many ages?
		List<String> ageValues = getTimePoints();

		DoubleMatrix<String, String> matrix = new DenseDoubleMatrix<String, String>(
				genes.size(), ageValues.size());

		matrix.setRowNames(new LinkedList<String>(genes));
		matrix.setColumnNames(ageValues);

		DoubleMatrix<String, String> resultSTD = new DenseDoubleMatrix<String, String>(
				genes.size(), ageValues.size());
		resultSTD.setRowNames(new LinkedList<String>(genes));
		resultSTD.setColumnNames(ageValues);

		Set<String> uniqueGenes = getUniqueGenes();
		log.info("Unique loaded genes:" + uniqueGenes.size());
		log.info("Genes not in set:" + Util.subtract(genes, uniqueGenes).size());
		log.info(Util.subtract(genes, uniqueGenes));
		genes.retainAll(uniqueGenes);

		// for a given age and a given gene, we have a list - search all?
		for (String gene : genes) {
			for (String age : ageValues) {
				List<Double> values = findData(gene, age, region, rank);
				double[] valueArray = ArrayUtils.toPrimitive(values
						.toArray(new Double[0]));
				DoubleArrayList expValuesDAL = new DoubleArrayList(valueArray);

				double mean = DescriptiveWithMissing.mean(expValuesDAL);
				double sampleStandardDeviation = Math
						.sqrt(DescriptiveWithMissing.sampleVariance(
								expValuesDAL, mean));

				matrix.setByKeys(gene, age, mean);
				resultSTD.setByKeys(gene, age, sampleStandardDeviation);
			}
		}

		ABAMSDataMatrix result = new ABAMSDataMatrix(matrix, "result",
				new CorrelationAdjacency(matrix));
		result = result.removeZeroRows();

		String filename = baseFolder + name + "." + type + ".size."
				+ result.rows() + ".txt";
		Util.writeRTable(filename, result);
		Util.writeRTable(baseFolder + name + "STD.size." + resultSTD.rows()
				+ ".txt", resultSTD);
		log.info("Files written:" + filename);
	}

	private List<String> getTimePoints() {
		Set<String> ageValuesSet = getUniqueValues(data.getColNames(), AGE_KEY);

		List<String> ageValues = new LinkedList<String>(ageValuesSet);
		Collections.sort(ageValues, new TimeSorter());
		return ageValues;
	}

	public List<Double> findData(String gene, String age, String region) {
		return findData(gene, age, region, false);
	}

	public List<Double> findData(String gene, String age, String region,
			boolean rank) {
		List<Double> result = new LinkedList<Double>();

		for (Map<String, String> col : data.getColNames()) {
			if (col.get(AGE_KEY).equals(age)) {
				if (region == null || col.get("structure_name").equals(region)) {

					List<Double> remaining = new LinkedList<Double>();
					double foundValue = Double.NaN;

					for (Map<String, String> row : data.getRowNames()) {
						if (row.get(getGeneSymbolKey()).equals(gene)) {
							foundValue = data.getByKeys(row, col);
							if (!rank)
								result.add(foundValue);
						} else if (rank) {
							remaining.add(data.getByKeys(row, col));
						}
					}
					if (rank) {
						double rankPos = 1;
						for (double otherValue : remaining) {
							if (otherValue > foundValue) {
								rankPos++;
							}
						}
						result.add(rankPos);
					}
				}
			}
		}
		return result;
	}

	public Map<String, List<Double>> findData(String gene, String age,
			Collection<String> regions) {
		Map<String, List<Double>> result = new HashMap<String, List<Double>>();
		for (Map<String, String> col : data.getColNames()) {
			if (col.get(AGE_KEY).equals(age)) {
				String region = col.get("structure_name");
				if (regions.contains(region)) {
					for (Map<String, String> row : data.getRowNames()) {
						if (row.get(getGeneSymbolKey()).equals(gene)) {
							List<Double> list = result.get(region);
							if (list == null) {
								list = new LinkedList<Double>();
								result.put(region, list);
							}
							list.add(data.getByKeys(row, col));
						}
					}
				}
			}
		}
		return result;
	}

	public DoubleMatrix<String, String> dump() {
		DoubleMatrix<String, String> result = new DenseDoubleMatrix<String, String>(
				data.asArray());
		List<String> colKeyOrder = new LinkedList<String>(data.getColName(0)
				.keySet());
		List<String> rowKeyOrder = new LinkedList<String>(data.getRowName(0)
				.keySet());

		for (Map<String, String> columnMap : data.getColNames()) {
			String columnString = "";
			for (String key : colKeyOrder) {
				columnString += "," + columnMap.get(key);
			}
			result.addColumnName(columnString.substring(1));
		}
		for (Map<String, String> rowMap : data.getRowNames()) {
			String rowString = "";
			for (String key : rowKeyOrder) {
				rowString += "," + rowMap.get(key);
			}
			result.addRowName(rowString.substring(1));
		}

		log.info("Column key:" + Arrays.toString(colKeyOrder.toArray()));
		log.info("Row key:" + Arrays.toString(rowKeyOrder.toArray()));
		return result;
	}

	public String getGeneSymbolKey() {
		return GENE_SYMBOL_KEY;
	}

	public DoubleMatrix<String, String> getRegionByAge(String gene, boolean rank) {
		// findData( gene, age, region );
		Set<String> regions = getUniqueRegions();

		Set<String> ages = getUniqueValues(data.getColNames(), AGE_KEY);

		log.info("Regions:" + regions.size() + " ages:" + ages.size());
		log.info(ages);

		DoubleMatrix<String, String> matrix = new DenseDoubleMatrix<String, String>(
				regions.size(), ages.size());
		matrix.setRowNames(new LinkedList<String>(regions));

		LinkedList<String> ageValues = new LinkedList<String>(ages);
		Collections.sort(ageValues, new TimeSorter());
		matrix.setColumnNames(ageValues);

		// does natural log before averaging
		for (String region : regions) {
			for (String age : ages) {
				List<Double> data = findData(gene, age, region, rank);
				double[] valueArray = ArrayUtils.toPrimitive(data
						.toArray(new Double[0]));
				if (!rank) {
					for (int i = 0; i < valueArray.length; i++) {
						// ln
						valueArray[i] = Math.log(valueArray[i]);
					}
				}
				DoubleArrayList expValuesDAL = new DoubleArrayList(valueArray);
				double mean = DescriptiveWithMissing.mean(expValuesDAL);

				matrix.setByKeys(region, age, mean);
			}
		}
		return matrix;
	}

	public DoubleMatrix<String, String> getGeneData(String gene) {
		// make a new matrix?
		DoubleMatrix<String, String> result;

		List<Map<String, String>> rowsWithGene = new LinkedList<Map<String, String>>();
		for (Map<String, String> row : data.getRowNames()) {
			log.info(getGeneSymbolKey());
			if (row.get(getGeneSymbolKey()).equals(gene)) {

				rowsWithGene.add(row);
			}
		}
		log.info("Rows/exons for gene:" + rowsWithGene.size());
		result = new DenseDoubleMatrix<String, String>(rowsWithGene.size(),
				data.columns());

		for (Map<String, String> row : rowsWithGene) {
			String exonInfo = "";
			if (type.startsWith("exon")) {
				exonInfo = "exon." + row.get("start") + ".to." + row.get("end")
						+ ".";
			}

			String rowName = exonInfo + row.get(getGeneSymbolKey());

			result.addRowName(rowName);
			for (int col = 0; col < data.columns(); col++) {
				result.set(result.getRowIndexByName(rowName), col,
						data.get(data.getRowIndexByName(row), col));
			}
		}

		// colapse column data
		for (Map<String, String> col : data.getColNames()) {
			String colName = col.get(AGE_KEY) + "." + col.get("gender") + "."
					+ col.get("structure_id") + "." + col.get("donor_name");
			result.addColumnName(colName);
		}
		return result;
	}

	/**
	 * Simple correlation between rows
	 * 
	 * @param gene
	 * @throws Exception
	 */
	public DoubleMatrix<String, String> correlate(String gene) throws Exception {
		Set<Map<String, String>> geneRows = getRowBySymbol(gene);
		log.info("Number of rows:" + geneRows.size());
		log.info("Warning: Using only one");
		Map<String, String> geneRow = geneRows.iterator().next();
		double[] dataRow = data.getRowByName(geneRow);
		return correlate(dataRow, gene);

	}

	public DoubleMatrix<String, String> correlateToTime() throws Exception {
		double[] dataRow = new double[data.columns()];
		// Vector for time
		for (int i = 0; i < dataRow.length; i++) {
			dataRow[i] = i;
		}
		return correlate(dataRow, "Developmental Intervals");
	}

	public DoubleMatrix<String, String> correlateViaTime(
			DoubleMatrix<String, String> a, DoubleMatrix<String, String> b) {
		List<String> ageValues = getTimePoints();

		DoubleMatrix<String, String> correlations = new DenseDoubleMatrix<String, String>(
				1, ageValues.size());
		correlations.addRowName("Spearman");

		correlations.setColumnNames(ageValues);

		for (String agePoint : ageValues) {
			double[] aVec = a.getColumnByName(agePoint);
			double[] bVec = b.getColumnByName(agePoint);
			double spearmanCor = Util.spearmanCorrel(aVec, bVec);

			correlations.setByKeys("Spearman", agePoint, spearmanCor);
		}
		return correlations;
	}

	// TODO: filter correlations by genes of interest
	public DoubleMatrix<String, String> correlate(double[] dataRow, String name)
			throws Exception {

		DoubleMatrix<String, String> correlations = new DenseDoubleMatrix<String, String>(
				getUniqueGenes().size(), 4);
		correlations.setRowNames(new LinkedList<String>(getUniqueGenes()));
		correlations.addColumnName("Spearman " + name);
		correlations.addColumnName("Pearson " + name);
		correlations.addColumnName("Spearman P " + name);
		correlations.addColumnName("Pearson P " + name);

		for (Map<String, String> rowToTest : data.getRowNames()) {
			double[] testRow = data.getRowByName(rowToTest);
			double spearmanCor = Util.spearmanCorrel(dataRow, testRow);
			double spearmanP = CorrelationStats.spearmanPvalue(spearmanCor,
					testRow.length) * data.rows();
			double pearsonCor = CorrelationStats.correl(dataRow, testRow);
			double pearsonP = CorrelationStats.pvalue(pearsonCor,
					testRow.length) * data.rows();
			String geneTested = rowToTest.get(getGeneSymbolKey());
			correlations.setByKeys(geneTested, "Spearman " + name, spearmanCor);
			correlations.setByKeys(geneTested, "Pearson " + name, pearsonCor);
			correlations.setByKeys(geneTested, "Spearman P " + name, spearmanP);
			correlations.setByKeys(geneTested, "Pearson P " + name, pearsonP);
		}
		return correlations;
	}

	public Set<String> getRandomCoronalSubset(int size) throws Exception {
		ImageSeriesInfoLoader infoLoader = new ImageSeriesInfoLoader();
		HomoloGeneLoader homoLoader = new HomoloGeneLoader();
		Set<String> result = infoLoader.getAllRowsWithCoronalGenes();
		RankedGeneListLoader loader = new RankedGeneListLoader(
				new LinkedList<String>(result), "x");
		loader.convertToNCBIIDs();
		Set<String> humanSymbols = homoLoader.getHumanSymbolFromMouseID(loader
				.getLines());
		humanSymbols = RandomChooser.chooseRandomSubset(size, humanSymbols);
		return humanSymbols;
	}

	public Set<String> getRandomSubsetGenes(int size) {
		Set<String> result = new HashSet<String>();
		RandomChooser.init(1);

		for (Map<String, String> row : RandomChooser.chooseRandomSubset(size,
				data.getRowNames())) {
			result.add(row.get(getGeneSymbolKey()));
		}
		return result;
	}

	public void makeAgeMaricesForPLoSSets(boolean rank) throws Exception {
		String base = SetupParameters.config.getString("abams.dataFolder")
				+ "rankedGenes/near final ammon/";
		String outFile = base
				+ "LOOGenesInOrder.out.partialcon.ammon.txt.374.0.016424.topGenes.txt";
		String inFile = base
				+ "LOOGenesInOrder.in.partialcon.ammon.txt.452.0.01985.topGenes.txt";
		String proxFile = base
				+ "LOOGenesInOrder.space.ammon.txt.420.0.018435.topGenes.txt";

		writeMatrixForMouseList(outFile, "out.partial", rank);
		writeMatrixForMouseList(inFile, "in.partial", rank);
		writeMatrixForMouseList(proxFile, "proximity", rank);

	}

	public void writeMatrixForMouseList(String rankedFile, String name,
			boolean rank) throws Exception {
		RankedGeneListLoader rankedLoader = new RankedGeneListLoader(rankedFile);
		rankedLoader.convertToNCBIIDs();
		List<String> mouseIDs = rankedLoader.getLines();
		HomoloGeneLoader homoLoader = new HomoloGeneLoader();
		Set<String> humanSymbols = homoLoader
				.getHumanSymbolFromMouseID(mouseIDs);
		log.info("Human genes:" + humanSymbols.size());
		writeAgeMatrices(humanSymbols, name, rank);
	}

	public void writeHumanSets(String filename) throws Exception {
		// List<String> rankedFile = FileTools.getLines( filename );
		RankedGeneListLoader rankedLoader = new RankedGeneListLoader(filename);
		rankedLoader.convertToNCBIIDs();
		List<String> mouseIDs = rankedLoader.getLines();

		HomoloGeneLoader homoLoader = new HomoloGeneLoader();
		File f = new File(filename);
		// f.getName();
		Set<String> humanSymbols = homoLoader
				.getHumanSymbolFromMouseID(mouseIDs);
		FileTools.stringsToFile(
				humanSymbols,
				"/home/leon/Downloads/BrainSpan/human gene lists/"
						+ f.getName());
	}

	public static void main(String[] args) throws Exception {
		BrainSpanGeneric loader;
		DoubleMatrix<String, String> result;

		loader = new BrainSpanGeneric("genes");
		List<String> elisaGenes = FileTools
				.getLines("/home/abams/Desktop/human/elisa/wntGenes.txt");
		// );loader.writeAgeMatrices(loader.getUniqueGenes(), "Elisa");

		String region;
		boolean rank = true;
		region = "orbital frontal cortex";
		loader.writeAgeMatrices(elisaGenes, "WntElisa.ranked." + region, region, rank);
		region = "dorsolateral prefrontal cortex";
		loader.writeAgeMatrices(elisaGenes, "WntElisa.ranked." + region, region, rank);
		region = "ventrolateral prefrontal cortex";
		loader.writeAgeMatrices(elisaGenes, "WntElisa.ranked." + region, region, rank);

		for (String region2 : loader.getUniqueRegions()) {
			System.out.println(region2);
		}

		System.exit(1);

		DoubleMatrix<String, String> snca = loader.getRegionByAge("SNCA", true);
		DoubleMatrix<String, String> cd45 = loader
				.getRegionByAge("HLA-F", true);

		// Util.writeRTable("/home/abams/Desktop/human/BrainSpan/SNCA/HLAF.tsv",
		// cd45);
		// System.exit(1);

		LinkedList<String> emailedGenes = new LinkedList<String>();
		emailedGenes.add("IFNGR1");
		emailedGenes.add("HLA-E");
		emailedGenes.add("IRF7");
		emailedGenes.add("GBP1");
		emailedGenes.add("HLA-B");
		emailedGenes.add("CD44");
		emailedGenes.add("HLA-A");
		emailedGenes.add("IRF1");
		emailedGenes.add("B2M");
		emailedGenes.add("IRF6");
		emailedGenes.add("HLA-F");
		emailedGenes.add("SP100");
		emailedGenes.add("HLA-DPB1");
		emailedGenes.add("OASL");
		emailedGenes.add("SOCS3");
		emailedGenes.add("GBP2");
		emailedGenes.add("SOCS1");
		emailedGenes.add("IFI30");

		Set<String> interferonGOgroup = new HashSet<String>();
		interferonGOgroup.add("IRF1");
		interferonGOgroup.add("JAK2");
		interferonGOgroup.add("SOCS3");
		interferonGOgroup.add("IRF6");
		interferonGOgroup.add("SOCS1");
		interferonGOgroup.add("JAK2");
		interferonGOgroup.add("PIAS1");
		interferonGOgroup.add("OAS1");
		interferonGOgroup.add("IFNG");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-A");
		interferonGOgroup.add("HLA-DRA");
		interferonGOgroup.add("HLA-DQA2");
		interferonGOgroup.add("HLA-DQA1");
		interferonGOgroup.add("HLA-DRB1");
		interferonGOgroup.add("HLA-DRB1");
		interferonGOgroup.add("HLA-DQB1");
		interferonGOgroup.add("MT2A");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-C");
		interferonGOgroup.add("HLA-DRB1");
		interferonGOgroup.add("HLA-A");
		interferonGOgroup.add("HLA-DPB1");
		interferonGOgroup.add("ICAM1");
		interferonGOgroup.add("HLA-DQB2");
		interferonGOgroup.add("HCK");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-C");
		interferonGOgroup.add("IRF1");
		interferonGOgroup.add("IRF1");
		interferonGOgroup.add("FCGR1A");
		interferonGOgroup.add("IFI30");
		interferonGOgroup.add("NCAM1");
		interferonGOgroup.add("HLA-E");
		interferonGOgroup.add("HLA-DRB1");
		interferonGOgroup.add("HLA-DRB1");
		interferonGOgroup.add("HLA-DRB4");
		interferonGOgroup.add("IRF2");
		interferonGOgroup.add("IFNGR1");
		interferonGOgroup.add("CD44");
		interferonGOgroup.add("HLA-G");
		interferonGOgroup.add("PTPN2");
		interferonGOgroup.add("PTPN1");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("VCAM1");
		interferonGOgroup.add("HLA-DPA1");
		interferonGOgroup.add("HLA-DRB1");
		interferonGOgroup.add("JAK1");
		interferonGOgroup.add("SP100");
		interferonGOgroup.add("SP100");
		interferonGOgroup.add("PTAFR");
		interferonGOgroup.add("PTPN6");
		interferonGOgroup.add("PML");
		interferonGOgroup.add("OAS2");
		interferonGOgroup.add("HLA-A");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-C");
		interferonGOgroup.add("HLA-C");
		interferonGOgroup.add("HLA-C");
		interferonGOgroup.add("HLA-C");
		interferonGOgroup.add("HLA-C");
		interferonGOgroup.add("HLA-C");
		interferonGOgroup.add("HLA-F");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("GBP1");
		interferonGOgroup.add("GBP2");
		interferonGOgroup.add("CIITA");
		interferonGOgroup.add("IFNGR2");
		interferonGOgroup.add("STAT1");
		interferonGOgroup.add("STAT1");
		interferonGOgroup.add("B2M");
		interferonGOgroup.add("SUMO1");
		interferonGOgroup.add("HLA-DRB3");
		interferonGOgroup.add("IRF9");
		interferonGOgroup.add("IRF8");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("PRKCD");
		interferonGOgroup.add("PTPN11");
		interferonGOgroup.add("HLA-C");
		interferonGOgroup.add("CAMK2B");
		interferonGOgroup.add("CAMK2G");
		interferonGOgroup.add("CAMK2D");
		interferonGOgroup.add("IRF5");
		interferonGOgroup.add("IRF3");
		interferonGOgroup.add("IRF4");
		interferonGOgroup.add("OASL");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-C");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-C");
		interferonGOgroup.add("HLA-C");
		interferonGOgroup.add("HLA-DRB1");
		interferonGOgroup.add("HLA-DRB1");
		interferonGOgroup.add("HLA-DRB5");
		interferonGOgroup.add("HLA-DRB1");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-DRB1");
		interferonGOgroup.add("FCGR1B");
		interferonGOgroup.add("IRF7");
		interferonGOgroup.add("HLA-B");
		interferonGOgroup.add("HLA-C");
		interferonGOgroup.add("HLA-DRB1");
		interferonGOgroup.add("HLA-DRB1");
		interferonGOgroup.add("HLA-C");
		interferonGOgroup.add("HLA-DRB1");
		interferonGOgroup.add("CAMK2A");
		interferonGOgroup.add("OAS3");
		interferonGOgroup.add("SOCS3");
		interferonGOgroup.add("SOCS1");
		interferonGOgroup.add("JAK2");
		interferonGOgroup.add("PIAS1");
		interferonGOgroup.add("IFNG");
		interferonGOgroup.add("HSP90AB1");
		interferonGOgroup.add("IFNGR1");
		interferonGOgroup.add("PTPN2");
		interferonGOgroup.add("PTPN1");
		interferonGOgroup.add("JAK1");
		interferonGOgroup.add("PTPN6");
		interferonGOgroup.add("IFNGR2");
		interferonGOgroup.add("STAT1");
		interferonGOgroup.add("SUMO1");
		interferonGOgroup.add("PTPN11");
		interferonGOgroup.add("CDC37");
		interferonGOgroup.add("PPARG");
		interferonGOgroup.add("NR1H2");
		interferonGOgroup.add("NR1H3");
		interferonGOgroup.add("DNAJA3");
		interferonGOgroup.add("DNAJA3");
		interferonGOgroup.add("HPX");
		interferonGOgroup.add("TXK");
		interferonGOgroup.add("MED1");
		interferonGOgroup.add("NLRC5");
		interferonGOgroup.add("NLRC5");

		String filename;
		createCorrelationMatrix(loader, snca, interferonGOgroup, "interferon");

		HashSet<String> ubList = new HashSet<String>();
		ubList.add("RNF4");
		ubList.add("RNF4");
		ubList.add("RNF4");
		ubList.add("WWP2");
		ubList.add("WWP2");
		ubList.add("WWP2");
		ubList.add("WWP2");
		ubList.add("PARK2");
		ubList.add("PARK2");
		ubList.add("RNF8");
		ubList.add("ARIH2");
		ubList.add("NEDD4");
		ubList.add("UBE2N");
		ubList.add("UBE2N");
		ubList.add("UBE2N");
		ubList.add("UBE2G1");
		ubList.add("UBE2B");
		ubList.add("RNF4");
		ubList.add("TRAF2");
		ubList.add("TRAF2");
		ubList.add("UBE2V1");
		ubList.add("UBE2V1");
		ubList.add("UBE2S");
		ubList.add("RNF168");
		ubList.add("CDKN2A");
		ubList.add("UBE2E3");
		ubList.add("ITCH");
		ubList.add("ITCH");
		ubList.add("ZFP91");
		ubList.add("UBE2E2");
		ubList.add("RNF5");
		ubList.add("TRIM56");
		ubList.add("TRIM5");
		ubList.add("UBE2T");
		ubList.add("STUB1");
		ubList.add("UBE2D4");
		ubList.add("TRAF6");
		ubList.add("TRAF6");
		ubList.add("TRAF6");
		ubList.add("RNF4");
		ubList.add("RNF4");
		ubList.add("RNF4");
		ubList.add("UBE2C");
		ubList.add("PARK2");
		ubList.add("MARCH6");
		ubList.add("RNF8");
		ubList.add("RNF8");
		ubList.add("RNF8");
		ubList.add("ARIH2");
		ubList.add("ARIH2");
		ubList.add("TNFAIP3");
		ubList.add("CDC34");
		ubList.add("UBE2A");
		ubList.add("UBE2D1");
		ubList.add("UBE2E1");
		ubList.add("TTC3");
		ubList.add("UBE2G2");
		ubList.add("UBE2D3");
		ubList.add("UBE2K");
		ubList.add("UBE2G1");
		ubList.add("UBE2H");
		ubList.add("UBE2D2");
		ubList.add("UBE2B");
		ubList.add("RNF4");
		ubList.add("UBE3A");
		ubList.add("RNF187");
		ubList.add("UBE2R2");
		ubList.add("HACE1");
		ubList.add("RNF152");
		ubList.add("ZNRF1");
		ubList.add("UBE2Q2");
		ubList.add("UBE2E3");
		ubList.add("PELI1");
		ubList.add("ITCH");
		ubList.add("UBE2E2");
		ubList.add("RNF5");
		ubList.add("UBE2T");
		ubList.add("TOPORS");
		ubList.add("RNF146");
		ubList.add("RNF216");
		ubList.add("RNF6");
		ubList.add("UBE2D4");
		ubList.add("UBE2S");
		ubList.add("ITCH");
		ubList.add("ITCH");
		ubList.add("UBE2T");
		ubList.add("UBE2D4");
		ubList.add("UBE2S");
		ubList.add("UBE2T");
		ubList.add("RNF6");
		ubList.add("UBE2D4");
		ubList.add("RNF8");
		ubList.add("RNF8");
		ubList.add("RNF168");
		ubList.add("RNF168");
		ubList.add("RNF168");
		ubList.add("RNF168");
		ubList.add("RNF4");
		ubList.add("RNF4");
		ubList.add("RNF4");
		ubList.add("UBE2C");
		ubList.add("UBE2C");
		ubList.add("UBE2C");
		ubList.add("CDC27");
		ubList.add("UBE2A");
		ubList.add("UBE2D3");
		ubList.add("UBE2H");
		ubList.add("UBE2B");
		ubList.add("UBE2L3");
		ubList.add("RNF4");
		ubList.add("CDC16");
		ubList.add("UBE2S");
		ubList.add("UBE2S");
		ubList.add("CDC26");
		ubList.add("UBE2E3");
		ubList.add("UBE2W");
		ubList.add("UBE2E2");
		ubList.add("ANAPC13");
		ubList.add("ANAPC1");
		ubList.add("UBE2T");
		ubList.add("ANAPC11");
		ubList.add("CDC23");
		ubList.add("ANAPC7");
		ubList.add("ANAPC5");
		ubList.add("ANAPC4");
		ubList.add("ANAPC2");
		ubList.add("FZR1");
		ubList.add("ANAPC10");
		ubList.add("UBE2D4");
		ubList.add("RNF4");
		ubList.add("RNF4");
		ubList.add("RNF4");
		ubList.add("BRCA1");
		ubList.add("BRCA1");
		ubList.add("RNF4");
		ubList.add("UBE2S");
		ubList.add("BARD1");
		ubList.add("BARD1");
		ubList.add("UBE2T");
		ubList.add("RNF6");
		ubList.add("UBE2D4");
		ubList.add("UBR5");
		ubList.add("TRIP12");
		ubList.add("OTUB1");
		ubList.add("UBOX5");
		ubList.add("PSMB11");
		ubList.add("HDAC6");
		ubList.add("UBE4A");
		ubList.add("HDAC6");
		ubList.add("HDAC6");
		ubList.add("MARCH1");
		ubList.add("HDAC6");
		ubList.add("HDAC6");
		ubList.add("HDAC6");
		ubList.add("G2E3");
		ubList.add("CHFR");
		ubList.add("CHFR");
		ubList.add("G2E3");
		ubList.add("G2E3");
		ubList.add("G2E3");
		ubList.add("G2E3");
		ubList.add("RNF111");
		ubList.add("FBXO22");
		ubList.add("FBXO22");
		ubList.add("FBXO22");
		ubList.add("PSMD11");
		ubList.add("PSMD12");
		ubList.add("PSMD9");
		ubList.add("PSMD14");
		ubList.add("PSMA7");
		ubList.add("AXIN1");
		ubList.add("PSMD3");
		ubList.add("PARK2");
		ubList.add("PARK2");
		ubList.add("PARK2");
		ubList.add("MGRN1");
		ubList.add("PSMD10");
		ubList.add("UBOX5");
		ubList.add("UBR5");
		ubList.add("UBE4B");
		ubList.add("TNKS");
		ubList.add("TNKS");
		ubList.add("ARIH2");
		ubList.add("ARIH2");
		ubList.add("UBB");
		ubList.add("UBC");
		ubList.add("BCL2");
		ubList.add("PSMC3");
		ubList.add("TRIM21");
		ubList.add("PSMB1");
		ubList.add("PSMA1");
		ubList.add("PSMA2");
		ubList.add("PSMA3");
		ubList.add("PSMA4");
		ubList.add("PSMB8");
		ubList.add("PSMB9");
		ubList.add("PSMA5");
		ubList.add("PSMB4");
		ubList.add("PSMB6");
		ubList.add("PSMB5");
		ubList.add("TPP2");
		ubList.add("PSMC2");
		ubList.add("PSMB10");
		ubList.add("PSMC4");
		ubList.add("PSMD8");
		ubList.add("CDC34");
		ubList.add("UBE2A");
		ubList.add("PSMB3");
		ubList.add("PSMB2");
		ubList.add("PSMD7");
		ubList.add("UBE2D1");
		ubList.add("UBE2D1");
		ubList.add("UBE2E1");
		ubList.add("UBE2E1");
		ubList.add("PSMD4");
		ubList.add("NPEPPS");
		ubList.add("PSMA6");
		ubList.add("UBE2D3");
		ubList.add("UBE2D3");
		ubList.add("PSME3");
		ubList.add("PSMC1");
		ubList.add("PSMC5");
		ubList.add("PSMC6");
		ubList.add("UBE2D2");
		ubList.add("RPS27A");
		ubList.add("UBA52");
		ubList.add("UBE2B");
		ubList.add("UBE2L3");
		ubList.add("UBE2L3");
		ubList.add("PSME1");
		ubList.add("TRIM32");
		ubList.add("PSMD2");
		ubList.add("ERCC8");
		ubList.add("MAP3K1");
		ubList.add("SKP2");
		ubList.add("PPIL2");
		ubList.add("UBE2V1");
		ubList.add("BIRC2");
		ubList.add("BLMH");
		ubList.add("UBE4A");
		ubList.add("PSME4");
		ubList.add("PSMD6");
		ubList.add("UBE3C");
		ubList.add("UBE2V2");
		ubList.add("PSMD5");
		ubList.add("MARCH8");
		ubList.add("BTRC");
		ubList.add("RNF20");
		ubList.add("LRSAM1");
		ubList.add("NHLRC1");
		ubList.add("RNF111");
		ubList.add("G2E3");
		ubList.add("HUWE1");
		ubList.add("DZIP3");
		ubList.add("CDKN2A");
		ubList.add("ZNF738");
		ubList.add("ZSWIM2");
		ubList.add("FBXO22");
		ubList.add("PSMA8");
		ubList.add("MARCH1");
		ubList.add("DDB2");
		ubList.add("PSMF1");
		ubList.add("RNF31");
		ubList.add("RNF31");
		ubList.add("CHFR");
		ubList.add("PSMB7");
		ubList.add("PSMD1");
		ubList.add("RUSC1");
		ubList.add("RBCK1");
		ubList.add("RBCK1");
		ubList.add("TNKS2");
		ubList.add("TNKS2");
		ubList.add("RNF41");
		ubList.add("RNF41");
		ubList.add("RNF167");
		ubList.add("SMURF1");
		ubList.add("MARCH5");
		ubList.add("DTL");
		ubList.add("DTL");
		ubList.add("HDAC6");
		ubList.add("LNPEP");
		ubList.add("FBXW11");
		ubList.add("FBXO4");
		ubList.add("FBXO4");
		ubList.add("AMFR");
		ubList.add("AMFR");
		ubList.add("PSME2");
		ubList.add("PRPF19");
		ubList.add("STUB1");
		ubList.add("STUB1");
		ubList.add("PSMD13");
		ubList.add("BTRC");
		ubList.add("TRAF6");
		ubList.add("TRAF6");
		ubList.add("RNF31");
		ubList.add("RNF31");
		ubList.add("RNF31");
		ubList.add("RBCK1");
		ubList.add("RBCK1");
		ubList.add("RBCK1");
		ubList.add("SHARPIN");
		ubList.add("SHARPIN");
		ubList.add("SHARPIN");
		// createCorrelationMatrix(loader, snca, ubList,
		// "protein polyubiquitination");
		// // DoubleMatrix<String, String> result =
		// loader.correlateViaTime(snca,
		// // cd45);
		// // String filename =
		// "/home/abams/Desktop/human/BrainSpan/SNCA/correlations/cd45.tsv";
		// Util.writeRTable(filename, result);
		System.exit(1);
		loader = new BrainSpanGeneric("blueprint");
		Set<String> caml = new HashSet<String>();
		caml.add("LMACSD1");
		caml.add("ZEB2");
		loader.writeRegionMatrix(caml, "LMACSD1", "3 mo");
		// loader.writeRegionMatrix( caml, "LMACSD1", "12 mo" );
		// loader.writeRegionMatrix( caml, "LMACSD1", "48 mo" );
		System.exit(1);

		loader = new BrainSpanGeneric("genes");

		BrainSpanGeneric lmd = new BrainSpanGeneric("genes");
		result = lmd.getRegionByAge("LMACSD1", false);

		filename = "/home/leon/Downloads/prenatalLMDmicroarray/LMACSD1/matrix.brainspan.csv";
		Util.writeRTable(filename, result);

		System.exit(1);
		log.info("Unique genes:" + loader.getUniqueGenes().size());

		String base = SetupParameters.config.getString("abams.dataFolder")
				+ "rankedGenes/near final ammon/";

		// loader.writeHumanSets( base +
		// "LOOGenesInOrder.out.partialcon.ammon.txt.374.0.016424.topGenes.txt"
		// );
		// loader.writeHumanSets( base +
		// "LOOGenesInOrder.in.partialcon.ammon.txt.452.0.01985.topGenes.txt" );
		// System.exit( 1 );

		String geneCor = "SLC20A2";
		// Util.writeRTable( "/home/leon/Downloads/SLCCorrelations.txt",
		// loader.correlate( geneCor ) );
		// Util.writeRTable( "/home/leon/Downloads/TimeCorrelations.txt",
		// loader.correlateToTime() );
		// System.exit( 1 );

		// Set<String> SZUp = new HashSet<String>( FileTools.getLines(
		// "/home/leon/Downloads/query_szUp.txt" ) );
		// loader.writeAgeMatrices( SZUp, "SZUp",
		// "Dorsolateral prefrontal cortex" );
		//
		// Set<String> SZDown = new HashSet<String>( FileTools.getLines(
		// "/home/leon/Downloads/query_szDown.txt" ) );
		// loader.writeAgeMatrices( SZDown, "SZDown",
		// "Dorsolateral prefrontal cortex" );
		// loader.writeAgeMatrices( loader.getRandomCoronalSubset( 300 ),
		// "RandomCoronal" );

		// Set<String> axonGuidance = new HashSet<String>(
		// FileTools.getLines( "/home/leon/Downloads/axon guidance genes.txt" )
		// );
		// loader.writeAgeMatrices( axonGuidance, "axonGuidance" );
		// System.exit( 1 );
		//
		// for ( String region : loader.getUniqueRegions() ) {
		// System.out.println( region );
		// }
		//
		// loader.makeAgeMaricesForPLoSSets();

		// System.exit( 1 );

		// // TODO Auto-generated method stub
		// BrainSpanLoader loader = new BrainSpanLoader();
		// String gene = "SLC20A2";
		// DoubleMatrix<String, String> geneMatrix = loader.getGeneData( gene );
		// Util.writeRTable( "/home/leon/Downloads/BrainSpan/" + gene + ".txt",
		// geneMatrix );

		// String gene = "SLC20A2";
		// DoubleMatrix<String, String> geneMatrix = loader.getGeneData( gene );
		// Util.writeRTable( "/home/leon/Downloads/BrainSpan/" + gene + "." +
		// loader.getType() + ".txt", geneMatrix );

		// Set<String> geneSet = new HashSet<String>();
		// geneSet.add( "SLC20A2" );
		//
		// for ( String gene : loader.getUniqueGenes() ) {
		// if ( gene.startsWith( "SEMA" ) ) geneSet.add( gene );
		// }
		//
		// loader.writeAgeMatrices( geneSet, "Sema" );
		//
		// List<String> NEgeneSet = FileTools.getLines(
		// "/home/leon/Downloads/NEHuman.txt" );
		// loader.writeAgeMatrices( NEgeneSet, "NE" );
		//
		// List<String> OEgeneSet = FileTools.getLines(
		// "/home/leon/Downloads/OEHuman.txt" );
		// loader.writeAgeMatrices( OEgeneSet, "OE" );

		// List<String> MeetaUp = new LinkedList<String>( new HashSet<String>(
		// FileTools.getLines( "/home/leon/Downloads/query_szUp (1).txt" ) ) );
		// loader.writeAgeMatrices( MeetaUp, "MeetaUp" );

		String file = "/home/leon/Downloads/query_szDown (3).txt";

		// file = "/home/leon/MaternalSet/maternalunder.txt";
		// List<String> MeetaDown = new LinkedList<String>( new HashSet<String>(
		// FileTools.getLines( file ) ) );
		// loader.writeAgeMatrices( MeetaDown, "maternalunder" );
		//
		// file = "/home/leon/MaternalSet/maternalover.txt";
		// MeetaDown = new LinkedList<String>( new HashSet<String>(
		// FileTools.getLines( file ) ) );
		// loader.writeAgeMatrices( MeetaDown, "maternalover" );

		// getRandomSubsetGenes( int size )
		List<String> MeetaDown = new LinkedList<String>(new HashSet<String>(
				loader.getRandomSubsetGenes(331)));
		loader.writeAgeMatrices(MeetaDown, "random.331");

		MeetaDown = new LinkedList<String>(new HashSet<String>(
				loader.getRandomSubsetGenes(162)));
		loader.writeAgeMatrices(MeetaDown, "random.162");

		System.exit(1);

		List<String> NEgeneSet = FileTools
				.getLines("/home/leon/Downloads/NEHuman.txt");
		loader.writeAgeMatrices(NEgeneSet, "NE");
		loader.writeRegionMatrix(NEgeneSet, "NE", "0 mo");
		loader.writeRegionMatrix(NEgeneSet, "NE", "48 mo");

		List<String> OEgeneSet = FileTools
				.getLines("/home/leon/Downloads/OEHuman.txt");
		loader.writeAgeMatrices(OEgeneSet, "OE");
		loader.writeRegionMatrix(OEgeneSet, "OE", "0 mo");
		loader.writeRegionMatrix(OEgeneSet, "OE", "48 mo");

	}

	public static void createCorrelationMatrix(BrainSpanGeneric loader,
			DoubleMatrix<String, String> snca, Set<String> interferonGOgroup,
			String name) throws Exception {
		DoubleMatrix<String, String> interferon = new DenseDoubleMatrix<String, String>(
				interferonGOgroup.size(), snca.columns());
		interferon.setColumnNames(snca.getColNames());
		interferon.setRowNames(new LinkedList<String>(interferonGOgroup));
		for (String testGene : interferonGOgroup) {
			DoubleMatrix<String, String> testGeneData = loader.getRegionByAge(
					testGene, true);
			DoubleMatrix<String, String> corResult = loader.correlateViaTime(
					snca, testGeneData);

			for (String colName : corResult.getColNames()) {
				interferon.setByKeys(testGene, colName,
						corResult.getByKeys("Spearman", colName));
			}

		}
		String filename = "/home/abams/Desktop/human/BrainSpan/SNCA/correlations/"
				+ name + ".tsv";
		Util.writeRTable(filename, interferon);
	}

}
