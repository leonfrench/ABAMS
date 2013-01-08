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

import ubic.BAMSandAllen.RankedGeneListLoader;
import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.AllenDataLoaders.ImageSeriesInfoLoader;
import ubic.BAMSandAllen.Homologene.HomoloGeneLoader;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.math.CorrelationStats;
import ubic.basecode.math.DescriptiveWithMissing;
import ubic.basecode.math.RandomChooser;
import ubic.basecode.util.FileTools;
import au.com.bytecode.opencsv.CSVReader;
import cern.colt.list.DoubleArrayList;

public class BrainSpanGeneric {
    private static Log log = LogFactory.getLog( BrainSpanGeneric.class.getName() );

    String baseFolder = "/home/leon/Downloads/BrainSpan/";

    DoubleMatrix<Map<String, String>, Map<String, String>> data;
    String type;

    public Set<Map<String, String>> getRowBySymbol( String geneSymbol ) {
        Set<Map<String, String>> result = new HashSet<Map<String, String>>();
        for ( Map<String, String> row : data.getRowNames() ) {
            if ( row.get( "gene_symbol" ).equals( geneSymbol ) ) {
                result.add( row );
            }
        }
        return result;
    }

    public BrainSpanGeneric( String type ) throws Exception {
        this.type = type;
        baseFolder += type + "/";

        // read in columns
        List<Map<String, String>> colMaps = getLines( baseFolder + type + "_columns_metadata.csv" );

        // read in rows
        List<Map<String, String>> rowMaps = getLines( baseFolder + type + "_rows_metadata.csv" );

        log.info( "rowMaps.size(), colMaps.size() " + rowMaps.size() + "x" + colMaps.size() );

        data = new DenseDoubleMatrix<Map<String, String>, Map<String, String>>( rowMaps.size(), colMaps.size() );

        data.setColumnNames( colMaps );

        log.info( "setting row maps" );
        data.setRowNames( rowMaps );
        log.info( "done rownames" );

        String dataFile = baseFolder + type + "_matrix.csv";
        CSVReader reader = new CSVReader( new FileReader( dataFile ) );
        String line[];
        int row = 0;
        while ( ( line = reader.readNext() ) != null ) {
            // log.info( line.length );
            if ( row % 10000 == 0 ) {
                log.info( "Reading row:" + row );
                // break;
            }
            // shift over by one because the first column is the Ensembl ID (probably should check that against the row
            // data
            for ( int i = 1; i < line.length; i++ ) {
                data.set( row, i - 1, Double.parseDouble( line[i] ) );
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

    public List<Map<String, String>> getLines( String filename ) throws Exception {
        log.info( "Reading:" + filename );
        CSVReader reader = new CSVReader( new FileReader( filename ) );
        List<Map<String, String>> result = new LinkedList<Map<String, String>>();
        List<String[]> lines = ( List<String[]> ) reader.readAll();
        // grab the header line
        String[] keys = lines.get( 0 );
        lines.remove( 0 );
        for ( String key : keys )
            log.info( key );
        for ( String[] line : lines ) {
            Map<String, String> lineMap = new HashMap<String, String>();
            for ( int i = 0; i < line.length; i++ ) {
                lineMap.put( keys[i], line[i] );
            }
            result.add( lineMap );
        }
        reader.close();

        return result;
    }

    public Set<String> getUniqueGenes() {
        return getUniqueValues( data.getRowNames(), "gene_symbol" );
    }

    public Set<String> getUniqueRegions() {
        return getUniqueValues( data.getColNames(), "structure_name" );
    }

    public Set<String> getUniqueValues( List<Map<String, String>> maps, String key ) {
        Set<String> result = new HashSet<String>();
        for ( Map<String, String> map : maps ) {
            result.add( map.get( key ) );
        }
        return result;
    }

    public void writeAgeMatrices( Collection<String> genes, String name ) throws Exception {
        writeAgeMatrices( genes, name, null );
    }

    public void writeAgeMatrices( Collection<String> genes, String name, String region ) throws Exception {
        // how many ages?
        Set<String> ageValuesSet = getUniqueValues( data.getColNames(), "age" );

        List<String> ageValues = new LinkedList<String>( ageValuesSet );
        Collections.sort( ageValues, new TimeSorter() );

        DoubleMatrix<String, String> result = new DenseDoubleMatrix<String, String>( genes.size(), ageValues.size() );
        result.setRowNames( new LinkedList<String>( genes ) );
        result.setColumnNames( ageValues );

        DoubleMatrix<String, String> resultSTD = new DenseDoubleMatrix<String, String>( genes.size(), ageValues.size() );
        resultSTD.setRowNames( new LinkedList<String>( genes ) );
        resultSTD.setColumnNames( ageValues );

        Set<String> uniqueGenes = getUniqueGenes();
        log.info( "Genes not in set:" + Util.subtract( genes, uniqueGenes ).size() );
        genes.retainAll( uniqueGenes );

        // for a given age and a given gene, we have a list - search all?
        for ( String gene : genes ) {
            for ( String age : ageValues ) {
                List<Double> values = findData( gene, age, region );
                double[] valueArray = ArrayUtils.toPrimitive( values.toArray( new Double[0] ) );
                DoubleArrayList expValuesDAL = new DoubleArrayList( valueArray );

                double mean = DescriptiveWithMissing.mean( expValuesDAL );
                double sampleStandardDeviation = Math
                        .sqrt( DescriptiveWithMissing.sampleVariance( expValuesDAL, mean ) );

                result.setByKeys( gene, age, mean );
                resultSTD.setByKeys( gene, age, sampleStandardDeviation );
            }
        }

        String filename = baseFolder + name + ".size." + genes.size() + ".txt";
        Util.writeRTable( filename, result );
        Util.writeRTable( baseFolder + name + "STD.size." + genes.size() + ".txt", resultSTD );
        log.info( "Files written:" + filename );
    }

    public List<Double> findData( String gene, String age, String region ) {
        List<Double> result = new LinkedList<Double>();
        for ( Map<String, String> col : data.getColNames() ) {
            if ( col.get( "age" ).equals( age ) ) {
                if ( region == null || col.get( "structure_name" ).equals( region ) ) {
                    for ( Map<String, String> row : data.getRowNames() ) {
                        if ( row.get( "gene_symbol" ).equals( gene ) ) {
                            result.add( data.getByKeys( row, col ) );
                        }
                    }
                }
            }
        }
        return result;
    }

    public DoubleMatrix<String, String> dump() {
        DoubleMatrix<String, String> result = new DenseDoubleMatrix<String, String>( data.asArray() );
        List<String> colKeyOrder = new LinkedList<String>( data.getColName( 0 ).keySet() );
        List<String> rowKeyOrder = new LinkedList<String>( data.getRowName( 0 ).keySet() );

        for ( Map<String, String> columnMap : data.getColNames() ) {
            String columnString = "";
            for ( String key : colKeyOrder ) {
                columnString += "," + columnMap.get( key );
            }
            result.addColumnName( columnString.substring( 1 ) );
        }
        for ( Map<String, String> rowMap : data.getRowNames() ) {
            String rowString = "";
            for ( String key : rowKeyOrder ) {
                rowString += "," + rowMap.get( key );
            }
            result.addRowName( rowString.substring( 1 ) );
        }

        log.info( "Column key:" + Arrays.toString( colKeyOrder.toArray() ) );
        log.info( "Row key:" + Arrays.toString( rowKeyOrder.toArray() ) );
        return result;
    }

    public DoubleMatrix<String, String> getGeneData( String gene ) {
        // make a new matrix?
        DoubleMatrix<String, String> result;

        List<Map<String, String>> rowsWithGene = new LinkedList<Map<String, String>>();
        for ( Map<String, String> row : data.getRowNames() ) {
            if ( row.get( "gene_symbol" ).equals( gene ) ) {
                rowsWithGene.add( row );
            }
        }
        log.info( "Rows/exons for gene:" + rowsWithGene.size() );
        result = new DenseDoubleMatrix<String, String>( rowsWithGene.size(), data.columns() );

        for ( Map<String, String> row : rowsWithGene ) {
            String exonInfo = "";
            if ( type.startsWith( "exon" ) ) {
                exonInfo = "exon." + row.get( "start" ) + ".to." + row.get( "end" ) + ".";
            }

            String rowName = exonInfo + row.get( "gene_symbol" );

            result.addRowName( rowName );
            for ( int col = 0; col < data.columns(); col++ ) {
                result.set( result.getRowIndexByName( rowName ), col, data.get( data.getRowIndexByName( row ), col ) );
            }
        }

        // colapse column data
        for ( Map<String, String> col : data.getColNames() ) {
            String colName = col.get( "age" ) + "." + col.get( "gender" ) + "." + col.get( "structure_id" ) + "."
                    + col.get( "donor_name" );
            result.addColumnName( colName );
        }
        return result;
    }

    /**
     * Simple correlation between rows
     * 
     * @param gene
     * @throws Exception
     */
    public DoubleMatrix<String, String> correlate( String gene ) throws Exception {
        Set<Map<String, String>> geneRows = getRowBySymbol( gene );
        log.info( "Number of rows:" + geneRows.size() );
        log.info( "Warning: Using only one" );
        Map<String, String> geneRow = geneRows.iterator().next();
        double[] dataRow = data.getRowByName( geneRow );
        return correlate( dataRow, gene );

    }

    public DoubleMatrix<String, String> correlateToTime() throws Exception {
        double[] dataRow = new double[data.columns()];
        // Vector for time
        for ( int i = 0; i < dataRow.length; i++ ) {
            dataRow[i] = i;
        }
        return correlate( dataRow, "Developmental Intervals" );
    }

    public DoubleMatrix<String, String> correlate( double[] dataRow, String name ) throws Exception {

        DoubleMatrix<String, String> correlations = new DenseDoubleMatrix<String, String>( getUniqueGenes().size(), 4 );
        correlations.setRowNames( new LinkedList<String>( getUniqueGenes() ) );
        correlations.addColumnName( "Spearman " + name );
        correlations.addColumnName( "Pearson " + name );
        correlations.addColumnName( "Spearman P " + name );
        correlations.addColumnName( "Pearson P " + name );

        for ( Map<String, String> rowToTest : data.getRowNames() ) {
            double[] testRow = data.getRowByName( rowToTest );
            double spearmanCor = Util.spearmanCorrel( dataRow, testRow );
            double spearmanP = CorrelationStats.spearmanPvalue( spearmanCor, testRow.length ) * data.rows();
            double pearsonCor = CorrelationStats.correl( dataRow, testRow );
            double pearsonP = CorrelationStats.pvalue( pearsonCor, testRow.length ) * data.rows();
            String geneTested = rowToTest.get( "gene_symbol" );
            correlations.setByKeys( geneTested, "Spearman " + name, spearmanCor );
            correlations.setByKeys( geneTested, "Pearson " + name, pearsonCor );
            correlations.setByKeys( geneTested, "Spearman P " + name, spearmanP );
            correlations.setByKeys( geneTested, "Pearson P " + name, pearsonP );
        }
        return correlations;
    }

    public Set<String> getRandomCoronalSubset( int size ) throws Exception {
        ImageSeriesInfoLoader infoLoader = new ImageSeriesInfoLoader();
        HomoloGeneLoader homoLoader = new HomoloGeneLoader();
        Set<String> result = infoLoader.getAllRowsWithCoronalGenes();
        RankedGeneListLoader loader = new RankedGeneListLoader( new LinkedList<String>( result ), "x" );
        loader.convertToNCBIIDs();
        Set<String> humanSymbols = homoLoader.getHumanSymbolFromMouseID( loader.getLines() );
        humanSymbols = RandomChooser.chooseRandomSubset( size, humanSymbols );
        return humanSymbols;
    }

    public Set<String> getRandomSubsetGenes( int size ) {
        Set<String> result = new HashSet<String>();
        RandomChooser.init( 1 );

        for ( Map<String, String> row : RandomChooser.chooseRandomSubset( size, data.getRowNames() ) ) {
            result.add( row.get( "gene_symbol" ) );
        }
        return result;
    }

    public void makeAgeMaricesForPLoSSets() throws Exception {
        String base = SetupParameters.config.getString( "abams.dataFolder" ) + "rankedGenes/near final ammon/";
        String outFile = base + "LOOGenesInOrder.out.partialcon.ammon.txt.374.0.016424.topGenes.txt";
        String inFile = base + "LOOGenesInOrder.in.partialcon.ammon.txt.452.0.01985.topGenes.txt";
        String proxFile = base + "LOOGenesInOrder.space.ammon.txt.420.0.018435.topGenes.txt";

        writeMatrixForMouseList( outFile, "out.partial" );
        writeMatrixForMouseList( inFile, "in.partial" );
        writeMatrixForMouseList( proxFile, "proximity" );

    }

    public void writeMatrixForMouseList( String rankedFile, String name ) throws Exception {
        RankedGeneListLoader rankedLoader = new RankedGeneListLoader( rankedFile );
        rankedLoader.convertToNCBIIDs();
        List<String> mouseIDs = rankedLoader.getLines();
        HomoloGeneLoader homoLoader = new HomoloGeneLoader();
        Set<String> humanSymbols = homoLoader.getHumanSymbolFromMouseID( mouseIDs );
        log.info( "Human genes:" + humanSymbols.size() );
        writeAgeMatrices( humanSymbols, name );
    }

    public void writeHumanSets( String filename ) throws Exception {
        // List<String> rankedFile = FileTools.getLines( filename );
        RankedGeneListLoader rankedLoader = new RankedGeneListLoader( filename );
        rankedLoader.convertToNCBIIDs();
        List<String> mouseIDs = rankedLoader.getLines();

        HomoloGeneLoader homoLoader = new HomoloGeneLoader();
        File f = new File( filename );
        // f.getName();
        Set<String> humanSymbols = homoLoader.getHumanSymbolFromMouseID( mouseIDs );
        FileTools.stringsToFile( humanSymbols, "/home/leon/Downloads/BrainSpan/human gene lists/" + f.getName() );
    }

    public static void main( String[] args ) throws Exception {
        BrainSpanGeneric loader = new BrainSpanGeneric( "genes" );

        String base = SetupParameters.config.getString( "abams.dataFolder" ) + "rankedGenes/near final ammon/";

        // loader.writeHumanSets( base + "LOOGenesInOrder.out.partialcon.ammon.txt.374.0.016424.topGenes.txt" );
        // loader.writeHumanSets( base + "LOOGenesInOrder.in.partialcon.ammon.txt.452.0.01985.topGenes.txt" );
        // System.exit( 1 );

        String geneCor = "SLC20A2";
        // Util.writeRTable( "/home/leon/Downloads/SLCCorrelations.txt", loader.correlate( geneCor ) );
        Util.writeRTable( "/home/leon/Downloads/TimeCorrelations.txt", loader.correlateToTime() );
        System.exit( 1 );

        // Set<String> SZUp = new HashSet<String>( FileTools.getLines( "/home/leon/Downloads/query_szUp.txt" ) );
        // loader.writeAgeMatrices( SZUp, "SZUp", "Dorsolateral prefrontal cortex" );
        //
        // Set<String> SZDown = new HashSet<String>( FileTools.getLines( "/home/leon/Downloads/query_szDown.txt" ) );
        // loader.writeAgeMatrices( SZDown, "SZDown", "Dorsolateral prefrontal cortex" );
        // loader.writeAgeMatrices( loader.getRandomCoronalSubset( 300 ), "RandomCoronal" );

        Set<String> axonGuidance = new HashSet<String>( FileTools
                .getLines( "/home/leon/Downloads/axon guidance genes.txt" ) );
        loader.writeAgeMatrices( axonGuidance, "axonGuidance" );
        System.exit( 1 );

        for ( String region : loader.getUniqueRegions() ) {
            System.out.println( region );
        }

        loader.makeAgeMaricesForPLoSSets();

        System.exit( 1 );

        // // TODO Auto-generated method stub
        // BrainSpanLoader loader = new BrainSpanLoader();
        // String gene = "SLC20A2";
        // DoubleMatrix<String, String> geneMatrix = loader.getGeneData( gene );
        // Util.writeRTable( "/home/leon/Downloads/BrainSpan/" + gene + ".txt", geneMatrix );

        // String gene = "SLC20A2";
        // DoubleMatrix<String, String> geneMatrix = loader.getGeneData( gene );
        // Util.writeRTable( "/home/leon/Downloads/BrainSpan/" + gene + "." + loader.getType() + ".txt", geneMatrix );

        Set<String> geneSet = new HashSet<String>();
        geneSet.add( "SLC20A2" );

        for ( String gene : loader.getUniqueGenes() ) {
            if ( gene.startsWith( "SEMA" ) ) geneSet.add( gene );
        }

        loader.writeAgeMatrices( geneSet, "Sema" );

        List<String> NEgeneSet = FileTools.getLines( "/home/leon/Downloads/NEHuman.txt" );
        loader.writeAgeMatrices( NEgeneSet, "NE" );

        List<String> OEgeneSet = FileTools.getLines( "/home/leon/Downloads/OEHuman.txt" );
        loader.writeAgeMatrices( OEgeneSet, "OE" );

    }

}
