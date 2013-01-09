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
package ubic.BAMSandAllen.human;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.ABAMSDataMatrix;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.Homologene.HomoloGeneLoader;
import ubic.BAMSandAllen.adjacency.CorrelationAdjacency;
import ubic.basecode.dataStructure.CountingMap;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.io.reader.DoubleMatrixReader;
import ubic.basecode.math.CorrelationStats;
import ubic.basecode.math.DescriptiveWithMissing;
import ubic.basecode.util.FileTools;
import au.com.bytecode.opencsv.CSVReader;
import cern.colt.list.DoubleArrayList;

public class LoadHumanCSV {
    private static Log log = LogFactory.getLog( LoadHumanCSV.class.getName() );
    ABAMSDataMatrix matrix;
    String ID;

    public LoadHumanCSV( String ID ) throws Exception {
        this( ID, false );
    }

    public LoadHumanCSV( String ID, boolean commonRegions ) throws Exception {
        this.ID = ID;
        try {
            DoubleMatrixReader reader = new DoubleMatrixReader();
//            reader.setTopLeft( false );

            String filename = "";
            if ( ID.equals( "9861" ) )
                filename = "C:\\Users\\leon\\Desktop\\Human array\\" + ID + "\\matrix.29192 x 323.txt";
            if ( ID.equals( "10021" ) )
                filename = "C:\\Users\\leon\\Desktop\\Human array\\" + ID + "\\matrix.29192 x 346.txt";

            DoubleMatrix<String, String> inputMatrix = reader.read( filename );
            matrix = new ABAMSDataMatrix( inputMatrix, "HumanArray", new CorrelationAdjacency( inputMatrix ) );

            LinkedList<String> newRowNames = new LinkedList<String>();
            for ( String row : matrix.getRowNames() ) {
                newRowNames.addLast( row.substring( 1, row.length() - 1 ) );
            }

            LinkedList<String> newColNames = new LinkedList<String>();
            for ( String row : matrix.getColNames() ) {
                newColNames.addLast( row.substring( 1, row.length() - 1 ) );
            }

            matrix.setColumnNames( newColNames );

            matrix.setRowNames( newRowNames );

            if ( commonRegions ) {
                matrix = matrix.retainColumns( FileTools
                        .getLines( "C:\\Users\\leon\\Desktop\\Human array\\commonRegions.txt" ) );
            }
        } catch ( Exception e ) {
            e.printStackTrace();
            System.exit( 1 );
            makeMatrix();
        }
    }

    public void setToRows( Set<String> rowsToKeep ) {
        matrix = slimToRows( rowsToKeep );
    }

    public ABAMSDataMatrix slimToRows( Set<String> rowsToKeep ) {
        // log.info( rowsToKeep );
        // log.info( "Rows given:" + rowsToKeep.size() );
        // log.info( "All rows:" + matrix.rows() );
        Set<String> toRemove = new HashSet<String>( matrix.getRowNames() );
        toRemove.removeAll( rowsToKeep );
        // log.info( "Rows to remove:" + toRemove.size() );
        return matrix.removeRows( toRemove );
    }

    public Map<String, String> getAcroToNameMap() throws Exception {
        Map<String, String> acroToName = new HashMap<String, String>();

        // FileReader fr = new FileReader( "C:\\Users\\leon\\Desktop\\Human array\\" + ID + "\\SampleAnnot.csv" );
        FileReader fr = new FileReader( "/home/leon/Desktop/Windows Desktop/Stuff from Japan/Human array/" + ID
                + "/SampleAnnot.csv" );

        CSVReader reader = new CSVReader( fr );
        // eat header
        reader.readNext();
        String line[];
        CountingMap<String> samplesPerRegion = new CountingMap<String>();
        while ( ( line = reader.readNext() ) != null ) {
            acroToName.put( line[4], line[5] );
            log.info( line[5] );
        }
        fr.close();

        return acroToName;
    }

    public void makeMatrix() throws Exception {
        LinkedList<String> regionNames = new LinkedList<String>();
        FileReader fr = new FileReader( "C:\\Users\\leon\\Desktop\\Human array\\" + ID + "\\SampleAnnot.csv" );
        CSVReader reader = new CSVReader( fr );
        // eat header
        reader.readNext();
        String line[];
        CountingMap<String> samplesPerRegion = new CountingMap<String>();
        while ( ( line = reader.readNext() ) != null ) {
            log.info( line[5] );
            regionNames.addLast( line[5] );
            samplesPerRegion.increment( line[5] );
        }

        fr.close();
        log.info( regionNames.size() );
        log.info( ( new HashSet<String>( regionNames ) ).size() );

        // load in probe to gene
        // CountingMap<Integer> probesPerGene= new CountingMap<Integer>();
        Map<String, String> probesToGene = new HashMap<String, String>();
        fr = new FileReader( "C:\\Users\\leon\\Desktop\\Human array\\" + ID + "\\Probes.csv" );
        reader = new CSVReader( fr );
        while ( ( line = reader.readNext() ) != null ) {
            String probe = line[0];
            String gene = line[3];
            probesToGene.put( probe, gene );
        }
        fr.close();

        LinkedList<String> geneNames = new LinkedList<String>( new HashSet<String>( probesToGene.values() ) );
        LinkedList<String> regionNamesInMatrix = new LinkedList<String>( new HashSet<String>( regionNames ) );
        DoubleMatrix<String, String> matrix = new DenseDoubleMatrix<String, String>( geneNames.size(),
                regionNamesInMatrix.size() );
        matrix.setRowNames( geneNames );
        matrix.setColumnNames( regionNamesInMatrix );

        CountingMap<String> probesUsed = new CountingMap<String>();

        fr = new FileReader( "C:\\Users\\leon\\Desktop\\Human array\\" + ID + "\\MicroarrayExpression.csv" );
        reader = new CSVReader( fr );
        int lineCount = 0;
        while ( ( line = reader.readNext() ) != null ) {
            if ( lineCount++ % 1000 == 0 ) log.info( lineCount );
            String probe = line[0];
            String gene = probesToGene.get( probe );
            probesUsed.increment( gene );
            for ( int i = 1; i < line.length; i++ ) {
                String regionName = regionNames.get( i - 1 );
                double currentValue = matrix.getByKeys( gene, regionName );
                matrix.setByKeys( gene, regionName, currentValue + Double.parseDouble( line[i] ) );
            }
        }
        fr.close();

        // divide by number of probes
        for ( String row : matrix.getRowNames() ) {
            int probes = probesUsed.get( row );
            for ( String region : matrix.getColNames() ) {
                double current = matrix.getByKeys( row, region );
                int samples = samplesPerRegion.get( region );
                current = current / ( double ) probes;
                current = current / ( double ) samples;
                matrix.setByKeys( row, region, current );
            }
        }

        log.info( matrix.rows() + " x " + matrix.columns() );
        Util.writeRTable( "C:\\Users\\leon\\Desktop\\Human array\\" + ID + "\\matrix." + matrix.rows() + " x "
                + matrix.columns() + ".txt", matrix );
    }

    /**
     * @param args
     */

    public ABAMSDataMatrix getRandomSubMatrix( int size, Random r ) {
        List<String> rownames = new LinkedList<String>( matrix.getRowNames() );
        Collections.shuffle( rownames, r );
        return slimToRows( new HashSet<String>( rownames.subList( 0, size ) ) );
    }

    public void test( int size1, int size2, double statistic, boolean coronal ) throws Exception {
        Random r = new Random( 1 );
        if ( coronal ) {
            reduceToCoronal();
        }
        int lesserHits = 0;
        DoubleArrayList sampleHistory = new DoubleArrayList();

        int greaterHits = 0;
        for ( int i = 0; i < 10000; i++ ) {
            ABAMSDataMatrix NEmatrix = getRandomSubMatrix( size1, r );
            ABAMSDataMatrix OEmatrix = getRandomSubMatrix( size2, r );

            boolean removeNan = true;
            double[] NEarray = Util.columnSums( NEmatrix, removeNan ).getRow( 0 );
            double[] OEarray = Util.columnSums( OEmatrix, removeNan ).getRow( 0 );

            // log.info( "Pearson:" + CorrelationStats.correl( NEarray, OEarray ) );
            double sampleStatistic = Util.spearmanCorrel( NEarray, OEarray );

            sampleHistory.add( sampleStatistic );

            if ( sampleStatistic < statistic ) {
                lesserHits++;
                log.info( "Sample " + "Spearman correlation" + ":" + sampleStatistic + " less than " + statistic );
            }
            if ( sampleStatistic > statistic ) {
                greaterHits++;
                // log.info( "Sample " + "Spearman correlation" + ":" + sampleStatistic + " greater than " + statistic
                // );
            }
        }
        double mean = DescriptiveWithMissing.mean( sampleHistory );
        String statisticString = "Spearman";
        log.info( "Sample " + statisticString + " greater than " + statistic + ", " + greaterHits + " times" );
        log.info( "Sample " + statisticString + " less than " + statistic + ", " + lesserHits + " times" );
        log.info( "average resample" + statisticString + " " + mean + "" );
        log.info( "standard deviation resample" + statisticString + " "
                + Math.sqrt( DescriptiveWithMissing.sampleVariance( sampleHistory, mean ) ) );

    }

    private void reduceToCoronal() throws IOException, Exception {
        String base = "C:\\Users\\leon\\Desktop\\Human array\\";

        List<String> coronalIDs = FileTools.getLines( base + "coronal.txt" );

        HomoloGeneLoader loader = new HomoloGeneLoader();
        log.info( coronalIDs.size() );

        Set<String> coronalHuman = loader.getHumanIDsFromMouseID( coronalIDs );

        log.info( coronalHuman.size() );
        setToRows( coronalHuman );
        log.info( "Moving to coronal genes only:" + matrix.rows() );
    }

    public static void missingRegions() throws Exception {
        String ID = "10021";
        Set<String> h10021 = getRegionsUsed( ID );
        ID = "9861";
        Set<String> h9861 = getRegionsUsed( ID );

        log.info( "10021 size:" + h10021.size() );
        log.info( "9861 size:" + h9861.size() );

        log.info( "98 intersect 100: " + Util.intersectSize( h9861, h10021 ) );
        log.info( "98 intersect 100: " + Util.intersect( h9861, h10021 ) );
        log.info( "98 - 100: " + Util.subtract( h9861, h10021 ).size() );
        log.info( "100 - 98: " + Util.subtract( h10021, h9861 ).size() );
        log.info( "100 + 98: " + Util.union( h10021, h9861 ).size() );

        log.info( "98 - 100: " + Util.subtract( h9861, h10021 ) );
        log.info( "100 - 98: " + Util.subtract( h10021, h9861 ) );

        String base = "C:\\Users\\leon\\Desktop\\Human array\\commonRegions.txt";
        FileTools.stringsToFile( ( Set<String> ) Util.intersect( h10021, h9861 ), base );

        base = "C:\\Users\\leon\\Desktop\\Human array\\unionRegions.txt";
        FileTools.stringsToFile( ( Set<String> ) Util.union( h10021, h9861 ), base );
    }

    public static void forJoao( String ID ) throws Exception {
        String base = "C:\\Users\\leon\\Desktop\\Human array\\";
        List<String> JoaoIDs = FileTools.getLines( base + "Joao.txt" );
        // String ID = "10021";
        // ID = "9861";
        LoadHumanCSV human = new LoadHumanCSV( ID );

        ABAMSDataMatrix toWrite = human.slimToRows( new HashSet<String>( JoaoIDs ) );
        Util.writeImage( base + "Joao." + ID + ".png", toWrite );
        Util.writeRTable( base + "Joao.R." + ID + ".txt", toWrite );

        Random r = new Random( 1 );
        toWrite = human.getRandomSubMatrix( JoaoIDs.size(), r );
        Util.writeRTable( base + "Joao.R.random." + ID + ".txt", toWrite );

        // get random

    }

    public ABAMSDataMatrix getMatrix() {
        return matrix;
    }



    public static void main( String[] args ) throws Exception {
        String base = "C:\\Users\\leon\\Desktop\\Human array\\";
        // TODO Auto-generated method stub
        // human.makeMatrix();
        // forJoao( "10021" );
        // forJoao( "9861" );
        // missingRegions();
        // System.exit( 1 );

        List<String> oeIDs = FileTools.getLines( base + "oe.txt" );
        List<String> neIDs = FileTools.getLines( base + "ne.txt" );
        // List<String> coronalIDs = FileTools.getLines( base + "coronal.txt" );

        HomoloGeneLoader loader = new HomoloGeneLoader();

        Set<String> OE = loader.getHumanIDsFromMouseID( oeIDs );
        Set<String> NE = loader.getHumanIDsFromMouseID( neIDs );
        FileTools.stringsToFile( OE, base + "OEHuman.txt" );
        FileTools.stringsToFile( NE, base + "NEHuman.txt" );
        Set<String> both = ( Set<String> ) Util.union( NE, OE );

        log.info( NE );
        log.info( OE );
        String ID = "10021";
        ID = "9861";
        boolean commonRegions = true;
        LoadHumanCSV human = new LoadHumanCSV( ID, commonRegions );
        ABAMSDataMatrix NEmatrix = human.slimToRows( NE );

        ABAMSDataMatrix OEmatrix = human.slimToRows( OE );

        boolean removeNan = true;
        double[] NEarray = Util.columnSums( NEmatrix, removeNan ).getRow( 0 );
        double[] OEarray = Util.columnSums( OEmatrix, removeNan ).getRow( 0 );

        log.info( "Pearson:" + CorrelationStats.correl( NEarray, OEarray ) );
        log.info( "Spearman:" + Util.spearmanCorrel( NEarray, OEarray ) );
        boolean coronal = true;
        human.test( NEmatrix.rows(), OEmatrix.rows(), Util.spearmanCorrel( NEarray, OEarray ), coronal );
        log.info( "regions used:" + NEarray.length );

        ABAMSDataMatrix bothmatrix = human.slimToRows( both );

        // LinkedList<String> newNames = new LinkedList<String>();
        // for ( String geneName : bothmatrix.getRowNames() ) {
        // if ( NE.contains( geneName ) )
        // newNames.addLast( geneName + "**NE**" );
        // else
        // newNames.addLast( geneName + "||OE||" );
        // }
        // bothmatrix.setRowNames( newNames );
        // Util.writeRTable( base + "bothMatrix.txt", bothmatrix );
    }

    public Set<String> getRegions() {
        return new HashSet<String>( matrix.getColNames() );
    }

    public static Set<String> getRegionsUsed( String ID ) throws Exception {
        LoadHumanCSV human = new LoadHumanCSV( ID );
        return human.getRegions();
    }
}
