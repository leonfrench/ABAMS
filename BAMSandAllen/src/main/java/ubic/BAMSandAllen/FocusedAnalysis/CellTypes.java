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
package ubic.BAMSandAllen.FocusedAnalysis;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes;
import ubic.BAMSandAllen.RankedGeneListLoader;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.MatrixPairs.AllenMatrixPair;
import ubic.BAMSandAllen.MatrixPairs.AllenMatrixPairFactory;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenExpressionMatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ExpressionMatrixPairFactory;
import ubic.BAMSandAllen.MatrixPairs.MatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenPartialExpressionMatrixPair.RegressMatrix;
import ubic.BAMSandAllen.geneFilters.PlaneRemoveFilter;
import ubic.BAMSandAllen.geneFilters.PlaneRemoveFilter.Plane;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.dataStructure.params.ParamKeeper;
import ubic.basecode.graphics.ColorMatrix;
import ubic.basecode.graphics.MatrixDisplay;
import ubic.basecode.util.FileTools;

public class CellTypes {
    private static Log log = LogFactory.getLog( CellTypes.class.getName() );
    public String outputDir = "/grp/java/workspace/BAMSandAllen/data/ABACellEnriched/";
    MatrixPair pair;
    Map<File, List<String>> files;
    Direction direction;
    boolean degree;

    public CellTypes( MatrixPair pair, Direction direction ) {
        this.pair = pair;
        this.direction = direction;
        files = new HashMap<File, List<String>>();
        // files.put( new File( "All rows" ), pair.getMatrixBDataRows() );
    }

    public void createVennFile() throws Exception {
        File genSetRowsFile = new File( outputDir + "GeneSetVennMasterRows." + getFileBaseName() + ".txt" );
        File genSetGenesFile = new File( outputDir + "GeneSetVennMasterGenes." + getFileBaseName() + ".txt" );
        genSetGenesFile.delete();
        List<File> fileList = new LinkedList<File>( files.keySet() );
        Collections.sort( fileList );
        for ( File file : fileList ) {
            List<String> rows = files.get( file );
            if ( file.getName().equals( "All rows" ) ) continue;
            Util.addToVennMasterFile( genSetRowsFile, rows, file.getName() );
            Util.addToVennMasterFile( genSetGenesFile, Util.getUniqueGenes( rows ), file.getName() );
        }
    }

    public static void testXcorrelations( int testIterations, Direction direction ) throws Exception {
        boolean logDistance = true;
        boolean removeNonExp = true;
        boolean degree = true;

        AllenMatrixPair pair = AllenMatrixPairFactory.getSpaceAndExpressionPair( direction, removeNonExp, logDistance );
        pair.applyMatrixBRowFilter( new PlaneRemoveFilter( Plane.SAGITTAL ) );
        Set<String> removeRows = new HashSet<String>();
        removeRows.add( "y" );
        removeRows.add( "z" );
        pair.removeMatrixADataRows( removeRows );
        pair.printDimensions();
        log.info( pair.getFlattenedCorrelation( false ) );
        log.info( "Spearman:" + pair.getFlattenedCorrelation( true ) );
        CellTypes celltypes = new CellTypes( pair, direction );
        celltypes.addFrontiersList();
        // celltypes.addOptimizedSets();
        celltypes.run( testIterations, degree );
    }

    public static void testVolumecorrelations( int testIterations, Direction direction ) throws Exception {
        boolean removeNonExp = true;
        boolean degree = true;

        AllenMatrixPair pair = AllenMatrixPairFactory.getVolumeAndExpressionPair( direction, removeNonExp );

        // AllenMatrixPair pair = AllenMatrixPairFactory.getSpaceAndExpressionPair( direction, removeNonExp, logDistance
        // );
        pair.applyMatrixBRowFilter( new PlaneRemoveFilter( Plane.SAGITTAL ) );

        pair.printDimensions();
        log.info( "Pearson:" + pair.getFlattenedCorrelation( false ) );
        log.info( "Spearman:" + pair.getFlattenedCorrelation( true ) );

        CellTypes celltypes = new CellTypes( pair, direction );
        celltypes.addFrontiersList();
        // celltypes.addOptimizedSets();
        celltypes.run( testIterations, degree );
    }

    private void addFrontiersList() throws Exception {
        File folder;
        // folder = new File( "/grp/java/workspace/BAMSandAllen/data/ABACellEnriched/" );
        // for ( File file : folder.listFiles() ) {
        // if ( file.getAbsoluteFile().toString().endsWith( ".txt" ) && !file.getName().startsWith( "GeneSet" ) ) {
        // addGeneSet( file );
        // }
        // }

        folder = new File( "/grp/java/workspace/BAMSandAllen/data/ABACellEnriched/patterns/" );
        for ( File file : folder.listFiles() ) {
            if ( file.getAbsoluteFile().toString().endsWith( ".txt" ) && !file.getName().startsWith( "GeneSet" ) ) {
                addRowSet( file );
            }
        }

    }

    private void addMainLists() throws Exception {

        // G2Cdbgenelists
        File folder = new File( "/grp/java/workspace/BAMSandAllen/data/G2Cdbgenelists/" );
        for ( File file : folder.listFiles() ) {
            addGeneSet( file );
        }

        folder = new File( "/grp/java/workspace/BAMSandAllen/data/ABACellEnriched/everything/" );
        for ( File file : folder.listFiles() ) {
            if ( file.getAbsoluteFile().toString().endsWith( ".txt" ) && !file.getName().startsWith( "GeneSet" ) ) {
                addGeneSet( file );
            }
        }

        folder = new File( "/home/leon/Desktop/disease genes/" );
        for ( File file : folder.listFiles() ) {
            addGeneSet( file );
        }
    }

    public void addOptimizedSets() throws Exception {
        String baseFileName = "/grp/java/workspace/BAMSandAllen/data/rankedGenes/near final ammon/";
        addRowSet( new File( baseFileName + "LOOGenesInOrder.in.partialcon.ammon.txt.452.0.01985.topGenes.txt" ) );
        addRowSet( new File( baseFileName + "LOOGenesInOrder.out.partialcon.ammon.txt.374.0.016424.topGenes.txt" ) );
        addRowSet( new File( baseFileName + "LOOGenesInOrder.space.ammon.txt.420.0.018435.topGenes.txt" ) );

    }

    public void addGeneSet( File file ) throws Exception {
        List<String> dataRows = Util.getRowNamesFromGeneNames( file.getAbsolutePath(), pair.getMatrixBDataRows() );
        files.put( file, dataRows );
    }

    public void addRowSet( File file ) throws Exception {
        List<String> rows = FileTools.getLines( file );
        files.put( file, rows );
    }

    public void run( int testIterations, boolean degree ) throws Exception {

        ParamKeeper results = new ParamKeeper();

        // AllenMatrixPair pair = AllenMatrixPairFactory.getSpaceAndExpressionPair( direction, useVirtual, logDistance
        // );

        for ( File fileObj : files.keySet() ) {

            List<String> dataRows = files.get( fileObj );

            try {
                boolean spearman = false;
                Map<String, String> result = pair.testDataRows( dataRows, testIterations, degree, spearman );
                spearman = true;
                result.putAll( pair.testDataRows( dataRows, testIterations, degree, spearman ) );
                result.put( "name", fileObj.getName() );
                if ( fileObj.exists() ) {
                    Set<String> originalGenes = new HashSet<String>( FileTools.getLines( fileObj ) );
                    result.put( "Unique lines in file", "" + originalGenes.size() );
                }
                result.put( "Unique genes used", "" + Util.getUniqueGenes( dataRows ).size() );
                results.addParamInstance( result );
            } catch ( Exception e ) {
                e.printStackTrace();
                log.info( "Error for:" + fileObj.toString() );
                // System.exit( 1 );
            }

        }

        results.writeExcel( outputDir + "results." + getFileBaseName() + ".degree." + degree + ".iterations."
                + testIterations + ".xls" );
    }

    public String getFileBaseName() {
        String matrixPairName = pair.getName();
        return matrixPairName + "." + direction.toString() + "." + pair.getMatrixB().rows() + "rows";
    }

    public void doGO() throws Exception {
        int minSize = 5;
        int maxSize = 300;

        boolean NCBI = true;

        for ( File file : files.keySet() ) {
            if ( file.getName().contains( "pattern" ) ) {
                log.info( file.toString() );
                List<String> allRows = pair.getMatrixB().getRowNames();
                RankedGeneListLoader rankedList = new RankedGeneListLoader( allRows, file.toString() );
                List<String> newLines = new LinkedList<String>( files.get( file ) );
                rankedList.addToFront( newLines );
                double threshold = newLines.size() / ( double ) allRows.size();
                String ermineJFile = "";
                if ( NCBI ) {
                    ermineJFile = rankedList.forErmineJNCBI();
                } else {
                    ermineJFile = rankedList.forErmineJ();
                }
                File ermineJFile2 = new File( ermineJFile );
                log.info( threshold );
                RankedGeneListLoader.runErmineJORA( threshold, ermineJFile2.toString(), NCBI, minSize, maxSize );
                ermineJFile2.renameTo( new File( ermineJFile2.toString() + "." + threshold + ".threshold" ) );
            }
        }

    }

    public static void main( String args[] ) throws Exception {
        // this code overlaps with EstrogenTests.java and inputfilegenefilter
        Direction direction = AnalyzeBAMSandAllenGenes.Direction.APPENDED;

        int testIterations = 100;
        boolean proximityControlled = false;
        boolean degree = false;
        boolean logDistance = true;
        boolean slow = false;
        RegressMatrix regressType = RegressMatrix.BOTH;
        boolean useVirtual = true;
        // CHECK remove non exp
        boolean removeNonExp = true;
        boolean coronal = true;

        testVolumecorrelations( testIterations, direction );
        System.exit( 1 );

        ConnectivityAndAllenExpressionMatrixPair pair;
        if ( proximityControlled ) {
            pair = ExpressionMatrixPairFactory.connectivityPartial( direction, slow, regressType, useVirtual,
                    removeNonExp, logDistance );
        } else {
            pair = ExpressionMatrixPairFactory.connectivityAndExpression( direction, useVirtual, removeNonExp );
        }
        if ( coronal ) pair.applyGeneFilter( new PlaneRemoveFilter( Plane.SAGITTAL ) );
        log.info( "Full correlation:" + pair.getCorrelation() );
        log.info( "Flattened correlation:" + pair.getFlattenedCorrelation( false ) );
        log.info( "Spearman Flattened correlation:" + pair.getFlattenedCorrelation( true ) );

        CellTypes celltypes = new CellTypes( pair, direction );
        // celltypes.addGeneSet( new File( "/grp/java/workspace/BAMSandAllen/data/G2Cdbgenelists/ColinsMousePSD.txt" )
        // );

        // celltypes.addMainLists();
        // celltypes.addOptimizedSets();
        celltypes.addFrontiersList();

        // List<String> rowsA = FileTools
        // .getLines( "/grp/java/workspace/BAMSandAllen/data/ABACellEnriched/patterns/top100-pattern-a.txt" );
        // List<String> rowsB = FileTools
        // .getLines( "/grp/java/workspace/BAMSandAllen/data/ABACellEnriched/patterns/top100-pattern-b.txt" );
        //
        // pair.testDataRows( rowsA, 0, true, true );
        // pair.testDataRows( rowsB, 0, true, true );

        // celltypes.createVennFile();

        // celltypes.doGO();
        celltypes.run( testIterations, degree );
        // testXcorrelations( testIterations, direction );

        pair.printDimensions();
        log.info( "Unique:" + Util.getUniqueGenes( pair.getMatrixBDataRows() ).size() );

    }
}
