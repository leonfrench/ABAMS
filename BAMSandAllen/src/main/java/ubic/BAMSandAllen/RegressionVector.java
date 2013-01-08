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
package ubic.BAMSandAllen;

import java.io.Serializable;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.AllenDataLoaders.AllenAtlasAnnotationLoader;
import ubic.BAMSandAllen.AllenDataLoaders.AllenCatalogMatrices2;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.ClassSelectors.BrainRegionClassSelector;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenDimensionsPair;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenExpressionMatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenNomenclaturePair;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenPartialExpressionMatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenSpacePair;
import ubic.BAMSandAllen.MatrixPairs.ExpressionMatrixPairFactory;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenPartialExpressionMatrixPair.RegressMatrix;
import ubic.BAMSandAllen.adjacency.CorrelationAdjacency;
import ubic.BAMSandAllen.adjacency.EuclidAdjacency;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.util.FileTools;

public class RegressionVector implements Serializable {
    private static Log log = LogFactory.getLog( RegressionVector.class.getName() );

    int spot;
    int cols;
    List<String> colNames;
    String filename;
    DoubleMatrix<String, String[]> dataPoints;

    public RegressionVector( int rows, int cols ) {
        this( rows, cols, true );
    }

    public RegressionVector( int rows, int cols, boolean adjacency ) {
        spot = 0;
        this.cols = cols;
        filename = "";
        int length;
        if ( adjacency ) {
            length = cols * ( cols - 1 ) / 2;
        } else {
            length = cols;
        }
        dataPoints = new DenseDoubleMatrix<String, String[]>( rows, length );
    }

    public DoubleMatrix<String, String> getRowAsMatrix( String rowName ) {
        // reverse of the add method
        DoubleMatrix<String, String> rowMatrix = new DenseDoubleMatrix<String, String>( cols, cols );
        rowMatrix.setColumnNames( colNames );
        rowMatrix.setRowNames( colNames );
        for ( String[] colPair : dataPoints.getColNames() ) {
            String a = colPair[0];
            String b = colPair[1];
            double value = dataPoints.getByKeys( rowName, colPair );
            rowMatrix.setByKeys( a, b, value );
            rowMatrix.setByKeys( b, a, value );
        }
        return rowMatrix;
    }

    /*
     * dangerous, assumes the double[] is in the right order! used for regression
     */
    public void add( String rowName, double[] entries ) {
        if ( dataPoints.getRowNames().contains( rowName ) ) {
            put( rowName, entries );
        } else {
            dataPoints.addRowName( rowName );
            put( rowName, entries );
            spot++;
        }
    }

    protected void put( String rowName, double[] entries ) {
        if ( entries.length != dataPoints.columns() ) throw new RuntimeException( "Not array of right length" );
        int rowIndex = dataPoints.getRowIndexByName( rowName );
        for ( int i = 0; i < entries.length; i++ ) {
            dataPoints.set( rowIndex, i, entries[i] );
        }
    }

    protected void put( String rowName, DoubleMatrix<String, String> adjacencyMatrix ) {
        if ( adjacencyMatrix.columns() != adjacencyMatrix.rows() ) throw new RuntimeException( "Not a square matrix" );
        if ( adjacencyMatrix.columns() != cols ) throw new RuntimeException( "Not matrix of right dimensions" );

        if ( spot == 0 ) {
            setupColNames( adjacencyMatrix.getColNames() );
        }

        for ( String[] namePair : dataPoints.getColNames() ) {
            double value = adjacencyMatrix.getByKeys( namePair[0], namePair[1] );
            dataPoints.setByKeys( rowName, namePair, value );
        }
    }

    public void add( String rowName, DoubleMatrix<String, String> adjacencyMatrix ) {
        if ( dataPoints.getRowNames().contains( rowName ) ) {
            put( rowName, adjacencyMatrix );
        } else {
            dataPoints.addRowName( rowName );
            put( rowName, adjacencyMatrix );
            filename += rowName + ".";
            spot++;
        }
    }

    public double[] getRow( String name ) {
        return dataPoints.getRowByName( name );
    }

    public void writeRTable() throws Exception {
        // write it out for R
        List<String> newColNames = new LinkedList<String>();
        for ( String[] colName : dataPoints.getColNames() ) {
            newColNames.add( colName[0] + ".AND." + colName[1] );
        }
        DoubleMatrix<String, String> writeMatrix = new DenseDoubleMatrix<String, String>( dataPoints.asArray() );
        writeMatrix.setRowNames( dataPoints.getRowNames() );
        writeMatrix.setColumnNames( newColNames );

        String filename = SetupParameters.getDataFolder();
        for ( String rowName : dataPoints.getRowNames() ) {
            if ( filename.length() > 150 ) break;
            filename += rowName + ".";
        }
        filename += "forR.txt";

        Util.writeRTable( filename, writeMatrix );
        log.info( "Wrote to: " + filename );

    }

    public void setupColNames( List<String> names ) {
        this.colNames = names;
        int size = names.size();
        for ( int i = 0; i < size; i++ ) {
            for ( int j = i + 1; j < size; j++ ) {
                String colName = names.get( i );
                String rowName = names.get( j );
                String[] x = { colName, rowName };
                dataPoints.addColumnName( x );
            }
        }
    }

    public static void main( String args[] ) throws Exception {
        forPatrick();
    }

    /**
     * @param args
     */
    public static void fromPaper() throws Exception {
        boolean squareMatrix = false;
        boolean virtualRegions = true;
        boolean logDistance = true;
        Direction direction = AnalyzeBAMSandAllenGenes.Direction.OUTGOING; // was anydirection
        Set<String> colNames = new HashSet<String>( ExpressionMatrixPairFactory.getUsedCols( squareMatrix, direction ) );

        // get expression
        // get connectivity
        boolean useVirtual = true;
        boolean removeNonExp = true;
        ConnectivityAndAllenExpressionMatrixPair expressionPair = ExpressionMatrixPairFactory
                .connectivityAndExpression( direction, useVirtual, removeNonExp );

        // get space in BAMS space
        ConnectivityAndAllenSpacePair spacePair = new ConnectivityAndAllenSpacePair( new BrainRegionClassSelector(),
                squareMatrix, colNames, direction, logDistance );
        // !! using virtual regions!
        spacePair.makeVirtualRegions();
        spacePair.run();

        // get nomenclature in BAMS space
        ConnectivityAndAllenNomenclaturePair nomenPair = new ConnectivityAndAllenNomenclaturePair(
                new BrainRegionClassSelector(), squareMatrix, false, colNames, direction );
        if ( virtualRegions ) nomenPair.makeVirtualRegions();
        nomenPair.run();

        log.info( "_____________________________" );

        RegressionVector tryOut = new RegressionVector( 9, colNames.size() );
        tryOut.add( "Expression", expressionPair.getAdjacencyB() );
        tryOut.add( "Distance", spacePair.getAdjacencyB() );
        tryOut.add( "Connectivity", spacePair.getAdjacencyA() );
        tryOut.add( "Nomenclature", nomenPair.getAdjacencyB() );
        boolean add = true;

        expressionPair.setToSingleGene( "Pcp2[77413702]", add );
        tryOut.add( "Pcp2", expressionPair.getAdjacencyB() );

        expressionPair.setToSingleGene( "Pgrmc1[797]", add );
        tryOut.add( "Pgrmc1", expressionPair.getAdjacencyB() );

        expressionPair.setToSingleGene( "Rgs3[1822]", add );
        tryOut.add( "Rgs3", expressionPair.getAdjacencyB() );
        tryOut.add( "Connections", expressionPair.getConnectionMatrix( true ) );

        boolean slow = false;
        RegressMatrix regressType = RegressMatrix.BOTH;
        expressionPair = ExpressionMatrixPairFactory.connectivityPartial( direction, slow, regressType, useVirtual,
                removeNonExp, logDistance );
        String filename = "/grp/java/workspace/BAMSandAllen/data/rankedGenes/near final ammon/";
        filename += "LOOGenesInOrder.out.partialcon.ammon.txt.374.0.016424.topGenes.txt";
        RankedGeneListLoader aLook = new RankedGeneListLoader( filename );
        expressionPair.setMatrixBDataRows( aLook.getLines() );
        log.info( expressionPair.getCorrelation() );
        tryOut.add( "TopExpression", expressionPair.getAdjacencyB() );

        tryOut.writeRTable();
    }

    public static void forPatrick() throws Exception {
        boolean squareMatrix = false;
        boolean virtualRegions = true;
        boolean logDistance = true;
        Direction direction = AnalyzeBAMSandAllenGenes.Direction.ANYDIRECTION; // was anydirection

        Set<String> colNames = new HashSet<String>( ExpressionMatrixPairFactory.getUsedCols( squareMatrix, direction ) );

        // get expression
        // get connectivity
        boolean useVirtual = true;
        boolean removeNonExp = true;
        ConnectivityAndAllenExpressionMatrixPair expressionPair = ExpressionMatrixPairFactory
                .connectivityAndExpression( direction, useVirtual, removeNonExp );

        // get space in BAMS space
        ConnectivityAndAllenSpacePair spacePair = new ConnectivityAndAllenSpacePair( new BrainRegionClassSelector(),
                squareMatrix, colNames, direction, logDistance );
        // !! using virtual regions!
        spacePair.makeVirtualRegions();
        spacePair.run();

        // get nomenclature in BAMS space
        ConnectivityAndAllenNomenclaturePair nomenPair = new ConnectivityAndAllenNomenclaturePair(
                new BrainRegionClassSelector(), squareMatrix, false, colNames, direction );
        if ( virtualRegions ) nomenPair.makeVirtualRegions();
        nomenPair.run();

        log.info( "_____________________________" );

        String patternAStr = "/grp/java/workspace/BAMSandAllen/data/ABACellEnriched/patterns/top100-pattern-a.txt";
        String patternBStr = "/grp/java/workspace/BAMSandAllen/data/ABACellEnriched/patterns/top100-pattern-b.txt";

        List<String> patternArows = FileTools.getLines( patternAStr );
        List<String> patternBrows = FileTools.getLines( patternBStr );

        RegressionVector tryOut = new RegressionVector( 11 + 2*(patternArows.size() + patternBrows.size()), colNames.size() );
        tryOut.add( "Expression", expressionPair.getAdjacencyB() );
        tryOut.add( "Distance", spacePair.getAdjacencyB() );
        tryOut.add( "Connectivity", spacePair.getAdjacencyA() );
        tryOut.add( "Nomenclature", nomenPair.getAdjacencyB() );
        boolean add = true;
        boolean subtract = false;

        // Mef2c[79567505]
        expressionPair.setToSingleGene( "Mef2c[79567505]", add );
        tryOut.add( "Mef2cSum", expressionPair.getAdjacencyB() );

        expressionPair.setToSingleGene( "Mef2c[79567505]", subtract );
        tryOut.add( "Mef2cDiff", expressionPair.getAdjacencyB() );

        // Nefh[74512048].subtract
        expressionPair.setToSingleGene( "Nefh[74512048]", add );
        tryOut.add( "NefhSum", expressionPair.getAdjacencyB() );

        expressionPair.setToSingleGene( "Nefh[74512048]", subtract );
        tryOut.add( "NefhDiff", expressionPair.getAdjacencyB() );

        expressionPair.setToSingleGene( "Pgrmc1[797]", add );
        tryOut.add( "Pgrmc1Sum", expressionPair.getAdjacencyB() );

        expressionPair.setToSingleGene( "Pgrmc1[797]", subtract );
        tryOut.add( "Pgrmc1Diff", expressionPair.getAdjacencyB() );

        patternArows.addAll( patternBrows );

        for ( String row : patternArows ) {
            expressionPair.setToSingleGene( row, add );
            tryOut.add( row + ".Sum", expressionPair.getAdjacencyB() );

            expressionPair.setToSingleGene( row, subtract );
            tryOut.add( row + ".Diff", expressionPair.getAdjacencyB() );
        }

        tryOut.add( "Connections", expressionPair.getConnectionMatrix( true ) );

        tryOut.writeRTable();
    }

    // public static void mainOld( String[] args ) throws Exception {
    //
    // AllenCatalogMatrices2 allenMatrices = new AllenCatalogMatrices2();
    // DoubleMatrix<String, String> energyMatrix = allenMatrices.getFromDisk( "NewEnergies" );
    // energyMatrix = Util.logMatrix( energyMatrix, Double.NaN );
    //
    // ABAMSDataMatrix energyABAMS = new ABAMSDataMatrix( energyMatrix, "NewEnergies", new CorrelationAdjacency(
    // energyMatrix ) );
    //
    // AllenAtlasAnnotationLoader spaceLoader = new AllenAtlasAnnotationLoader();
    // DoubleMatrix<String, String> spaceMatrix = spaceLoader.getCenterMatrix();
    // ABAMSDataMatrix spaceABAMS = new ABAMSDataMatrix( spaceMatrix, "Space", new EuclidAdjacency() );
    //
    // // DoubleMatrix<String, String> dimsMatrix = spaceLoader.getDimensionsMatrix();
    // // // eculid makes it significant, use volume then!
    // // ABAMSDataMatrix dimsABAMS = new ABAMSDataMatrix( dimsMatrix, "Dimensions", new VolumeDiffAdjacency() );
    // //
    // // StructureCatalogLoader loader = new StructureCatalogLoader();
    // // DoubleMatrix<String, String> nomenclatureMatrix = loader.getNomenclatureMatrix();
    // // ABAMSDataMatrix nomenclatureABAMS = new ABAMSDataMatrix( nomenclatureMatrix, "Nomenclature",
    // // new DotProductAdjacency() );
    //
    // Direction direction = AnalyzeBAMSandAllenGenes.Direction.OUTGOING;
    //
    // Set<String> colNames = new HashSet<String>( ExpressionMatrixPairFactory.getUsedCols( squareMatrix, direction ) );
    //
    // // ConnectivityAndAllenExpressionMatrixPair expressionPair = new ConnectivityAndAllenExpressionMatrixPair(
    // // new BrainRegionClassSelector(), true, false, squareMatrix, Double.NaN, "NewEnergies", direction );
    // // System.exit( 1 );
    // // forR.writeRMatrices();
    // // expressionPair.run();
    //
    // ConnectivityAndAllenSpacePair spacePair = new ConnectivityAndAllenSpacePair( new BrainRegionClassSelector(),
    // squareMatrix, colNames, direction );
    // spacePair.run();
    //
    // ConnectivityAndAllenNomenclaturePair x = new ConnectivityAndAllenNomenclaturePair(
    // new BrainRegionClassSelector(), squareMatrix, false, colNames, direction );
    // x.run();
    //
    // ConnectivityAndAllenDimensionsPair dims = new ConnectivityAndAllenDimensionsPair(
    // new BrainRegionClassSelector(), squareMatrix, colNames, direction );
    // dims.run();
    // dims.writeImages();
    //
    // log.info( "_____________________________" );
    //
    // RegressionVector tryOut = new RegressionVector( 7, colNames.size() );
    // tryOut.add( "Expression", expressionPair.getAdjacencyB() );
    // tryOut.add( "Distance", spacePair.getAdjacencyB() );
    // tryOut.add( "Connectivity", spacePair.getAdjacencyA() );
    // // spacePair.
    // tryOut.add( "Nomenclature", x.getAdjacencyB() );
    // tryOut.add( "Dimensions", dims.getAdjacencyB() );
    // tryOut.add( "Connections", expressionPair.getConnectionMatrix( true ) );
    //
    // String filename = "/grp/java/workspace/BAMSandAllen/data/rankedGenes/near final ammon/";
    // filename += "LOOGenesInOrder.out.partialcon.ammon.txt.374.0.016424.topGenes.txt";
    // RankedGeneListLoader aLook = new RankedGeneListLoader( filename );
    // expressionPair.setMatrixBDataRows( aLook.getLines() );
    // log.info( expressionPair.getCorrelation() );
    // tryOut.add( "TopExpression", expressionPair.getAdjacencyB() );
    //
    // // expressionPair.setDataRows( rows
    //
    // tryOut.writeRTable();
    //
    // }
}
