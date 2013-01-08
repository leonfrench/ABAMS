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
package ubic.BAMSandAllen.MatrixPairs;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.lang.time.StopWatch;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import flanagan.analysis.Regression;

import ubic.BAMSandAllen.ABAMSDataMatrix;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes;
import ubic.BAMSandAllen.RankedGeneListLoader;
import ubic.BAMSandAllen.RegressionVector;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.ClassSelectors.BrainRegionClassSelector;
import ubic.BAMSandAllen.adjacency.CorrelationAdjacency;
import ubic.BAMSandAllen.adjacency.PartialCorrelationCompute;
import ubic.BAMSandAllen.geneFilters.EphrinGeneFilter;
import ubic.BAMSandAllen.geneFilters.EstrogenGeneFilter;
import ubic.BAMSandAllen.geneFilters.GeneFilter;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.math.CorrelationStats;

public class ConnectivityAndAllenPartialExpressionMatrixPair extends ConnectivityAndAllenExpressionMatrixPair {
    // enum for which matrix is regressed on, the traditional partial mantel test uses both
    public enum RegressMatrix {
        CONNECTIVITY, EXPRESSION, BOTH
    };

    private static Log log = LogFactory.getLog( ConnectivityAndAllenPartialExpressionMatrixPair.class.getName() );
    // this is the explainatory matrix
    ABAMSDataMatrix regressMatrix;
    // this is the targetMatrix (either matrixA or matrixB)
    RegressMatrix regressType;

    public ConnectivityAndAllenPartialExpressionMatrixPair( BrainRegionClassSelector selector, boolean doLog,
            boolean onePlusLog, boolean square, double zeroReplacement, String matrixName, Direction direction,
            ABAMSDataMatrix regressMatrix ) throws Exception {
        // assume regression is explaining matrixB(Expression), not connectivity
        this( selector, doLog, onePlusLog, square, zeroReplacement, matrixName, direction, regressMatrix,
                RegressMatrix.EXPRESSION );
    }

    public ConnectivityAndAllenPartialExpressionMatrixPair( BrainRegionClassSelector selector, boolean doLog,
            boolean onePlusLog, boolean square, double zeroReplacement, String matrixName, Direction direction,
            ABAMSDataMatrix regressMatrix, RegressMatrix regressType ) throws Exception {
        super( selector, doLog, onePlusLog, square, zeroReplacement, matrixName, direction );
        this.regressMatrix = regressMatrix;
        this.regressType = regressType;
    }

    public ConnectivityAndAllenPartialExpressionMatrixPair( String connectivityMatrixFile, boolean doLog,
            boolean onePlusLog, double zeroReplacement, String matrixName, ABAMSDataMatrix regressMatrix,
            RegressMatrix regressType ) throws Exception {
        super( connectivityMatrixFile, doLog, onePlusLog, zeroReplacement, matrixName );
        this.regressMatrix = regressMatrix;
        this.regressType = regressType;
    }

    /**
     * This returns which matrices will be used in the partial regression
     * 
     * @return
     */
    public Set<ABAMSDataMatrix> getTargetMatrix() {
        Set<ABAMSDataMatrix> result = new HashSet<ABAMSDataMatrix>();
        if ( regressType == RegressMatrix.CONNECTIVITY || regressType == RegressMatrix.BOTH ) result.add( matrixA );
        if ( regressType == RegressMatrix.EXPRESSION || regressType == RegressMatrix.BOTH ) result.add( matrixB );
        return result;
    }

    public RegressMatrix getRegressType() {
        return regressType;
    }

    // below should be refactored
    public void setRegressionComputeMatrixA( boolean rCompute ) {
        ( ( PartialCorrelationCompute ) matrixA.getAdjacencyCompute() ).setComputeRegression( rCompute );
    }

    public void setRegressionComputeMatrixB( boolean rCompute ) {
        ( ( PartialCorrelationCompute ) matrixB.getAdjacencyCompute() ).setComputeRegression( rCompute );
    }

    public boolean getRegressionComputeMatrixB()   {
        return ( ( PartialCorrelationCompute ) matrixB.getAdjacencyCompute() ).getComputeRegression();
    }

    public void setTrianglesMatrixA( RegressionVector triangles ) {
        ( ( PartialCorrelationCompute ) matrixA.getAdjacencyCompute() ).setTriangles( triangles );
    }

    public void setTrianglesMatrixB( RegressionVector triangles ) {
        ( ( PartialCorrelationCompute ) matrixB.getAdjacencyCompute() ).setTriangles( triangles );
    }

    public RegressionVector getTrianglesMatrixA() {
        return ( ( PartialCorrelationCompute ) matrixA.getAdjacencyCompute() ).getTriangles();
    }

    public RegressionVector getTrianglesMatrixB() {
        return ( ( PartialCorrelationCompute ) matrixB.getAdjacencyCompute() ).getTriangles();
    }

    public double run() throws Exception {
        clean();
        sameSpace();

        // assume regressmatrix is in ABAMS region space
        setupAdjacencyComputes();
        log.info( "After same space and renaming" );

        printDimensions();
        double correlation = getCorrelation();

        log.info( "Correlation:" + correlation );
        return correlation;

    }

    private void setupAdjacencyComputes() throws Exception {
        regressMatrix = regressMatrix.retainColumns( getAllenDataColNames() );
        for ( ABAMSDataMatrix targetMatrix : getTargetMatrix() ) {
            PartialCorrelationCompute adjacencyCompute = new PartialCorrelationCompute( targetMatrix, regressMatrix,
                    new CorrelationAdjacency( targetMatrix ) );
            targetMatrix.setAdjacencyCompute( adjacencyCompute );
            targetMatrix.setName( targetMatrix.getName() + ".Residuals(" + regressMatrix.getName() + ")" );

            // for display and R purposes
            DoubleMatrix<String, String> explainAdjacenty = adjacencyCompute.getExplainMatrix().getAdjacency();
            Util.writeRTable(
                    baseName + "." + adjacencyCompute.getExplainMatrix().getName() + ".Explain.Adjacency.txt",
                    explainAdjacenty );

        }
    }

    // public void removeMatrixBDataRows( Collection<String> rows ) throws Exception {
    // super.removeMatrixBDataRows( rows );
    // setupAdjacencyComputes();
    // }

    public void printDimensions() {
        super.printDimensions();
        log.info( matrixB.getAdjacency().rows() + " x " + matrixB.getAdjacency().columns() );
    }

    public double getFlattenedCorrelation( boolean spearman ) throws Exception {
        // are they in the right order?
        boolean removeNan = true;
        DoubleMatrix<String, String> matrixASums = Util.columnSums( matrixA, removeNan );
        DoubleMatrix<String, String> matrixBSums = Util.columnSums( matrixB, removeNan );
        if ( !matrixASums.getColNames().equals( matrixB.getColNames() ) )
            throw new RuntimeException( "Error column names do not match" );

        if ( !regressMatrix.getColNames().equals( matrixASums.getColNames() ) )
            throw new RuntimeException( "Error x matrix names don't match" );

        double[] xVec = regressMatrix.getRowByName( "x" );
        double[] matrixAVec = matrixASums.getRow( 0 );
        double[] matrixBVec = matrixBSums.getRow( 0 );

        int length = xVec.length;
        log.info( "X lenght:" + length );
        double matrixAVecControlled[] = matrixAVec;
        double matrixBVecControlled[] = matrixBVec;

        if ( regressType.equals( RegressMatrix.BOTH ) || regressType.equals( RegressMatrix.CONNECTIVITY ) ) {
            Regression regression = new Regression( xVec, matrixAVec );
            regression.linear();
            regression.toString();

            matrixAVecControlled = regression.getResiduals();
        }

        if ( regressType.equals( RegressMatrix.BOTH ) || regressType.equals( RegressMatrix.EXPRESSION ) ) {
            Regression regression = new Regression( xVec, matrixBVec );
            regression.linear();
            regression.toString();

            matrixBVecControlled = regression.getResiduals();
        }

        int rows = 5;

        // // write out for testing
        // DenseDoubleMatrix<String, String> dataPoints = new DenseDoubleMatrix<String, String>( rows, length );
        // dataPoints.setColumnNames( matrixA.getColNames() );
        // addToMatrix( "x", xVec, dataPoints );
        // addToMatrix( "degree", matrixAVec, dataPoints );
        // addToMatrix( "exp", matrixBVec, dataPoints );
        // addToMatrix( "degreeCon", matrixAVecControlled, dataPoints );
        // addToMatrix( "expCon", matrixBVecControlled, dataPoints );
        //
        // Util.writeRTable( "/grp/java/workspace/BAMSandAllen/data/ABACellEnriched/patterns/flattened.ap.."
        // + matrixB.rows() + "exp.txt", dataPoints );

        if ( !spearman ) {
            return CorrelationStats.correl( matrixAVecControlled, matrixBVecControlled );
        } else {
            return Util.spearmanCorrel( matrixAVecControlled, matrixBVecControlled );
        }
    }

    // ugly, should be worked into regressionvector as superclass or something (regression vector is for adjacencies)
    public void addToMatrix( String rowName, double[] values, DenseDoubleMatrix<String, String> matrix ) {
        matrix.addRowName( rowName );
        int r = matrix.getRowIndexByName( rowName );
        for ( int i = 0; i < values.length; i++ ) {
            matrix.set( r, i, values[i] );
        }
    }

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        // AllenAtlasAnnotationLoader spaceLoader = new AllenAtlasAnnotationLoader();
        // DoubleMatrix<String, String> spaceMatrix = spaceLoader.getCenterMatrix();
        // ABAMSDataMatrix spaceABAMS = new ABAMSDataMatrix( spaceMatrix, "Space", new EuclidAdjacency() );
        Direction direction = AnalyzeBAMSandAllenGenes.Direction.INCOMING;
        boolean makeVirtualRegions = true;

        boolean squareMatrix = false;
        boolean logDistance = true;
        boolean slow = false;
        RegressMatrix regressType = RegressMatrix.CONNECTIVITY;
        boolean useVirtual = true;
        boolean removeNonExp = true;

        // ConnectivityAndAllenPartialExpressionMatrixPair pair = ExpressionMatrixPairFactory.connectivityPartial(
        // direction, slow, regressType, useVirtual, removeNonExp, logDistance );
        MatrixPair pair = AllenMatrixPairFactory.getSpaceAndExpressionPair( direction, removeNonExp, logDistance );
        log.info( pair.getCorrelation() );

        pair.writeRMatrices();

        // pair.test( 1000 );

        // System.exit( 1 );

        String filename = "/grp/java/workspace/BAMSandAllen/data/rankedGenes/near final ammon/";
        // filename += "LOOGenesInOrder.in.partialcon.ammon.txt.452.0.01985.topGenes.txt";
        // filename += "LOOGenesInOrder.out.partialcon.ammon.txt.374.0.016424.topGenes.txt";
        filename += "LOOGenesInOrder.space.ammon.txt.420.0.018435.topGenes.txt";
        RankedGeneListLoader aLook = new RankedGeneListLoader( filename );
        pair.setMatrixBDataRows( aLook.getLines() );
        log.info( pair.getCorrelation() );
        pair.test( 1000 );
        System.exit( 1 );

        String division = "Hindbrain";
        // divisions.add( "Interbrain" );
        // divisions.add( "Midbrain" );
        // divisions.add( "Cerebrum" );
        // divisions.add( "Cerebellum" );

        // for ( String rowName : pair.getMatrixBDataRows() ) {
        // log.info( rowName + ":" );
        // pair.setToSingleGene( rowName );
        // pair.getCorrelation();
        // double pval = pair.test( 1000, false );
        // if ( pval < 0.05 ) {
        // log.info( "LOW:" + rowName );
        // }
        // }

        System.exit( 1 );

        // pair.applyGeneFilter( new EphrinGeneFilter() );
        // pair.applyGeneFilter( new EstrogenGeneFilter() );
        // log.info( pair.getCorrelation() );
        // log.info( pair.test( 1000, false ) );

        // log.info( pair.getMatrixBDataRows().size() );

        // log.info( pair.getCorrelation() );
        // pair.test( 1000 );

        // log.info( pair.testOldHypoth() );

        // this pair is used to convert the space pair region names to BAMS names
        ConnectivityAndAllenSpacePair spacePair = new ConnectivityAndAllenSpacePair( new BrainRegionClassSelector(),
                squareMatrix, null, direction, logDistance );
        if ( makeVirtualRegions ) spacePair.makeVirtualRegions();
        spacePair.run();

        ConnectivityAndAllenPartialExpressionMatrixPair forR = new ConnectivityAndAllenPartialExpressionMatrixPair(
                new BrainRegionClassSelector(), true, false, squareMatrix, Double.NaN, "NewEnergies", direction,
                spacePair.getMatrixB(), regressType );

        // forR = new ConnectivityAndAllenExpressionMatrixPair( new BrainRegionClassSelector(), true, false,
        // squareMatrix,
        // Double.NaN, "NewEnergies", AnalyzeBAMSandAllenGenes.Direction.OUTGOING );

        // forR.writeRMatrices();
        if ( makeVirtualRegions ) forR.makeVirtualRegions();
        log.info( forR.run() );
        log.info( forR.getCorrelation() );

        // forR.runAllenStyle();
        // forR.writeRMatrices();
        System.exit( 1 );

        // forR.runAllenStyle();

        String rowName = forR.getMatrixBDataRows().get( 1 );

        StopWatch timer = new StopWatch();
        timer.start();
        log.info( forR.correWithoutMatrixBDataRow( rowName ) );
        log.info( forR.correWithoutMatrixBDataRow( rowName ) );
        log.info( forR.correWithoutMatrixBDataRow( rowName ) );
        log.info( forR.correWithoutMatrixBDataRow( rowName ) );
        log.info( forR.correWithoutMatrixBDataRow( rowName ) );
        log.info( forR.correWithoutMatrixBDataRow( rowName ) );
        log.info( forR.correWithoutMatrixBDataRow( rowName ) );
        log.info( forR.correWithoutMatrixBDataRow( rowName ) );
        log.info( timer.getTime() );

        forR.setRegressionComputeMatrixB( false );
        timer = new StopWatch();
        timer.start();
        log.info( forR.correWithoutMatrixBDataRow( rowName ) );
        log.info( forR.correWithoutMatrixBDataRow( rowName ) );
        log.info( forR.correWithoutMatrixBDataRow( rowName ) );
        log.info( forR.correWithoutMatrixBDataRow( rowName ) );
        log.info( forR.correWithoutMatrixBDataRow( rowName ) );
        log.info( forR.correWithoutMatrixBDataRow( rowName ) );
        log.info( forR.correWithoutMatrixBDataRow( rowName ) );
        log.info( forR.correWithoutMatrixBDataRow( rowName ) );
        log.info( "After shutting down regression compute:" + timer.getTime() );

    }
}
