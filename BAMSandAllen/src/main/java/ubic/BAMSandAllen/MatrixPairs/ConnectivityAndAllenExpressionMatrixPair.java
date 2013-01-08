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

import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.StringTokenizer;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import com.hp.hpl.jena.util.FileUtils;

import ubic.BAMSandAllen.ABAMSDataMatrix;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes;
import ubic.BAMSandAllen.RankedGeneListLoader;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.AllenDataLoaders.AllenCatalogMatrices2;
import ubic.BAMSandAllen.AllenDataLoaders.ImageSeriesInfoLoader;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.ClassSelectors.BrainRegionClassSelector;
import ubic.BAMSandAllen.FocusedAnalysis.ExploreRegionNames;
import ubic.BAMSandAllen.adjacency.CorrelationAdjacency;
import ubic.BAMSandAllen.adjacency.SingleRowAdjacency;
import ubic.BAMSandAllen.geneFilters.GeneFilter;
import ubic.basecode.dataStructure.StringToStringSetMap;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.math.CorrelationStats;
import ubic.basecode.util.FileTools;

public class ConnectivityAndAllenExpressionMatrixPair extends ConnectivityAndAllenDataPair {
    private static Log log = LogFactory.getLog( ConnectivityAndAllenExpressionMatrixPair.class.getName() );
    String matrixName;
    boolean doLog;
    double zeroReplacement;
    boolean onePlusLog;

    public ConnectivityAndAllenExpressionMatrixPair( BrainRegionClassSelector selector, boolean doLog,
            boolean onePlusLog, boolean square, double zeroReplacement, String matrixName, Direction direction )
            throws Exception {
        super( selector, square, direction );
        loadExpression( doLog, onePlusLog, zeroReplacement, matrixName );
    }

    /**
     * For loading literature matrices.
     * 
     * @throws Exception
     */
    public ConnectivityAndAllenExpressionMatrixPair( String filename, boolean doLog, boolean onePlusLog,
            double zeroReplacement, String matrixName ) throws Exception {
        super( filename );
        loadExpression( doLog, onePlusLog, zeroReplacement, matrixName );
    }

    private void loadExpression( boolean doLog, boolean onePlusLog, double zeroReplacement, String matrixName ) {
        this.matrixName = matrixName;
        this.zeroReplacement = zeroReplacement;
        this.doLog = doLog;
        this.onePlusLog = onePlusLog;

        // load the data here
        try {
            AllenCatalogMatrices2 allenMatrices = new AllenCatalogMatrices2();
            DoubleMatrix<String, String> dataMatrix = allenMatrices.getFromDisk( matrixName );
            matrixB = new ABAMSDataMatrix( dataMatrix, matrixName, new CorrelationAdjacency( dataMatrix ) );
        } catch ( Exception e ) {
            e.printStackTrace();
            System.exit( 1 );
        }

        matrixB = matrixB.removeZeroColumns();

        if ( doLog ) {
            if ( onePlusLog ) {
                matrixB = matrixB.replaceMatrix( Util.log1pMatrix( matrixB ) );
            } else {
                matrixB = matrixB.replaceMatrix( Util.logMatrix( matrixB, zeroReplacement ) );
            }
        }

        log.info( "got Matrix B - Expression level" );
    }

    public void setToSingleGene( String rowName ) {
        boolean add = false;
        setToSingleGene( rowName, add );
    }

    public void setToSingleGene( String rowName, boolean add ) {
        matrixB.setAdjacencyCompute( new SingleRowAdjacency( matrixB, rowName, add ) );

    }

    public void applyGeneFilter( GeneFilter filter ) throws Exception {
        applyMatrixBRowFilter( filter );
    }

    public void mergeExpressionTest() {
        List<String> rows = getMatrixBDataRows();
        StringToStringSetMap geneToRows = new StringToStringSetMap();

        for ( String rowName : rows ) {
            String gene = rowName.substring( 0, rowName.indexOf( "[" ) );
            geneToRows.put( gene, rowName );
        }
        // remove the null gene from ABA data, I suspect they are QC fails
        log.info( "Null genes:" + geneToRows.getSize( "null" ) );
        geneToRows.remove( "null" );
        log.info( "Genes:" + geneToRows.size() + " rows:" + rows.size() );
        int totalcount = 0;
        int morethan1 = 0;
        double total = 0;
        for ( String gene : geneToRows.keySet() ) {
            // correlate?
            Set<String> dataRows = geneToRows.get( gene );
            if ( dataRows.size() > 1 ) morethan1++;
            for ( String dataRow : dataRows ) {
                double[] row1 = matrixB.getRowByName( dataRow );
                for ( String dataRow2 : dataRows ) {
                    if ( !dataRow.equals( dataRow2 ) ) {
                        double[] row2 = matrixB.getRowByName( dataRow2 );
                        log.info( gene + " " + dataRow2 + "-" + dataRow + " size:" + dataRows.size() + " cor:"
                                + CorrelationStats.correl( row1, row2 ) );
                        if ( !Double.isNaN( CorrelationStats.correl( row1, row2 ) ) ) {
                            total += CorrelationStats.correl( row1, row2 );
                            totalcount++;
                        }
                    }
                }

            }
        }
        log.info( "Average correlation:" + ( total / ( double ) totalcount ) );
        log.info( "1+:" + morethan1 );
        // go through all rows, find gene name matches
        // merge
        // check correlation

    }

    // should be in a driver class
    public static void usingOutgoinRegions() throws Exception {
        boolean squareMatrix = false;
        Direction direction = AnalyzeBAMSandAllenGenes.Direction.INCOMING;
        ConnectivityAndAllenExpressionMatrixPair forR = new ConnectivityAndAllenExpressionMatrixPair(
                new BrainRegionClassSelector(), true, false, squareMatrix, Double.NaN, "NewEnergies", direction );

        forR.removeAllenCols( ExpressionMatrixPairFactory.getUsedCols( squareMatrix, direction.OUTGOING ) );
        forR.run();

        log.info( "using outgoing regions:" + forR.getCorrelation() );
        forR.test( 1000 );

    }

    public double getDegreeCorrelation( String rowName ) {
        double[] exp = matrixB.getRowByName( rowName );
        DoubleMatrix<String, String> connectionDegrees = Util.columnSums( matrixA );
        double[] degrees = connectionDegrees.getRow( 0 );
        return CorrelationStats.correl( degrees, exp );
    }

    public double getRankDegreeCorrelation( String rowName ) {
        double[] exp = matrixB.getRowByName( rowName );
        DoubleMatrix<String, String> connectionDegrees = Util.columnSums( matrixA );
        double[] degrees = connectionDegrees.getRow( 0 );
        return Util.spearmanCorrel( degrees, exp );
    }

    public void printMappingRelations() {
        super.printMappingRelations();
        log.info( "Virtual regions:" + virtualRegions.size() );
    }

    public boolean isVirtualRegion( String colName ) {
        return virtualRegions.contains( colName );
    }

    public Set<String> getUniqueGenes() {
        Set<String> result = new HashSet<String>();
        for ( String rowName : matrixB.getRowNames() ) {
            result.add( ImageSeriesInfoLoader.getGeneNameFromRowName( rowName ) );
        }
        return result;
    }


    public static void main( String[] args ) throws Exception {

        Direction direction;
        direction = AnalyzeBAMSandAllenGenes.Direction.INCOMING;

        // usingOutgoinRegions();
        boolean useVirtual = true;
        boolean removeNonExp = true;
        boolean square = false;
        ConnectivityAndAllenExpressionMatrixPair forR;
        forR = ExpressionMatrixPairFactory.connectivityAndExpression( direction, useVirtual, removeNonExp );
        // forR.testZeroes();
        // forR.writeImages();

        System.exit( 1 );

        forR.removeZeroConnectionRows();
        log.info( forR.getCorrelation() );
        forR.printDimensions();

        String filename = "/grp/java/workspace/BAMSandAllen/data/rankedGenes/near final ammon/";
        filename += "LOOGenesInOrder.out.partialcon.ammon.txt.374.0.016424.topGenes.txt";

        RankedGeneListLoader aLook = new RankedGeneListLoader( filename );
        forR.setMatrixBDataRows( aLook.getLines() );
        log.info( forR.getCorrelation() );
        forR.test( 1000 );

        forR.printDimensions();
        forR.writeRMatrices();
        System.exit( 1 );

        direction = AnalyzeBAMSandAllenGenes.Direction.OUTGOING;
        ConnectivityAndAllenExpressionMatrixPair forR2;
        forR2 = ExpressionMatrixPairFactory.connectivityAndExpression( direction, useVirtual, removeNonExp );
        // forR.writeImages();
        forR2.removeZeroConnectionRows();
        log.info( forR2.getCorrelation() );
        forR2.printDimensions();

        log
                .info( "Intersect:"
                        + Util.intersectSize( forR2.getMatrixA().getRowNames(), forR.getMatrixA().getRowNames() ) );
        log.info( "Union:" + Util.union( forR2.getMatrixA().getRowNames(), forR.getMatrixA().getRowNames() ).size() );

        System.exit( 1 );

        boolean run = false;
        forR = ExpressionMatrixPairFactory.connectivityAndExpression( direction, useVirtual, removeNonExp, run, square );
        // log.info(forR.get)
        // log.info( forR.getUniqueGenes().size() );
        // forR.printDimensions();
        // forR.printConnectionInfo();
        // System.exit( 1 );

        // Hindbrain
        // Interbrain
        // Midbrain
        // Cerebrum
        // Cerebellum

        String focusRegion = "Hindbrain";

        forR.setDivision( focusRegion );

        double sig = forR.test( 1000 );
        log.info( focusRegion + ":" + forR.getCorrelation() );
        log.info( focusRegion + " P:" + sig );
        log.info( focusRegion + " P<0.01:" + ( sig < 0.01 ) );

        // forR.test( 10000 );
        // System.exit( 1 );
        //
        // log.info( forR.getCorrelation() );
        // log.info( forR.test( 1000 ) );
        // forR.exploreSelfConnections();
        // log.info( forR.getCorrelation() );
        // log.info( forR.test( 1000 ) );
        //
        // log.info( "Diff:" + forR.bedStriaDiff( false, false ) );

        // log.info( forR.bedStriaDiff( true ) );
        // System.exit( 1 );

        // String filename = SetupParameters.getDataFolder() + "ABANonexpressed.txt";
        // List<String> rows = Util.getRowNamesFromGeneNames( filename, forR.getMatrixBDataRows() );

        // RankedGeneListLoader loader = new RankedGeneListLoader( forR.getMatrixBDataRows(),
        // SetupParameters.getDataFolder()
        // + "AllTopTen" );
        // TopTenInfo topTen = new TopTenInfo( loader, forR );
        // topTen.writeExpressionInfo();
        // System.exit( 1 );

        // direction = AnalyzeBAMSandAllenGenes.Direction.OUTGOING;
        // forR = ExpressionMatrixPairFactory.connectivityAndExpression( direction, useVirtual, removeNonExp );
        //
        // log.info( "out:" + forR.getCorrelation() );
        // log.info( "in:" + in );
        // log.info( "out:" + forR.test( 1000 ) );
        //
        // Set<String> genes = forR.getUniqueGenes();
        // log.info( "genes:" + genes.size() );
        //
        // System.exit( 1 );

        // forR.applyGeneFilter( new EstrogenGeneFilter() );
        // forR.applyGeneFilter( new EphrinGeneFilter() );
        //
        // log.info( forR.getCorrelation() );
        // // log.info( forR.getMatrixBDataRows() );
        // log.info( forR.getMatrixBDataRows().size() );
        // log.info( "genes:" + Util.getUniqueGenes( forR.getMatrixBDataRows() ) );
        // forR.test( 1000 );
        // forR.applyGeneFilter( new NaNGeneFilter( 50 ) );
        // forR.printDimensions();
        // forR.writeRMatrices();
        // log.info( forR.getCorrelation() );
        // forR.applyGeneFilter( new NaNGeneFilter() );
        // forR.printDimensions();
        // log.info( forR.getCorrelation() );

        // forR.applyGeneFilter( new NaNGeneFilter(33) );
        // forR.printDimensions();
        // log.info( forR.getCorrelation() );
        // forR.applyGeneFilter( new PlaneRemoveFilter( PlaneRemoveFilter.Plane.SAGITTAL ) );
        // forR.printDimensions();
        // log.info( forR.getCorrelation() );

        // System.exit( 1 );

        // usingOutgoinRegions();
        //
        //
        // // forR.writeRMatrices();
        // // forR.runAllenStyle();
        //
        // // log.info( forR.getCorrelation() );
        // // // log.info( "Removing coronal:" + forR.removeCoronalRows() );
        // // log.info( forR.getCorrelation() );
        //
        // System.exit( 1 );
        //
        // log.info( "Bed diff:" + forR.bedStriaDiff() );

        // List<String> ubiq = forR.getGeneRowNamesFromFile( SetupParameters.getDataFolder() + "ABAUbiquitous.txt" );
        // List<String> nonExp = forR.getGeneRowNamesFromFile( SetupParameters.getDataFolder() + "ABANonexpressed.txt"
        // );
        // List<String> both = new LinkedList<String>();
        // both.addAll( ubiq );
        // both.addAll( nonExp );
        // log.info( "Ubiquitous:" + ubiq.size() );
        // log.info( "Non expressors:" + nonExp.size() );
        // log.info( "Both:" + both.size() );
        //
        // // both.addAll( ubiq );
        // // both.addAll( nonExp );
        // log.info( forR.correlationReducedDataMatrix( ubiq ) );
        // log.info( forR.correlationReducedDataMatrix( nonExp ) );
        // log.info( forR.correlationReducedDataMatrix( both ) );
        // forR.writeRMatrices();
        //
        // forR.removeMatrixBDataRows( both );
        // log.info( forR.getCorrelation() );
        //
        // // log.info( forR.correlationReducedDataMatrix( forR.getLowVarianceDataRows( 2.5 ) ) );
        // // log.info( forR.correlationReducedDataMatrix( forR.getLowVarianceDataRows( 3 ) ) );
        //
        // // forR.test( 1000 );
        //
        // log.info( "Ttest Pvalue=" + forR.testOldHypoth() );
    }
}
