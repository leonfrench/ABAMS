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

import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.AllenDataLoaders.AllenCatalogMatrices;
import ubic.BAMSandAllen.AllenDataLoaders.StructureCatalogLoader;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.ClassSelectors.BrainRegionClassSelector;
import ubic.BAMSandAllen.ClassSelectors.StructureDataSelector;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.graphics.MatrixDisplay;
import ubic.basecode.math.CorrelationStats;

import com.hp.hpl.jena.ontology.OntClass;

public class StructureCatalogAnalyze extends AnalyzeBAMSandAllenGenes {
    private static Log log = LogFactory.getLog( StructureCatalogAnalyze.class.getName() );

    public StructureCatalogAnalyze() {
        // ensures the BAMS regions we deal with, have a mapping to allen and one of those mappings has expression
        // information
        super( new StructureDataSelector() );
    }

    public StructureCatalogAnalyze( BrainRegionClassSelector selector ) {
        super( selector );
    }

    @Deprecated
    public void analyze() throws Exception {
        // need correlation matrix between all regions for -
        // connectivity
        // gene expression
        Direction direction = Direction.INCOMING;

        AllenCatalogMatrices allenMatrices = new AllenCatalogMatrices();
        DoubleMatrix<String, String> connectionMatrix;
        connectionMatrix = makeConnectionMatrix( direction );
        int BAMSRegionCount = connectionMatrix.columns();
        Set<String> BAMSregionNames = new HashSet<String>( connectionMatrix.getColNames() );

        DoubleMatrix<String, String> allExpressionMatrix;
        allExpressionMatrix = allenMatrices.getEnergies();

        // setup expression for BAMS regions, average when required
        StructureCatalogLoader allenCatalog = new StructureCatalogLoader();

        DoubleMatrix<String, String> expressionMatrix = new DenseDoubleMatrix<String, String>( allExpressionMatrix
                .rows(), BAMSregionNames.size() );
        expressionMatrix.setColumnNames( connectionMatrix.getColNames() );
        expressionMatrix.setRowNames( allExpressionMatrix.getRowNames() );
        log.info( "Expression matrix is " + expressionMatrix.columns() + " by " + expressionMatrix.rows() );
        // log.info( expressionMatrix.toString() );

        for ( String BAMSregion : BAMSregionNames ) {
            Set<String> allenRegions = allenCatalog.getAllenMappedRegions( BAMSregion );
            // keep Allen regions that have expression information
            allenRegions.retainAll( allExpressionMatrix.getColNames() );

            // so now these allen regions need to be moved into the BAMS region space
            for ( String allenRegion : allenRegions ) {
                int allColumnIndex = allExpressionMatrix.getColIndexByName( allenRegion );
                double[] expressionColumn = allExpressionMatrix.getColumn( allColumnIndex );

                for ( int row = 0; row < expressionColumn.length; row++ ) {
                    int columnIndex = expressionMatrix.getColIndexByName( BAMSregion );
                    double current = expressionMatrix.get( row, columnIndex );
                    // average things out
                    double value = current + expressionColumn[row] / ( ( double ) allenRegions.size() );
                    expressionMatrix.set( row, columnIndex, value );
                }
            }
        }
        // print the matrix?
        MatrixDisplay matDisplay;
        // MatrixDisplay matDisplay = new MatrixDisplay( expressionMatrix );
        // matDisplay.setLabelsVisible( true );
        // matDisplay.saveImage( SetupParameters.getDataFolder() + "expressionMatrixCatalog.png" );

        // for (String regionA)
        // allenCatalog.getAllenMappedRegions()

        // Set<String> zeroColumns = findZeroColumns( connectionMatrix );
        // zeroColumns.addAll( findZeroColumns( expressionMatrix ) );

        log.info( "Conn Zeroes:" + Util.findZeroColumns( connectionMatrix ).size() );
        log.info( "Exp Zeroes:" + Util.findZeroColumns( expressionMatrix ).size() );

        // two square matrices - size is based on BAMS mapped terms but Allen could be used (how to merge connection
        // profiles?)
        DoubleMatrix<String, String> connectionCorrelations = Util.correlateColumns( connectionMatrix, false );
        matDisplay = new MatrixDisplay( connectionCorrelations );
        matDisplay.setLabelsVisible( true );
        matDisplay.saveImage( SetupParameters.getDataFolder() + "connectionCorrelation.png" );

        DoubleMatrix<String, String> expressionCorrelations = Util.correlateColumns( expressionMatrix, false );
        matDisplay = new MatrixDisplay( expressionCorrelations );
        matDisplay.setLabelsVisible( true );
        matDisplay.saveImage( SetupParameters.getDataFolder() + "expressionCorrelation.png" );

        double[] expVecTri = Util.getTriangle( expressionCorrelations );
        double[] conVecTri = Util.getTriangle( connectionCorrelations );

        // for ( double d : expVecTri )
        // System.out.print( " E" + d );
        // System.out.println();
        // System.out.println();
        // for ( double d : conVecTri )
        // System.out.print( " C" + d );
        log.info( conVecTri[conVecTri.length - 1] );
        log.info( expVecTri[conVecTri.length - 1] );
        log.info( "Size:" + conVecTri.length );
        log.info( "Size:" + expVecTri.length );

        log.info( CorrelationStats.correl( getTriangle( expressionMatrix ), getTriangle( connectionMatrix ) ) );
        // convert to vector?
        // correlate vectors?
        List<String> names = new LinkedList<String>( expressionMatrix.getColNames() );
        Random random = new Random( 1 );
        log.info( expressionMatrix.asArray().length );
        log.info( expressionMatrix.asArray()[0].length );
        // System.exit(1);
        double real = CorrelationStats.correl( getTriangle( expressionMatrix ), getTriangle( connectionMatrix ) );
        int lowerScore = 0;
        for ( int i = 0; i < 1000; i++ ) {
            if ( i % 100 == 0 ) log.info( "Shuffle " + i );
            Collections.shuffle( names, random );
            DoubleMatrix<String, String> shuffled = new DenseDoubleMatrix<String, String>( expressionMatrix.rows(),
                    expressionMatrix.columns() );
            shuffled.setColumnNames( names );
            // ugly
            for ( int row = 0; row < shuffled.rows(); row++ ) {
                for ( int column = 0; column < shuffled.columns(); column++ ) {
                    String name = names.get( column );
                    int originalLocation = expressionMatrix.getColIndexByName( name );
                    shuffled.set( row, column, expressionMatrix.get( row, originalLocation ) );
                }
            }

            // log.info( names.get( 0 ) );
            // log.info( expressionMatrix.getColName( 0 ) );

            double shuffledResult = CorrelationStats.correl( getTriangle( shuffled ), getTriangle( connectionMatrix ) );
            if ( Math.abs( shuffledResult ) >= real ) {
                log.info( shuffledResult );
                lowerScore++;
            }
        }
        log.info( lowerScore + " random sets have correlation with absolute value above " + real );

    }

    public DoubleMatrix<String, String> shuffle( DoubleMatrix<String, String> input, Random rgen ) {
        double[][] array = input.asArray();
        // --- Shuffle by exchanging each element randomly
        // From http://leepoint.net/notes-java/algorithms/random/random-shuffling.html
        for ( int i = 0; i < array.length; i++ ) {
            int randomPosition = rgen.nextInt( array.length );
            double[] col = array[i];
            array[i] = array[randomPosition];
            array[randomPosition] = col;
        }
        return null;
    }

    public double[] getTriangle( DoubleMatrix<String, String> matrix ) {
        DoubleMatrix<String, String> correlations = Util.correlateColumns( matrix, false );
        double[] vecTri = Util.getTriangle( correlations );
        return vecTri;
    }

    public static void main( String[] args ) throws Exception {
        Direction direction = Direction.INCOMING;

        StructureCatalogAnalyze catalog = new StructureCatalogAnalyze();
        // catalog.propagateConnections();
        // instead of propigating, load a previous model
        catalog.readModel( SetupParameters.getDataFolder() + "Propigated.rdf" );

        DoubleMatrix<String, String> connectionMatrix;
        connectionMatrix = catalog.makeConnectionMatrix( direction );

        MatrixDisplay matDisplay = new MatrixDisplay( connectionMatrix );
        matDisplay.setLabelsVisible( true );
        matDisplay.saveImage( SetupParameters.getDataFolder() + "connectionMatrixCatalog.png" );

        // catalog.showMappings();

        catalog.analyze();

    }

    public String classToString( OntClass ontClass ) {
        return BAMSData.convertClassRegionToString( ontClass );
    }

}
