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
import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.StructureCatalogAnalyze;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.AllenDataLoaders.StructureCatalogLoader;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.BAMSDataLoaders.BAMSDataLoader;
import ubic.BAMSandAllen.BAMSDataLoaders.BAMSFromNIFDataLoader;
import ubic.BAMSandAllen.ClassSelectors.BrainRegionClassSelector;
import ubic.BAMSandAllen.FocusedAnalysis.ExploreRegionNames;
import ubic.BAMSandAllen.adjacency.CorrelationAdjacency;
import ubic.BAMSandAllen.adjacency.IdentityAdjacency;
import ubic.basecode.dataStructure.StringToStringSetMap;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.io.reader.DoubleMatrixReader;

/**
 * This class has special methods for converting Allen brain regin names to BAMS region names.
 * 
 * @author leon
 */
public abstract class ConnectivityAndAllenDataPair extends MatrixPair {
    private static Log log = LogFactory.getLog( ConnectivityAndAllenDataPair.class.getName() );

    StructureCatalogLoader allenCatalog;
    boolean squareConnectivity;
    List<String> virtualRegions;
    Direction direction;

    // matrixA is Connections
    // matrixB is Expression
    public ConnectivityAndAllenDataPair( BrainRegionClassSelector selector, boolean squareConnectivity,
            Direction direction ) throws Exception {
        this( selector, squareConnectivity, direction, false );
    }

    public ConnectivityAndAllenDataPair( String matrixFilename ) throws Exception {
        allenCatalog = new StructureCatalogLoader();

        direction = direction.ANYDIRECTION;

        DoubleMatrixReader matrixReader = new DoubleMatrixReader();
//        matrixReader.setTopLeft( false );
        DoubleMatrix<String, String> connectionMatrix = matrixReader.read( matrixFilename );

        if ( !connectionMatrix.getRowNames().equals( connectionMatrix.getColNames() ) ) {
            throw new RuntimeException( "error row and col names in diff order" );
        }
        // Remove quotes
        // assumes they are in the same order!
        LinkedList<String> newRowNames = new LinkedList<String>();
        for ( String row : connectionMatrix.getRowNames() ) {
            newRowNames.addLast( row.substring( 1, row.length() - 1 ) );
        }

        connectionMatrix.setRowNames( newRowNames );
        connectionMatrix.setColumnNames( newRowNames );

        List<String> regions = connectionMatrix.getColNames();

        // WARNING, removes occurance counts!!
        // set large propigated values to one
        for ( String row : regions ) {
            for ( String col : regions ) {
                double connectionMatrixValue = connectionMatrix.getByKeys( row, col );
                if ( connectionMatrixValue > 0 ) {
                    connectionMatrix.setByKeys( row, col, 1d );
                } else {
                    connectionMatrix.setByKeys( row, col, 0d );
                }
            }
        }

        matrixA = new ABAMSDataMatrix( connectionMatrix, "ConnectivityLoadedMatrix", new CorrelationAdjacency(
                connectionMatrix ) );

        printConnectionInfo();

        matrixA = matrixA.removeZeroColumns();

        log.info( "got Matrix A - connections" );
        virtualRegions = new LinkedList<String>();
    }

    public ConnectivityAndAllenDataPair( BrainRegionClassSelector selector, boolean squareConnectivity,
            Direction direction, boolean NIFMatrix ) throws Exception {
        allenCatalog = new StructureCatalogLoader();
        this.squareConnectivity = squareConnectivity;
        this.direction = direction;

        if ( squareConnectivity && direction.equals( Direction.APPENDED ) ) {
            throw new RuntimeException( "Error square matrix incompatable with appended connectivity" );
        }
        if ( !NIFMatrix ) {
            StructureCatalogAnalyze forMatrix = new StructureCatalogAnalyze( selector );
            // could be an option here for propigated or non-propigated connection matrix
            forMatrix.readModel( SetupParameters.getDataFolder() + "Propigated.rdf" );

            DoubleMatrix<String, String> dataMatrix = forMatrix.makeConnectionMatrix( direction );

            // matrixA = new ABAMSDataMatrix( dataMatrix, "Connectivity", new ChiIndexAdjacency() );
            // matrixA = new ABAMSDataMatrix( dataMatrix, "Connectivity", new IdentityAdjacency( dataMatrix ) );
            matrixA = new ABAMSDataMatrix( dataMatrix, "Connectivity", new CorrelationAdjacency( dataMatrix ) );
        } else {
            BAMSFromNIFDataLoader connectionLoader = new BAMSFromNIFDataLoader();
            boolean skipFibers = true;
            DoubleMatrix<String, String> dataMatrix = connectionLoader.getBAMSMatrix( direction, false, false,
                    skipFibers );
            matrixA = new ABAMSDataMatrix( dataMatrix, "Connectivity", new CorrelationAdjacency( dataMatrix ) );
        }

        log.info( "Removing bed Nuclei Stria" );
        removeBedNucleiStria();

        //
        // matrixA = matrixA.removeRows( Util.findZeroRows( matrixA ) );

        printConnectionInfo();
        //
        matrixA = matrixA.removeZeroColumns();

        log.info( "got Matrix A - connections" );
        virtualRegions = new LinkedList<String>();
    }

    public void printConnectionInfo() {
        Set<String> zeroes = Util.findZeroColumns( matrixA );
        Util.findZeroRows( matrixA );

        log.info( "Connection Zero rows(" + Util.findZeroRows( matrixA ).size() + ")" ); // :" + zeroes );
        log.info( "Connection Nonzero rows(" + ( matrixA.rows() - Util.findZeroRows( matrixA ).size() ) + ")" ); // :" +
        log.info( "Connection Zero cols(" + zeroes.size() + ")" ); // :" + zeroes );
        log.info( "Total rows(" + matrixA.rows() + ")" );

    }

    // in a bad spot, should be in constructor
    public void makeDirect() {
        matrixA = new ABAMSDataMatrix( matrixA, "Connectivity", new IdentityAdjacency( matrixA ) );

    }

    public void removeBedNucleiStria() {
        boolean inverse = false;
        boolean leaveCols = false;
        removeBedNucleiStria( inverse, leaveCols );
    }

    /**
     * removes all or retains only the connections to and from the Bed nuclei of the stria terminalis
     * 
     * @param inverse true to retain the connections, false to remove
     * @param leaveCols true to leave the columns of the matrix untouched
     */
    public void removeBedNucleiStria( boolean inverse, boolean leaveCols ) {
        BAMSDataLoader BAMSData = new BAMSDataLoader();
        boolean indirect = true;
        Set<String> bedNuclei = BAMSData.getChildren( "Bed nuclei of the stria terminalis", indirect );
        bedNuclei.add( "Bed nuclei of the stria terminalis" );

        log.info( "Bed nuclei:" + bedNuclei );
        log.info( "Bed nuclei terms:" + bedNuclei.size() );

        List<String> connectionRows = new LinkedList<String>( matrixA.getRowNames() );
        if ( inverse ) {
            connectionRows.removeAll( bedNuclei );
        } else {
            connectionRows.retainAll( bedNuclei );
        }
        log.info( "Bed nuclei rows removed:" + connectionRows.size() );
        matrixA = matrixA.removeRows( connectionRows );

        if ( !leaveCols ) {
            Set<String> connectionCols = new HashSet<String>( matrixA.getColNames() );
            if ( inverse ) {
                // do nothing
            } else {
                connectionCols.retainAll( bedNuclei );
                log.info( "Bed nuclei cols removed:" + connectionCols.size() );
                matrixA = matrixA.removeColumns( connectionCols );
                if ( isInSameSpace ) {
                    // remove from matrix B too
                    log.info( "MatrixB before:" + matrixB.columns() );
                    matrixB = matrixB.removeColumns( connectionCols );
                    log.info( "MatrixB after:" + matrixB.columns() );
                } else {
                    // if it's not in the same space then it will fail to map it to a BAMS region, and will be removed
                }
            }
        }
    }

    public void exploreSelfConnections() {
        // create self connections?
        double sum = 0;
        double count = 0;
        for ( String region : matrixA.getColNames() ) {
            if ( matrixA.getRowNames().contains( region ) ) {
                // log.info( region + " SELF:" + matrixA.getByKeys( region, region ) );

                sum += matrixA.getByKeys( region, region );
                count++;
            }
            // dataMatrix.setByKeys( region, region, 1d );
        }
        log.info( "Self connections:" + sum );
        log.info( "Possible self connections:" + count );
    }

    public void removeZeroConnectionRows() {
        matrixA = matrixA.removeRows( Util.findZeroRows( matrixA ) );
    }

    public List<String> getAllenDataColNames() {
        return new LinkedList<String>( matrixB.getColNames() );
    }

    public void printDimensions() {
        super.printDimensions();
        log.info( "Connections:" + getConnectionCount() );
    }

    /*
     * allows correlation for reduced expression/space matrix
     */
    public double correlationReducedDataMatrix( Collection<String> removeRows ) throws Exception {
        boolean fast = true;
        Set<String> removeRowsSet = new HashSet<String>( removeRows );
        ABAMSDataMatrix matrixBforCor = matrixB.removeRows( removeRowsSet );

        // log.info( matrixBName + " = " + matrixBforCor.rows() + " x " + matrixBforCor.columns() );
        return getCorrelation( matrixBforCor, fast );
    }

    public double getConnectionCount() {
        return Util.zSum( matrixA );
    }

    /**
     * Converts connectivity region name to an Allen name
     */
    public Set<String> convertANametoB( String BAMSregion ) {
        // the region has been modified to create a virtual region (the virtual mapping is 1-1)
        if ( virtualRegions.contains( BAMSregion ) ) {
            Set<String> result = new HashSet<String>();
            result.add( BAMSregion );
            return result;
        }

        Set<String> result = allenCatalog.getAllenMappedRegions( BAMSregion );

        if ( result != null ) {
            result.retainAll( allenCatalog.getLeafs() );
            // if it is in ABAMS space for both matrices, then we can't filter using this
            if ( !isInSameSpace ) result.retainAll( matrixB.getColNames() );
        }
        return result;
    }

    public Set<String> convertBNametoA( String allenRegion ) {
        boolean useVirtual = true;
        return convertBNametoA( allenRegion, useVirtual );
    }

    /**
     * Converts allen name to connectivity region name
     */
    public Set<String> convertBNametoA( String allenRegion, boolean useVirtual ) {
        if ( !allenCatalog.getLeafs().contains( allenRegion ) ) return null;

        // the region has been modified to create a virtual region
        if ( useVirtual && virtualRegions.contains( allenRegion ) ) {
            Set<String> result = new HashSet<String>();
            result.add( allenRegion );
            return result;
        }
        Set<String> result = allenCatalog.getBAMSMappedRegions( allenRegion );
        result.retainAll( matrixA.getColNames() );
        return result;
    }

    public boolean hasVirtualRegions() {
        return !virtualRegions.isEmpty();
    }

    public void makeVirtualConnectivityRegionRow( String virtualRegion ) {

    }

    public void makeVirtualConnectivityRegionColumn( String bName ) {
        // if ( squareConnectivity ) log.warn( "Virtual regions incompatable with square connectivity" );
        // put in a datastructure for convertAnametoBname
        // take the cols and OR them together
        Set<String> aNames = convertBNametoA( bName );
        // only those we have data for
        aNames.retainAll( matrixA.getColNames() );

        log.info( "Making virtual region for:" + matrixB.getName() + " " + bName + " -> " + aNames );

        virtualRegions.add( bName );
        double[] newColValues = new double[matrixA.rows()];

        for ( String aName : aNames ) {
            double[] matrixColumn = matrixA.getColumnByName( aName );
            for ( int row = 0; row < matrixColumn.length; row++ ) {
                // logical OR operation
                if ( matrixColumn[row] == 1 ) {
                    newColValues[row] = 1;
                }
            }
        }

        // remove the three cols and add a new one - use the Allen Name
        matrixA = matrixA.removeColumns( aNames );
        matrixA = matrixA.addColumn( bName, newColValues );

    }

    public void makeVirtualRegions() {
        // for each many B->many mappings
        for ( String bName : matrixB.getColNames() ) {
            Set<String> aNames = convertBNametoA( bName );
            if ( aNames != null ) {
                if ( aNames.size() > 1 ) {
                    makeVirtualConnectivityRegionColumn( bName );
                }
            }
        }

        // deal with rows
        if ( !direction.equals( Direction.APPENDED ) ) {
            DoubleMatrix<String, String> additionalRows = new DenseDoubleMatrix<String, String>( virtualRegions.size(),
                    matrixA.columns() );
            additionalRows.setRowNames( virtualRegions );
            additionalRows.setColumnNames( matrixA.getColNames() );
            for ( String virtualRegion : virtualRegions ) {
                for ( String column : additionalRows.getColNames() ) {
                    // set self loops to zero
                    if ( column.equals( virtualRegion ) ) continue;
                    // this transposes the column for the virtual region into a row
                    // since the columns no contain virtual regions we must not look those up in the rows

                    double value = 0;
                    // its a normal region
                    if ( !virtualRegions.contains( column ) ) {
                        value = matrixA.getByKeys( column, virtualRegion );
                    } else { // we have to figure out the connectivity based on the components
                        // get the regions for this virtual region and set OR the connections to this column
                        log.info( virtualRegion + " -> " + column );
                        // get all the subregions and check for connectivity
                        boolean useVirtual = false;
                        Set<String> aNames = convertBNametoA( column, useVirtual );
                        for ( String aName : aNames ) {
                            if ( matrixA.getByKeys( aName, virtualRegion ) == 1d ) {
                                value = 1;
                            }
                        }
                    }
                    additionalRows.setByKeys( virtualRegion, column, value );
                    // additionalRows.setByKeys( virtualRegion, column, 1d );

                }

            }

            matrixA = matrixA.addRows( additionalRows );
        }
    }

    /*
     * Needs to be done twice
     */
    public void removeConnectionZeroes() {
        log.info( "Removing " + Util.findZeroColumns( matrixA ).size() + " columns that have only zero values" );
        matrixA = matrixA.removeZeroColumns();
    }

    public void squareConnectionMatrix() {
        boolean removeZeroLines = true;
        squareConnectionMatrix( removeZeroLines );
    }

    public void squareConnectionMatrix( boolean removeZeroLines ) {
        // and also remove zeroes
        if ( removeZeroLines ) matrixA = removeZeroLines( matrixA );
        // make the connection matrix square
        matrixA = matrixA.getSquare();
        // and also remove zeroes
        if ( removeZeroLines ) matrixA = removeZeroLines( matrixA );

    }

    public void clean() {
        slimMatrices();
        testZeroes();
        printDimensions();
        removeConnectionZeroes();
        if ( squareConnectivity ) squareConnectionMatrix();
        if ( squareConnectivity ) squareConnectionMatrix();
        printMappingRelations();
        log.info( "After slimming" );

        printDimensions();
        testZeroes();
    }

    public double run() throws Exception {
        clean();
        sameSpace();
        log.info( "After same space" );
        printDimensions();
        double correlation = getCorrelation();
        log.info( "Correlation:" + correlation );
        return correlation;
    }

    public void removeDataRow( String row ) {
        List<String> rows = new LinkedList<String>();
        rows.add( row );
        removeMatrixBDataRows( rows );
    }

    public void orderDataRows( List<String> rows ) {
        matrixB = matrixB.orderRows( rows );
    }

    public double bedStriaDiff( boolean inverse, boolean leaveCols ) throws Exception {
        double before = getCorrelation( false );
        printDimensions();
        log.info( "Before:" + before );
        removeBedNucleiStria( inverse, leaveCols );
        printDimensions();
        double after = getCorrelation( false );
        log.info( "After:" + after );
        return before - after;
    }

    public DoubleMatrix<String, String> getConnectionMatrix( boolean square ) {
        if ( square ) {
            return matrixA.getSquare();
        } else
            return matrixA.copy();
    }

    public double runAllenStyle() throws Exception {
        log.info( "___________________________________________________________________________" );
        double correl = getCorrelation( false );
        printDimensions();
        log.info( matrixA.getName() + " and " + matrixB.getName() + " correlation is " + correl );
        testBoth( 1000, false );
        // writeImages();
        log.info( "___________________________________________________________________________" );
        return correl;
    }

    public static void main( String[] args ) throws Exception {
    }

    public void shuffleDataCols( int seed ) {
        matrixB = matrixB.shuffleColData( seed );
    }

    public void shuffleConnectivityCols( int seed ) {
        matrixA = matrixA.shuffleColData( seed );
    }

    public void removeAllenCols( Collection<String> colNames ) {
        // matrixA = matrixA.slimMatrix( colNames );
        matrixB = matrixB.retainColumns( colNames );
    }

    public List<Double> runJackKnife( boolean returnPvalues ) throws Exception {
        List<Double> resultCor = new LinkedList<Double>();
        List<Double> resultSig = new LinkedList<Double>();
        for ( String region : matrixA.getColNames() ) {
            double[] matrixAValues = matrixA.getColByName( region );
            double[] matrixBValues = matrixB.getColByName( region );

            matrixA = matrixA.removeColumn( region );
            matrixB = matrixB.removeColumn( region );

            resultCor.add( getCorrelation() );
            if ( returnPvalues ) {
                resultSig.add( test( 1000 ) );
            }
            log.info( getCorrelation() );

            // put it back
            matrixA = matrixA.addColumn( region, matrixAValues );
            matrixB = matrixB.addColumn( region, matrixBValues );
        }
        if ( returnPvalues ) return resultSig;
        return resultCor;
    }

    public Direction getDirection() {
        return direction;
    }

    public boolean isSquareConnectivity() {
        return squareConnectivity;
    }

    public void setDivision( String focusRegion ) throws Exception {
        ExploreRegionNames explore = new ExploreRegionNames( this );
        StringToStringSetMap parents = explore.getParents();
        Set<String> ROIs = parents.get( focusRegion );
        // some may have no exp
        ROIs.retainAll( getAllenDataColNames() );
        log.info( ROIs.size() );
        removeAllenCols( ROIs );
        run();
    }

}
