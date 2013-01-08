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

import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.Set;

import ubic.BAMSandAllen.adjacency.AdjacencyCompute;
import ubic.BAMSandAllen.geneFilters.GeneFilter;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;

public class ABAMSDataMatrix extends DenseDoubleMatrix<String, String> {

    String name;
    AdjacencyCompute adjacency;

    public ABAMSDataMatrix( DoubleMatrix<String, String> matrix, String name, AdjacencyCompute adjacency ) {
        super( matrix.asArray() );
        this.setRowNames( matrix.getRowNames() );
        this.setColumnNames( matrix.getColNames() );
        this.name = name;
        this.adjacency = adjacency;
        adjacency.setMatrix( this, true );
    }

    public ABAMSDataMatrix( int rows, int cols, String name, AdjacencyCompute adjacency ) {
        super( rows, cols );
        this.name = name;
        this.adjacency = adjacency;
    }

    public void setAdjacencyCompute( AdjacencyCompute adjacency ) {
        this.adjacency = adjacency;
    }

    public DoubleMatrix<String, String> getAdjacency() {
        return adjacency.getAdjacency();
    }

    public ABAMSDataMatrix applyRowFilter( GeneFilter filter ) throws Exception {
        List<String> removeRows = filter.getRowsToRemove( this );
        log.info( "Removed:" + removeRows.size() + " Name:" + filter.getName() );
        return removeRows( removeRows );
    }

    public double[] getTriangle() {
        return Util.getTriangle( adjacency.getAdjacency() );
    }

    public double[] getTriangleWithoutRow( String removed ) {
        return Util.getTriangle( adjacency.getAdjacency( removed ) );
    }

    public String getDimensionString() {
        return getName() + " = " + rows() + " x " + columns();

    }

    public AdjacencyCompute getAdjacencyCompute() {
        return adjacency;
    }

    public String getName() {
        return name;
    }

    public void setName( String name ) {
        this.name = name;
    }

    /*
     * Returns a new object with the same name and adjacency function
     */
    public ABAMSDataMatrix replaceMatrix( DoubleMatrix<String, String> newMatrix ) {
        ABAMSDataMatrix result = getReplacement( newMatrix );
        return result;
    }

    public ABAMSDataMatrix removeZeroColumns() {
        DoubleMatrix<String, String> newMatrix = removeColumns( Util.findZeroColumns( this ) );
        return replaceMatrix( newMatrix );
    }

    public ABAMSDataMatrix removeZeroRows() {
        return removeRows( Util.findZeroRows( this ) );
    }

    public ABAMSDataMatrix removeColumn( String removeColNames ) {
        Set<String> temp = new HashSet<String>();
        temp.add( removeColNames );
        return removeColumns( temp );
    }

    public ABAMSDataMatrix removeColumns( Collection<String> removeColNames ) {
        Set<String> colNames = new HashSet<String>( getColNames() );
        colNames.removeAll( removeColNames );
        return retainColumns( colNames );
    }

    /**
     * Slims matrix columns to only those in the input parameter.
     * 
     * @param colNames
     * @return
     */
    public ABAMSDataMatrix retainColumns( Collection<String> colNames ) {
        DoubleMatrix<String, String> output = new DenseDoubleMatrix<String, String>( rows(), colNames.size() );
        output.setRowNames( getRowNames() );
        output.setColumnNames( new LinkedList<String>( colNames ) );
        for ( String col : colNames ) {
            int colIndexOutput = output.getColIndexByName( col );
            int colIndexInput = getColIndexByName( col );
            for ( int row = 0; row < rows(); row++ ) {
                output.set( row, colIndexOutput, get( row, colIndexInput ) );
            }
        }
        return replaceMatrix( output );
    }

    /**
     * Adds rows to the matrix. If the same row already exists it is overwritten.
     * 
     * @param
     * @param
     * @return the new matrix
     */
    public ABAMSDataMatrix addRows( DoubleMatrix<String, String> newRows ) {
        Set<String> newRowNames = new HashSet<String>( getRowNames() );
        newRowNames.addAll( newRows.getRowNames() );
        DoubleMatrix<String, String> output = new DenseDoubleMatrix<String, String>( newRowNames.size(), columns() );
        output.setRowNames( new LinkedList<String>( newRowNames ) );
        output.setColumnNames( getColNames() );

        for ( String rowName : newRowNames ) {
            for ( String colName : output.getColNames() ) {
                double value;
                if ( newRows.containsRowName( rowName ) ) {
                    value = newRows.getByKeys( rowName, colName );
                } else {
                    value = getByKeys( rowName, colName );
                }
                output.setByKeys( rowName, colName, value );
            }
        }

        return replaceMatrix( output );
    }

    /**
     * Adds a column to the matrix
     * 
     * @param
     * @param
     * @return the new matrix
     */
    public ABAMSDataMatrix addColumn( String name, double[] values ) {
        DoubleMatrix<String, String> output = new DenseDoubleMatrix<String, String>( rows(), columns() + 1 );
        output.setRowNames( getRowNames() );
        List<String> newColnames = new LinkedList<String>( getColNames() );
        newColnames.add( name );
        output.setColumnNames( newColnames );

        for ( String col : getColNames() ) {
            int colIndexOutput = output.getColIndexByName( col );
            int colIndexInput = getColIndexByName( col );
            for ( int row = 0; row < rows(); row++ ) {
                output.set( row, colIndexOutput, get( row, colIndexInput ) );
            }
        }

        int colIndexOutput = output.getColIndexByName( name );
        for ( int row = 0; row < rows(); row++ ) {
            output.set( row, colIndexOutput, values[row] );
        }

        return replaceMatrix( output );
    }

    /*
     * assume the triangle compute is cloned
     */
    protected ABAMSDataMatrix getReplacement( DoubleMatrix<String, String> newMatrix ) {
        return getReplacement( newMatrix, true );
    }

    protected ABAMSDataMatrix getReplacement( DoubleMatrix<String, String> newMatrix, boolean clone ) {
        AdjacencyCompute newAdj = adjacency;
        if ( clone ) newAdj = adjacency.clone();
        return new ABAMSDataMatrix( newMatrix, name, newAdj );
    }

    public ABAMSDataMatrix copy() {
        AdjacencyCompute newAdj = adjacency;
        newAdj = adjacency.clone();
        return new ABAMSDataMatrix( super.copy(), name, newAdj );
    }

    /*
     * Beware: does not clone the adjacency object
     */
    public ABAMSDataMatrix copyWithoutRow( String rowNameRemove ) {
        adjacency.getAdjacency( rowNameRemove, false );

        ABAMSDataMatrix returnval = new ABAMSDataMatrix( this.rows() - 1, this.columns(), name, adjacency );
        returnval.setColumnNames( this.getColNames() );

        int skip = 0;
        for ( int i = 0, n = this.rows(); i < n; i++ ) {
            String rowName = this.getRowName( i );
            assert rowName != null : "Row " + i + " has null name";
            if ( rowName.equals( rowNameRemove ) ) {
                skip = 1;
                continue;
            }
            returnval.setRowName( rowName, i - skip );
            for ( int j = 0, m = this.columns(); j < m; j++ ) {
                returnval.set( i - skip, j, this.get( i, j ) );
            }
        }
        return returnval;
    }

    /**
     * @param args
     */
    public static void main( String[] args ) {
        // TODO Auto-generated method stub

    }

    public ABAMSDataMatrix retainRows( Collection<String> retainRows ) {
        return removeRows( ( Collection<String> ) ( Util.subtract( getRowNames(), retainRows ) ), true );
    }

    public ABAMSDataMatrix removeRows( Collection<String> removeRows ) {
        return removeRows( removeRows, true );
    }

    public ABAMSDataMatrix removeRows( Collection<String> removeRows, boolean clone ) {
        // if ( !clone ) {
        // // tell the triangle about the removed rows
        // boolean addBack = false;
        // // update the triangle, don't put the row back
        // for ( String remove : removeRows ) {
        // adjacency.getAdjacency( remove, addBack );
        // }
        // }

        ABAMSDataMatrix result = getReplacement( Util.removeRows( this, removeRows ), clone );
        boolean safe = false;
        result.adjacency.setMatrix( result, safe );
        return result;
    }

    public ABAMSDataMatrix shuffleColData( int seed ) {
        // maintains the original column name ordered and shuffles the data
        List<String> shuffledColNames = new LinkedList<String>( getColNames() );
        Random random = new Random( seed );
        Collections.shuffle( shuffledColNames, random );

        DoubleMatrix<String, String> output = new DenseDoubleMatrix<String, String>( rows(), columns() );
        output.setColumnNames( getColNames() );
        output.setRowNames( getRowNames() );

        for ( String col : getColNames() ) {
            int shuffledColIndex = shuffledColNames.indexOf( col );
            int oldColIndex = getColIndexByName( col );
            for ( int row = 0; row < rows(); row++ ) {
                output.set( row, oldColIndex, get( row, shuffledColIndex ) );
            }
        }

        boolean clone = true;
        ABAMSDataMatrix result = getReplacement( output, clone );
        return result;
    }

    public ABAMSDataMatrix orderCols( List<String> colOrder ) {
        DoubleMatrix<String, String> output = new DenseDoubleMatrix<String, String>( rows(), colOrder.size() );
        output.setRowNames( getRowNames() );
        int index = 0;
        for ( String col : colOrder ) {
            output.setColumnName( col, index++ );
            int newColIndex = output.getColIndexByName( col );
            int oldColIndex = getColIndexByName( col );
            for ( int row = 0; row < output.rows(); row++ ) {
                output.set( row, newColIndex, get( row, oldColIndex ) );
            }
        }

        boolean clone = true;
        ABAMSDataMatrix result = getReplacement( output, clone );
        boolean safe = false;
        result.adjacency.setMatrix( result, safe );
        return result;
    }

    public ABAMSDataMatrix orderRows( List<String> rowOrder ) {
        DoubleMatrix<String, String> output = new DenseDoubleMatrix<String, String>( rowOrder.size(), columns() );
        output.setColumnNames( getColNames() );
        int index = 0;
        for ( String row : rowOrder ) {
            output.setRowName( row, index++ );
            int newRowIndex = output.getRowIndexByName( row );
            int oldRowIndex = getRowIndexByName( row );
            for ( int col = 0; col < output.columns(); col++ ) {
                output.set( newRowIndex, col, get( oldRowIndex, col ) );
            }
        }

        boolean clone = true;
        ABAMSDataMatrix result = getReplacement( output, clone );
        boolean safe = false;
        result.adjacency.setMatrix( result, safe );
        return result;
    }

    public ABAMSDataMatrix getSquare() {
        Set<String> colNames = new HashSet<String>( getColNames() );
        Set<String> rowNames = new HashSet<String>( getRowNames() );
        rowNames.removeAll( colNames );
        ABAMSDataMatrix result = removeRows( rowNames, true );
        result = orderRows( getColNames() );
        return result;
    }

}
